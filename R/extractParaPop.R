#' Extract parameter matrix for a new covariate data frame
#'
#' This is the main function for extracting parameter matrices.
#'
#' ## Parallelization options
#' If "parallel" is used then \code{mcmapply} is called from the \code{parallel} package;
#' if "biocparallel" is used, then \code{bpmapply} is called from the \code{BiocParallel} package;
#' if "future.apply" is used, then \code{future_mapply} is called from the \code{future.apply} package;
#' if "pbmcapply" is used, then \code{pbmcmapply} is called from the \code{pbmcapply} package.
#'
#' @param sce a SingleCellExperiment object.
#' @param assay_use a string scalar specifying the slot to use in input \code{sce}.
#'     The default is "counts".
#' @param marginal_list a list of named features, each with the fitted object
#'     and other variables as output from [fitMarginalPop()].
#' @param n_cores a positive integer value (greater or equal to 1) to specify the
#'     number of CPU cores used in parallelization. The default is 2.
#' @param family_use a string scalar or vector of marginal distribution used.
#' @param new_covariate a cell-by-covariate data frame obtained in the list output from
#'     [constructDataPop()]. It must have a corr_group variable.
#' @param new_eqtl_geno_list a list of eQTL genotype data frames for each gene
#'     to be simulated.  If using same list as in [fitMarginalPop()],
#'     then the in those samples
#' @param indiv_colname a string scalar of the sample ID variable in cell covariate
#'     of \code{sce}. The default is "indiv".
#' @param snp_colname a string scalar for the SNP variable in \code{eqtlgeno_df}
#'     used in [constructDataPop()]. The default is "snp_id".
#' @param loc_colname a string scalar for the last column of eQTL annotation in
#'     \code{eqtlgeno_df}. The default is "POS".
#' @param parallelization a string scalar specifying the parallelization backend
#'     used when extracting parameters. Must be one of "parallel", "future.apply",
#'     "biocparallel", or "pbmcapply". The default value is "parallel". See details.
#' @param BPPARAM a BiocParallelParam class object (from \code{BiocParallel} R package)
#'     that must be specified when using \code{parallelization = "biocparallel"}. Either
#'     \code{BiocParallel::SnowParam()} or \code{BiocParallel::MulticoreParam()}
#'     can be used to initialize, depending on the operating system. BPPARAM is
#'     not used in other parallelization options. The default is NULL.
#' @param future.seed a logical or an integer (of length one or seven), or a list
#'     of length(X) with pre-generated random seeds that can be specified when using
#'     \code{parallelization = "future.apply"}. See \code{future.apply::future_eapply}
#'     documentation for more details on its usage. future.seed is not used in
#'     other parallelization options. The default is FALSE.
#' @param data_maxsize a positive numeric value used to set max marginal_list size
#'     in GiB increments. Used only when \code{parallelization = "future.apply"}.
#'     The default is 1.
#' @param data a cell-by-covariate data frame used to fit marginal models in
#'     [fitMarginalPop()]. It must have a corr_group variable.
#' @param ... additional arguments passed to internal functions.
#'
#' @return a list of mean, sigma, and zero parameter cell by feature matrices:
#' \describe{
#'      \item{\code{mean_mat}}{a cell by feature matrix containing the conditional
#'      mean values.}
#'      \item{\code{sigma_mat}}{a cell by feature matrix containing the gene specific
#'      dispersion values.}
#'      \item{\code{zero_mat}}{a cell by feature matrix containing the gene specific
#'      zero probability values (for zip and zinb models).}
#' }
#'
#' @export
#'
#' @examples
#' NULL
extractParaPop <- function(sce,
                           assay_use = "counts",
                           marginal_list,
                           n_cores = 2L,
                           family_use,
                           new_covariate,
                           new_eqtl_geno_list,
                           indiv_colname = "indiv",
                           snp_colname = "snp_id",
                           loc_colname = "POS",
                           parallelization = c("pbmcapply", "future.apply",
                                               "parallel", "biocparallel"),
                           BPPARAM = NULL,
                           future.seed = FALSE,
                           data_maxsize = 1,
                           data,
                           ...) {

    # backward compat.
    if(parallelization == "mcmapply") {
        parallelization <- "parallel"
        message("'mcmapply' option is deprecated and will be removed in a future update.",
                "Please use 'parallel' instead.")
    } else if(parallelization == "pbmcmapply") {
        parallelization <- "pbmcapply"
        message("'pbmcmapply' option is deprecated and will be removed in a future update.",
                "Please use 'pbmcapply' instead.")
    } else if(parallelization == "bpmapply") {
        parallelization <- "biocparallel"
        message("'bpmapply' option is deprecated and will be removed in a future update.",
                "Please use 'biocparallel' instead.")
    }

    # checks
    assertthat::assert_that(assertthat::has_name(new_covariate, colnames(data)))

    checkIndividuals(eqtlgeno_list = new_eqtl_geno_list,
                     covariate_df = new_covariate,
                     eqtlgeno_name = "new_eqtl_geno_list",
                     covariate_name = "new_covariate",
                     loc_colname = loc_colname,
                     indiv_colname = indiv_colname)

    # check for new indiv
    has_newindiv <- !checkVectorEqual(unique(data[[indiv_colname]]),
                                      unique(new_covariate[[indiv_colname]]),
                                      ignore_order = TRUE)

    parallelization <- match.arg(parallelization)
    more_args <- list(...)

    removed_cell_list <- lapply(marginal_list, function(x) { x$removed_cell })
    marginal_list <- lapply(marginal_list, function(x) { x$fit })


    # find gene whose marginal is fitted
    qc_gene_idx <- which(!sapply(marginal_list, is.null))

    ## TODO: insert checks for family_use

    # check if new covariate same as train data
    data_tmp <- data[, colnames(new_covariate), drop = FALSE]
    if(identical(data_tmp, new_covariate)){
        same_cellcov <- TRUE
    } else {
        same_cellcov <- FALSE
    }

    count_mat <- SummarizedExperiment::assay(sce, assay_use)

    mat_function <- function(x, y) {
        # outputs a cell-by-3 column matrix of mean, sigma, zero parameters
        # x : feature index
        # y : model family

        fit <- marginal_list[[x]]
        removed_cell <- removed_cell_list[[x]]

        # extract response vector of cells' counts
        data$gene <- count_mat[x, ]

        total_cells <- dim(new_covariate)[1]
        cell_names <- rownames(new_covariate)

        # construct design matrix
        res_list <- constructDesignMatrix(response_vec = NULL,
                                          cellcov_df = new_covariate,
                                          eqtlgeno_df = new_eqtl_geno_list[[x]],
                                          loc_colname = loc_colname,
                                          snp_colname = snp_colname,
                                          indiv_colname = indiv_colname,
                                          filter_snps = FALSE,
                                          cleanup = FALSE)

        new_covariate <- res_list[["dmat_df"]]

        # TODO: test this code; currently skipped since removed_cell is NA for all genes
        if (length(removed_cell) > 0 && !any(is.na(removed_cell))) {

            if (is.null(new_covariate)) {
                data <- data[-removed_cell, ]

            } else {

                # extract string vector of covariates from marginal fit
                if(methods::is(fit, "gamlss")) {
                    all_covariates <- all.vars(fit$mu.formula)[-1]

                    # TODO: need to test this
                    # } else if(methods::is(fit, "glmmTMB")) {
                    #     all.vars(base::formula(fit, fixed.only = FALSE, component = "cond"))[-1]  # same as fit$call$formula

                } else {
                    all_covariates <- all.vars(fit$formula)[-1]

                }


                # find which features contain zero counts across all cells in new_covariate
                remove_idx <- lapply(
                    all_covariates, function(x) {   # TODO: figure out what this lapply does

                        curr_x <- tapply(data$gene, data[, x], sum)   # TODO: figure out what this tapply does
                        zero_group <- which(curr_x == 0)

                        if (length(zero_group) == 0) {
                            return(NA)

                        } else {
                            type <- names(curr_x)[zero_group]

                            return(which(new_covariate[, x] %in% type))
                        }
                    })

                # returns vector of indices to remove cells
                remove_cell_idx <- unlist(remove_idx)
                remove_cell_idx <- unique(stats::na.omit(remove_cell_idx))

                if(length(remove_cell_idx) > 0) {
                    new_covariate <- new_covariate[-remove_cell_idx, ]
                }
            }
        }  # end of removed_cell if statement


        # compute parameter vectors using S3 function
        param_list <- calcParaVectors(fit = fit,
                                      family_use = y,
                                      new_covariate = new_covariate,
                                      data = data,   # used in gamlss S3 method
                                      indiv_colname = indiv_colname,
                                      has_newindiv = has_newindiv,
                                      same_cellcov = same_cellcov,
                                      convert_theta = TRUE)

        mean_vec <- param_list$mean_vec
        theta_vec <- param_list$theta_vec
        zero_vec <- param_list$zero_vec

        names(zero_vec) <- names(mean_vec)

        # account for removed cells
        if (length(mean_vec) < total_cells) {

            full_means <- rep(0, total_cells)
            names(full_means) <- cell_names

            full_means[names(mean_vec)] <- mean_vec  # replace means for non-removed cells

            full_theta <- rep(NA, total_cells)
            names(full_theta) <- cell_names

            full_zero <- rep(NA, total_cells)
            names(full_zero) <- cell_names

            if (is.null(names(theta_vec))) {
                if (length(theta_vec) == length(mean_vec)) {
                    names(theta_vec) <- names(mean_vec)  # for gamlss case

                } else {
                    names(theta_vec) <- cell_names  # for gam case
                    # TODO: check for glmmTMB case
                }
            }

            full_theta[names(theta_vec)] <- theta_vec
            full_zero[names(zero_vec)] <- zero_vec

            mean_vec <- full_means
            theta_vec <- full_theta
            zero_vec <- full_zero
        }

        para_mat <- cbind(mean_vec, theta_vec, zero_vec)
        rownames(para_mat) <- cell_names

        return(para_mat)
    }  # end of mat_function()


    ## parallelize across features
    paraFunc <- parallel::mcmapply

    if (parallelization == "biocparallel") {
        paraFunc <- BiocParallel::bpmapply
        if (.Platform$OS.type == "unix") {
            BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores)
        } else if (.Platform$OS.type == "windows") {
            BPPARAM <- BiocParallel::SnowParam(workers = n_cores)
        }
    } else if (parallelization == "future.apply"){
        paraFunc <- future.apply::future_mapply
    } else if (parallelization == "pbmcapply") {
        paraFunc <- pbmcapply::pbmcmapply
    }

    if (parallelization == "biocparallel") {
        # BPPARAM$workers <- n_cores  # old code
        mat <- suppressMessages(paraFunc(mat_function,
                                         x = seq_len(dim(sce)[1])[qc_gene_idx],
                                         y = family_use,
                                         BPPARAM = BPPARAM,
                                         SIMPLIFY = FALSE))
    } else if (parallelization == "future.apply") {

        old_opt <- options("future.globals.maxSize")
        on.exit(options(old_opt), add = TRUE)
        options(future.globals.maxSize = data_maxsize * 1024^3)

        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)
        future::plan(future::multisession, workers = n_cores)

        mat <- suppressMessages(paraFunc(mat_function,
                                         x = seq_len(dim(sce)[1])[qc_gene_idx],
                                         y = family_use,
                                         future.seed = future.seed,
                                         SIMPLIFY = FALSE))
    } else {  # parallel or pbmcapply
        mat <- suppressMessages(paraFunc(mat_function,
                                         x = seq_len(dim(sce)[1])[qc_gene_idx],
                                         y = family_use,
                                         SIMPLIFY = FALSE,
                                         mc.cores = n_cores))
    }


    # create cell-by-feature mean, sigma and zero parameter matrices
    mean_mat <- sapply(mat, function(x) x[, 1])
    sigma_mat <- sapply(mat, function(x) x[, 2])  # theta renamed as sigma for scDesign3 convention
    zero_mat <- sapply(mat, function(x) x[, 3])

    #
    if (length(qc_gene_idx) < dim(sce)[1]) {  # changed in scDesign3 v0.99.7

        colnames(mean_mat) <- colnames(sigma_mat) <- colnames(zero_mat) <- rownames(sce)[qc_gene_idx]

        zeros <- matrix(0, nrow = dim(mean_mat)[1],   # from scDesign3 v0.99.7
                        ncol = dim(sce)[1] - length(qc_gene_idx))

        na_mat <- matrix(NA, nrow = dim(mean_mat)[1],
                         ncol = dim(sce)[1] - length(qc_gene_idx))

        rownames(zeros) <- rownames(mean_mat)   # from scDesign3 v0.99.7
        colnames(zeros) <- rownames(sce)[-qc_gene_idx]

        rownames(na_mat) <- rownames(mean_mat)
        colnames(na_mat) <- rownames(sce)[-qc_gene_idx]

        # pad columns at end for filtered features
        mean_mat <- cbind(mean_mat, na_mat)
        sigma_mat <- cbind(sigma_mat, na_mat)
        zero_mat <- cbind(zero_mat, na_mat)

        mean_mat <- mean_mat[, rownames(sce)]
        sigma_mat <- sigma_mat[, rownames(sce)]
        zero_mat <- zero_mat[, rownames(sce)]

    } else {
        colnames(mean_mat) <- colnames(sigma_mat) <- colnames(zero_mat) <- rownames(sce)
    }

    zero_mat <- Matrix::Matrix(zero_mat, sparse = TRUE)

    return(list(mean_mat = mean_mat,
                sigma_mat = sigma_mat,
                zero_mat = zero_mat))
}


#' Generic function to compute model parameter vectors
#'
#' A S3 generic function for computing the mean, theta, zero parameter vectors with
#' covariates for a feature, given the parametric family of the marginal model.
#'
#' It also uses either train or new data, and has option to compute parameter vectors
#' for new individuals options. Note that there is randomness introduced when
#' computing the parameter vectors for new individuals, as there are random effects.
#'
#' @inheritParams extractParaPop
#' @param fit a fitted object in the marginal_list.
#' @param ... additional arguments passed to calcParaVectors S3 method functions.
#'
#' @export
calcParaVectors <- function(fit,
                            family_use,
                            new_covariate,
                            data,
                            ...) {
    UseMethod("calcParaVectors")
}



#' A calcParaVectors Method for gamlss Objects
#'
#' @inheritParams calcParaVectors
#'
#' @exportS3Method
calcParaVectors.gamlss <- function(fit,
                                   family_use,
                                   new_covariate,
                                   data,
                                   ...) {

    more_args <- list(...)
    has_newindiv <- if(!is.null(more_args$has_newindiv)) { more_args$has_newindiv }
    same_cellcov <- if(!is.null(more_args$same_cellcov)) { more_args$same_cellcov }
    indiv_colname <- if(!is.null(more_args$indiv_colname)) { more_args$indiv_colname }
    convert_theta <- if(!is.null(more_args$convert_theta)) { more_args$convert_theta }

    mean_vec <- stats::predict(fit,
                               type = "response",
                               what = "mu",  # what option used in gamlss
                               newdata = new_covariate,
                               data = data)

    total_cells <- length(mean_vec)

    if (family_use == "poisson" | family_use == "binomial") {
        theta_vec <- rep(NA, total_cells)

    } else if (family_use == "gaussian") {
        theta_vec <- stats::predict(fit,
                                    type = "response",
                                    what = "sigma",
                                    newdata = new_covariate,
                                    data = data)

    } else if (family_use == "nb") {
        theta_vec <- stats::predict(fit,
                                    type = "response",
                                    what = "sigma",
                                    newdata = new_covariate,
                                    data = data)

    } else if (family_use == "zip") {
        theta_vec <- rep(NA, total_cells)

        zero_vec <- stats::predict(fit,
                                   type = "response",
                                   what = "sigma",
                                   newdata = new_covariate,
                                   data = data)

    } else if (family_use == "zinb") {
        theta_vec <- stats::predict(fit,
                                    type = "response",
                                    what = "sigma",
                                    newdata = new_covariate,
                                    data = data)

        zero_vec <- stats::predict(fit,
                                   type = "response",
                                   what = "nu",
                                   newdata = new_covariate,
                                   data = data)

    } else {
        stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb, zip, or zinb!")
    }

    return(list(mean_vec = mean_vec,
                theta_vec = theta_vec,
                zero_vec = zero_vec))
}


#' A calcParaVectors Method for glmmTMB Objects
#'
#' @inheritParams calcParaVectors
#'
#' @exportS3Method
calcParaVectors.glmmTMB <- function(fit,
                                    family_use,
                                    new_covariate,
                                    ...) {
    # fit : a glmmTMB fit object
    # family_use : a string scalar specifying the model parametric family used
    # new_covariate : a data frame containing the design matrix for which to compute
    #       the parameter vector
    # has_newindiv : a logical scalar for whether there are new individuals in the
    #       new_covariate that are not in the fit train data
    # same_cellcov : a logical scalar for whether the new_covariate is the same cell
    #       covariate as the fit train data. If true, then the same individuals exist
    #       as in the fit object
    # indiv_colname : a string scalar specifying the variable containing the individual
    #       IDs
    # convert_theta : a logical scalar for whether to convert theta parameter to the
    #       parameterization used in gamlss' ZINBI or NBI (type I)

    more_args <- list(...)
    has_newindiv <- if(!is.null(more_args$has_newindiv)) { more_args$has_newindiv }
    same_cellcov <- if(!is.null(more_args$same_cellcov)) { more_args$same_cellcov }
    indiv_colname <- if(!is.null(more_args$indiv_colname)) { more_args$indiv_colname }
    convert_theta <- if(!is.null(more_args$convert_theta)) { more_args$convert_theta }

    # family_use <- stats::family(fit)$family[1]
    link_type <- stats::family(fit)$link

    # if (grepl("nbinom2", family_use)) {
    #     family_use <- "nb"  # variance parameterization: var = mu + mu^2 / theta
    # }

    # detect zero-inflation
    # zi_form <- paste(as.character(fit[["modelInfo"]][["allForm"]][["ziformula"]]), collapse = "")
    # has_zi <- !is.null(zi_form) && !identical(zi_form, "~0")

    # re-define based on zero-inflation
    # family_use <- switch(
    #     family_use,
    #     "nb" = if(has_zi) { "zinb" } else { "nb" },
    #     "poisson" = if(has_zi) { "zip" } else { "poisson" },
    #     "gaussian" = "gaussian",
    #     "binomial" = "binomial",
    #     "lognormal" = "lognormal"
    # )

    if(has_newindiv) {
        mean_vec <- calcMeanVector.glmmTMB(fit = fit,
                                           new_covariate = new_covariate,
                                           has_newindiv = has_newindiv,
                                           same_cellcov = same_cellcov,
                                           indiv_colname = indiv_colname,
                                           link_type = link_type)
    } else {
        mean_vec <- calcMeanVector.glmmTMB(fit = fit,
                                           new_covariate = new_covariate,
                                           has_newindiv = has_newindiv,
                                           same_cellcov = same_cellcov,
                                           indiv_colname = indiv_colname,
                                           link_type = link_type)
    }

    total_cells <- length(mean_vec)

    # add barcode as rownames
    # names(mean_vec) <- rownames(fit$frame)  # REMOVED

    # extract dispersion and zero-inflat. prob scalar parameters
    sigma_val <- glmmTMB::sigma(fit)
    ziprob_val <- getZiProb.glmmTMB(fit)

    theta_vec <- switch(
        family_use,
        "nb" = rep(if(convert_theta) { 1 / sigma_val }    # convert for gamlss.dist::qNBI
                   else { sigma_val },
                   total_cells),
        "zinb" = rep(if(convert_theta) { 1 / sigma_val }  # convert for gamlss.dist::qZINBI
                     else { sigma_val },
                     total_cells),
        "poisson" = rep(NA, total_cells),
        "zip" = rep(NA, total_cells),
        "gaussian" = rep(sigma_val, total_cells),   # used as std. dev. in Gaussian
        "lognormal" = rep(sigma_val, total_cells),  # used as std. dev. in lognormal
        "binomial" = rep(NA, total_cells),
        stop("Distribution of glmmTMB must be one of 'nb', 'poisson', 'gaussian', 'binomial', 'zinb', 'zip', or 'lognormal'!")
        )

    # extract zero-inflation probability
    zero_vec <- switch(
        family_use,
        "nb" = rep(0, total_cells),
        "poisson" = rep(0, total_cells),
        "gaussian" = rep(0, total_cells),
        "binomial" = rep(0, total_cells),
        "zinb" = rep(ziprob_val, total_cells),
        "zip" = rep(ziprob_val, total_cells),
        "lognormal" = rep(0, total_cells),
        stop("Distribution of glmmTMB must be one of 'nb', 'poisson', 'gaussian', 'binomial', 'zinb', 'zip', or 'lognormal'!")
        )

    return(list("mean_vec" = mean_vec,
                "theta_vec" = theta_vec,
                "zero_vec" = zero_vec))
}


#' A calcParaVectors Method for gam (mgcv package) Objects
#'
#' @inheritParams calcParaVectors
#'
#' @exportS3Method
calcParaVectors.gam <- function(fit,
                                family_use,
                                new_covariate,
                                ...) {

    more_args <- list(...)
    has_newindiv <- if(!is.null(more_args$has_newindiv)) { more_args$has_newindiv }
    same_cellcov <- if(!is.null(more_args$same_cellcov)) { more_args$same_cellcov }
    indiv_colname <- if(!is.null(more_args$indiv_colname)) { more_args$indiv_colname }
    convert_theta <- if(!is.null(more_args$convert_theta)) { more_args$convert_theta }

    family_use <- stats::family(fit)$family[1]
    if (grepl("Negative Binomial", family_use)) {
        family_use <- "nb"
    }

    if (same_cellcov) {  # no new covariate

        mean_vec <- stats::predict(fit, type = "response")
        total_cells <- length(mean_vec)

        theta_vec <- switch(
            family_use,
            "nb" = rep(if(convert_theta) { 1 / fit$family$getTheta(TRUE) }    # convert for gamlss.dist::qNBI
                       else { fit$family$getTheta(TRUE) },
                       total_cells),
            "poisson" = rep(NA, total_cells),
            "gaussian" = rep(sqrt(fit$sig2), total_cells),   # this theta_vec is used for sigma_vec
            "binomial" = rep(NA, total_cells),
            stop("Distribution of gam must be one of gaussian, binomial, poisson, nb!")
        )

    } else {  # has new covariate

        mean_vec <- stats::predict(fit,
                                   type = "response",
                                   newdata = new_covariate)
        total_cells <- length(mean_vec)

        theta_vec <- switch(
            family_use,
            "nb" = rep(if(convert_theta) { 1 / fit$family$getTheta(TRUE) }    # convert for gamlss.dist::qNBI
                       else { fit$family$getTheta(TRUE) },
                       total_cells),
            "poisson" = rep(NA, total_cells),
            "gaussian" = stats::predict(fit,
                                        type = "response",
                                        what = "sigma",
                                        newdata = new_covariate),
            "binomial" = rep(NA, total_cells),
            stop("Distribution of gam must be one of gaussian, binomial, poisson, nb!")
        )
    }

    # extract zero-inflation probability
    zero_vec <- switch(
        family_use,
        "nb" = rep(0, total_cells),
        "poisson" = rep(0, total_cells),
        "gaussian" = rep(0, total_cells),
        "binomial" = rep(0, total_cells),
        stop("Distribution of gam must be one of gaussian, binomial, poisson, nb!")
    )

    return(list(mean_vec = mean_vec,
                theta_vec = theta_vec,
                zero_vec = zero_vec))
}


## computes a mean vector for cells using a glmmTMB fit
calcMeanVector.glmmTMB <- function(fit,
                                   new_covariate = NULL,
                                   has_newindiv = FALSE,
                                   same_cellcov = TRUE,
                                   indiv_colname = "indiv",
                                   link_type) {
    # fit : a glmmTMB fit object
    # new_covariate : a data frame that has the cell covariates for a feature and
    #       includes the eQTL genotypes for all the individuals
    # same_cellcov : a logical scalar for whether new_covariate has same cell
    #       covariates as "data" (ie. covariate)
    # has_newindiv : a logical scalar for whether there are new (unseen) individuals
    # link_type : a string scalar that specifies the link function used in the fit object

    if(same_cellcov) {  # no new covariate
        mean_vec <- stats::predict(fit, type = "conditional")  # uses fit$frame

    } else {  # has new covariate

        if (has_newindiv) {

            # get sigma param and individuals
            rand_sd <- getRandEffVarCov.glmmTMB(fit = fit,
                                                indiv_colname = indiv_colname)
            indivs <- unique(new_covariate[[indiv_colname]])

            # simulate random effects per indiv
            indiv_rand <- data.frame(indivs,
                                     stats::rnorm(length(indivs),
                                                  mean = 0,
                                                  sd = rand_sd))
            colnames(indiv_rand) <- c(indiv_colname, "sim_randeff")

            # Note: new random effects are not currently simulated despite using re.form = NULL
            newindiv_df <- data.frame(new_covariate[[indiv_colname]],
                                      stats::predict(fit,
                                                     type = "link",  # on g(mu) scale (eg. log scale)
                                                     newdata = new_covariate,
                                                     allow.new.levels = TRUE,
                                                     re.form = ~0))  # set random effects to zero
            colnames(newindiv_df) <- c(indiv_colname, "pop_condmean")

            newindiv_df <- dplyr::left_join(newindiv_df, indiv_rand,
                                            by = indiv_colname) %>%
                dplyr::mutate(
                    linear_pred = !!rlang::sym("pop_condmean") + !!rlang::sym("sim_randeff"),
                    mean_vec = invLinkFunc(.data$linear_pred, link_type)
                    )

            mean_vec <- newindiv_df$mean_vec

        } else {  # no new indiv
            mean_vec <- stats::predict(fit,
                                       type = "conditional",
                                       newdata = new_covariate,
                                       allow.new.levels = FALSE,
                                       re.form = NULL)  # use fitted random effects for same indivs
        }
    }

    return(mean_vec)
}


## helper function to apply inverse link function
invLinkFunc <- function(eta, link_type) {

    switch(link_type,
           "identity" = eta,
           "log"      = base::exp(eta),
           "logit"    = stats::plogis(eta),  # exp(eta) / (1 + exp(eta))
           "probit"   = stats::pnorm(eta),   # \Phi(eta)
           "inverse"  = 1 / eta,
           "sqrt"     = eta^2,
           stop(sprintf("Unsupported link function: %s!\nThe link_type must be one of 'identity', 'log', 'logit', 'probit', 'inverse', or 'sqrt'.",
                        link_type))
           )
}


## helper function to extract zero-inflation parameter from a glmmTMB object
getZiProb.glmmTMB <- function(fit) {

    zi <- glmmTMB::fixef(fit)$zi

    if(length(zi) == 0L) return(0)
    if(!identical(names(zi), "(Intercept)")) {
        stop("Currently the zero-inflation specification only supports an intercept-only model. Please check!")
    }

    return(stats::plogis(unname(zi)))
}


## helper function to extract variance-covariance matrix from a glmmTMB object
getRandEffVarCov.glmmTMB <- function(fit,
                                     indiv_colname) {

    vc <- glmmTMB::VarCorr(fit)$cond[[indiv_colname]]

    if(is.null(vc) || length(vc) == 0L) return(0)

    term_names <- colnames(vc)
    if(is.null(term_names)) {
        term_names <- rownames(vc)
    }

    if(length(term_names) != 1L || !identical(term_names, "(Intercept)")) {
        stop("Currently only a random-intercept model is supported. Please check!")
    }

    return(sqrt(vc[1, 1]))
}
