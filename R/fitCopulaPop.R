# TODO: missing convert_u function

#' Fits copula for input
#'
#' This is the main function for fitting a copula.
#'
#' ## Parallelization options
#' If "parallel" is used then \code{mcmapply} is called from the \code{parallel} package; if
#' "biocparallel" is used, then \code{bpmapply} is called from the \code{BiocParallel} package; if
#' "future.apply" is used, then \code{future_mapply} is called from the \code{future.apply} package;
#' if "pbmcapply" is used, then \code{pbmcmapply} is called from the \code{pbmcapply} package.
#'
#' @param sce a SingleCellExperiment object.
#' @param assay_use a string scalar specifying the slot to use in input \code{sce}.
#'     The default is "counts".
#' @param input_data a cell-by-covariate data frame obtained in the list output from
#'     \code{\link{constructDataPop}}. It must have a corr_group variable.
#' @param marginal_list a list of named features, each with the fitted object
#'     and other variables as output from \code{\link{fitMarginalPop}}.
#' @param family_use a string scalar to specify model fitting used.
#' @param copula a string value to specify the type of Copula fitting to use.
#'     Currently, only Gaussian copula is supported. The default is "gaussian".
#' @param DT a logic value to perform the distributional transformation. The default
#'     is \code{TRUE}.
#' @param pseudo_obs a logic value. If TRUE, use the empirical quantiles instead
#'     of theoretical quantiles for fitting copula. The default is FALSE.
#'     (not currently implemented)
#' @param epsilon a numeric value close to 0 to specify tolerance for avoiding
#'     0 or 1 quantiles. The default value is 1e-06.
#' @param family_set a string or a string vector of the bivariate copula families.
#'     Default is c("gaussian", "indep").
#' @param important_feature a string vector of genes.
#' @param n_cores positive integer value (greater or equal to 1) to specify the
#'     number of CPU cores used in parallelization. The default is 2.
#' @param parallelization a string value specifying the parallelization backend
#'     used during copula fitting. Must be one of "parallel", "future.apply",
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
#' @param ... additional arguments passed to internal functions.
#'
#' @return outputs a list with following elements:
#' \describe{
#'      \item{\code{model_aic}}{total model AIC value.}
#'      \item{\code{model_bic}}{total model BIC value.}
#'      \item{\code{copula_list}}{a list of fitted copulas for each \code{corr_group}.}
#'      \item{\code{important_features}}{a string vector of genes.}
#'      \item{\code{na_marginal_aic}}{a named string vector of genes that did not
#'          have AIC in the marginal models.}
#'      \item{\code{na_marginal_bic}}{a named string vector of genes that did not
#'          have BIC in the marginal models.}
#' }
#'
#' @import glmmTMB
#' @export
#'
#' @examples
#' NULL
fitCopulaPop <- function(sce,
                         assay_use,
                         input_data,   # cell covariate from construct_data w/ corr_group variable
                         # empirical_quantile = FALSE,  option not implemented at moment
                         # new_covariate = NULL,  # removed in scDesign3 v0.99.7
                         marginal_list,
                         family_use,
                         copula = c("gaussian", "vine"),
                         DT = TRUE,  # distributional transformation
                         pseudo_obs = FALSE,
                         epsilon = 1e-06,  # tolerance for avoiding 0 or 1 quantiles
                         family_set = c("gaussian", "indep"),
                         important_feature = "all",
                         n_cores = 2L,
                         parallelization = c("pbmcapply", "future.apply",
                                             "parallel", "biocparallel"),
                         BPPARAM = NULL,
                         future.seed = FALSE,
                         data_maxsize = 1,
                         ...) {

    # backward compat.
    if(parallelization == "mcmapply") {
        parallelization <- "parallel"
    } else if(parallelization == "pbmcmapply") {
        parallelization <- "pbmcapply"
    } else if(parallelization == "bpmapply") {
        parallelization <- "biocparallel"
    }

    copula <- match.arg(copula)
    parallelization <- match.arg(parallelization)
    more_args <- list(...)

    if(important_feature == "all") {
        important_feature <- rep(TRUE, dim(sce)[1])
    }

    marginals <- lapply(marginal_list, function(x) { x$fit })

    # find gene whose marginal is fitted
    qc_gene_idx <- which(!is.na(marginals))



    # NOTE: AIC and BIC computation tmp not implemented
    # NOTE: empirical_quantile chunk code tmp not implemented
    # NOTE: if ind chunk code tmp not implemented

    group_index <- unique(input_data$corr_group)
    ind <- group_index[1] == "ind"


    if (copula == "gaussian") {
        message("Convert Residuals to Multivariate Gaussian")
        newmat <- convert_n_pop(sce = sce[qc_gene_idx, ],
                                assay_use = assay_use,
                                marginal_list = marginal_list[qc_gene_idx],
                                data = input_data,
                                DT = DT,
                                pseudo_obs = pseudo_obs,
                                n_cores = n_cores,
                                family_use = family_use,
                                epsilon = epsilon,
                                parallelization = parallelization,
                                BPPARAM = BPPARAM,  # added in scDesign3 v0.99.7
                                future.seed = future.seed,
                                data_maxsize = data_maxsize)
        message("Converting End")

    } else {  # not implemented
        message("Convert Residuals to Uniform")
        newmat <- convert_u(sce = sce[qc_gene_idx, ],  # convert_u not implemented
                            assay_use = assay_use,
                            marginal_list = marginal_list[qc_gene_idx],
                            data = input_data,
                            DT = DT,
                            pseudo_obs = pseudo_obs,
                            family_use = family_use,
                            n_cores = n_cores,
                            epsilon = epsilon,
                            parallelization = parallelization,
                            BPPARAM = BPPARAM)  # added in scDesign3 v0.99.7
        message("Converting End")
    }

    ## select important genes
    if (is.vector(important_feature) & methods::is(important_feature, "logical")) {
        if (length(important_feature) != dim(sce)[1]) {
            stop("The important_feature should either be 'auto' or a logical vector with the length equals to the number of genes in the input data")
        }
    } else {
        if (important_feature == "auto") {
            gene_zero_prop <- apply(as.matrix(SummarizedExperiment::assay(sce, assay_use)), 1,
                                    function(y) { sum(y < 1e-05) / dim(sce)[2] })
            # TODO: add option for sparse calculation if sparse sce matrix is used

            important_feature <- gene_zero_prop < 0.8  # default zero prop in scDesign2
            names(important_feature) <- rownames(sce)

        } else {
            stop("The important_feature should either be 'auto' or a logical vector with the length equals to the number of genes in the input data")
        }
    }

    important_feature <- important_feature[qc_gene_idx]

    # group_index <- unique(input_data$corr_group)  # removed in scDesign3 v0.99.7
    # ind <- group_index[1] == "ind"                # removed in scDesign3 v0.99.7

    corr_group <- as.data.frame(input_data$corr_group)
    colnames(corr_group) <- "corr_group"

    ## removed in scDesign3 v0.99.7
    # if (is.null(new_covariate)) {
    #     new_corr_group <- NULL
    # } else {
    #     new_corr_group <- as.data.frame(new_covariate$corr_group)
    #     colnames(new_corr_group) <- "corr_group"
    # }
    ##

    newmvn.list <- lapply(group_index, function(x,
                                                sce,
                                                newmat,
                                                corr_group,
                                                # new_corr_group,  # removed in scDesign3 v0.99.7
                                                ind,
                                                n_cores,
                                                important_feature) {
        # compute covariance matrix for each corr_group

        message(paste0("Copula group ", x, " starts"))
        curr_index <- which(corr_group[, 1] == x)

        ## removed in scDesign3 v0.99.7
        # if (is.null(new_covariate)) {
        #     curr_ncell <- length(curr_index)
        #     curr_ncell_idx <- curr_index

        # } else {
        #     curr_ncell <- length(which(new_corr_group[, 1] == x))
        #     curr_ncell_idx <- paste0("Cell", which(new_corr_group[, 1] == x))

        # }
        ##

        # skipped empirical_quantile code

        if (copula == "gaussian") {
            curr_mat <- newmat[curr_index, , drop = FALSE]
            cor.mat <- cal_cor(curr_mat, important_feature = important_feature,
                               if.sparse = FALSE, lambda = 0.05, tol = 1e-08,
                               ind = ind)
            model_aic <- cal_aic(norm.mat = newmat, cor.mat = cor.mat,
                                 ind = ind)
            model_bic <- cal_bic(norm.mat = newmat, cor.mat = cor.mat,
                                 ind = ind)

        } else if (copula == "vine") {
            if (!ind) {
                message("Vine Copula Estimation Starts")
                start <- Sys.time()
                curr_mat <- newmat[curr_index, , drop = FALSE]
                curr_mat <- curr_mat[, which(important_feature)]
                vine.fit <- rvinecopulib::vinecop(data = curr_mat,
                                                  family_set = family_set, show_trace = FALSE,
                                                  par_method = "mle", cores = n_cores)
                end <- Sys.time()
                print(end - start)

                message("Vine Copula Estimation Ends")
                model_aic <- stats::AIC(vine.fit)
                model_bic <- stats::BIC(vine.fit)
                cor.mat <- vine.fit

            } else {

                model_aic <- 0
                model_bic <- 0
                cor.mat <- NULL
            }
        } else {

            stop("Copula must be one of 'vine' or 'gaussian'")
        }    # end of copula if statement

        return(list("model_aic" = model_aic,
                    "model_bic" = model_bic,
                    # "copula_list" = copula_list,              # added in scDesign3 v0.99.7
                    # "important_feature" = important_feature,  # added in scDesign3 v0.99.7
                    "cor.mat" = cor.mat)  # cor.mat not returned in scDesign v0.99.7
        )
    },
    sce = sce,
    newmat = newmat,
    ind = ind,
    n_cores = n_cores,
    corr_group = corr_group,
    # new_corr_group = new_corr_group,   # removed in scDesign3 v0.99.7
    important_feature = important_feature)  # end of newmvn.list <- lapply()

    # compute AIC/BIC
    marginal.aic <- aggregateMarginalAIC(model_list = marginals,
                                         indx_vec = qc_gene_idx,
                                         type = "AIC")
    marginal.bic <- aggregateMarginalAIC(model_list = marginals,
                                         indx_vec = qc_gene_idx,
                                         type = "BIC")

    copula.aic <- sum(sapply(newmvn.list, function(x) x$model_aic))
    # marginal.aic <- sum(sapply(marginals[qc_gene_idx], stats::AIC))

    copula.bic <- sum(sapply(newmvn.list, function(x) x$model_bic))
    # marginal.bic <- sum(sapply(marginals[qc_gene_idx], stats::BIC))

    # model_aic <- c(marginal.aic, copula.aic, marginal.aic + copula.aic)
    model_aic <- c(marginal.aic$total, copula.aic, marginal.aic$total + copula.aic)
    names(model_aic) <- c("aic.marginal", "aic.copula", "aic.total")

    # model_bic <- c(marginal.bic, copula.bic, marginal.bic + copula.bic)
    model_bic <- c(marginal.bic$total, copula.bic, marginal.bic$total + copula.bic)
    names(model_bic) <- c("bic.marginal", "bic.copula", "bic.total")

    copula_list <- lapply(newmvn.list, function(x) x$cor.mat)
    names(copula_list) <- group_index

    return(list(model_aic = model_aic,
                model_bic = model_bic,
                copula_list = copula_list,
                important_feature = important_feature,
                na_marginal_aic = marginal.aic$na_vec,
                na_marginal_bic = marginal.bic$na_vec))
}



## Convert marginal distributions to standard normals
convert_n_pop <- function(sce,
                          assay_use,
                          marginal_list,
                          data,
                          DT = TRUE,
                          pseudo_obs = FALSE,
                          epsilon = 1e-06,
                          family_use,
                          n_cores,
                          parallelization,
                          BPPARAM,   # Note: BPPARAM option added in scDesign v0.99.7
                          future.seed,
                          data_maxsize = 1,
                          ...) {

    more_args <- list(...)

    ## Extract gene-by-cell count matrix
    # count_mat <- t(as.matrix(SummarizedExperiment::assay(sce, assay_use)))  # old code
    count_mat <- SummarizedExperiment::assay(sce, assay_use)

    # get removed cells
    removed_cell_list <- lapply(marginal_list, function(x) { x$removed_cell })
    marginal_list <- lapply(marginal_list, function(x) { x$fit })

    # ncell <- dim(count_mat)[1]
    ncell <- dim(count_mat)[2]  # total cells


    mat_function <- function(x, y) {
        # x : feature index
        # y : model family

        fit <- marginal_list[[x]]
        removed_cell <- removed_cell_list[[x]]

        if(length(removed_cell) > 0 && !any(is.na(removed_cell))) {
            data <- data[-removed_cell, ]
        }

        if(methods::is(fit, "gamlss")) {  # gamlss currently NOT IMPLEMENTED
            # mean_vec <- stats::predict(fit, type = "response", what = "mu", data = data)

            # code chunk skipped for rest of gamlss models

        } else if(methods::is(fit, "glmmTMB")) {
            # if input is from glmmTMB

            # model type for marginal fit of y response
            y <- stats::family(fit)$family[1]

            if (grepl("nbinom2", y)) {  # var = mu + mu^2 / phi
                y <- "nb"
            }

            mean_vec <- stats::predict(fit, type = "response")  # TODO: use option newdata

            # add barcode as rownames
            # names(mean_vec) <- rownames(fit$frame)  # REMOVED

            # extract dispersion parameter
            if (y == "poisson" | y == "binomial") {
                theta_vec <- rep(NA, length(mean_vec))

            } else if (y == "gaussian") {
                theta <- glmmTMB::sigma(fit)
                theta_vec <- rep(theta, length(mean_vec))  # used as sigma for Gaussian

            } else if (y == "nb") {
                theta <- glmmTMB::sigma(fit)
                theta_vec <- rep(theta, length(mean_vec))

            } else {
                stop("Distribution of glmmTMB must be one of gaussian, poisson, binomial or nb!")
            }

        } else {  # if input is from mgcv

            # model type for marginal fit of y response
            y <- stats::family(fit)$family[1]

            if (grepl("Negative Binomial", y)) {
                y <- "nb"
            }

            mean_vec <- stats::predict(fit, type = "response")

            # add barcode as rownames
            # names(mean_vec) <- rownames(fit$frame)  # REMOVED

            # extract dispersion parameter
            if (y == "poisson" | y == "binomial") {
                theta_vec <- rep(NA, length(mean_vec))

            } else if (y == "gaussian") {
                theta_vec <- rep(sqrt(fit$sig2), length(mean_vec))  # called the theta_vec but actually used as sigma_vec for Gaussian

            } else if (y == "nb") {
                theta <- fit$family$getTheta(TRUE)
                theta_vec <- rep(theta, length(mean_vec))
            } else {
                stop("Distribution of mgcv must be one of gaussian, poisson or nb!")
            }
        }  # end of mgcv


        ## Mean for Each Cell
        # Y <- count_mat[names(mean_vec), x]   # REVISED
        Y <- count_mat[x, ]

        ## Frame
        if (!exists("zero_vec")) {
            zero_vec <- 0
        }

        family_frame <- cbind(Y, mean_vec, theta_vec, zero_vec)

        if (y == "binomial") {
            pvec <- apply(family_frame, 1, function(x) {
                stats::pbinom(x[1], prob = x[2], size = 1)
            })
        } else if (y == "poisson") {
            pvec <- apply(family_frame, 1, function(x) {
                stats::ppois(x[1], lambda = x[2])
            })
        } else if (y == "gaussian") {
            pvec <- apply(family_frame, 1, function(x) {
                gamlss.dist::pNO(x[1], mu = x[2], sigma = abs(x[3]))
            })
        } else if (y == "nb") {
            pvec <- apply(family_frame, 1, function(x) {
                stats::pnbinom(x[1], mu = x[2], size = x[3])
            })
        } else if (y == "zip") {
            pvec <- apply(family_frame, 1, function(x) {
                gamlss.dist::pZIP(x[1], mu = x[2], sigma = abs(x[4]))
            })
        } else if (y == "zinb") {
            pvec <- apply(family_frame, 1, function(x) {
                gamlss.dist::pZINBI(x[1],
                                    mu = x[2],
                                    sigma = abs(x[3]),
                                    nu = x[4])
            })
        } else {
            stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb, zip or zinb!")
        }


        ## CHECK ABOUT THE FIRST PARAM!!!!!
        if (DT) {  # distribution transform
            if (y == "poisson") {
                pvec2 <- apply(family_frame, 1, function(x) {
                    stats::ppois(x[1] - 1, lambda = x[2]) * as.integer(x[1] > 0)
                })
            } else if (y == "gaussian" | y == "binomial") {
                ## Gaussian is continuous, thus do not need DT.
                ## Binomial seems to be weird to have DT.
                message("Continuous gaussian does not need DT.")
                pvec2 <- pvec

            } else if (y == "nb") {

                pvec2 <- apply(family_frame, 1, function(x) {
                    stats::pnbinom(x[1] - 1, mu = x[2], size = x[3]) * as.integer(x[1] > 0)
                })
            } else if (y == "zip") {
                pvec2 <- apply(family_frame, 1, function(x) {
                    ifelse(x[1] > 0, gamlss.dist::pZIP(x[1] - 1, mu = x[2], sigma = abs(x[4])), 0)
                })
            } else if (y == "zinb") {
                pvec2 <- apply(family_frame, 1, function(x) {
                    ifelse(x[1] > 0,
                           gamlss.dist::pZINBI(x[1] - 1,
                                               mu = x[2],
                                               sigma = abs(x[3]),
                                               nu = x[4]
                                               ),
                           0)
                })
            } else {
                stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb, zip or zinb!")
            }

            u1 <- pvec
            u2 <- pvec2

            v <- stats::runif(length(mean_vec))

            ## Random mapping
            r <- u1 * v + u2 * (1 - v)
        } else {
            r <- pvec
        }

        if (length(r) < dim(sce)[2]) {
            new_r <- rep(1, dim(sce)[2])
            names(new_r) <- colnames(sce)
            new_r[names(r)] <- r
            r <- new_r
        }

        ## Avoid Inf
        idx_adjust <- which(1 - r < epsilon)
        r[idx_adjust] <- r[idx_adjust] - epsilon

        idx_adjust <- which(r < epsilon)
        r[idx_adjust] <- r[idx_adjust] + epsilon

        if (pseudo_obs) {
            r <- rvinecopulib::pseudo_obs(r)
        }

        stats::qnorm(r,
                     mean = 0,
                     sd = 1,
                     lower.tail = TRUE,
                     log.p = FALSE)

    }  # end of mat_function()

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
        mat <- paraFunc(mat_function,
                        x = seq_len(dim(sce)[1]),
                        y = family_use,
                        SIMPLIFY = TRUE,
                        BPPARAM = BPPARAM)
    } else if (parallelization == "future.apply") {

        old_opt <- options("future.globals.maxSize")
        on.exit(options(old_opt), add = TRUE)
        options(future.globals.maxSize = data_maxsize * 1000 * 1024^2)

        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)
        future::plan(future::multisession, workers = n_cores)

        mat <- paraFunc(mat_function,
                        x = seq_len(dim(sce)[1]),
                        y = family_use,
                        future.seed = future.seed,
                        SIMPLIFY = TRUE)
    } else {  # parallel or pbmcapply
        mat <- paraFunc(mat_function,
                        x = seq_len(dim(sce)[1]),
                        y = family_use,
                        SIMPLIFY = TRUE,
                        mc.cores = n_cores)
    }
    colnames(mat) <- rownames(sce)
    rownames(mat) <- colnames(sce)

    ## Remove inf
    mat[is.infinite(mat)] <- NA

    ## Use mean to replace the missing value
    for (i in 1:ncol(mat)) {
        mat[is.na(mat[, i]), i] <- mean(mat[, i], na.rm = TRUE)
    }

    # output quantile matrix
    return(mat)
}


## Calculate the correlation matrix
# Note: sparse cor estimation not implemented
cal_cor <- function(norm.mat,
                    important_feature,
                    if.sparse = FALSE,
                    lambda = 0.05,
                    tol = 1e-08,
                    ind = FALSE) {

    if (ind) {  # uncorrelated genes
        cor.mat <- diag(rep(1, dim(norm.mat)[2]))

        return(cor.mat)
    } else {

        cor.mat <- diag(rep(1, dim(norm.mat)[2]))
        rownames(cor.mat) <- colnames(norm.mat)
        colnames(cor.mat) <- colnames(norm.mat)

        important.mat <- norm.mat[, which(important_feature)]
        important_cor.mat <- correlation(important.mat)

        #s_d <- apply(norm.mat, 2, stats::sd)
        s_d <- matrixStats::colSds(important.mat, na.rm = TRUE)

        if (any(0 == s_d)) {  # zero std dev
            important_cor.mat[is.na(important_cor.mat)] <- 0
        }

        cor.mat[rownames(important_cor.mat), colnames(important_cor.mat)] <- important_cor.mat
    }

    n <- dim(cor.mat)[1]

    cor.mat
}


# TODO: Try fix issue when cor.mat is not pos. definite and dmvnorm() returns -Inf
## Calculate Gaussian copula AIC
cal_aic <- function(norm.mat,
                    cor.mat,
                    ind) {
    if (ind) {
        copula.aic <- 0
    } else {
        copula.nop <- (as.integer(sum(cor.mat != 0)) - dim(cor.mat)[1]) / 2

        copula.aic <- -2 * (sum(mvtnorm::dmvnorm(x = norm.mat,
                                                 mean = rep(0, dim(cor.mat)[1]),
                                                 sigma = cor.mat,
                                                 log = TRUE)
        ) - sum(rowSums(stats::dnorm(norm.mat, log = TRUE)))) + 2 * copula.nop
    }

    copula.aic
}


## Calculate Gaussian copula BIC
cal_bic <- function(norm.mat,
                    cor.mat,
                    ind) {

    n_obs <- dim(norm.mat)[1]
    if (ind) {
        copula.bic <- 0
    } else {
        copula.nop <- (as.integer(sum(cor.mat != 0)) - dim(cor.mat)[1]) / 2

        copula.bic <- -2 * (sum(mvtnorm::dmvnorm(x = norm.mat,
                                                 mean = rep(0, dim(cor.mat)[1]),
                                                 sigma = cor.mat,
                                                 log = TRUE)
        ) - sum(rowSums(stats::dnorm(norm.mat, log = TRUE)))) + log(n_obs) * copula.nop
    }

    copula.bic
}


## Similar to the cora function from "Rfast" but uses different functions to calculate column means and row sums
correlation <- function(x) {
    mat <- t(x) - matrixStats::colMeans2(x)
    mat <- mat / sqrt(matrixStats::rowSums2(mat^2))
    tcrossprod(mat)  # is this from Matrix package?
}


## aggregates computed information criterion (AIC or BIC) for a list of marginal fitted models
aggregateMarginalAIC <- function(model_list,
                               indx_vec,
                               type = c("AIC", "BIC"),
                               na_rm = TRUE) {
    # model_list : a list of fitted models
    # indx_vec : a vector of indices to select elements from model_list
    # na_rm : a logical value for whether to remove marginal models with NA
    #         computed AIC or BIC. The default is TRUE.

    type <- match.arg(type)

    if (type == "AIC") {
        computeFunc <- stats::AIC
    } else if (type == "BIC") {
        computeFunc <- stats::BIC
    } else {
        NULL
    }

    value_vec <- sapply(model_list[indx_vec], computeFunc)

    na_vec <- which(is.na(value_vec))

    if(length(na_vec) == 0) {
        na_vec <- NULL
    }

    total <- sum(value_vec, na.rm = na_rm)

    return(list("total" = total,
                "na_vec" = na_vec))
}
