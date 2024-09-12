# TODO: test removed_cell code chunk
# TODO: try not using lapply for removed_cell_list and marginal_list

#' Extract parameter matrix for new covariate df
#'
#' This is the main function.
#'
#' @param sce add later
#' @param assay_use add later
#' @param marginal_list add later
#' @param n_cores add later
#' @param family_use a string scalar or vector of marginal distribution used.
#' @param new_covariate a cell-by-feature covariate dataframe (from construct_data.R) plus corr_group.
#' @param new_eqtl_geno_list a list of eQTL genotype dataframes for each gene (to be predicted).
#' @param indiv_colname add later
#' @param snp_colname add later
#' @param loc_colname add later
#' @param parallelization add later
#' @param BPPARAM add later
#' @param data a cell-by-feature covariate dataframe (from construct_data.R) plus corr_group. Used only in gamlss fits.
#'
#' @return a list of mean, sigma, and zero parameter cell by feature matrices:
#' @export
#'
#' @examples
#' NULL
extractParaPop <- function(sce,
                           assay_use = "counts",
                           marginal_list,
                           n_cores,
                           family_use,
                           new_covariate,
                           new_eqtl_geno_list,
                           indiv_colname = "indiv",
                           snp_colname = "snp_id",
                           loc_colname = "POS",
                           parallelization = "mcmapply",
                           BPPARAM = NULL,
                           data) {

    removed_cell_list <- lapply(marginal_list, function(x) { x$removed_cell })
    marginal_list <- lapply(marginal_list, function(x) { x$fit })

    # find gene whose marginal is fitted
    qc_gene_idx <- which(!is.na(marginal_list))

    ## TODO: insert checks for family_use and new_covariate input

    mat_function <- function(x, y) {
        # returns a cell-by-1 vector of mean, sigma, zero
        # x : feature index
        # y : model family

        fit <- marginal_list[[x]]
        removed_cell <- removed_cell_list[[x]]

        # extract vector of cells' counts stored to data
        data$gene <- SummarizedExperiment::assay(sce, assay_use)[x, ]

        if (!is.null(new_covariate)) {  # has new covariates
            total_cells <- dim(new_covariate)[1]
            cell_names <- rownames(new_covariate)

            ### previous code chunk
            # construct the new_covariate
            # eqtl_geno_df <- new_eqtl_geno_list[[x]]

            # convert to sample by using indiv label
            # geno_df <- eqtl_geno_df %>%
            #     dplyr::select(c(!!rlang::sym(snp_colname), new_covariate[[indiv_colname]])) %>%
            #     tidyr::pivot_longer(., cols = -snp_colname,
            #                         values_to = "genotype",
            #                         names_to = indiv_colname) %>%
            #     tidyr::pivot_wider(., names_from = snp_colname, values_from = "genotype")


            # construct design matrix df for marginal fitting
            # new_covariate <- as.data.frame(new_covariate) %>%
            #     dplyr::select(-!!rlang::sym(indiv_colname), dplyr::everything()) %>%
            #     dplyr::left_join(., geno_df, by = indiv_colname) %>%
            #     dplyr::mutate(!!rlang::sym(indiv_colname) := as.factor(!!rlang::sym(indiv_colname)))

            # rownames(new_covariate) <- cell_names
            ###

            new_covariate_list <- constructDesignMatrix(response_vec = NULL,
                                                        cellcov_df = new_covariate,
                                                        eqtlgeno_df = new_eqtl_geno_list[[x]],
                                                        loc_colname = loc_colname,
                                                        snp_colname = snp_colname,
                                                        indiv_colname = indiv_colname,
                                                        filter_snps = FALSE,
                                                        cleanup = FALSE)

            new_covariate <- new_covariate_list[["dmat_df"]]
        }

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
                remove_idx <- lapply(all_covariates,
                                     function(x) {   # TODO: figure out what this lapply does

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


        # use S3 generic function to compute parameter vectors
        param_list <- calcParaVectors(fit = fit,
                                      family_use = y,
                                      new_covariate = new_covariate,
                                      data = data,   # only used in gamlss S3 method
                                      total_cells = total_cells)

        mean_vec <- param_list$mean_vec
        theta_vec <- param_list$theta_vec

        if ("zero_vec" %in% names(param_list)) {
            zero_vec <- param_list$zero_vec
        } else {
            zero_vec <- rep(0L, length(mean_vec))
        }
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
    paraFunc <- parallel::mcmapply  # set parallelization function

    if (parallelization == "bpmapply") {
        paraFunc <- BiocParallel::bpmapply
    }
    if (parallelization == "pbmcmapply") {
        paraFunc <- pbmcapply::pbmcmapply
    }

    if (parallelization == "bpmapply") {
        BPPARAM$workers <- n_cores

        # loop thru each feature using mat_function()
        mat <- suppressMessages(paraFunc(mat_function,
                                         x = seq_len(dim(sce)[1])[qc_gene_idx],
                                         y = family_use,
                                         BPPARAM = BPPARAM,
                                         SIMPLIFY = FALSE))
    } else {  # run mcmapply or pbmcmapply
        mat <- suppressMessages(paraFunc(mat_function,
                                         x = seq_len(dim(sce)[1])[qc_gene_idx],
                                         y = family_use,
                                         SIMPLIFY = FALSE,
                                         mc.cores = n_cores))
    }


    # create cell-by-feature mean, sigma and zero parameter matrices
    mean_mat <- sapply(mat, function(x) x[, 1])
    sigma_mat <- sapply(mat, function(x) x[, 2])
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
#' A S3 generic function for computing model parameters for with or without new
#'     covariate of a feature
#'
#' @param fit add later
#' @param family_use add later
#' @param new_covariate add later
#' @param total_cells add later
#' @param data add later
#' @param ... add later
#'
#' @export
calcParaVectors <- function(fit,
                            family_use,
                            new_covariate,
                            total_cells,
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
                                   total_cells,
                                   data,
                                   ...) {

    mean_vec <- stats::predict(fit,
                               type = "response",
                               what = "mu",  # what option used in gamlss
                               newdata = new_covariate,
                               data = data)

    if (family_use == "poisson" | family_use == "binomial") {
        theta_vec <- rep(NA, total_cells)

    } else if (family_use == "gaussian") {
        theta_vec <- stats::predict(fit,
                                    type = "response",  # "response" is on the scale of the response variable
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
                                   what = "sigma",  # should be "nu" instead?
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
                                    total_cells,
                                    ...) {

    family_use <- stats::family(fit)$family[1]
    if (grepl("nbinom2", family_use)) {
        family_use <- "nb"  # variance parameterization: var = mu + mu^2 / phi
    }

    if (is.null(new_covariate)) {  # no new covariate

        mean_vec <- stats::predict(fit, type = "response")  # use train data stored in fit$frame

        if (family_use == "poisson" | family_use == "binomial") {
            theta_vec <- rep(NA, total_cells)

        } else if (family_use == "gaussian") {
            theta <- glmmTMB::sigma(fit)  # used as sigma for Gaussian
            theta_vec <- rep(theta, total_cells)

        } else if (family_use == "nb") {
            theta <- glmmTMB::sigma(fit)
            theta_vec <- 1 / rep(theta, total_cells)  # convert to NB type I in gamlss.dist::qNBI

        } else {
            stop("Distribution of glmmTMB must be one of gaussian, binomial, poisson, nb!")
        }

    } else {  # has new covariate

        mean_vec <- stats::predict(fit,
                                   type = "response",
                                   newdata = new_covariate)
        # TODO: modify to allow new levels in new_eqtl_geno

        if (family_use == "poisson" | family_use == "binomial") {
            theta_vec <- rep(NA, total_cells)

        } else if (family_use == "gaussian") {
            theta <- glmmTMB::sigma(fit)
            theta_vec <- rep(theta, total_cells)

        } else if (family_use == "nb") {
            theta <- glmmTMB::sigma(fit)
            theta_vec <- 1 / rep(theta, total_cells)  # convert to NB type I in gamlss.dist::qNBI

        } else {
            stop("Distribution of glmmTMB must be one of gaussian, binomial, poisson, nb!")
        }
    }

    return(list(mean_vec = mean_vec,
                theta_vec = theta_vec))
}


#' A calcParaVectors Method for gam (mgcv package) Objects
#'
#' @inheritParams calcParaVectors
#'
#' @exportS3Method
calcParaVectors.gam <- function(fit,
                                family_use,
                                new_covariate,
                                total_cells,
                                ...) {

    family_use <- stats::family(fit)$family[1]
    if (grepl("Negative Binomial", family_use)) {
        family_use <- "nb"
    }

    if (is.null(new_covariate)) {  # no new covariate

        mean_vec <- stats::predict(fit, type = "response")

        if (family_use == "poisson" | family_use == "binomial") {
            theta_vec <- rep(NA, total_cells)

        } else if (family_use == "gaussian") {
            theta_vec <- rep(sqrt(fit$sig2), total_cells)  # this theta_vec is used for sigma_vec

        } else if (family_use == "nb") {
            theta <- fit$family$getTheta(TRUE)
            theta_vec <- 1 / rep(theta, total_cells)  # convert to NB type I in gamlss.dist::qNBI

        } else {
            stop("Distribution of gam must be one of gaussian, binomial, poisson, nb!")
        }

    } else {  # has new covariate

        mean_vec <- stats::predict(fit,
                                   type = "response",
                                   newdata = new_covariate)

        if (family_use == "poisson" | family_use == "binomial") {
            theta_vec <- rep(NA, total_cells)

        } else if (family_use == "gaussian") {
            theta_vec <- stats::predict(fit,
                                        type = "response",
                                        what = "sigma",
                                        newdata = new_covariate)

        } else if (family_use == "nb") {
            theta <- fit$family$getTheta(TRUE)
            theta_vec <- 1 / rep(theta, total_cells)  # convert to NB type I in gamlss.dist::qNBI

        } else {
            stop("Distribution of gam must be one of gaussian, binomial, poisson, or nb!")
        }
    }

    return(list(mean_vec = mean_vec,
                theta_vec = theta_vec))
}
