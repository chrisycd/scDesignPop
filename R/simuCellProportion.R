#' Simulate cell proportions using genotype principal components and population-level covariates
#'
#' Function fits a multinomial regression model using the input genotype principal
#' components (PCs) and other population-level covariates as training data and outputs
#' simulated cell proportions.
#'
#' @param sce a SingleCellExperiment object.
#' @param genoPC a data frame of individual by genotype principal components for \code{sce} input.
#'     The first column must be the variable for individual with same name as \code{indiv_colname}.
#' @param new_genoPC a data frame of individual by genotype principal components for
#'     simulated individuals. The first column must be the variable for individual
#'     with same name as \code{indiv_colname}, followed by "PC1", "PC2", etc.
#' @param new_othercov a data frame of the test data containing same additional covariates
#'     as in colData of \code{sce}.
#' @param PCnum an integer scalar specifying the number of principal components used
#'     in multinomial regression.
#' @param cov_colnames an optional string vector or scalar for the variable names to include
#'     in the cell proportion model. Variables must exist in both \code{new_othercov} and
#'     colData of \code{sce}.
#' @param indiv_colname a string scalar to specify the variable in \code{sce} containing
#'     individuals.
#' @param cellstate_colname a string scalar to specify the variable in \code{sce} containing
#'     cell states (ie. cell types).
#' @param cn_model_family a string scalar to specify the model family used for total
#'     cell modeling. Currently only 'lognormal' from \link[MASS]{fitdistr} is supported.
#' @param cn_meanlog a numeric scalar for the mean parameter (on log scale) of the total
#'     cell number model. When \code{cn_meanlog = NULL}, the parameter is estimated from
#'     input data.
#' @param cn_sdlog a numeric scalar for the standard deviation parameter (on log scale)
#'     of the total cell number model. When \code{n_sdlog = NULL}, the parameter is
#'     estimated from input data.
#' @param cp_model_family a string scalar to specify the model family used for cell
#'     proportion modeling. Currently only 'MN' from \link[MGLM]{dist} is supported.
#' @param cp_intercept a logical scalar for whether to include an intercept in the
#'     cell proportion model.
#' @param ... additional optional arguments.
#'
#' @return outputs a list with following elements:
#' \describe{
#'      \item{\code{simu_cov}}{a cell-by-covariate data frame of simulated cell types and corresponding individual.}
#'      \item{\code{cp_simu_df}}{a cell type-by-covariate data frame summarizing the simulate cell proportions, total cell numbers, and cells per cell types.}
#'      \item{\code{cp_modelfit}}{a fitted model object for the cell proportion model.}
#'      \item{\code{cn_modelfit}}{a fitted model object for the cell number model.}
#' }
#'
#' @export
#'
#' @examples
#' NULL
simuCellProportion <- function(sce,
                               genoPC,
                               new_genoPC,
                               new_othercov,
                               PCnum = 5L,
                               cov_colnames = NULL,
                               indiv_colname = "indiv",
                               cellstate_colname = "cell_type",
                               cn_model_family = "lognormal",
                               cn_meanlog = NULL,
                               cn_sdlog = NULL,
                               cp_model_family = "MN",
                               cp_intercept = TRUE,
                               ...) {

    # TODO: make genoPC, new_genoPC, and PCnum so multinomial can be fit without
    #       additional covariates

    opt_args <- list(...)

    # coerce to dataframes
    genoPC <- as.data.frame(genoPC)
    new_genoPC <- as.data.frame(new_genoPC)
    new_othercov <- as.data.frame(new_othercov)

    sce_colnames <- colnames(SummarizedExperiment::colData(sce))

    # checks
    assertthat::assert_that(assertthat::has_name(new_othercov, indiv_colname))

    stopifnot("Input new_genoPC does not have first column named as indiv_colname. Please check input!" =
                  (colnames(new_genoPC)[1] == indiv_colname))

    stopifnot("Input genoPC does not have first column named as indiv_colname. Please check input!" =
                  (colnames(genoPC)[1] == indiv_colname))

    stopifnot("Input new_genoPC and new_othercov should have same individuals. Please check input!" =
                  checkVectorEqual(new_genoPC[[indiv_colname]], new_othercov[[indiv_colname]], ignore_order = FALSE))

    stopifnot("Input sce does not have matching indiv_colname, cellstate_colname, or cov_colnames variables. Please check input!" =
                  checkVectorContain(c(indiv_colname, cellstate_colname, cov_colnames), sce_colnames))

    stopifnot("Input genoPC and new_genoPC should have same column variables. Please check input!" =
                  checkVectorEqual(colnames(genoPC), colnames(new_genoPC), ignore_order = FALSE))

    stopifnot("Input new_othercov does not have all the column variables specified by cov_colnames. Please check input!" =
                  checkVectorContain(cov_colnames, colnames(new_othercov)))


    # extract cell type and indiv dataframe from sce input
    frame <- as.data.frame(SummarizedExperiment::colData(sce)[, c(cellstate_colname, indiv_colname)])

    # extract population-level covariates from sce input
    orig_othercov <- unique(SummarizedExperiment::colData(sce)[, c(indiv_colname, cov_colnames), drop = FALSE])
    orig_othercov <- as.data.frame(orig_othercov)


    # fit cell proportion model using train data
    cp_modelfit <- fitCellPropModel(frame = frame,
                                    genoPC = genoPC,
                                    PCnum = PCnum,
                                    cp_model_family = cp_model_family,
                                    othercov = orig_othercov,
                                    indiv_colname = indiv_colname,
                                    cellstate_colname = cellstate_colname,
                                    intercept = cp_intercept)

    # simulate cell prop for test data
    X_new <- simuCellProp(genoPC = new_genoPC,
                          PCnum = PCnum,
                          othercov = new_othercov,
                          indiv_colname = indiv_colname,
                          intercept = cp_intercept)

    cp_simu <- MGLM::predict(cp_modelfit, newdata = X_new)

    # create df for simulated cell prop
    cp_simu_df <- data.frame(as.vector(cp_simu),
                                 rep(rownames(cp_simu), dim(cp_simu)[2]),
                                 rep(colnames(cp_simu), each = dim(cp_simu)[1]),
                                 "simulated")

    colnames(cp_simu_df) <- c("prop", indiv_colname, cellstate_colname, "group")


    # fit lognormal model on train data
    cn_modelfit <- fitCellTotalModel(frame = frame,
                                     cn_model_family = cn_model_family)

    # check for cell num input parameters
    if(is.null(cn_meanlog)) {
        cn_meanlog <- cn_modelfit[["estimate"]][["meanlog"]]
    }

    if(is.null(cn_sdlog)) {
        cn_sdlog <- cn_modelfit[["estimate"]][["sdlog"]]
    }

    # simulate total cells per indiv in test data
    cn_simu <- stats::rlnorm(dim(cp_simu)[1],  # nindiv
                                 meanlog = cn_meanlog,
                                 sdlog = cn_sdlog)
    # TODO: add check if any cn_simu is below threshold or negative, and display message then resample

    # add cell numbers per cell type
    cp_simu_df[["celltotalnum"]] <- rep(cn_simu, dim(cp_simu)[2])
    cp_simu_df[["cellnum"]] <- round(cp_simu_df[["celltotalnum"]] * cp_simu_df[["prop"]])

    # expand simulated cell types for each indiv by cellnum value
    simu_cov <- cp_simu_df[rep(1:nrow(cp_simu_df), cp_simu_df[["cellnum"]]), ] %>%
        dplyr::select(c(!!rlang::sym(cellstate_colname), !!rlang::sym(indiv_colname)))
    rownames(simu_cov) <- paste0("simcell", 1:nrow(simu_cov), "_", simu_cov[[indiv_colname]])
    # TODO: use cell index per indiv and cell type

    # convert all character columns to factor variables
    ischar_cols <- base::sapply(simu_cov, is.character)
    simu_cov[ischar_cols] <- base::lapply(simu_cov[ischar_cols], as.factor)


    return(list("simu_cov" = simu_cov,
                "cp_simu_df" = cp_simu_df,
                "cp_modelfit" = cp_modelfit,
                "cn_modelfit" = cn_modelfit))
}


#### helper functions

# fit cell proportion model
fitCellPropModel <- function(frame,
                             genoPC,
                             PCnum,
                             cp_model_family,
                             othercov,
                             indiv_colname,
                             cellstate_colname,
                             intercept,
                             ...) {

    opt_args <- list(...)

    Y <- unclass(as.matrix(t(table(frame))))
    X <- constructCellPropDesignMatrix(frame,
                                       genoPC,
                                       PCnum,
                                       othercov,
                                       indiv_colname,
                                       cellstate_colname,
                                       intercept)
    X <- X[rownames(Y), ]


    # create formula and design matrix
    formula <- paste("~", paste(colnames(X), collapse = " + "))

    if(intercept) {
        formula <- stats::as.formula(formula)
    } else {
        formula <- stats::as.formula(paste(formula, "-1"))
    }
    # TODO: fix no intercept issue

    X <- stats::model.matrix(formula, data = X)

    # fit multinomial regression
    return(MGLM::MGLMreg.fit(Y = Y, X = X, dist = cp_model_family))
}


# combine frame with principal components
constructCellPropDesignMatrix <- function(frame,
                                          genoPC,
                                          PCnum,
                                          othercov,
                                          indiv_colname,
                                          cellstate_colname,
                                          intercept,
                                          ...) {

    opt_args <- list(...)

    # merge with geno pcs and other covariates to create indiv by covariate design matrix
    X <- dplyr::distinct(frame, !!rlang::sym(indiv_colname),
                         !!rlang::sym(cellstate_colname), .keep_all = FALSE) %>%
        dplyr::left_join(genoPC[, 1:(1 + PCnum)], by = indiv_colname) %>%
        # TODO: remove hard-coding by 1
        dplyr::left_join(othercov, by = indiv_colname) %>%
        dplyr::select(-!!rlang::sym(cellstate_colname)) %>%
        dplyr::distinct(indiv, .keep_all = TRUE) %>%
        tibble::column_to_rownames(var = indiv_colname)

    # X <- data.matrix(X)

    # if(intercept) {
    #     X <- cbind(rep(1, nrow(X)), X)
    # }

    return(X)
}


# simulate cell proportion
simuCellProp <- function(genoPC,
                         PCnum,
                         othercov,
                         indiv_colname,
                         intercept) {

    X <- genoPC[, 1:(1 + PCnum)] %>%
        dplyr::left_join(othercov, by = indiv_colname) %>%
        tibble::column_to_rownames(var = indiv_colname)

    X <- data.matrix(X)
    # X <- stats::model.matrix(X)

    # adding intercept column
    if(intercept) {
        X <- cbind(rep(1, nrow(X)), X)
    }

    return(X)
}


# fit model for total cells across individuals
fitCellTotalModel <- function(frame,
                              cn_model_family) {

    Y_tot <- rowSums(unclass(as.matrix(t(table(frame)))))
    cellnum_fit <- MASS::fitdistr(Y_tot, densfun = cn_model_family)

    return(cellnum_fit)
}
