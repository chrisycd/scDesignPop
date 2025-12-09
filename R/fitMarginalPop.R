
#' Fit marginal models for every feature
#'
#' Fits a specified parametric model using various input parameters.
#'
#' @param data_list a list of input data.
#' @param mean_formula a string scalar to specify the mean formula, including
#'     random effects (if any) but without SNP genotypes or SNP genotype interaction effects.
#' @param model_family a string scalar to specify model fitting used.
#' @param interact_colnames a string scalar or vector for the variable names that
#'     have first-order interaction with SNP genotypes.
#' @param parallelization a string scalar specifying the type of parallelization
#'     used during marginal fitting. Must be one of either "future.apply",
#'     or "pbmcmapply". The default value is "pbmcmapply".
#' @param n_threads positive integer value (greater or equal to 1) to specify the
#'     number of CPU threads used in parallelization. The default is 2.
#' @param loc_colname a string scalar for column name of SNP position variable.
#' @param snp_colname a string scalar for column name of SNP id variable.
#' @param cellstate_colname a string scalar for column name of cell state (ie. cell type).
#' @param indiv_colname a string scalar for column name of individuals (samples).
#' @param filter_snps a logical scalar for whether to filter out SNP covariates
#'     with either low-variance or with only 1 distinct genotype (ie. all 1's)
#'     prior to fitting the model.
#' @param snpvar_thres a numeric scalar (between 0 to 1) used to filter out SNPs
#'     whose variance of genotypes across samples are below this threshold.
#'     Used together when \code{filter_snps = TRUE}.
#' @param force_formula a logical scalar for whether to bypass model parsimony check.
#'     If \code{force_formula = TRUE}, interaction terms whose covariates are not
#'     main effects in the model would be permitted. Results in error if
#'     \code{force_formula = FALSE} and \code{length(geno_interact_names) > 0}.
#' @param data_maxsize a positive numeric value used to set max data_list size
#'     in GiB increments. Used only when \code{parallelization = "future.apply"}.
#' @param keep_cellnames a logical scalar for whether to keep cell barcode names.
#'     If \code{keep_cellnames = TRUE}, the memory will be larger. The default is FALSE.
#'
#' @return a list of named features, each containing a list with the following items:
#' \describe{
#'      \item{\code{fit}}{a \code{glmmTMB} fit object.}
#'      \item{\code{time}}{a numeric scalar for the elapsed time to fit the given feature.}
#'      \item{\code{snp_cov}}{a string scalar or vector of SNP ids used in the
#'          fit for given feature.}
#'      \item{\code{model_attr}}{a list of attributes extracted from each model. NOT
#'          currently implemented.}
#'      \item{\code{removed_cell}}{a string scalar or vector of the cell names
#'          removed due to low-variance (NOT currently implemented).}
#' }
#' @export
#'
#' @examples
#' NULL
fitMarginalPop <- function(data_list,
                           mean_formula,
                           model_family = "nb",
                           interact_colnames = NULL,
                           parallelization = "pbmcapply",
                           n_threads = 2L,
                           loc_colname = "POS",
                           snp_colname = "snp_id",
                           cellstate_colname = "cell_type",
                           indiv_colname = "indiv",
                           filter_snps = TRUE,
                           snpvar_thres = 0,
                           force_formula = FALSE,
                           data_maxsize = 1,
                           keep_cellnames = FALSE) {

    # check cell covariates
    assertthat::assert_that(assertthat::has_name(data_list[["covariate"]],
                                                 c(interact_colnames,
                                                   indiv_colname,
                                                   cellstate_colname)))

    assertthat::assert_that(assertthat::has_name(data_list[["new_covariate"]],
                                                 c(interact_colnames,
                                                   indiv_colname,
                                                   cellstate_colname)))

    if(is.null(data_list[["eqtl_geno_list"]])) {
        has_eqtl <- FALSE
    } else {
        has_eqtl <- TRUE
    }

    # check features
    important_features <- data_list[["important_features"]]
    noeqtl_features <- data_list[["noeqtl_features"]]  # ADDED
    sce_features <- colnames(data_list[["count_mat"]])
    cell_names <- rownames(data_list[["count_mat"]])

    if(has_eqtl) {
        eqtl_features <- names(data_list[["eqtl_geno_list"]])

        stopifnot("Features in count_mat and eqtl_geno_list do not match. Please check input!" =
                      checkVectorEqual(sce_features, eqtl_features, ignore_order = TRUE))

        stopifnot("Features in important_features and eqtl_geno_list do not match. Please check input!" =
                      checkVectorEqual(important_features, eqtl_features, ignore_order = TRUE))
    }


    stopifnot("Features in important_features and count_mat do not match. Please check input!" =
                  checkVectorEqual(important_features, sce_features, ignore_order = TRUE))

    # check eQTL geno
    if(has_eqtl) {
        eqtl_colnames <- lapply(data_list[["eqtl_geno_list"]], function(x) colnames(x)) %>%
            Reduce(intersect, .)

        stopifnot("loc_colname, snp_colname, or cellstate_colname is missing in eqtl_geno_list. Please check input!" =
                      checkVectorContain(c(loc_colname, snp_colname, cellstate_colname), eqtl_colnames))
    }


    # remove cell barcodes
    if(!keep_cellnames) {
        rownames(data_list[["count_mat"]]) <- NULL
        rownames(data_list[["covariate"]]) <- NULL
    }

    ## run marginal fitting

    if(parallelization == "pbmcapply") {
        marginal_list <- pbmcapply::pbmclapply(
            sce_features,
            function(x) {
                fitModel(feature_name = x,
                         response_vec = data_list[["count_mat"]][, x],
                         eqtlgeno_df = data_list[["eqtl_geno_list"]][[x]],
                         cellcov_df = data_list[["covariate"]],
                         mu_formula = mean_formula,
                         interact_colnames = interact_colnames,
                         model_family = model_family,
                         loc_colname = loc_colname,
                         snp_colname = snp_colname,
                         cellstate_colname = cellstate_colname,
                         indiv_colname = indiv_colname,
                         filter_snps = filter_snps,
                         snpvar_thres = snpvar_thres,
                         force_formula = force_formula)
                }, mc.cores = n_threads)

    } else if(parallelization == "future.apply") {

        old_opt <- options("future.globals.maxSize")
        on.exit(options(old_opt), add = TRUE)
        options(future.globals.maxSize = data_maxsize * 1000 * 1024^2)  # set max size for data_list

        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)
        future::plan(future::multisession, workers = n_threads)

        marginal_list <- future.apply::future_lapply(
            sce_features,
            function(x) {
                fitModel(feature_name = x,
                         response_vec = data_list[["count_mat"]][, x],
                         eqtlgeno_df = data_list[["eqtl_geno_list"]][[x]],
                         cellcov_df = data_list[["covariate"]],
                         mu_formula = mean_formula,
                         interact_colnames = interact_colnames,
                         model_family = model_family,
                         loc_colname = loc_colname,
                         snp_colname = snp_colname,
                         cellstate_colname = cellstate_colname,
                         indiv_colname = indiv_colname,
                         filter_snps = filter_snps,
                         snpvar_thres = snpvar_thres,
                         force_formula = force_formula)
                }
            )
    }
    # TODO: add progress bar to future.apply (see https://github.com/HenrikBengtsson/future.apply/issues/34#issuecomment-549011124)

    # store gene ids as list names
    names(marginal_list) <- sce_features

    #
    # marginal_list[["cell_names"]] <- cell_names

    # clean up
    rm(data_list)
    gc()

    return(marginal_list)
}


#' Fit a Marginal Model
#'
#' Fits a specified parametric model for a feature using a response variable, and
#'     eQTL genotype and cell covariates as explanatory variables.
#'
#' @inheritParams fitMarginalPop
#' @param feature_name a string scalar of a feature's name (ie. gene id).
#' @param response_vec a vector of values for the response variable.
#' @param cellcov_df a cell-by-covariates dataframe containing the covariates
#'     (explanatory variables) for all cells.
#' @param eqtlgeno_df a SNP-by-sample genotype dataframe containing a feature's
#'     eQTL annotations and SNP genotypes (explanatory variables) for all samples
#'     (ie. individuals).
#' @param mu_formula a string scalar to specify the mean formula, including random
#'     effects (if any) but without SNP genotypes or SNP genotype interaction effects.
#'
#' @return a list containing the fitted model object, elapsed time, SNP ids of covariates,
#'     and removed cells (NOT currently implemented.)
#' @export
#'
#' @examples
#' NULL
fitModel <- function(feature_name,
                     response_vec,
                     cellcov_df,
                     eqtlgeno_df,
                     mu_formula,
                     model_family = "nb",
                     interact_colnames = NULL,
                     loc_colname = "POS",
                     snp_colname = "snp_id",
                     cellstate_colname = "cell_type",
                     indiv_colname = "indiv",
                     filter_snps = TRUE,
                     snpvar_thres = 0,
                     force_formula = FALSE) {
    # TODO: fix warning msg for using external vector in selections was deprecated in tidyselect
    #       use `all_of()` or `any_of()` instead

    if(is.null(eqtlgeno_df)) {
        has_eqtl <- FALSE
    } else {
        has_eqtl <- TRUE
    }

    # checks
    if(has_eqtl) {
        assertthat::assert_that(assertthat::has_name(eqtlgeno_df,
                                                     c(loc_colname, snp_colname, cellstate_colname)))
    }

    assertthat::assert_that(assertthat::has_name(cellcov_df,
                                                 c(interact_colnames, indiv_colname, cellstate_colname)))

    stopifnot("response_vec must be either numeric or factor. Please check input!" =
                  (class(response_vec) == "numeric" || class(response_vec) == "factor"))


    removed_cell <- NA  # TODO: add later

    # create design matrix and filtered SNP covariate vector
    res_list <- constructDesignMatrix(response_vec = response_vec,
                                      cellcov_df = cellcov_df,
                                      eqtlgeno_df = eqtlgeno_df,
                                      loc_colname = loc_colname,
                                      snp_colname = snp_colname,
                                      indiv_colname = indiv_colname,
                                      filter_snps = filter_snps,
                                      snpvar_thres = snpvar_thres,
                                      cleanup = TRUE)

    snp_cov <- res_list[["snp_cov"]]

    # create model formula
    model_formula <- constructFormula(model_formula = mu_formula,
                                      interact_colnames = interact_colnames,
                                      snp_cov = snp_cov,
                                      force_formula = force_formula)

    # fit marginal
    start_time <- Sys.time()

    glmmTMB.fit <- tryCatch({
        # code from https://stackoverflow.com/questions/68084740/r-trycatch-but-retain-the-expression-result-in-the-case-of-a-warning
        expr = {
            withCallingHandlers(
                expr = {

                    if(model_family == "nb") {  # GLMM NB
                        model <- glmmTMB::glmmTMB(formula = model_formula,
                                                  data = res_list[["dmat_df"]],
                                                  family = glmmTMB::nbinom2,
                                                  ziformula = ~0)
                    } else if(model_family == "poisson") {  # GLMM Poisson
                        model <- glmmTMB::glmmTMB(formula = model_formula,
                                                  data = res_list[["dmat_df"]],
                                                  family = stats::poisson,
                                                  ziformula = ~0)
                    } else if(model_family == "gaussian") {  # LMM
                        model <- glmmTMB::glmmTMB(formula = model_formula,
                                                  data = res_list[["dmat_df"]],
                                                  family = stats::gaussian,
                                                  ziformula = ~0)
                    } else if(model_family == "binomial") {
                        model <- glmmTMB::glmmTMB(formula = model_formula,
                                                  data = res_list[["dmat_df"]],
                                                  family = stats::binomial,
                                                  ziformula = ~0)
                    } else {
                        stop("The model_family option is not valid.\n
                             Must be one of 'nb', 'poisson', 'gaussian', or 'binomial'.")
                    }

                    # Note: use ziformula = ~1 for zero-inflation

                    # model[["frame"]] <- c()  # remove training dataframe
                    model_attr <- list("class" = attr(model[["frame"]], "class"),
                                       "terms" = attr(model[["frame"]], "terms"),
                                       "names" = attr(model[["frame"]], "names"))
                    # model[["frame"]] <- model[["frame"]] %>%    # retain response vec and snp covariate df
                    #     dplyr::select(c("response", dplyr::all_of(snp_cov)))

                    model
                },
                # If expression throws a warning, record diagnostics without halting,
                # so as to store the result of the expression.
                warning = function(w){

                    message(sprintf("glmmTMB fit issues on feature %s: %s \n", feature_name, base::conditionMessage(w)))

                    # parent <- parent.env(environment())
                    # parent$diag <- w

                }
            )
        }

    }, error = function(e) {

        message(sprintf("glmmTMB fit fails on feature %s: %s \n", feature_name, base::conditionMessage(e)))
        NULL

    }, silent = FALSE)

    end_time <- Sys.time()
    elap_time <- as.numeric(end_time - start_time)

    # clean-up
    rm(res_list, cellcov_df, response_vec)
    gc()

    return(list("fit" = glmmTMB.fit,
                "time" = elap_time,
                "snp_cov" = snp_cov,  # snp covariates used in model fitting
                "model_attr" = model_attr,
                "removed_cell" = removed_cell)
           )
}



#' Construct a Design Matrix Dataframe
#'
#' A helper function that constructs a design matrix using a given feature's
#'     expression vector for every cell, eQTL genotype dataframe (optional),
#'     and cell covariate dataframe.
#'
#' When \code{response_vec} is provided, the function constructs a full design
#'     matrix dataframe (ie. with response variable), whereas only a covariate
#'     dataframe is constructed when \code{response_vec = NULL}.
#'
#' @inheritParams fitMarginalPop
#' @inheritParams fitModel
#' @param cleanup a logical scalar for whether to clean up variables after
#'     constructing \code{dmat_df}.
#'
#' @return a list containing the following:
#' \describe{
#'      \item{\code{dmat_df}}{a dataframe of the design matrix containing all covariates
#'      (both cell covariates and eQTL genotype covariates) for a given feature.}
#'      \item{\code{snp_cov}}{a string scalar or vector of SNP ids in the design matrix.}
#' }
#' @export
#'
#' @examples
#' NULL
constructDesignMatrix <- function(response_vec,
                                  cellcov_df,
                                  eqtlgeno_df,
                                  loc_colname = "POS",
                                  snp_colname = "snp_id",
                                  indiv_colname = "indiv",
                                  filter_snps = TRUE,
                                  snpvar_thres = 0,
                                  cleanup = TRUE) {

    if(!is.null(eqtlgeno_df)) {

        # construct a sample by snp df using indiv label
        firstsnp_indx <- which(colnames(eqtlgeno_df) == loc_colname) + 1L  # last column preceding genotypes
        lastsnp_indx <- dim(eqtlgeno_df)[2]

        geno_df <- eqtlgeno_df[, firstsnp_indx:lastsnp_indx] %>%
            dplyr::bind_cols(!!as.name(snp_colname) := eqtlgeno_df[[snp_colname]], .) %>%
            # eqtlgeno_df %>%   # alternate code w/o using loc_colname variable
            # dplyr::select(tidyselect::all_of(
            #     c(snp_colname, as.character(unique(new_covariate[[indiv_colname]]))
            #     ))) %>%
            tidyr::pivot_longer(., cols = -snp_colname,
                                values_to = "genotype",
                                names_to = indiv_colname) %>%
            tidyr::pivot_wider(., names_from = snp_colname, values_from = "genotype")

        snp_cov <- eqtlgeno_df[[snp_colname]]

        # filter out low variance SNPs or single unique genotype SNPs (ie. all 1's)
        if(filter_snps) {

            if_zerovar <- geno_df[, base::setdiff(colnames(geno_df), indiv_colname)] %>%
                sapply(., stats::var) <= snpvar_thres

            if_oneunique <- geno_df[, base::setdiff(colnames(geno_df), indiv_colname)] %>%
                sapply(., function(x) length(unique(x))) == 1

            exclude_snps <- snp_cov[which(if_zerovar | if_oneunique)]

            # remove excluded SNP variables and update remaining SNP covariates
            if(length(exclude_snps) != 0L) {
                geno_df <- geno_df %>%
                    dplyr::select(-dplyr::all_of(exclude_snps))  # TODO: check this

                snp_cov <- snp_cov[!(if_zerovar | if_oneunique)]
            }
        }
    } else {   # w/o eqtl
        snp_cov <- NULL
    }

    # construct design matrix df for marginal fitting
    # dmat_df <- data.frame("response" = response_vec) %>%
    #     dplyr::bind_cols(., cellcov_df) %>%
    dmat_df <- cellcov_df %>%  # keeps row names
        { if(!is.null(eqtlgeno_df)) dplyr::left_join(., geno_df, by = indiv_colname) else . } %>%  # check eqtlgeno
        dplyr::mutate(dplyr::across(tidyselect::where(is.character), as.factor))
    # dplyr::mutate(!!rlang::sym(indiv_colname) := as.factor(!!rlang::sym(indiv_colname)))  # old code


    # add response variable to design matrix
    if(!is.null(response_vec)) {
        dmat_df <- dmat_df %>%
            dplyr::bind_cols(data.frame("response" = response_vec), .)
    }


    # clean up
    if(cleanup) {
        rm(cellcov_df)
        gc()
    }

    return(list("dmat_df" = dmat_df,
                "snp_cov" = snp_cov))
}



#' Constructs a Model Formula
#'
#' A helper function that specifies the formula for the marginal model.
#'
#' @inheritParams fitMarginalPop
#' @param model_formula a string scalar to specify the model formula, including
#'     random effects (if any) but without SNP genotypes or SNP genotype interaction effects.
#' @param snp_cov a string scalar or vector that specifies the SNP covariates used in the model.
#'
#' @return a formula object
#' @export
#'
#' @examples
#' NULL
constructFormula <- function(model_formula,
                             interact_colnames = NULL,
                             snp_cov = NULL,
                             force_formula = FALSE) {

    # check model parsimony
    if(!force_formula && (length(interact_colnames) > 0L)) {
        if(!all(base::sapply(interact_colnames, base::grepl, model_formula))) {
            stop(sprintf("Not all interact_colnames variables are included in the model_formula. Please check!"))
        }
    }

    # add SNP covariates to model formula
    if(length(snp_cov) > 0) {  # if not empty

        snp_cov <- paste0("`", snp_cov, "`")  # add backticks if snp has colon
        snp_terms <- paste(snp_cov, collapse = " + ")

        model_formula <- paste0(model_formula, " + ", snp_terms)
    }

    # add SNP interaction terms
    if(!is.null(interact_colnames) && (length(snp_cov) > 0)) {

        interact_terms <- expand.grid(snp_cov, interact_colnames, sep = ":") %>%  # all interaction combinations
            base::with(., paste(Var1, Var2, sep = ":")) %>%   # concat together
            paste(., collapse = " + ")  # combine strings
        # TODO: resolve notes for Var1, Var2 in devtools::check()

        model_formula <- paste0(model_formula, " + ", interact_terms)
    }

    model_formula <- stats::as.formula(paste("response ~ ", model_formula))

    return(model_formula)
}

