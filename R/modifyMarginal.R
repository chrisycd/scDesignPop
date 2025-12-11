#' Modify marginal models
#'
#' Modify the marginal model parameters of genes based on user inputs for
#'     cell-type-specific eQTLs.
#'
#' @param marginal_list A list of marginal model objects.
#' @param eqtlgeno_list A list of eqtl genotypes.
#' @param mean_log2fc A numeric scalar or vector for the log2 fold-change parameter to
#'     increase or decrease the conditional mean at genotype 1 \eqn{\mu_{1}} in
#'     a cell type. Default is \code{mean_log2fc = 0} (no parameters are modified
#'     and uses estimated parameters from the fitted marginal model).
#' @param eqtl_log2fc A numeric scalar or vector for the log2 fold-change parameter to
#'     increase or decrease the slope of eQTL effect in a celltype. The eQTL slope
#'     is defined as the difference between the conditional mean at genotype 1
#'     and genotype 0 (\eqn{\mu_{1}} - \eqn{\mu_{0}}). Default is
#'     \code{eqtl_log2fc = mean_log2fc} (eQTL slope is scaled the same as the
#'     conditional mean log2 fold-change).
#' @param mean_baseline A numeric scalar or vector to specify the minimum conditional mean
#'     at genotype 1 \eqn{\mu_{1}}.  If \code{mean_baseline_only = FALSE},
#'     then the conditional mean will be the maximum of the fitted (estimated from
#'     marginal model) and the \code{mean_baseline} value.  Otherwise, the
#'     conditional mean will be set to the \code{mean_baseline} value.  Default
#'     value is \code{NULL}.
#' @param eqtl_baseline A numeric scalar or vector to specify the minimum eQTL slope
#'     between genotype 1 and 0 (\eqn{\mu_{1}} - \eqn{\mu_{0}}).  If
#'     \code{eqtl_baseline_only = FALSE}, then the eQTL slope will be the
#'     maximum of the slope of fitted (estimated from marginal model) and
#'     the \code{eqtl_baseline} value.  Otherwise, the eQTL slope will be set to
#'     the \code{eqtl_baseline} value.  Default value is \code{NULL}.
#' @param mean_baseline_only A logical scalar or vector to force the conditional mean (in
#'     linear prediction) at genotype 1 \eqn{\mu_{1}}. Default is \code{FALSE}.
#' @param eqtl_baseline_only A logical scalar or vector to force the eQTL slope between
#'     genotype 1 and 0 (\eqn{\mu_{1}} - \eqn{\mu_{0}}). Default is \code{FALSE}.
#' @param features A scalar or vector of features (ie. genes) to apply the modifications.
#' @param debug A logical for whether to output a \code{mod_list} list in addition
#'     to \code{marginal_list}.
#' @inheritParams modifyModelPara
#'
#' @return A list of marginal models similar to \code{marginal_list} input.
#'     When \code{debug} = TRUE, the output have a \code{mod_list} list containing
#'     intermediate objects in addition to \code{marginal_list}.
#' @export
#'
#' @examples
#' NULL
modifyMarginalModels <- function(marginal_list,
                                 eqtlgeno_list,
                                 features,
                                 celltype,
                                 neg_ctrl = FALSE,
                                 mean_log2fc = 0,
                                 eqtl_log2fc = mean_log2fc,
                                 eqtl_reverse = FALSE,
                                 mean_baseline = NULL,
                                 eqtl_baseline = NULL,
                                 mean_baseline_only = FALSE,
                                 eqtl_baseline_only = FALSE,
                                 disp_scaling = "linear",
                                 celltype_colname = "cell_type",
                                 snp_colname = "snp_id",
                                 verbose = TRUE,
                                 debug = FALSE,
                                 ...) {

    opt_args <- list(...)

    allfeat_names <- names(marginal_list)
    n_feat <- length(features)

    # checks
    stopifnot("Please check features exist in marginal_list." =
                  (checkVectorContain(features, allfeat_names)))

    # create same length parameters as features
    neg_ctrl <- rep(neg_ctrl, length.out = n_feat)
    mean_log2fc <- rep(mean_log2fc, length.out = n_feat)
    eqtl_log2fc <- rep(eqtl_log2fc, length.out = n_feat)
    eqtl_reverse <- rep(eqtl_reverse, length.out = n_feat)
    mean_baseline <-  if(!is.null(mean_baseline)) rep(mean_baseline, length.out = n_feat)
    eqtl_baseline <- if(!is.null(eqtl_baseline)) rep(eqtl_baseline, length.out = n_feat)
    mean_baseline_only <- rep(mean_baseline_only, length.out = n_feat)
    eqtl_baseline_only <- rep(eqtl_baseline_only, length.out = n_feat)


    # iterate thru genes to modify parameters
    mod_list <- lapply(1:n_feat, FUN = function(idx) {

        if(verbose) {
            message(sprintf("Modifying parameters for %s in %s celltype...",
                            features[idx], celltype))
        }

        # for handling NULL
        mean_baseline_val <- if(is.null(mean_baseline)) NULL else mean_baseline[idx]
        eqtl_baseline_val <- if(is.null(eqtl_baseline)) NULL else eqtl_baseline[idx]

        modifyModelPara(model_obj = marginal_list[[features[idx]]][["fit"]],
                        eqtlgeno = eqtlgeno_list[[features[idx]]],
                        celltype = celltype,
                        neg_ctrl = neg_ctrl[idx],
                        mean_log2fc = mean_log2fc[idx],
                        eqtl_log2fc = eqtl_log2fc[idx],
                        eqtl_reverse = eqtl_reverse[idx],
                        mean_baseline = mean_baseline_val,
                        eqtl_baseline = eqtl_baseline_val,
                        mean_baseline_only = mean_baseline_only[idx],
                        eqtl_baseline_only = eqtl_baseline_only[idx],
                        celltype_colname = celltype_colname,
                        snp_colname = snp_colname,
                        verbose = verbose,
                        debug = debug,
                        ...)
    })
    names(mod_list) <- features


    # update models in marginal list
    for(i in features) {
        marginal_list[[i]][["fit"]] <- mod_list[[i]][["mod_new"]]
    }

    if(debug) {
        return(list("marginal_list" = marginal_list,
                    "mod_list" = mod_list))
    } else {

        mod_list_trunc <- lapply(features, FUN = function(feat) {
            mod_list[[feat]][["coef_new"]]
        })
        names(mod_list_trunc) <- features

        return(list("marginal_list" = marginal_list,
                    "coef_new" = mod_list_trunc))
    }
}


# TODO: handle vector of snps (multi-snps)
# TODO: handle marginal model w/o celltype SNP interactions

#' Modify parameters of a glmmTMB model object
#'
#' @param model_obj A marginal model object for a gene.
#' @param eqtlgeno A dataframe with eQTL annotations and samples' genotype for a gene.
#' @param celltype A string to specify the cell type in which to make the modification.
#' @param neg_ctrl A logical value for whether to set a negative control eQTL (ie. a
#'     non-eGene).  This option sets the conditional means to be identical across
#'     genotypes (0, 1, 2).  If \code{neg_ctrl = TRUE}, the \code{mean_log2fc} option
#'     will still be applied if set, but eqtl_log2fc will be overidden and have
#'     no impact.  Default is \code{FALSE}.
#' @param mean_log2fc A numeric value for the log2 fold-change parameter to
#'     increase or decrease the conditional mean at genotype 1 \eqn{\mu_{1}} in
#'     a cell type. Default is \code{mean_log2fc = 0} (no parameters are modified
#'     and uses estimated parameters from the fitted marginal model).
#' @param eqtl_log2fc A numeric value for the log2 fold-change parameter to
#'     increase or decrease the slope of eQTL effect in a celltype. The eQTL slope
#'     is defined as the difference between the conditional mean at genotype 1
#'     and genotype 0 (\eqn{\mu_{1}} - \eqn{\mu_{0}}). Default is
#'     \code{eqtl_log2fc = mean_log2fc} (eQTL slope is scaled the same as the
#'     conditional mean log2 fold-change).
#' @param eqtl_reverse A logical value to determine whether the eQTL slope trends
#'     in the reverse direction (TRUE) or same (FALSE). Default is \code{FALSE}.
#' @param mean_baseline A numeric value to specify the minimum conditional mean
#'     at genotype 1 \eqn{\mu_{1}}.  If \code{mean_baseline_only = FALSE},
#'     then the conditional mean will be the maximum of the fitted (estimated from
#'     marginal model) and the \code{mean_baseline} value.  Otherwise, the
#'     conditional mean will be set to the \code{mean_baseline} value.  Default
#'     value is \code{NULL}.
#' @param eqtl_baseline A numeric value to specify the minimum eQTL slope
#'     between genotype 1 and 0 (\eqn{\mu_{1}} - \eqn{\mu_{0}}).  If
#'     \code{eqtl_baseline_only = FALSE}, then the eQTL slope will be the
#'     maximum of the slope of fitted (estimated from marginal model) and
#'     the \code{eqtl_baseline} value.  Otherwise, the eQTL slope will be set to
#'     the \code{eqtl_baseline} value.  Default value is \code{NULL}.
#' @param mean_baseline_only A logical value to force the conditional mean (in
#'     linear prediction) at genotype 1 \eqn{\mu_{1}}. Default is \code{FALSE}.
#' @param eqtl_baseline_only A logical value to force the eQTL slope between
#'     genotype 1 and 0 (\eqn{\mu_{1}} - \eqn{\mu_{0}}). Default is \code{FALSE}.
#' @param disp_scaling A string value to specify the dispersion-mean scaling for
#'     certain parametric models. Current options are either \code{"linear"},
#'     \code{"quadratic"}, or \code{"none"}. (NOTE: currently only applicable to
#'     the negative binomial model.)
#' @param celltype_colname A string for cell type variable name.
#' @param snp_colname A string for SNP id variable name.
#' @param verbose A logical value for whether to output messages related to
#'     modified parameters. Default is \code{TRUE}.
#' @param debug A logical value for whether to output intermediate objects used
#'     for debugging purposes. Default is \code{FALSE}.
#' @param log_tol A numeric value used as tolerance in log computation. Default
#'     value is \eqn{1e-4}.
#' @param ... Additional options.
#'
#' @return A list of dataframe of coefficients, model objects, and optional outputs
#'     if debugging is enabled.
#' @export
#'
#' @examples
#' NULL
modifyModelPara <- function(model_obj,
                            eqtlgeno,
                            celltype,
                            neg_ctrl = FALSE,
                            mean_log2fc = 0,
                            eqtl_log2fc = mean_log2fc,
                            eqtl_reverse = FALSE,
                            mean_baseline = NULL,
                            eqtl_baseline = NULL,
                            mean_baseline_only = FALSE,
                            eqtl_baseline_only = FALSE,
                            disp_scaling = "linear",
                            celltype_colname = "cell_type",
                            snp_colname = "snp_id",
                            verbose = TRUE,
                            debug = FALSE,
                            log_tol = 1e-4,
                            ...) {
    # Note: currently works for single-SNP model with celltype SNP interaction effect
    #       for log link glmmTMB models.

    opt_args <- list(...)


    # extract model object, link function type, and snps
    snps <- eqtlgeno[[snp_colname]]

    mod_orig <- model_obj

    # extract model family and disp. parameter
    family_use <- stats::family(mod_orig)$family[1]
    phi_orig <- glmmTMB::sigma(mod_orig)

    if(grepl("nbinom2", family_use)) {
        family_use <- "nb"
    }


    if(methods::is(mod_orig, "glmmTMB")) {
        link_func <- mod_orig[["modelInfo"]][["family"]][["link"]]
    } else {
        stop(sprintf("Please ensure model object is glmmTMB."))
    }

    if(stringr::str_detect(snps, ":")) {
        interact_char <- ":`"
    } else {
        interact_char <- ":"
    }

    # caution messages
    if(abs(mean_log2fc) >= 5) {
        message(sprintf("Note: mean fold change is %s of the original conditional mean.",
                        2^mean_log2fc))
    }
    if(abs(eqtl_log2fc) >= 5) {
        message(sprintf("Note: eqtl fold change is %s of the original eQTL slope.",
                        2^eqtl_log2fc))
    }


    if(link_func == "log") {

        coef <- as.data.frame(summary(mod_orig)[["coefficients"]][["cond"]])

        # classify each term of model
        coef <- tibble::rownames_to_column(coef, "term") %>%
            dplyr::mutate(term = stringr::str_remove(term, celltype_colname),
                          celltype = stringr::str_detect(term, celltype),
                          snp = stringr::str_detect(term, snps),
                          interaction = stringr::str_detect(term, interact_char),
                          intercept = stringr::str_detect(term, "Intercept"))

        # check
        stopifnot("Please make sure cell type is not the baseline cell type used in contrast." =
                      (celltype %in% coef$term))

        # get parameter values and indices
        intcpt_val <- coef[coef$intercept == 1, "Estimate"]
        intcpt_idx <- which(coef[["Estimate"]] == intcpt_val)

        ct_val <- coef[coef$celltype == 1 & coef$snp == 0 & coef$interaction == 0, "Estimate"]
        ct_idx <- which(coef[["Estimate"]] == ct_val)

        snp_val <- coef[coef$celltype == 0 & coef$snp == 1 & coef$interaction == 0, "Estimate"]
        snp_idx <- which(coef[["Estimate"]] == snp_val)

        int_val <- coef[coef$celltype == 1 & coef$snp == 1 & coef$interaction == 1, "Estimate"]
        int_idx <- which(coef[["Estimate"]] == int_val)

        # initialize parameter vector for cell type
        paravec <- c(intcpt_val, ct_val, snp_val, int_val)
        names(paravec) <- c("intcpt", "ct", "snp", "int")

        paravec_new <- paravec

        # design matrix with intercept, SNP, and covariates
        design_mat <- matrix(1L, nrow = 3L, ncol = 4L,
                             dimnames = list(c("geno0", "geno1", "geno2"),
                                             c("intcpt", "ct", "snp", "int"))
                             )
        design_mat["geno0", 3:4] <- 0L
        design_mat["geno1", 3:4] <- 1L
        design_mat["geno2", 3:4] <- 2L

        # compute linear predictors at 3 genos
        predvec <- as.vector(design_mat %*% paravec)
        names(predvec) <- c("geno0", "geno1", "geno2")

        # compute nominal eqtl slope and means
        meanvec <- exp(predvec)  # conditional mean on exp scale
        eqtl <- meanvec["geno1"] - meanvec["geno0"]
        eqtl_sign <- sign(eqtl)

        if(neg_ctrl) {

            # adjust mean at geno 1
            mean1_new <- meanvec["geno1"] * 2^mean_log2fc

            if(!is.null(mean_baseline)) {

                mean_baseline <- mean_baseline * 2^mean_log2fc

                if(mean_baseline_only) {
                    mean1_new <- mean_baseline
                } else {
                    mean1_new <- ifelse(mean1_new < mean_baseline, mean_baseline, mean1_new)
                }
            }

            # solve for celltype main effect
            paravec_new["ct"] <- log(mean1_new) - paravec["intcpt"]

            # solve for celltype snp interaction effect
            paravec_new["int"] <- -paravec["snp"]

        } else {

            # adjust mean at geno 1 and eqtl slope
            mean1_new <- meanvec["geno1"] * 2^mean_log2fc

            # check if use mean baseline
            if(!is.null(mean_baseline)) {

                mean_baseline <- mean_baseline * 2^mean_log2fc

                if(mean_baseline_only) {
                    mean1_new <- mean_baseline
                } else {
                    mean1_new <- ifelse(mean1_new < mean_baseline, mean_baseline, mean1_new)
                }
            }

            # check if reverse eQTL direction
            eqtl_new <- ifelse(eqtl_reverse, -eqtl, eqtl) * 2^eqtl_log2fc

            # check if use eqtl baseline
            if(!is.null(eqtl_baseline)) {

                eqtl_baseline <- ifelse(eqtl_reverse, -eqtl_baseline, eqtl_baseline) *
                    2^eqtl_log2fc

                if(eqtl_baseline_only) {
                    eqtl_new <- ifelse(eqtl_sign >= 0, 1, -1) * eqtl_baseline
                } else {
                    eqtl_new <- ifelse(eqtl_new < eqtl_baseline,
                                       ifelse(eqtl_sign >= 0, 1, -1) * eqtl_baseline,
                                       eqtl_new)
                }
            }

            # solve for celltype main effect
            mean1eqtl_new_diff <- ifelse(mean1_new - eqtl_new <= 0,
                                         log_tol,
                                         mean1_new - eqtl_new)

            paravec_new["ct"] <- log(mean1eqtl_new_diff) - paravec["intcpt"]

            # solve for celltype snp interaction effect
            paravec_new["int"] <- log(mean1_new) - paravec["intcpt"] -
                                    paravec_new["ct"] - paravec["snp"]

        }

        # update adjusted eqtl slope and means
        predvec_new <- as.vector(design_mat %*% paravec_new)
        names(predvec_new) <- c("geno0", "geno1", "geno2")

        meanvec_new <- exp(predvec_new)

        # update new parameters
        coef_new <- coef
        coef_new[ct_idx, "Estimate"] <- paravec_new["ct"]
        coef_new[int_idx, "Estimate"] <- paravec_new["int"]

        # update parameters in model object
        mod_new <- mod_orig

        mod_new[["fit"]][["par"]][ct_idx] <- paravec_new["ct"]
        mod_new[["fit"]][["parfull"]][ct_idx] <- paravec_new["ct"]

        mod_new[["fit"]][["par"]][int_idx] <- paravec_new["int"]
        mod_new[["fit"]][["parfull"]][int_idx] <- paravec_new["int"]

        # update disp. parameter in model
        phi_new <- phi_orig

        if(!is.null(phi_orig) && disp_scaling == "linear") {

            if(family_use == "nb") {

                # linear scale to mu^2 / phi
                phi_new <- meanvec_new[["geno1"]] * phi_orig / meanvec[["geno1"]]

                mod_new[["fit"]][["par"]][["betad"]] <- log(phi_new)
                mod_new[["fit"]][["parfull"]]["betad"] <- log(phi_new)
            }
        } else if(!is.null(phi_orig) && disp_scaling == "quadratic") {

            if(family_use == "nb") {

                # linear scale to mu^2 / phi
                phi_new <- meanvec_new[["geno1"]]^2 * phi_orig / meanvec[["geno1"]]^2

                mod_new[["fit"]][["par"]][["betad"]] <- log(phi_new)
                mod_new[["fit"]][["parfull"]]["betad"] <- log(phi_new)
            }
        } else if(!is.null(phi_orig) && disp_scaling == "none") {
            NULL
        }

        # output messages
        if(verbose) {
            message(sprintf("celltype effect: %.5f ===> new value: %.5f",
                            paravec["ct"], paravec_new["ct"]))
            message(sprintf("interaction effect: %.5f ===> new value: %.5f",
                            paravec["int"], paravec_new["int"]))
            if(neg_ctrl) {
                message(sprintf("<< Setting conditional means to be negative controls (no eQTL effect) >>"))
            }
            message(sprintf("conditional means at geno 0, 1, 2: %.5f, %.5f, %.5f ===> new values: %.5f, %.5f, %.5f",
                            meanvec["geno0"], meanvec["geno1"], meanvec["geno2"],
                            meanvec_new["geno0"], meanvec_new["geno1"], meanvec_new["geno2"]))
            if(eqtl_reverse) {
                message(sprintf("<< Setting the eQTL slope in reverse direction >>"))
            }
            message(sprintf("eQTL slope between geno 1 and 0: %.5f ===> new value: %.5f",
                            meanvec["geno1"] - meanvec["geno0"],
                            meanvec_new["geno1"] - meanvec_new["geno0"]))
            if(!is.null(phi_new)) {
                message(sprintf("Phi parameter: %.5f ===> new value: %.5f\n",
                                phi_orig, phi_new))
            }
        }
    } else {
        stop(sprintf("Please check link function is log."))
    }

    if(debug) {
        return(list("coef_new" = coef_new,
                    "mod_new" = mod_new,
                    ## diagnosis outputs
                    "coef" = coef,
                    "mod_orig" = mod_orig,
                    "means_orig" = meanvec,
                    "means_new" = meanvec_new,
                    "snps" = snps,
                    "paravec" = paravec,
                    "paravec_new" = paravec_new,
                    "eqtlslope" = meanvec["geno1"] - meanvec["geno0"],
                    "eqtlslope_new" = meanvec_new["geno1"] - meanvec_new["geno0"]
                    ))
    } else {
        return(list("coef_new" = coef_new,
                    "mod_new" = mod_new
                    ))
    }
}
