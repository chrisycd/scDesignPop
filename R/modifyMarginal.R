#' Modify marginal models
#'
#' Iterate through genes and modify marginal model parameters based on user inputs
#'     for celltype-specific eQTLs.
#'
#' @param marginal_list A list of marginal model objects.
#' @param eqtlgeno_list A list of eqtl genotypes.
#' @param eqtl_loc A non-negative numeric scalar or vector  used as multiplicative factor to
#'     increase or decrease the celltype eQTL mean effect size at genotype 1. The
#'     default is eqtl_loc = NULL, where no effect size changes are made.
#' @param eqtl_scale A non-negative numeric scalar or vector used as multiplicative factor
#'     to adjust the slope of celltype eQTL effect size. If NULL then no effect
#'     size changes are made.  The default is same as \code{eqtl_loc}.
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
                                 eqtl_loc = NULL, # eqtl mean fold change
                                 eqtl_scale = eqtl_loc,  # eqtl diff fold change
                                 cellstate_colname = "cell_type",
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
    if(!is.null(eqtl_loc)) {
        eqtl_loc <- rep(eqtl_loc, length.out = n_feat)
    }

    if(!is.null(eqtl_scale)) {
        eqtl_scale <- rep(eqtl_scale, length.out = n_feat)
    }

    # iterate thru genes to modify parameters
    mod_list <- lapply(1:n_feat, FUN = function(idx) {

        if(verbose) {
            message(sprintf("Modifying parameters for %s in %s celltype...",
                            features[idx], celltype))
        }

        modifyModelPara(model_obj = marginal_list[[features[idx]]][["fit"]],
                        eqtlgeno = eqtlgeno_list[[features[idx]]],
                        celltype = celltype,
                        neg_ctrl = neg_ctrl,
                        eqtl_loc = eqtl_loc[idx], # eqtl mean fold change
                        eqtl_scale = eqtl_scale[idx],  # eqtl diff fold change
                        direction = NULL,  # not implemented
                        cellstate_colname = cellstate_colname,
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
        return(marginal_list)
    }
}


# TODO: use log fold change for eqtl_mean and eqtl_slope options
# TODO: handle vector of snps (multi-snps)
# TODO: handle marginal model w/o celltype SNP interactions

#' Modify parameters of a glmmTMB model object
#'
#' @param model_obj A marginal model object for a gene.
#' @param eqtlgeno A dataframe of eqtl genotypes for a gene.
#' @param celltype A string to specify the cell type.
#' @param neg_ctrl A logical value for whether to set negative control eQTL.
#'     If neg_ctrl = TRUE, eqtl_loc option can still be used, while eqtl_scale
#'     will have no impact. This option sets the conditional means to be the
#'     same across genotypes (0, 1, 2).
#' @param eqtl_loc A non-negative numeric value used as multiplicative factor to
#'     increase or decrease the celltype eQTL mean effect size at genotype 1. The
#'     default is eqtl_loc = NULL, where no effect size changes are made.
#' @param eqtl_scale A non-negative numeric value used as multiplicative factor
#'     to adjust the slope of celltype eQTL effect size. If NULL then no effect
#'     size changes are made.
#' @param direction A logical value to determine whether the eQTL slope is in
#'     nominal or reverse direction.
#' @param cellstate_colname A string for cell state variable name.
#' @param snp_colname A string for SNP id variable name.
#' @param verbose A logical value for whether to output messages related to
#'     modified parameters. The default is TRUE.
#' @param debug A logical value for whether to output intermediate objects used
#'     for debugging purposes. The default is FALSE.
#' @param ... Additional options.
#'
#' @return A list of dataframe of coefficients, model objects, and optional .
#' @export
#'
#' @examples
#' NULL
modifyModelPara <- function(model_obj,
                            eqtlgeno,
                            celltype,
                            neg_ctrl = FALSE,
                            eqtl_loc = NULL, # eqtl mean shift relative to genotype 1
                            eqtl_scale = eqtl_loc,  # eqtl mean slope change relative to genotype 0
                            direction = NULL,  # not implemented
                            cellstate_colname = "cell_type",
                            snp_colname = "snp_id",
                            verbose = TRUE,
                            debug = FALSE,
                            ...) {
    # Note: currently works for single-SNP model with celltype SNP interaction effect
    #       for log link glmmTMB models.

    opt_args <- list(...)


    # extract model object, link function type, and snps
    snps <- eqtlgeno[[snp_colname]]

    mod_orig <- model_obj

    if(methods::is(mod_orig, "glmmTMB")) {
        link <- mod_orig[["modelInfo"]][["family"]][["link"]]
    } else {
        stop(sprintf("Please ensure model object is glmmTMB."))
    }

    if(stringr::str_detect(snps, ":")) {
        interact_char <- ":`"
    } else {
        interact_char <- ":"
    }


    if(link == "log") {

        coef <- as.data.frame(summary(mod_orig)[["coefficients"]][["cond"]])

        # classify each term of model
        coef <- tibble::rownames_to_column(coef, "term") %>%
            dplyr::mutate(term = stringr::str_remove(term, cellstate_colname),
                          celltype = stringr::str_detect(term, celltype),
                          snp = stringr::str_detect(term, snps),
                          interaction = stringr::str_detect(term, interact_char),
                          intercept = stringr::str_detect(term, "Intercept"))

        # check
        stopifnot("Please make sure cell type is not the baseline cell type used in contrast." =
                      (celltype %in% coef$term))

        # get parameter values and indices
        intcpt_val <- coef[coef$intercept == 1, ][["Estimate"]]
        intcpt_idx <- which(coef[["Estimate"]] == intcpt_val)

        ct_val <- coef[coef$celltype == 1 & coef$snp == 0 & coef$interaction == 0, ][["Estimate"]]
        ct_idx <- which(coef[["Estimate"]] == ct_val)

        snp_val <- coef[coef$celltype == 0 & coef$snp == 1 & coef$interaction == 0, ][["Estimate"]]
        snp_idx <- which(coef[["Estimate"]] == snp_val)

        int_val <- coef[coef$celltype == 1 & coef$snp == 1 & coef$interaction == 1, ][["Estimate"]]
        int_idx <- which(coef[["Estimate"]] == int_val)


        # compute nominal mean at genotypes
        mean_0 <- exp(intcpt_val + ct_val + snp_val * 0 + int_val * 0)
        mean_1 <- exp(intcpt_val + ct_val + snp_val * 1 + int_val * 1)
        mean_2 <- exp(intcpt_val + ct_val + snp_val * 2 + int_val * 2)
        # formula: exp(intercept + celltype main eff + snp main eff + celltype snp interaction effect)


        if(neg_ctrl) {

            # adjust mean and solve for new params
            mean_1_new <- mean_1 * ifelse(is.null(eqtl_loc), 1, eqtl_loc)
            mean_0_new <- mean_1_new

            # solve for celltype main eff
            ct_new <- log(mean_0_new) - intcpt_val

            # solve for celltype snp interaction effect
            int_new <- -snp_val

        } else {

            # adjust mean and solve for new params
            mean_1_new <- mean_1 * ifelse(is.null(eqtl_loc), 1, eqtl_loc)
            mean_0_new <- mean_0 * ifelse(is.null(eqtl_scale), 1, eqtl_scale)

            # solve for celltype main eff
            ct_new <- log(mean_0_new) - intcpt_val

            # solve for celltype snp interaction effect
            int_new <- log(mean_1_new) - intcpt_val - ct_new - snp_val * 1
        }

        # update new parameters
        coef_new <- coef
        coef_new[ct_idx, ][["Estimate"]] <- ct_new
        coef_new[int_idx, ][["Estimate"]] <- int_new

        # compute updated mean at genotypes
        mean_0_calc <- coef_new[intcpt_idx, ][["Estimate"]] +
            coef_new[ct_idx, ][["Estimate"]] +       # celltype main eff
            coef_new[snp_idx, ][["Estimate"]] * 0 +  # snp main effect
            coef_new[int_idx, ][["Estimate"]] * 0    # celltype snp interaction effect
        mean_0_calc <- exp(mean_0_calc)

        mean_1_calc <- coef_new[intcpt_idx, ][["Estimate"]] +
            coef_new[ct_idx, ][["Estimate"]] +
            coef_new[snp_idx, ][["Estimate"]] * 1 +
            coef_new[int_idx, ][["Estimate"]] * 1

        mean_1_calc <- exp(mean_1_calc)

        mean_2_calc <- coef_new[intcpt_idx, ][["Estimate"]] +
            coef_new[ct_idx, ][["Estimate"]] +
            coef_new[snp_idx, ][["Estimate"]] * 2 +
            coef_new[int_idx, ][["Estimate"]] * 2

        mean_2_calc <- exp(mean_2_calc)


        # update parameters in model object
        mod_new <- mod_orig

        mod_new[["fit"]][["par"]][ct_idx] <- ct_new
        mod_new[["fit"]][["parfull"]][ct_idx] <- ct_new

        mod_new[["fit"]][["par"]][int_idx] <- int_new
        mod_new[["fit"]][["parfull"]][int_new] <- int_new

        if(verbose) {
            message(sprintf("celltype effect: %.4f ===> new value: %.4f",
                            ct_val, ct_new))
            message(sprintf("interaction effect: %.4f ===> new value: %.4f",
                            int_val, int_new))
            if(neg_ctrl) {
                message(sprintf("<< Setting conditional means to be negative controls (no eQTL effect) >>"))
            }
            message(sprintf("conditional means at 0,1,2 geno: %.4f, %.4f, %.4f ===> new values: %.4f, %.4f, %.4f\n",
                            mean_0, mean_1, mean_2, mean_0_calc, mean_1_calc, mean_2_calc))
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
                    "means_orig" = c(mean_0, mean_1, mean_2),
                    "snps" = snps,
                    "ct_new" = ct_new,
                    "int_new" = int_new,
                    "mean_1_new" = mean_1_new,
                    "mean_0_new" = mean_0_new,
                    "means_calc" = c(mean_0_calc, mean_1_calc, mean_2_calc)
        ))
    } else {
        return(list("coef_new" = coef_new,
                    "mod_new" = mod_new
        ))
    }
}
