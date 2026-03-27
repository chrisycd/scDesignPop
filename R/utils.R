################# Helper functions #################

#' Check if multiple vectors have same elements
#'
#' This is an internal helper function to check if two or more vectors have same
#'     elements (with their order considered as an option).
#'
#' @param ... arguments of two or more vectors.
#' @param ignore_order logical scalar to disregard the ordering of elements in
#'     input vectors. Default is TRUE.
#'
#' @return a logical scalar
#' @export
#'
#' @examples
#' vec1 <- c("cherry", "apple", "banana", "cherry")
#' vec2 <- c("cherry", "apple", "cherry", "banana")
#' vec3 <- c("banana", "cherry", "apple", "cherry")
#'
#' checkVectorEqual(vec1, vec2, vec3, ignore_order = TRUE)  # returns TRUE
#' checkVectorEqual(vec1, vec2, vec3, ignore_order = FALSE)  # returns FALSE
checkVectorEqual <- function(..., ignore_order = TRUE) {

    vec_list <- list(...)

    if(length(vec_list) < 2) {
        stop(sprintf("Input must be 2 or more vectors!"))
    }

    first_vec <- vec_list[[1]]

    # compare all vectors
    if(ignore_order) {
        res <- all(base::sapply(vec_list[-1], function(x) base::setequal(first_vec, x)))
    } else {
        res <- all(base::sapply(vec_list[-1], function(x) base::identical(first_vec, x)))
    }

    return(res)
}


#' Check membership of first vector compared to other vectors
#'
#' This is an internal helper function to check if elements of the first vector are
#'     subsets of one or more other vectors. Ordering of elements is ignored.
#'
#' @param ... arguments of two or more vectors.
#' @param ignore_dups logical scalar to disregard duplicate elements during comparison.
#'     Default is FALSE.
#'
#' @return a logical scalar
#' @export
#'
#' @examples
#' vec1 <- c("cherry", "apple", "cherry")
#' vec2 <- c("apple", "cherry", "banana", "cherry")
#' vec3 <- c("banana", "kiwi", "apple", "cherry")
#' vec4 <- c("banana", "apple", "apple", "cherry")
#'
#' checkVectorContain(vec1, vec2, vec3, ignore_dups = TRUE)   # returns TRUE
#' checkVectorContain(vec1, vec2, vec4, ignore_dups = FALSE)  # returns FALSE
checkVectorContain <- function(..., ignore_dups = FALSE) {

    vec_list <- list(...)

    if(length(vec_list) < 2) {
        stop(sprintf("Input must be 2 or more vectors!"))
    }

    first_vec <- vec_list[[1]]

    # compare membership of first vector to all others
    if(ignore_dups) {
        res <- all(base::sapply(vec_list[-1], function(x) base::is.element(first_vec, x)))
    } else {
        res <- all(base::sapply(vec_list[-1], function(x) length(base::setdiff(first_vec, x)) == 0))
    }

    return(res)
}


#' Check individuals in eQTL-genotype and covariate input
#'
#' This is an internal helper function to check that all data frames in \code{eqtlgeno_list}
#' contain identical individual ID columns, and that these individuals match those
#' in \code{covariate} data frame.
#'
#' @param eqtlgeno_list A list of eQTL genotype data frames.
#' @param covariate_df A data frame of covariates.
#' @param eqtlgeno_name A string scalar specifying name of \code{eqtlgeno_list}.
#' @param covariate_name A string scalar specifying name of \code{covariate_df}.
#' @param loc_colname A string scalar specifying the column name immediately before
#'     the individuals' genotype columns. Default is \code{"POS"}.
#' @param indiv_colname A string scalar specifying the individual ID variable in
#'     \code{covariate}. Default is \code{"indiv"}.
#'
#' @return Invisibly returns \code{NULL}. Stops if validation fails.
#'
#' @keywords internal
#' @noRd
checkIndividuals <- function(eqtlgeno_list,
                             covariate_df,
                             eqtlgeno_name = "eqtlgeno_list",
                             covariate_name = "covariate",
                             loc_colname = "POS",
                             indiv_colname = "indiv") {

    indivs_list <- lapply(eqtlgeno_list, function(df) {
        firstgeno_idx <- match(loc_colname, colnames(df)) + 1L
        eqtl_indivs <- colnames(df)[firstgeno_idx:ncol(df)]

        return(eqtl_indivs)
    })

    indiv_ref <- indivs_list[[1]]
    all_same <- all(vapply(indivs_list,
                           function(x) { identical(x, indiv_ref) },
                           logical(1))
    )

    if(!all_same) {
        stop(sprintf("Not all data frames have identical individual ID columns in %s. Please check input!",
                     eqtlgeno_name))
    }

    if(!checkVectorEqual(indiv_ref,
                         unique(covariate_df[[indiv_colname]]),
                         ignore_order = TRUE)) {
        stop(sprintf("The individual IDs in %s and %s do not match. Please check input!",
                     eqtlgeno_name, covariate_name))
    }

    invisible(NULL)
}


#' Check eQTL-geno data frame for duplicate SNPs but different genotypes for a feature
#'
#' This is an internal helper function to check an eQTL genotype data frame input.
#' It checks whether the same feature-SNP pair appears in multiple rows with different
#' genotype values across the genotype columns. Such duplicated keys with
#' inconsistent genotypes are not allowed.
#'
#' @param df A data frame containing eQTL annotations and genotype columns.
#' @param df_name A string scalar giving the dataset label to include in
#'     error messages.
#' @param feature_colname A string scalar specifying the feature column name.
#'     Default is \code{"gene_id"}.
#' @param snp_colname A string scalar specifying the SNP column name. Default is
#'     \code{"snp_id"}.
#' @param firstgeno_idx An integer scalar giving the first genotype column index
#'     in \code{df}.
#' @param lastgeno_idx An integer scalar giving the last genotype column index
#'     in \code{df}.
#'
#' @return Invisibly returns \code{TRUE} if no duplicated feature-SNP pairs with
#'     inconsistent genotypes are found. Otherwise, the function stops with an
#'     error.
#'
#' @details
#' The function first identifies duplicated feature-SNP keys using
#' \code{feature_colname} and \code{snp_colname}. For each duplicated key, it
#' then compares the genotype values across columns \code{firstgeno_idx} through
#' \code{lastgeno_idx}. If any duplicated key has non-identical genotype rows,
#' the input is considered invalid.
#'
#' Empty or \code{NULL} inputs are treated as valid and return invisibly.
#'
#' @keywords internal
#' @noRd
checkNonduplicateEqtlGeno <- function(df,
                                      df_name,
                                      feature_colname = "gene_id",
                                      snp_colname = "snp_id",
                                      firstgeno_idx,
                                      lastgeno_idx) {

    if(is.null(df) || nrow(df) == 0L) return(invisible(TRUE))

    dup_keys <- df %>%
        dplyr::count(.data[[feature_colname]],
                     .data[[snp_colname]]) %>%
        dplyr::filter(.data$n > 1L)

    if(nrow(dup_keys) == 0L) return(invisible(TRUE))

    bad <- dup_keys %>%
        dplyr::rowwise() %>%
        dplyr::filter({
            sub_df <- df[df[[feature_colname]] == .data[[feature_colname]] &
                             df[[snp_colname]] == .data[[snp_colname]], ]

            geno_mat <- as.matrix(sub_df[, firstgeno_idx:lastgeno_idx, drop = FALSE])
            length(unique(apply(geno_mat, 1, paste, collapse = "|"))) > 1L}) %>%
        dplyr::ungroup()

    if(nrow(bad) > 0L) {
        stop(sprintf(
            "Different genotypes found for the same SNP ID (e.g. %s, %s). Please check %s input!",
            bad[[feature_colname]][1], bad[[snp_colname]][1], df_name)
            )
    }

    invisible(TRUE)
}


# Fix global variable notes
utils::globalVariables(c("."))
utils::globalVariables(c(":="))
