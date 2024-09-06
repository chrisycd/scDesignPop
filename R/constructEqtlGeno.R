#' Construct eQTL annotation and genotype dataframe
#'
#' This function creates a combined dataframe using input eQTL annotations, and sample genotype dataframe.
#'
#' @param eqtl_annot_df eQTL annotation dataframe
#' @param geno_df eQTL by sample genotype dataframe
#' @param sampid_vec sample id vector
#' @param name string or integer scalar to identify input
#' @param feature_colname feature variable name
#' @param cellstate_colname cell state variable name
#' @param snp_colname SNP id variable name
#' @param loc_colname SNP position variable name
#'
#' @return a dataframe with eQTL annotations and genotype of SNPs
#' @export
#'
#' @examples
#' NULL
constructEqtlGeno <- function(eqtl_annot_df,
                              geno_df,
                              sampid_vec,
                              name = NULL,
                              feature_colname = "gene_name",
                              cellstate_colname = "cell_type",
                              snp_colname = "snp_id",
                              loc_colname = "POS") {
    # TODO: add checks and filtering for sample id in geno_df
    # TODO: add assertthat checks for CHR columns

    assertthat::assert_that(assertthat::has_name(eqtl_annot_df,
                                                 c(feature_colname, snp_colname,
                                                   loc_colname, cellstate_colname)))
    assertthat::assert_that(assertthat::has_name(geno_df, snp_colname))

    if(is.null(name)) {
        name <- ""
    } else {
        name <- sprintf(" on %s", name)
    }

    # convert to additive dosage 0, 1, or 2 (number of alternative alleles)
    geno_dat <- geno_df %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::all_of(sampid_vec),
                                    .fns = ~ dplyr::case_when(. == "0|0" | . == "0/0" ~ 0L,
                                                              . == "0|1" | . == "1|0" | . == "0/1" | . == "1/0" ~ 1L,
                                                              . == "1|1" | . == "1/1" ~ 2L,
                                                              . == ".|." | . == "./." ~ NA_integer_,
                                                              TRUE ~ NaN)))


    # check number of samples with missing genotypes
    na_count <- sapply(geno_dat, function(x) sum(length(which(is.na(x)))))
    if(any(na_count[-1] > 0)) {
        message(sprintf("Sample with missing genotypes%s! Consider filtering before proceeding.", name))
    }

    # merge eQTL annotations with SNP genotype info
    eqtl_geno <- eqtl_annot_df %>%
        dplyr::relocate(!!rlang::sym(loc_colname), .after = tidyselect::last_col()) %>%  # move POS to last column (used later)
        dplyr::left_join(., geno_dat, by = snp_colname, multiple = "first") %>%
        dplyr::arrange(!!rlang::sym(feature_colname)) %>%
        dplyr::ungroup(.)
    # Note: duplicate eSNPs IDs in geno_dat will use first one with multiple = "first" option

    # filter out rows with NA genotype and output message
    firstcol_indx <- which(colnames(eqtl_geno) == loc_colname) + 1L
    lastcol_indx <- dim(eqtl_geno)[2]

    initial_rows <- dim(eqtl_geno)[1]

    eqtl_geno <- eqtl_geno %>%
        dplyr::filter(base::rowSums(dplyr::across(colnames(.)[firstcol_indx]:colnames(.)[lastcol_indx], ~ !is.na(.))) > 0)

    final_rows <- dim(eqtl_geno)[1]

    filter_row_count <- initial_rows - final_rows

    if(filter_row_count > 0) {
        message(sprintf("%s unmatched rows in eqtl_geno filtered out%s.", filter_row_count, name))
    }

    return(eqtl_geno)
}
