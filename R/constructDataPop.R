
#' Construct a list of input data
#'
#' Function extracts an expression matrix, cell covariates, and filters the SNPs
#'     in the eQTL genotype dataframe.
#'
#' @param sce a SingleCellExperiment object
#' @param eqtlgeno_df a dataframe with eQTL annotations and SNP genotypes for each gene
#' @param new_eqtlgeno_df a dataframe with eQTL annotations and SNP genotypes for each
#'     gene in new individuals
#' @param new_covariate a cell covariate dataframe for which to simulate
#' @param overlap_features an optional string vector to filter for features
#' @param sampid_vec an optional string vector to filter for sample ids
#' @param slot_name a string scalar to specify the input SCE slot used
#' @param snp_model a string scalar to specify the type of SNP model used
#' @param cellstate_colname a string scalar for the cell state variable in eqtlgeno_df
#'     and cell covariate of SCE
#' @param feature_colname a string scalar for the feature variable in eqtlgeno_df
#' @param snp_colname a string scalar for the SNP variable in eqtlgeno_df
#' @param loc_colname a string scalar for the last column of eQTL annotation in eqtlgeno_df
#' @param chrom_colname a string scalar for the chromosome variable in eqtlgeno_df
#' @param indiv_colname a string scalar for the sample ID variable in cell covariate of SCE
#' @param prune_thres numerical value between 0 and 1 to threshold the SNPs
#' @param ct_copula a logical scalar for whether to fit the copula by cell type
#'
#' @return outputs a list with following elements:
#' \describe{
#'      \item{\code{count_mat}}{a cell-by-gene matrix of response values.}
#'      \item{\code{covariate}}{a cell-by-covariate data frame used for fit marginal.}
#'      \item{\code{new_covariate}}{an optional cell-by-covariate data frame used for prediction.}
#'      \item{\code{important_features}}{a string vector of gene ids.}
#'      \item{\code{eqtl_geno_list}}{a list of eQTL genotype dataframes for each gene
#'          (for fit marginal).}
#'      \item{\code{new_eqtl_geno_list}}{a optional list of eQTL genotype dataframes
#'          for each gene (for prediction on new individuals).}
#'      \item{\code{filtered_gene}}{string vector of features QC filtered.}
#' }
#'
#' @export
#'
#' @examples
#' NULL
constructDataPop <- function(sce,
                             eqtlgeno_df,
                             new_eqtlgeno_df = NULL,
                             new_covariate = NULL,
                             overlap_features = NULL,
                             sampid_vec = NULL,
                             ct_copula = TRUE,
                             slot_name = "counts",
                             snp_model = c("single", "multi"),
                             cellstate_colname = "cell_type",
                             feature_colname = "gene_id",
                             snp_colname = "snp_id",
                             loc_colname = "POS",
                             chrom_colname = "CHR",
                             indiv_colname = "indiv",
                             prune_thres = 0.9) {
    # TODO: add new_eqtlgeno_df : an optional list of eQTL genotype dataframes for each gene (for prediction)

    snp_model <- match.arg(snp_model)

    if(!is.null(new_eqtlgeno_df)) {
        has_newindiv <- TRUE
    } else {
        has_newindiv <- FALSE
    }

    # checks
    assertthat::assert_that(assertthat::has_name(eqtlgeno_df,
                                                 c(cellstate_colname, feature_colname,
                                                   snp_colname, chrom_colname, loc_colname)))

    if(has_newindiv) {
        assertthat::assert_that(assertthat::has_name(new_eqtlgeno_df,
                                                     c(cellstate_colname, feature_colname,
                                                       snp_colname, chrom_colname, loc_colname)))
    }

    sce_features <- rownames(sce)

    firstcol_indx <- which(colnames(eqtlgeno_df) == loc_colname) + 1L  # used as last column b4 genotype
    lastcol_indx <- dim(eqtlgeno_df)[2]
    eqtl_indivs <- colnames(eqtlgeno_df)[firstcol_indx:lastcol_indx]

    ## check unique cell names and feature names
    stopifnot("Please make sure there are no duplicate cell names in SingleCellExperiment object input." =
                  (length(colnames(sce)) == length(unique(colnames(sce)))))

    stopifnot("Please make sure there are no duplicate feature names in SingleCellExperiment object input." =
                  (length(sce_features) == length(unique(sce_features))))

    stopifnot("Please make sure there are no duplicate sample IDs in eQTL geno input." =
                  (length(eqtl_indivs) == length(unique(eqtl_indivs))))

    if(has_newindiv) {
        new_firstcol_indx <- which(colnames(new_eqtlgeno_df) == loc_colname) + 1L  # used as last column b4 genotype
        new_lastcol_indx <- dim(new_eqtlgeno_df)[2]
        new_eqtl_indivs <- colnames(new_eqtlgeno_df)[new_firstcol_indx:new_lastcol_indx]

        stopifnot("Please make sure there are no duplicate sample IDs in new eQTL geno input." =
                      (length(new_eqtl_indivs) == length(unique(new_eqtl_indivs))))
    }

    ## filter and check features
    if(!is.null(overlap_features)) {

        stopifnot("Please make sure there are no duplicate feature names in overlap_features." =
                      (length(overlap_features) == length(unique(overlap_features))))

        stopifnot("Please make sure all feature names in overlap_features exist in SingleCellExperiment object input." =
                      checkVectorContain(overlap_features, sce_features))

        # filter sce inputs
        sce <- sce[overlap_features, ]

        # filter eQTL geno input
        eqtlgeno_df <- eqtlgeno_df %>%
            dplyr::filter(!!rlang::sym(feature_colname) %in% overlap_features)

        if(has_newindiv) {
            new_eqtlgeno_df <- new_eqtlgeno_df %>%
                dplyr::filter(!!rlang::sym(feature_colname) %in% overlap_features)
        }
    }

    eqtl_features <- unique(eqtlgeno_df[[feature_colname]])
    sce_features <- rownames(sce)  # update

    if(has_newindiv) {
        new_eqtl_features <- unique(new_eqtlgeno_df[[feature_colname]])

        stopifnot("The features in eqtlgeno_df, new_eqtlgeno_df, sce do not match. Please check input!" =
                      checkVectorEqual(eqtl_features, new_eqtl_features, sce_features,
                                       ignore_order = TRUE))
    } else {
        new_eqtl_features <- NULL

        stopifnot("The features in eqtlgeno_df and sce do not match. Please check input!" =
                      checkVectorEqual(eqtl_features, sce_features, ignore_order = TRUE))
    }


    ## check and filter indiv ids
    if(!is.null(sampid_vec)) {  # filter for sample ids

        stopifnot("Please make sure there are no duplicate sample IDs in sampid_vec." =
                      (length(sampid_vec) == length(unique(sampid_vec))))

        stopifnot("Please make sure all sample IDs in sampid_vec exist in eqtlgeno_df input." =
                      checkVectorContain(sampid_vec, eqtl_indivs))

        # filter eQTL genotype
        eqtlgeno_df <- eqtlgeno_df[, c(colnames(eqtlgeno_df)[1:(firstcol_indx - 1)], sampid_vec)]

        # filter sce inputs
        sce <- which((SingleCellExperiment::colData(sce)[[indiv_colname]] %>%   # TODO: try as.character(sce[[indiv_colname]]) %in% sampid_vec
                          as.character(.)) %in% sampid_vec) %>%
            sce[, .]

    }

    # update
    firstcol_indx <- which(colnames(eqtlgeno_df) == loc_colname) + 1L  # used as last column b4 genotype
    lastcol_indx <- dim(eqtlgeno_df)[2]
    eqtl_indivs <- colnames(eqtlgeno_df)[firstcol_indx:lastcol_indx]

    sce_indivs <- extractFromSCE(sce,
                                 slot_name = slot_name,
                                 indiv_colname = indiv_colname,
                                 sc_indiv = TRUE)[["sc_indiv"]]


    stopifnot("The distinct sample IDs in SingleCellExperiment input does not match columns of eQTL geno input!" =
                  checkVectorEqual(sce_indivs, eqtl_indivs, ignore_order = TRUE))


    ## filter SNP genotypes for each gene in eqtlgeno_df
    important_features <- sce_features

    message("Constructing eqtlgeno list...")
    eqtl_geno_list <- pbapply::pblapply(important_features, function(g) {  # pbapply shows progress bar

        eqtl_tmp <- which(eqtlgeno_df[[feature_colname]] == g) %>%
            eqtlgeno_df[., ]

        geno_tmp <- eqtl_tmp[, firstcol_indx:lastcol_indx]

        if(snp_model == "single") {  # single SNP model

            # choose highest variance across indiv
            var_tmp <- geno_tmp %>%
                base::apply(., MARGIN = 1, FUN = stats::var)

            eqtl_tmp <- eqtl_tmp %>%
                dplyr::mutate("var" = var_tmp) %>%
                dplyr::relocate("var", .before = !!rlang::sym(loc_colname))

            max_var_indx <- base::which.max(var_tmp)

            return(eqtl_tmp[max_var_indx, ])
            # assign(g, eqtl_tmp[max_var_indx, ])
            # TODO: add feature name to output list

        } else if(snp_model == "multi") {  # multi-SNP model

            # prune SNPs then extract dataframe per gene
            snp_indx <- geno_tmp %>%
                as.matrix(.) %>%
                pruneSnp(., ld_threshold = prune_thres)

            eqtl_tmp <- eqtl_tmp %>%
                dplyr::slice(snp_indx)

            return(eqtl_tmp)

        } else {
            stop("Option for snp_model must be either 'single' or 'multi'!")
        }

    })
    names(eqtl_geno_list) <- important_features

    ## filter SNP genotypes for each gene in new_eqtlgeno_df
    if(has_newindiv) {
        message("New indiv eqtlgeno input detected. Constructing new_eqtlgeno_list...")

        new_eqtl_geno_list <- pbapply::pblapply(important_features, function(g) {

            ref_eqtl_tmp <- eqtl_geno_list[[g]]

            eqtl_tmp <- new_eqtlgeno_df[which(new_eqtlgeno_df[[feature_colname]] == g), ] %>%
                dplyr::filter(!!rlang::sym(cellstate_colname) %in% ref_eqtl_tmp[[cellstate_colname]],
                              !!rlang::sym(snp_colname) %in% ref_eqtl_tmp[[snp_colname]])

            return(eqtl_tmp)
        })
        names(new_eqtl_geno_list) <- important_features

    } else {
        new_eqtl_geno_list <- NULL
    }


    # create data_list
    count_mat <- extractFromSCE(sce,
                                slot_name = slot_name,
                                indiv_colname = indiv_colname,
                                expr_mat = TRUE)[["expr_mat"]]
    count_mat <- t(as.matrix(count_mat)) # %>%
    # Matrix::Matrix(., sparse = TRUE)  # TODO: fix later


    covariate_df <- extractFromSCE(sce,
                                   slot_name = slot_name,
                                   indiv_colname = indiv_colname,
                                   sc_cov  = TRUE)[["sc_cov"]]


    ## QC features based on expression
    qc_features <- apply(count_mat, MARGIN = 2, function(x){
        return(length(which(x < 1e-5 & x > -1e-5)) > length(x) - 2)
    })

    if(length(which(qc_features)) == 0){
        filtered_gene <- NULL
    } else {
        filtered_gene <- names(which(qc_features))
    }


    # add corr group for cell type copula
    covariate_df <- covariate_df %>%
        { if(ct_copula) dplyr::mutate(., corr_group = as.integer(!!rlang::sym(cellstate_colname)))
            else dplyr::mutate(., corr_group = 1L) }

    if(!is.null(new_covariate)) {
        new_covariate <- new_covariate %>%
            { if(ct_copula) dplyr::mutate(., corr_group = as.integer(!!rlang::sym(cellstate_colname)))
                else dplyr::mutate(., corr_group = 1L) }
    }

    data_list <- list("count_mat" = count_mat,
                      "covariate" = covariate_df,
                      "new_covariate" = new_covariate,
                      "important_features" = important_features,
                      "eqtl_geno_list" = eqtl_geno_list,
                      "new_eqtl_geno_list" = new_eqtl_geno_list,
                      "filtered_gene" = filtered_gene)

    return(data_list)
}



calculateSnpLD <- function(geno_mat) {
    # FUN: helper function in pruneSnp() used to compute correlation matrix with
    # Pearson correlation as measure and outputs first vector column
    # geno_mat : SNP by sample matrix with additive genotype

    # TODO: fix warning for zero std output

    out <- stats::cor(t(geno_mat), method = "pearson")

    # output column vector of absolute value of correlation
    as.vector(abs(out[1, ]))
}



pruneSnp <- function(geno_mat, ld_threshold = 0.9) {
    # FUN : greedily prunes SNPs and outputs row indices for SNPs to keep
    # geno_mat : snp by sample matrix with additive genotype
    # ld_threshold : numeric scalar (between 0 and 1) used to filter out SNPs
    #                whose absolute Pearson correlation value exceeds this

    selected_snps <- 1:nrow(geno_mat)

    # intialize vector to keep snps
    to_keep <- c()
    i <- 1

    while(i <= length(selected_snps)) {

        current_snp <- selected_snps[1]
        remaining_snps <- setdiff(selected_snps, current_snp)

        if(is.null(remaining_snps) || length(remaining_snps) == 0) {
            to_keep <- c(to_keep, current_snp)
            break
        }

        # compute LD between current SNP and remaining SNPs
        ld_values <- calculateSnpLD(geno_mat[c(current_snp, remaining_snps), , drop = FALSE])

        # filter SNPs based on LD threshold
        to_remove <- remaining_snps[ld_values[-1] > ld_threshold]  # exclude current SNP from check

        # Remove the SNPs that exceed the LD threshold
        selected_snps <- setdiff(selected_snps, to_remove)

        to_keep <- c(to_keep, current_snp)

        # Remove current SNP from set
        selected_snps <- setdiff(selected_snps, current_snp)
    }

    return(to_keep)
}
