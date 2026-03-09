#' Creates pseudobulk expression grouped by eQTL genotypes.
#'
#' Computes pseudobulk based on normalization and aggregation options, and combines
#' with eQTL genotypes to output dataframe of pseudobulk expression vs genotypes.
#'
#' @param sce_list A list of named SingleCellExperiment (SCE) objects. Each SCE
#'     object must have the \code{indiv_colname} and \code{celltype_colname}
#'     columns in its colData. If no names are provided in the list, the objects
#'     will be named "sce_1", "sce_2", etc, by default.
#' @param eqtlgeno A dataframe of eQTL with genotypes for samples. Must have the
#'     columns corresponding to \code{feature_colname}, \code{snp_colname},
#'     \code{celltype_colname}, and \code{loc_colname}. The \code{loc_colname}
#'     variable must be the last column that precedes the genotypes of the
#'     individuals (samples).
#' @param feature_sel A string scalar for the feature name (eg. a gene id) to
#'     select. Must exist in rowname of every SCE object in \code{sce_list}.
#' @param celltype_sel A string scalar or vector of the cell states (eg. cell types)
#'     to select. Must be one of the cell types in the \code{celltype_colname} of
#'     every SCE object in \code{sce_list}.
#' @param eqtl_indx An integer scalar or vector to select the eQTL from \code{eqtlgeno}
#'     for a given celltype and feature. Only used if \code{eqtl_snp = NULL}. The default is 1.
#'     The values correspond to \code{celltype_sel} and are recycled if the length of
#'     elements is shorter.
#' @param eqtl_snp A string scalar or vector to select the eQTL from \code{eqtlgeno}.
#'     Overrides \code{eqtl_indx} option if both \code{eqtl_indx} and \code{eqtl_snp}
#'     are specified. Can select any SNP in the \code{eqtlgeno} (even from another
#'     feature). Must exist in the \code{snp_colname} column of \code{eqtlgeno}.
#'     The values correspond to \code{celltype_sel} and are recycled if the length of
#'     elements is shorter.
#' @param normalize_type A string scalar for normalization method applied to
#'     \code{sce_list}.
#' @param aggregate_type A string scalar for pseudobulk aggregation applied to
#'     \code{sce_list} after subsetting by \code{celltype_sel}.
#' @param slot_name A string scalar for the SCE slot used. The default is "counts".
#' @param rescale A logical scalar for whether to rescale the response values
#'     based on maximum expression value.
#' @param overwrite A logical value on whether to overwrite the "logcounts" slot in
#'     input sce object. Must be set to \code{overwrite = TRUE} to overwrite. The
#'     default is FALSE.
#' @param feature_colname A string scalar for column name of feature variable
#'     (i.e. gene). The default is "gene_id".
#' @param snp_colname A string scalar for column name of the SNP variable. The
#'     default is "snp_id".
#' @param indiv_colname A string scalar for column name of individuals (samples).
#'     The default is "indiv".
#' @param celltype_colname A string scalar for column name of cell type. The
#'     default is "cell_type".
#' @param loc_colname A string scalar for column name of SNP position variable.
#'     The default is "POS".
#' @param if_plot A logical scalar specifying whether to return a pseudobulk plot.
#'     The default is TRUE.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return Outputs a dataframe with following columns:
#' \describe{
#'      \item{\code{indiv}}{Name of individual's sample ID.}
#'      \item{\code{response}}{Pseudobulk expression of given feature.}
#'      \item{\code{genotype}}{Genotype values of individual for selected eQTL SNP.}
#'      \item{\code{feature}}{Name of selected feature.}
#'      \item{\code{celltype}}{Name of selected cell type.}
#'      \item{\code{snp}}{Name of selected eQTL SNP.}
#'      \item{\code{sce_name}}{Name of the SingleCellExperiment object.}
#'      \item{\code{normalization}}{Type of normalization method used.}
#'      \item{\code{aggregate_type}}{Type of pseudobulk aggregation used.}
#'      \item{\code{slot_used}}{Slot used from input SingleCellExperiment object.}
#' }
#'
#' @export
#'
#' @examples
#' NULL
createPbulkExprGeno <- function(sce_list,
                                eqtlgeno,
                                feature_sel,
                                celltype_sel = "allct",
                                eqtl_indx = 1L,
                                eqtl_snp = NULL,
                                normalize_type = c("log1p", "seurat", "scran", "none"),
                                aggregate_type = c("mean", "sum"),
                                slot_name = c("counts", "logcounts"),
                                rescale = FALSE,
                                overwrite = FALSE,
                                feature_colname = "gene_id",
                                snp_colname = "snp_id",
                                indiv_colname = "indiv",
                                celltype_colname = "cell_type",
                                loc_colname = "POS",
                                if_plot = TRUE,
                                ...) {
    # Note: uses normalizeSCE, createPbulk, getCommonCelltypes, subsetSCEbyCelltype,
    #       recycleTo, getCommonFeatures functions

    normalize_type <- match.arg(normalize_type)
    aggregate_type <- match.arg(aggregate_type)
    slot_name <- match.arg(slot_name)

    # check
    assertthat::assert_that(
        assertthat::has_name(eqtlgeno,
                             c(feature_colname, snp_colname,
                               celltype_colname, loc_colname))
    )

    # check normalization and slot used
    if(normalize_type %in% c("scran", "seurat", "log1p") && slot_name == "counts") {
        stop("normalize_type and slot_name option mismatch. Please check!")
    }

    # check SCE names
    sce_names <- names(sce_list)
    if(is.null(sce_names)) { sce_names <- paste0("sce_", 1:length(sce_list)) }

    celltypes <- getCommonCelltypes(sce_list = sce_list,
                                    celltype_colname = celltype_colname)

    # check cell types
    celltypes_only <- celltype_sel[grep("allct", celltype_sel, invert = TRUE)]

    stopifnot("celltype_sel must be present in sce_list. Please check!" =
                  all(celltypes_only %in% celltypes))

    features <- getCommonFeatures(sce_list)

    # check
    stopifnot("feature_sel must be present in sce_list. Please check!" =
                  (feature_sel %in% features))

    # indices to subset geno columns
    firstcol_indx <- which(colnames(eqtlgeno) == loc_colname) + 1L
    lastcol_indx <- dim(eqtlgeno)[2]

    ## step 1. Normalize each SCE list object
    sce_list <- base::lapply(sce_list,
                             normalizeSCE,
                             method = normalize_type,
                             overwrite = overwrite,
                             ...)

    logic_vec <- sapply(sce_list, FUN = function(sce) {
        assay_names <- SummarizedExperiment::assayNames(sce)
        return(slot_name %in% assay_names)
    })

    if(!all(logic_vec)) {
        stop(sprintf("%s assay is not found in some SCE in sce_list. Please check!", slot_name))
    }

    ## step 2. create eqtl vectors matched to same length as input cell type selection
    n_ct <- length(celltype_sel)
    if(!is.null(eqtl_snp)) {

        if(!is.null(eqtl_indx)) {
            message("eqtl_snp input used to override eqtl_indx...")
        }

        eqtl_snp <- recycleTo(eqtl_snp, n_ct, "eqtl_snp")

        eqtl_snp_list <- sapply(eqtl_snp, list, USE.NAMES = FALSE)
        names(eqtl_snp_list) <- celltype_sel

        eqtl_idx_list <- NULL

    } else if(!is.null(eqtl_indx) && is.null(eqtl_snp)) {

        eqtl_indx <- recycleTo(eqtl_indx, n_ct, "eqtl_indx")

        eqtl_idx_list <- sapply(eqtl_indx, list, USE.NAMES = FALSE)
        names(eqtl_idx_list) <- celltype_sel

        eqtl_snp_list <- NULL

    } else {
        stop("Either 'eqtl_indx' or 'eqtl_snp' option must be specified. Please check input!")
    }

    # iterate thru cell types and apply to sce_list
    res_list <- lapply(celltype_sel, FUN = function(ct) {

        ## step 3. subset by cell type
        sce_list <- base::lapply(sce_list,
                                 FUN = subsetSCEbyCelltype,
                                 ct = ct,
                                 celltype_colname = celltype_colname)

        ## step 4. aggregate to create pseudobulk-by-indiv list
        pbulk_list <- base::lapply(sce_list, FUN = function(sce) {
            mat <- as.data.frame(createPbulk(sce = sce,
                                             slot_name = slot_name,
                                             aggregate_type = aggregate_type))
            return(mat)
        })
        names(pbulk_list) <- sce_names

        # get indiv colnames
        indivname_list <- lapply(pbulk_list, FUN = colnames)
        names(indivname_list) <- sce_names

        ## step 5. get eQTL of selected gene and cell type
        if(is.null(eqtl_snp_list)) {

            eqtl_sub <- dplyr::filter(eqtlgeno,
                                      !!rlang::sym(feature_colname) == feature_sel)

            if(ct != "allct") {
                eqtl_sub <- dplyr::filter(eqtl_sub,
                                          !!rlang::sym(celltype_colname) == ct)
            }

            eqtls <- unique(eqtl_sub[[snp_colname]])

            if(length(eqtls) == 0) {
                stop(sprintf("There are no eQTLs for selected feature. Please check!"))
            }

            message(sprintf("%s eQTLs available in %s celltype for %s feature:\n  %s \nSelected %s snp using 'eqtl_indx=%s'.",
                            length(eqtls),
                            ct,
                            feature_sel,
                            paste(eqtls, collapse = ", "),
                            eqtls[eqtl_idx_list[[ct]]],
                            eqtl_idx_list[[ct]]
            ))

            eqtl_sel <- eqtl_sub[eqtl_idx_list[[ct]], ]

        } else if(is.null(eqtl_idx_list)) {

            eqtl_sub <- dplyr::filter(eqtlgeno,
                                      !!rlang::sym(snp_colname) == eqtl_snp_list[[ct]])

            if(nrow(eqtl_sub) == 0) {
                stop(sprintf("eQTL SNP must be present in eqtlgeno. Please check!"))
            }

            cts <- unique(eqtl_sub[[celltype_colname]])

            message(sprintf("%s cell types use this SNP:\n  %s \nSelected %s snp for %s celltype, %s feature.",
                            length(cts),
                            paste(cts, collapse = ", "),
                            eqtl_snp_list[[ct]],
                            ct,
                            feature_sel
            ))

            eqtl_sel <- dplyr::slice_head(eqtl_sub)
            # Note: feature in eqtl_sel does not need to be same feature_sel
        }

        ## step 6. merge genotypes of selected indivs and pseudobulk expression values
        res_list <- lapply(1:length(sce_list), FUN = function(i) {

            geno_select <- tidyr::pivot_longer(eqtl_sel[, indivname_list[[i]]],
                                               cols = tidyselect::everything(),
                                               values_to = "genotype",
                                               names_to = indiv_colname) %>%
                dplyr::mutate(genotype = as.factor(.data$genotype))

            eqtl_annot <- eqtl_sel[, 1:(firstcol_indx - 1)]

            # create a df based on cell type and pseudobulk expression
            response <- tidyr::pivot_longer(pbulk_list[[i]][feature_sel, ],
                                            tidyselect::everything(),
                                            names_to = indiv_colname,
                                            values_to = "response")

            # merge response with genotype and add other column info
            res_df <- dplyr::left_join(response,
                                       geno_select,
                                       by = indiv_colname) %>%
                dplyr::mutate(!!rlang::sym(feature_colname) := feature_sel,
                              !!rlang::sym(celltype_colname) := ct,
                              !!rlang::sym(snp_colname) := eqtl_annot[[snp_colname]],
                              sce_name = sce_names[i],
                              normalization = normalize_type,
                              aggregate_type = aggregate_type,
                              slot_used = slot_name)

            # re-scale response
            if(rescale) {
                res_df <- dplyr::mutate(res_df,
                                        response = response / max(response))
            }

            return(res_df)
        })

        return(do.call(rbind, res_list))  # combine all SCEs
    })

    res_df <- do.call(rbind, res_list)  # combine all cell types

    if(if_plot) {
        p_pbulk <- plotPbulkGeno(res_df = res_df,
                               feature_sel = feature_sel,
                               celltype_colname = celltype_colname,
                               snp_colname = snp_colname,
                               ...
                               )
        # print(p_pbulk)

        return(list("res_df" = res_df,
                    "p_pbulk" = p_pbulk))

    } else {
        return(list("res_df" = res_df))
    }
}



#' Plot pseudobulk expression vs. genotype
#'
#' Creates a faceted violin plot of pseudobulk expression values across genotype
#' groups, with optional boxplots and mean summary points.
#'
#' @param res_df A dataframe
#' @param show_boxplot A logical scalar for whether to overlay boxplots.
#' @param show_meanstat A logical scalar for whether to display the mean marked
#'     by a red dot.
#' @param show_snp A logical scalar for whether to display selected SNP for each
#'     y-axis facet.
#' @param ... Additional arguments passed to internal functions.
#' @inheritParams createPbulkExprGeno
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' NULL
plotPbulkGeno <- function(res_df,
                          celltype_colname,
                          snp_colname,
                          show_boxplot = TRUE,
                          show_meanstat = TRUE,
                          show_snp = FALSE,
                          ...) {

    # include SNP info in y-axis facet
    if(show_snp) {
        res_df <- dplyr::mutate(res_df,
                                !!rlang::sym(celltype_colname) :=
                                    paste0(!!rlang::sym(celltype_colname), "\n",
                                           !!rlang::sym(snp_colname)))
    }

    p_pbulk <- ggplot2::ggplot(res_df,
                               ggplot2::aes(x = .data$genotype,
                                            y = .data$response)
                               ) +
        ggplot2::geom_violin(ggplot2::aes(fill = .data$genotype),
                             trim = TRUE,  # cutoff at extremes
                             linewidth = 1.0,
                             position = ggplot2::position_dodge(1),
                             scale = "width") +
        ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(celltype_colname)),
                            cols = ggplot2::vars(.data$sce_name),
                            scales = "free_y",
                            drop = FALSE)
    if(show_boxplot) {
        p_pbulk <- p_pbulk +
            ggplot2::geom_boxplot(ggplot2::aes(group = .data$genotype),
                                  width = 0.2,
                                  position = ggplot2::position_dodge(1))
    }

    if(show_meanstat) {
        p_pbulk <- p_pbulk +
            ggplot2::stat_summary(fun = mean,
                                  geom = "point",
                                  shape = 16,
                                  size = 2,
                                  color = "red",
                                  fill = "red")
    }

    return(p_pbulk)
}


################# Helper Functions #################

# aggregates single-cell data into pseudobulk expression
createPbulk <- function(sce,
                        aggregate_type = c("sum", "mean"),
                        slot_name = "counts",
                        indiv_colname = "indiv") {
    # sce : a SingleCellExperiment object
    # aggregate_type : aggregation option. Options are either "sum" or "mean"
    # slot_name : string scalar of the assay name to extract from sce_obj
    # indiv_colname : a string scalar for column name of individuals (samples)

    aggregate_type <- match.arg(aggregate_type)

    # check
    assertthat::assert_that(assertthat::has_name(SingleCellExperiment::colData(sce),
                                                 indiv_colname))

    # aggregate counts by individual
    indivNames <- as.character(unique(SingleCellExperiment::colData(sce)[[indiv_colname]]))

    pseudo_bulk <- c()

    for(i in indivNames) {
        logic_vec <- SingleCellExperiment::colData(sce)[[indiv_colname]] == i  # columns matching indiv

        if(sum(logic_vec) > 1){  # if >1 cells exist
            if(aggregate_type == "sum") {  # aggregate raw counts
                countAggregateIndiv <- SummarizedExperiment::assay(sce, slot_name)[, logic_vec] %>%
                    Matrix::rowSums()
            }

            if(aggregate_type == "mean") {  # aggregate normalized counts
                countAggregateIndiv <- SummarizedExperiment::assay(sce, slot_name)[, logic_vec] %>%
                    Matrix::rowMeans()
            }

        } else {  # if only 1 cell exists
            countAggregateIndiv <- SummarizedExperiment::assay(sce, slot_name)[, logic_vec]
        }

        pseudo_bulk <- cbind(pseudo_bulk, countAggregateIndiv)
    }
    colnames(pseudo_bulk) <- indivNames

    return(pseudo_bulk)
}


# get overlap cell types
getCommonCelltypes <- function(sce_list,
                               celltype_colname) {
    ct_list <- lapply(sce_list, function(sce) {
        as.character(unique(SummarizedExperiment::colData(sce)[[celltype_colname]]))
    })
    return(base::Reduce(base::intersect, ct_list))
}

# subset SingleCellExperiment object by cell type
subsetSCEbyCelltype <- function(sce,
                                ct,
                                celltype_colname) {

    if(ct == "allct") {
        return(sce)
    } else {
        return(sce[, sce[[celltype_colname]] == ct])
    }
}

# get intersect features
getCommonFeatures <- function(sce_list) {
    base::Reduce(base::intersect, lapply(sce_list, FUN = rownames))
}

# recycles input x to have same length given by y_n
recycleTo <- function(x, y_n, what) {
    if(length(x) == 0) {
        stop(what, " must have length >= 1.")
    } else if(length(x) < y_n) {
        message(what, " length is shorter than celltype_sel... \nRecycling elements to match length of celltype_sel.")
    } else if(length(x) > y_n) {
        message(what, " length is longer than celltype_sel... \nTruncating elements to match length of celltype_sel.")
    }
    return(rep_len(x, y_n))
}
