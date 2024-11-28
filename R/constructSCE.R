#' Construct a SingleCellExperiment object with specified covariates.
#'
#' @param data_obj a Seurat, SingleCellExperiment, or expression matrix (gene-by-cell)
#'     used to construct the new SingleCellExperiment object.
#' @param cellcov_df a dataframe containing covariates for all cells. Row names must
#'     have cell names that match \code{data_obj}. Only used when \code{data_obj}
#'     is a matrix object.
#' @param featcov_df a dataframe containing covariates for all features (ie genes).
#'     Row names must have feature names that match \code{data_obj}. Only used when
#'     \code{data_obj} is a matrix object.
#' @param cellcov_names a string vector of all columns to extract from cell covariates.
#' @param sampid_vec an optional string vector of sample ids (individuals). Default
#'     is NULL. If provided, it must match sample ids found in \code{data_obj}.
#' @param overlap_features an optional string vector for features to filter \code{data_obj}.
#'     Default is NULL. If provided, all features must be present in the \code{data_obj}.
#' @param assay_name a string scalar of name of the Seurat object slot of interest.
#'     Only used if \code{data_obj} is a Seurat object.
#' @param slot_name a string scalar for the type of assay data (either 'logcounts'
#'     or 'counts'). Default is 'counts'. Only used if \code{data_obj} is a Seurat or
#'     SingleCellExperiment object.
#' @param sce_name a string scalar for the type of assay data for output SingleCellExperiment.
#' @param indiv_colname a string scalar for the sample id colname. Must be present
#'     in \code{cellcov_df} or \code{data_obj}. Default is 'indiv'.
#' @param cellstate_colname a string scalar for the cell state colname. Must be present
#'     in \code{cellcov_df} or \code{data_obj}. Default is 'cell_type'.
#' @param factor_colnames an optional string vector or scalar for the columns to
#'     coerce into factor variables. Default is NULL. If provided, it must be present
#'     in \code{cellcov_df} or \code{data_obj}.
#' @param cellcov_renames an optional string vector to rename the cell covariates.
#'     Must be same length as \code{cellcov_names}. Default is same as \code{cellcov_names}.
#'
#' @return a SingleCellExperiment object
#' @export
#'
#' @examples
#' NULL
constructSCE <- function(data_obj,
                         cellcov_df = NULL,
                         featcov_df = NULL,
                         cellcov_names,
                         sampid_vec,
                         overlap_features = NULL,
                         assay_name = "RNA",
                         slot_name = "counts",
                         sce_name = slot_name,
                         indiv_colname = "indiv",
                         cellstate_colname = "cell_type",
                         factor_colnames = NULL,
                         cellcov_renames = cellcov_names) {

    # checks
    stopifnot("cellcov_renames length is not same as cellcov_names. Please check input!" =
                  (length(cellcov_renames) == length(cellcov_names)))

    stopifnot("indiv_colname is not included in cellcov_names. Please check input!" =
                  (indiv_colname %in% cellcov_names))

    stopifnot("cellstate_colname is not included in cellcov_names. Please check input!" =
                  (cellstate_colname %in% cellcov_names))

    stopifnot("factor_colnames is not included in cellcov_names. Please check input!" =
                  (factor_colnames %in% cellcov_names))

    # extract column types, indiv ids and feature names
    if(methods::is(data_obj, "Seurat")) {

        seurat_res <- extractFromSeurat(data_obj,
                                        features = NULL,
                                        assay_name = assay_name,
                                        slot_name = slot_name,
                                        col_class = TRUE,
                                        feat_names = TRUE,
                                        sc_indiv = TRUE,
                                        indiv_colname = indiv_colname)

        col_class <- seurat_res[["col_class"]]
        sc_indiv <- seurat_res[["sc_indiv"]]
        feat_names <- seurat_res[["feat_names"]]

    } else if(methods::is(data_obj, "SingleCellExperiment")) {

        sce_res <- extractFromSCE(data_obj,
                                  features = NULL,
                                  slot_name = slot_name,
                                  col_class = TRUE,
                                  feat_names = TRUE,
                                  sc_indiv = TRUE,
                                  indiv_colname = indiv_colname)

        col_class <- sce_res[["col_class"]]
        sc_indiv <- sce_res[["sc_indiv"]]
        feat_names <- sce_res[["feat_names"]]

    } else if(methods::is(data_obj, "matrix") || methods::is(data_obj, "dgCMatrix")) {

        col_class <- base::sapply(cellcov_df, class)

        sc_indiv <- cellcov_df[[indiv_colname]] %>%
            unique(.) %>%
            as.character(.)

        feat_names <- rownames(data_obj)

    } else {
        stop(sprintf("data_obj must be one of Seurat, SingleCellExperiment, or matrix type. Please check input!"))
    }

    # checks
    stopifnot("cellstate_colname is not found. Please check input!" =
                  cellstate_colname %in% names(col_class))

    stopifnot("indiv_colname is not found. Please check input!" =
                  indiv_colname %in% names(col_class))

    if(!is.null(sampid_vec)) {
        stopifnot("Sample ids in data_obj and sampid_vec are not same. Please check input!" =
                      checkVectorEqual(sampid_vec, sc_indiv, ignore_order = TRUE))
    }


    # check all features in overlap_features exist
    if(sum(feat_names %in% overlap_features) != length(overlap_features)) {
        stop(sprintf("Number of features in overlap_features file are not same as in Seurat object. Please check input!"))
    }

    # extract expression matrix, cell covariates, and feature covariates
    if(methods::is(data_obj, "Seurat")) {

        seurat_list <- extractFromSeurat(seurat_obj = data_obj,
                                         features = overlap_features,
                                         assay_name = assay_name,
                                         slot_name = slot_name,
                                         expr_mat = TRUE,
                                         sc_cov = TRUE,
                                         feat_cov = TRUE)
        # Note: currently works on Seurat v4 objects

        data_obj <- seurat_list[["expr_mat"]]
        cellcov_df <- seurat_list[["sc_cov"]]
        featcov_df <- seurat_list[["feat_cov"]]

    } else if(methods::is(data_obj, "SingleCellExperiment")) {

        sce_list <- extractFromSCE(sce_obj = data_obj,
                                   features = overlap_features,
                                   slot_name = slot_name,
                                   expr_mat = TRUE,
                                   sc_cov = TRUE,
                                   feat_cov = TRUE)

        data_obj <- sce_list[["expr_mat"]]
        cellcov_df <- sce_list[["sc_cov"]]
        featcov_df <- sce_list[["feat_cov"]]

    } else if(methods::is(data_obj, "matrix") || methods::is(data_obj, "dgCMatrix")) {

        # checks
        stopifnot("cellcov_df is missing. Please check input!" = !is.null(cellcov_df))
        stopifnot("featcov_df is missing. Please check input!" = !is.null(featcov_df))

        stopifnot("Cell ids in data_obj and cellcov_df must be the same. Please check input!" =
                      checkVectorEqual(rownames(cellcov_df), colnames(data_obj), ignore_order = FALSE))

        stopifnot("Feature ids in data_obj and featcov_df must be the same. Please check input!" =
                      checkVectorEqual(rownames(featcov_df), rownames(data_obj), ignore_order = FALSE))
    }

    # subset columns in cell covariates
    cellcov_df <- cellcov_df %>%
        dplyr::select(dplyr::all_of(cellcov_names))


    # explicitly convert columns to factor variables
    if(!is.null(factor_colnames)) {
        cellcov_df[factor_colnames] <- base::lapply(cellcov_df[factor_colnames], function(x) {
            as.character(x) %>%
                as.factor(.)
        })
    }

    # convert remaining character columns to factor variables
    ischar_cols <- base::sapply(cellcov_df, is.character)
    cellcov_df[ischar_cols] <- base::lapply(cellcov_df[ischar_cols], as.factor)

    # rename cell column names
    if(!identical(cellcov_renames, cellcov_names)) {
        cellcov_df <- cellcov_df %>%
            dplyr::rename_with(~ cellcov_renames, tidyselect::everything())
    }

    # construct new sce object
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(data_obj),
                                                      colData = cellcov_df)

    SummarizedExperiment::assayNames(sce) <- sce_name

    # remove non-sparse count matrix in metadata
    sce@rowRanges@elementMetadata@listData <- list()

    # add feature covariate as metadata
    SummarizedExperiment::rowData(sce) <- featcov_df

    return(sce)
}


#' Extract data from Seurat object
#'
#' @param seurat_obj a Seurat object
#' @param features string vector to filter for features
#' @param assay_name string scalar to specify the name of input assay data
#' @param slot_name string scalar to specify the type of input assay data
#' @param indiv_colname string scalar to specify the column that has the indiv ids
#' @param ... other options to specify extracted outputs. Currently available:
#' \itemize{
#'      \item \code{col_class} : logical value
#'      \item \code{feat_names} : logical value
#'      \item \code{cell_names} : logical value
#'      \item \code{sc_indiv} : logical value
#'      \item \code{expr_mat} : logical value
#'      \item \code{sc_cov} : logical value
#'      \item \code{feat_cov} : logical value
#' }
#'
#' @return a list of extracted data with following optional outputs:
#' \itemize{
#'      \item \code{col_class} : a vector of the column classes for cell covariate variables
#'      \item \code{feat_names} : a vector of feature names
#'      \item \code{cell_names} : a vector of cell names
#'      \item \code{sc_indiv} : a vector of distinct sample ids extracted from cell covariate using \code{indiv_colname}
#'      \item \code{expr_mat} : an expression matrix extracted using \code{assay_name} and \code{slot_name} options
#'      \item \code{sc_cov} : a dataframe for cell covariates
#'      \item \code{feat_cov} : a dataframe for feature covariates
#' }
#' @export
#'
#' @examples
#' NULL
extractFromSeurat <- function(seurat_obj,
                              features = NULL,
                              assay_name = "RNA",
                              slot_name = "counts",
                              indiv_colname = "indiv",
                              ...) {

    opt_list <- list(...)

    # filter for features
    if(!is.null(features)) {
        seurat_obj <- seurat_obj[features, ]
    }

    res_list <- list()

    # extract cell cov column types, feature names, cell names, and indiv ids
    if(isTRUE(opt_list[["col_class"]])) {
        res_list[["col_class"]] <- base::sapply(seurat_obj@meta.data, class)
    }

    if(isTRUE(opt_list[["feat_names"]])) {
        res_list[["feat_names"]] <- rownames(seurat_obj)
    }

    if(isTRUE(opt_list[["cell_names"]])) {
        res_list[["cell_names"]] <- colnames(seurat_obj)
    }

    if(isTRUE(opt_list[["sc_indiv"]])) {
        res_list[["sc_indiv"]] <- seurat_obj@meta.data[[indiv_colname]] %>%
            unique(.) %>%
            as.character(.)
    }

    # extract gene-by-cell expression matrix
    if(isTRUE(opt_list[["expr_mat"]])) {
        if(slot_name == "logcounts") {
            res_list[["expr_mat"]] <- SeuratObject::GetAssayData(seurat_obj,
                                                                 slot = "data",
                                                                 assay = assay_name)
        } else if(slot_name == "counts") {
            res_list[["expr_mat"]] <- SeuratObject::GetAssayData(seurat_obj,
                                                                 slot = "counts",
                                                                 assay = assay_name)
        } else if(slot_name == "data") {
            res_list[["expr_mat"]] <- SeuratObject::GetAssayData(seurat_obj,
                                                                 slot = "data",
                                                                 assay = assay_name)
        } else {
            stop("slot_name option must be either 'logcounts', 'counts', or 'data' and available in Seurat object. Please check!")
        }

    }

    # extract cell and feature covariates
    if(isTRUE(opt_list[["sc_cov"]])) {
        res_list[["sc_cov"]] <- seurat_obj@meta.data
    }

    if(isTRUE(opt_list[["feat_cov"]])) {
        res_list[["feat_cov"]] <- seurat_obj@assays[[assay_name]]@meta.features
    }

    res_list
}


#' Extract data from SingleCellExperiment object
#'
#' @param sce_obj a SingleCellExperiment object
#' @param features string vector to filter for features
#' @param slot_name string scalar to specify the type of input assay data
#' @param indiv_colname string scalar to specify the column that has the indiv ids
#' @param ... other options to specify extracted outputs
#'
#' @return a list of extracted data
#' @export
#'
#' @examples
#' NULL
extractFromSCE <- function(sce_obj,
                           features = NULL,
                           slot_name = "counts",
                           indiv_colname = "indiv",
                           ...) {

    opt_list <- list(...)

    # filter for features
    if(!is.null(features)) {
        sce_obj <- sce_obj[features, ]
    }

    res_list <- list()

    # extract cell cov column types, feature names, cell names, and indiv ids
    if(isTRUE(opt_list[["col_class"]])) {
        res_list[["col_class"]] <- base::sapply(SingleCellExperiment::colData(sce_obj), class)
    }

    if(isTRUE(opt_list[["feat_names"]])) {
        res_list[["feat_names"]] <- rownames(sce_obj)
    }

    if(isTRUE(opt_list[["cell_names"]])) {
        res_list[["cell_names"]] <- colnames(sce_obj)
    }

    if(isTRUE(opt_list[["sc_indiv"]])) {
        res_list[["sc_indiv"]] <- SingleCellExperiment::colData(sce_obj)[[indiv_colname]] %>%
            unique(.) %>%
            as.character(.)
    }

    # extract gene-by-cell expression matrix
    if(isTRUE(opt_list[["expr_mat"]])) {
        if(slot_name == "logcounts") {
            res_list[["expr_mat"]] <- SingleCellExperiment::logcounts(sce_obj)
        } else if(slot_name == "counts") {
            res_list[["expr_mat"]] <- SingleCellExperiment::counts(sce_obj)
        } else if(slot_name == "normcounts") {
            res_list[["expr_mat"]] <- SingleCellExperiment::normcounts(sce_obj)
        } else {
            stop("slot_name option must be either 'logcounts', 'counts', or 'normcounts' and available in SingleCellExperiment object. Please check input!")
        }
    }

    # extract cell and feature covariates
    if(isTRUE(opt_list[["sc_cov"]])) {
        res_list[["sc_cov"]] <- SingleCellExperiment::colData(sce_obj) %>%
            as.data.frame(.)
    }

    if(isTRUE(opt_list[["feat_cov"]])) {
        res_list[["feat_cov"]] <- SingleCellExperiment::rowData(sce_obj) %>%
            as.data.frame(.)
    }

    res_list
}
