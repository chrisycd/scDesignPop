#' Dimensionality reduction and visualization for population-scale data
#'
#' This function takes a reference sce and a list of new sces,
#' performs the dimensionality reduction on the reference data,
#' projects the synthetic datasets on the same low dimensional
#' space (PCA and UMAP) of the reference sce, and then visualize the results.
#'
#' @param ref_sce a reference sce object which synthetic sce objects will
#' be projected to.
#' @param sce_list a list of synthetic sce objects.
#' @param name_vec a string vector specifiying the names of each dataset.
#' The length should be \code{length(sce_list) + 1}, where the first name is
#' for \code{ref_sce}.
#' @param assay_use a string scalar which indicates the assay you will use
#' in the sce. Default is 'logcounts'.
#' @param n_pc an integer specifying the number of PCs.
#' @param pc_umap a boolean value specifying whether using PCs as the input of UMAP.
#' Default is TRUE.
#' @param center a boolean value specifying whether centering the data before PCA.
#' Default is TRUE.
#' @param scale. a boolean value specifying whether scaling the data before PCA.
#' Default is TRUE.
#' @param if_plot a boolean value specifying whether returning the plot.
#' If FALSE, return the reduced dimensions of each dataset.
#' @param shape_by a string scalar which indicates the column in \code{colData}
#' used for shape.
#' @param color_by a string scalar which indicates the column in \code{colData}
#' used for color.
#' @param point_size an numeric scalar specifying of the point size in the
#' final plot. Default is 1.
#'
#' @return A data frame of the reduced dimensions (both PCA and UMAP) or a list
#' contains both the data frame and two ggplot2 object of PCA plot and UMAP plot.
#'
#' @import ggplot2
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom irlba irlba
#' @importFrom uwot umap
#' @importFrom viridis viridis
#'
#' @export
#'
#' @examples
#' NULL
plotReducedDimPop <- function (
    ref_sce, sce_list, name_vec, assay_use = "logcounts",
    pc_umap = TRUE, n_pc = 50, center = TRUE, scale. = TRUE,
    if_plot = TRUE, shape_by = NULL, color_by, point_size = 1)
{
  Method <- NULL
  stopifnot(length(name_vec) == (length(sce_list) + 1))
  ref_sce <- ref_sce[,order(SummarizedExperiment::colData(ref_sce)[,color_by])]
  mat_ref <- t(as.matrix(SummarizedExperiment::assay(ref_sce,
                                                     assay_use)))
  if (sum(matrixStats::colVars(mat_ref) == 0) > 0) {
    stop("The ref dataset contains 0 variance features. Please remove them.")
  }
  if (!is.null(shape_by)) {
    if (!(shape_by %in% colnames(SummarizedExperiment::colData(ref_sce)))) {
      stop("The shape_by in not in your ref_sce's colData. Please double check the variable name for shape_by.")
    }
    shape_by_check <- sapply(sce_list, function(x) {
      shape_by %in% colnames(SummarizedExperiment::colData(x))
    })
    if (!all(shape_by_check)) {
      stop("The shape_by in not in your sce_list's colData. Please double check the variable name for shape_by.")
    }
  }
  mat_list <- lapply(sce_list, function(x) {
    mat <- t(as.matrix(SummarizedExperiment::assay(x, assay_use)))
    mat
  })

  RhpcBLASctl::blas_set_num_threads(4)
  ref_pca_fit <- irlba::prcomp_irlba(mat_ref, center = center,
                                     scale. = scale., n = n_pc)
  ref_pca <- ref_pca_fit$x
  if (pc_umap) {
    # ref_umap_fit <- umap::umap(ref_pca_fit$x)
    uwot_model <- uwot::umap(
      ref_pca_fit$x,
      n_neighbors = 15,
      min_dist    = 0.1,
      ret_model   = TRUE,       # return full model for transforms
      n_threads   = 4           # parallelism for fit
    )
  }
  else {
    # ref_umap_fit <- umap::umap(mat_ref)
    uwot_model <- uwot::umap(
      mat_ref,
      n_neighbors = 15,
      min_dist    = 0.1,
      ret_model   = TRUE,       # return full model for transforms
      n_threads   = 4           # parallelism for fit
    )
  }
  ref_umap <- uwot_model$embedding
  colnames(ref_umap) <- c("UMAP1", "UMAP2")
  SingleCellExperiment::reducedDim(ref_sce, "PCA") <- ref_pca
  SingleCellExperiment::reducedDim(ref_sce, "UMAP") <- ref_umap
  sce_list <- lapply(sce_list, function(x) {
    x <- x[,order(colData(x)[,color_by])]
    mat <- t(as.matrix(SummarizedExperiment::assay(x, assay_use)))
    SingleCellExperiment::reducedDim(x, "PCA") <- stats::predict(ref_pca_fit,
                                                                 newdata = mat)
    if (pc_umap) {
      # res <- stats::predict(object = ref_umap_fit, data = SingleCellExperiment::reducedDim(x,"PCA"))
      res <- uwot::umap_transform(
        SingleCellExperiment::reducedDim(x, "PCA"),
        model     = uwot_model,
        n_threads = 4            # parallelism for transform
      )
    }
    else {
      # res <- stats::predict(object = ref_umap_fit, data = mat)
      res <- uwot::umap_transform(
        mat,
        model     = uwot_model,
        n_threads = 4            # parallelism for transform
      )
    }
    colnames(res) <- c("UMAP1", "UMAP2")
    SingleCellExperiment::reducedDim(x, "UMAP") <- res
    return(x)
  })
  sce_list_new <- c(list(ref_sce), sce_list)
  names(sce_list_new) <- name_vec
  rd_list <- lapply(sce_list_new, function(x) {
    rd <- tibble::as_tibble(SummarizedExperiment::colData(x))
    if (is.null(shape_by)) {
      rd <- dplyr::select(rd, color_by)
    }
    else {
      rd <- dplyr::select(rd, c(color_by, shape_by))
    }
    rd_pca <- tibble::as_tibble(SingleCellExperiment::reducedDim(x,
                                                                 "PCA"))
    rd_umap <- tibble::as_tibble(SingleCellExperiment::reducedDim(x,
                                                                  "UMAP"))
    rd <- dplyr::bind_cols(rd, rd_pca)
    rd <- dplyr::bind_cols(rd, rd_umap)
    rd
  })
  names(rd_list) <- names(sce_list_new)
  rd_tbl <- dplyr::bind_rows(rd_list, .id = "Method")
  rd_tbl <- dplyr::mutate(rd_tbl, Method = factor(Method,
                                                  levels = name_vec))

  ## switch back
  RhpcBLASctl::blas_set_num_threads(1)
  ##
  if (if_plot) {
    p_pca <- ggplot(rd_tbl, aes_string(x = "PC1", y = "PC2",
                                       color = color_by)) + geom_point(alpha = 0.8, size = point_size,
                                                                       aes_string(shape = shape_by)) + facet_wrap(~Method,
                                                                                                                  nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom")
    p_umap <- ggplot(rd_tbl, aes_string(x = "UMAP1", y = "UMAP2",
                                        color = color_by)) + geom_point(alpha = 0.8, size = point_size,
                                                                        aes_string(shape = shape_by)) + facet_wrap(~Method,
                                                                                                                   nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom")
    if (is.numeric(unlist(rd_tbl[, color_by]))) {
      p_pca <- p_pca + viridis::scale_color_viridis()
      p_umap <- p_umap + viridis::scale_color_viridis()
    }
    else {
      p_pca <- p_pca + guides(color = guide_legend(override.aes = list(size = 2,
                                                                       alpha = 1)))
      p_umap <- p_umap + guides(color = guide_legend(override.aes = list(size = 2,
                                                                         alpha = 1)))
    }
    return(list(p_pca = p_pca, p_umap = p_umap, rd_tbl = rd_tbl))
  }
  else {
    return(rd_tbl)
  }
}
