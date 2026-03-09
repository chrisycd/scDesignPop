#' Normalize data in a SingleCellExperiment object
#'
#' Normalizes the data in a SingleCellExperiment object with specified method
#' and returns a SingleCellExperiment object.
#'
#' ## Normalization options
#' If 'seurat' is used, the default options from \code{NormalizeData} from the
#' \code{Seurat} package is used. Seurat v4 is currently supported;
#' If 'scran' is used, the default options from \code{computePooledFactors},
#' followed by \code{logNormCounts} from the \code{scuttle} package is used.
#' If 'log1p' is used, the \code{log1p} function from base R is used. In general,
#' the computational time from fast to slow is 'log1p', 'seurat', and 'scran'.
#'
#' @param sce A SingleCellExperiment object.
#' @param method A string value specifying the normalization method. Must be one
#'     of 'log1p', seurat', 'scran', or 'none'. The default is 'log1p'.
#'     See Normalization details.
#' @param overwrite A logical value on whether to overwrite the "logcounts" slot in
#'     input sce object. Must be set to \code{overwrite = TRUE} to overwrite. The
#'     default is FALSE.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A SingleCellExperiment object
#' @export
#'
#' @examples
#' data("example_sce")
#' example_sce <- normalizeSCE(example_sce, method = "log1p", overwrite = TRUE)
normalizeSCE <- function(sce,
                         method = c("log1p", "seurat", "scran", "none"),
                         overwrite = FALSE,
                         ...) {

    method <- match.arg(method)

    has_logcount <- "logcounts" %in% names(SummarizedExperiment::assays(sce))

    if(has_logcount && !overwrite) {
        warning("The 'logcounts' slot already exists! Set overwrite = TRUE to overwrite.")
        return(sce)
    }

    if(method == "scran") {
        sce <- scranNormalize(sce,
                              ...)
    } else if(method == "seurat") {  # current works for Seurat 4.4
        sce <- seuratNormalize(sce,
                               ...)
    } else if(method == "log1p") {
        sce <- log1pNormalize(sce)
    } else if(method == "none") {
        sce <- sce
    } else {
        warning("Options for method is either 'log1p', 'seurat', 'scran', or 'none'!")
    }

    return(sce)
}


# helper function to perform scran's default size factor normalization
scranNormalize <- function(sce,
                           ...) {

    more_args <- list(...)

    message("Computing size factors for pooled cells...")

    # compute size factors across pooled cells
    sce <- scuttle::computePooledFactors(sce, positive = TRUE)   # same as scran::computeSumFactors() function

    message("Completed computing size factors!")

    # log normalize counts
    sce <- scuttle::logNormCounts(sce)

    return(sce)
}

# helper function to perform Seurat's default log1pCP10k normalization
seuratNormalize <- function(sce,
                            ...) {

    more_args <- list(...)

    scale.factor <- if(!is.null(more_args$scale.factor)) {
        more_args$scale.factor
    } else {
        10000L
    }
    normalization.method <- if(!is.null(more_args$normalization.method)) {
        more_args$normalization.method
    } else {
        "LogNormalize"
    }

    seurat <- Seurat::as.Seurat(sce, data = NULL)

    seurat <- Seurat::NormalizeData(seurat,
                                    scale.factor = scale.factor,
                                    normalization.method = normalization.method)

    return(Seurat::as.SingleCellExperiment(seurat))
}

# helper function to perform log1p normalization
log1pNormalize <- function(sce) {
    SingleCellExperiment::logcounts(sce) <- base::log1p(SingleCellExperiment::counts(sce))

    return(sce)
}
