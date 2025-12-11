#' Example single-cell RNA-seq data (multi–cell type)
#'
#' A small \code{SingleCellExperiment} object used as example input for
#' scDesignPop functions. It contains a subset of genes and cells from the
#' OneK1K cohort with associated individual and cell-type metadata.
#'
#' @format A \code{SingleCellExperiment} object with 1,000 genes and 7,811 cells.
#' \describe{
#'   \item{assays}{
#'     \itemize{
#'       \item \code{counts}: raw UMI count matrix (genes x cells).
#'       \item \code{logcounts}: log-normalized expression values.
#'     }
#'   }
#'   \item{rowData}{(empty in this toy example).}
#'   \item{colData}{A \code{DataFrame} with 4 columns:
#'     \describe{
#'       \item{\code{indiv}}{Individual identifier (e.g., \code{"SAMP1"}, \code{"SAMP2"}).}
#'       \item{\code{cell_type}}{Cell-type annotation (e.g., T cells, B cells, NK cells).}
#'       \item{\code{sex}}{Sex of the donor.}
#'       \item{\code{age}}{Age of the donor.}
#'     }
#'   }
#'   \item{reducedDimNames}{None in this example.}
#' }
#'
#' @usage data("example_sce")
#'
#' @examples
#' data("example_sce")
#' example_sce
#' table(SingleCellExperiment::colData(example_sce)$cell_type)
"example_sce"

#' Example genotype data for cell-type-specific eQTLs
#'
#' A tibble containing example genotype data for lead cis-eQTL SNPs across
#' multiple cell types and genes in OneK1K cohort. Each row corresponds to a cell-type–gene–SNP
#' combination, and each sample column stores genotype dosages (0/1/2) for one
#' individual. No eQTL effect size results are included. Genotypes are all permuted.
#'
#' @format A tibble with 2,826 rows and 45 columns:
#' \describe{
#'   \item{\code{cell_type}}{Cell type (e.g., \code{"cd4nc"}, \code{"cd8nc"}, \code{"nk"}, \code{"bin"}, \code{"bmem"}).}
#'   \item{\code{gene_id}}{Ensembl gene ID (e.g., \code{"ENSG00000023902"}).}
#'   \item{\code{snp_id}}{SNP identifier in \code{CHR:POS} format (e.g., \code{"1:150123456"}).}
#'   \item{\code{CHR}}{Chromosome number.}
#'   \item{\code{POS}}{Genomic position (base-pair coordinate).}
#'   \item{\code{SAMP1}, \code{SAMP2}, \dots, \code{SAMP40}}{Genotype dosage for each individual
#'         (typically encoded as 0, 1, or 2).}
#' }
#'
#' @usage data("example_eqtlgeno")
#'
#' @examples
#' data("example_eqtlgeno")
#' example_eqtlgeno
#' dplyr::count(example_eqtlgeno, cell_type)
"example_eqtlgeno"

#' Example B cell single-cell RNA-seq data with pseudotime
#'
#' A \code{SingleCellExperiment} object containing a subset of B cells from the
#' OneK1K cohort with associated individual metadata and a
#' pseudotime trajectory. This dataset is useful for illustrating dynamic
#' eQTL modeling along a continuous trajectory.
#'
#' @format A \code{SingleCellExperiment} object with 817 genes and 3,726 cells.
#' \describe{
#'   \item{assays}{
#'     \itemize{
#'       \item \code{counts}: raw UMI count matrix.
#'       \item \code{norm}: normalized expression values.
#'       \item \code{logcounts}: log-transformed normalized expression values.
#'     }
#'   }
#'   \item{rowData}{A \code{DataFrame} with 2 columns:
#'     \describe{
#'       \item{\code{GeneSymbol}}{Gene symbol.}
#'       \item{\code{features}}{Feature identifier (e.g., Ensembl ID).}
#'     }
#'   }
#'   \item{colData}{A \code{DataFrame} with 5 columns:
#'     \describe{
#'       \item{\code{cell_type}}{Cell-type annotation within the B-cell compartment
#'         (e.g., immature/naive B cells, memory B cells).}
#'       \item{\code{indiv}}{Individual identifier.}
#'       \item{\code{sex}}{Sex of the donor.}
#'       \item{\code{age}}{Age of the donor.}
#'       \item{\code{slingPseudotime_1}}{Inferred pseudotime along a trajectory
#'         (e.g., from immature to memory B cells).}
#'     }
#'   }
#'   \item{reducedDimNames}{\code{"PHATE"}: PHATE embedding used for trajectory
#'         inference and visualization.}
#' }
#'
#' @usage data("example_sce_Bcell")
#'
#' @examples
#' data("example_sce_Bcell")
#' example_sce_Bcell
#' summary(SingleCellExperiment::colData(example_sce_Bcell)$slingPseudotime_1)
"example_sce_Bcell"

#' Example B cell eQTL genotype data
#'
#' A tibble containing example genotype data for selected cis-eQTL SNPs
#' for the selected B cells from OneK1K. This
#' dataset can be used together with \code{example_sce_Bcell} to show
#' dynamic or cell-type–specific eQTL modeling.
#' No eQTL effect size results are included. Genotypes are all permuted.
#'
#' @format A tibble with 2,406 rows and 106 columns:
#' \describe{
#'   \item{\code{cell_type}}{Cell type (e.g., \code{"cd4nc"}, \code{"cd8nc"}, \code{"nk"}, \code{"cd4et"}, etc.).}
#'   \item{\code{gene_name}}{Gene symbol (e.g., \code{"PLEKHO1"}, \code{"TNFRSF1B"}, \code{"UTS2"}).}
#'   \item{\code{gene_id}}{Ensembl gene ID.}
#'   \item{\code{snp_id}}{SNP identifier in \code{CHR:POS} format.}
#'   \item{\code{CHR}}{Chromosome number.}
#'   \item{\code{POS}}{Genomic position (base-pair coordinate).}
#'   \item{\code{SAMP1}, \code{SAMP2}, \dots}{Genotype dosage (0/1/2) for each
#'     individual. The remaining columns \code{SAMPk} store genotypes for all
#'     samples used in the example.}
#' }
#'
#' @usage data("example_eqtlgeno_Bcell")
#'
#' @examples
#' data("example_eqtlgeno_Bcell")
#' example_eqtlgeno_Bcell
#' dplyr::count(example_eqtlgeno_Bcell, gene_name)
"example_eqtlgeno_Bcell"

#' Example list of scDesignPop's marginal models
#'
#' A named list of fitted marginal models for selected gene–SNP
#' pairs used in the examples of scDesignPop. Each element corresponds to one
#' gene and contains the fitted model object, the SNP used as covariate, and
#' additional model attributes.
#'
#' This dataset is useful for demonstrating how to directly perform power
#' analysis using the fitted marginal models of the scDesignPop.
#'
#' @format A named list. Each element corresponds to one gene (e.g.,
#'   \code{"ENSG00000163221"}, \code{"ENSG00000135218"}) and is a list with
#'   the following components:
#' \describe{
#'   \item{\code{fit}}{A fitted generalized linear mixed model, typically
#'     a negative binomial mixed model (e.g., from \pkg{glmmTMB}) with
#'     random intercepts for individuals and fixed effects for cell type,
#'     SNP genotype, and SNP–cell-type interactions.}
#'   \item{\code{time}}{Computation time (in seconds) required to fit the model.}
#'   \item{\code{snp_cov}}{Character string giving the SNP covariate used in
#'     the model (e.g., \code{"1:153337943"}).}
#'   \item{\code{model_attr}}{A list of auxiliary model attributes, including
#'     the \code{terms} object, variable names, and data classes used for
#'     constructing the model matrix.}
#'   \item{\code{removed_cell}}{Information about removed cells (if any)
#'     during model fitting; \code{NA} if no cells were removed.}
#' }
#'

#'
#' @usage data("marginal_list_sel")
#'
#' @examples
#' data("marginal_list_sel")
#' names(marginal_list_sel)
#' str(marginal_list_sel[[1]]$fit)
"marginal_list_sel"

#' Genotype principal components for training individuals
#'
#' Principal component (PC) scores computed from genotype data for a set of
#' training individuals. These PCs are typically used as covariates to account
#' for population structure in the cell type proportion modeling.
#'
#' @format A tibble with 40 rows and 31 columns:
#' \describe{
#'   \item{\code{indiv}}{Individual identifier (e.g., \code{"SAMP1"}, \code{"SAMP2"}).}
#'   \item{\code{PC1}, \code{PC2}, \dots, \code{PC30}}{Principal component scores
#'     summarizing genotype variation across individuals.}
#' }
#'
#' @usage data("example_genopc_train")
#'
#' @examples
#' data("example_genopc_train")
#' head(example_genopc_train)
"example_genopc_train"

#' Genotype principal components for new individuals
#'
#' Principal component (PC) scores computed from genotype data for a set of
#' new individuals. These individuals are projected into the same PC space as
#' the training set to be used in the cell type proportion modeling.
#'
#' @format A tibble with 982 rows and 31 columns:
#' \describe{
#'   \item{\code{indiv}}{Individual identifier (e.g., \code{"NEW_SAMP1"}).}
#'   \item{\code{PC1}, \code{PC2}, \dots, \code{PC30}}{Principal component scores
#'     for new individuals projected into the PC space derived from training individuals.}
#' }
#'
#' @usage data("example_genopc_new")
#'
#' @examples
#' data("example_genopc_new")
#' head(example_genopc_new)
"example_genopc_new"



