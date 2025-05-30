% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructDataPop.R
\name{constructDataPop}
\alias{constructDataPop}
\title{Construct a list of input data}
\usage{
constructDataPop(
  sce,
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
  prune_thres = 0.9
)
}
\arguments{
\item{sce}{a SingleCellExperiment object.}

\item{eqtlgeno_df}{a dataframe with eQTL annotations and SNP genotypes for each gene.}

\item{new_eqtlgeno_df}{a dataframe with eQTL annotations and SNP genotypes for each
gene in new individuals. The default is NULL.}

\item{new_covariate}{a cell covariate dataframe for which to simulate.
The default is NULL.}

\item{overlap_features}{an optional string vector to filter for features (ie. genes).
The default is NULL.}

\item{sampid_vec}{an optional string vector to filter for sample ids.
The default is NULL.}

\item{ct_copula}{a logical scalar for whether to fit the copula by cell state
variable specified by \code{cellstate_colname} option. The default is TRUE.}

\item{slot_name}{a string scalar specifying the slot to use in input \code{sce}.
The default is "counts".}

\item{snp_model}{a string scalar specifying the type of SNP model used. Options
are either "single" for single-SNP, or "multi" for multi-SNP.}

\item{cellstate_colname}{a string scalar specifying the cell state variable in
\code{eqtlgeno_df} and cell covariate of \code{sce} object.
The default is "cell_type".}

\item{feature_colname}{a string scalar specifying the feature variable (ie. genes)
in \code{eqtlgeno_df}. The default is "gene_id".}

\item{snp_colname}{a string scalar for the SNP variable in \code{eqtlgeno_df}.
The default is "snp_id".}

\item{loc_colname}{a string scalar for the last column of eQTL annotation in
\code{eqtlgeno_df}. The default is "POS".}

\item{chrom_colname}{a string scalar of the chromosome variable in \code{eqtlgeno_df}.
The default is "CHR".}

\item{indiv_colname}{a string scalar of the sample ID variable in cell covariate
of \code{sce}. The default is "indiv".}

\item{prune_thres}{a numerical value between 0 and 1 used to threshold the pairwise
correlations of eQTLs' genotypes for each feature. The default value is 0.9.}
}
\value{
outputs a list with following elements:
\describe{
\item{\code{count_mat}}{a cell-by-gene matrix of response values.}
\item{\code{covariate}}{a cell-by-covariate data frame used for fit marginal.}
\item{\code{new_covariate}}{an optional cell-by-covariate data frame used for prediction.}
\item{\code{important_features}}{a string vector of gene ids.}
\item{\code{eqtl_geno_list}}{a list of eQTL genotype dataframes for each gene
(for fit marginal).}
\item{\code{new_eqtl_geno_list}}{a optional list of eQTL genotype dataframes
for each gene (for new individual simulation).}
\item{\code{filtered_gene}}{string vector of features QC filtered.}
}
}
\description{
Function extracts an expression matrix, cell covariates, and filters the SNPs
in the eQTL genotype dataframe.
}
\examples{
NULL
}
