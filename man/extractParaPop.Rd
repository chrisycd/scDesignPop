% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractParaPop.R
\name{extractParaPop}
\alias{extractParaPop}
\title{Extract parameter matrix for new covariate df}
\usage{
extractParaPop(
  sce,
  assay_use = "counts",
  marginal_list,
  n_cores,
  family_use,
  new_covariate,
  new_eqtl_geno_list,
  indiv_colname = "indiv",
  snp_colname = "snp_id",
  loc_colname = "POS",
  parallelization = "mcmapply",
  BPPARAM = NULL,
  data
)
}
\arguments{
\item{sce}{add later}

\item{assay_use}{add later}

\item{marginal_list}{add later}

\item{n_cores}{add later}

\item{family_use}{a string scalar or vector of marginal distribution used.}

\item{new_covariate}{a cell-by-feature covariate dataframe (from construct_data.R) plus corr_group.}

\item{new_eqtl_geno_list}{a list of eQTL genotype dataframes for each gene (to be predicted).}

\item{indiv_colname}{add later}

\item{snp_colname}{add later}

\item{loc_colname}{add later}

\item{parallelization}{add later}

\item{BPPARAM}{add later}

\item{data}{a cell-by-feature covariate dataframe (from construct_data.R) plus corr_group. Used only in gamlss fits.}
}
\value{
a list of mean, sigma, and zero parameter cell by feature matrices:
}
\description{
This is the main function.
}
\examples{
NULL
}
