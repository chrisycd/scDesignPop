#' Simulate new data
#'
#' \code{simuNewPop} generates new simulated data based on fitted marginal and copula models.
#' This function is adapted from simu_new function in scDesign3 v0.99.7
#'
#'
#' The function takes the new covariate (if use) from \code{\link{constructDataPop}},
#' parameter matrices from \code{\link{extractParaPop}} and multivariate Unifs
#' from \code{\link{fitCopulaPop}}.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assay_use A string which indicates the assay you will use in the sce.
#'                  Default is 'counts'.
#' @param mean_mat A cell by feature matrix of the mean parameter.
#' @param sigma_mat A cell by feature matrix of the sigma parameter.
#' @param zero_mat A cell by feature matrix of the zero-inflation parameter.
#' @param quantile_mat A cell by feature matrix of the multivariate quantile.
#' @param copula_list A list of copulas for generating the multivariate quantile
#'                    matrix. If provided, the \code{quantile_mat} must be NULL.
#' @param n_cores An integer. The number of cores to use.
#' @param fastmvn An logical variable. If TRUE, the sampling of multivariate Gaussian is done
#'                by \code{mvnfast}, otherwise by \code{mvtnorm}. Default is FALSE.
#' @param family_use A string of the marginal distribution. Must be one of
#'                   'poisson', 'nb', or 'gaussian'.
#' @param nonnegative A logical variable. If TRUE, values < 0 in the synthetic
#'                    data will be converted to 0. Default is TRUE (since the expression
#'                    matrix is nonnegative).
#' @param nonzerovar A logical variable. If TRUE, for any gene with zero variance,
#'                   a cell will be replaced with 1. This is designed for avoiding
#'                   potential errors, for example, PCA.
#' @param input_data A input count matrix.
#' @param new_covariate A data.frame which contains covariates of targeted simulated
#'                      data from  \code{\link{constructDataPop}}.
#' @param important_feature important_feature A string or vector which indicates whether a gene will be used in correlation estimation or not. If this is a string, then
#' this string must be either "all" (using all genes) or "auto", which indicates that the genes will be automatically selected based on the proportion of zero expression across cells
#' for each gene. Gene with zero proportion greater than 0.8 will be excluded form gene-gene correlation estimation. If this is a vector, then this should
#' be a logical vector with length equal to the number of genes in \code{sce}. \code{TRUE} in the logical vector means the corresponding gene will be included in
#' gene-gene correlation estimation and \code{FALSE} in the logical vector means the corresponding gene will be excluded from the gene-gene correlation estimation.
#' The default value for is "all".
#' @param parallelization A string indicating the specific parallelization function to use.
#' Must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which corresponds to the parallelization function in the package
#' \code{parallel},\code{BiocParallel}, and \code{pbmcapply} respectively. The default value is 'mcmapply'.
#' @param BPPARAM A \code{MulticoreParam} object or NULL. When the parameter parallelization = 'mcmapply' or 'pbmcmapply',
#' this parameter must be NULL. When the parameter parallelization = 'bpmapply',  this parameter must be one of the
#' \code{MulticoreParam} object offered by the package 'BiocParallel. The default value is NULL.
#' @param filtered_gene A vector or NULL which contains genes that are excluded in the marginal and copula fitting
#' @param mean_limit A numeric scalar to filter genes which has cells that exceed
#'     the limit in the \code{mean_mat}. The default value is 1e15. This is to
#'     avoid features that have extremely high and unreasonable means.
#' @param debug A logical scalar for whether to return a list of variables in
#'     addition to simulated count matrix. The default is FALSE.
#'
#' @return A feature by cell matrix of the new simulated count (expression) matrix or sparse matrix.
#' @export simuNewPop
#'
#' @examples
#' NULL
simuNewPop <- function(sce,
                       assay_use = "counts",
                       mean_mat,
                       sigma_mat,
                       zero_mat,
                       quantile_mat = NULL,
                       copula_list,
                       n_cores,
                       fastmvn = FALSE,
                       family_use,
                       nonnegative = TRUE,
                       nonzerovar = FALSE,
                       input_data,
                       new_covariate,
                       important_feature = "all",
                       parallelization = "mcmapply",
                       BPPARAM = NULL,
                       filtered_gene,
                       mean_limit = 1e15,
                       debug = FALSE) {

    stopifnot("Features in sce and mean_mat do not match. Please check input!" !=
                  checkVectorEqual(rownames(sce), colnames(mean_mat), ignore_order = FALSE))

    if(!is.null(quantile_mat) & !is.null(copula_list)) {
        stop("You can only provide either the quantile_mat or the copula_list!")
    }

    # check if user inputted new covariates
    data_temp <- input_data[, colnames(new_covariate), drop = FALSE]

    if(identical(data_temp, new_covariate)){
        new_covariate <- NULL
    }

    qc_gene_idx <- which(!rownames(sce) %in% filtered_gene)

    # filter genes exceeding mean value
    mean_lim_idx <- which(apply(mean_mat, 2, function(col_vec) any(col_vec >= mean_limit)))
    qc_gene_idx <- base::setdiff(qc_gene_idx, mean_lim_idx)

    if(length(family_use) != 1){
        family_use <- family_use[qc_gene_idx]
    }

    if(!is.null(quantile_mat)) {
        message("Multivariate quantile matrix is provided")

        newmvn_full <- matrix(NA, nrow = dim(quantile_mat)[1], ncol = dim(sce)[1])
        rownames(newmvn_full) <- rownames(quantile_mat)
        colnames(newmvn_full) <- rownames(sce)
        newmvn_full[rownames(quantile_mat), colnames(quantile_mat)] <- quantile_mat
        quantile_mat <- as.matrix(newmvn_full)

    } else {
        message("Use Copula to sample a multivariate quantile matrix")

        group_index <- unique(input_data$corr_group)
        corr_group <- as.data.frame(input_data$corr_group)
        colnames(corr_group) <- "corr_group"
        ngene <- length(qc_gene_idx)

        if(is.null(new_covariate)) {
            new_corr_group <- NULL
        } else {
            new_corr_group <- as.data.frame(new_covariate$corr_group)
            colnames(new_corr_group) <- "corr_group"
        }

        ind <- group_index[1] == "ind"

        newmvn.list <-
            lapply(group_index, function(x,
                                         sce,
                                         corr_group,
                                         new_corr_group,
                                         ind,
                                         n_cores,
                                         copula_list) {
                message(paste0("Sample Copula group ", x, " starts"))
                curr_index <- which(corr_group[, 1] == x)

                if (is.null(new_covariate)) {
                    curr_ncell <- length(curr_index)
                    curr_ncell_idx <- curr_index
                } else {
                    curr_ncell <- length(which(new_corr_group[, 1] == x))
                    curr_ncell_idx <-which(new_corr_group[, 1] == x)
                    #paste0("Cell", which(new_corr_group[, 1] == x))
                }
                cor.mat <- copula_list[[x]]

                if(curr_ncell == 0) {
                    new_mvu <- NULL
                } else {
                    if (methods::is(cor.mat, "matrix")) {
                        #message(paste0("Group ", group_index, " Start"))

                        #message("Sample MVN")
                        #sample from mvn for important genes only
                        corr_gene_idx <- apply(cor.mat, 2, function(x) length(which(x < 1e-5)) != length(x)-1)
                        corr_gene <- colnames(cor.mat)[which(corr_gene_idx)]
                        new_mvn_important <- sampleMVN(n = curr_ncell,
                                                       Sigma = cor.mat[corr_gene, corr_gene],
                                                       n_cores = n_cores,
                                                       fastmvn = fastmvn)
                        colnames(new_mvn_important) <- corr_gene
                        #message("MVN Sampling End")
                        ind_gene <- colnames(cor.mat)[which(corr_gene_idx == FALSE)]

                        if(length(ind_gene) > 0) {
                            new_mvn_non_important <- lapply(ind_gene, function(x) return(stats::rnorm(n = curr_ncell)))
                            new_mvn_non_important_mat <- do.call("cbind",new_mvn_non_important)
                            colnames(new_mvn_non_important_mat) <- ind_gene

                            mvnrvq <- apply(new_mvn_non_important_mat, 2, stats::pnorm)
                            new_mvu <- cbind(new_mvn_important, mvnrvq)
                            new_mvu <- new_mvu[,colnames(cor.mat)]
                        } else {
                            new_mvu <- new_mvn_important
                        }

                        rownames(new_mvu) <- curr_ncell_idx
                    } else if (methods::is(cor.mat, "vinecop")) {
                        new_mvu <- matrix(0, nrow = curr_ncell, ncol = ngene)
                        #message("Sampling Vine Copula Starts")
                        mvu <- rvinecopulib::rvinecop(
                            curr_ncell,
                            vine = cor.mat,
                            cores = n_cores,
                            qrng = TRUE
                        )

                        new_mvu[, which(important_feature)] <- mvu
                        if(length(which(important_feature)) != ngene){
                            cor.mat <- diag(rep(1, length(which(!important_feature))))
                            mvu2 <- sampleMVN(n = curr_ncell,
                                              Sigma = cor.mat,
                                              n_cores = n_cores,
                                              fastmvn = fastmvn)
                            new_mvu[, which(!important_feature)] <- mvu2
                        }
                        #message("Sampling Vine Copula Ends")
                        rownames(new_mvu) <- curr_ncell_idx
                    } else if(ind) {
                        "Use independent copula (random Unif)."
                        new_mvu <-
                            matrix(data = stats::runif(curr_ncell * ngene),
                                   nrow = curr_ncell)
                        rownames(new_mvu) <- curr_ncell_idx
                    } else{
                        stop("Copula must be one from 'vine' or 'gaussian', or assume gene-gene is independent")
                    }
                }

                return(list(new_mvu = new_mvu))
                }, sce = sce,
                   ind = ind,
                   n_cores = n_cores,
                   corr_group = corr_group,
                   new_corr_group = new_corr_group,
                   copula_list = copula_list)

        newmvn <-
            do.call(rbind, lapply(newmvn.list, function(x) {
                x$new_mvu[, qc_gene_idx]  # use filtered genes
                }))

        newmvn[as.numeric(rownames(newmvn)), ] <- newmvn
        rownames(newmvn) <- as.character(1:dim(newmvn)[1])
        colnames(newmvn) <- rownames(sce)[qc_gene_idx]

        newmvn_full <- matrix(NA, nrow = dim(newmvn)[1], ncol = dim(sce)[1])

        rownames(newmvn_full) <- rownames(newmvn)
        colnames(newmvn_full) <- rownames(sce)
        newmvn_full[rownames(newmvn), colnames(newmvn)] <- newmvn

        quantile_mat <- as.matrix(newmvn_full)
    }

    mat_function <- function(x, y) {

        idx <- which(mean_mat[, x] != 0)
        para_mat <- cbind(mean_mat[idx, x], sigma_mat[idx, x], quantile_mat[idx, x], zero_mat[idx, x])

        if (y == "binomial") {
            qfvec <- stats::qbinom(p = para_mat[, 3], prob = para_mat[, 1], size = 1)
        } else if (y == "poisson") {

            qfvec <- stats::qpois(p = para_mat[, 3], lambda = para_mat[, 1])
        } else if (y == "gaussian") {
            qfvec <-
                gamlss.dist::qNO(p = para_mat[, 3],
                                 mu = para_mat[, 1],
                                 sigma = abs(para_mat[, 2]))
        } else if (y == "nb") {

            ## avoids excess qNBI compute time when NB mean exceeds 1e8
            mu_scaled <- para_mat[, 1]
            is_scaled <- mu_scaled > 1e+8  # logic vec for scaling condition
            k <- ifelse(is_scaled, mu_scaled / 1e+4, 1)  # dynamic k based on ratio

            # apply scaling
            mu_scaled <- mu_scaled / k

            qfvec <-
                gamlss.dist::qNBI(p = para_mat[, 3],
                                  mu = mu_scaled,  # use modified mu vec
                                  sigma = para_mat[, 2])

            # scale quantiles back to original
            qfvec[is_scaled] <- round(qfvec[is_scaled] * k[is_scaled])
        } else if (y == "zip") {
            qfvec <-
                gamlss.dist::qZIP(p = para_mat[, 3],
                                  mu = para_mat[, 1],
                                  sigma = ifelse(para_mat[, 4] != 0, para_mat[, 4],  2.2e-16))  ## Avoid zero zero-inflated prob
        } else if (y == "zinb") {

            qfvec <-
                gamlss.dist::qZINBI(p = para_mat[, 3],
                                    mu = para_mat[, 1],
                                    sigma = para_mat[, 2],
                                    nu = ifelse(para_mat[, 4] != 0, para_mat[, 4],  2.2e-16))
        } else {
            stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
        }

        #message(paste0("Gene ", x , " End!"))

        r <- as.vector(qfvec)
        if(length(r) < total_cells){
            new_r <- rep(0, total_cells)
            new_r[idx] <- r
            names(new_r) <- cell_names
            r <- new_r
        }

        r
    }


    ## New count
    paraFunc <- parallel::mcmapply

    if(parallelization == "bpmapply"){
        paraFunc <- BiocParallel::bpmapply
    }
    if(parallelization == "pbmcmapply"){
        paraFunc <- pbmcapply::pbmcmapply
    }

    if(is.null(new_covariate)){
        total_cells <- dim(sce)[2]
        cell_names <- colnames(sce)
    }else{
        total_cells <- dim(new_covariate)[1]
        cell_names <- rownames(new_covariate)
    }

    if(parallelization == "bpmapply"){
        BPPARAM$workers <- n_cores
        mat <-  paraFunc(mat_function, x = seq_len(dim(sce)[1])[qc_gene_idx], y = family_use, SIMPLIFY = TRUE, BPPARAM = BPPARAM)
    }else{
        mat <- paraFunc(mat_function, x = seq_len(dim(sce)[1])[qc_gene_idx], y = family_use, SIMPLIFY = TRUE
                        , mc.cores = n_cores
        )
    }
    new_count <- mat #simplify2array(mat)
    rownames(new_count) <- cell_names
    colnames(new_count) <- rownames(sce)[qc_gene_idx]

    if(length(qc_gene_idx) < dim(sce)[1]){
        temp_count <- matrix(0, total_cells, dim(sce)[1])
        rownames(temp_count) <- cell_names
        colnames(temp_count) <- rownames(sce)
        temp_count[rownames(new_count),colnames(new_count)] <- new_count
        new_count <- temp_count
    }
    new_count <- as.matrix(t(new_count))

    if(nonnegative) new_count[new_count < 0] <- 0

    if(nonzerovar) {
        row_vars <- matrixStats::rowVars(new_count[qc_gene_idx,])
        if(sum(row_vars == 0) > 0) {
            message("Some genes have zero variance. Replace a random one with 1.")
            row_vars_index <- which(row_vars == 0)
            col_index <- seq_len(dim(new_count)[2])
            for(i in row_vars_index) {
                new_count[i, sample(col_index, 1)] <- 1
            }
        }
    }

    if(methods::is(SummarizedExperiment::assay(sce, assay_use), "sparseMatrix")){
        new_count<- Matrix::Matrix(new_count, sparse = TRUE)
    }

    if(debug) {
        return(list("new_count" = new_count,
                    "qc_gene_idx" = qc_gene_idx,
                    "mean_lim_idx" = mean_lim_idx))
    } else {
        return(new_count)
    }
}


### Sample MVN based on cor
sampleMVN <- function(n,
                      Sigma,
                      n_cores = n_cores,
                      fastmvn = fastmvn) {
    if(fastmvn) {
        mvnrv <- mvnfast::rmvn(n, mu = rep(0, dim(Sigma)[1]), sigma = Sigma, ncores = n_cores)
    } else {
        mvnrv <-
            rmvnorm(n, mean = rep(0, dim(Sigma)[1]), sigma = Sigma, checkSymmetry = FALSE, method="eigen")
    }
    mvnrvq <- apply(mvnrv, 2, stats::pnorm)

    return(mvnrvq)
}


## fix integer overflow issue when # of gene* # of cell is too larger in rnorm(n * ncol(sigma))
## This function comes from package Mvnorm.
rmvnorm <- function(n,
                    mean = rep(0, nrow(sigma)),
                    sigma = diag(length(mean)),
                    method = c("eigen", "svd", "chol"),
                    pre0.9_9994 = FALSE,
                    checkSymmetry = TRUE) {

    if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                                      check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma))
        stop("mean and sigma have non-conforming size")

    method <- match.arg(method)

    R <- if(method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))){
            warning("sigma is numerically not positive semidefinite")
        }
        ## ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
        ## faster for large  nrow(sigma):
        t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
    }
    else if(method == "svd"){
        s. <- svd(sigma)
        if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
            warning("sigma is numerically not positive semidefinite")
        }
        t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
    }
    else if(method == "chol"){
        R <- chol(sigma, pivot = TRUE)
        R[, order(attr(R, "pivot"))]
    }

    retval <- matrix(stats::rnorm(as.double(n) * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
}
