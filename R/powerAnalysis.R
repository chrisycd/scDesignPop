#' Construct a model formula for power analysis
#'
#' @param fm a stats::formula object from the full marginal model.
#' @param method a character object specifying the method that will be analyzed for power.
#' (Options: nb,poisson,gaussian,pseudoBulkLinear).
#' @param snpid a character object contains snpid.
#' @param type_specific a character object contains the name of the covariate that the analysis is specific on.
#'
#' @return a stats::formula object for power analysis in each specified type.
#' @export
#'
#' @examples
#' NULL
constructPAFormula <- function(fm,
                               method=c("nb","poisson","gaussian","pseudoBulkLinear"),
                               snpid=NULL,
                               type_specific=NULL){
  # checks
  method = match.arg(method)
  stopifnot("snpid is not included in the formula. Please check input!" =
                (stringr::str_detect(as.character(fm)[3],snpid)))
  stopifnot("type_specific is not included in the formula. Please check input!" =
                (stringr::str_detect(as.character(fm)[3],type_specific)))

  # model formula
  fm <- as.character(fm)
  # separate all terms from the rightside of the formula
  fm_rs <- fm[3]
  fm_rspart <- unlist(stringr::str_split(fm_rs,pattern = " \\+ "))
  # remove all effects related to the type_specific parameter
  fm_rspart1 <- fm_rspart[!stringr::str_detect(fm_rspart,type_specific)]
  # remove all genotype effects
  fm_rspart2 <- fm_rspart1[!stringr::str_detect(fm_rspart1,"`")]
  # add specified snp's genotype effect
  fm_rspart3 <- c(fm_rspart2,paste0("`",snpid,"`"))

  if(method=="nb"||method=="poisson"||method=="gaussian"){
    fm_rs_re <- paste(fm_rspart3,collapse = " + ")
    fm[3] <- fm_rs_re
    fm <- paste(fm[2],fm[1],fm[3]," ")
  }
  else if(method=="pseudoBulkLinear"){
    # remove individual random effect
    fm_rspart4 <- fm_rspart3[!stringr::str_detect(fm_rspart3,"\\|")]
    fm_rs_re <- paste(fm_rspart4,collapse = " + ")
    fm[3] <- fm_rs_re
    fm <- paste(fm[2],fm[1],fm[3]," ")
  }

  # construct the formula
  model_formula <- stats::as.formula(fm)
  return(model_formula)
}

#' Fit a marginal model for power analysis
#'
#' @param df a data frame object contains the design matrix.
#' @param model_formula a stats::formula object contains the model formula for power analysis.
#' @param idx a numeric value recording the serial number of the simulation
#' @param method a character object specifying the method that will be analyzed for power.
#' (Options: nb,poisson,gaussian,pseudoBulkLinear).
#' @param snpid a character object contains snpid.
#'
#' @return a fitted stats::model object.
#' @export
#'
#' @examples
#' NULL
fitPAModel <- function(df,
                       model_formula,
                       idx,
                       method=c("nb","poisson","gaussian","pseudoBulkLinear"),
                       snpid) {
  # checks
  method = match.arg(method)
  stopifnot("snpid is not included in the formula. Please check input!" =
                  (stringr::str_detect(as.character(model_formula)[3],snpid)))

  if (method=="nb"){
    result <- tryCatch({

      expr = {
        withCallingHandlers(
          expr = {
            model <- glmmTMB::glmmTMB(formula = model_formula,
                                      data = df,
                                      family = glmmTMB::nbinom2,
                                      ziformula = ~0)
            # Note: use ziformula = ~1 for zero-inflation
            result <- glmmTMB::fixef(model)$cond[paste0("`",snpid,"`")]
          },
          # If expression throws a warning, record diagnostics without halting,
          # so as to store the result of the expression.
          warning = function(w){

            message(sprintf("glmmTMB fit issues on i=%s: %s \n", idx, base::conditionMessage(w)))

            # parent <- parent.env(environment())
            # parent$diag <- w

          }
        )
      }

    }, error = function(e) {

      message(sprintf("glmmTMB fit fails on i=%s: %s \n", idx, base::conditionMessage(e)))

      NA

    }, silent = FALSE)

  }else if (method=="poisson"){
      result <- tryCatch({

          expr = {
              withCallingHandlers(
                  expr = {
                      model <- glmmTMB::glmmTMB(formula = model_formula,
                                                data = df,
                                                family = stats::poisson())
                      # Note: use ziformula = ~1 for zero-inflation
                      result <- glmmTMB::fixef(model)$cond[paste0("`",snpid,"`")]
                  },
                  # If expression throws a warning, record diagnostics without halting,
                  # so as to store the result of the expression.
                  warning = function(w){

                      message(sprintf("glmmTMB fit issues on i=%s: %s \n", idx, base::conditionMessage(w)))

                      # parent <- parent.env(environment())
                      # parent$diag <- w

                  }
              )
          }

      }, error = function(e) {

          message(sprintf("glmmTMB fit fails on i=%s: %s \n", idx, base::conditionMessage(e)))

          NA

      }, silent = FALSE)

  }else if (method=="gaussian"){
    # # generate pseudo-bulk
    # df <- df %>% dplyr::group_by(indiv) %>%
    #   dplyr::mutate(response=sum(response)) %>%
    #   dplyr::distinct()
    #colnames(df)[length(colnames(df))] <- "genotype"

    result <- tryCatch({

      expr = {
        withCallingHandlers(
          expr = {
            # model <- nlme::lme(fixed = response ~ genotype,random = ~1|indiv,data = df,
            #                    control = nlme::lmeControl(opt = "optim"))
            # result <- summary(model)$coefficients$fixed[2]

            model <- glmmTMB::glmmTMB(formula = model_formula,
                                      data = df,
                                      family = "gaussian")
            result <- glmmTMB::fixef(model)$cond[paste0("`",snpid,"`")]
          },
          # If expression throws a warning, record diagnostics without halting,
          # so as to store the result of the expression.
          warning = function(w){

            message(sprintf("glmmTMB fit issues on i=%s: %s \n", idx, base::conditionMessage(w)))

            # parent <- parent.env(environment())
            # parent$diag <- w

          }
        )
      }

    }, error = function(e) {

      message(sprintf("glmmTMB fit fails on i=%s: %s \n", idx, base::conditionMessage(e)))

      NA

    }, silent = FALSE)

  }else if (method=="pseudoBulkLinear"){
      # generate pseudo-bulk
      df <- df %>% dplyr::group_by(!!rlang::sym("indiv")) %>%
          dplyr::mutate("response"=sum(!!rlang::sym("response"))) %>%
          dplyr::distinct()

      result <- tryCatch({

          expr = {
              withCallingHandlers(
                  expr = {
                      model <- stats::lm(formula = model_formula,data = df)
                      result <- summary(model)$coefficients[paste0("`",snpid,"`"),1]
                  },
                  # If expression throws a warning, record diagnostics without halting,
                  # so as to store the result of the expression.
                  warning = function(w){

                      message(sprintf("lm fit issues on i=%s: %s \n", idx, base::conditionMessage(w)))

                      # parent <- parent.env(environment())
                      # parent$diag <- w

                  }
              )
          }

      }, error = function(e) {

          message(sprintf("lm fit fails on i=%s: %s \n", idx, base::conditionMessage(e)))

          NA

      }, silent = FALSE)

  }



  return(result)
}


#' Simulate new design matrix for power analysis
#'
#' @param fit a fitted stats::model object.
#' @param df_sel a data frame contains the new design matrix.
#' @param nindiv_total a vector contains the number of individuals for each genotype.
#' @param model a character showing the model types of the full marginal model.
#' @param snpid a character object contains snpid.
#' @param nindiv a numeric value showing the number of individuals that user wants to simulate.
#' @param ncell a numeric value showing the number of cells per each individual that user wants to simulate.
#' @param type_specific a character object contains the name of the covariate that the analysis is specific on.
#'
#' @return a new data frame contains the design matrix with simulated response.
#' @export
#'
#' @examples
#' NULL
simulatePADesignMatrix <- function(fit,
                                   df_sel,
                                   nindiv_total,
                                   model=c("nb","poisson","gaussian"),
                                   snpid,
                                   nindiv,
                                   ncell,
                                   type_specific){
    # checks
    stopifnot("nindiv_total doesn't contain three numbers of individuals for the three genotypes. Please check input!" =
                  (length(nindiv_total) == 3))
    model = match.arg(model)
    stopifnot("snpid is not included in the df_sel data.frame. Please check input!" =
                  (snpid%in%colnames(df_sel)))
    stopifnot("nindiv is not a integer. Please check input!" =
                  (abs(nindiv - round(nindiv)) < .Machine$double.eps^0.5))
    stopifnot("ncell is not a integer. Please check input!" =
                  (abs(ncell - round(ncell)) < .Machine$double.eps^0.5))

    # extract genotypes
    geno <- df_sel[,snpid]
    names(geno) <- df_sel$indiv
    geno <- geno[unique(df_sel$indiv)]
    geno0 <- geno[which(geno==0)]
    geno1 <- geno[which(geno==1)]
    geno2 <- geno[which(geno==2)]

    # sampling until all combination presented
    # avoid no points in a subgroup
    flag1=F
    while(flag1==F){
        rand_indiv <- c(sample(names(geno0),nindiv_total[1],replace = T),
                        sample(names(geno1),nindiv_total[2],replace = T),
                        sample(names(geno2),nindiv_total[3],replace = T))
        df_tmp <- df_sel[which(df_sel$indiv%in%rand_indiv),]
        df_tmp <- df_tmp[,setdiff(colnames(df_tmp),c("response","indiv"))]
        if(length(unique(df_tmp))>=length(colnames(df_tmp))){
            flag1=T
        }
    }

    # avoid no enough different response values
    flag2=F
    while(flag2==F){
        rand_cellindex <- lapply(rand_indiv,function(indiv){
            return(sample(x = rep(which(df_sel$indiv==indiv),2),size = ncell,replace = T))
        }) # in case only one cells, causing program truffles

        # construct new covariates data with new individuals
        df_news <- do.call(rbind,lapply(1:length(rand_cellindex),function(n){
            cellindex <- rand_cellindex[[n]]
            df_new <- df_sel[cellindex,]
            df_new$indiv <- paste0("NewIndiv",n)
            return(df_new)
        }))
        if(model=="nb" || model=="poisson"){
            flag2=T
        }else if(model=="gaussian"){
            if(length(unique(df_news$response))>=4){
                flag2=T
            }
        }
    }

    df_news$indiv=factor(df_news$indiv,
                         levels=paste0("NewIndiv",1:length(rand_cellindex)))

    # use the following line to generate mean estimates without random effect
    colnames(df_news)[colnames(df_news)=="type"] <- type_specific
    response_new <- stats::predict(fit,type = "response",newdata=df_news,
                                   allow.new.levels = TRUE)
    # manually add random effect
    # simulate
    rand_sd <- sqrt(as.numeric(summary(fit)$varcor$cond$indiv))
    newindiv_effect <- stats::rnorm(nindiv,mean = 0,sd = rand_sd)

    if(model=="gaussian"){
        response_new <- response_new+rep(newindiv_effect,each=ncell)
        df_news$response <- stats::rnorm(n=length(response_new),mean=response_new,
                                  sd = glmmTMB::sigma(fit))
    }else if(model=="nb"){
        response_new <- exp(log(response_new)+rep(newindiv_effect,each=ncell))
        df_news$response <- stats::rnbinom(n=length(response_new),mu=response_new,
                                  size = glmmTMB::sigma(fit))
    }else if(model=="poisson"){
        response_new <- exp(log(response_new)+rep(newindiv_effect,each=ncell))
        df_news$response <- stats::rpois(n=length(response_new),lambda=response_new)
    }

    return(df_news)
}


#' Perform a power analysis
#'
#' @param marginal_list the output of function fitMarginalPop().
#' @param marginal_model a character showing the model types of the full marginal model.
#' @param refit_formula the formula used to refit the marginal full model if user wants to. Default is null.
#' @param geneid a character object contains geneid.
#' @param snpid a character object contains snpid.
#' @param type_specific a character object contains the name of the covariate that the analysis is specific on.
#' @param type_vector a vector object contains the specified type/level names of the covariate.
#' @param method a character object specifying the method that will be analyzed for power.
#' (Options: nb,poisson,gaussian,pseudoBulkLinear).
#' @param nindivs a vector of numeric values showing the numbers of individuals that user wants to simulate.
#' @param ncells a vector of numeric values showing the numbers of cells per each individual that user wants to simulate.
#' @param nPool a vector of numeric values showing how many pools of sequencing has been performed.
#' @param nIndivPerPool a numerical value showing how many individuals are sequenced in one pool.
#' @param nCellPerPool a vector of numeric values showing how many cells are sequenced in one pool.
#' @param alpha the p value threshold for rejecting the H0 hypothesis.
#' @param nsims number of simulations for calculating the power. This parameter will affect the resolution of the power value.
#' @param ncores number of CPU cores user wants to use.
#'
#' @return a list of named features, each containing a list with the following items:
#' \describe{
#'      \item{\code{intercept}}{the intercept for the genotype effect in the specified type/level.}
#'      \item{\code{slope}}{the slope for the genotype effect in the specified type/level.}
#'      \item{\code{power}}{a data frame contains power values in different parameter settings.}
#'      \item{\code{data}}{a data frame contains both the H1 and H0 genotype effect estimates in different parameter settings and simulation times.}
#' }
#' @export
#'
#' @examples
#' NULL
powerAnalysis <- function(marginal_list,
                          marginal_model,
                          refit_formula = NULL,
                          geneid,
                          snpid,
                          type_specific,
                          type_vector,
                          method = c("nb","poisson","gaussian","pseudoBulkLinear"),
                          nindivs = NULL,
                          ncells = NULL,
                          nPool = NULL,
                          nIndivPerPool = NULL,
                          nCellPerPool = NULL,
                          alpha=0.05,
                          nsims=100,
                          ncores=1){
  # checks
  stopifnot("marginal_model is not available. Please check input!" =
                 (marginal_model %in% c("nb","poisson","gaussian")))
  stopifnot("geneid is not included in the marginal_list. Please check input!" =
                (geneid%in%names(marginal_list)))
  stopifnot("snpid is not included in the selected gene's marginal model. Please check input!" =
                 (snpid%in%colnames(marginal_list[[geneid]]$fit$frame)))
  stopifnot("type_specific is not included in the selected gene's marginal model. Please check input!" =
                (type_specific%in%colnames(marginal_list[[geneid]]$fit$frame)))
  stopifnot("type_vector is not included in the selected gene's marginal model. Please check input!" =
                (type_vector%in%unique(marginal_list[[geneid]]$fit$frame[,type_specific])))
  method = match.arg(method)
  stopifnot("alpha input is not between 0 and 1. Please check input!" =
                (alpha >=0 && alpha <=1))
  stopifnot("nsims input is not a integer. Please check input!" =
                (abs(nsims - round(nsims)) < .Machine$double.eps^0.5))
  stopifnot("ncores input is not a integer. Please check input!" =
                (abs(ncores - round(ncores)) < .Machine$double.eps^0.5))

  # # set seed
  # set.seed(111)

  # count combinations
  if(is.null(nindivs) && is.null(ncells)){
      # checks when nindivs and ncells are not used
      stopifnot("nPool input is not a vector of integer values. Please check input!" =
                    (mean(abs(nPool - round(nPool)) < .Machine$double.eps^0.5)==1))
      stopifnot("nIndivPerPool input is not a vector of integer values. Please check input!" =
                    (mean(abs(nIndivPerPool - round(nIndivPerPool)) < .Machine$double.eps^0.5)==1))
      stopifnot("nCellPerPool input is not a vector of integer values. Please check input!" =
                    (mean(abs(nCellPerPool - round(nCellPerPool)) < .Machine$double.eps^0.5)==1))

      nindivs <- merge(nPool,nIndivPerPool)
      nindivs$z <- nindivs$x * nindivs$y
      nindivs <- sort(unique(nindivs$z))

      ncells <- merge(nCellPerPool,nIndivPerPool)
      ncells$z <- ncells$x / ncells$y
      ncells <- sort(unique(ncells$z))

  }else{
      stopifnot("nindivs input is not a vector of integer values. Please check input!" =
                    (mean(abs(nindivs - round(nindivs)) < .Machine$double.eps^0.5)==1))
      stopifnot("ncells input is not a vector of integer values. Please check input!" =
                    (mean(abs(ncells - round(ncells)) < .Machine$double.eps^0.5)==1))
  }

  message(paste("Number of individuals:",paste(nindivs,collapse = ",")))
  message(paste("Number of cells per individual:",paste(ncells,collapse = ",")))
  nindiv_cell_dat <- merge(nindivs,ncells)

  # select gene
  #df <- df_construct(data_list,geneid)
  df <- marginal_list[[geneid]]$fit$frame

  # user-specified full model refitting
  if(!is.null(refit_formula)){
    message("Refitting the whole marginal data with input geneid...")
    refit_formula <- stats::as.formula(refit_formula)
    if(marginal_model=="nb"){
        fit <- glmmTMB::glmmTMB(formula = refit_formula,
                                data = df,
                                family = glmmTMB::nbinom2,
                                ziformula = ~0)
    }else if(marginal_model=="poisson"){
        fit <- glmmTMB::glmmTMB(formula = refit_formula,
                                data = df,
                                family = stats::poisson())
    }else if(marginal_model=="gaussian"){
        fit <- glmmTMB::glmmTMB(formula = refit_formula,
                                data = df,
                                family = "gaussian")
    }

    # build the formula using snpid
    model_formula <- constructPAFormula(fm=refit_formula,method = method,snpid = snpid,type_specific = type_specific)

  }else{
    fit <- marginal_list[[geneid]]$fit
    fm <- marginal_list[[geneid]][["fit"]][["call"]][["formula"]]

    # build the formula using snpid
    model_formula <- constructPAFormula(fm=fm,method = method,snpid = snpid,type_specific = type_specific)
  }

  # null model
  fit0 <- fit
  snp_index <- which(stringr::str_detect(names(glmmTMB::fixef(fit)$cond),snpid))

  fit0$fit$par[snp_index] <- 0 # change the corresponding sizes to zeros
  fit0$fit$parfull[snp_index] <- 0
  fit0$sdr$par.fixed[snp_index] <- 0

  # main
  res <- list()
  for(type in type_vector){

    #TODO: construct the PAModel based on user specified input

    message(paste0(type_specific,": ",type))

    # collect intercept and slope of genotype effects for celltypes
    # TODO: for other type_specific order
    ct_index <- 1:length(type_vector)
    # d_index <- ct_index[length(ct_index)]+1
    # a_index <- d_index[length(d_index)]+1
    snp_cov <-glmmTMB::fixef(fit)$cond[snp_index]

    colnames(df)[colnames(df)==type_specific] <- "type"
    if(which(levels(df$type)==type)==1){
      intercept <- as.numeric(glmmTMB::fixef(fit)$cond[1])
    }else{
      intercept <- sum(glmmTMB::fixef(fit)$cond[c(1,which(levels(df$type)==type))])
    }
    if(which(levels(df$type)==type)==1){
      slope <- as.numeric(snp_cov[1])
    }else if(mean(stringr::str_detect(names(snp_cov),paste0(type_specific,type))) > 0){
      slope <- sum(snp_cov[1],snp_cov[which(stringr::str_detect(names(snp_cov),paste0(type_specific,type)))])
    }else{
      slope <- as.numeric(snp_cov[1])
    }
    print(intercept)
    print(slope)

    # select type
    df_sel <- df[which(df$type==type),]

    output <- list()
    powers <- list()
    for(i in 1:length(nindiv_cell_dat[,1])){
      # fetch nindiv and ncell
      nindiv <- nindiv_cell_dat[i,1]
      ncell <- nindiv_cell_dat[i,2]
      message(paste("Loop:",i," nindiv:",nindiv," ncell:",ncell))

      # power analysis
      # control MAFs
      geno <- df_sel[,snpid]
      names(geno) <- df_sel$indiv
      geno <- geno[unique(df_sel$indiv)]
      geno0 <- geno[which(geno==0)]
      geno1 <- geno[which(geno==1)]
      geno2 <- geno[which(geno==2)]
      l_tot <- length(geno)
      p0 <- length(geno0)/l_tot
      p1 <- length(geno1)/l_tot
      p2 <- length(geno2)/l_tot
      nindiv0 <- round(nindiv * p0)
      nindiv1 <- round(nindiv * p1)
      nindiv2 <- round(nindiv * p2)
      nindiv_total <- c(nindiv0,nindiv1,nindiv2)
      if(sum(nindiv_total)>nindiv){
        index <- which.max(c(nindiv0,nindiv1,nindiv2))
        nindiv_total[index] <- nindiv_total[index] - 1
      }else if(sum(nindiv_total)<nindiv){
        index <- which.max(c(nindiv0,nindiv1,nindiv2))
        nindiv_total[index] <- nindiv_total[index] + 1
      }

      # Calculation
      message("Computing the statistical power...")
      result <- unlist(pbmcapply::pbmclapply(1:nsims,function(x){

        # allow more than 100 individuals
        # H1 data
        df_news <- simulatePADesignMatrix(fit = fit,
                                          df_sel = df_sel,
                                          nindiv_total = nindiv_total,
                                          model = marginal_model,
                                          snpid = snpid,
                                          nindiv = nindiv,
                                          ncell = ncell,
                                          type_specific = type_specific)
        ######
        # H0 data
        df0_news <- simulatePADesignMatrix(fit = fit0,
                                           df_sel = df_sel,
                                           nindiv_total = nindiv_total,
                                           model = marginal_model,
                                           snpid = snpid,
                                           nindiv = nindiv,
                                           ncell = ncell,
                                           type_specific = type_specific)

        # refit on H0 and H1 data
        # record the effect coefficients
        stat1 <- fitPAModel(df = df_news,model_formula = model_formula,idx = x,
                            method = method,snpid = snpid)
        stat0 <- fitPAModel(df = df0_news,model_formula = model_formula,idx = x,
                            method = method,snpid = snpid)

        return(c(as.numeric(stat1),as.numeric(stat0)))
      },mc.cores = ncores))

      stat1s <- result[seq(1,length(result),2)]
      stat0s <- result[seq(2,length(result),2)]

      if(slope < 0){

        power <- mean(stats::quantile(na.omit(stat0s),alpha) > na.omit(stat1s))

      }else if(slope > 0){

        power <- mean(stats::quantile(na.omit(stat0s),1-alpha) < na.omit(stat1s))

      }else{
        message("True eQTL effect size is zero.")
      }
      powers <- c(powers, list(data.frame(power=power,
                                          nindiv=nindiv,
                                          ncell=ncell)))
      output <- c(output,list(data.frame(stats=c(stat1s,stat0s),
                                         group=c(rep("stat1",nsims),rep("stat0",nsims)),
                                         nindiv=c(rep(nindiv,2*nsims)),
                                         ncell=c(rep(ncell,2*nsims)))))
    }
    powers <- do.call(rbind,powers)
    output <- do.call(rbind,output)
    res_tmp <- list(intercept = intercept,slope = slope,powers=powers,data=output)
    res <- c(res,list(res_tmp))
  }
  names(res) <- type_vector
  return(res)
}

#' Calculate a bootstrap confidence interval for each power
#'
#' @param res the output of function powerAnalysis()
#' @param types a vector object contains the specified type/level names of the covariate.
#' @param nindivs a vector of numeric values showing the numbers of individuals that user wants to simulate.
#' @param ncells a vector of numeric values showing the numbers of cells per each individual that user wants to simulate.
#' @param snp_number the number of SNPs for multiple testing correction.
#' @param gene_number the number of genes for multiple testing correction.
#' @param alpha the p value threshold for rejecting the H0 hypothesis.
#' @param nsim number of simulations for calculating the Bootstrap CI.
#' @param conf Bootstrap CI interval.
#'
#' @return a data frame contains average power with standard deviations in each parameter settings.
#' @export
#'
#' @examples
#' NULL
powerCICalculation <- function(res,
                               types,
                               nindivs,
                               ncells,
                               snp_number=10,
                               gene_number=800,
                               alpha=0.05,
                               nsim=1000,
                               conf=0.05){
  # check
  stopifnot("types is included in the res object. Please check input!" =
                 (types %in% names(res)))
  stopifnot("nindivs input is not a vector of numeric values. Please check input!" =
                (class(nindivs)=="numeric"))
  stopifnot("nindivs input is not a vector of integer values. Please check input!" =
                (mean(abs(nindivs - round(nindivs)) < .Machine$double.eps^0.5)==1))
  stopifnot("ncells input is not a vector of integer values. Please check input!" =
                (mean(abs(ncells - round(ncells)) < .Machine$double.eps^0.5)==1))
  stopifnot("snp_number input is not a vector of integer values. Please check input!" =
                (abs(snp_number - round(snp_number)) < .Machine$double.eps^0.5))
  stopifnot("gene_number input is not a vector of integer values. Please check input!" =
                (abs(gene_number - round(gene_number)) < .Machine$double.eps^0.5))
  stopifnot("alpha input is not between 0 and 1. Please check input!" =
                (alpha >=0 && alpha <=1))
  stopifnot("nsim input is not a integer. Please check input!" =
                (abs(nsim - round(nsim)) < .Machine$double.eps^0.5))
  stopifnot("conf input is not between 0 and 1. Please check input!" =
                (conf >=0 && conf <=1))

  # set.seed(111)
  res_new <- pbapply::pblapply(types,function(type){
      data <- res[[type]]$data
      intercept <- res[[type]]$intercept
      slope <- res[[type]]$slope
      powers <- res[[type]]$powers

      powerCIs <- list()
      for(nindiv in nindivs){
          for(ncell in ncells){
              stats1 <- data[which(data$nindiv==nindiv & data$ncell==ncell & data$group=="stat1"),"stats"]
              stats0 <- data[which(data$nindiv==nindiv & data$ncell==ncell & data$group=="stat0"),"stats"]
              power <- powers[which(powers$nindiv==nindiv & powers$ncell==ncell),"power"]

              ps <- c()
              for(i in 1:nsim){
                  stats1_tmp <- sample(stats1,length(stats1),replace = T)
                  stats0_tmp <- sample(stats0,length(stats0),replace = T)

                  if(slope < 0){

                      p <- mean(stats::quantile(na.omit(stats0_tmp),alpha/(snp_number*gene_number)) > na.omit(stats1_tmp))

                  }else if(slope > 0){

                      p <- mean(stats::quantile(na.omit(stats0_tmp),1-(alpha/(snp_number*gene_number))) < na.omit(stats1_tmp))

                  }
                  ps <- c(ps,p)
              }
              powerCI <- data.frame(power=power,
                                    nindiv=nindiv,
                                    ncell=ncell,
                                    mean=mean(ps),
                                    sd=stats::sd(ps),
                                    ci1=stats::quantile(ps,conf/2),
                                    ci2=stats::quantile(ps,(1-(conf/2))),
                                    intercept=intercept,
                                    slope=slope)
              powerCIs <- c(powerCIs,list(powerCI))
          }
      }

      powerCIs <- do.call(rbind,powerCIs)
      res_new_tmp <- res[[type]]
      res_new_tmp$powers <- powerCIs

      return(res_new_tmp)
  })

  return(res_new)
}



#' The wrapper function for power analysis
#'
#' @param marginal_list the output of function fitMarginalPop().
#' @param marginal_model a character showing the model types of the full marginal model.
#' @param refit_formula the formula used to refit the marginal full model if user wants to. Default is null.
#' @param geneid a character object contains geneid.
#' @param snpid a character object contains snpid.
#' @param type_specific a character object contains the name of the covariate that the analysis is specific on.
#' @param type_vector a vector object contains the specified type/level names of the covariate.
#' @param methods a vector of character objects specifying the methods that will be analyzed for power.
#' (Options: nb,poisson,gaussian,pseudoBulkLinear).
#' @param nindivs a vector of numeric values showing the numbers of individuals that user wants to simulate.
#' @param ncells a vector of numeric values showing the numbers of cells per each individual that user wants to simulate.
#' @param nPool a vector of numeric values showing how many pools of sequencing has been performed.
#' @param nIndivPerPool a numerical value showing how many individuals are sequenced in one pool.
#' @param nCellPerPool a vector of numeric values showing how many cells are sequenced in one pool.
#' @param alpha the p value threshold for rejecting the H0 hypothesis.
#' @param power_nsim a number of simulations for calculating the power. This parameter will affect the resolution of the power value.
#' @param snp_number the number of SNPs for multiple testing correction.
#' @param gene_number the number of genes for multiple testing correction.
#' @param CI_nsim number of simulations for calculating the Bootstrap CI.
#' @param CI_conf Bootstrap CI interval.
#' @param ncores number of CPU cores user wants to use.
#'
#' @return a data frame contains power analysis result in different parameter settings.
#' @export
#'
#' @examples
#' NULL
runPowerAnalysis <- function(marginal_list,
                             marginal_model,
                             refit_formula = NULL,
                             geneid,
                             snpid,
                             type_specific,
                             type_vector,
                             methods,
                             nindivs = NULL,
                             ncells = NULL,
                             nPool = NULL,
                             nIndivPerPool = NULL,
                             nCellPerPool = NULL,
                             alpha=0.05,
                             power_nsim=100,
                             snp_number = 10,
                             gene_number = 800,
                             CI_nsim=1000,
                             CI_conf=0.05,
                             ncores=1){
    output <- list()
    for(method in methods){
        message(paste("Performing simulations for method",method))
        res <- powerAnalysis(marginal_list = marginal_list,
                             marginal_model = marginal_model,
                             refit_formula = refit_formula,
                             geneid = geneid,
                             snpid = snpid,
                             type_specific = type_specific,
                             type_vector = type_vector,
                             method = method,
                             nindivs = nindivs,
                             ncells = ncells,
                             nPool = nPool,
                             nIndivPerPool = nIndivPerPool,
                             nCellPerPool = nCellPerPool,
                             alpha=alpha,
                             nsims=power_nsim,
                             ncores=ncores)


        message(paste("Calculating confidence intervals for method",method))
        res_CI <- powerCICalculation(res = res,
                                     types = type_vector,
                                     nindivs = nindivs,
                                     ncells = ncells,
                                     snp_number = snp_number,
                                     gene_number = gene_number,
                                     alpha = alpha,
                                     nsim = CI_nsim,
                                     conf = CI_conf)

        power_data_CI <- res_CI[[1]]$powers

        if(length(res_CI)>1){
            for(i in 1:(length(type_vector)-1)){
                power_data_CI <- rbind(power_data_CI,res_CI[[i+1]]$powers)
            }
        }

        power_data_CI$celltype <- c(rep(type_vector,each=length(res_CI[[1]]$powers[,1])))
        power_data_CI$ncell <- factor(power_data_CI$ncell)

        if(method=="nb"){
            method_name <- "NB mixed"
        }else if(method=="poisson"){
            method_name <- "Poisson mixed"
        }else if(method=="gaussian"){
            method_name <- "Linear mixed"
        }else if(method=="pseudoBulkLinear"){
            method_name <- "Pseudobulk linear"
        }else{
            stop(paste("Cannot identify input method:",method))
        }
        power_data_CI$method <- method_name
        output <- c(output,list(power_data_CI))
    }

    output <- do.call(rbind,output)
    output$nindiv <- factor(output$nindiv)
    output$method <- factor(output$method,levels = c("NB mixed","Poisson mixed","Linear mixed","Pseudobulk linear"))
    output$celltype <- factor(output$celltype)

    message("Completed!")

    return(output)
}


#' Visualize the power analysis result
#'
#' @param power_result a data frame contains power analysis result in different parameter settings.
#' @param celltypes a vector of cell types that will be selected for visualization.
#' @param x_axis a character that specifies the x axis. Default is the number of individuals.
#' @param y_axis a character that specifies the y axis. Default is the number of cells per individual.
#' @param col_group a character that specifies the color groups. Default is the eQTL model.
#' @param geneid a character object contains geneid to be part of the plot title.
#' @param snpid a character object contains snpid to be part of the plot title.
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
#' NULL
visualizePowerResult <- function(power_result,
                                 celltypes,
                                 x_axis="nindiv",
                                 y_axis="ncell",
                                 col_group="method",
                                 geneid,
                                 snpid){
    # select cell types
    power_result <- power_result[which(power_result$celltype%in%celltypes),]
    # add slope
    power_result$celltype <- as.character(power_result$celltype)

    for(i in 1:length(celltypes)){
        celltype = celltypes[i]
        power_result$celltype[which(power_result$celltype==celltype)] <- paste0(celltype,
                                                                                " (ES:",
                                                                                round(unique(as.numeric(power_result[which(power_result$celltype==celltype),"slope"])),3),
                                                                                ")")
    }
    power_result$celltype <- factor(power_result$celltype)

    p <- ggplot(power_result)
    if(x_axis=="nindiv" & y_axis=="ncell" & col_group=="method"){
        p <- p +
            geom_point(aes(nindiv,mean,color=method,shape=method),size=1,position = position_dodge(width = 0.5))+
            geom_line(aes(nindiv,mean,color=method,group=method),position = position_dodge(width = 0.5))+
            #geom_text(aes(nindiv,power,label=power),nudge_x = 10,nudge_y = 0.02)+
            geom_errorbar(aes(x=nindiv, y=mean, ymin = mean-sd,ymax=mean+sd,color=method),linewidth = 0.5,
                          width = 0.45,alpha=0.8,position = position_dodge(width = 0.5))+
            facet_grid(ncell~celltype,axes = "all",axis.labels = "margins")+
            labs(x="Number of individuals",y="Power",color="EQTL model",shape="EQTL model")+
            scale_y_continuous(
                sec.axis = sec_axis(~., name = "Number of cells per individual")
            )
    }else if(x_axis=="ncell" & y_axis=="nindiv" & col_group=="method"){
        p <- p +
            geom_point(aes(ncell,mean,color=method,shape=method),size=1,position = position_dodge(width = 0.5))+
            geom_line(aes(ncell,mean,color=method,group=method),position = position_dodge(width = 0.5))+
            #geom_text(aes(nindiv,power,label=power),nudge_x = 10,nudge_y = 0.02)+
            geom_errorbar(aes(x=ncell, y=mean, ymin = mean-sd,ymax=mean+sd,color=method),linewidth = 0.5,
                          width = 0.45,alpha=0.8,position = position_dodge(width = 0.5))+
            facet_grid(nindiv~celltype,axes = "all",axis.labels = "margins")+
            labs(x="Number of cells per individual",y="Power",color="EQTL model",shape="EQTL model")+
            scale_y_continuous(
                sec.axis = sec_axis(~., name = "Number ofindividuals")
            )
    }else{
        stop("Parameter combination unavailable.")
    }

    p <- p + theme_classic()+
        theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
        theme(
            panel.spacing = unit(1, "lines"),
            strip.background = element_blank(),
            strip.placement = "outside",
            axis.line.x = element_line(color = "black"),
            legend.title.align = 0
        )+geom_hline(yintercept = 0.8,colour = "red",linetype = "dashed")+
        ggtitle(paste(geneid,"-",snpid),subtitle = "Cell type")+
        theme(
            plot.subtitle = element_text(hjust = 0.5),
            axis.ticks.y.right = element_blank(),
            axis.text.y.right = element_blank()
        )

    return(p)
}
