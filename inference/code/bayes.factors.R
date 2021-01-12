
### bayes.factors.R
### Nathan Wikle

### Calculate Bayes Factors for multiple models, using Importance Sampling
###   (posterior fits used for importance distribution). These functions only work
###   on output from mcmc.odesim.2.0.R updated on July 18th.

################################################################################
###### 1. Functions
################################################################################

library(matrixStats)

### A. Grab samples from posterior output
sample.grab <- function(S, directory, burnin = 0){
  # Input: 
  #   S: number of samples to grab from posterior output.
  #   directory: location of posterior output.
  #   burnin: number of burn-in samples to leave out.
  # Output:
  #   A sample of S log-likelihood evaluations saved with posterior output.
  
  # grab mcmc output
  mcmc.output <- mcmc.params(out.folder = directory, burnin = burnin)
  inds <- dim(mcmc.output$beta)

  # subsample the mcmc output
  subsamples <- cbind(sample(1:inds[1], size = S),
                        sample(1:inds[3], size = S, replace = T))
  
  # store results as a list
  res <- list()
  
  res$betas <- t(sapply(1:S, function(x){mcmc.output$beta[subsamples[x,1],,subsamples[x,2]]}))
  res$ode <- t(sapply(1:S, function(x){mcmc.output$ode[subsamples[x,1],,subsamples[x,2]]}))
  res$rr <- t(sapply(1:S, function(x){mcmc.output$rr[subsamples[x,1],,subsamples[x,2]]}))
  res$lik.params <- t(sapply(1:S, function(x){mcmc.output$lik.params[subsamples[x,1],,subsamples[x,2]]}))
  res$loglik <- unlist(sapply(1:S, function(x){mcmc.output$loglik[subsamples[x,1],1,subsamples[x,2]]}))
  
  res$days <- mcmc.output$days
  res$Z.beta <- mcmc.output$Z.beta
  res$Z.rr <- mcmc.output$Z.rr
  colnames(res$ode) <- mcmc.output$col.names
  res$df <- mcmc.output$df
  res$loc <- mcmc.output$loc
  res$const <- mcmc.output$const
  res$introday <- mcmc.output$introday
  
  return(res)
}

mcmc.params <- function(out.folder, burnin = 0){
  # Input:
  #   out.folder: location of posterior output.
  #   burnin: number of burn-in samples to leave out.
  # Output:
  #   All MCMC parameters saved after burnin.
  
  # save output files in a list
  files.list <- list.files(out.folder,"*.Rdata")
  files.list <- files.list[order(nchar(files.list), files.list)]
  max.iter=length(files.list)
  
  # load mcmc output
  load(paste(out.folder, files.list[1], sep=""))
  data <- out
  
  # initialize parameters
  n.chains <- length(data)
  n.mcmc <- nrow(data[[1]]$beta)
  n.beta <- ncol(data[[1]]$beta)
  n.ode.params <- ncol(data[[1]]$ode.params)
  n.rr.params <- ncol(data[[1]]$rr.params)
  n.lik.params <- ncol(data[[1]]$lik.params)
  spline <- data[[1]]$spline
  df <- data[[1]]$df
  
  # data structures
  beta.chains=array(NA,dim=c(n.mcmc*max.iter,n.beta,n.chains))
  ode.chains=array(NA,dim=c(n.mcmc*max.iter,n.ode.params,n.chains))
  rr.chains=array(NA,dim=c(n.mcmc*max.iter,n.rr.params,n.chains))
  lik.chains=array(NA,dim=c(n.mcmc*max.iter,n.lik.params,n.chains))
  lik.vals=array(NA, c(n.mcmc*max.iter,n.lik.params,n.chains))
  
  for(iter in 1:max.iter){
    ## load in data and call it "data"
    load(paste(out.folder,files.list[iter],sep=""))
    data=out
    ## get betas and ode.params
    for(i in 1:n.chains){
      beta.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$beta
      ode.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$ode.params
      rr.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$rr.params
      lik.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$lik.params
      lik.vals[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$loglik
    }
  }
  
  # return list of values
  results <- list()
  if (burnin > 0){
    results$beta <- beta.chains[-c(1:burnin),,]
    results$ode <- ode.chains[-c(1:burnin),,]
    results$rr <- rr.chains[-c(1:burnin),,]
    results$lik.params <- lik.chains[-c(1:burnin),,]
    results$loglik <- lik.vals[-c(1:burnin),,]
  } else {
    results$beta <- beta.chains
    results$ode <- ode.chains
    results$rr <- rr.chains
    results$lik.params <- lik.chains
    results$loglik <- lik.vals
  }
  
  library(fda)
  df=data[[1]]$df
  results$days <- df$daynum
  results$Z.beta <- eval.basis(61:max(results$days),data[[1]]$spline.beta)
  results$Z.rr <- eval.basis(61:max(results$days),data[[1]]$spline.rr)
  results$col.names <- colnames(data[[1]]$ode.params)
  results$df <- df
  results$loc <- data[[1]]$loc
  results$const <- data[[1]]$const.params
  results$introday <- data[[1]]$introday
  
  return(results)
}


### B. Calculate IS estimate of p(D | M)
IS.estimate <- function(ll){
  # Input:
  #   ll: sample of log-likelihood evaluations.
  # Output:
  #   IS estimate of log(p(D|M)), and estimate of Monte Carlo SE.
  
  # number of samples
  S <- length(ll)
  # log(p(D|M))
  log.marg.dens <- log(S) - logSumExp(-ll)
  # Monte Carlo SE estimate
  log.MCSE <- log.marg.dens + logSumExp(c(logSumExp(-2*ll) - 2 * logSumExp(-ll), log(S))) / 2
    
  # return log density and log st. error
  res <- list(log.marg.dens = log.marg.dens,
              log.MCSE = log.MCSE)
  return(res)
}


psis.estimate <- function(ll){
  
}
















