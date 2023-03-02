mcmcMH <- function(target, init.theta, proposal.sd = NULL,
                   n.iterations, covmat = NULL,
                   limits=list(lower = NULL, upper = NULL),
                   adapt.size.start = NULL, adapt.size.cooling = 0.99,
                   adapt.shape.start = NULL, adapt.shape.stop = NULL,
                   print.info.every = n.iterations/100,
                   verbose = FALSE, max.scaling.sd = 50) {
  
  # initialise theta
  theta.current <- init.theta
  theta.propose <- init.theta
  
  # extract theta of gaussian proposal
  covmat.proposal <- covmat
  lower.proposal <- limits$lower
  upper.proposal <- limits$upper
  
  # reorder vector and matrix by names, set to default if necessary
  theta.names <- names(init.theta)
  if (!is.null(proposal.sd) && is.null(names(proposal.sd))) {
    names(proposal.sd) <- theta.names
  }
  
  if (is.null(covmat.proposal)) {
    if (is.null(proposal.sd)) {
      proposal.sd <- init.theta/10
    }
    covmat.proposal <-
      matrix(diag(proposal.sd[theta.names]^2, nrow = length(theta.names)),
             nrow = length(theta.names),
             dimnames = list(theta.names, theta.names))
  } else {
    covmat.proposal <- covmat.proposal[theta.names,theta.names]
  }
  
  if (is.null(lower.proposal)) {
    lower.proposal <- init.theta
    lower.proposal[] <- -Inf
  } else {
    lower.proposal <- lower.proposal[theta.names]
  }
  
  if (is.null(upper.proposal)) {
    upper.proposal <- init.theta
    upper.proposal[] <- Inf
  } else {
    upper.proposal <- upper.proposal[theta.names]
  }
  
  # covmat init
  covmat.proposal.init <- covmat.proposal
  
  adapting.size <- FALSE # will be set to TRUE once we start
  # adapting the size
  
  adapting.shape <- 0  # will be set to the iteration at which
  # adaptation starts
  
  # find estimated theta
  theta.estimated.names <- names(which(diag(covmat.proposal) > 0))
  
  # evaluate target at theta init
  target.theta.current <- target(theta.current)
  
  if (!is.null(print.info.every)) {
    message(Sys.time(), ", Init: ", printNamedVector(theta.current[theta.estimated.names]),
            ", target: ", target.theta.current)
  }
  
  # trace
  trace <- matrix(ncol=length(theta.current)+1, nrow=n.iterations, 0)
  colnames(trace) <- c(theta.estimated.names, "log.density")
  
  # acceptance rate
  acceptance.rate <- 0
  
  # scaling factor for covmat size
  scaling.sd  <- 1
  
  # scaling multiplier
  scaling.multiplier <- 1
  
  # empirical covariance matrix (0 everywhere initially)
  covmat.empirical <- covmat.proposal
  covmat.empirical[,] <- 0
  
  # empirical mean vector
  theta.mean <- theta.current
  
  # if print.info.every is null never print info
  if (is.null(print.info.every)) {
    print.info.every <- n.iterations + 1
  }
  
  start_iteration_time <- Sys.time()
  
  for (i.iteration in seq_len(n.iterations)) {
    
    # adaptive step
    if (!is.null(adapt.size.start) && i.iteration >= adapt.size.start &&
        (is.null(adapt.shape.start) || acceptance.rate*i.iteration < adapt.shape.start)) {
      if (!adapting.size) {
        message("\n---> Start adapting size of covariance matrix")
        adapting.size <- TRUE
      }
      # adapt size of covmat until we get enough accepted jumps
      scaling.multiplier <- exp(adapt.size.cooling^(i.iteration-adapt.size.start) * (acceptance.rate - 0.234))
      scaling.sd <- scaling.sd * scaling.multiplier
      scaling.sd <- min(c(scaling.sd,max.scaling.sd))
      # only scale if it doesn't reduce the covariance matrix to 0
      covmat.proposal.new <- scaling.sd^2*covmat.proposal.init
      if (!(any(diag(covmat.proposal.new)[theta.estimated.names] <
                .Machine$double.eps))) {
        covmat.proposal <- covmat.proposal.new
      }
      
    } else if (!is.null(adapt.shape.start) &&
               acceptance.rate*i.iteration >= adapt.shape.start &&
               (adapting.shape == 0 || is.null(adapt.shape.stop) ||
                i.iteration < adapting.shape + adapt.shape.stop)) {
      if (!adapting.shape) {
        message("\n---> Start adapting shape of covariance matrix")
        # flush.console()
        adapting.shape <- i.iteration
      }
      
      ## adapt shape of covmat using optimal scaling factor for multivariate target distributions
      scaling.sd <- 2.38/sqrt(length(theta.estimated.names))
      
      covmat.proposal <- scaling.sd^2 * covmat.empirical
    } else if (adapting.shape > 0) {
      message("\n---> Stop adapting shape of covariance matrix")
      adapting.shape <- -1
    }
    
    # print info
    if (i.iteration %% ceiling(print.info.every) == 0) {
      message(Sys.time(), ", Iteration: ",i.iteration,"/", n.iterations,
              ", acceptance rate: ",
              sprintf("%.3f",acceptance.rate), appendLF=FALSE)
      if (!is.null(adapt.size.start) || !is.null(adapt.shape.start)) {
        message(", scaling.sd: ", sprintf("%.3f", scaling.sd),
                ", scaling.multiplier: ", sprintf("%.3f", scaling.multiplier),
                appendLF=FALSE)
      }
      message(", state: ",(printNamedVector(theta.current)))
      message(", logdensity: ", target.theta.current)
    }
    
    # propose another parameter set
    if (any(diag(covmat.proposal)[theta.estimated.names] <
            .Machine$double.eps)) {
      print(covmat.proposal[theta.estimated.names,theta.estimated.names])
      stop("non-positive definite covmat",call.=FALSE)
    }
    if (length(theta.estimated.names) > 0) {
      theta.propose[theta.estimated.names] <-
        as.vector(rtmvnorm(1,
                           mean =
                             theta.current[theta.estimated.names],
                           sigma =
                             covmat.proposal[theta.estimated.names,theta.estimated.names],
                           lower =
                             lower.proposal[theta.estimated.names],
                           upper = upper.proposal[theta.estimated.names]))
    }
    
    # evaluate posterior of proposed parameter
    target.theta.propose <- target(theta.propose)
    # if return value is a vector, set log.density and trace
    
    if (!any(is.finite(target.theta.propose))) { # GK: changed to include "any" as getting error with list
      # if posterior is 0 then do not compute anything else and don't accept
      log.acceptance <- -Inf
      
    }else{
      
      # compute Metropolis-Hastings ratio (acceptance probability)
      log.acceptance <- target.theta.propose - target.theta.current
      log.acceptance <- log.acceptance +
        dtmvnorm(x = theta.current[theta.estimated.names],
                 mean =
                   theta.propose[theta.estimated.names],
                 sigma =
                   covmat.proposal[theta.estimated.names,
                                   theta.estimated.names],
                 lower =
                   lower.proposal[theta.estimated.names],
                 upper =
                   upper.proposal[theta.estimated.names],
                 log = TRUE)
      log.acceptance <- log.acceptance -
        dtmvnorm(x = theta.propose[theta.estimated.names],
                 mean = theta.current[theta.estimated.names],
                 sigma =
                   covmat.proposal[theta.estimated.names,
                                   theta.estimated.names],
                 lower =
                   lower.proposal[theta.estimated.names],
                 upper =
                   upper.proposal[theta.estimated.names],
                 log = TRUE)
      
    }
    
    if (verbose) {
      message("Propose: ", theta.propose[theta.estimated.names],
              ", target: ", target.theta.propose,
              ", acc prob: ", exp(log.acceptance), ", ",
              appendLF = FALSE)
    }
    
    if (is.accepted <- (log(runif (1)) < log.acceptance)) {
      # accept proposed parameter set
      theta.current <- theta.propose
      target.theta.current <- target.theta.propose
      if (verbose) {
        message("accepted")
      }
    } else if (verbose) {
      message("rejected")
    }
    trace[i.iteration, ] <- c(theta.current, target.theta.current)
    
    # update acceptance rate
    if (i.iteration == 1) {
      acceptance.rate <- is.accepted
    } else {
      acceptance.rate <- acceptance.rate +
        (is.accepted - acceptance.rate) / i.iteration
    }
    
    # update empirical covariance matrix
    if (adapting.shape >= 0) {
      tmp <- updateCovmat(covmat.empirical, theta.mean,
                          theta.current, i.iteration)
      covmat.empirical <- tmp$covmat
      theta.mean <- tmp$theta.mean
    }
    
  }
  
  return(list(trace = trace,
              acceptance.rate = acceptance.rate,
              covmat.empirical = covmat.empirical))
}


#'Print named vector
#'
#'Print named vector with format specified by \code{fmt} (2 decimal places by default).
#' @param x named vector
#' @inheritParams base::sprintf
#' @inheritParams base::paste
#' @export
#' @seealso \code{\link[base]{sprintf}}
#' @keywords internal
printNamedVector <- function(x, fmt="%.2f", sep=" | ") {
  
  paste(paste(names(x),sprintf(fmt,x),sep=" = "),collapse=sep)
  
}

#' Simulate forward a stochastic model
#'
#' This function uses the function \code{\link[adaptivetau]{ssa.adaptivetau}} to simulate the model and returns the trajectories in a valid format for the class \code{\link{fitmodel}}.
#' @param theta named vector of model parameters.
#' @param init.state named vector of initial state of the model.
#' @param times time sequence for which state of the model is wanted; the first value of times must be the initial time.
#' @inheritParams adaptivetau::ssa.adaptivetau
#' @export
#' @import adaptivetau
#' @return a data.frame of dimension \code{length(times)x(length(init.state)+1)} with column names equal to \code{c("time",names(init.state))}.
simulateModelStochastic <- function(theta,init.state,times,transitions,rateFunc) {
  
  
  stoch <- as.data.frame(ssa.adaptivetau(init.state,transitions,rateFunc,theta,tf=diff(range(times))))
  
  # rescale time as absolute value
  stoch$time <- stoch$time + min(times)
  
  # interpolate
  traj <- cbind(time=times,apply(stoch[,-1],2,function(col){approx(x=stoch[,1],y=col,xout=times,method="constant")$y}))
  
  return(as.data.frame(traj))
  
}



#'Simulate several replicate of the model
#'
#'Simulate several replicate of a fitmodel using its function simulate
#' @param times vector of times at which you want to observe the states of the model.
#' @param n number of replicated simulations.
#' @param observation logical, if \code{TRUE} simulated observation are generated by \code{rTrajObs}.
#' @inheritParams testFitmodel
#' @export
#' @import plyr
#' @return a data.frame of dimension \code{[nxlength(times)]x[length(init.state)+2]} with column names equal to \code{c("replicate","time",names(init.state))}.
simulateModelReplicates <- function(fitmodel,theta, init.state, times, n, observation=FALSE) {
  
  stopifnot(inherits(fitmodel,"fitmodel"),n>0)
  
  if(observation && is.null(fitmodel$dPointObs)){
    stop("Can't generate observation as ",sQuote("fitmodel")," doesn't have a ",sQuote("dPointObs")," function.")
  }
  
  rep <- as.list(1:n)
  names(rep) <- rep
  
  if (n > 1) {
    progress = "text"
  } else {
    progress = "none"
  }
  
  traj.rep <- ldply(rep,function(x) {
    
    if(observation){
      traj <- rTrajObs(fitmodel, theta, init.state, times)
    } else {
      traj <- fitmodel$simulate(theta,init.state,times)
    }
    
    return(traj)
    
  },.progress=progress,.id="replicate")
  
  return(traj.rep)
}


#'Simulate model until extinction
#'
#'Return final state at extinction
#' @param extinct character vetor. Simulations stop when all these state are extinct.
#' @param time.init numeric. Start time of simulation.
#' @param time.step numeric. Time step at which extinction is checked
#' @inheritParams testFitmodel
#' @inheritParams simulateModelReplicates
#' @inheritParams particleFilter
#' @export
#' @import plyr parallel doParallel
simulateFinalStateAtExtinction <- function(fitmodel, theta, init.state, extinct=NULL ,time.init=0, time.step=100, n=100, observation=FALSE, n.cores = 1) {
  
  stopifnot(inherits(fitmodel,"fitmodel"),n>0)
  
  if(observation && is.null(fitmodel$dPointObs)){
    stop("Can't generate observation as ",sQuote("fitmodel")," doesn't have a ",sQuote("dPointObs")," function.")
  }
  
  if(is.null(n.cores)){
    n.cores <- detectCores()
  }
  
  if(n.cores > 1){
    registerDoParallel(cores=n.cores)
  }
  
  rep <- as.list(1:n)
  names(rep) <- rep
  
  if (n > 1 && n.cores==1) {
    progress = "text"
  } else {
    progress = "none"
  }
  
  times <- c(time.init, time.step)
  
  final.state.rep <- ldply(rep,function(x) {
    
    if(observation){
      traj <- rTrajObs(fitmodel, theta, init.state, times)
    } else {
      traj <- fitmodel$simulate(theta,init.state,times)
    }
    
    current.state <- unlist(traj[nrow(traj),fitmodel$state.names])
    current.time <- last(traj$time)
    
    while(any(current.state[extinct]>=0.5)){
      
      times <- times + current.time
      
      if(observation){
        traj <- rTrajObs(fitmodel, theta, current.state, times)
      } else {
        traj <- fitmodel$simulate(theta, current.state,times)
      }
      
      current.state <- unlist(traj[nrow(traj),fitmodel$state.names])
      current.time <- last(traj$time)
    }
    
    return(data.frame(t(c(time=current.time,current.state))))
    
  },.progress=progress,.id="replicate",.parallel=(n.cores > 1),.paropts=list(.inorder=FALSE))
  
  return(final.state.rep)
}


#'Update covariance matrix
#'
#'Update covariance matrix using a stable one-pass algorithm. This is much more efficient than using \code{\link{cov}} on the full data.
#' @param covmat covariance matrix at iteration \code{i-1}. Must be numeric, symmetrical and named.
#' @param theta.mean mean vector at iteration \code{i-1}. Must be numeric and named.
#' @param theta vector of new value at iteration \code{i}. Must be numeric and named.
#' @param i current iteration.
#' @references \url{http://en.wikipedia.org/wiki/Algorithms\%5Ffor\%5Fcalculating\%5Fvariance#Covariance}
#' @export
#' @keywords internal
#' @return A list of two elements
#' \itemize{
#' \item \code{covmat} update covariance matrix
#' \item \code{theta.mean} updated mean vector
#' }
updateCovmat <- function(covmat,theta.mean,theta,i) {
  
  if(is.null(names(theta))){
    stop("Argument ",sQuote("theta")," must be named.",.call=FALSE)
  }
  if(is.null(names(theta.mean))){
    stop("Argument ",sQuote("theta.mean")," must be named.",.call=FALSE)
  }
  if(is.null(rownames(covmat))){
    stop("Argument ",sQuote("covmat")," must have named rows.",.call=FALSE)
  }
  if(is.null(colnames(covmat))){
    stop("Argument ",sQuote("covmat")," must have named columns.",.call=FALSE)
  }
  
  covmat <- covmat[names(theta),names(theta)]
  theta.mean <- theta.mean[names(theta)]
  
  residual <- as.vector(theta-theta.mean)
  covmat <- (covmat*(i-1)+(i-1)/i*residual%*%t(residual))/i
  theta.mean <- theta.mean + residual/i
  
  return(list(covmat=covmat,theta.mean=theta.mean))
}



#'Burn and thin MCMC chain
#'
#'Return a burned and thined trace of the chain.
#' @param trace either a \code{data.frame} or a \code{list} of \code{data.frame} with all variables in column, as outputed by \code{\link{mcmcMH}}. Accept also an \code{mcmc} or \code{mcmc.list} object.
#' @param burn proportion of the chain to burn.
#' @param thin number of samples to discard per sample that is being kept
#' @export
#' @import coda
#' @return an object with the same format as \code{trace} (\code{data.frame} or \code{list} of \code{data.frame} or \code{mcmc} object or \code{mcmc.list} object)
burnAndThin <- function(trace, burn = 0, thin = 0) {
  
  convert_to_mcmc <- FALSE
  
  if(class(trace)=="mcmc"){
    convert_to_mcmc <- TRUE
    trace <- as.data.frame(trace)
  } else if(class(trace)=="mcmc.list"){
    convert_to_mcmc <- TRUE
    trace <- as.list(trace)
  }
  
  if(is.data.frame(trace) || is.matrix(trace)){
    
    # remove burn
    if (burn > 0) {
      trace <- trace[-(1:burn), ]
    }
    # thin
    trace <- trace[seq(1, nrow(trace), thin + 1), ]
    
    if(convert_to_mcmc){
      trace <- mcmc(trace)
    }
    
  } else {
    
    trace <- lapply(trace, function(x) {
      
      # remove burn
      if (burn > 0) {
        x <- x[-(1:burn), ]
      }
      # thin
      x <- x[seq(1, nrow(x), thin + 1), ]
      
      if(convert_to_mcmc){
        x <- mcmc(x)
      }
      
      return(x)
    }) 
    
    if(convert_to_mcmc){
      trace <- mcmc.list(trace)            
    }
  }
  
  return(trace)
}


#'Distance weighted by number of oscillations
#'
#'This positive distance is the mean squared differences between \code{x} and the \code{y}, divided by the square of the number of times the \code{x} oscillates around the \code{y} (see note below for illustration).
#' @param x,y numeric vectors of the same length.
#' @note To illustrate this distance, suppose we observed a time series \code{y = c(1,3,5,7,5,3,1)} and we have two simulated time series \code{x1 = (3,5,7,9,7,5,3)} and \code{x2 = (3,5,3,5,7,5,3)}; \code{x1} is consistently above \code{y} and \code{x2} oscillates around \code{y}. While the squared differences are the same, we obtain \eqn{d(y, x1) = 4} and \eqn{d(y, x2) = 1.3}.
#' @export
distanceOscillation <- function(x, y) {
  
  # check x and y have same length
  if(length(x)!=length(y)){
    stop(sQuote("x")," and ",sQuote("y")," must be vector of the same length")
  }
  
  # 1 + number of times x oscillates around y
  n.oscillations <- 1+length(which(diff((x-y)>0)!=0))
  
  dist <- sum((x-y)^2)/(length(x)*n.oscillations)
  
  return(dist)
}


#'Export trace in Tracer format
#'
#'Print \code{trace} in a \code{file} that can be read by the software Tracer.
#' @param trace a \code{data.frame} with one column per estimated parameter, as returned by \code{\link{burnAndThin}}
#' @inheritParams utils::write.table
#' @note Tracer is a program for analysing the trace files generated by Bayesian MCMC runs. It can be dowloaded at \url{http://tree.bio.ed.ac.uk/software/tracer}.
#' @export
#' @seealso burnAndThin
#' @keywords internal
export2Tracer <- function(trace, file) {
  
  if(is.mcmc(trace)){
    trace <- as.data.frame(trace)
  }
  
  if(!"iteration"%in%names(trace)){
    trace$iteration <- (1:nrow(trace) - 1)
  }
  
  trace <- trace[c("iteration",setdiff(names(trace),c("iteration","weight")))]
  write.table(trace,file=file,quote=FALSE,row.names=FALSE,sep="\t")
  
}