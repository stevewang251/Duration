# ---------------------------------------------------------------------------- #
#
# Confidence interval for extinction duration
# from Confidence intervals for the duration of a mass extinction
# by Steve Wang, Phil Everson, Brendan McVeigh, Aaron Zimmerman, Heidi Wong
# Paleobiology 38:2, 2012
# updated version: allows for non-uniform preservation
# additional code by Madison Shoraka, Melissa Zavez, Kevin Choi, Eric Zhang
#
# ---------------------------------------------------------------------------- #



# ------------------------- Function definitions ------------------------------ #


# read in code for adaptive beta method CIs
source("~/Documents/Research/Adaptive range extensions/ abm38.3g.R")

# this is a test

# Draw range chart
plotrangechart <- function(data, ci)  {
# data: dataset with taxa in columns
# ci: vector of either theta-hat values or CI upper limits for all taxa
  ord <- order(data[1,])
  data <- data[,ord]           # sort data by highest find
  ci <- ci[ord]
  ntaxa <- dim(data)[2]
  maxfinds <- dim(data)[1]
  ymax <- max(data[1,])
  plot( NULL, pch="", xlab="", ylab="", xlim=c(0, ntaxa+1), 
        ylim=c(-.3, max(ci)), bty="L", xaxt="n" )
  for(taxon in 1:ntaxa)  {
    # draw taxon lines
    lines( c(taxon,taxon), c(0,data[1,taxon]), lwd=.9, col="darkgray" )
    # plot fossil horizons
    points(rep(taxon,maxfinds), data[,taxon], pch=16, cex=.65)
    # plot confidence interval tops
    points(taxon, ci[taxon], pch="^", cex=.65, col="red")
  }
}


# Simulate a sample from a reflected beta distribution with the given lambda
rrefbeta <- function(n, lambda)  {
  if(lambda<=0)  { 
    return(rbeta(n, 1, 1-lambda))
  }  else  return(rbeta(n, 1+lambda, 1))
}


# Simulate the gap between the highest fossil find and the true extinction time 
#   assumes theta = 1 (i.e. gaps are on the unit interval) and lambda <= 0 
simulategap <- function(lambda, n) {
  # lambda and n must be vectors of the same length (does not check for inconsistent input)
  # assumes all lambda values <= 0 (does not check)
  g <- function(x, lambda, n) (1-(1-x)^(1/n))^(1/(1-lambda))
  x <- runif(length(n))
  return(g(x, lambda, n))
}


# Simulate a dataset from thetas having true delta equal to ‘delta’, 
#   and calculate the d value from the dataset.
# Simulates only the highest finds, not the entire dataset (to save time)
simdvalues <- function(data, thetahat, ci, ymaxtaxon, y, ymax, n, ntaxa, delta, nsims, lambda, DEBUG=0)  {
  # thetahat: vector of estimated extinction times
  # ymaxtaxon: the taxon with the highest y value
  # y: vector of highest finds for all taxa
  # ymax: highest of the highest finds
  # n: vector of sample sizes for all taxa
  # ntaxa: number of taxa
  # delta: desired delta under which to simulate fossil finds
  # nsims: number of sets of simulations to run
  # lambda: vector of lambda values for all taxa
  # note that we pass in many of these values to save time instead of recalculating them
  
  simd <-  rep(NA, nsims)           # vector of simulated values
  if(DEBUG)   cat(" \n* -------------------  DELTA:", (round(delta)), " ------------------- * \n")   

  # set which taxa will have the highest and lowest thetas based on their CIs
  taxonhi <- which.max(ci)
  taxonlo <- which.min(ci)
  
  for(sim in 1:nsims)  {     

    # Choose lowest and highest theta values
    # first randomly choose gaps between highest fossil and true theta for taxonhi
#    gap <- simulategap(lambda[taxonhi], n[taxonhi]) * thetahat[taxonhi]
    numsimgaps <- 100
    gap <- quantile( simulategap(rep(lambda[taxonhi],numsimgaps), rep(n[taxonhi],numsimgaps)) * thetahat[taxonhi], .04)
    thetahi <- y[taxonhi] + gap

    # set lowest theta based on thetahi and delta
    thetalo <- thetahi - delta
    # if thetalo falls below a fossil, push thetalo and thetahi up accordingly
    if(thetalo <= y[taxonlo])  {
      thetalo <- y[taxonlo]
      thetahi <- thetalo + delta
    }

    # save lowest and highest values of theta and randomly choose intermediate values
    theta <- runif(ntaxa, thetalo, thetahi)
    theta[taxonlo] <- thetalo
    theta[taxonhi] <- thetahi
    
    if(DEBUG & sim<=5)    {  cat("\nthetas   ");  print(round((theta)))  }
    if(DEBUG & sim<=5)    {  cat(round(thetahi));  cat("  ");  cat(round(thetalo));  cat("\n")  }

    # simulate dataset from this set of thetas 
    simy <- rep(NA,ntaxa)        # vector of simulated highest finds
    for(i in 1:ntaxa)
      simy[i] <- max(rrefbeta(n[i], lambda[i])) * theta[i]
    
    # save d from this dataset
    simd[sim] <- max(simy) - min(simy)
  }
  return(simd)
}


# Draw the distribution of simulated d values on an existing plot
plotDvalues <- function(delta, nsims, simd, whichStep, color)  {
# assumes plot window has already been opened by deltaCI function
  points(rep(delta,nsims), simd, col=gray(.95), cex=.2)
  lines(rep(delta,2), quantile(simd,c(.25,.75)), col=color, lwd=4) 
  lines(rep(delta,2), quantile(simd,c(.05,.95)), col=color, lwd=2) 
  mtext(whichStep, side=3, line=0, at=delta, col=color, cex=.6)
}



# Calculate a CI for the range of extinction durations consistent with the data
deltaCI47 <- function(x, conflevel=.9, PLOT=0, ylim=NULL, xlim=NULL, nsims1=100, nsims2=100, 
             mindiff=2, abm.prAlpha=0, abm.prBeta=1, DEBUG=0)  {
  # x: dataset; assumes each taxon (column) is sorted high to low and padded with NAs
  # conflevel: CI confidence level
  # PLOT: make plot?
  # ylim, xlim: axis limits of plot
  # nsims1: no. of simulations to run at each value of delta to get sampling dist. for lower endpoint
  # nsims2: no. of simulations to run at each value of delta to get sampling dist. for upper endpoint
  # mindiff: how close the binary search steps need to be to halt

  
  # ----------------- Initialize and get dataset parameters ----------------- #

  ntaxa <- ncol(x)                      # number of taxa
  n <- rep(NA, ntaxa)                   # vector of sample sizes for each taxon
  for(i in 1:ntaxa) 
    n[i] <- length(na.omit(x[,i]))      # get sample sizes for each taxon
  y <- x[1,]                            # vector of highest fossil finds for each taxon
  ymax <- max(y)
  ymin <- min(y)
  ymaxtaxon <- which.max(y)             # taxon with highest find (assumed unique)
  obsd <- ymax - ymin                   # observed d value
  thetahat <- rep(NA, ntaxa)            # vector of estimated thetas for each taxon
  ci <- rep(NA, ntaxa)                  # vector of CIs for thetas for each taxon
  lambda <- rep(NA, ntaxa)              # vector of estimated lambdas for each taxon
  lambdavars <- rep(NA, ntaxa)          # vector of estimated variances for the lambdas of each taxon
  if(is.null(xlim))    xlim <- obsd*3   # default x-axis limits 
  if(is.null(ylim))    ylim <- ymax*2   # default y-axis limits 
  
  
  # calculate range extensions on each taxon using Adaptive Beta Method (Wang et al 2016)
  for(i in 1:ntaxa) { 
    indivconflevel = .9
    temp <- abm(na.omit(x[,i]), conf=indivconflevel, prAlpha=abm.prAlpha, prBeta=abm.prBeta)
    thetahat[i] <- temp[1]              # point estimate of extinction time
    ci[i] <- temp[3]                    # upper endpoint of CI for extinction time
    lambda[i] <- temp[4]                # estimated lambda value for each taxon
    lambdavars[i] <- temp[5]            # estimated variance of lambda values for each taxon
  }

  
  # set up plot
  if(PLOT)  {
    par(mfrow=c(2,1), mai=c(1,1,.2,.5))
    plotrangechart(x, thetahat)
    plot(0,0, type="n", xlim=c(0,xlim), ylim=c(0,ylim), xlab="delta ", ylab="simulated d values")    
    abline(h=obsd, lty=3, col=gray(.8))
    abline(v=obsd, lty=3, col=gray(.8))
  }
  
  
  # ------------- Shrink lambdas towards the mean of all lambdas ------------- #

  # calculate mean and variance of the lambdas 
  meanLambdas <- mean(lambda)
  varLambdas <- var(lambda)
  
  # calculate scaling factor to shrink lambdas towards the grand mean
  b <- lambdavars / (lambdavars + rep(varLambdas, ntaxa))
  
  # shrink the lambdas
  lambdaShrunken <- b*meanLambdas + (1-b)*lambda
  lambdaShrunken[lambdaShrunken>0] <- 0    # lambdas assumed <= 0; otherwise, set to 0

# ignore this
#lambdaShrunken <- (lambdaShrunken + 1)*2 + -1
#lambdaShrunken <- pmin(lambdaShrunken, 0)


  # make shrinkage plot
  if(PLOT)  {
    pdf("results abm lambda_hats.pdf", w=8.5, h=11)
    plot(lambda, rep(1,ntaxa), ylim = c(-1, 2))
    points(lambdaShrunken, rep(0, ntaxa))
    for(i in 1:ntaxa) 
      segments(lambda[i], 1, lambdaShrunken[i],0)
    dev.off()
  }
  

  # ------------- Calculate lower endpoint of confidence interval ------------- #

  # loop starting at currdelta = 0; low and high are the current bounds of the binary search
  low <- 0
  high <- 1.2*obsd          # assumes obsd must be in the CI, so lower endpoint <= obsd
  currdelta <- low          # the value of delta we are currently checking
  difference <- high - low
  whichStep <- 0            # indicates how many steps have we taken
  
  # loop starting at currdelta = 0, iterate until we converge on the lower endpoint
  while(difference > mindiff)  {
  
    # simulate d values and determine cutoff for acceptance
    whichStep <- whichStep + 1
    simd <- simdvalues(x, thetahat, ci, ymaxtaxon, y, ymax, n, ntaxa, currdelta, nsims1, lambdaShrunken, DEBUG)
    temp <- quantile(simd, c( (1-conflevel)/2,  1-(1-conflevel)/2 ))
    lowercutoff <- temp[1]                              # lower-tail 1-sided rejection region
    uppercutoff <- temp[2];  rm(temp)                   # upper-tail 1-sided rejection region
    
    # check if our observed d value is accepted (i.e., consistent w/ the simulated d values)
    if (obsd <= uppercutoff)  {                         # if so, search to the left
      high <- currdelta    # make the current delta the new upper bound, keep the old lower bound
      if(PLOT)  plotDvalues(currdelta, nsims1, simd, whichStep, "red") 
    } 
    else {    # if not accepted (observed d value not consistent with simulated values), search to the right
      low <- currdelta     # make the current delta the new lower bound, keep the old upper bound
      if(PLOT)  plotDvalues(currdelta, nsims1, simd, whichStep, "pink") 
    } 
    
    # set the new delta value and update the difference between the lower and upper search bounds
    currdelta <- (low + high)/2 
    difference = high - low     

  }  # while 
  # when we exit this while loop, we have converged to the lower edge of the acceptance region
  
  # save lower endpoint
  deltamin <- currdelta
  if(PLOT)    points(deltamin,0, pch=6, col=gray(.7))

  
  # ------------- Calculate upper endpoint of confidence interval ------------- #

  low <- deltamin + mindiff
  high <- ymax * (-min(lambda) + 1)      # upper search limit 
  currdelta <- low
  difference <- high - low
  whichStep <- 0
  
  # loop starting at currdelta = obsd, iterate until we converge on the upper endpoint
  while(difference > mindiff)  {  

    # simulate d values and determine cutoff for acceptance
    whichStep <- whichStep + 1    
    simd <- simdvalues(x, thetahat, ci, ymaxtaxon, y, ymax, n, ntaxa, currdelta, nsims2, lambdaShrunken, DEBUG)
    temp <- quantile(simd, c( (1-conflevel)/2,  1-(1-conflevel)/2 ))
    lowercutoff <- temp[1]                              # lower-tail 1-sided rejection region
    uppercutoff <- temp[2];  rm(temp)                   # upper-tail 1-sided rejection region

    # check if our observed d value is accepted (i.e., consistent w/ the simulated d values)
    if (lowercutoff <= obsd)  {        # if so, search to the right
      low <- currdelta      # make the current delta the new lower bound, keep the old upper bound
      if(PLOT)  plotDvalues(currdelta, nsims2, simd, whichStep, "blue") 
    } 
    else {    # if not accepted (observed d value not consistent with simulated values), search to the left
      high <- currdelta     # make the current delta the new upper bound, keep the old lower bound
      if(PLOT)  plotDvalues(currdelta, nsims2, simd, whichStep, "steel blue 2") 
    }

    # set the new delta value and update the difference between the lower and upper search bounds
    currdelta <- (low + high)/2
    difference = high - low

  }  # while
  # when we exit this while loop, we have converged to the upper edge of the acceptance region
  
  # save upper endpoint
  deltamax <- currdelta
  if(PLOT)    points(deltamax,0, pch=6, col=gray(.7))

    
  # -------- Return CI lower and upper bound, and observed value of d -------- #

  deltamin <- min(deltamin, obsd)
  deltamax <- max(deltamax, obsd)
  return(c(deltamin,deltamax, obsd, mean(lambdaShrunken)))

}


