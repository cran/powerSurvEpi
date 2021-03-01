# created on Feb. 20, 2021
#  (1) power and sample size for score test of conditional logistic regression


###########
# Lachin2008, Section 3.3, Formula (38)
###########
# N - number of sets (total number of subjects = N * (nD + nH))
# power - power of the test
# OR - odds ratio = exp(theta), which is regression coefficient of the exposure variable
# pE - population prevalence of exposure
# nD - number of cases per set
# nH - number of controls per set
# R2 - coefficient of determination of the exposure variable and other covariates

# binary covariate
powerConLogistic.bin = function(N = NULL, 
                            power = 0.8, 
                            OR, 
                            pE, 
                            nD, 
                            nH,
                            R2 = 0,
                            alpha = 0.05, 
                            nTests = 1,
                            OR.low = 1.01,
                            OR.upp = 100
)
{
  alpha2 = alpha/2
  za = qnorm(1-alpha2/nTests)
 
  if(is.null(power) == TRUE && is.null(N) == FALSE && is.null(OR) == FALSE)
  {
    theta = log(OR)
    corePart = theta^2*pE*(1-pE)*(1-R2)*nD*nH/(nD+nH)
 
    # calculate power
    power = pnorm(sqrt(N*corePart) - za)

    return(power)
  } else if(is.null(power) == FALSE && is.null(N) == TRUE && is.null(OR) == FALSE) {
    theta = log(OR)
    corePart = theta^2*pE*(1-pE)*(1-R2)*nD*nH/(nD+nH)

    # calculate sample size
    zb = qnorm(power)
    numer = (zb + za)^2
    denom = corePart
    N = numer/denom
    
    return(N)
  } else if(is.null(power) == FALSE && is.null(N) == FALSE && is.null(OR) == TRUE) {

    # find minimum detectable OR
    OR = powerConLogistic.bin.est.OR(
      N = N, 
      power = power, 
      pE = pE, 
      nD = nD, 
      nH = nH,
      R2 = R2,
      alpha = alpha, 
      nTests = nTests,
      OR.low = OR.low,
      OR.upp = OR.upp
    )

    return(OR)

  } else {
    stop("one and only one of power, N, and OR should be null!")
  }
  
}

# difference between desired power and estimated power given theta = log(OR)
diffPower.ConLogistic.bin.OR = function(
  theta,
  N, 
  power = 0.8, 
  pE, 
  nD, 
  nH,
  R2 = 0,
  alpha=0.05, 
  nTests=1
)
{
  alpha2 = alpha/2
  za = qnorm(1-alpha2/nTests)
  corePart = theta^2*pE*(1-pE)*(1-R2)*nD*nH/(nD+nH)
  
  # calculate power
  power.est = pnorm(sqrt(N*corePart) - za)

  diff = power - power.est

  return(diff)
}

powerConLogistic.bin.est.OR = function(
  N, 
  power = 0.8, 
  pE, 
  nD, 
  nH,
  R2 = 0,
  alpha=0.05, 
  nTests=1,
  OR.low = 1.01,
  OR.upp = 100
)
{
  theta.low = log(OR.low)
  theta.upp = log(OR.upp)

  aa = try(res <- uniroot(f = diffPower.ConLogistic.bin.OR, interval = c(theta.low, theta.upp), 
    N = N,
    power = power,
    pE = pE,
    nD = nD,
    nH = nH,
    R2 = R2,
    alpha = alpha,
    nTests = nTests
  ), silent = FALSE)

  res.aa = attr(aa, which = "class")
  if(is.null(res.aa) == FALSE)
  {
    OR = NA
    cat("Error occurs when searhing optimal OR!\nPlease adjust range [OR.low, OR.upp]!\n")
  } else {
    OR = exp(res$roo)
  }

  return(OR)
}


####
###########
# Lachin2008, Section 3.2.1, Formula (24) and (25)
###########

# continuous covariate

# N - number of sets (total number of subjects = N * (nD + nH))
# power - power of the test
# OR - odds ratio = exp(theta), which is regression coefficient of the exposure variable
# sigma - standard deviation of the continuous covariate
# nD - number of cases per set
# nH - number of controls per set
# R2 - coefficient of determination of the exposure variable and other covariates

powerConLogistic.con = function(N = NULL, 
                                power = 0.8, 
                                OR, 
                                sigma, 
                                nD, 
                                nH,
                                R2 = 0,
                                alpha = 0.05, 
                                nTests = 1,
                                OR.low = 1.01,
                                OR.upp = 100
)
{
  alpha2 = alpha/2
  za = qnorm(1-alpha2/nTests)
  n = nD + nH
  nchoosed = nchoosek(n, nD)
  
  if(is.null(power) == TRUE && is.null(N) == FALSE && is.null(OR) == FALSE)
  {
    theta = log(OR)
    corePart = theta^2*sigma^2*nD*(1-1/nchoosed)*(1-R2)
    # calculate power
    power = pnorm(sqrt(N*corePart) - za)
    
    return(power)
  } else if(is.null(power) == FALSE && is.null(N) == TRUE && is.null(OR) == FALSE) {

    theta = log(OR)
    corePart = theta^2*sigma^2*nD*(1-1/nchoosed)*(1-R2)

    # calculate sample size
    zb = qnorm(power)
    numer = (zb + za)^2
    denom = corePart
    N = numer/denom
    
    return(N)
  } else if(is.null(power) == FALSE && is.null(N) == FALSE && is.null(OR) == TRUE) {

    # find minimum detectable OR
    OR = powerConLogistic.con.est.OR(
      N = N, 
      power = power, 
      sigma = sigma, 
      nD = nD, 
      nH = nH,
      R2 = R2,
      alpha = alpha, 
      nTests = nTests,
      OR.low = OR.low,
      OR.upp = OR.upp
    )

    return(OR)

  } else {
    stop("one and only one of power, N, and OR should be null!")
  }
  
}

# difference between desired power and estimated power given theta = log(OR)
diffPower.ConLogistic.con.OR = function(
  theta,
  N, 
  power = 0.8, 
  sigma, 
  nD, 
  nH,
  R2 = 0,
  alpha=0.05, 
  nTests=1
)
{
  alpha2 = alpha/2
  za = qnorm(1-alpha2/nTests)
  n = nD + nH
  nchoosed = nchoosek(n, nD)
 
  corePart = theta^2*sigma^2*nD*(1-1/nchoosed)*(1-R2)
  # calculate power
  power.est = pnorm(sqrt(N*corePart) - za)
 
  diff = power - power.est

  return(diff)
}

powerConLogistic.con.est.OR = function(
  N, 
  power = 0.8, 
  sigma, 
  nD, 
  nH,
  R2 = 0,
  alpha=0.05, 
  nTests=1,
  OR.low = 1.01,
  OR.upp = 100
)
{
  theta.low = log(OR.low)
  theta.upp = log(OR.upp)

  aa = try(res <- uniroot(f = diffPower.ConLogistic.con.OR, interval = c(theta.low, theta.upp), 
    N = N,
    power = power,
    sigma = sigma,
    nD = nD,
    nH = nH,
    R2 = R2,
    alpha = alpha,
    nTests = nTests
  ), silent = FALSE)

  res.aa = attr(aa, which = "class")
  if(is.null(res.aa) == FALSE)
  {
    OR = NA
    cat("Error occurs when searhing optimal OR!\nPlease adjust range [OR.low, OR.upp]!\n")
  } else {
    OR = exp(res$roo)
  }

  return(OR)
}


######

####
###########
# Lachin2008, Section 3.2.2, Formula (33) and (34)
###########


# continuous covariate
# power as a function of mean difference of continuous covariate and nD=1

# N - number of sets (total number of subjects = N * (nD + nH))
# power - power of the test
# delta - difference between mean of cases and mean of controls of continous covariate
# sigma - standard deviation of the continuous covariate
# nD = 1 - number of cases per set
# nH - number of controls per set
# R2 - coefficient of determination of the exposure variable and other covariates

powerConLogistic.con2 = function(N = NULL, 
                                power = 0.8, 
                                delta, 
                                sigma, 
                                nH,
                                R2 = 0,
                                alpha=0.05, 
                                nTests=1)
{
  # nD = 1
  alpha2 = alpha/2
  za = qnorm(1-alpha2/nTests)
  corePart = nH*delta^2*(1-R2)/((1+nH)*sigma^2)
  
  if(is.null(power) == TRUE && is.null(N) == FALSE)
  {
    # calculate power
    power = pnorm(sqrt(N*corePart) - za)
    
    return(power)
  } else if(is.null(power) == FALSE && is.null(N) == TRUE) {
    # calculate sample size
    zb = qnorm(power)
    numer = (zb + za)^2
    denom = corePart
    N = numer/denom
    
    return(N)
  } else {
    stop("one and only one of power and N should be null!")
  }
  
}

