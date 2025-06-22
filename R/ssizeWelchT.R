
ttFuncWelchT <- function(n1,
                         ratioN2toN1,
                         meanDiff,
                         sd1,
                         sd2,
                         power,
                         alpha,
                         minN1 = 3)
{

  n2 <-  ceiling(n1*ratioN2toN1)

  ttpower <- powerWelchT(
    n1 = n1,  
    n2 = n2, 
    meanDiff = meanDiff, 
    sd1 = sd1, 
    sd2 = sd2, 
    alpha = alpha
  )
  
  diff = (ttpower - power)^2
  
  
  return(diff)
}

ssizeWelchT <- function(
    ratioN2toN1, # ratio of sample size for group 2 to sample size for group 1
    meanDiff, # mean difference between 2 groups
    sd1, # SD of group 1
    sd2, # SD of group 2
    power = 0.8, # power
    alpha = 0.05, # type I error rate
    minN1 = 3 # minimu possible sample size for group 1
)
{
  
  # Use optim with lower bound minN1
  result <- optim(par = minN1+1, 
                  fn = ttFuncWelchT, 
                  method = "L-BFGS-B", 
                  lower = minN1,
                  ratioN2toN1 = ratioN2toN1,
                  meanDiff = meanDiff,
                  sd1 = sd1,
                  sd2 = sd2,
                  power = power,
                  alpha = alpha,
                  minN1 = minN1)
  

  n1 = ceiling(result$par)    # Optimal x
  n2 = ceiling(n1 * ratioN2toN1)

  ttpower <- powerWelchT(
    n1 = n1,  
    n2 = n2, 
    meanDiff = meanDiff, 
    sd1 = sd1, 
    sd2 = sd2, 
    alpha = alpha
  )
  
  if(ttpower < power)
  {
    n1 = n1 + 1
  }
  
  
  res = list(n1 = n1, n2 = n2)
  
  return(res)
}
