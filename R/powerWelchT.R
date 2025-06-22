# calculate power of two-sided 2-sample t test with unequal variances and unequal sample sizes

powerWelchT <- function(
    n1, # sample size for group 1 
    n2, # sample size for group 2
    meanDiff, # mean difference between 2 groups
    sd1, # SD of group 1
    sd2, # SD of group 2
    alpha = 0.05 # type I error rate
)
{
  # the noncentrality parameter
  lambda <- abs(meanDiff) / sqrt((sd1^2 / n1) + (sd2^2 / n2))
  
  # Satterthwaite approximation of the degrees of freedom
  numer <- ((sd1^2 / n1) + (sd2^2 / n2))^2
  denom <- ((sd1^2 / n1)^2) / (n1 - 1) + ((sd2^2 / n2)^2) / (n2 - 1)
  df <- numer / denom
  
  # upper percentile of t distri with df
  cutoff <- qt(1 - alpha/2, df)
  
  power <- 1 - pt(cutoff, df, ncp = lambda) + pt(-cutoff, df, ncp = lambda)
  
  return(power)
}
