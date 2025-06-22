\name{ssizeWelchT}
\alias{ssizeWelchT}
\title{
Sample Size Calculation For Two-Sided Two Sample T Test With Unequal Variances And Unequal Sample Sizes
}
\description{
Sample size calculation for two-sided two sample t test with unequal variances and unequal sample sizes
}
\usage{
ssizeWelchT(
  ratioN2toN1, 
  meanDiff, 
  sd1, 
  sd2, 
  power = 0.8, 
  alpha = 0.05, 
  minN1 = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ratioN2toN1}{
numeric. The ratio of sample size for group 2 to sample size for group 1
}
  \item{meanDiff}{
mean difference between 2 groups
}
  \item{sd1}{
standard deviation of group 1
}
  \item{sd2}{
standard deviation of group 2
}
  \item{power}{
power
}
  \item{alpha}{
Type I error rate
}
  \item{minN1}{
minimum sample size for group 1
}
}
\details{

The power formula is
\deqn{
power = Pr\left(|T| > t_{1-\alpha/2, \nu} | T \sim t_{\nu, \lambda}\right),
}
where \eqn{\lambda} is the noncentrality parameter of the t distribution with degree of freedom 
\eqn{\nu}. \eqn{t_{1-\alpha/2, \nu}} is the upper \eqn{100\alpha/2} percentile of the t distribution
with degree of freedom \eqn{\nu}. \eqn{\alpha} is the significance level.
The noncentrality parameter \eqn{\lambda} is defined as
\deqn{
\lambda = \frac{|\mu_1 - \mu_2|}{\sqrt{\frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2}}}.
}

The degree \eqn{\nu} of freedom is the Satterthwaite approximation and is defined as
\deqn{
\nu = \frac{\left(\frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2}\right)^2}{
\frac{\left(\frac{\sigma_1^2}{n_1}\right)^2}{n_1-1}
+
\frac{\left(\frac{\sigma_2^2}{n_2}\right)^2}{n_2-1}
}
}

}
\value{
A list with 2 elements
  \item{n1}{sample size for group 1}
  \item{n2}{sample size for group 2}
}

\examples{
ssizeWelchT(
    ratioN2toN1=30/64, # ratio of sample size for group 2 to sample size for group 1
    meanDiff = 1, # mean difference between 2 groups
    sd1 = 2, # SD of group 1
    sd2 = 1, # SD of group 2
    power = 0.8918191, # power
    alpha = 0.05, # type I error rate
    minN1 = 3 # minimu possible sample size for group 1
)
# n1 = 64 and n2 = 30
}
\keyword{ survival }
\keyword{ design }
