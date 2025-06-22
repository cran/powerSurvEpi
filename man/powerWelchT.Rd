\name{powerWelchT}
\alias{powerWelchT}
\title{
Power of Two-Sided Two Sample T Test With Unequal Variances And Unequal Sample Sizes
}
\description{
Power of two-sided 2 sample t test with unequal variances and unequal sample sizes.
}
\usage{
powerWelchT(
  n1, 
  n2, 
  meanDiff, 
  sd1, 
  sd2, 
  alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n1}{
sample size for group 1
}
  \item{n2}{
sample size for group 2
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
  \item{alpha}{
Type I error rate
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
power
}

\examples{
powerWelchT(
    n1 = 64, # sample size for group 1 
    n2 = 30, # sample size for group 2
    meanDiff = 1, # mean difference between 2 groups
    sd1 = 2, # SD of group 1
    sd2 = 1, # SD of group 2
    alpha = 0.05 # type I error rate
)
# 0.8918191

}
\keyword{ survival }
\keyword{ design }