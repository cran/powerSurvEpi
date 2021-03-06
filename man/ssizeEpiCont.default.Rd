\name{ssizeEpiCont.default}
\alias{ssizeEpiCont.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Sample Size Calculation for Cox Proportional Hazards Regression with Nonbinary Covariates for Epidemiological Studies
}
\description{
Sample size calculation for Cox proportional hazards regression with nonbinary covariates for Epidemiological Studies. 
}
\usage{
ssizeEpiCont.default(power, 
		     theta, 
		     sigma2, 
		     psi, 
		     rho2, 
		     alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{numeric. postulated power. }
  \item{theta}{numeric. postulated hazard ratio. }
  \item{sigma2}{numeric. variance of the covariate of interest. }
  \item{psi}{numeric. proportion of subjects died of the disease of interest. }
  \item{rho2}{numeric. square of the multiple correlation coefficient between the covariate of interest and other covariates. }
  \item{alpha}{numeric. type I error rate. }
}
\details{
This is an implementation of the sample size calculation formula 
derived by Hsieh and Lavori (2000)
for the following Cox proportional hazards regression in the epidemiological
studies: 
\deqn{h(t|x_1, \boldsymbol{x}_2)=h_0(t)\exp(\beta_1 x_1+\boldsymbol{\beta}_2 
\boldsymbol{x}_2,}
where the covariate \eqn{X_1} is a nonbinary variable and 
\eqn{\boldsymbol{X}_2} is a vector of other covariates.
 
Suppose we want to check if  
the hazard ratio of the main effect \eqn{X_1=1} to \eqn{X_1=0} is equal to 
\eqn{1} or is equal to \eqn{\exp(\beta_1)=\theta}.
Given the type I error rate \eqn{\alpha} for a two-sided test, the total 
number of subjects required to achieve a sample size of \eqn{1-\beta} is
\deqn{n=\frac{\left(z_{1-\alpha/2}+z_{1-\beta}\right)^2}{
[\log(\theta)]^2 \sigma^2 \psi (1-\rho^2)
},}
where \eqn{z_{a}} is the \eqn{100 a}-th percentile of the standard normal distribution, \eqn{\sigma^2=Var(X_1)}, \eqn{\psi} is the proportion of subjects died of
the disease of interest, and \eqn{\rho} is the multiple correlation coefficient
of the following linear regression:
\deqn{x_1=b_0+\boldsymbol{b}^T\boldsymbol{x}_2.} 
That is, \eqn{\rho^2=R^2}, where \eqn{R^2} is the proportion of variance
explained by the regression of \eqn{X_1} on the vector of covriates
\eqn{\boldsymbol{X}_2}.
}
\value{
The total number of subjects required.
}
\references{ 
Hsieh F.Y. and Lavori P.W. (2000).
Sample-size calculation for the Cox proportional hazards regression model with nonbinary covariates.
\emph{Controlled Clinical Trials}. 21:552-560.

}
\note{ 
(1) Hsieh and Lavori (2000) assumed one-sided test, while this implementation assumed two-sided test.
(2) The formula can be used to calculate
ssize for a randomized trial study by setting \code{rho2=0}.
}
\seealso{ \code{\link{ssizeEpiCont}} }
\examples{
  # example in the EXAMPLE section (page 557) of Hsieh and Lavori (2000).
  # Hsieh and Lavori (2000) assumed one-sided test, 
  # while this implementation assumed two-sided test. 
  # Hence alpha=0.1 here (two-sided test) will correspond
  # to alpha=0.05 of one-sided test in Hsieh and Lavori's (2000) example.
  ssizeEpiCont.default(power = 0.806, 
		       theta = exp(1), 
		       sigma2 = 0.3126^2, 
                       psi = 0.738, 
		       rho2 = 0.1837, 
		       alpha = 0.1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ design }% __ONLY ONE__ keyword per line
