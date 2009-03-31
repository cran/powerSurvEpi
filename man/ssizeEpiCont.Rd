\name{ssizeEpiCont}
\alias{ssizeEpiCont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Sample Size Calculation for Cox Proportional Hazards Regression with nonbinary covariates for Epidemiological Studies
}
\description{
Sample size calculation for Cox proportional hazards regression with nonbinary covariates for Epidemiological Studies. 
}
\usage{
ssizeEpiCont(formula, dat, X1, failureFlag, power, theta, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ a formula object relating the covariate of interest
    to other covariates to calculate the multiple correlation coefficient. The
    variables in formula must be in the data frame \code{dat}.
  }
  \item{dat}{ a \code{nPilot} by \code{p} data frame representing the pilot 
    data set, where \code{nPilot} is the number of subjects in the pilot study and
    the \code{p} (\eqn{>1}) columns contains the covariate of interest and other
    covariates.
  }
  \item{X1}{the covariate of interest.}
  \item{failureFlag}{a \code{nPilot} by 1 vector of indicators indicating if a subject is 
    failure (\code{failureFlag=1}) or alive (\code{failureFlag=0}). 
  }
  \item{power}{ postulated power. }
  \item{theta}{ postulated hazard ratio. }
  \item{alpha}{ type I error rate. }
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
where \eqn{\sigma^2=Var(X_1)}, \eqn{\psi} is the proportion of subjects died of
the disease of interest, and \eqn{\rho} is the multiple correlation coefficient
of the following linear regression:
\deqn{x_1=b_0+\boldsymbol{b}^T\boldsymbol{x}_2.} 
That is, \eqn{\rho^2=R^2}, where \eqn{R^2} is the proportion of variance
explained by the regression of \eqn{X_1} on the vector of covriates
\eqn{\boldsymbol{X}_2}.


\eqn{rho^2}, \eqn{\sigma^2}, and \eqn{\psi} will be estimated from a pilot study.
}
\value{
  \item{n}{the total number of subjects required.}
  \item{rho2}{square of the correlation between \eqn{X_1} and \eqn{X_2}.}
  \item{sigma2}{ variance of the covariate of interest. }
  \item{psi}{ proportion of subjects died of the disease of interest. }

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
\seealso{ \code{\link{ssizeEpiCont.default}} }
\examples{
  # generate a toy pilot data set
  set.seed(123456)
  X1 <- rnorm(100, mean = 0, sd = 0.3126)
  X2 <- sample(c(0, 1), 100, replace = TRUE)
  failureFlag <- sample(c(0, 1), 100, prob = c(0.25, 0.75), replace = TRUE)
  dat <- data.frame(X1 = X1, X2 = X2, failureFlag = failureFlag)

  ssizeEpiCont(formula = X1 ~ X2, dat = dat, X1 = X1, failureFlag = failureFlag, 
    power = 0.806, theta = exp(1), alpha = 0.05)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ design }% __ONLY ONE__ keyword per line
