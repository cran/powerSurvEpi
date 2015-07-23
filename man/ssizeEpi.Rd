\name{ssizeEpi}
\alias{ssizeEpi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sample Size Calculation for Cox Proportional Hazards Regression}
\description{
Sample size calculation for Cox proportional hazards regression with two covariates for Epidemiological Studies. The covariate of interest should be a binary variable. The other covariate can be either binary or non-binary. The formula takes
into account competing risks and the correlation between the two covariates.
}
\usage{
ssizeEpi(X1, X2, failureFlag, power, theta, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X1}{ a \code{nPilot} by 1 vector, where \code{nPilot} is the number of subjects 
    in the pilot data set. This vector records the values of the covariate of 
    interest for the \code{nPilot} subjects in the pilot study. \code{X1} should 
    be binary and take only two possible values: zero and one.
  }
  \item{X2}{ a \code{nPilot} by 1 vector, where \code{nPilot} is the number of subjects
    in the pilot study. This vector records the values of the second covariate 
    for the \code{nPilot} subjects in the pilot study. \code{X2} can be binary or 
    non-binary.
  }
  \item{failureFlag}{a \code{nPilot} by 1 vector of indicators indicating if a subject is 
    failure (\code{failureFlag=1}) or alive (\code{failureFlag=0}). 
  }
  \item{power}{ postulated power. }
  \item{theta}{ postulated hazard ratio. }
  \item{alpha}{ type I error rate. }
}
\details{
This is an implementation of the sample size formula 
derived by Latouche et al. (2004)
for the following Cox proportional hazards regression in the epidemiological
studies:
\deqn{h(t|x_1, x_2)=h_0(t)\exp(\beta_1 x_1+\beta_2 x_2),}
where the covariate \eqn{X_1} is of our interest. The covariate \eqn{X_1} has to be
a binary variable taking two possible values: zero and one, while the 
covariate \eqn{X_2} can be binary or continuous.

Suppose we want to check if the hazard of \eqn{X_1=1} is equal to 
the hazard of \eqn{X_1=0} or not. Equivalently, we want to check if 
the hazard ratio of \eqn{X_1=1} to \eqn{X_1=0} is equal to \eqn{1} 
or is equal to \eqn{\exp(\beta_1)=\theta}.
Given the type I error rate \eqn{\alpha} for a two-sided test, the total 
number of subjects required to achieve a power of \eqn{1-\beta} is
\deqn{n=\frac{\left(z_{1-\alpha/2}+z_{1-\beta}\right)^2}{
[\log(\theta)]^2 p (1-p) \psi (1-\rho^2)},}
where \eqn{\psi} is the proportion of subjects died of
the disease of interest, and 
\deqn{\rho=corr(X_1, X_2)=(p_1-p_0)\times\sqrt{\frac{q(1-q)}{p(1-p)}},} 
and
\eqn{p=Pr(X_1=1)}, \eqn{q=Pr(X_2=1)}, \eqn{p_0=Pr(X_1=1|X_2=0)},
and \eqn{p_1=Pr(X_1=1 | X_2=1)}.

\eqn{p}, \eqn{\rho^2}, and \eqn{\psi} will be estimated from a pilot study.

}
\value{
  \item{n}{the total number of subjects required.}
  \item{p}{the proportion that \eqn{X_1} takes value one.}
  \item{rho2}{square of the correlation between \eqn{X_1} and \eqn{X_2}.}
  \item{psi}{ proportion of subjects died of the disease of interest. }

}
\references{ 
Schoenfeld DA. (1983).
Sample-size formula for the proportional-hazards regression model.
\emph{Biometrics}. 39:499-503.

Latouche A., Porcher R. and Chevret S. (2004).
Sample size formula for proportional hazards modelling of competing risks.
\emph{Statistics in Medicine}. 23:3263-3274.
}
\note{ 
(1) The calculated sample size will be round up to an integer.  

(2) The formula can be used to calculate
sample size required for a randomized trial study by setting \code{rho2=0}.

(3) When \code{rho2=0}, the formula derived by Latouche et al. (2004)
looks the same as that derived by Schoenfeld (1983). Latouche et al. (2004) pointed out that in this situation, the interpretations are different hence
the two formulae are actually different. In Latouched et al. (2004), the
hazard ratio \eqn{\exp(\beta_1)=\theta} measures the difference of effect of a covariate 
at two different levels on  cause-specific hazard for a particular failure, 
while in Schoenfeld (1983), the hazard ratio \eqn{\theta} measures 
the difference of effect on subdistribution hazard.
}
\seealso{ \code{\link{ssizeEpi.default}} }
\examples{
  # generate a toy pilot data set
  X1 <- c(rep(1, 39), rep(0, 61))
  set.seed(123456)
  X2 <- sample(c(0, 1), 100, replace = TRUE)
  failureFlag <- sample(c(0, 1), 100, prob = c(0.5, 0.5), replace = TRUE)

  ssizeEpi(X1 = X1, X2 = X2, failureFlag = failureFlag, 
    power = 0.80, theta = 2, alpha = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ design }% __ONLY ONE__ keyword per line
