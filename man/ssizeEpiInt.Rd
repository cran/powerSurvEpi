\name{ssizeEpiInt}
\alias{ssizeEpiInt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Sample Size Calculation Testing Interaction Effect for Cox Proportional Hazards Regression with two covariates for Epidemiological Studies (Both covariates should be binary)
}
\description{
Sample size calculation testing interaction effect for Cox proportional hazards regression with two covariates for Epidemiological Studies. Both covariates should be  binary variables. The formula takes
into account the correlation between the two covariates.
}
\usage{
ssizeEpiInt(X1, X2, failureFlag, power, theta, alpha = 0.05)
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
    for the \code{nPilot} subjects in the pilot study. 
    \code{X2} should be binary and take only two possible values: zero and one.
  }
  \item{failureFlag}{a \code{nPilot} by 1 vector of indicators indicating if a subject is 
    failure (\code{failureFlag=1}) or alive (\code{failureFlag=0}). 
  }
  \item{power}{ postulated power. }
  \item{theta}{ postulated hazard ratio. }
  \item{alpha}{ type I error rate. }
}
\details{
This is an implementation of the sample size calculation formula 
derived by Schmoor et al. (2000) for
the following Cox proportional hazards regression in the epidemoilogical
studies: 
\deqn{h(t|x_1, x_2)=h_0(t)\exp(\beta_1 x_1+\beta_2 x_2 + \gamma (x_1 x_2)),}
where both covariates \eqn{X_1} and \eqn{X_2} are binary variables.

Suppose we want to check if  
the hazard ratio of the interaction effect \eqn{X_1 X_2=1} to \eqn{X_1 X_2=0} is equal to \eqn{1} 
or is equal to \eqn{\exp(\gamma)=\theta}.
Given the type I error rate \eqn{\alpha} for a two-sided test, the total
number of subjects required to achieve the desired power \eqn{1-\beta} is:
\deqn{n=\frac{\left(z_{1-\alpha/2}+z_{1-\beta}\right)^2 G}{
[\log(\theta)]^2 \psi (1-p) p (1-\rho^2)
},}
where \eqn{\psi} is the proportion of subjects died of
the disease of interest, and 
\deqn{\rho=corr(X_1, X_2)=(p_1-p_0)\times\sqrt{\frac{q(1-q)}{p(1-p)}},} 
and
\eqn{p=Pr(X_1=1)}, \eqn{q=Pr(X_2=1)}, \eqn{p_0=Pr(X_1=1|X_2=0)},
and \eqn{p_1=Pr(X_1=1 | X_2=1)}, and
\deqn{G=\frac{[(1-q)(1-p_0)p_0+q(1-p_1)p_1]^2}{(1-q)q (1-p_0)p_0 (1-p_1) p_1},}
and
\eqn{p0=Pr(X_1=1 | X_2=0)=myc/(mya+myc)},
\eqn{p1=Pr(X_1=1 | X_2=1)=myd/(myb+myd)},
\eqn{p=Pr(X_1=1)=(myc+myd)/n},
\eqn{q=Pr(X_2=1)=(myb+myd)/n},
\eqn{n=mya+myb+myc+myd}.


\eqn{p_{00}=Pr(X_1=0,\mbox{and}, X_2=0)}, 
\eqn{p_{01}=Pr(X_1=0,\mbox{and}, X_2=1)}, 
\eqn{p_{10}=Pr(X_1=1,\mbox{and}, X_2=0)}, 
\eqn{p_{11}=Pr(X_1=1,\mbox{and}, X_2=1)}. 

\eqn{p_{00}}, \eqn{p_{01}}, \eqn{p_{10}}, \eqn{p_{11}}, and \eqn{\psi} will be
estimated from the pilot data.
}
\value{
  \item{n}{the total number of subjects required.}
  \item{p}{estimated \eqn{Pr(X_1=1)}}
  \item{q}{estimated \eqn{Pr(X_2=1)}}
  \item{p0}{estimated \eqn{Pr(X_1=1 | X_2=0)}}
  \item{p1}{estimated \eqn{Pr(X_1=1 | X_2=1)}}
  \item{rho2}{square of the estimated \eqn{corr(X_1, X_2)}}
  \item{G}{a factor adjusting the sample size. The sample size needed to
detect an effect of a prognostic factor with given error probabilities has
to be multiplied by the factor \code{G} when an interaction of the
same magnitude is to be detected.}
  \item{mya}{estimated number of subjects taking values \eqn{X_1=0} and 
    \eqn{X_2=0}.}
  \item{myb}{estimated number of subjects taking values \eqn{X_1=0} and 
    \eqn{X_2=1}.}
  \item{myc}{estimated number of subjects taking values \eqn{X_1=1} and 
    \eqn{X_2=0}.}
  \item{myd}{estimated number of subjects taking values \eqn{X_1=1} and 
    \eqn{X_2=1}.}
  \item{psi}{ proportion of subjects died of the disease of interest. }
}
\references{ 
Schmoor C., Sauerbrei W., and Schumacher M. (2000).
Sample size considerations for the evaluation of prognostic factors in survival analysis.
\emph{Statistics in Medicine}. 19:441-452.

}
\seealso{ \code{\link{ssizeEpiInt.default0}}, \code{\link{ssizeEpiInt2}} }
\examples{
  # generate a toy pilot data set
  X1 <- c(rep(1, 39), rep(0, 61))
  set.seed(123456)
  X2 <- sample(c(0, 1), 100, replace = TRUE)
  failureFlag <- sample(c(0, 1), 100, prob = c(0.25, 0.75), replace = TRUE)

  ssizeEpiInt(X1, X2, failureFlag, power = 0.88, theta = 3, alpha = 0.05)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ design }% __ONLY ONE__ keyword per line
