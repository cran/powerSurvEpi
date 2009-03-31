\name{powerEpiInt2}
\alias{powerEpiInt2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Power Calculation Testing Interaction Effect for Cox Proportional Hazards Regression with two covariates for Epidemiological Studies (Both covariates should be binary)
}
\description{
Power calculation testing interaction effect for Cox proportional hazards regression with two covariates for Epidemiological Studies. Both covariates should be  binary variables. The formula takes
into account the correlation between the two covariates.
}
\usage{
powerEpiInt2(n, theta, psi, mya, myb, myc, myd, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ total number of subjects. }
  \item{theta}{ postulated hazard ratio. }
  \item{psi}{ proportion of subjects died of the disease of interest. }
  \item{mya}{ number of subjects taking values \eqn{X_1=0} and 
\eqn{X_2=0} obtained from a pilot study.}
  \item{myb}{ number of subjects taking values \eqn{X_1=0} and 
\eqn{X_2=1} obtained from a pilot study.}
  \item{myc}{ number of subjects taking values \eqn{X_1=1} and 
\eqn{X_2=0} obtained from a pilot study.}
  \item{myd}{ proportion of subjects taking values \eqn{X_1=1} and 
\eqn{X_2=1} obtained from a pilot study.}
  \item{alpha}{ type I error rate. }
}
\details{
This is an implementation of the power calculation formula 
derived by Schmoor et al. (2000) for
the following Cox proportional hazards regression in the epidemiological
studies: 
\deqn{h(t|x_1, x_2)=h_0(t)\exp(\beta_1 x_1+\beta_2 x_2 + \gamma (x_1 x_2)),}
where both covariates \eqn{X_1} and \eqn{X_2} are binary variables.

Suppose we want to check if  
the hazard ratio of the interaction effect \eqn{X_1 X_2=1} to \eqn{X_1 X_2=0} is equal to \eqn{1} 
or is equal to \eqn{\exp(\gamma)=\theta}.
Given the type I error rate \eqn{\alpha} for a two-sided test, the power
required to detect a hazard ratio as small as \eqn{\exp(\gamma)=\theta} is
\deqn{power=\Phi\left(-z_{1-\alpha/2}+\sqrt{\frac{n}{G}[\log(\theta)]^2 p (1-p) \psi (1-\rho^2)}\right),}
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
\eqn{p=Pr(X_1=1)=(myc+myd)/n_{obs}},
\eqn{q=Pr(X_2=1)=(myb+myd)/n_{obs}},
\eqn{n_{obs}=mya+myb+myc+myd}.


\eqn{p_{00}=Pr(X_1=0,\mbox{and}, X_2=0)}, 
\eqn{p_{01}=Pr(X_1=0,\mbox{and}, X_2=1)}, 
\eqn{p_{10}=Pr(X_1=1,\mbox{and}, X_2=0)}, 
\eqn{p_{11}=Pr(X_1=1,\mbox{and}, X_2=1)}. 
}
\value{
The power of the test.
}
\references{ 
Schmoor C., Sauerbrei W., and Schumacher M. (2000).
Sample size considerations for the evaluation of prognostic factors in survival analysis.
\emph{Statistics in Medicine}. 19:441-452.

}
\seealso{ \code{\link{powerEpiInt.default0}}, \code{\link{powerEpiInt.default1}} }
\examples{
  # Example at the end of Section 4 of Schmoor et al. (2000).
  # mya, myb, myc, and myd are obtained from Table III on page 448
  # of Schmoor et al. (2000).
  powerEpiInt2(n = 184, theta = 3, psi = 139 / 184,
    mya = 50, myb = 21, myc = 78, myd = 35, alpha = 0.05)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ design }% __ONLY ONE__ keyword per line
