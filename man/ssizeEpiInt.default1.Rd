\name{ssizeEpiInt.default1}
\alias{ssizeEpiInt.default1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Sample Size Calculation Testing Interaction Effect for Cox Proportional Hazards Regression
}
\description{
Sample size calculation testing interaction effect for Cox proportional hazards regression with two covariates for Epidemiological Studies. Both covariates should be  binary variables. The formula takes
into account the correlation between the two covariates.
}
\usage{
ssizeEpiInt.default1(power, 
		     theta, 
		     psi, 
		     p00, 
		     p01, 
		     p10, 
		     p11, 
		     alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{numeric. postulated power. }
  \item{theta}{numeric. postulated hazard ratio. }
  \item{psi}{numeric. proportion of subjects died of the disease of interest. }
  \item{p00}{numeric. proportion of subjects taking values \eqn{X_1=0} and 
\eqn{X_2=0}, i.e., \eqn{p_{00}=Pr(X_1=0,\mbox{and}, X_2=0)}.}
  \item{p01}{numeric. proportion of subjects taking values \eqn{X_1=0} and 
\eqn{X_2=1}, i.e., \eqn{p_{01}=Pr(X_1=0,\mbox{and}, X_2=1)}.}
  \item{p10}{numeric. proportion of subjects taking values \eqn{X_1=1} and 
\eqn{X_2=0}, i.e., \eqn{p_{10}=Pr(X_1=1,\mbox{and}, X_2=0)}.}
  \item{p11}{numeric. proportion of subjects taking values \eqn{X_1=1} and 
\eqn{X_2=1}, i.e., \eqn{p_{11}=Pr(X_1=1,\mbox{and}, X_2=1)}.}
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
number of subjects required to achieve a power of \eqn{1-\beta} is
\deqn{n=\frac{\left(z_{1-\alpha/2}+z_{1-\beta}\right)^2\delta}{[\log(\theta)]^2 \psi},}
where \eqn{z_{a}} is the \eqn{100 a}-th percentile of the standard normal distribution, 
\eqn{\psi} is the proportion of subjects died of
the disease of interest, 
\deqn{\delta=\frac{1}{p_{00}}+\frac{1}{p_{01}}+\frac{1}{p_{10}}
+\frac{1}{p_{11}},}
and
\eqn{p_{00}=Pr(X_1=0,\mbox{and}, X_2=0)}, 
\eqn{p_{01}=Pr(X_1=0,\mbox{and}, X_2=1)}, 
\eqn{p_{10}=Pr(X_1=1,\mbox{and}, X_2=0)}, 
\eqn{p_{11}=Pr(X_1=1,\mbox{and}, X_2=1)}. 
}
\value{
The ssize of the test.
}
\references{ 
Schmoor C., Sauerbrei W., and Schumacher M. (2000).
Sample size considerations for the evaluation of prognostic factors in survival analysis.
\emph{Statistics in Medicine}. 19:441-452.

}
\seealso{ \code{\link{ssizeEpiInt.default0}}, \code{\link{ssizeEpiInt2}} }
\examples{
  # Example at the end of Section 4 of Schmoor et al. (2000).
  # p00, p01, p10, and p11 are calculated based on Table III on page 448
  # of Schmoor et al. (2000).
  ssizeEpiInt.default1(power = 0.8227, 
		       theta = 3, 
		       psi = 139 / 184,
                       p00 = 50/184, 
		       p01 = 21 / 184, 
		       p10 = 78 / 184, 
		       p11 = 35 / 184,
                       alpha = 0.05)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ design }% __ONLY ONE__ keyword per line
