\name{ssizeCT.default}
\alias{ssizeCT.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sample Size Calculation in the Analysis of Survival Data for Clinical Trials }
\description{
Sample size calculation for the Comparison of Survival Curves Between Two Groups under
 the Cox Proportional-Hazards Model for clinical trials. 
}
\usage{
ssizeCT.default(power, k, pE, pC, RR, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{ power to detect the magnitude of the hazard ratio as small as that specified by \code{RR}.}
  \item{k}{ ratio of participants in group E (experimental group) compared to group C (control group).}
  \item{pE}{ probability of failure in group E (experimental group) over the maximum time period of the study (t years). }
  \item{pC}{ probability of failure in group C (control group) over the maximum time period of the study (t years). }
  \item{RR}{ postulated hazard ratio. }
  \item{alpha}{ type I error rate. }
}
\details{
This is an implementation of the sample size calculation method described in Section 14.12 (page 807)
of Rosner (2006). The method was proposed by Freedman (1982).

Suppose we want to compare the survival curves between an experimental group (\eqn{E}) and 
a control group (\eqn{C}) in a clinical trial with a maximum follow-up of \eqn{t} years. 
The Cox proportional hazards regression model is assumed to have the form:
\deqn{h(t|X_1)=h_0(t)\exp(\beta_1 X_1).}
Let \eqn{n_E} be the number of participants in the \eqn{E} group
and \eqn{n_C} be the number of participants in the \eqn{C} group.  
We wish to test the hypothesis \eqn{H0: RR=1} versus \eqn{H1: RR} not equal to 1, 
where \eqn{RR=\exp(\beta_1)=}underlying hazard ratio
for the \eqn{E} group versus the \eqn{C} group. Let \eqn{RR} be the postulated hazard ratio,
\eqn{\alpha} be the significance level. Assume that the test is a two-sided test.
If the ratio of participants in group 
E compared to group C \eqn{= n_E/n_C=k}, then the number of participants needed in each group to 
achieve a power of \eqn{1-\beta} is
\deqn{n_E=\frac{m k}{k p_E + p_C}, n_C=\frac{m}{k p_E + p_C}} 
where \deqn{m=\frac{1}{k}\left(\frac{k RR + 1}{RR - 1}\right)^2\left(
z_{1-\alpha/2}+z_{1-\beta}
\right)^2,} 
and \eqn{z_{1-\alpha/2}} 
is the \eqn{100 (1-\alpha/2)} percentile of
the standard normal distribution \eqn{N(0, 1)}.

}
\value{
  A two-element vector. The first element is \eqn{n_E} and the second 
    element is \eqn{n_C}.
}
\references{ 
  Freedman, L.S. (1982). 
  Tables of the number of patients required in clinical trials using the log-rank test. 
  \emph{Statistics in Medicine}. 1: 121-129

  Rosner B. (2006).
  \emph{Fundamentals of Biostatistics}. (6-th edition). Thomson Brooks/Cole.
}
\note{ 
(1) The sample size formula assumes that the central-limit theorem is valid and hence is appropriate for large samples.
(2) \eqn{n_E} and \eqn{n_C} will be rounded up to integers.
}
\seealso{ \code{\link{ssizeCT}} }
\examples{
  # Example 14.42 in Rosner B. Fundamentals of Biostatistics. 
  # (6-th edition). (2006) page 809
  ssizeCT.default(power = 0.8, k = 1, pE = 0.3707, pC = 0.4890, 
    RR = 0.7, alpha = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ design }% __ONLY ONE__ keyword per line
