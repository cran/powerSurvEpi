\name{powerCT}
\alias{powerCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Power Calculation in the Analysis of Survival Data for Clinical Trials }
\description{
Power calculation for the Comparison of Survival Curves Between Two Groups under
 the Cox Proportional-Hazards Model for clinical trials. Some parameters will be estimated based on a pilot data set.
}
\usage{
powerCT(formula, 
	dat, 
	nE, 
	nC, 
	RR, 
	alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ A formula object, e.g. \code{Surv(time, status) ~ x}, where 
    \code{time} is a vector of survival/censoring time, \code{status} is a vector of 
    censoring indicator, \code{x} is the group indicator,
    which is a factor object in R and takes only two possible values (\code{C} for control group 
    and \code{E} for experimental group). See also the documentation of the function \code{survfit} in the 
    library \code{survival}.}
  \item{dat}{ a data frame representing the pilot data set and containing at least 3 columns: (1) survival/censoring time; (2) censoring indicator;
    (3) group indicator which is a factor object in R and takes only two possible values (\code{C} for control group 
    and \code{E} for experimental group). }
  \item{nE}{integer. number of participants in the experimental group. }
  \item{nC}{integer. number of participants in the control group. }
  \item{RR}{numeric. postulated hazard ratio.}
  \item{alpha}{numeric. type I error rate. }
}
\details{
This is an implementation of the power calculation method described in Section 14.12 (page 807)
of Rosner (2006). The method was proposed by Freedman (1982).

The movitation of this function is that some times we do not have information about \eqn{m} or \eqn{p_E}
and \eqn{p_C} available, but we have a pilot data set that can be used to estimate \eqn{p_E} and \eqn{p_C} hence
\eqn{m}, where  \eqn{m=n_E p_E + n_C p_C} is the expected total number of events over both groups, \eqn{n_E}
and \eqn{n_C} are numbers of participants in group E (experimental group) and group C (control group), respectively.
\eqn{p_E} is the probability of failure in group E (experimental group) over the maximum time period of the study (t years).  \eqn{p_C} is the probability of failure in group C (control group) over the maximum time period of the study 
(t years). 

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
E compared to group C \eqn{= n_E/n_C=k}, then the power of the test is
\deqn{power=\Phi(\sqrt{k*m}*|RR-1|/(k*RR+1)-z_{1-\alpha/2}),} 
where \deqn{m=n_E p_E+n_C p_C,} 
and \eqn{z_{1-\alpha/2}} 
is the \eqn{100 (1-\alpha/2)}-th percentile of
the standard normal distribution \eqn{N(0, 1)}, \eqn{\Phi} is the cumulative distribution function (CDF) 
of \eqn{N(0, 1)}.

\eqn{p_C} and \eqn{p_E} can be calculated from the following formulaes:
\deqn{p_C=\sum_{i=1}^{t}D_i, p_E=\sum_{i=1}^{t}E_i,}
where \eqn{D_i=\lambda_i A_i C_i}, \eqn{E_i=RR\lambda_i B_i C_i}, 
\eqn{A_i=\prod_{j=0}^{i-1}(1-\lambda_j)}, \eqn{B_i=\prod_{j=0}^{i-1}(1-RR\lambda_j)}, 
     \eqn{C_i=\prod_{j=0}^{i-1}(1-\delta_j)}.  And
\eqn{\lambda_i} is the probability of failure at time \code{i} among participants in the 
       control group, given that a participant has survived to time \eqn{i-1} and is not censored at time \eqn{i-1},
       i.e., the approximate hazard time \eqn{i} in the control group, \eqn{i=1,...,t}; 
     \eqn{RRlambda_i} is the probability of failure at time \code{i} among participants in the 
       experimental group, given that a participant has survived to time \eqn{i-1} and is not censored at time \eqn{i-1},
       i.e., the approximate hazard time \eqn{i} in the experimental group, \eqn{i=1,...,t}; 
     \eqn{delta} is the prbability that a participant is censored at time \eqn{i} given that he was 
       followed up to time \eqn{i} and has not failed, \eqn{i=0, 1, ..., t}, which is assumed the same in each group. 
}
\value{
  \item{mat.lambda }{a matrix with 9 columns and \code{nTimes+1} rows, where \code{nTimes} is the number of
     observed time points for the control group in the data set. The 9 columns are 
     (1) \code{time} - observed time point for the control group; 
     (2) \code{lambda};  
     (3) \code{RRlambda}; 
     (4) \code{delta}; 
     (5) \code{A}; 
     (6) \code{B}; 
     (7) \code{C}; 
     (8) \code{D}; 
     (9) \code{E}. 
     Please refer to the Details section for the definitions of elements
of these quantities. See also Table 14.24 on page 809 of Rosner (2006). 
  }
  \item{mat.event }{a matrix with 5 columns and \code{nTimes+1} rows, where \code{nTimes} is the number of
     observed time points for control group in the data set. The 5 columns are 
     (1) \code{time} - observed time point for the control group; 
     (2) \code{nEvent.C} - number of events in the control group at each time point;
     (3) \code{nCensored.C} - number of censorings in the control group at each time point;
     (4) \code{nSurvive.C} - number of alived in the control group at each time point;
     (5) \code{nRisk.C} - number of participants at risk in the control group at each time point.
     Please refer to Table 14.12 on page 787 of Rosner (2006). 
  }
  \item{pC}{ estimated probability of failure in group C (control group) over the maximum time period of the study (t years). }
  \item{pE}{ estimated probability of failure in group E (experimental group) over the maximum time period of the study (t years). }
  \item{power}{the power of the test.}
}
\references{ 
  Freedman, L.S. (1982). 
  Tables of the number of patients required in clinical trials using the log-rank test. 
  \emph{Statistics in Medicine}. 1: 121-129

  Rosner B. (2006).
  \emph{Fundamentals of Biostatistics}. (6-th edition). Thomson Brooks/Cole.
}
\note{ 
(1) The estimates of \eqn{RRlambda_i=RR*\lambda_i}. That is, \code{RRlambda} is not directly estimated based on data
from the experimental group;
(2) The power formula assumes that the central-limit theorem is valid and hence is appropriate for large samples.
}
\seealso{ \code{\link{powerCT.default0}}, 
\code{\link{powerCT.default}} }
\examples{
  # Example 14.42 in Rosner B. Fundamentals of Biostatistics. 
  # (6-th edition). (2006) page 809

  library(survival)

  data(Oph)
  res <- powerCT(formula = Surv(times, status) ~ group, 
		 dat = Oph, 
                 nE = 200, 
		 nC = 200, 
		 RR = 0.7, 
		 alpha = 0.05)

  # Table 14.24 on page 809 of Rosner (2006)
  print(round(res$mat.lambda, 4))

  # Table 14.12 on page 787 of Rosner (2006)
  print(round(res$mat.event, 4))

  # the power
  print(round(res$power, 2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ design }% __ONLY ONE__ keyword per line
