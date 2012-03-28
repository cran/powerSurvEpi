\name{power.stratify}
\alias{power.stratify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power calculation for survival analysis with binary predictor and exponential survival function
}
\description{
Power calculation for survival analysis with binary predictor and exponential survival function.
}
\usage{
power.stratify(
    n, 
    timeUnit, 
    gVec, 
    PVec, 
    HR, 
    lambda0Vec, 
    power.ini = 0.8, 
    power.low = 0.001, 
    power.upp = 0.999, 
    alpha = 0.05, 
    verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Scalar. Sample size.
}
  \item{timeUnit}{
Scalar. Total study length.
}
  \item{gVec}{
m by 1 vector. The s-th element is the proportion of the total sample size for 
the s-th stratum, where m is the number of strata.
}
  \item{PVec}{
m by 1 vector. The s-th element is the proportion of subjects in treatment group 1 for the s-th stratum, where m is the number of strata.
}
  \item{HR}{
Scalar. Hazard ratio (Ratio of the hazard for treatment group 1 to the hazard for treatment group 0, i.e. reference group).
}
  \item{lambda0Vec}{
m by 1 vector. The s-th element is the hazard for treatment group 0 (i.e., reference group) in the s-th stratum.
}
  \item{power.ini}{
Scalar. Initial power estimate.
}
  \item{power.low}{
Scalar. Lower bound for power.
}
  \item{power.upp}{
Scalar. Upper bound for power.
}
  \item{alpha}{
Scalar. Type I error rate.
}
  \item{verbose}{
Logical. Indicating if intermediate results will be output or not.
}
}
\details{
We assume (1) there is only one predictor and no covariates in the survival model 
(exponential survival function); (2) there are \code{m} strata; (3) the predictor \code{x}
is a binary variable indicating treatment group 1 (\eqn{x=1}) or treatment group 0
(\eqn{x=0}); (3) the treatment effect is constant over time (proportional hazards);
(4) the hazard ratio is the same in all strata, and (5) the data will be analyzed by 
the stratified log rank test.

The sample size formula is Formula (1) on page 801 of Palta M and Amini SB (1985):
\deqn{
n=(Z_{\alpha}+Z_{\beta})^2/\mu^2
}{
n=(Z_{\alpha}+Z_{\beta})^2/\mu^2,
} 
where \eqn{\alpha} is the Type I error rate,
\eqn{\beta} is the Type II error rate (power\eqn{=1-\beta}),
\eqn{Z_{\alpha}} is the \eqn{100(1-\alpha)} percentile of standard normal distribution, and
\deqn{
\mu=\log(\delta)\sqrt{ \sum_{s=1}^{m} g_s P_s (1 - P_s) V_s }
}{
\mu=\log(\delta)\sqrt{ \sum_{s=1}^{m} g_s P_s (1 - P_s) V_s },
}
and
\deqn{
V_s=P_s\left[1-\frac{1}{\lambda_{1s}} \left\{
\exp\left[-\lambda_{1s}(T-1)\right]
-\exp(-\lambda_{1s}T)
\right\}
\right]
+(1-P_s)\left[
1-\frac{1}{\lambda_{2s}}
\left\{
\exp\left[-\lambda_{2s}(T-1)\right]
-\exp(-\lambda_{2s}T
\right\}
\right]
}{
V_s=P_s\left[1-\frac{1}{\lambda_{1s}} \left\{
\exp\left[-\lambda_{1s}(T-1)\right]
-\exp(-\lambda_{1s}T)
\right\}
\right]
+(1-P_s)\left[
1-\frac{1}{\lambda_{2s}}
\left\{
\exp\left[-\lambda_{2s}(T-1)\right]
-\exp(-\lambda_{2s}T
\right\}
\right].
}
In the above formulas, \eqn{m} is the number of strata,
\eqn{T} is the total study length, \eqn{\delta} is the hazard ratio,
\eqn{g_s} is the proportion of the total sample size in stratum \eqn{s},
\eqn{P_s} is the proportion of stratum \eqn{s}, which is in treatment group 1,
and \eqn{\lambda_{is}} is the hazard for the \eqn{i}-th treatment group in 
stratum \eqn{s}.


}
\value{
A list of 2 elments.
\item{power }{Estimated power}
\item{res.optim }{Object returned by funciton \code{optim}. We used numerical optimization method to calculate power based on sample size calculation formula.}
}
\references{
Palta M and Amini SB. (1985). 
Consideration of covariates and stratification in sample size determination for survival time studies. 
\emph{Journal of Chronic Diseases}. 38(9):801-809.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ssize.stratify}}
}
\examples{

# example on page 803 of Palta M and Amini SB. (1985). 
res.power <- power.stratify(
  n = 146, 
  timeUnit = 1.25, 
  gVec = c(0.5, 0.5),
  PVec = c(0.5, 0.5), 
  HR = 1 / 1.91, 
  lambda0Vec = c(2.303, 1.139),
  power.ini = 0.8, 
  power.low = 0.001, 
  power.upp = 0.999,
  alpha = 0.05, 
  verbose = TRUE
  )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
