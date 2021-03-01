\name{powerConLogistic.bin}
\alias{powerConLogistic.bin}
\title{
Sample Size Calculation for Conditional Logistic Regression with Binary Covariate
}
\description{
Sample Size Calculation for Conditional Logistic Regression with Binary Covariate, such as matched logistic regression or nested case-control study.
}
\usage{
powerConLogistic.bin(
  N = NULL, 
  power = 0.8, 
  OR, 
  pE, 
  nD, 
  nH, 
  R2 = 0, 
  alpha = 0.05, 
  nTests = 1,
  OR.low = 1.01,
  OR.upp = 100
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
integer. Number of sets. Each set contains \code{nD} cases and \code{nH} controls.  
}
  \item{power}{
numeric. Power of the test for if the exposure variable is associated with the risk of diseases
}
  \item{OR}{
numeric. Odds ratio \eqn{=exp(\theta)}, where \eqn{\theta} is the 
regression coefficient of the exposure variable.
}
  \item{pE}{
numeric. Population prevalence of exposure.
}
  \item{nD}{
integer. Number of cases per set.
}
  \item{nH}{
integer. Number of controls per set.
}
  \item{R2}{
numeric. Coefficient of determination of the exposure variable and other covariates
}
  \item{alpha}{
numeric. family-wise type I error rate.
}
  \item{nTests}{
integer. Number of tests.
}
  \item{OR.low}{
numeric. Lower bound of odds ratio. Only used when \code{OR=NULL} and \code{power} and \code{N} are not equal to NULL.
}
  \item{OR.upp}{
numeric. Upper bound of odds ratio. Only used when \code{OR=NULL} and \code{power} and \code{N} are not equal to NULL.
}
}
\details{
The power and sample size calculation formulas are provided by Lachin (2008, Section 3.3, Formula (38))

\deqn{power = \Phi\left( \sqrt{N c} - z_{\alpha/(2 nTests)}\right)}
and
\deqn{N = (z_{power} + z_{\alpha/(2 nTests)})^2/ c}
where \eqn{\Phi} is the cumulative distribution function of the 
standard normal distribution \eqn{N(0, 1)}, \eqn{z_{a}} is the upper \eqn{100 a}-th
percentile of \eqn{N(0, 1)},
\deqn{c = \theta^2 pE (1-pE) (1-R^2)nD*nH/(nD+nH)}
and \eqn{R^2} is the coefficient of determination for linear regression
linking the exposure with other covariates.
}
\value{
If the inputs \code{is.null(N) = TRUE} and \code{is.null(power) = FALSE}, 
then the function returns the number \code{N} of sets.

If the inputs \code{is.null(N) = FALSE} and \code{is.null(power) = TRUE}, 
then the function returns the power. 

Otherwise, an error message is output.
}
\references{
Lachin, JM
Sample Size Evaluation for a Multiply Matched Case-Control Study Using the Score Test From a Conditional Logistic (Discrete Cox PH) Regression Model.
Stat Med. 2008 27(14): 2509-2523
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

# estimate power
power = powerConLogistic.bin(
  N = 59, 
  power = NULL, 
  OR = 3.5, 
  pE = 0.15, 
  nD = 1, 
  nH = 2, 
  R2 = 0, 
  alpha = 0.05, 
  nTests = 1)

print(power) # 0.80

# estimate N (number of sets)
N = powerConLogistic.bin(
  N = NULL, 
  power = 0.80, 
  OR = 3.5, 
  pE = 0.15, 
  nD = 1, 
  nH = 2, 
  R2 = 0, 
  alpha = 0.05, 
  nTests = 1)

print(ceiling(N)) # 59

# estimate OR
OR = powerConLogistic.bin(
  N = 59, 
  power = 0.80, 
  OR = NULL, 
  pE = 0.15, 
  nD = 1, 
  nH = 2, 
  R2 = 0, 
  alpha = 0.05, 
  nTests = 1,
  OR.low = 1.01,
  OR.upp = 100)

print(OR) # 3.49



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ method }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
