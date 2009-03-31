\name{Oph}
\alias{Oph}
\docType{data}
\title{Ophthalmology Data}
\description{
The Ophthalmology data set is described in Example 14.41 on page 807 in Rosner (2006).
}
\usage{data(Oph)}
\format{
  A data frame with 354 observations on the following 3 variables.
  \describe{
    \item{\code{times}}{a numeric vector recording the survival/censoring time for each event/censoring.} 
    \item{\code{status}}{a numeric vector recording if a observed time is event time (\code{status=1}) 
      or censoring time (\code{status=0}).}
    \item{\code{group}}{a factor with levels \code{C} (indicating control group) and \code{E} (indicating 
      experimental group).}
  }
}
\details{
This data set was from a clinical trial (Berson et al., 1993) conducted to test the efficacy of different 
vitamin supplements in preventing visual loss in patients with retinitis pigmentosa. Rosner (2006) used the
data from this clinical trial to illustrate the analysis of survival data (Sections 14.9-14.12 of Rosner (2006)).

The data set consists of two groups of participants: 
(1) the experimental group (i.e., group E in which participants receiving 15,000 IU of vitamin A per day)
and (2) the control group (i.e., group C in which participants receiving 75 IU of vitamin A per day).

The participants were enrolled over a 2-year period (1984-1987) and followed for a maximum of 6 years. The
follow-up was terminated in September 1991. Some participants dropped out of the study before September 1991 
and had not failed. Dropouts were due to death, other diseases, or side effects possibly due to the 
study medications, or unwillingness to comply (take study medications). There are 6 time points (at 1st year, 
2nd year, 3rd year, 4th year, 5-th year, and 6-th year) in this data set.

Rosner (2006, page 786) defined
the participants who do not reach a disease endpoint during their period of follow-up as censored observations.
A participant has been censored at time \eqn{t} if the participant has been followed up to time \eqn{t} and has
not failed. Noninformative censoring is assumed. That is, participants who are censored have the same 
underlying survival curve after their censoring time as patients who are not censored.
}
\source{
Created based on Table 14.12 on page 787 of Rosner (2006).
}
\references{
  Berson, E.L., Rosner, B., Sandberg, M.A., Hayes, K.C., Nicholson, B.W., Weigel-DiFranco, C., and Willett, W.C. (1993).
  A randomized trial of vitamin A and vitamin E supplementation for retinitis pigmentosa.
  \emph{Archives of Ophthalmology}. 111:761-772.

  Rosner B. (2006).
  \emph{Fundamentals of Biostatistics}. (6-th edition). Thomson Brooks/Cole.
}
\examples{
data(Oph)
}
\keyword{datasets}
