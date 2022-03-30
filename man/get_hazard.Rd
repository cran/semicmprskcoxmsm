\name{get_hazard}
\alias{get_hazard}

\title{
Compute the (Cumulative) Baseline Hazard from Cox Model
}

\description{
Compute the Breslow type baseline hazard and cumulative baseline hazard at each event time from a Cox model.
}

\usage{
get_hazard(fit)
}

\arguments{

\item{fit}{
The results of a \code{coxph} fit.
}

}

\value{

A list contains two dataframes.

\item{Lambda}{
See also \code{basehaz}, returns the Breslow type cumulative baseline hazard.
}

\item{lambda}{
Returns the Breslow type baseline hazard.
}


}

\details{

See also \code{basehaz}, we only extract the estimated baseline hazard and baseline cumulative hazard from the results of a \code{coxph} fit.

}

\keyword{Baseline hazard}

\seealso{ \code{basehaz} }
