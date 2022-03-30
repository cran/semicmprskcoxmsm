\name{initial_lambda_em}
\alias{initial_lambda_em}

\title{
Compute the Initial (Cumulative) Baseline Hazard From the MSM Illness-death Model
}

\description{
Compute the Breslow type baseline hazard and cumulative baseline hazard at each event time from the MSM illness-death model.
}

\usage{
initial_lambda_em (OUT)
}

\arguments{

\item{OUT}{
The results of a \code{initial_fit_em_weights} fit.
}

}

\value{

A list contains six dataframes: including baseline hazard and cumulative baseline hazard for non-terminal event, terminal event without non-terminal event, and terminal event following non-terminal event.


}

\details{

See also \code{get_hazard}

}

\keyword{Baseline hazard}

\seealso{ \code{get_hazard} }
