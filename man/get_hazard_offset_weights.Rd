\name{get_hazard_offset_weights}
\alias{get_hazard_offset_weights}

\title{
Compute the (Cumulative) Baseline Hazard from Cox Model with Offsets
}

\description{
Compute the Breslow type baseline hazard and cumulative baseline hazard at each event time from a weighted Cox model with offsets.
}

\usage{
get_hazard_offset_weights(fit,data,time1= NULL,time2,w)
}

\arguments{

\item{fit}{
The results of a weighted \code{coxph} fit.
}

\item{data}{
The original data for fitting the weighted Cox model.
}

\item{time1}{
The default is \code{NULL}. For left truncation data, which refers to transition rate for terminal event following non-terminal events, this argument is the time to non-terminal event.
}

\item{time2}{
For right censored data, this is the event time or censoring time. For left truncation data, the argument is the time to terminal event or the censoring time.
}

\item{w}{
IP weights.
}

}

\value{

A list contains two dataframes.

\item{Lambda}{
See also \code{get_hazard}, returns a step function for cumulative baseline hazard.
}

\item{lambda}{
Returns a dataframe for baseline hazard.
}

\item{cum_base_haz}{
Returns a dataframe for cumulative baseline hazard.
}

}

\details{

See also \code{get_hazard}, handles the offset term in \code{coxph} for predicting the baseline hazard.

}

\concept{Baseline hazard}
\concept{Offset}

\seealso{ \code{get_hazard}, \code{basehaz} }
