\name{initial_fit_em_weights}
\alias{initial_fit_em_weights}

\title{
Fit the MSM Cox Model with IP Weights
}

\description{
Fit the MSM cox model with IPW as the initial value for EM algorithm to fit the illness-death general Markov model
}

\usage{
initial_fit_em_weights(data,X1,X2,event1,event2,w,Trt)
}

\arguments{

\item{data}{
The dataset, includes non-terminal events, terminal events as well as event indicator.
}

\item{X1}{
Time to non-terminal event, could be censored by terminal event or lost to follow up.
}

\item{X2}{
Time to terminal event, could be censored by lost to follow up.
}

\item{event1}{
Event indicator for non-terminal event.
}

\item{event2}{
Event indicator for terminal event.
}

\item{w}{
IP weights.
}

\item{Trt}{
Treatment variable.
}

}

\value{

A list of objects from \code{survival} package:

\item{event1}{
An object of class \code{Surv} for non-terminal event.
}

\item{event2}{
An object of class \code{Surv} for terminal event without non-terminal event.
}

\item{event3}{
An object of class \code{Surv} for terminal event following non-terminal event.
}

\item{fit1}{
An object of class \code{coxph} representing the fit for time to non-terminal event. See \code{coxph.object} for details.
}

\item{fit2}{
An object of class \code{coxph} representing the fit for time to terminal event without non-terminal event.
}

\item{fit3}{
An object of class \code{coxph} representing the fit for time to terminal event following non-terminal event.
}


}

\details{

As initial values we use for \eqn{\beta_j}, \eqn{j=1, 2, 3}, the estimates from IP weighted Cox regression without the offsets, i.e. from the usual Markov model.

}

\concept{Initial Value}
\concept{EM Algorithm}
\concept{General Markov Model}

\seealso{ \code{Surv}, \code{coxph} }

