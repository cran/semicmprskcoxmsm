\name{OUT_em_weights}
\alias{OUT_em_weights}

\title{
Initial Value For Fitting the General Markov Model
}

\description{
Compute the initial value for fitting the MSM illness-death general Markov model using EM type algorithm
}

\usage{
OUT_em_weights(data,X1,X2,event1,event2,w,Trt)
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

A list of vectors and dataframes:

\item{beta1}{
Initial value for \eqn{\beta_1}, the coefficient for the non-terminal event model.
}

\item{beta2}{
Initial value for \eqn{\beta_2}, the coefficient for the terminal event without non-terminal event model.
}

\item{beta3}{
Initial value for \eqn{\beta_3}, the coefficient for the terminal event following non-terminal event model.
}

\item{lambda1}{
Initial value for \eqn{ \lambda_{01} }, the estimated baseline hazard for the non-terminal event model.
}

\item{lambda2}{
Initial value for \eqn{ \lambda_{02} }, the estimated baseline hazard for the terminal event without non-terminal event model.
}

\item{lambda3}{
Initial value for \eqn{ \lambda_{03} }, the estimated baseline hazard for the terminal event following non-terminal event model.
}

\item{Lambda1}{
Initial value for \eqn{ \Lambda_{01} }, the estimated cumulative baseline hazard for the non-terminal event model.
}

\item{Lambda2}{
Initial value for \eqn{ \Lambda_{02} }, the estimated cumulative baseline hazard for the terminal event without non-terminal event model.
}

\item{Lambda3}{
Initial value for \eqn{ \Lambda_{03} }, the estimated cumulative baseline hazard for the terminal event following non-terminal event model.
}

\item{event1}{
An object of class \code{Surv} for non-terminal event.
}

\item{event2}{
An object of class \code{Surv} for terminal event without non-terminal event.
}

\item{event3}{
An object of class \code{Surv} for terminal event following non-terminal event.
}

}

\details{

See \code{usual_illness_death_weight}

}

\concept{EM algorithm}
\concept{General Markov Model}


\seealso{ \code{usual_illness_death_weight} }

\examples{
n <- 500
set.seed(1234)
Cens = runif(n,0.7,0.9)
set.seed(1234)
OUT1 <- sim_cox_msm_semicmrsk(beta1 = 1,beta2 = 1,beta3 = 0.5,
                              sigma_2 = 1,
                              alpha0 = 0.5, alpha1 = 0.1, alpha2 = -0.1, alpha3 = -0.2,
                              n=n, Cens = Cens)
data_test <- OUT1$data0

## Get the PS weights
vars <- c("Z1","Z2","Z3")
ps1 <- doPS(data = data_test,
            Trt = "A",
            Trt.name = 1,
            VARS. = vars,
            logistic = TRUE,w= NULL)
w <- ps1$Data$ipw_ate_stab

### Get the initial value
EM_initial <- OUT_em_weights(data = data_test,
                             X1 = "X1",
                             X2 = "X2",
                             event1 = "delta1",
                             event2 = "delta2",
                             w = w,
                             Trt = "A")
}

