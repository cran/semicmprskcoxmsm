\name{cif_est_usual}
\alias{cif_est_usual}

\title{
Estimating Three Cumulative Incidence Functions Using the Usual Markov Model
}

\description{
\code{cif_est_usual} estimates the cumulative incidence function (CIF, i.e.risk) based on the MSM illness-death usual Markov model.
}

\usage{
cif_est_usual(data,X1,X2,event1,event2,w,Trt,
              t1_star = t1_star)
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

\item{t1_star}{
Fixed non-terminal event time for estimating CIF function for terminal event following the non-terminal event.
}


}

\value{

Returns a table containing the estimated CIF for the event of interest for control and treated group.

}

\details{

After estimating the parameters in the illness-death model \eqn{\lambda_{j}^a} using IPW, we could estimate the corresponding CIF:

\deqn{
    \hat{P}(T_1^a<t,\delta_1^a=1) = \int_{0}^{t} \hat{S}^a(u) d\hat{\Lambda}_{1}^a(u),
}

\deqn{
    \hat{P}(T_2^a<t,\delta_1^a=0,\delta_2^a=1) = \int_{0}^{t} \hat{S}^a(u) d\hat{\Lambda}_{2}^a(u),
}

and

\deqn{
    \hat{P}(T_2^a<t_2 \mid T_1^a<t_1, T_2^a>t_1) = 1- e^{- \int_{t_1}^{t_2} d \hat{\Lambda}_{12}^a(u) },
}

where \eqn{\hat{S}^a} is the estimated overall survial function for joint \eqn{T_1^a, T_2^a}, \eqn{ \hat{S}^a(u) = e^{-\hat{\Lambda}_{1}^a(u)} - \hat{\Lambda}_{2}^a(u) }. We obtain three hazards by fitting the MSM illness-death model \eqn{ \hat\Lambda_{j}^a(u) = \hat\Lambda_{0j}(u)e^{\hat\beta_j*a} }, \eqn{ \hat\Lambda_{12}^a(u) = \hat\Lambda_{03}(u)e^{\hat\beta_3*a} }, and \eqn{ \hat\Lambda_{0j}(u) } is a Breslow-type estimator of the baseline cumulative hazard.

}

\concept{Cumulative Incidence Function}
\concept{Multi-state Model}

\references{
Meira-Machado, Luis and Sestelo, Marta (2019). \dQuote{Estimation in the progressive illness-death model: A nonexhaustive review,}
\emph{Biometrical Journal} 61(2), 245--263.
}
