\name{conditional_cif_b}
\alias{conditional_cif_b}

\title{
Estimating Three Conditional Cumulative Incidence Functions Using the General Markov Model Conditional on Random Effect
}

\description{
\code{conditional_cif_b} estimates the cumulative incidence function based on the MSM illness-death general Markov model conditional on the fixed random effect b.
}

\usage{
conditional_cif_b(res1,
                  t1_star,
                  b)
}

\arguments{

\item{res1}{
The output from \code{em_illness_death_phmm_weight}, the general Markov model result.
}

\item{t1_star}{
Fixed non-terminal event time for estimating CIF function for terminal event following the non-terminal event.
}

\item{b}{
Fixed random effect value.
}


}

\value{

\item{a1}{
The step function for estimated CIF conditional on b for time to non-terminal event for control group.
}

\item{b1}{
The step function for estimated CIF conditional on b for time to non-terminal event for treated group.
}

\item{a2}{
The step function for estimated CIF conditional on b for time to terminal event without non-terminal event for control group.
}

\item{b2}{
The step function for estimated CIF conditional on b for time to terminal event without non-terminal event for treated group.
}

\item{a3}{
The step function for estimated CIF conditional on b for time to terminal event following non-terminal event by t1_start for control group.
}

\item{b3}{
The step function for estimated CIF conditional on b for time to terminal event without non-terminal event by t1_start for treated group.
}

\item{cif.1}{
A data frame with time and estimated CIF conditional on b if is treated or controlled for time to non-terminal event.
}

\item{cif.2}{
A data frame with time and estimated CIF conditional on b if is treated or controlled for time to terminal event without non-terminal event.
}

\item{cif.3}{
A data frame with time and estimated CIF conditional on b if is treated or controlled for time to terminal event without non-terminal event by t1_start.
}

}

\details{

Similar as \code{cif_est_usual}, after estimating the parameters in the illness-death model \eqn{\lambda_{j}^a} using IPW, we could estimate the corresponding conditional CIF under fixed b:

\deqn{
    \hat{P}(T_1^a<t,\delta_1^a=1 \mid b) = \int_{0}^{t} \hat{S}^a(u \mid b) d\hat{\Lambda}_{1}^a(u \mid b ),
}

\deqn{
    \hat{P}(T_2^a<t,\delta_1^a=0,\delta_2^a=1 \mid b) = \int_{0}^{t} \hat{S}^a(u \mid b) d\hat{\Lambda}_{2}^a(u \mid b),
}

and

\deqn{
    \hat{P}(T_2^a<t_2 \mid T_1^a<t_1, T_2^a>t_1 \mid b) = 1- e^{- \int_{t_1}^{t_2} d \hat{\Lambda}_{12}^a(u \mid b) },
}

where \eqn{\hat{S}^a} is the estimated overall survial function for joint \eqn{T_1^a, T_2^a}, \eqn{ \hat{S}^a(u) = e^{-\hat{\Lambda}_{1}^a(u)} - \hat{\Lambda}_{2}^a(u) }. We obtain three hazards by fitting the MSM illness-death model \eqn{ \hat\Lambda_{j}^a(u) = \hat\Lambda_{0j}(u)e^{\hat\beta_j*a} }, \eqn{ \hat\Lambda_{12}^a(u) = \hat\Lambda_{03}(u)e^{\hat\beta_3*a} }, and \eqn{ \hat\Lambda_{0j}(u) } is a Breslow-type estimator of the baseline cumulative hazard.

where \eqn{ S(t \mid b;a) = \exp[- \int_0^{t} \{ \lambda_{01} (u)e^{\beta_1a + b} + \lambda_{02} (u )e^{\beta_2a + b} \} d u ] = \exp \{- e^{\beta_1a + b} \Lambda_{01}(t) - e^{\beta_2a + b}  \Lambda_{02} (t ) \} }

}

\concept{Conditional Cumulative Incidence Function}
\concept{Multi-state Model}
\concept{Frailty}


\seealso{ \code{cif_est_usual} }
