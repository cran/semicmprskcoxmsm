\name{individual_RR_RD}
\alias{individual_RR_RD}

\title{
Estimating Three Individual Risk Difference and Risk Ratio Using the General Markov Model Conditional on Predicted Random Effect
}

\description{
\code{individual_RR_RD} estimates the individual risk difference and risk ratio based on the MSM illness-death general Markov model conditional on predicted random effect for each data point at a fixed time point.
}

\usage{
individual_RR_RD(dat1,res1,t1_star ,t)
}

\arguments{

\item{dat1}{
The dataset, includes non-terminal events, terminal events as well as event indicator.
}

\item{res1}{
The output from \code{em_illness_death_phmm_weight}, the general Markov model result, the result data includes the predicted random effect.
}

\item{t1_star}{
Fixed non-terminal event time for estimating risk difference/ratio for terminal event following the non-terminal event.
}

\item{t}{
Fixed time point of interest to compare the individual risk difference / ratio.
}


}

\value{

Returns a data frame that includes the individual risk difference / ratio for three type of events.

}

\details{

Similar as \code{cif_est_usual}, after estimating the parameters in the illness-death model \eqn{\lambda_{j}^a} using IPW, we could estimate the corresponding conditional CIF under the predicted b:

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

The frailty term, or equivalently, the random effect b represents the unobserved heterogeneity among the individuals. As such, the above conditional risk represents individual risk, and the risk contrasts the individual risk contrasts. We therefore have the individual risk difference (IRD) and the individual risk ratio (IRR).

Under the random effects model, for \eqn{i = 1,2,...,n}, the predicted random effect is \eqn{\hat{b}_i = E(b_i \mid O_i, \hat{\theta})}. We then obtain the predicted IRD and the predicted IRR.

}

\concept{Individual Risk Difference}
\concept{Individual Risk Ratio}
\concept{Frailty}
\concept{Multi state Model}


