\name{usual_illness_death_weight}
\alias{usual_illness_death_weight}

\title{
Fit MSM Illness-death Usual Markov Model For Semi-competing Risks Data
}

\description{
Fit the marginal structural three-state illness-death model with Cox representation and IP weights for semi-competing risks data. Inference under this model can be carried out using estimating equations with IP weights.
}

\usage{
usual_illness_death_weight(data,X1,X2,event1,event2,w,Trt)
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

A list of values and dataframes:

\item{beta1}{
Estimated \eqn{\beta_1}, the coefficient for the non-terminal event model.
}

\item{beta2}{
Estimated \eqn{\beta_2}, the coefficient for the terminal event without non-terminal event model.
}

\item{beta3}{
Estimated \eqn{\beta_3}, the coefficient for the terminal event following non-terminal event model.
}

\item{sd_beta1}{
Model based standard error for \eqn{\beta_1}.
}

\item{sd_beta2}{
Model based standard error for  \eqn{\beta_2}.
}

\item{sd_beta3}{
Model based standard error for  \eqn{\beta_3}.
}


\item{Lambda01}{
See also \code{get_hazard}. List of two dataframes for estimated \eqn{ \Lambda_{01} } and \eqn{ \lambda_{01} }, the estimated (cumulative) baseline hazard for the non-terminal event model.
}

\item{Lambda02}{
List of two dataframes for estimated \eqn{ \Lambda_{02} } and \eqn{ \lambda_{02} }, the estimated (cumulative) baseline hazard for the terminal event without non-terminal event model.
}

\item{Lambda03}{
List of two dataframes for estimated \eqn{ \Lambda_{03} } and \eqn{ \lambda_{03} }, the estimated (cumulative) baseline hazard for the terminal event following non-terminal event model.
}

}

\details{

Let \eqn{T_1}, \eqn{T_2} be the time to non-terminal event and terminal event, A be the treatment assignment. We postulate the semi-parametric Cox models for three transition rates in marginal structural illness-death model:

\deqn{
    \lambda_{1}(t_1 ; a) = \lambda_{01}(t)e^{\beta_1 a}, t_1 > 0 ;
}

\deqn{
    \lambda_{2}(t_2 ; a) = \lambda_{02}(t)e^{\beta_2 a}, t_2 > 0 ;
}

and

\deqn{
    \lambda_{12}(t_2 \mid t_1 ; a) = \lambda_{03}(t_2)e^{\beta_3 a}, 0 < t_1 <  t_2 .
}

The coefficients as well as Breslow type baseline hazards can be estimated by fitting the IP weights Cox proportional hazards models. Meanwhile, if we assume the estimated weights as known, then the robust sandwich variance estimator can be used to obtain the estimated variance.

The usual Markov model is also the same as the initial value for the general Markov model.

}

\concept{Usual Markov model}
\concept{Semi-competing risks}
\concept{Illness-death model}

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
            logistic = TRUE,w=NULL)
w <- ps1$Data$ipw_ate_stab

### Fit the Usual Markov model
res1 <- usual_illness_death_weight(data = data_test,
                                   X1 = "X1",
                                   X2 = "X2",
                                   event1 = "delta1",
                                   event2 = "delta2",
                                   w = w,
                                   Trt = "A")
print(paste("The estimated value for beta1 is:", round(res1$beta1,5) ) )

}

