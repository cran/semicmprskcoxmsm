\name{em_illness_death_phmm_weight}
\alias{em_illness_death_phmm_weight}

\title{
Using EM Type Algorithm for MSM Illness-death General Markov Model
}

\description{
Under the general Markov illness-death model, with normal frailty term which is a latent variable. We use the EM type algorithm to estimate the coefficient in the MSM illness-death general Markov model.
}

\usage{
em_illness_death_phmm_weight(data,X1,X2,event1,event2,w,Trt,
                            EM_initial,sigma_2_0)
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

\item{EM_initial}{
Initial value for the EM algorithm, the output of \code{OUT_em_weights}.
}

\item{sigma_2_0}{
Initial value for \eqn{\sigma^2}, the variance of zero-mean normal frailty, usually starts with 1.
}

}

\value{

\item{beta1}{
The EM sequence for estimating \eqn{\beta_1} at each iteration.
}

\item{beta2}{
The EM sequence for estimating \eqn{\beta_2} at each iteration.
}

\item{beta3}{
The EM sequence for estimating \eqn{\beta_3} at each iteration.
}

\item{Lambda01}{
List of two dataframes for estimated \eqn{ \Lambda_{01} } and \eqn{ \lambda_{01} } when EM converges.
}

\item{Lambda02}{
List of two dataframes for estimated \eqn{ \Lambda_{02} } and \eqn{ \lambda_{02} } when EM converges.
}

\item{Lambda03}{
List of two dataframes for estimated \eqn{ \Lambda_{03} } and \eqn{ \lambda_{03} } when EM converges.
}

\item{sigma_2}{
The EM sequence for estimating \eqn{\sigma^2} at each iteration.
}

\item{loglik}{
The EM sequence for log-likelihood at each iteration.
}

\item{em.n}{
Number of EM steps to converge.
}

\item{data}{
Data after running the EM.
}

}

\details{

Similar as the usual Markov model. We postulate the semi-parametric Cox models with a frailty term for three transition rates in marginal structural illness-death model:

\deqn{
    \lambda_{1}(t_1 ; a) = \lambda_{01}(t)e^{\beta_1 a + b}, t_1 > 0 ;
}

\deqn{
    \lambda_{2}(t_2 ; a) = \lambda_{02}(t)e^{\beta_2 a + b}, t_2 > 0 ;
}

and

\deqn{
    \lambda_{12}(t_2 \mid t_1 ; a) = \lambda_{03}(t_2)e^{\beta_3 a + b}, 0 < t_1 <  t_2 ,
}

where \eqn{b \sim N(0,1)}. Since b is not observed in the data, we use the IP weighted EM type algorithm to estimate all the parameters in the MSM illness-death general Markov model.

}

\concept{General Markov model}
\concept{Semi-competing risks}
\concept{EM Algorithm}

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

### Fit the General Markov model
EM_initial <- OUT_em_weights(data = data_test,
                             X1 = "X1",
                             X2 = "X2",
                             event1 = "delta1",
                             event2 = "delta2",
                             w = w,
                             Trt = "A")

res1 <- em_illness_death_phmm_weight(data = data_test,
                                     X1 = "X1",
                                     X2 = "X2",
                                     event1 = "delta1",
                                     event2 = "delta2",
                                     w = w,
                                     Trt = "A",
                                     EM_initial = EM_initial,
                                     sigma_2_0 = 2)

print(paste("The estimated value for beta1 is:", round(res1$beta1[res1$em.n],5) ) )

}

