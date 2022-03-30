\name{sim_cox_msm_semicmrsk}
\alias{sim_cox_msm_semicmrsk}

\title{
Simulating Semi-competing Risks with Right-censored Survival Data under Marginal Structural Illness-death Cox Model
}

\description{
The function to simulate semi-competing risk with right-censored survival data under marginal structural illness-death Cox model.
}

\usage{
sim_cox_msm_semicmrsk(beta1,beta2,beta3,sigma_2,
        alpha0,alpha1,alpha2,alpha3,
        n,Cens)
}

\arguments{

\item{beta1}{
True value of \eqn{\beta_1}  in the illness-death model.
}

\item{beta2}{
True value of \eqn{\beta_2}  in the illness-death model.
}

\item{beta3}{
True value of \eqn{\beta_3} in the illness-death model.
}

\item{sigma_2}{
True value of variance of normal frailty \eqn{\sigma^2} in the illness-death model, if \eqn{\sigma^2} = 0, then there is no frailty term.
}

\item{alpha0}{
True value of \eqn{\alpha_0}  in the propensity score model.
}

\item{alpha1}{
True value of \eqn{\alpha_1}  in the propensity score model.
}

\item{alpha2}{
True value of \eqn{\alpha_2}  in the propensity score model.
}

\item{alpha3}{
True value of \eqn{\alpha_3}  in the propensity score model.
}

\item{n}{
Sample size.
}

\item{Cens}{
Censoring distribution.
}

}

\value{

Returns a data frame that contains time to non-terminal event, T1, terminal event, T2 and censoring time C with their event indicator, delta1 and delta2. Three covariates Z1, Z2, Z3, and treatment assignment A are also included.

}

\details{

We simulate data followed by Xu(2010) to generate semi-competing risk data under illness-death model, where we have baseline hazard \eqn{\lambda_{01}(t) = \lambda_{02}(t) = 2exp(-t)I(0 \le t \le 3) + 2exp(-3)I(t \ge 3)}, and \eqn{\lambda_{03}(t) = 2\lambda_{01}(t)}.

We also have the propensity score model to generate treatment assignment \eqn{P_A = logit^{-1}(\alpha_0 + \alpha_1 Z_1 + \alpha_2 Z_2 + \alpha_3 Z3)}.

}

\concept{Illness-death model}
\concept{Cox model}
\concept{Semi-competing risk}

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
}
