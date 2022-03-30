\name{var_em_illness_death_phmm}
\alias{var_em_illness_death_phmm}

\title{
Variance of parameters in MSM Illness-death General Markov Model
}

\description{
Use bootstrap to obtain the variance estimator for parameters in MSM illness-death general markov model.
}

\usage{
var_em_illness_death_phmm(data,sigma_2_0,VARS.)
}

\arguments{

\item{data}{
The output dataset from \code{em_illness_death_phmm_weight}.
}

\item{sigma_2_0}{
Initial value for \eqn{\sigma^2}, the variance of zero-mean normal frailty, usually starts with 1.
}

\item{VARS.}{
Confounder sets.
}

}

\value{

List of bootstrap SE for all the parameters in the general Markov model

}

\details{
See \code{em_illness_death_phmm_weight}. In each bootstrap, the propensity score model needs to be re-fitted, and fit the MSM illness-death general markov model with new IP weights.
}

\concept{General Markov model}
\concept{Bootstrap Variance}
