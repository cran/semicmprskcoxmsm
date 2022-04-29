\name{bayesian_boot_irrd}
\alias{bayesian_boot_irrd}

\title{
Obtaining Bayesian Bootstrap Sample for Individual Risk Difference and Risk Ratio.
}

\description{
\code{bayesian_boot_irrd} provides the bootstrap sample for individual risk difference and risk ratio, it can be used for further inferences.
}

\usage{
bayesian_boot_irrd(dat2,B,sigma_2_0, EM_initial, varlist, t1_star,t)
}

\arguments{

\item{dat2}{
The dataset, includes non-terminal events, terminal events as well as event indicator.
}

\item{B}{
Number of bootstraps that the user want to run, typically we use B = 500.
}

\item{sigma_2_0}{
Initial value for sigma_2 for the general Markov model
}

\item{EM_initial}{
Initial value for the EM algorithm, the output of \code{OUT_em_weights}.
}

\item{varlist}{
Confounder list for the propensity score model.
}

\item{t1_star}{
Fixed non-terminal event time for estimating risk difference/ratio for terminal event following the non-terminal event.
}

\item{t}{
Fixed time point of interest to compare the individual risk difference / ratio.
}


}

\value{

\item{RD1_boot}{
A n times B matrix as the Bayesian bootstrap sample for each data point. The sample is for individual risk difference for time to non-terminal event at time t.
}

\item{RD2_boot}{
A n times B matrix as the Bayesian bootstrap sample for each data point. The sample is for individual risk difference for time to terminal event without non-terminal event at time t.
}

\item{RD3_boot}{
A n times B matrix as the Bayesian bootstrap sample for each data point. The sample is for individual risk difference for time to terminal event following non-terminal event by t1_start at time t.
}

\item{RR1_boot}{
A n times B matrix as the Bayesian bootstrap sample for each data point. The sample is for individual risk ratio for time to non-terminal event at time t.
}

\item{RR2_boot}{
A n times B matrix as the Bayesian bootstrap sample for each data point. The sample is for individual risk ratio for time to terminal event without non-terminal event at time t.
}

\item{RR3_boot}{
A n times B matrix as the Bayesian bootstrap sample for each data point. The sample is for individual risk ratio for time to terminal event following non-terminal event by t1_start at time t.
}



}

\details{

For each  bootstrap sample:

1. Generate \eqn{n} standard exponential (mean and variance 1) random variates : \eqn{u_1, u_2,..., u_n};

2. The weights for the Bayesian bootstrap are: \eqn{w_{i}^{boot} = u_i / \bar{u}}, where \eqn{\bar{u} = n^{-1}\sum_{i=1}^{n} u_i};

3. Calculate the propensity score and IP weights \eqn{w_{i}^{IPW}} based on Bayesian bootstrap weighted data, and assigned the weights for fitting the MSM general Markov model as \eqn{w_i = w_{i}^{boot} * w_{i}^{IPW}}.

4. After obtaining \eqn{\hat{\theta}} and \eqn{\hat{b}_i}, for each individual i, calculate the IRR and IRD by plugging \eqn{\hat{\theta}, \hat{b}_i} and a=0, a=1 separately at time t.

The 95\% prediction intervals (PI) cam be obtained by the normal approximation using bootstrap standard error.

}

\concept{Individual Risk Difference}
\concept{Individual Risk Ratio}
\concept{Frailty}
\concept{Multi state Model}
\concept{Bootstrap sample}


