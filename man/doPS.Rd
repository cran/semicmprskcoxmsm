\name{doPS}
\alias{doPS}

\title{
Generate the Inverse Probability Treatment Weights
}

\description{
\code{doPS} calculates the unstabilized and stabilized inverse probability treatment weights (IPW) for average treatment effect using propensity score. The propensity score is calculated by \code{twang} package using the boosted logistic regression.
}

\usage{
doPS(data,Trt,Trt.name,VARS.,logistic = FALSE)
}

\arguments{

\item{data}{
The dataset, includes treatment assignment as well as covariates.
}

\item{Trt}{
The name of the treatment variable in the dataset.
}

\item{Trt.name}{
The treated group name of the treatment variable in the dataset.
}

\item{VARS.}{
The vector of the name of potential confounding variables in the dataset.
}

\item{logistic}{
A logical value indicating whether use logistic regression (TRUE) or non-parametric boosted tree (FALSE).
}

}

\value{

doPS returns an object of class "PS".
An object of class "PS" is a list containing the following components:

\item{Data}{A new dataset which excludes all the missing value on the potential confounders from input data, add the propensity score and IPW into the new dataset.
 \describe{
   \item{ps_ate}{The estimated propensity scores with estimand of interest as ATE.}
   \item{ipw_ate_unstab}{Unstabilized ipw calculated from ps_ate.}
   \item{ipw_ate_stab}{Stabilized ipw calculated from ps_ate.}
  }}

\item{ps}{an object of class \code{ps}, See the help for \link[twang]{ps} for details of the ps class.
}

}

\details{

The treatment variable should only contain 2 levels of treatment, and one should be viewed as treated group and another is control group.

For stabilized weights:

For the treated individuals, we assign the IPW: w = Pr(T=1)/Pr(T=1|X=x), for control individuals, the stabilized weight is: w = (1-Pr(T=1))/(1-Pr(T=1|X=x)).

}

\concept{Propensity Score}
\concept{Inverse Probability Weighting}

\seealso{ \code{\link[twang]{ps}} }

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
            logistic = TRUE)
w <- ps1$Data$ipw_ate_stab
}
