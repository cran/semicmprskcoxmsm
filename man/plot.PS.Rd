\name{plot.PS}
\alias{plot.PS}

\title{
Plotting Histogram of Propensity Score and Balancing Plot for Covariates in the Propensity Score Model
}

\description{
Displays a the histogram plots for the propensity score, stratified by treated and control group and a graph of standardized mean difference of potential confounders before and after weigthing.
}

\usage{
\method{plot}{PS}(x,...)
}

\arguments{

\item{x}{The results of \code{doPS} function.}

\item{...}{the other arguments you want to put in the built-in plot function.}

}

\details{
Only available when logistic = FALSE in doPS. The standardized mean difference (SMD), defined as the (weighted) treatment group mean minus the (weighted) control group mean divided by the (weighted) pooled sample (treatment and control) standard deviation. SMD between -0.1 and 0.1 typically indicates good balance.
}

\value{
Histogram of propensity score and balancing plot for covariates in the propensity score model corresponding to the output from \code{doPS}.
}

\keyword{Standardized Mean Difference}
\keyword{Histogram}

\seealso{ \code{\link[twang]{bal.table}} }




