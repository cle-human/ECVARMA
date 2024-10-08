% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{ecvarma.define}
\alias{ecvarma.define}
\title{Defining ECVARMA model}
\usage{
ecvarma.define(
  A0 = NULL,
  alphas,
  betas,
  AR = NULL,
  MA = NULL,
  Sigma,
  include = c("none", "SRconst", "LRconst")
)
}
\arguments{
\item{A0}{a matrix of parameters before the dependend variables (kxk)}

\item{alphas}{a matrix of adjustment parameters (k x r)}

\item{betas}{a matrix of cointegration parameters (r x k) including betas that are normalized to 1 or 0. One more column must be added to the end if a constant in the long-run equilibrium is desired.}

\item{AR}{a matrix of the parameters before p lagged first differences (k x kp) One more column must be added to the beginning if a constant in the short-run dynamics is desired.}

\item{MA}{a matrix of the parameters before q lagged error terms (k x kq)}

\item{Sigma}{the non-aggregated variance-covariance matrix of the error terms}

\item{include}{allows to add a constant in the short-run dynamics or the long-run equilibrium.}
}
\value{
a list (ECVARMA class)
}
\description{
Function to define an ECVARMA process, which can be used to find the aggregated model or to simulate data. \deqn{A_0\Delta Y_t=\alpha\beta^{'} Y_{t-1}+ \mu +\sum_{i=1}^{p} \Gamma _i  \Delta Y_{t-i}+u_t+\sum_{j=1}^{q} M _j  u_{t-j};\; \beta^{'}= (I_r\; \beta\; \beta_0);\; u_t\sim N(0,\Sigma)}, where \eqn{Y_t} is the time series data in levels, \eqn{\Delta Y_t} the first difference of the time series data and \eqn{u_t} the error term.
}
