
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ECVARMA

<!-- badges: start -->
<!-- badges: end -->

The ECVARMA package is a R package which allows to analyse error
correction processes with moving average dynamic. It includes three
estimation functions for EC-VARMA models in the following
specifications: final moving average (FMA) form, diagonal moving average
(DMA) from and the scalar component model (SCM). Custom restrictions can
be used as well. Additionally, information criteria for lag order
identification of the FMA form and for a non-parametric identification
of the cointegration rank are included too. Further, there are functions
to aggregate a model temporally and contemporanously.

## Installation

You can install the development version of ECVARMA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cle-human/ECVARMA")
```

## Example

In the following example, an EC-VARMA model is specified, simulated,
estimated and aggregated. The package works around ECVARMA objects.
These are the output of most functions and the input if parameter values
are needed. In this example, the SCM specification is used to estimate
the model. Here, we know the correct scalar components, but otherwise
the function “SCMid” from the “MTS” package can be used to determine
them.

``` r
library(ECVARMA)

#specification of the example 
specified_model<-ecvarma.define(A0 = NULL,
                         alphas = matrix(c(-0.16,-0.04,0.11,0.1,-0.08,0.02),3,2),
                         betas = cbind(diag(2),-1),
                         AR = matrix(c(-0.24,0.67,0.36,0.12,0.07,0.13,0.06,-0.17,0.03),3,3),
                         MA = matrix(c(0.85,0.49,0.55,0.26,0.1,-0.4,0.26,0,0.43),3,3),
                         Sigma = diag(3),
                         include = "none")
#simulation of example data
example_data<-ecvarma.sim(ecvarma = specified_model, #an ECVARMA object 
                     sample.size = 1000,
                     n = 100,
                     seed = 10)

#find cointegration rank
coint.criterion(example_data)

#find scalar components
#MTS::SCMid(example_data)

#estimation
estimated_model<-ecvarma.scm(input = example_data,
                             sc = matrix(c(2,1,
                                           2,1,
                                           2,1),3,2,byrow = TRUE), #scalar components
                             r=2, #cointegration rank
                             n=10, #number of estimation iterations
                             include = "none") #no constant included in the model
  
#temporal aggregation
aggregated_model<-ecvarma.agg.M(ecvarma = specified_model,
                                l=3, #aggregation period length 3
                                f="bop") #beginning of period values
```
