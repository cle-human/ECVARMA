#' Defining ECVARMA model
#' @description Function to define an ECVARMA process, which can be used to find the aggregated model or to simulate data. \deqn{A_0\Delta Y_t=\alpha\beta^{'} Y_{t-1}+ \mu +\sum_{i=1}^{p} \Gamma _i  \Delta Y_{t-i}+u_t+\sum_{j=1}^{q} M _j  u_{t-j};\; \beta^{'}= (I_r\; \beta\; \beta_0);\; u_t\sim N(0,\Sigma)}, where \eqn{Y_t} is the time series data in levels, \eqn{\Delta Y_t} the first difference of the time series data and \eqn{u_t} the error term.
#' @param A0 a matrix of parameters before the dependend variables (kxk)
#' @param alphas a matrix of adjustment parameters (k x r)
#' @param betas a matrix of cointegration parameters (r x k) including betas that are normalized to 1 or 0. One more column must be added to the end if a constant in the long-run equilibrium is desired.
#' @param AR a matrix of the parameters before p lagged first differences (k x kp) One more column must be added to the beginning if a constant in the short-run dynamics is desired.
#' @param MA a matrix of the parameters before q lagged error terms (k x kq)
#' @param Sigma the non-aggregated variance-covariance matrix of the error terms
#' @param include allows to add a constant in the short-run dynamics or the long-run equilibrium.
#' @return a list (ECVARMA class)
#' @export
ecvarma.define<-function(A0=NULL,alphas,betas,AR=NULL,MA=NULL,Sigma,include=c("none","SRconst","LRconst")){
  k<-nrow(alphas)
  if (is.null(A0)) {
    A0<-diag(k)
  }
  if (is.null(AR)) {
    p<-0
  }else{
    p<-(ncol(AR)-ifelse(include=="SRconst",1,0))/k
  }
  if (is.null(MA)) {
    q<-0
  }else{
    q<-(ncol(MA))/k
  }
  if (ncol(betas)<nrow(betas)) {
    betas<-t(betas)
  }
  r<-nrow(betas)
  out<-list(
    output="defined",
    model.specification=list(
      k = k,
      p = p,
      q = q,
      r = r,
      include = include
    ),
    A0 = A0,
    alphas = alphas,
    betas = betas,
    AR = AR,
    MA = MA,
    Sigma = Sigma
  )
  class(out)<-"ECVARMA"
  attr(out,"determination")<-"defined"
  return(out)
}

#' Simulation of ECVARMA process
#' @description A function to simulate an ECVARMA process.
#' @param ecvarma an ECVARMA object, created by ecvarma.define, an ECVARMA aggregation function or an estimated ECVARMA process
#' @param sample.size final sample size of the simulated data set
#' @param n number of omitted initial observations (burn in phase), minimal value is the biggest lag order of the model
#' @param seed value for set.seed function
#' @return a matrix of size (T x k)
#' @export

ecvarma.sim<-function(ecvarma,sample.size,n,seed=0){
  k<-ecvarma[["model.specification"]][["k"]]
  p<-ecvarma[["model.specification"]][["p"]]
  q<-ecvarma[["model.specification"]][["q"]]
  r<-ecvarma[["model.specification"]][["r"]]
  A0<-ecvarma$A0
  alphas<-ecvarma$alphas
  betas<-ecvarma$betas
  AR<-ecvarma$AR
  MA<-ecvarma$MA
  Sigma<-ecvarma$Sigma
  include<-ecvarma[["model.specification"]][["include"]]
  if (!all(A0==diag(k))) {
    A0i<-solve(A0)
    alphas<-A0i%*%alphas
    AR<-A0i%*%AR
    MA<-A0i%*%MA
  }

  n<-max(n,p,q)
  set.seed(seed)
  errors<-mvtnorm::rmvnorm(sample.size+n,sigma = Sigma)
  P<-matrix(0,sample.size+n,k)
  del.P<-matrix(0,sample.size+n,k)
  nn<-max(p,q,1)
  for (i in (nn+1):(sample.size+n)) {
    if (include=="LRconst") {
      del.P[i,]<-alphas%*%betas%*%c(P[i-1,],1)
    }else{
      del.P[i,]<-alphas%*%betas%*%P[i-1,]
    }
    if (include=="SRconst") {
      del.P[i,]<-del.P[i,]+AR[,1]
    }
    if (p>0) {
      for (j in 1:p-1) {
        del.P[i,]<-del.P[i,]+AR[,1:k+ifelse(include=="SRconst",1,0)+k*j]%*%del.P[i-1-j,]
      }
    }
    del.P[i,]<-del.P[i,]+errors[i,]
    if (q>0) {
      for (j in 1:q-1) {
        del.P[i,]<-del.P[i,]+MA[,1:k+k*j]%*%errors[i-1-j,]
      }
    }
    P[i,]<-del.P[i,]+P[i-1,]
  }
  P<-P[-(1:n),]
  return(P)
}
