#wrapper functions
#' ECVARMA estimation with SCM restrictions
#' @description Estimation algorithm for: \deqn{A_0\Delta Y_t=\alpha\beta^{'} Y_{t-1}+ \mu +\sum_{i=1}^{p} \Gamma _i  \Delta Y_{t-i}+u_t+\sum_{j=1}^{q} M _j  u_{t-j};\; \beta^{'}= (I_r\; \beta_{free}\; \beta_0);\; u_t\sim N(0,\Sigma)}, where \eqn{Y_t} is the time series data in levels, \eqn{\Delta Y_t} the first difference of the time series data and \eqn{u_t} the error term. The cointegration parameters for the first r variables are normalized to an identity matrix. Only \eqn{\beta_{free}} and \eqn{\beta_0} are part of the estimation.
#'
#' @param input a matrix of time series data in levels (the bigger dimension is considered to be t)
#' @param sc a matrix of scalar components, first column AR order including lag for error correction, second column MA order
#' @param r number of cointegrating vectors, only for 0<r<k
#' @param n number of interations performed with the optimization algorithm
#' @param include inclusion of a constant in the short run dynamics, cointegration relationship or not at all
#' @param initial_values a vector of first set of all parameter values \eqn{vec(\Sigma,A_0,\beta_{free},\beta_0,\mu,\alpha,\Gamma,M)}, \eqn{\beta_0} and \eqn{\mu} only if a constant in the cointegration relationship or in the short run dynamics is included respectively. With "auto", these are determined automatically, but in rare cases these might not work with the iteration algorithm.
#' @param h number of lags in level VAR estimation to determine first values of the error terms, only relevant when initial_values is set to "auto"
#' @param Newton option to use either Newton method or nlminb to estimate parameters
#' @param p.values whether p-values should be calculated, only possible if last iteration didn't fail
#'
#' @returns an ECVARMA object
#' @references Athanasopoulos, G., Poskitt, D. S., Vahid, F., & Yao, W. (2016). Determination of Long‐run and Short‐run Dynamics in EC‐VARMA Models via Canonical Correlations. \emph{Journal of Applied Econometrics}, 31(6), 1100-1119.
#' @references Yap, S. F., & Reinsel, G. C. (1995). Estimation and testing for unit roots in a partially nonstationary vector autoregressive moving average model.  \emph{Journal of the American Statistical Association}, 90(429), 253-267.
#' @references Kascha, C., & Trenkler, C. (2015). Simple identification and specification of cointegrated VARMA models. \emph{Journal of Applied Econometrics}, 30(4), 675-702.
#' @export
ecvarma.scm<-function(input,sc,r,n=5,include=c("none","SRconst","LRconst"),initial_values="auto",h="auto",Newton=TRUE,p.values=TRUE){
  #preparation
  input<-as.matrix(input)
  if(dim(input)[1]>dim(input)[2]){input<-t(input)}
  k<-dim(input)[1]
  T<-dim(input)[2]
  d<-k-r
  if (h=="auto") {
    h<-floor(max(log(T),tsDyn::lags.select(t(input))[["AIC_min"]][2]))
  }
  #if(R==c(NULL,"scm","fma","own"))stop("Choose restrictions")
  #if(all(include==c("none","SRconst","LRconst")))include="none"
  if(!(include%in% list("none","LRconst","SRconst")))stop("include is neither \"none\", \"LRconst\" nor \"SRconst\"")
  if (r==k)stop("r=k isn't implemented")
  if (r==0)stop("r=0 isn't implemented, use ecvarma with suitable restrictions")
  if (any(dim(sc)!=c(k,2)))stop("sc matrix has false dimensions")
  #restrictions
  res<-ecvarma.res(R="SCM",r=r,include=include,sc,p=NULL,q=NULL,k=k,d=d)
  R<-res$R
  R0<-res$R0
  R.f<-res$R.f
  R0.f<-res$R0.f
  p<-max(1,max(sc[,1])-1)
  q<-max(1,max(sc[,2]))
  #initial value estimation
  if (all(initial_values=="auto")) {
    initial_values<-ECVARMA.ini(input=input,p=p,r=r,q=q,h=h,k=k,T=T,R.f=R.f,R0.f=R0.f,include=include)
    del<-initial_values$del
    sig0<-initial_values$sig
  }else{
    if (length(initial_values)!=(dim(R)[1]+k^2)) stop(paste("initial_values must be a vector of length ",dim(R)[1]+k^2))
    sig0<-matrix(initial_values[1:(k^2)],k,k)
    del<-as.matrix(c(initial_values[-(1:k^2)]),dim(R)[1],1)
  }

  #iteration
  est_results<-ECVARMA.iter(input=input,del=del,sig0=sig0,k=k,T=T,p=p,r=r,q=q,n=n,R=R,R0=R0,include=include,Newton=Newton)

  #p-values
  if (p.values) {
    p.values<-tryCatch(expr = {ECVARMA.pvalues(input=input,T=T,k=k,d=d,del.r=as.vector(est_results[["del.r"]]),M=est_results[["M"]],bet=est_results[["bet"]],p=p,r=r,q=q,R=R,include=include,Sig=est_results[["Sig"]],u=est_results[["errors"]])},error=function(e){
      p.values<-list()
      #p.values$del.r<-"calculation failed, because no further iteration was possible"
      p.values$std.error<-"calculation failed, because no further iteration was possible"
      p.values$p.values<-"calculation failed, because no further iteration was possible"
      return(p.values)})
    #function(e){return(message="no p-values could be calculated, because no further iteration was possible")})
  }else{
    p.values<-list()
    #p.values$del.r<-"not calculated"
    p.values$std.error<-"not calculated"
    p.values$p.values<-"not calculated"
  }

  #output
  out<-list(
    output="estimated",
    model.specification=list(
      k = k,
      p = p,
      q = q,
      r = r,
      include = include
    ),
    A0 = est_results[["A0"]],
    alphas = est_results[["alp"]],
    betas = cbind(diag(r),est_results[["bet"]]),
    AR = est_results[["G"]],
    MA = est_results[["M"]],
    Sigma = est_results[["Sig"]],
    estimation = list(
      sample.size=T,
      specification.scheme = "scalar component method",
      scalar.components = sc,
      R = R,
      R0 = R0,
      del.r = as.vector(est_results[["del.r"]]),
      del.r2 = p.values$del.r,
      std.errors = p.values$std.error,
      p.values = p.values$p.values,
      error.terms = est_results[["errors"]],
      initital.values = initial_values,
      interations = "no iteration values with nleqslv"
    )
  )
  if (Newton) {
    out[["estimation"]][["interations"]] = est_results[["iterations"]]
  }
  class(out)<-"ECVARMA"
  attr(out,"determination")<-"estimated"
  return(out)

  #results<-list()
  #results$est_results<-est_results
  #results$p.values<-p.values
  #return(results)
}


#' ECVARMA estimation with FMA restrictions
#' @description Estimation algorithm for: \deqn{A_0\Delta Y_t=\alpha\beta^{'} Y_{t-1}+ \mu +\sum_{i=1}^{p} \Gamma _i  \Delta Y_{t-i}+u_t+\sum_{j=1}^{q} M _j  u_{t-j};\; \beta^{'}= (I_r\; \beta_{free}\; \beta_0);\; u_t\sim N(0,\Sigma)}, where \eqn{Y_t} is the time series data in levels, \eqn{\Delta Y_t} the first difference of the time series data and \eqn{u_t} the error term. The cointegration parameters for the first r variables are normalized to an identity matrix. Only \eqn{\beta_{free}} and \eqn{\beta_0} are part of the estimation. The \eqn{M_j}'s are restricted to scalars, i.e. the MA dynamic is identical for each error series and each dependent variable only uses lagged error terms from its own error series.
#'
#' @param input a matrix of time series data in levels (the bigger dimension is considered to be t)
#' @param p number of lagged first differences
#' @param r number of cointegrating vectors, only for 0<r<k
#' @param q number of lagged error terms
#' @param n number of interations performed with the optimization algorithm
#' @param include inclusion of a constant in the short run dynamics, cointegration relationship or not at all
#' @param initial_values a vector of first set of all parameter values \eqn{vec(\Sigma,A_0,\beta_{free},\beta_0,\mu,\alpha,\Gamma,M)}, \eqn{\beta_0} and \eqn{\mu} only if a constant in the cointegration relationship or in the short run dynamics is included respectively. With "auto", these are determined automatically, but in rare cases these might not work with the iteration algorithm.
#' @param h number of lags in level VAR estimation to determine first values of the error terms, only relevant when initial_values is set to "auto"
#' @param Newton option to use either Newton method or nlminb to estimate parameters
#' @param p.values whether p-values should be calculated, only possible if last iteration didn't fail
#'
#' @return an ECVARMA object
#' @references Yap, S. F., & Reinsel, G. C. (1995). Estimation and testing for unit roots in a partially nonstationary vector autoregressive moving average model.  \emph{Journal of the American Statistical Association}, 90(429), 253-267.
#' @references Kascha, C., & Trenkler, C. (2015). Simple identification and specification of cointegrated VARMA models. \emph{Journal of Applied Econometrics}, 30(4), 675-702.
#' @export
ecvarma.fma<-function(input,p,r,q,n=5,include=c("none","SRconst","LRconst"),initial_values="auto",h="auto",Newton=TRUE,p.values=TRUE){
  #preparation
  input<-as.matrix(input)
  if(dim(input)[1]>dim(input)[2]){input<-t(input)}
  k<-dim(input)[1]
  T<-dim(input)[2]
  d<-k-r
  if (h=="auto") {
    h<-floor(max(log(T),tsDyn::lags.select(t(input))[["AIC_min"]][2]))
  }
  if(!(include%in% list("none","LRconst","SRconst")))stop("include is neither \"none\", \"LRconst\" nor \"SRconst\"")
  #restrictions
  res<-ecvarma.res(R="FMA",r=r,include=include,sc=NULL,p=p,q=q,k=k,d=d)
  R<-res$R
  R0<-res$R0
  R.f<-res$R.f
  R0.f<-res$R0.f
  #initial value estimation
  if (all(initial_values=="auto")) {
    initial_values<-ECVARMA.ini(input=input,p=p,r=r,q=q,h=h,k=k,T=T,R.f=R.f,R0.f=R0.f,include=include)
    del<-initial_values$del
    sig0<-initial_values$sig
  }else{
    if (length(initial_values)!=(dim(R)[1]+k^2)) stop(paste("initial_values must be a vector of length",dim(R)[1]))
    sig0<-matrix(initial_values[1:9],k,k)
    del<-as.matrix(c(initial_values[-(1:k^2)]),dim(R)[1],1)
  }

  #iteration
  est_results<-ECVARMA.iter(input=input,del=del,sig0=sig0,k=k,T=T,p=p,r=r,q=q,n=n,R=R,R0=R0,include=include,Newton=Newton)

  #p-values
  if (p.values) {
    p.values<-tryCatch(expr = {ECVARMA.pvalues(input=input,T=T,k=k,d=d,del.r=as.vector(est_results[["del.r"]]),M=est_results[["M"]],bet=est_results[["bet"]],p=p,r=r,q=q,R=R,include=include,Sig=est_results[["Sig"]],u=est_results[["errors"]])},error=function(e){
      p.values<-list()
      #p.values$del.r<-"calculation failed, because no further iteration was possible"
      p.values$std.error<-"calculation failed, because no further iteration was possible"
      p.values$p.values<-"calculation failed, because no further iteration was possible"
      return(p.values)})
    #function(e){return(message="no p-values could be calculated, because no further iteration was possible")})
  }else{
    p.values<-list()
    #p.values$del.r<-"not calculated"
    p.values$std.error<-"not calculated"
    p.values$p.values<-"not calculated"
  }

  #output
  out<-list(
    output="estimated",
    model.specification=list(
      k = k,
      p = p,
      q = q,
      r = r,
      include = include
    ),
    A0 = est_results[["A0"]],
    alphas = est_results[["alp"]],
    betas = cbind(diag(r),est_results[["bet"]]),
    AR = est_results[["G"]],
    MA = est_results[["M"]],
    Sigma = est_results[["Sig"]],
    estimation = list(
      sample.size=T,
      specification.scheme = "final moving average",
      R = R,
      R0 = R0,
      del.r = as.vector(est_results[["del.r"]]),
      std.errors = p.values$std.error,
      p.values = p.values$p.values,
      error.terms = est_results[["errors"]],
      initital.values = initial_values,
      interations = "no iteration values with nleqslv"
    )
  )
  if (Newton) {
    out[["estimation"]][["interations"]] = est_results[["iterations"]]
  }
  class(out)<-"ECVARMA"
  attr(out,"determination")<-"estimated"
  return(out)
  # results<-list()
  #
  # results$est_results<-est_results
  # results$p.values<-p.values
  #
  # return(results)
}


#' ECVARMA estimation with DMA restrictions
#' @description Estimation algorithm for: \deqn{A_0\Delta Y_t=\alpha\beta^{'} Y_{t-1}+ \mu +\sum_{i=1}^{p} \Gamma _i  \Delta Y_{t-i}+u_t+\sum_{j=1}^{q} M _j  u_{t-j};\; \beta^{'}= (I_r\; \beta_{free}\; \beta_0);\; u_t\sim N(0,\Sigma)}, where \eqn{Y_t} is the time series data in levels, \eqn{\Delta Y_t} the first difference of the time series data and \eqn{u_t} the error term. The cointegration parameters for the first r variables are normalized to an identity matrix. Only \eqn{\beta_{free}} and \eqn{\beta_0} are part of the estimation. The \eqn{M_j}'s are restricted to diagonal matrices, i.e. the MA dynamic isn't identical for each error series, but each dependent variable only uses lagged error terms from its own error series.
#'
#' @param input a matrix of time series data in levels (the bigger dimension is considered to be t)
#' @param p number of lagged first differences
#' @param r number of cointegrating vectors, only for 0<r<k
#' @param q number of lagged error terms
#' @param n number of interations performed with the optimization algorithm
#' @param include inclusion of a constant in the short run dynamics, cointegration relationship or not at all
#' @param initial_values a vector of first set of all parameter values \eqn{vec(\Sigma,A_0,\beta_{free},\beta_0,\mu,\alpha,\Gamma,M)}, \eqn{\beta_0} and \eqn{\mu} only if a constant in the cointegration relationship or in the short run dynamics is included respectively. With "auto", these are determined automatically, but in rare cases these might not work with the iteration algorithm.
#' @param h number of lags in level VAR estimation to determine first values of the error terms, only relevant when initial_values is set to "auto"
#' @param Newton option to use either Newton method or nlminb to estimate parameters
#' @param p.values whether p-values should be calculated, only possible if last iteration didn't fail
#'
#' @return an ECVARMA object
#' @references Yap, S. F., & Reinsel, G. C. (1995). Estimation and testing for unit roots in a partially nonstationary vector autoregressive moving average model.  \emph{Journal of the American Statistical Association}, 90(429), 253-267.
#' @references Kascha, C., & Trenkler, C. (2015). Simple identification and specification of cointegrated VARMA models. \emph{Journal of Applied Econometrics}, 30(4), 675-702.
#' @export
ecvarma.dma<-function(input,p,r,q,n=5,include=C("none","SRconst","LRconst"),initial_values="auto",h="auto",Newton=TRUE,p.values=TRUE){
  #preparation
  input<-as.matrix(input)
  if(dim(input)[1]>dim(input)[2]){input<-t(input)}
  k<-dim(input)[1]
  T<-dim(input)[2]
  d<-k-r
  if (h=="auto") {
    h<-floor(max(log(T),tsDyn::lags.select(t(input))[["AIC_min"]][2]))
  }
  if(!(include%in% list("none","LRconst","SRconst")))stop("include is neither \"none\", \"LRconst\" nor \"SRconst\"")
  #restrictions
  res<-ecvarma.res(R="DMA",r=r,include=include,sc=NULL,p=p,q=q,k=k,d=d)
  R<-res$R
  R0<-res$R0
  R.f<-res$R.f
  R0.f<-res$R0.f
  #initial value estimation
  if (all(initial_values=="auto")) {
    initial_values<-ECVARMA.ini(input=input,p=p,r=r,q=q,h=h,k=k,T=T,R.f=R.f,R0.f=R0.f,include=include)
    del<-initial_values$del
    sig0<-initial_values$sig
  }else{
    if (length(initial_values)!=(dim(R)[1]+k^2)) stop(paste("initial_values must be a vector of length",dim(R)[1]))
    sig0<-matrix(initial_values[1:9],k,k)
    del<-as.matrix(c(initial_values[-(1:k^2)]),dim(R)[1],1)
  }

  #iteration
  est_results<-ECVARMA.iter(input=input,del=del,sig0=sig0,k=k,T=T,p=p,r=r,q=q,n=n,R=R,R0=R0,include=include,Newton=Newton)

  #p-values
  if (p.values) {
    p.values<-tryCatch(expr = {ECVARMA.pvalues(input=input,T=T,k=k,d=d,del.r=as.vector(est_results[["del.r"]]),M=est_results[["M"]],bet=est_results[["bet"]],p=p,r=r,q=q,R=R,include=include,Sig=est_results[["Sig"]],u=est_results[["errors"]])},error=function(e){
      p.values<-list()
      #p.values$del.r<-"calculation failed, because no further iteration was possible"
      p.values$std.error<-"calculation failed, because no further iteration was possible"
      p.values$p.values<-"calculation failed, because no further iteration was possible"
      return(p.values)})
    #function(e){return(message="no p-values could be calculated, because no further iteration was possible")})
  }else{
    p.values<-list()
    #p.values$del.r<-"not calculated"
    p.values$std.error<-"not calculated"
    p.values$p.values<-"not calculated"
  }

  #output
  out<-list(
    output="estimated",
    model.specification=list(
      k = k,
      p = p,
      q = q,
      r = r,
      include = include
    ),
    A0 = est_results[["A0"]],
    alphas = est_results[["alp"]],
    betas = cbind(diag(r),est_results[["bet"]]),
    AR = est_results[["G"]],
    MA = est_results[["M"]],
    Sigma = est_results[["Sig"]],
    estimation = list(
      sample.size=T,
      specification.scheme = "diagaonal moving average",
      R = R,
      R0 = R0,
      del.r = as.vector(est_results[["del.r"]]),
      std.errors = p.values$std.error,
      p.values = p.values$p.values,
      error.terms = est_results[["errors"]],
      initital.values = initial_values,
      interations = "no iteration values with nleqslv"
    )
  )
  if (Newton) {
    out[["estimation"]][["interations"]] = est_results[["iterations"]]
  }
  class(out)<-"ECVARMA"
  attr(out,"determination")<-"estimated"
  return(out)
  # results<-list()
  #
  # results$est_results<-est_results
  # results$p.values<-p.values
  #
  # return(results)
}



#' ECVARMA estimation with restrictions
#' @description Estimation algorithm for: \deqn{A_0\Delta Y_t=\alpha\beta^{'} Y_{t-1}+ \mu +\sum_{i=1}^{p} \Gamma _i  \Delta Y_{t-i}+u_t+\sum_{j=1}^{q} M _j  u_{t-j};\; \beta^{'}= (I_r\; \beta_{free}\; \beta_0);\; u_t\sim N(0,\Sigma)}, where \eqn{Y_t} is the time series data in levels, \eqn{\Delta Y_t} the first difference of the time series data and \eqn{u_t} the error term. The cointegration parameters for the first r variables are normalized to an identity matrix. Only \eqn{\beta_{free}} and \eqn{\beta_0} are part of the estimation.
#'
#' @param input a matrix of time series data in levels (the bigger dimension is considered to be t)
#' @param p number of lagged first differences
#' @param r number of cointegrating vectors, only for r<k, for r=0 use R to fix the alphas to 0
#' @param q number of lagged error terms
#' @param n number of interations performed with the optimization algorithm
#' @param h number of variable lags for first values of the error terms
#' @param R parameter restriction matrix in \eqn{vec(A_0,\beta_{free},\beta_0,\mu,\alpha,\Gamma,M)=R\delta+R_0}. \eqn{\delta} are the parameters that are actually estimated.
#' @param R0 a vector of constant values in \eqn{vec(A_0,\beta_{free},\beta_0,\mu,\alpha,\Gamma,M)=R\delta+R_0}
#' @param include inclusion of a constant in the short run dynamics, cointegration relationship or not at all
#' @param initial_values a vector of first set of all parameter values \eqn{vec(\Sigma,A_0,\beta_{free},\beta_0,\mu,\alpha,\Gamma,M)}, \eqn{\beta_0} and \eqn{\mu} only if a constant in the cointegration relationship or in the short run dynamics is included respectively. With "auto", these are determined automatically, but in rare cases these might not work with the iteration algorithm.
#' @param h number of lags in level VAR estimation to determine first values of the error terms, only relevant when initial_values is set to "auto"
#' @param Newton option to use either Newton method or nlminb to estimate parameters
#' @param p.values whether p-values should be calculated, only possible if last iteration didn't fail
#'
#' @return an ECVARMA object
#' @references Yap, S. F., & Reinsel, G. C. (1995). Estimation and testing for unit roots in a partially nonstationary vector autoregressive moving average model.  \emph{Journal of the American Statistical Association}, 90(429), 253-267.
#' @references Kascha, C., & Trenkler, C. (2015). Simple identification and specification of cointegrated VARMA models. \emph{Journal of Applied Econometrics}, 30(4), 675-702.
#' @export
ecvarma<-function(input,p,r,q,n=5,include=c("none","SRconst","LRconst"),R,R0=NULL,initial_values="auto",h="auto",Newton=TRUE,p.values=TRUE){
  #preparation
  input<-as.matrix(input)
  if(dim(input)[1]>dim(input)[2]){input<-t(input)}
  k<-dim(input)[1]
  T<-dim(input)[2]
  d<-k-r
  if (h=="auto") {
    h<-floor(max(log(T),tsDyn::lags.select(t(input))[["AIC_min"]][2]))
  }
  if(!(include%in% list("none","LRconst","SRconst")))stop("include is neither \"none\", \"LRconst\" nor \"SRconst\"")
  if (p==0) stop("p=0 isn't implemented, use R to set AR lags to zero")
  if (q==0) stop("q=0 isn't implemented, use R to set MA lags to zero")
  #restrictions for full estimation
  Rpq<-R[-(1:(r*d+r*k+k^2+ifelse(include=="LRconst",r,0))),]
  for(g in dim(Rpq)[2]:1){
    if(all(Rpq[,g]==0)){
      Rpq<-Rpq[,-g]
    }
  }
  if (all(R[(r*d+1):(r*d+r*k)+k^2,]==0)) {
    R.f<-rbind(cbind(diag(0,k^2+ifelse(include=="LRconst",k,0)),matrix(0,k*k+ifelse(include=="LRconst",k,0),dim(Rpq)[2])),cbind(matrix(0,dim(Rpq)[1],k*k+ifelse(include=="LRconst",k,0)),Rpq))
  }else{
    R.f<-rbind(cbind(diag(1,k^2+ifelse(include=="LRconst",k,0)),matrix(0,k*k+ifelse(include=="LRconst",k,0),dim(Rpq)[2])),cbind(matrix(0,dim(Rpq)[1],k*k+ifelse(include=="LRconst",k,0)),Rpq))
  }
  R.f<-rbind(cbind(R[1:k^2,1:k^2],matrix(0,k^2,dim(R.f)[2])),cbind(matrix(0,dim(R.f)[1],k^2),R.f))
  for(g in dim(R.f)[2]:1){
    if(all(R.f[,g]==0)){
      R.f<-R.f[,-g]
    }
  }
  if (all(is.null(R0))) {
    R0<-matrix(0,dim(R)[1],1)
  }
  R0.f<-R0[-(1:(r*d+r*k+ifelse(include=="LRconst",r,0))+k^2),]
  R0.f<-c(R0.f[1:k^2],matrix(0,k^2+ifelse(include=="LRconst",k,0),1),R0.f[-(1:k^2)])

  #initial value estimation
  if (all(initial_values=="auto")) {
    initial_values<-ECVARMA.ini(input=input,p=p,r=r,q=q,h=h,k=k,T=T,R.f=R.f,R0.f=R0.f,include=include)
    del<-initial_values$del
    sig0<-initial_values$sig
  }else{
    if (length(initial_values)!=(dim(R)[1]+k^2)) stop(paste("initial_values must be a vector of length",dim(R)[1]))
    sig0<-matrix(initial_values[1:9],k,k)
    del<-as.matrix(c(initial_values[-(1:k^2)]),dim(R)[1],1)
  }

  #iteration
  est_results<-ECVARMA.iter(input=input,del=del,sig0=sig0,k=k,T=T,p=p,r=r,q=q,n=n,R=R,R0=R0,include=include,Newton=Newton)


  #p-values
  if (p.values) {
    p.values<-tryCatch(expr = {ECVARMA.pvalues(input=input,T=T,k=k,d=d,del.r=as.vector(est_results[["del.r"]]),M=est_results[["M"]],bet=est_results[["bet"]],p=p,r=r,q=q,R=R,include=include,Sig=est_results[["Sig"]],u=est_results[["errors"]])},error=function(e){
      p.values<-list()
      #p.values$del.r<-"calculation failed, because no further iteration was possible"
      p.values$std.error<-"calculation failed, because no further iteration was possible"
      p.values$p.values<-"calculation failed, because no further iteration was possible"
      return(p.values)})
#function(e){return(message="no p-values could be calculated, because no further iteration was possible")})
  }else{
    p.values<-list()
    #p.values$del.r<-"not calculated"
    p.values$std.error<-"not calculated"
    p.values$p.values<-"not calculated"
  }
  #output
  out<-list(
    output="estimated",
    model.specification=list(
      k = k,
      p = p,
      q = q,
      r = r,
      include = include
    ),
    A0 = est_results[["A0"]],
    alphas = est_results[["alp"]],
    betas = cbind(diag(r),est_results[["bet"]]),
    AR = est_results[["G"]],
    MA = est_results[["M"]],
    Sigma = est_results[["Sig"]],
    estimation = list(
      sample.size=T,
      specification.scheme = "-",
      R = R,
      R0 = R0,
      del.r = as.vector(est_results[["del.r"]]),
      del.r2 = p.values$del.r,
      std.errors = p.values[["std.error"]],
      p.values = p.values[["p.values"]],
      error.terms = est_results[["errors"]],
      initital.values = initial_values,
      interations = "no iteration values with nleqslv"
    )
  )
  if (Newton) {
    out[["estimation"]][["interations"]] = est_results[["iterations"]]
  }

  class(out)<-"ECVARMA"
  attr(out,"determination")<-"estimated"
  return(out)
  # results<-list()
  #
  # results$est_results<-est_results
  # results$p.values<-p.values
  #
  # return(results)
}

#ecvarma.re
