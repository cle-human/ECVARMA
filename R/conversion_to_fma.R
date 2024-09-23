#' Conversion of ECVARMA to FMA specification
#' @description
#' Conversion of an ECVARMA model to its FMA specification
#' @param ecvarma an ECVARMA object, created by ecvarma.define, an ECVARMA aggregation function or an estimated ECVARMA process
#' @return an ECVARMA object
#' @export
ecvarma.toFMA<-function(ecvarma){ #para=result of theoretical aggregation function
  if (attr(ecvarma,"determination")=="converted to FMA") {
    stop("ECVARMA is already in FMA form.")
  }
  if (attr(ecvarma,"determination")=="estimated"&ecvarma[["estimation"]][["specification.scheme"]]=="final moving average") {
    stop("ECVARMA is already in FMA form.")
  }
  k<-ecvarma[["model.specification"]][["k"]]
  p<-ecvarma[["model.specification"]][["p"]]
  q<-ecvarma[["model.specification"]][["q"]]
  r<-ecvarma[["model.specification"]][["r"]]
  alphas<-ecvarma$alphas
  betas<-ecvarma$betas
  AR<-ecvarma$AR
  MA<-ecvarma$MA
  sig<-ecvarma$Sigma
  include<-ecvarma[["model.specification"]][["include"]]
  A0<-ecvarma$A0
  if (!all(A0==diag(k))) {
    A0i<-solve(A0)
    alphas<-A0i%*%alphas
    AR<-A0i%*%AR
    MA<-A0i%*%MA
  }
  if (include=="LRconst") {
    betas<-betas[,-(k+1)]
  }
  if (include=="SRconst") {
    SRconst<-AR[,1]
    AR<-AR[,-1]
  }
  if(is.null(sig)){sig<-diag(1,k)}
  #if(is.null(AR)){p<-0}else{p <- dim(AR)[2]/k} #p-1 lags
  p.AR<-p+1
  #if(is.null(MA)){q<-0}else{q <- dim(MA)[2]/k}

  #transform VECM to level VAR
  pi.mat <- alphas %*% (betas)
  Amat <- matrix(NA, nrow = k, ncol = (p.AR) * (k))
  if (p!=0) {
    Amat[, 1:k] <- pi.mat + diag(k) + AR[,(1:k)]
    if (p > 1){
      for (i in 1:(p - 1)){
        Amat[, (1:k) + k * i] <- AR[,(1:k) + k * (i)]-AR[,(1:k) + k * (i-1)]
      }
    }
    Amat[,(1:k) + k * p] <- -AR[,(1:k) + k * (p-1)]
  }else {
    Amat[, 1:k] <- pi.mat + diag(k)
  }
  #y<-para
  #k<-nrow(para[[4]])
  B_L<-polyMatrix(c(diag(1,k),MA),k,k,q)
  B_L.adj<-adjoi(B_L)
  MA.fma<-t(as.vector(polyMatrix::det(B_L))[-1])%x%diag(1,k)
  pureAR.fma<-B_L.adj%*%polyMatrix(c(diag(1,k),-Amat),k,k,p+1)
  pureAR.fma<- -pureAR.fma@coef[,-(1:k)]
  p.fma<-ncol(pureAR.fma)/k-1

  AR.fma<-matrix(0,k,k)
  if (p.fma!=0) {
    AR.fma<-matrix(0,k,k*p.fma)
    for (i in 1:p.fma) {
      for (j in (i):(p.fma)) {
        AR.fma[,(1:k)+(i-1)*k]<-AR.fma[,(1:k)+(i-1)*k]+pureAR.fma[,(1:k)+j*k]
      }
    }
    AR.fma<- -AR.fma
  }

  pi.fma<- (diag(1,k))
  for (i in (1:(p.fma+1))-1) {
    pi.fma<-pi.fma-as.matrix(pureAR.fma[,(i*k+1):(i*k+k)])#-
  }
  pi.fma<- -pi.fma
  alpha.fma<-pi.fma[,1:r]
  if (include=="SRconst") {
    SRconst.fma<-B_L%*%SRconst
    SRconst.fma<-colSums(SRconst.fma@coef)
    AR.fma = cbind(SRconst.fma,AR.fma)
  }
  #results<-list(AR.fma,pi.fma)
  #return(results)

  #pi.fma<-ARtoEC(pureAR.fma,k,ncol(pureAR.fma)/k-1)[[2]]
  #AR.fma<-ARtoEC(pureAR.fma,k,ncol(pureAR.fma)/k-1)[[1]]
  out<-list(
    output="respecified",
    model.specification=list(
      k = k,
      p = p.fma,
      q = q*k,
      r = r,
      include = include
    ),
    alphas = alpha.fma,
    betas = ecvarma$betas,
    AR = AR.fma,
    MA = MA.fma,
    Sigma = sig,
    conversion = list(
      pi = pi.fma,
      pureAR.fma = pureAR.fma
    )
  )
  class(out)<-"ECVARMA"
  attr(out,"determination")<-"converted to FMA"
  return(out)
}
#
# ARtoEC<-function(AR.agg,k,lag){
#   lag.AR<-lag+1
#   dot.AR<-matrix(0,k,k)
#   if (lag!=0) {
#     dot.AR<-matrix(0,k,k*lag)
#     for (i in 1:lag) {
#       for (j in (i):(lag.AR-1)) {
#         dot.AR[,(1:k)+(i-1)*k]<-dot.AR[,(1:k)+(i-1)*k]-AR.agg[,(1:k)+j*k]
#       }
#     }
#   }
#
#   dot.pi<- (diag(1,k))
#   for (i in (1:(lag.AR))-1) {
#     dot.pi<-dot.pi+as.matrix(AR.agg[,(i*k+1):(i*k+k)])#-
#   }
#   dot.pi<- -dot.pi
#   results<-list(dot.AR,dot.pi)
#   return(results)
# }

determinant.poly<-function(x){
  if (is.polynomial(x)) {
    return(x)
  }else if(is.polyMatrix(x)){
    y<-polyMatrix::det(x)
    return(y)
  }else{stop("no polynom")}
}

adjoi<-function(x){
  y<-x
  for(i in 1:nrow(x)){
    for (j in 1:ncol(x)) {
      y[i,j]<-determinant.poly(x[-i,-j])*(-1)^(i+j)
    }
  }
  y<-polyMatrix::t(y)
  return(y)
}
