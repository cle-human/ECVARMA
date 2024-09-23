#' Method for contemporaneous aggregation
#' @description A function to derive the contemporaneously aggregated error correction model.
#'
#' @param ecvarma an ECVARMA object, created by ecvarma.define, an ECVARMA aggregation function or an estimated ECVARMA process
#' @param l the number of omitted variables
#' @param Fmat a matrix to change order of variables or create contemporaneous averages, must have full rank (k x k)
#' @return an ECVARMA object of the aggregated model
#' @export
#' @import polyMatrix
ecvarma.agg.con<-function(ecvarma,l,Fmat=NULL){

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
    LRconst<-betas[,k+1]
    betas<-betas[,-(k+1)]
  }else{
    LRconst<-NULL
  }
  if (include=="SRconst") {
    SRconst<-AR[,1]
    AR<-AR[,-1]
  }else{
    SRconst<-NULL
  }

  k-l->m
  if (m<1)stop("l cannot be larger than k")
  if(is.null(sig)){sig<-diag(1,k)}
  if(is.null(AR)){lag<-0}else{lag <- dim(AR)[2]/k} #p-1 lags
  lag.AR<-lag+1
  if(is.null(MA)){lag.MA<-0}else{lag.MA <- dim(MA)[2]/k}
  if(is.null(Fmat)){Fmat<-diag(1,k)}
  if(qr(Fmat)$rank<k)stop("F-Matrix must be of full rank.")
  lag.b<-ceiling(l*lag.AR/m)
  #g<-(lag.AR)
  #k<-l
  #s<-lag.MA
  #pi.mat <- matrix(round(rnorm(k^2,0,0.1),2),k,k)
  pi.mat <- alphas %*% betas
  Amat <- matrix(NA, nrow = k, ncol = (lag.AR) * (k))
  if (is.null(AR) == FALSE) {
    Amat[, 1:k] <- pi.mat + diag(k) + AR[,(1:k)]
    if (lag > 1){
      for (i in 1:(lag - 1)){
        Amat[, (1:k) + k * i] <- AR[,(1:k) + k * (i)]-AR[,(1:k) + k * (i-1)]
      }
    }
    Amat[,(1:k) + k * lag] <- -AR[,(1:k) + k * (lag-1)]
  }else {
    Amat[, 1:k] <- pi.mat + diag(k)
  }

  A_L<-Fmat%mult%polyMatrix::polyMatrix(cbind(diag(1,k),-Amat),k,k,lag.AR)%mult%solve(Fmat)
  norm<-solve(A_L@coef[,1:k])
  if (is.null(MA)) {
    M_L<-Fmat%mult%polyMatrix::polyMatrix(diag(1,k),k,k,0)
  }else{
    M_L<-Fmat%mult%polyMatrix::polyMatrix(cbind(diag(1,k),MA),k,k,dim(MA)[2])
  }
  M_L<-norm%mult%M_L
  A_L<-norm%mult%A_L
  #AR.agg<-A_L
  #MA.agg<-M_L

  A11_L<-A_L[1:(m),1:(m)]     #mxm
  A12_L<-A_L[1:(m),(m+1):k]   #mxl
  A21_L<-A_L[(m+1):k,1:(m)]   #lxm
  A22_L<-A_L[(m+1):k,(m+1):k] #lxl

  lag.A12<-polyMatrix::degree(A12_L)
  lag.A22<-polyMatrix::degree(A22_L)

  #A11<-A11_L@coef
  A12<-A12_L@coef
  if (lag.A12<lag.AR) {
    A12<-cbind(A12,matrix(0,dim(A12)[1],dim(A12)[2]*(lag.AR-lag.A12)))
  }
  #A21<-A21_L@coef
  if (l==1) {
    A22<-matrix(A22_L,l,l*(lag.AR+1))
  }else{A22<-A22_L@coef}

  if (lag.b==1) {
    lhs<-matrix(0,l*lag.AR,m)
    for (i in 1:lag.AR) {
      lhs[1:l+l*(i-1),1:m]<-t(A12[,1:l+l*i])
    }
    lhs<-lhs%x%diag(1,m)
  }else{
    lhs<-matrix(0,l*(lag.AR+lag.b-1),m*(lag.b)+l*(lag.b-1)) #ein Zeilenrank zu wenig
    for (i in 1:(lag.AR+lag.b-1)) {
      for (j in 1:lag.b) {
        lhs[1:l+((i-1)*l),1:m+((j-1)*m)]<-t(
          if(1+l*i-l*(j-1)<1){
            matrix(0,m,l)}else{
              if(l+l*(i+lag.b-3)-l*(j-1)>(lag.AR+lag.b-1)*l-1){
                matrix(0,m,l)}else{
                  A12[,1:l+l*i-l*(j-1)]}
            })
      }
      for (j in 1:(lag.b-1)) {
        lhs[1:l+((i-1)*l),1:l+((j-1)*l+lag.b*m)]<-t(if(1+l*(i-1)-l*(j-1)<1){matrix(0,l,l)}else{if(l+l*(i+lag.b-3)-l*(j)>(lag.AR+lag.b-1)*l-1){matrix(0,l,l)}else{-A22[,1:l+l*(i-1)-l*(j-1)]}})
      }
    }
    lhs<-lhs%x%diag(1,m)
  }
  rhs<-matrix(0,m*l*(lag.AR+lag.b-1),1)
  for (i in 1:(lag.AR+lag.b-1)) {
    rhs[1:(m*l)+(i-1)*m*l,]<-as.vector((A12[,1:l+l]%*%as.matrix(if(lag.AR<i){matrix(0,l,l)}else{A22[,1:l+l*i]}))-if(lag.AR<(i+1)){matrix(0,m,l)}else{A12[,1:l+l*(i+1)]})
  }
  Bmatrices<-MASS::ginv(lhs)%*%rhs
  B1_L<-polyMatrix::polyMatrix(c(diag(1,m),Bmatrices[1:(m^2*lag.b)]),m,m,lag.b)
  A12.1<-A12[,(l)+(1:(l))]
  if (lag.b==1) {
    B2_L<-polyMatrix::polyMatrix(c(matrix(0,m,l),A12.1),m,l,lag.b)
  }else{
    B2_L<-polyMatrix::polyMatrix(c(matrix(0,m,l),A12.1,Bmatrices[(m^2*lag.b+1):(m^2*lag.b+m*l*(lag.b-1))]),m,l,lag.b)
  }

  if (l==1) {
    AR.agg<-B1_L%mult%A11_L-cbind(B2_L,0)%mult%rbind(A21_L,0)
  }else{
    AR.agg<-B1_L%mult%A11_L-B2_L%mult%A21_L
  }

  AR.agg<-AR.agg[-(k+1-(l:1)),-(k+1-(l:1))]
  AR.agg<--AR.agg@coef[,-(1:m)]
  lag.AR<-dim(AR.agg)[2]/m
  lag<-lag.AR-1
  dot.AR<-matrix(0,m,m)
  if (lag!=0) {
    dot.AR<-matrix(0,m,m*lag)
    for (i in 1:lag) {
      for (j in (i):(lag.AR-1)) {
        dot.AR[,(1:m)+(i-1)*m]<-dot.AR[,(1:m)+(i-1)*m]-AR.agg[,(1:m)+j*m]
      }
    }
  }

  dot.pi<- (diag(1,m))
  for (i in (1:(lag.AR))-1) {
    dot.pi<-dot.pi-as.matrix(AR.agg[,(i*m+1):(i*m+m)])#-
  }
  dot.pi<- -dot.pi

  if (is.null(SRconst)) {
    SRconst<-matrix(0,k,1)
  }
  SRconst<-Fmat%*%SRconst
  SR.const.agg.l<-B1_L%mult%matrix(SRconst[1:m],m,1)-B2_L%mult%matrix(SRconst[(m+1):k],l,1)
  SR.const.agg<-matrix(0,m,1)
  for (i in 1:max(1,polyMatrix::degree(SR.const.agg.l)+1)) {
    SR.const.agg<-SR.const.agg+SR.const.agg.l@coef[,i]
  }
  if (is.null(LRconst)) {
    LRconst<-matrix(0,r,1)
  }
  Pi.const<-alphas%*%matrix(LRconst,r,1)
  Pi.const<-Fmat%*%Pi.const
  LR.const.agg.l<-B1_L%mult%matrix(Pi.const[1:m,],m,1)-B2_L%mult%matrix(Pi.const[(m+1):k,1],l,1)
  LR.const.agg<-matrix(0,m,1)
  for (i in 1:max(1,polyMatrix::degree(LR.const.agg.l)+1)) {
    LR.const.agg<-LR.const.agg+LR.const.agg.l@coef[,i]
  }


  #MA
  M11_L<-M_L[1:m,1:m]
  M12_L<-M_L[1:m,(m+1):k]
  M21_L<-M_L[(m+1):k,1:m]
  M22_L<-M_L[(m+1):k,(m+1):k]
  MA.fra<-cbind(B1_L%mult%M11_L-B2_L%mult%M21_L,B1_L%mult%M12_L-B2_L%mult%M22_L)

  N<-MA.fra@coef
  lag.N<-dim(N)[2]/k-1

  #aggregierte error cov
  sig.fra<-matrix(0,m,m)
  for (i in 0:(dim(N)[2]/k-1)){#(lag.AR*l-lag.AR+lag.MA)) {
    sig.fra<-sig.fra+N[,(1:k)+i*k]%*%sig%*%t(N[,(1:k)+i*k])
  }

  if (lag.N>0) {
    #error auto cov
    autocov<-matrix(0,m,m*lag.N)

    N.fra.cov<-cbind(N,matrix(0,m,k*lag.N))
    for (i in 1:lag.N) {
      N.fra.cov<-rbind(N.fra.cov,cbind(matrix(0,m,i*k),N,matrix(0,m,(lag.N-i)*k)))
      for (j in (i):(dim(N.fra.cov)[2]/k-1-i)) {
        autocov[,1:m+(m*(i-1))]<-autocov[,1:m+(m*(i-1))]+N.fra.cov[1:m,(1:k)+j*k]%*%sig%*%t(N.fra.cov[1:m+i*m,(1:k)+j*k])
      }
    }

    for (i in lag.N:1) {
      if (all(autocov[,1:m+(i-1)*m]==0)) {
        autocov<-autocov[,-(1:m+(i-1)*m)]
        lag.N<-lag.N-1
      }else{break}
    }
    MA.function<-function(x){
      y<-matrix(0,m,m+m*lag.N)
      sig.agg<-matrix(x[1:m^2],m,m)
      MA.cof.agg<-matrix(x[(m^2+1):length(x)],m,m*lag.N)
      #error Cov
      y[,1:m]<-sig.agg
      for (i in 1:lag.N-1) {
        y[,1:m]<-y[,1:m]+MA.cof.agg[,1:m+m*i]%*%sig.agg%*%t(MA.cof.agg[,1:m+m*i])
      }
      y[,1:m]<-y[,1:m]-sig.fra
      #error autocovs
      autocov.agg<-matrix(0,m,m*lag.N)
      MA.fra.cov.agg<-cbind(diag(1,m,m),MA.cof.agg,matrix(0,m,m*lag.N))
      for (i in 1:lag.N) {
        MA.fra.cov.agg<-rbind(MA.fra.cov.agg,cbind(matrix(0,m,i*m),diag(1,m),MA.cof.agg,matrix(0,m,(lag.N-i)*m)))
        for (j in (i):(dim(MA.fra.cov.agg)[2]/m-1-i)) {
          autocov.agg[,1:m+(m*(i-1))]<-autocov.agg[,1:m+(m*(i-1))]+MA.fra.cov.agg[1:m,(1:m)+j*m]%*%sig.agg%*%t(MA.fra.cov.agg[1:m+i*m,(1:m)+j*m])
        }
      }
      y[,(m+1):(m+m*lag.N)]<-autocov.agg-autocov
      y
    }
    xstart<-c(sig[1:m,1:m],matrix(0,m,m*lag.N))
    MA.sol<-nleqslv::nleqslv(xstart,MA.function,control = list(maxit=200))
    sig.agg<-matrix(MA.sol[["x"]][1:m^2],m,m)
    MA.agg<-matrix(MA.sol[["x"]][(m^2+1):(m^2*lag.N+m^2)],m,m*lag.N)
  }else{
    sig.agg<-sig.fra
    MA.agg<-matrix(0,m,m)
  }
  r<-qr(dot.pi)$rank
  if (all(SR.const.agg!=0)) {
    include<-"unrestricted constant"
    AR = cbind(SR.const.agg,dot.AR)
    betas = solve(dot.pi[1:r,1:r])%*%dot.pi[1:r,] #r x (k-r)
  }else if (all(LR.const.agg!=0)) {
    include<-"restricted constant"
    AR = dot.AR
    betas = solve(dot.pi[1:r,1:r])%*%cbind(matrix(dot.pi[1:r,],r,m),LR.const.agg[1:r,]) #r x (k-r)
  }else{
    include <- "none"
    AR = dot.AR
    betas = solve(dot.pi[1:r,1:r])%*%dot.pi[1:r,]
  }

  out<-list(
    output="aggregated",
    model.specification=list(
      k = k,
      p = lag.AR,
      q = lag.N,
      r = r,
      include = include
    ),
    alphas = dot.pi[,1:r],
    betas = betas,
    AR = AR,
    MA = MA.agg,
    Sigma = sig.agg,
    aggregation = list(
      l = l,
      Fmat = Fmat,
      method = "Hoffmann",
      VAR.rep = AR.agg,
      pi.agg = dot.pi,
      B1_L = B1_L,
      B2_L = B2_L
    )
  )
  class(out)<-"ECVARMA"
  attr(out,"determination")<-"aggregated"
  return(out)
}

