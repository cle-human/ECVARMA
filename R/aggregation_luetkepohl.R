#' Lütkepohl's method for temporal aggregation
#' @description A function to derive the temporally aggregated model using Lütkepohl's method. The aggregated model is in final equation form.
#'
#' @param alphas a matrix of adjustment parameters (k x r)
#' @param betas a matrix of cointegration parameters (r x k) including betas that are normalized to 1 or 0
#' @param AR a matrix of the parameters before p lagged first differences (k x kp)
#' @param l the aggregation period length
#' @param f the used aggregation method: average aggregation - avg, beginning of period sampling - bop, end of period sampling - eop, or an aggregation matrix f (k x kl)
#' @param sig the non-aggregated variance-covariance matrix of the error terms
#' @return agg.coeff is a list of the aggregated AR parameters, MA parameters and the variance-covariance matrix of the aggregated error terms
#' @export

VECMagg.L<-function(alphas,betas,AR=NULL,l,f=c("avg","bop","eop","specify own F-matrix"),sig=NULL){
  k <- dim(alphas)[1]
  if(is.null(sig) == T){sig<-diag(1,k)}
  if(is.null(AR)){lag<-0}else{lag <- dim(AR)[2]/k} #p-1 lags
  r <- dim(alphas)[2]

  if (is.matrix(f)) {
    Fmat<-f
  }else{
    Fmat<-switch(f,
                 avg = matrix(rep(diag(1,k),times=l),k,k*l)*(1/l),
                 bop = cbind(diag(1,k),matrix(0,k,k*(l-1))),
                 eop = cbind(matrix(0,k,k*(l-1)),diag(1,k)) )
  }

  pi.mat <- alphas %*% t(betas)
  Amat <- matrix(NA, nrow = k, ncol = (lag+1) * (k))
  if (is.null(AR) == F) {
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

  if (l==1) {
    AR.agg<-polyMatrix::polyMatrix(cbind(diag(1,k),-Amat),k,k,lag+1)

    dot.AR<-matrix(0,k,k*lag)
    for (i in 1:lag) {
      for (j in (i+1):(lag+1)) {
        dot.AR[,(1:k)+(i-1)*k]<-dot.AR[,(1:k)+(i-1)*k]-AR.agg@coef[,(1:k)+j*k]
      }
    }

    dot.pi<- (diag(1,k))
    for (i in 1:(lag+1)) {
      dot.pi<-dot.pi+as.matrix(AR.agg@coef[,(i*k+1):(i*k+k)])#-
    }
    dot.pi<- -dot.pi

    dot.ec<- dot.pi[,1]

    sig.agg<-sig
    agg.coeff<-list(AR.agg,dot.pi,dot.AR,sig.agg)

  }else{
    #lag<-dim(Amat)[2]/k-1
    lag.fra<-ceiling((lag+1)/l)
    #Amat

    colmatlhs<-diag(k*l)
    for (i in 1:l){#rows
      for (j in 1:l) {#columns
        if ((i>j)&((lag+1)>=(i-j))){
          for (ii in 1:k) {
            for (jj in 1:k) {
              colmatlhs[(k*(i-1)+1)+ii-1,(k*(j-1)+1)+jj-1] <- -Amat[ii,(i-j-1)*k+jj]
            }
          }
        }
      }
    }

    colmatrhs<-matrix(0,k*l,k*l*lag.fra) #lag matrizenreihenfolge im framework getauscht
    #testp<-polyMatrix::polyMatrix(1,k*l,k*l*lag.fra)
    for (i in 1:l){#rows
      for (j in 1:(l*lag.fra)) {#columns
        for (ii in 1:k) {
          for (jj in 1:k) {
            if ((i+(l*lag.fra-j))<=(lag+1)){
              colmatrhs[(k*(i-1))+ii,(k*(j-1))+jj] <- Amat[ii,((l*lag.fra-j)+i-1)*k+jj]
            }
          }
        }
      }
    }
    #Inverse links multipliziert
    colmatma<-solve(colmatlhs)%*%diag(k*l)
    colmatrhs<-solve(colmatlhs)%*%colmatrhs
    #Lag operator
    #x<-polyMatrix::polyMatrix(c(0,1),1,1,1)
    colmatrhs<-polyMatrix::polyMatrix(colmatrhs,k*l,k*l*lag.fra,0)
    for (ii in 1:(k*l)) {
      for (jj in 1:(k*l*lag.fra)) {
        i<-lag.fra-ceiling(jj/k/l)+1
        colmatrhs[ii,jj]<-colmatrhs[ii,jj]*polyMatrix::polyMatrix(c(seq(0,0,length.out=i),1),1,1,i)

      }
    }
    AR.fra<-diag(k*l)
    for (i in 1:lag.fra) {
      AR.fra<-AR.fra-colmatrhs[1:(k*l),(1:(k*l))+(i-1)*k*l]
    }

    # cfunc<-function(x){
    #   y<-numeric(k*(k*l-k)+1)
    #   c<-matrix(x,k*l,k*l-k) #ortho
    #   y[1:(k*(k*l-k))]<-as.vector((Fmat)%*%c)#k x kl   kl x kl-k
    #
    #   qr(c)$rank-k*l-k->y[k*(k*l-k)+1]
    #   y
    # }
    # c.start<-numeric(k*l*(k*l-k))
    # c.sol<-nleqslv::nleqslv(c.start,cfunc)

    C<-matrix(0,k*l-k,k*l)
    n<-1
    for (i in 1:(k*l)) {
      if (all(Fmat[,i]==0)) {
        C[n,i]<-1
        n<-n+1
      }
    }
    n<-1
    F.rearrage<-rbind((Fmat),C)
    for (i in 1:(k*l)) {
      if (all(F.rearrage[i,]==0)) {
        F.rearrage[i,n]<-1
        n<-n+1
      }
    }

    adj<-function(x){
      n <- nrow(x)
      if (nrow(x) != ncol(x)) {
        stop("Square matrix is expected")
      }
      y<-polyMatrix::polyMatrix(0,n,n)
      for (i in 1:n) {
        for (j in 1:n) {
          y[i,j]<- (-1)^(i+j)*det(polyMatrix::polyMatrix(x[-i,-j],n-1,n-1))
        }
      }
      return(y)
    }

    AR.fra.re<-F.rearrage%*%AR.fra%*%solve(F.rearrage)
    AR.agg<-det(polyMatrix::polyMatrix(AR.fra.re[(k+1):(k*l),(k+1):(k*l)],k*l-k,k*l-k))*AR.fra.re[(1):(k),(1):(k)]-AR.fra.re[(1):(k),(k+1):(k*l)]%*%adj(polyMatrix::polyMatrix(AR.fra.re[(k+1):(k*l),(k+1):(k*l)],k*l-k,k*l-k))%*%AR.fra.re[(k+1):(k*l),(1):(k)]

    # dot.AR<-matrix(0,k,k*lag)
    # for (i in 1:lag) {
    #   for (j in (i+1):(lag+1)) {
    #     dot.AR[,(1:k)+(i-1)*k]<-dot.AR[,(1:k)+(i-1)*k]-AR.agg@coef[,(1:k)+j*k]
    #   }
    # }
    #
    # dot.pi<- (diag(1,k))
    # for (i in 1:(lag+1)) {
    #   dot.pi<-dot.pi+as.matrix(AR.agg@coef[,(i*k+1):(i*k+k)])#-
    # }
    # dot.pi<- -dot.pi
    #
    # dot.ec<- dot.pi[,1]


    MA.fra.re<-F.rearrage%*%colmatma#%*%solve(F.rearrage)

    MA.fra.u1<-det(polyMatrix::polyMatrix(AR.fra.re[(k+1):(k*l),(k+1):(k*l)],k*l-k,k*l-k))*polyMatrix::polyMatrix(MA.fra.re[(1):(k),(1):(k)],k,k)-AR.fra.re[(1):(k),(k+1):(k*l)]%*%adj(polyMatrix::polyMatrix(AR.fra.re[(k+1):(k*l),(k+1):(k*l)],k*l-k,k*l-k))%*%MA.fra.re[(k+1):(k*l),(1):(k)]
    MA.fra.u2<-det(polyMatrix::polyMatrix(AR.fra.re[(k+1):(k*l),(k+1):(k*l)],k*l-k,k*l-k))*polyMatrix::polyMatrix(MA.fra.re[(1):(k),(k+1):(k*l)],k,k*l-k)-AR.fra.re[(1):(k),(k+1):(k*l)]%*%adj(polyMatrix::polyMatrix(AR.fra.re[(k+1):(k*l),(k+1):(k*l)],k*l-k,k*l-k))%*%MA.fra.re[(k+1):(k*l),(k+1):(k*l)]

    MA.fra<-cbind(MA.fra.u1,MA.fra.u2)

    #bereinige auf extrem kleine Lags
    for (i in 1:length(MA.fra@coef)) {
      if(abs(MA.fra@coef[i])<1e-6){MA.fra@coef[i]<-0}
    }

    #Lag-Operator ausschreiben
    MA.fra.no.lag.notation<-MA.fra@coef[,1:(l*k)]
    if (degree(MA.fra)>0) {
      for (i in 1:degree(MA.fra)) {
        MA.fra.no.lag.notation<-cbind(MA.fra@coef[,1:(l*k)+i*l*k],MA.fra.no.lag.notation)
      }
    }

    #?bersch?ssige NULL-Matrizen elimienieren
    for (i in 1:((degree(MA.fra)+1)*l)) {
      if(all(MA.fra.no.lag.notation[,1:k+(i-1)*k]==0)){MA.fra.no.lag.notation<-MA.fra.no.lag.notation[,-(1:k+(i-1)*k)]}else{break}
    }
    for (i in (dim(MA.fra.no.lag.notation)[2]/k):1) {
      if(all(MA.fra.no.lag.notation[,1:k+(i-1)*k]==0)){MA.fra.no.lag.notation<-MA.fra.no.lag.notation[,-(1:k+(i-1)*k)]}else{break}
    }
    #Lag der MA bestimmen
    Lag.MA.agg<-ceiling(dim(MA.fra.no.lag.notation)[2]/k/l)-1

    #aggregierte error cov
    sig.fra<-matrix(0,k,k)
    for (i in 0:(dim(MA.fra.no.lag.notation)[2]/k-1)) {
      sig.fra<-sig.fra+MA.fra.no.lag.notation[,(1:k)+i*k]%*%sig%*%t(MA.fra.no.lag.notation[,(1:k)+i*k])
    }


    if (Lag.MA.agg>0) {
      #error auto cov
      autocov<-matrix(0,k,k*Lag.MA.agg)
      MA.fra.cov<-cbind(matrix(0,k,(Lag.MA.agg+1)*k*l-dim(MA.fra.no.lag.notation)[2]+Lag.MA.agg*k*l),MA.fra.no.lag.notation)
      for (i in 1:Lag.MA.agg) {
        MA.fra.cov<-rbind(MA.fra.cov,cbind(matrix(0,k,dim(MA.fra.cov)[2]-dim(MA.fra.no.lag.notation)[2]-i*k*l),MA.fra.no.lag.notation,matrix(0,k,i*k*l)))
        for (j in (l*i):(dim(MA.fra.cov)[2]/k-1-l*i)) {
          autocov[,1:k+(k*(i-1))]<-autocov[,1:k+(k*(i-1))]+MA.fra.cov[1:k,(1:k)+j*k]%*%sig%*%t(MA.fra.cov[1:k+i*k,(1:k)+j*k])
        }
      }
      MA.function<-function(x){  #nur sig wird optimiert MA nicht?
        y<-matrix(0,k,k+k*Lag.MA.agg)
        sig.agg<-matrix(x[1:k^2],k,k)
        MA.cof.agg<-matrix(x[(k^2+1):length(x)],k,k*Lag.MA.agg)
        #error Cov
        y[,1:k]<-sig.agg
        for (i in 1:Lag.MA.agg-1) {
          y[,1:k]<-y[,1:k]+MA.cof.agg[,1:k+k*i]%*%sig.agg%*%t(MA.cof.agg[,1:k+k*i])
        }
        y[,1:k]<-y[,1:k]-sig.fra
        #error autocovs
        autocov.agg<-matrix(0,k,k*Lag.MA.agg)
        MA.fra.cov.agg<-cbind(diag(1,k,k),MA.cof.agg,matrix(0,k,k*Lag.MA.agg))
        for (i in 1:Lag.MA.agg) {
          MA.fra.cov.agg<-rbind(MA.fra.cov.agg,cbind(matrix(0,k,i*k),diag(1,k),MA.cof.agg,matrix(0,k,(Lag.MA.agg-i)*k)))
          for (j in (i):(dim(MA.fra.cov.agg)[2]/k-1-i)) {
            autocov.agg[,1:k+(k*(i-1))]<-autocov.agg[,1:k+(k*(i-1))]+MA.fra.cov.agg[1:k,(1:k)+j*k]%*%sig.agg%*%t(MA.fra.cov.agg[1:k+i*k,(1:k)+j*k])
          }
        }
        y[,(k+1):(k+k*Lag.MA.agg)]<-autocov.agg-autocov
        y
      }

      xstart<-c(l*sig,matrix(0,k,k*Lag.MA.agg))
      #xstart<-c(sig.fra,autocov)
      #x<-c(sig.avg.mean[l,],MA.cof.avg.mean[l,])
      MA.sol<-nleqslv::nleqslv(xstart,MA.function)
      sig.agg<-matrix(MA.sol[["x"]][1:k^2],k,k)
      MA.agg<-matrix(MA.sol[["x"]][(k^2+1):(k^2*Lag.MA.agg+k^2)],k,k*Lag.MA.agg)
      agg.coeff<-list(AR.agg,MA.agg,sig.agg)#dot.pi,dot.AR,
    }else{
      sig.agg<-sig.fra
      agg.coeff<-list(AR.agg,sig.agg)#dot.pi,dot.AR,
    }
  }
  out<-list(
    output="aggregated",
    model.specification=list(
      k = k,
      p = (lag+1)*k,
      q = Lag.MA.agg,
      r = 0,
      include = "none"
    ),
    alphas = NULL,
    betas = NULL,
    AR = NULL,
    MA = MA.agg,
    Sigma = sig.agg,
    aggregation = list(
      l = l,
      scheme = f,
      method = "Lütkepohl",
      VAR.rep = AR.agg
    )
  )
  class(out)<-"ECVARMA"
  attr(out,"determination")<-"aggregated"
  return(out)
}
