#' temporal aggregation with Marcelino's method (1999)
#' @description A function to derive the temporally aggregated error correction model.
#'
#' @param ecvarma an ECVARMA object, created by ecvarma.define, an ECVARMA aggregation function or an estimated ECVARMA process
#' @param l the aggregation period length, i.e. if the original frequency is weekly and l=2, the new frequency will be biweekly
#' @param f the used aggregation method: average aggregation - avg, beginning of period sampling - bop
#' @param start.values for the Newton method to find the aggregated error variance-covariance matrix and MA parameter (k x k+k*q_agg), with "auto" these are generated automatically
#' @return an ECVARMA object of the aggregated model
#' @export

ecvarma.agg.M<-function(ecvarma,l,f=c("avg","bop"),start.values="auto"){
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
  if(is.null(AR)){lag<-0}else{lag <- dim(AR)[2]/k} #p-1 lags
  lag.AR<-lag+1
  if(is.null(MA)){lag.MA<-0}else{lag.MA <- dim(MA)[2]/k}

  #transform VECM to level VAR
  pi.mat <- alphas %*% (betas)
  Amat <- matrix(NA, nrow = k, ncol = (lag.AR) * (k))
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

  #aggregate level VAR
  Gv<-cbind(--Amat,matrix(0,k,(lag.AR)*l*k-(lag.AR)*k))

  Gm<-matrix(0,(lag.AR)*l*k-(lag.AR)*k,(lag.AR)*l*k)
  n<-0
  for (i in 1:((lag.AR)*l-(lag.AR))) {
    Gm[(1:k)+(i-1)*k,]<-cbind(matrix(0,k,n*k),-diag(1,k),--Amat,matrix(0,k,(lag.AR)*l*k-n*k-k-dim(Amat)[2]))
    n<-n+1
  }
  Gvmin<-Gv
  Gmmin<-Gm
  for (i in lag.AR:1) {
    Gvmin<-Gvmin[,-((1:k)+k*i*l-k)]
    Gmmin<-Gmmin[,-((1:k)+k*i*l-k)]
  }

  AR.agg<-(Gvmin%*%solve(Gmmin)%*%Gm-Gv)
  #n<-0
  for (n in 1:lag.AR-1) {
    AR.agg<-AR.agg[,-(1:((l-1)*k)+n*k)]
    #n<-n+1
  }
  AR.agg<- -1*AR.agg
  A_L<-polyMatrix::polyMatrix(cbind(diag(1,k),-Amat),k,k,dim(Amat)[2]/k)
  #transform level VAR to VECM
  dot.AR<-matrix(0,k,k)
  if (lag!=0) {
    dot.AR<-matrix(0,k,k*lag)
    for (i in 1:lag) {
      for (j in (i):(lag.AR-1)) {
        dot.AR[,(1:k)+(i-1)*k]<-dot.AR[,(1:k)+(i-1)*k]-AR.agg[,(1:k)+j*k]
      }
    }
  }

  dot.pi<- (diag(1,k))
  for (i in (1:(lag.AR))-1) {
    dot.pi<-dot.pi-as.matrix(AR.agg[,(i*k+1):(i*k+k)])
  }
  dot.pi<- -dot.pi

  #derive MA dynamic
  Bv<-Gvmin%*%solve(Gmmin)
  B_L<-polyMatrix::polyMatrix(cbind(diag(1,k),-Bv),k,k,dim(Bv)[2]/k)

  if(is.null(MA)){
    S_L<-diag(1,k)
  }else{
    S_L<-polyMatrix::polyMatrix(cbind(diag(1,k),MA),k,k,dim(MA)[2]/k)
  }
  if(f=="avg"){
    S_L<-S_L%mult%polyMatrix::polyMatrix(diag(1/l,k),k,k,l-1)
    lag.MA<-lag.MA+l-1
  }

  N_L<-B_L%mult%S_L

  N<-N_L@coef
  lag.N<-dim(N)[2]/k-1

  #MA lag agg
  if(f=="avg"){
    if (lag.AR>lag.MA) {
      q<-0
      co<-0
      while (co<2) {
        co<-0
        if (q*l<lag.AR-lag.MA+1) {
          co<-co+1
        }else{q<-q+1}
        if (lag.AR-lag.MA+1<=(q+1)*l) {
          co<-co+1
        }else{q<-q+1}
      }
      lag.MA.agg<-lag.AR-q
    }
    if (lag.AR==lag.MA) {
      lag.MA.agg<-lag.AR
    }
    if (lag.AR<lag.MA) {
      q<-0
      co<-0
      while (co<2) {
        co<-0
        if (q*l<=lag.MA-lag.AR-1) {
          co<-co+1
        }else{q<-q+1}
        if (lag.MA-lag.AR-1<(q+1)*l) {
          co<-co+1
        }else{q<-q+1}
      }
      lag.MA.agg<-lag.AR+q+1
    }
  }else{
    if (lag.AR>lag.MA) {
      q<-0
      co<-0
      while (co<2) {
        co<-0
        if (q*l<lag.AR-lag.MA) {
          co<-co+1
        }else{q<-q+1}
        if (lag.AR-lag.MA<=(q+1)*l) {
          co<-co+1
        }else{q<-q+1}
      }
      lag.MA.agg<-lag.AR-1-q
    }
    if (lag.AR==lag.MA) {
      lag.MA.agg<-lag.AR
    }
    if (lag.AR<lag.MA) {
      q<-0
      co<-0
      while (co<2) {
        co<-0
        if (q*l<=lag.MA-lag.AR) {
          co<-co+1
        }else{q<-q+1}
        if (lag.MA-lag.AR<(q+1)*l) {
          co<-co+1
        }else{q<-q+1}
      }
      lag.MA.agg<-lag.AR+q
    }
  }

  #aggregated error cov
  sig.fra<-matrix(0,k,k)
  for (i in 0:(dim(N)[2]/k-1)){#(lag.AR*l-lag.AR+lag.MA)) {
    sig.fra<-sig.fra+N[,(1:k)+i*k]%*%sig%*%t(N[,(1:k)+i*k])
  }

  if (lag.MA.agg>0) {
    #error auto cov
    autocov<-matrix(0,k,k*lag.MA.agg)

    N.fra.cov<-cbind(N,matrix(0,k,k*lag.MA.agg*l))
    for (i in 1:lag.MA.agg) {
      N.fra.cov<-rbind(N.fra.cov,cbind(matrix(0,k,i*k*l),N,matrix(0,k,(lag.MA.agg-i)*k*l)))
      for (j in (l*i):(dim(N.fra.cov)[2]/k-1-i*l)) {
        autocov[,1:k+(k*(i-1))]<-autocov[,1:k+(k*(i-1))]+N.fra.cov[1:k,(1:k)+j*k]%*%sig%*%t(N.fra.cov[1:k+i*k,(1:k)+j*k])
      }
    }
    for (i in lag.MA.agg:1) {
      if (all(autocov[,1:k+(i-1)*k]==0)) {
        autocov<-autocov[,-(1:k+(i-1)*k)]
        lag.MA.agg<-lag.MA.agg-1
      }else{break}
    }

    #derive new aggregated cov and MA dynamic
    MA.function<-function(x){
      y<-matrix(0,k,k+k*lag.MA.agg)
      sig.agg<-matrix(x[1:k^2],k,k)
      MA.cof.agg<-matrix(x[(k^2+1):length(x)],k,k*lag.MA.agg)
      #error Cov
      y[,1:k]<-sig.agg
      for (i in 1:lag.MA.agg-1) {
        y[,1:k]<-y[,1:k]+MA.cof.agg[,1:k+k*i]%*%sig.agg%*%t(MA.cof.agg[,1:k+k*i])
      }
      y[,1:k]<-y[,1:k]-sig.fra
      #error autocovs
      autocov.agg<-matrix(0,k,k*lag.MA.agg)
      MA.fra.cov.agg<-cbind(diag(1,k,k),MA.cof.agg,matrix(0,k,k*lag.MA.agg))
      for (i in 1:lag.MA.agg) {
        MA.fra.cov.agg<-rbind(MA.fra.cov.agg,cbind(matrix(0,k,i*k),diag(1,k),MA.cof.agg,matrix(0,k,(lag.MA.agg-i)*k)))
        for (j in (i):(dim(MA.fra.cov.agg)[2]/k-1-i)) {
          autocov.agg[,1:k+(k*(i-1))]<-autocov.agg[,1:k+(k*(i-1))]+MA.fra.cov.agg[1:k,(1:k)+j*k]%*%sig.agg%*%t(MA.fra.cov.agg[1:k+i*k,(1:k)+j*k])
        }
      }
      y[,(k+1):(k+k*lag.MA.agg)]<-autocov.agg-autocov
      y
    }
    if (all(start.values=="auto")) {
      xstart<-c(l*sig,matrix(0+ifelse(f=="avg",.2,0),k,k*lag.MA.agg))
    }else{
      if (nrow(start.values)!=k) {
        stop(paste0("rows of start.values not ",k))
      }
      if (ncol(start.values)!=(k+k*lag.MA.agg)) {
        stop(paste0("columns of start.values not ",k," x ",k*lag.MA.agg))
      }
      xstart<-start.values
    }

    MA.sol<-nleqslv::nleqslv(xstart,MA.function,control = list(maxit=200))
    sig.agg<-matrix(MA.sol[["x"]][1:k^2],k,k)
    MA.agg<-matrix(MA.sol[["x"]][(k^2+1):(k^2*lag.MA.agg+k^2)],k,k*lag.MA.agg)
    #agg.coeff<-list(pi.agg=dot.pi,AR.agg=dot.AR,MA.agg=MA.agg,sig.agg=sig.agg,level.AR.agg=AR.agg)


  }else{
    sig.agg<-sig.fra
    MA.agg<-matrix(0,k,k)
    #agg.coeff<-list(pi.agg=dot.pi,AR.agg=dot.AR,MA.agg=MA.agg,sig.agg=sig.agg,level.AR.agg=AR.agg)
  }

  if (include =="SRconst") {
    SRconst.agg<-B_L%mult%SRconst
    SRconst.agg<-rowSums(SRconst.agg@coef)
    dot.AR<-cbind(SRconst.agg,dot.AR)
  }
  out<-list(
    output="aggregated",
    model.specification=list(
      k = k,
      p = lag.AR-1,
      q = lag.MA.agg,
      r = r,
      include = include
    ),
    alphas = dot.pi[,1:r],
    betas = ecvarma$betas,
    AR = dot.AR,
    MA = MA.agg,
    Sigma = sig.agg,
    aggregation = list(
      l = l,
      scheme = f,
      method = "Marcellino",
      VAR.rep = AR.agg,
      pi.agg = dot.pi,
      B_L = B_L
    )
  )
  class(out)<-"ECVARMA"
  attr(out,"determination")<-"aggregated"
  return(out)
}
