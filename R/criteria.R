#' Cointegration criterion of Athanasopoulos et al. (2016)
#' @description an information criterion to determine the cointegration rank of a data set
#' @param input a matrix of time series data (the bigger dimension is considered to be the sample size)
#' @param detrend whether the data should be detrended, TRUE is recommended
#' @return rank of cointegration as numeric
#' @export
coint.criterion<-function(input,detrend=TRUE){
  input<-as.matrix(input)
  if (dim(input)[1]<dim(input)[2]) {
    input<-t(input)
  }
  t<-dim(input)[1]
  k<-dim(input)[2]
  Pt<-log(t)
  if(detrend){
    x<-cbind(rep(1,t),1:t)
    y<-input-x%*%(solve(crossprod(x))%*%t(x)%*%input)
    yt<-y[-1,]
    y1<-y[-t,]
  }else{
    yt<-input[-1,]
    y1<-input[-t,]
  }

  lam<-sort(c(cancor(yt,y1)[["cor"]])^2,decreasing = FALSE)#sort(c(cancor(yt^2,y1^2)[["cor"]]),decreasing = TRUE)
  if (lam[k]<=(1-(log(t)/t)^.5)) {
    cat("cointegration rank = ",k)
    invisible(k)
  }else{
  Lam<-numeric(k) #1:k instead of 0:k-1
  zeta<-numeric(k)
  for (rho in 0:(k-1)) {#i-1 is rho, i in k:1
    Lam[rho+1]<-(k-rho)^-1*sum(lam[(rho+1):k])/(prod(lam[(rho+1):k])^(1/(k-rho)))
    zeta[rho+1]<-t*(k-rho)*log(Lam[rho+1])+rho*(2*k-rho+1)*Pt/2
    #Lam[i-k+1]<-1/i*sum(lam[1:i])/(prod(lam[1:i])^(1/i))
  }
  # zeta<-numeric(k)
  # for (i in 1:k) {
  #   zeta[i]<-t*(k-i+1)*log(Lam[i])+(i-1)*(2*k-i)*Pt/2
  # }
  coin.rank<-which.min(zeta)-1
  cat("cointegration rank = ",coin.rank)
  invisible(coin.rank)
  }
}

#' FMA lag criterion by Kascha and Trenkler (2015)
#' @description
#' An information criterion to specify the number of AR and MA lags of an (EC)VARMA model in final moving average specification
#' @param input a matrix of time series data (the bigger dimension is considered to be the sample size)
#' @param h an integer that equals the maximum lags tested for the AR and MA dynamic, "auto" uses the recommendation of Kascha and Trenkler
#' @param v a parameter for the calculation of the lag criterion, must be greater than 0, Kascha and Trenkler use a value of 0.5
#' @return a list with the recommended number of AR and MA lags, and a matrix with the information criteria values, AR lags = rows (first AR lag reserved for error correction, total number of AR lags is then reduced by one), MA lags = columns
#' @export
lag.select.fma<-function(input,h="auto",v=0.5){
  input<-as.matrix(input)
  if (dim(input)[1]<dim(input)[2]) {
    input<-t(input)
  }
  if (h=="auto") {
    h<-floor(log(dim(input)[1])^1.25)
  }
  pmax<-h
  qmax<-h

  #1.Long VAR
  #if(h<max(pmax,qmax))stop("h must be bigger than pmax and qmax")
  if(dim(input)[1]>dim(input)[2]){input<-t(input)}
  k<-dim(input)[1]
  T<-dim(input)[2]
  #input of dimension kxt
  #mean centering of inputs
  for(i in 1:k){
    input[i,]<-input[i,]-mean(input[i,])
  }

  Y<-input[,-(1:h)]				#t columns			#k x t-h
  Z<-t(embed(t(input),h+1))[-(1:k),]				#kh x t-h
  u<-Y-Y%*%t(Z)%*%MASS::ginv(tcrossprod(Z))%*%Z	#k x t-h
  #2.Information criterion
  s<-max(pmax,qmax)+h
  Sig<-list()
  DP<-matrix(NA,pmax+1,qmax+1)
  for (p in 0:pmax){
    Sig[[p+1]]<-list()
    for (q in 0:qmax){
      if (p==0 & q!=0) {
        Z2<-t(embed(t(u),q+1))[-(1:k),]		#kp+kq x t-h
        Y<-t(embed(t(input),p+1))[(1:k),-(1:(h-p+q))]
        u2<-Y-Y%*%t(Z2)%*%MASS::ginv(tcrossprod(Z2))%*%Z2
        Sig.m<-(1/(T-s))*tcrossprod(u2)
        #restriction
        R<-rbind(cbind(diag(1,k*(k*p)),matrix(0,k*(k*p),q)),cbind(matrix(0,k*k*q,k*(k*p)),diag(1,q)%x%as.vector(diag(1,k))))
        gam<-MASS::ginv(t(R)%*%(tcrossprod(Z2)%x%MASS::ginv(Sig.m))%*%R)%*%t(R)%*%(Z2%x%MASS::ginv(Sig.m))%*%as.vector(Y)
        qmat<-diag(gam[k^2*p+1],k)
        if(q>1){
          for (o in 1:(q-1)) {
            qmat<-cbind(qmat,diag(gam[k^2*p+1+o],k))
          }
        }
        u3<-Y-qmat%*%Z2
        Sig.m<-(1/(T-s))*tcrossprod(u3)
      }else if (p!=0 & q==0) {
        Z2<-t(embed(t(input),p+1))[-(1:k),-(1:(h-p+q))]		#kp+kq x t-h
        Y<-t(embed(t(input),p+1))[(1:k),-(1:(h-p+q))]
        u2<-Y-Y%*%t(Z2)%*%MASS::ginv(tcrossprod(Z2))%*%Z2
        Sig.m<-(1/(T-s))*tcrossprod(u2)
        #restriction
        R<-rbind(cbind(diag(1,k*(k*p)),matrix(0,k*(k*p),q)),cbind(matrix(0,k*k*q,k*(k*p)),diag(1,q)%x%as.vector(diag(1,k))))
        gam<-MASS::ginv(t(R)%*%(tcrossprod(Z2)%x%MASS::ginv(Sig.m))%*%R)%*%t(R)%*%(Z2%x%MASS::ginv(Sig.m))%*%as.vector(Y)
        bet<-matrix(gam[1:(k^2*p)],k,p*k)
        u3<-Y-bet%*%Z2
        Sig.m<-(1/(T-s))*tcrossprod(u3)
      }else if (p==0 & q==0) {
        Z2<-0		#kp+kq x t-h
        u3<-Y
        Sig.m<-(1/(T-s))*tcrossprod(u3)
      }else{
        Z2<-rbind(t(embed(t(input),p+1))[-(1:k),-(1:(h-p+q))],t(embed(t(u),q+1))[-(1:k),])		#kp+kq x t-h
        Y<-t(embed(t(input),p+1))[(1:k),-(1:(h-p+q))]
        u2<-Y-Y%*%t(Z2)%*%MASS::ginv(tcrossprod(Z2))%*%Z2
        Sig.m<-(1/(T-s))*tcrossprod(u2)
        #restriction
        R<-rbind(cbind(diag(1,k*(k*p)),matrix(0,k*(k*p),q)),cbind(matrix(0,k*k*q,k*(k*p)),diag(1,q)%x%as.vector(diag(1,k))))
        gam<-MASS::ginv(t(R)%*%(tcrossprod(Z2)%x%MASS::ginv(Sig.m))%*%R)%*%t(R)%*%(Z2%x%MASS::ginv(Sig.m))%*%as.vector(Y)
        bet<-matrix(gam[1:(k^2*p)],k,p*k)
        qmat<-diag(gam[k^2*p+1],k)
        if(q>1){
          for (o in 1:(q-1)) {
            qmat<-cbind(qmat,diag(gam[k^2*p+1+o],k))
          }
        }
        bet<-cbind(bet,qmat)
        u3<-Y-bet%*%Z2
        Sig.m<-(1/(T-s))*tcrossprod(u3)
      }

      Sig[[p+1]][[q+1]]<-list(Sig.m)
      DP[p+1,q+1]<-log(det(Sig.m))+(p*k+q)*(log(T-s))^(1+v)/(T-s)
    }
  }
  rownames(DP)<-paste0(0:p)
  colnames(DP)<-paste0(0:q)
  minDP<-as.vector(which(DP==min(DP),arr.ind = TRUE))-1
  results<-list(lag.recommendation=minDP,DP.values=DP)
  cat("rows = AR order, columns = MA order \n")
  print(DP)
  return(results)
}
