ECVARMA.iter<-function(input,del,sig0,k,T,p,r,q,n=5,R,R0,include,Newton){
  d<-k-r
  Sig<-sig0
  del0<-del
  #3.Iteration
  del.r<-MASS::ginv(R)%*%del

  Y<-input																						#k x T
  dY<-Y[,-1]-Y[,-dim(Y)[2]]		                      	#k x T-1
  Y<-Y[,-dim(Y)[2]]																		#k x T-1
  if(include=="LRconst"){
    Y<-rbind(Y,1)                                     #k+1 x T-1
  }

  results<-list()

  #start iteration

  fit<-function(del.r){
    del<-R%*%del.r+R0
    #del to coefficient matrices
    A0<-diag(k)-matrix(del[1:k^2],k,k)
    #A0i<-solve(A0)
    G<-matrix(del[1:(ifelse(include=="SRconst",k,0)+(p)*k^2)+k^2+(k*r)+r*(k-r)+ifelse(include=="LRconst",r,0)],k,k*(p)+ifelse(include=="SRconst",1,0))
    M<-matrix(del[(length(del)-k^2*q+1):(length(del))],k,k*q)
    alp<-matrix(del[1:(k*r)+r*(k-r)+ifelse(include=="LRconst",r,0)+k^2],k,r)
    bet<-matrix(del[1:(r*(k-r)+ifelse(include=="LRconst",r,0))+k^2],r,k-r+ifelse(include=="LRconst",1,0))

    #error
    u<-matrix(0,k,T-1)

    for (i in 1:(T-1)) {
      u[,i]<-A0%*%dY[,i]-alp%*%(cbind(diag(1,r),(bet)))%*%as.matrix(Y[,i]) #t(bet)
      if (include=="SRconst") {
        u[,i]<-u[,i]-G[,1]
      }
      for (g in 1:(p)) {
        if(i-g<1){
          u[,i]<-u[,i]
        }else{
          u[,i]<-u[,i]-G[,(1:k)+k*(g-1)+ifelse(include=="SRconst",1,0)]%*%as.matrix(dY[,i-g])##const
        }
      }

      for (g in 1:q) {
        if(i-g<1){
          u[,i]<-u[,i]
        }else{
          u[,i]<-u[,i]-M[,(1:k)+k*(g-1)]%*%as.matrix(u[,i-g])###
        }
      }
    }
    U<-matrix(0,T-1,r+(p)*k+q*k+ifelse(include=="SRconst",1,0)) #U trans
    for (i in 1:(T-1)) {
      U[i,1:r]<-t(cbind(diag(1,r),bet)%*%as.matrix(Y[,i]))  #t((r x k) x (kx1)) => 1xr ###const
      if (include=="SRconst") {
        U[i,r+1]<-1
      }
      for (g in 1:(p)) {
        if(i-g<1){
          U[i,(r+1+k*(g-1)):(r+g*k)+ifelse(include=="SRconst",1,0)]<-0
        }else{
          U[i,(r+1+k*(g-1)):(r+g*k)+ifelse(include=="SRconst",1,0)]<-as.matrix(dY[,i-g])
        }
      }

      for (g in 1:q) {
        if(i-g<1){
          U[i,((r+1+k*(g-1)):(r+g*k))+(p)*k+ifelse(include=="SRconst",1,0)]<-0
        }else{
          U[i,((r+1+k*(g-1)):(r+g*k))+(p)*k+ifelse(include=="SRconst",1,0)]<- as.matrix(u[,i-g]) ####
        }
      }
    }

    #Error cov
    Sig<-1/(T-1)*tcrossprod(u)
    if (Newton==FALSE) {
      return(det(Sig))
    }else{
      tryCatch(
        expr ={
          Xs<-array(0,dim = c(T-1,k,dim(R)[2]))
          for (i in 1:(T-1)) {
            Xs[i,,]<- -cbind(t(as.matrix(dY[,i]))%x%diag(1,k),t(cbind(matrix(0,d,r),diag(d))%*%as.matrix(Y[-(k+1),i]))%x%(alp),if(include=="LRconst"){t(diag(d)%*%as.matrix(Y[(k+1),i]))%x%(alp)},t(as.matrix(U[i,]))%x%diag(k))%*%R

            #-(t(as.matrix(dY[,i]))%x%diag(1,k))%*%R[(1:k^2),]+cbind(t(cbind(matrix(0,d,r),diag(d))%*%as.matrix(Y[-(k+1),i]))%x%(alp),if(include=="LRconst"){t(diag(d)%*%as.matrix(Y[(k+1),i]))%x%(alp)},t(U[i,])%x%diag(k))%*%R[-(1:k^2),] #%x%A0?

            #-(t(A0i%*%(cbind(alp,G,M)%*%U[i,]))%x%A0i)%*%R[(1:k^2),]+cbind(t(cbind(matrix(0,d,r),diag(1,d))%*%as.matrix(Y[-(k+1),i]))%x%(A0i%*%alp),if(include=="LRconst"){t(diag(1,d)%*%as.matrix(Y[(k+1),i]))%x%(A0i%*%alp)},t(U[i,])%x%A0i)%*%R[-(1:k^2),]

            for (g in 1:q) {
              if(i-g<1){
                Xs[i,,]<-Xs[i,,]
              }else{
                Xs[i,,]<-Xs[i,,]-M[,(1:k)+(g-1)*k]%*%Xs[i-g,,]#A0i%*%
              }
            }
          }

          XsXssum<-matrix(0,dim(R)[2],dim(R)[2])#r*k+(p+p)*k^2
          Xsusum<-matrix(0,dim(R)[2],1)
          for(j in (1):(T-1)){
            XsXssum<-XsXssum+t(Xs[j,,])%*%solve(Sig)%*%Xs[j,,]
            Xsusum<-Xsusum+t(as.vector(u[,j])%*%solve(Sig)%*%Xs[j,,])
          }

          del.r<-del.r-solve(XsXssum)%*%Xsusum
          del<-R%*%del.r+R0
          ##############

          #collection of values of each iteration
          if (any(abs(del[-(1:(k^2+d*r+r*ifelse(include=="LRconst",1,0))),])>10)) {
            return(list(message="Iteration failed",del=del0,Sig=sig0,det.Sig=det(sig0)))
          }
          results<-list(message="Iteration successful",del.r=del.r,Sig=Sig,det.Sig=det(Sig))
          return(results)
        },
        error=function(e){
          return(list(message="Iteration failed",del=del,Sig=Sig,det.Sig=det(Sig)))
        }
      )
    }
  }
  results$initial.values<-del.r
  results$initial.sig<-Sig

  if(Newton){
    Sig.0<-abs(det(Sig))
    nn<-1
    #results$iterations<-list()
    while(n>=nn){
      sol<-fit(del.r)
      results[["iterations"]][[nn]]<-sol
      if (results[["iterations"]][[nn]][["message"]]=="Iteration failed") {
        break
      }

      #Sig<-results[["iterations"]][[nn]][["message"]]
      del.r<-sol$del.r
      if (Sig.0-abs(results[["iterations"]][[nn]][["det.Sig"]])<0.000001&nn>=2) {
        break
      }else{
        Sig.0<-abs(results[["iterations"]][[nn]][["det.Sig"]])
      }
      nn<-nn+1
    }
  }else{
    sol<-stats::nlminb(start=del.r,objective=fit,control = list(iter.max=n))
    results$sol<-sol
    del.r<-sol$par
  }

  #error terms and covariance matrix
  del<-R%*%del.r+R0
  #del to coefficient matrices
  A0<-diag(k)-matrix(del[1:k^2],k,k)
  #A0i<-solve(A0)
  G<-matrix(del[1:(ifelse(include=="SRconst",k,0)+(p)*k^2)+k^2+(k*r)+r*(k-r)+ifelse(include=="LRconst",r,0)],k,k*(p)+ifelse(include=="SRconst",1,0))
  M<-matrix(del[(length(del)-k^2*q+1):(length(del))],k,k*q)
  alp<-matrix(del[1:(k*r)+r*(k-r)+ifelse(include=="LRconst",r,0)+k^2],k,r)
  bet<-matrix(del[1:(r*(k-r)+ifelse(include=="LRconst",r,0))+k^2],r,k-r+ifelse(include=="LRconst",1,0))

  u<-matrix(0,k,T-1)

  for (i in 1:(T-1)) {
    u[,i]<-A0%*%dY[,i]-alp%*%(cbind(diag(1,r),(bet)))%*%as.matrix(Y[,i]) #t(bet)
    if (include=="SRconst") {
      u[,i]<-u[,i]-G[,1]
    }
    for (g in 1:(p)) {
      if(i-g<1){
        u[,i]<-u[,i]
      }else{
        u[,i]<-u[,i]-G[,(1:k)+k*(g-1)+ifelse(include=="SRconst",1,0)]%*%as.matrix(dY[,i-g])##const
      }
    }

    for (g in 1:q) {
      if(i-g<1){
        u[,i]<-u[,i]
      }else{
        u[,i]<-u[,i]-M[,(1:k)+k*(g-1)]%*%as.matrix(u[,i-g])###
      }
    }
  }

  #Error cov
  Sig<-1/(T-1)*tcrossprod(u)

  #results$del<-del
  results$del.r<-del.r
  results$A0<-A0
  results$bet<-bet
  results$alp<-alp
  results$G<-G
  results$M<-M
  results$Sig<-Sig
  results$errors<-u
  return(results)
}
