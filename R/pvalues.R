#p-values
ECVARMA.pvalues<-function(input,T,k,d,del.r,M,bet,p,r,q,R,include,Sig,u){
  #Error cov
  Y<-input																						#k x T
  dY<-Y[,-1]-Y[,-dim(Y)[2]]		                      	#k x T-1
  Y<-Y[,-dim(Y)[2]]																		#k x T-1
  if(include=="LRconst"){
    Y<-rbind(Y,1)                                     #k+1 x T-1
  }

  iSig<-solve(Sig)

  U<-matrix(0,T-1,r+p*k+q*k+ifelse(include=="SRconst",1,0)) #U trans
  for (i in 1:(T-1)) {
    U[i,1:r]<-t(cbind(diag(1,r),bet)%*%as.matrix(Y[,i]))  #t((r x k) x (kx1)) => 1xr ###const
    if (include=="SRconst") {
      U[i,r+1]<-1
    }
    for (g in 1:p) {
      if(i-g<1){
        U[i,(r+1+k*(g-1)):(r+g*k)+ifelse(include=="SRconst",1,0)]<-0
      }else{
        U[i,(r+1+k*(g-1)):(r+g*k)+ifelse(include=="SRconst",1,0)]<-as.matrix(dY[,i-g])
      }
    }

    for (g in 1:q) {
      if(i-g<1){
        U[i,((r+1+k*(g-1)):(r+g*k))+p*k+ifelse(include=="SRconst",1,0)]<-0
      }else{
        U[i,((r+1+k*(g-1)):(r+g*k))+p*k+ifelse(include=="SRconst",1,0)]<- as.matrix(u[,i-g]) ####
      }
    }
  }

  R.minbet<-R[-(1:(r*d+ifelse(include=="LRconst",r,0))+k^2),]
  for(g in dim(R.minbet)[2]:1){
    if(all(R.minbet[,g]==0)){
      R.minbet<-R.minbet[,-g]
    }
  }

  cfv<-function(R){
    zero.values<-apply(R,2,function(col)any(col!=0))
    sum(zero.values)
  }
  free.A0<-cfv(R[1:k^2,])
  free.beta<-cfv(R[1:(r*d+ifelse(include=="LRconst",r,0))+k^2,])
  if (free.beta>0) {
    del.r<-del.r[-(1:free.beta+free.A0)]
  }

  #del.r<-MASS::ginv(R.minbet)%*%del

  Us<-array(0,dim = c(T-1,k,dim(R.minbet)[2]))
  for (i in 1:(T-1)) {
    Us[i,,]<- -cbind(t(as.matrix(dY[,i]))%x%diag(-1,k),t(U[i,])%x%diag(k))%*%R.minbet
    for (g in 1:q) {
      if(i-g<1){
        Us[i,,]<-Us[i,,]
      }else{
        Us[i,,]<-Us[i,,]-M[,(1:k)+(g-1)*k]%*%Us[i-g,,]#A0i%*%
      }
    }
  }

  UsUssum<-matrix(0,dim(R.minbet)[2],dim(R.minbet)[2])

  for(j in (1):(T-1)){
    UsUssum<-UsUssum+t(Us[j,,])%*%solve(Sig)%*%Us[j,,]
  }

  V<-UsUssum #kein /T ? #p 479
  iV<-solve(V)
  std.error<-numeric(dim(R.minbet)[2])
  p.values<-numeric(dim(R.minbet)[2])
  for (i in 1:dim(R.minbet)[2]){
    std.error[i]<-iV[i,i]^0.5
    p.values[i]<-pt( abs(del.r[i]/std.error[i]),T-dim(R)[2],lower.tail = FALSE)*2
  }

  results<-list()
  results$del.r<-del.r
  results$std.error<-std.error #only for non-fixed parameters
  results$p.values<-p.values
  return(results)
}
