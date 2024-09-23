#initial value estimation
ECVARMA.ini<-function(input,p,r,q,h,k,T,R.f,R0.f,include){
  #1.Long VAR
  Y<-input-matrix(rowMeans(input[,-1]-input[,-T]),k,T)%*%diag(1:(T))	#t columns
  Y<-Y-matrix(rowMeans(Y),k,T)
  Z<-t(embed(t(Y),h+1))[-(1:k),]						#kh x t-h
  Y<-Y[,-(1:h)]                                 #k x t-h
  u<-Y-Y%*%t(Z)%*%solve(tcrossprod(Z))%*%Z		  #k x t-h
  Sig<-tcrossprod(u)/(T-h)
  #iSig<-solve(Sig)
  iSig<-diag(k)

  #2.First values

  Y <-input																												#k x T
  dY<-Y[,-1]-Y[,-dim(Y)[2]] 																			#k x T-1
  Z<-rbind(t(embed(t(dY[,-(1:(h-p+q-1))]),p+1))[-(1:k),],t(embed(t(u),q+1))[-(1:k),]) 	#kp+kq x T-h-q-1
  Y <-Y[,-c(1:(h+q-1),dim(Y)[2])]																	      	#k x T-h-q-1
  dY<-dY[,-(1:(h+q-1))]																				      			#k x T-h-q-1
  #Z<-rbind(t(embed(t(dY[,-(1:(h))]),p+1))[-(1:k),],t(embed(t(u),q+1))[-(1:k),]) 	#kp+kq x T-h-p-1
  if (include=="SRconst") {
    Z<-rbind(1,Z)
  }
  #Y <-Y[,-c(1:(h+p-1),dim(Y)[2])]																	      	              #k x T-h-p-1
  if (include=="LRconst") {
    Y<-rbind(Y,1)
  }

  X<-array(0,dim = c(T-h+q-1,k,dim(R.f)[2]))   #k*(k+k*p+k*q)   +k if const
  const<-array(0,dim = c(T-h+q-1,k,1)) #k x 1
  i<-1

  for(i in 1:(T-h-q-1)){
    X[i,,]<-cbind(t(dY[,i])%x%diag(1,k),t(Y[,i])%x%diag(1,k),t(Z[,i])%x%diag(1,k))%*%R.f   #k x number of free parameter
    const[i,,]<-cbind(t(dY[,i])%x%diag(1,k),t(Y[,i])%x%diag(1,k),t(Z[,i])%x%diag(1,k))%*%R0.f  #k x 1
  }

  gam<-numeric(length=dim(R.f)[2]) #k^2*(p+p+1)
  XXsum<-matrix(0,dim(R.f)[2],dim(R.f)[2])
  Xysum<-matrix(0,dim(R.f)[2],1)
  for (i in 1:(T-h-q-1)){
    XXsum<-XXsum+t(X[i,,])%*%iSig%*%X[i,,]
    Xysum<-Xysum+t(X[i,,])%*%iSig%*%(dY[,i]-const[i,,])
  }

  gam<-solve(XXsum)%*%Xysum  #Pi sind die ersten k^2(+k) Elemente
  para<-R.f%*%gam+R0.f

  A0<-matrix(para[1:(k^2)],k,k)#-diag(1,k)
  pi.mat<-matrix(para[1:(k^2+ifelse(include=="LRconst",k,0))+k^2],k,k+ifelse(include=="LRconst",1,0))
  alp<-pi.mat[,1:r]
  M<-matrix(para[(length(para)-k^2*q+1):(length(para))],k,k*q)
  G<-matrix(para[(1):(ifelse(include=="SRconst",k,0)+(p)*k^2)+2*k^2+ifelse(include=="LRconst",k,0)],k,k*(p)+ifelse(include=="SRconst",1,0))
  if (all(pi.mat==0)) {
    bet<-matrix(-1,r,k-r)
  }else if(all(M==0)){
    bet<-MASS::ginv(t(alp)%*%MASS::ginv(diag(1,q)%x%Sig)%*%alp)%*%(t(alp)%*%MASS::ginv(diag(1,q)%x%Sig)%*%pi.mat[,(r+1):(k+ifelse(include=="LRconst",1,0))])
  }else{
    bet<-MASS::ginv(t(alp)%*%MASS::ginv(M%*%(diag(1,q)%x%Sig)%*%t(M))%*%alp)%*%(t(alp)%*%MASS::ginv(M%*%(diag(1,q)%x%Sig)%*%t(M))%*%pi.mat[,(r+1):(k+ifelse(include=="LRconst",1,0))])
  }

  del<-c(as.vector(A0),as.vector((bet)),as.vector(alp),as.vector(G),as.vector(M))
  results<-list()
  results$del<-del
  results$sig<-Sig
  return(results)
}
