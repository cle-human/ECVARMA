#restriction function
ecvarma.res<-function(R,r,include,sc,p,q,k,d){
  # 0<r<k
  if (R=="FMA") {
    R<-diag(k^2+r*d+ifelse(include=="LRconst",r,0)+k*r+ifelse(include=="SRconst",k,0)+max(p,1)*k^2)
    Rcol<-dim(R)[2]
    if (p==0) {
      R<-R[,-((Rcol-k^2+1):Rcol)]
    }
    R<-R[,-(1:k^2)]
    for (i in 1:(k^2*max(1,q))) {
      R<-rbind(R,0)
    }
    if (q>0) {
      R<-cbind(R,matrix(0,dim(R)[1],q))
      R[(dim(R)[1]-k^2*q+1):dim(R)[1],(dim(R)[2]-q+1):dim(R)[2]]<-diag(1,q)%x%as.vector(diag(1,k))
    }
    R.f<-R[-(1:(r*d+ifelse(include=="LRconst",r,0)+k*r+k^2)),-(1:(r*d+ifelse(include=="LRconst",r,0)+k*r))]
    R.f<-cbind(rbind(matrix(0,k^2,k^2+ifelse(include=="LRconst",k,0)),diag(k^2+ifelse(include=="LRconst",k,0)),matrix(0,ifelse(include=="SRconst",k,0)+max(p,1)*k^2+max(q,1)*k^2,k^2+ifelse(include=="LRconst",k,0))),rbind(matrix(0,k^2+k^2+ifelse(include=="LRconst",k,0),dim(R.f)[2]),R.f))

    R0<-matrix(0,dim(R)[1],1)
    R0.f<-matrix(0,dim(R.f)[1],1)
    #return(R)
  } else if (R=="DMA") {
    R<-diag(k^2+r*d+ifelse(include=="LRconst",r,0)+k*r+ifelse(include=="SRconst",k,0)+max(p,1)*k^2)
    Rcol<-dim(R)[2]
    if (p==0) {
      R<-R[,-((Rcol-k^2+1):Rcol)]
    }
    R<-R[,-(1:k^2)]
    for (i in 1:(k^2*max(1,q))) {
      R<-rbind(R,0)
    }
    if (q>0) {
      R<-cbind(R,matrix(0,dim(R)[1],q*k))
      dia.mat<-matrix(0,k^2,k)#diag MA
      n<-1
      for (i in 1:k) {
        dia.mat[n,i]<-1
        n<-n+k+1
      }
      R[(dim(R)[1]-k^2*q+1):dim(R)[1],(dim(R)[2]-q*k+1):dim(R)[2]]<-diag(1,q)%x%dia.mat
    }
    R.f<-R[-(1:(r*d+ifelse(include=="LRconst",r,0)+k*r+k^2)),-(1:(r*d+ifelse(include=="LRconst",r,0)+k*r))]
    R.f<-cbind(rbind(matrix(0,k^2,k^2+ifelse(include=="LRconst",k,0)),diag(k^2+ifelse(include=="LRconst",k,0)),matrix(0,ifelse(include=="SRconst",k,0)+max(p,1)*k^2+max(q,1)*k^2,k^2+ifelse(include=="LRconst",k,0))),rbind(matrix(0,k^2+k^2+ifelse(include=="LRconst",k,0),dim(R.f)[2]),R.f))

    R0<-matrix(0,dim(R)[1],1)
    R0.f<-matrix(0,dim(R.f)[1],1)
    #return(R)
  }else if (R=="SCM") { #sc 1.col p, 2.col q
    p<-max(max(sc[,1])-1,1)
    p.i<-pmax(sc[,1]-1,0)
    ec.i<-pmin(sc[,1],1)
    q<-max(max(sc[,2]),1)
    q.i<-sc[,2]

    A<-matrix(1,k,k)-diag(k)
    alp<-matrix(1,k,r)
    bet<-matrix(1,r,d+ifelse(include=="LRconst",1,0))
    G<-matrix(1,k,k*p+ifelse(include=="SRconst",1,0))
    M<-matrix(1,k,k*q)

    for (i in 1:k) {
      if (ec.i[i]==0) {
        alp[i,]<-0
      }
      if (p.i[i]<p) {
        G[i,(k*p):(p.i[i]*k+1)+ifelse(include=="SRconst",1,0)]<-0
      }
      if (q.i[i]<q) {
        M[i,(k*q):(sc[i,2]*k+1)]<-0
      }
    }
    #add 1. res and A
    for (i in 1:k) {
      for (j in 1:k) {
        if ((sc[i,1]>sc[j,1])&(sc[i,2]>sc[j,2])) {
          M[i,j+((1:min(sc[i,1]-sc[j,1],sc[i,2]-sc[j,2])-1)*k)]<-0
        }
        if ((sc[i,1]>=sc[j,1])&(sc[i,2]>=sc[j,2])) {
          A[i,j]<-0
        }
      }
    }
    R<-diag(c(A,bet,alp,G,M))
    for (i in dim(R)[2]:1) {
      if (all(R[,i]==0)) {
        R<-R[,-i]
      }
    }
    R.f<-diag(c(A,sign(alp%*%cbind(diag(r),bet)),G,M))
    for (i in dim(R.f)[2]:1) {
      if (all(R.f[,i]==0)) {
        R.f<-R.f[,-i]
      }
    }
    R0<-matrix(0,dim(R)[1],1)
    R0.f<-matrix(0,dim(R.f)[1],1)
  }

  result<-list(R=R,R0=R0,R.f=R.f,R0.f=R0.f)
  return(result)
}
