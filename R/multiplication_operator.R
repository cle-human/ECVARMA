setGeneric("%mult%", function(x, y) {
  standardGeneric("%mult%")
})

setMethod("%mult%",signature ( x="polynomial",y="polyMatrix"),function(x,y){
  return(x*y)
})

setMethod("%mult%",signature ( x="polyMatrix",y="polynomial"),function(x,y){
  return(x*y)
})

setMethod("%mult%",signature ( x="polynomial",y="polynomial"),function(x,y){
  return(polyMatrix(x*y))
})

setMethod("%mult%",signature ( x="polyMatrix",y="polyMatrix"),function(x,y){
  if (ncol(x)!=nrow(y)) stop("ncol(x)!=nrow(y)")
  z<-polyMatrix(0,nrow(x),ncol(y),0)
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(y)) {
      for (ii in 1:ncol(x)) {
        z[i,j]<-z[i,j]+x[i,ii]*y[ii,j]
      }
    }
  }
  return(z)
})

setMethod("%mult%",signature ( x="matrix",y="polyMatrix"),function(x,y){
  if (ncol(x)!=nrow(y)) stop("ncol(x)!=nrow(y)")
  z<-polyMatrix(0,nrow(x),ncol(y),0)
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(y)) {
      for (ii in 1:ncol(x)) {
        z[i,j]<-z[i,j]+x[i,ii]*y[ii,j]
      }
    }
  }
  return(z)
})
setMethod("%mult%",signature ( x="polyMatrix",y="matrix"),function(x,y){
  if (ncol(x)!=nrow(y)) stop("ncol(x)!=nrow(y)")
  z<-polyMatrix(0,nrow(x),ncol(y),0)
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(y)) {
      for (ii in 1:ncol(x)) {
        z[i,j]<-z[i,j]+x[i,ii]*y[ii,j]
      }
    }
  }
  return(z)
})

setMethod("%mult%",signature ( x="numeric",y="polyMatrix"),function(x,y){
  return(x*y)
})
setMethod("%mult%",signature ( x="polyMatrix",y="numeric"),function(x,y){
  return(x*y)
})

setMethod("%mult%",signature ( x="numeric",y="matrix"),function(x,y){
  return(polyMatrix(x*y,nrow(y),ncol(y),0))
})
setMethod("%mult%",signature ( x="matrix",y="numeric"),function(x,y){
  return(polyMatrix(x*y,nrow(x),ncol(x),0))
})

setMethod("%mult%",signature ( x="polynomial",y="matrix"),function(x,y){
  z<-y
  for (i in 1:nrow(y)) {
    for (j in 1:ncol(y)) {
      z[i,j]<-x*y[i,j]
    }
  }
  return(polyMatrix(z,nrow(y),ncol(y),length(x)-1))
})
setMethod("%mult%",signature ( x="matrix",y="polynomial"),function(x,y){
  z<-x
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      z[i,j]<-y*x[i,j]
    }
  }
  return(polyMatrix(z,nrow(x),ncol(x),length(y)-1))
})

setMethod("%mult%",signature ( x="polynomial",y="numeric"),function(x,y){
  return(polyMatrix(x*y))
})

setMethod("%mult%",signature ( x="numeric",y="polynomial"),function(x,y){
  return(polyMatrix(x*y))
})

setMethod("%mult%",signature ( x="numeric",y="numeric"),function(x,y){
  return(polyMatrix(x*y,1,1,0))
})

setMethod("%mult%",signature ( x="matrix",y="matrix"),function(x,y){
  return(polyMatrix(x%*%y,nrow(x),ncol(y),0))
})
