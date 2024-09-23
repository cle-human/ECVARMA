# #test ECVARMA estimations
# #library(mvtnorm)
# #sample size
# x<-480*10#*10
# #aggregation period length
# y<-4
# #number of samples
# s<-20#10000
#
# alphas<-matrix(c(-0.08,0.09),2,1)
# betas<-matrix(c(1,-1.03),2,1)
# pi.mat<-alphas%*%t(betas)
# AR<-matrix(c(.72,.29,.06,.29,-.28,0,-.01,.16),2,4)
# sig.non.agg=matrix(c(5.21,0.85,0.85,4.83),2,2)
#
# # test.ecvarma<-ecvarma.define(A0=NULL,alphas,betas,AR,MA=NULL,Sigma=sig.non.agg,include = "none")
# # test.agg<-ecvarma.agg.M(test.ecvarma, 4, "bop")
# # test.agg.fma<-ecvarma.toFMA(test.agg)
# #
# # test.fma.sim<-ecvarma.sim(test.agg.fma,10000,100,seed = 10)
# # test.fma.est<-ecvarma.fma(as.matrix(test.fma.sim),1,1,2,include = "none")
#
# # test.data<-ecvarma.sim(test.ecvarma,x,20)
# #set.seed(s+100)
# e<-t(mvtnorm::rmvnorm(x+100,mean = rep(0,2),sigma = sig.non.agg))
# dp<-matrix(0,2,x+100)
# p<-matrix(0,2,x+100)
#
# for(i in 1:(x+100-2)){
#   dp[,i+2]<-pi.mat%*%p[,i+1]+t(t(e[,i+2]))+AR[,1:2]%*%dp[,i+1]+AR[,1:2+2]%*%dp[,i]
#   p[,i+2]<-p[,i+1]+dp[,i+2]
# }
# p<-p[,-(1:100)]
# dp<-dp[,-(1:100)]
#
# input<-p
# R<-diag(4+2+2+8+4)
# R<-R[,-c(1:4,17:20)]
#
# test.est<-ecvarma(input = input,p=2,r=1,q=1,n=5,h="auto",include = "LRconst",R=R,R0=NULL,initial_values = "auto",p.values = FALSE)
#
#
# # alphas<-matrix(round(rnorm(k,0,0.1),2),k,1)
# # alphas<-matrix(c(-.1,.1,.1,rep(0,l)),k,1)
# # #alphas<-matrix(c(-.1,.1,.1,.1,rep(0,1)),k,1)
# # alphas<-matrix(c(-.1,.1,.1,0,0,.1,-.1,.1,.1,0),k,2)
# # # if (alphas[1,]>0)alphas[1,]<-alphas[1,]*-1
# # # for (i in 2:k) {
# # #   if (alphas[i,]<0)alphas[i,]<-alphas[i,]*-1
# # # }
# # betas<-matrix(c(1,rep(-1,k-1)),k,1)
# # betas<-matrix(c(1,rep(-1,m-1),rep(0,l)),k,1)
# # betas<-matrix(c(1,rep(-1,m-2),rep(0,l+1),0,1,0,-1,0),k,2)
# # AR<-matrix(round(rnorm(k^2,0,0.1),2),k,k)
# # LRconst<-matrix(c(1,0),2,1)
# # SRconst<-matrix(rep(1,5),5,1)
# # VECMagg.con(alphas,betas,AR=AR,SRconst=SRconst,LRconst=NULL,MA=NULL,l=2,Fmat=NULL,sig=diag(k))
# # #AR<-matrix(0,k,k)
# # Fmat<-matrix(c(1/3,1/3,1/3,0,0,
# #                0,0,0,1/2,1/2,
# #                0,1,0,0,0,
# #                0,0,1,0,0,
# #                0,0,0,0,1),5,5,byrow = T)
# #VECMagg.con(alphas,betas,AR=AR,SRconst=SRconst,LRconst=NULL,MA=NULL,l=2,Fmat=Fmat,sig=diag(k))
#
# # solve(matrix(c(.5,.5,0,0,0,0,
# #                0,0,.25,.25,.25,.25,
# #                0,1,0,0,0,0,
# #                0,0,0,1,0,0,
# #                0,0,0,0,1,0,
# #                0,0,0,0,0,1
# #                ),6,6,byrow = T))
#
# A0<-matrix(c(1,0,0,
#              .2,1,0,
#              .1,-.05,1),3,3,byrow = T)
# alphas<-matrix(c(-0.08,0.09,.1,.18,-.09,.17),3,2)
# betas<-matrix(c(1,0,0,1,-1.03,-1.05,-5,-10),2,4)
# #pi.mat<-alphas%*%t(betas)
# AR<-matrix(c(-0.21, 0, 0, 0.23, 0, 0, -0.3, 0, 0),3,3)
# MA<-matrix(c(-0.04, 0.17, 0, -0.3, -0.09, 0, 0, 0.4, 0),3,3)
# sig.non.agg=diag(3)
#
# test.process<-ecvarma.define(A0=A0,alphas,betas,AR,MA,Sigma=sig.non.agg,include = "LRconst")
# test.data<-ecvarma.sim(test.process,500,100,seed = 201010234)
#
# MTS::SCMid(test.data)
#
# sc=matrix(c(2,1,
#             1,1,
#             1,0),3,2,byrow = T)
#
# test.est<-ecvarma.scm(test.data,sc,r=2,include = "LRconst",Newton = F,n=20)
#
# A0<-diag(3)
# alphas<-matrix(c(-0.08,0.09,.1,.18,-.09,.17),3,2)
# betas<-matrix(c(1,0,0,1,-1.03,-1.05,-5,-10),2,4)
# #pi.mat<-alphas%*%t(betas)
# AR<-matrix(c(-0.16, 0.02, 0.07, -0.01, -0.41, -0.15, 0.25, 0.06, -0.15),3,3)
# MA<-matrix(c(-0.1, 0, 0, 0, -0.1, 0, 0, 0, 0.1),3,3)
#
# test.process<-ecvarma.define(A0=A0,alphas,betas,AR,MA,Sigma=sig.non.agg,include = "LRconst")
# test.data<-ecvarma.sim(test.process,500,100,seed = 201331)
#
# test.est<-ecvarma.fma(input=test.data,r=2,p=1,q=1,include="LRconst")
