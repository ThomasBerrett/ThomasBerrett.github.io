############ Proposed Methods #############

library(mvtnorm)
library(latex2exp)
library(copula)
library(mice)

# Direct estimator
thetaDirect=function(X,Y,X1,X2,M,h){
  lambda1=length(X1)/length(Y); lambda2=length(X2)/length(Y)
  
  alpha1=0; alpha2=0
  if(M>1){
    for(m in 1:(M-1)){
      alpha1new=(lambda1/(1+lambda1))*ksmooth(X[,1],Y-alpha2,x.points=X[,1],bandwidth = h, kernel="normal")$y[rank(X[,1])]
      alpha2=(lambda2/(1+lambda2))*ksmooth(X[,2],Y-alpha1,x.points=X[,2],bandwidth = h, kernel="normal")$y[rank(X[,2])]
      alpha1=alpha1new
    }
  }
  
  alpha1comp=(lambda1/(1+lambda1))*ksmooth(X[,1],Y-alpha2,x.points=X[,1],bandwidth = h, kernel="normal")$y[rank(X[,1])]
  alpha1miss=(lambda1/(1+lambda1))*ksmooth(X[,1],Y-alpha2,x.points=X1,bandwidth = h, kernel="normal")$y[rank(X1)]
  alpha2comp=(lambda2/(1+lambda2))*ksmooth(X[,2],Y-alpha1,x.points=X[,2],bandwidth = h, kernel="normal")$y[rank(X[,2])]
  alpha2miss=(lambda2/(1+lambda2))*ksmooth(X[,2],Y-alpha1,x.points=X2,bandwidth = h, kernel="normal")$y[rank(X2)]
  
  return(mean(Y)-mean(alpha1comp)-mean(alpha2comp)+mean(alpha1miss)+mean(alpha2miss))
}


# Cross-fit estimator (half of)
thetaf=function(X,Y,X1,X2,Xt,Yt,M,h){
  lambda1=length(X1)/(length(Y)+length(Yt)); lambda2=length(X2)/(length(Y)+length(Yt))
  n=length(Y)
  
  alpha1=0; alpha2=0
  if(M>1){
    for(m in 1:(M-1)){
      alpha1new=(lambda1/(1+lambda1))*ksmooth(X[((m-1)*n/M+1):(m*n/M),1],Y[((m-1)*n/M+1):(m*n/M)]-alpha2,
                                              x.points=X[(m*n/M+1):((m+1)*n/M),1],bandwidth = h, kernel="normal")$y[rank(X[(m*n/M+1):((m+1)*n/M),1])]
      alpha2=(lambda2/(1+lambda2))*ksmooth(X[((m-1)*n/M+1):(m*n/M),2],Y[((m-1)*n/M+1):(m*n/M)]-alpha1,
                                           x.points=X[(m*n/M+1):((m+1)*n/M),2],bandwidth = h, kernel="normal")$y[rank(X[(m*n/M+1):((m+1)*n/M),2])]
      alpha1=alpha1new
    }
  }
  
  alpha1comp=(lambda1/(1+lambda1))*ksmooth(X[((M-1)*n/M+1):n,1],Y[((M-1)*n/M+1):n]-alpha2,
                                           x.points=Xt[,1],bandwidth = h, kernel="normal")$y[rank(Xt[,1])]
  alpha1miss=(lambda1/(1+lambda1))*ksmooth(X[((M-1)*n/M+1):n,1],Y[((M-1)*n/M+1):n]-alpha2,
                                           x.points=X1,bandwidth = h, kernel="normal")$y[rank(X1)]
  alpha2comp=(lambda2/(1+lambda2))*ksmooth(X[((M-1)*n/M+1):n,2],Y[((M-1)*n/M+1):n]-alpha1,
                                           x.points=Xt[,2],bandwidth = h, kernel="normal")$y[rank(Xt[,2])]
  alpha2miss=(lambda2/(1+lambda2))*ksmooth(X[((M-1)*n/M+1):n,2],Y[((M-1)*n/M+1):n]-alpha1,
                                           x.points=X2,bandwidth = h, kernel="normal")$y[rank(X2)]
  
  return(mean(Y)-mean(alpha1comp,na.rm=TRUE)-mean(alpha2comp,na.rm=TRUE)+mean(alpha1miss,na.rm=TRUE)+mean(alpha2miss,na.rm=TRUE))
}

## Cross-fit Estimator
thetac=function(X,Y,X1,X2,M,h){
  n=length(Y)
  theta1=thetaf(X[1:(n/2),],Y[1:(n/2)],X1,X2,X[(n/2+1):n,],Y[(n/2+1):n],M,h)
  theta2=thetaf(X[(n/2+1):n,],Y[(n/2+1):n],X1,X2,X[1:(n/2),],Y[1:(n/2)],M,h)
  return(0.5*theta1+0.5*theta2)
}

########### Setting (i) ############

n=600; lambda1=10; lambda2=10
sigY=0.3; Rho=0.1*(0:10)
h=0.1
Nrep=10000

MSE0=rep(0,length(Rho)); Var0=rep(0,length(Rho))
MSE1=rep(0,length(Rho)); Var1=rep(0,length(Rho))
MSE2=rep(0,length(Rho)); Var2=rep(0,length(Rho))
MSE3=rep(0,length(Rho)); Var3=rep(0,length(Rho))
MSE4=rep(0,length(Rho)); Var4=rep(0,length(Rho))
MSE1c=rep(0,length(Rho)); Var1c=rep(0,length(Rho))
MSE2c=rep(0,length(Rho)); Var2c=rep(0,length(Rho))
MSE3c=rep(0,length(Rho)); Var3c=rep(0,length(Rho))
MSE4c=rep(0,length(Rho)); Var4c=rep(0,length(Rho))
for(i in 1:length(Rho)){
  print(i)
  Theta0=rep(0,Nrep)
  Theta1=rep(0,Nrep)
  Theta2=rep(0,Nrep)
  Theta3=rep(0,Nrep)
  Theta4=rep(0,Nrep)
  Theta1c=rep(0,Nrep)
  Theta2c=rep(0,Nrep)
  Theta3c=rep(0,Nrep)
  Theta4c=rep(0,Nrep)
  for(nrep in 1:Nrep){
    X=pnorm(rmvnorm(n,sigma=matrix(c(1,-Rho[i],-Rho[i],1),2,2)))
    Y=X[,1]-X[,2]+rnorm(n,sd=sigY)
    
    X1=runif(n*lambda1); X2=runif(n*lambda2)
    
    Theta0[nrep]=mean(Y) 
    Theta1[nrep]=thetaDirect(X,Y,X1,X2,1,h)
    Theta2[nrep]=thetaDirect(X,Y,X1,X2,2,h)
    Theta3[nrep]=thetaDirect(X,Y,X1,X2,3,h)
    Theta4[nrep]=thetaDirect(X,Y,X1,X2,4,h)
    Theta1c[nrep]=thetac(X,Y,X1,X2,1,h)
    Theta2c[nrep]=thetac(X,Y,X1,X2,2,h)
    Theta3c[nrep]=thetac(X,Y,X1,X2,3,h)
    Theta4c[nrep]=thetac(X,Y,X1,X2,4,h)
    
    thetaTRUE=0
    
  }
  MSE0[i]=mean((Theta0-thetaTRUE)^2); MSE1[i]=mean((Theta1-thetaTRUE)^2); MSE2[i]=mean((Theta2-thetaTRUE)^2); MSE3[i]=mean((Theta3-thetaTRUE)^2); MSE4[i]=mean((Theta4-thetaTRUE)^2)
  MSE1c[i]=mean((Theta1c-thetaTRUE)^2); MSE2c[i]=mean((Theta2c-thetaTRUE)^2); MSE3c[i]=mean((Theta3c-thetaTRUE)^2); MSE4c[i]=mean((Theta4c-thetaTRUE)^2)
  Var0[i]=(mean((Theta0-thetaTRUE)^4)-MSE0[i]^2)/Nrep; Var1[i]=(mean((Theta1-thetaTRUE)^4)-MSE1[i]^2)/Nrep; Var2[i]=(mean((Theta2-thetaTRUE)^4)-MSE2[i]^2)/Nrep
  Var3[i]=(mean((Theta3-thetaTRUE)^4)-MSE3[i]^2)/Nrep; Var4[i]=(mean((Theta4-thetaTRUE)^4)-MSE4[i]^2)/Nrep
  Var1c[i]=(mean((Theta1c-thetaTRUE)^4)-MSE1c[i]^2)/Nrep; Var2c[i]=(mean((Theta2c-thetaTRUE)^4)-MSE2c[i]^2)/Nrep
  Var3c[i]=(mean((Theta3c-thetaTRUE)^4)-MSE3c[i]^2)/Nrep; Var4c[i]=(mean((Theta4c-thetaTRUE)^4)-MSE4c[i]^2)/Nrep
}
plot(Rho,n*MSE0,type="l",ylim=c(0,0.5),ylab='n*MSE',xlab=TeX(r"($\rho$)"))
lines(Rho,n*MSE1,col=2)
lines(Rho,n*MSE2,col=3)
lines(Rho,n*MSE3,col=4)
lines(Rho,n*MSE4,col=5)
lines(Rho,n*MSE1c,col=2,lty=2)
lines(Rho,n*MSE2c,col=3,lty=2)
lines(Rho,n*MSE3c,col=4,lty=2)
lines(Rho,n*MSE4c,col=5,lty=2)
lines(Rho,n*MSEmice,lty=2)
for(i in 1:length(Rho)){
  points(Rho[i],n*MSE0[i],pch=1)
  segments(Rho[i],n*MSE0[i]-3*n*sqrt(Var0[i]),Rho[i],n*MSE0[i]+3*n*sqrt(Var0[i]))
  points(Rho[i],n*MSE1[i],pch=1,col=2)
  segments(Rho[i],n*MSE1[i]-3*n*sqrt(Var1[i]),Rho[i],n*MSE1[i]+3*n*sqrt(Var1[i]),col=2)
  points(Rho[i],n*MSE2[i],pch=1,col=3)
  segments(Rho[i],n*MSE2[i]-3*n*sqrt(Var2[i]),Rho[i],n*MSE2[i]+3*n*sqrt(Var2[i]),col=3)
  points(Rho[i],n*MSE3[i],pch=1,col=4)
  segments(Rho[i],n*MSE3[i]-3*n*sqrt(Var3[i]),Rho[i],n*MSE3[i]+3*n*sqrt(Var3[i]),col=4)
  points(Rho[i],n*MSE4[i],pch=1,col=5)
  segments(Rho[i],n*MSE4[i]-3*n*sqrt(Var4[i]),Rho[i],n*MSE4[i]+3*n*sqrt(Var4[i]),col=5)
  points(Rho[i],n*MSE1c[i],pch=2,col=2)
  segments(Rho[i],n*MSE1c[i]-3*n*sqrt(Var1c[i]),Rho[i],n*MSE1c[i]+3*n*sqrt(Var1c[i]),col=2)
  points(Rho[i],n*MSE2c[i],pch=2,col=3)
  segments(Rho[i],n*MSE2c[i]-3*n*sqrt(Var2c[i]),Rho[i],n*MSE2c[i]+3*n*sqrt(Var2c[i]),col=3)
  points(Rho[i],n*MSE3c[i],pch=2,col=4)
  segments(Rho[i],n*MSE3c[i]-3*n*sqrt(Var3c[i]),Rho[i],n*MSE3c[i]+3*n*sqrt(Var3c[i]),col=4)
  points(Rho[i],n*MSE4c[i],pch=2,col=5)
  segments(Rho[i],n*MSE4c[i]-3*n*sqrt(Var4c[i]),Rho[i],n*MSE4c[i]+3*n*sqrt(Var4c[i]),col=5)
  points(Rho[i],n*MSEmice[i],pch=2)
  segments(Rho[i],n*MSEmice[i]-3*n*sqrt(Varmice[i]),Rho[i],n*MSEmice[i]+3*n*sqrt(Varmice[i]))
}
legend(0,0.5,legend=c("CC",TeX(r"($M = 1$)"),TeX(r"($M = 2$)"),TeX(r"($M = 3$)"),TeX(r"($M = 4$)"),
                      TeX(r"(Direct)"),TeX(r"(Cross fit)")),
       col=c("black","red","green","blue","cyan","black","black"), lty=c(1,1,1,1,1,1,2))

Results=rbind(MSE0,MSE1,MSE2,MSE3,MSE4,MSE1c,MSE2c,MSE3c,MSE4c,
              Var0,Var1,Var2,Var3,Var4,Var1c,Var2c,Var3c,Var4c)
save(Results, file="Gaussian600.RData")

n=600; Rho=0.1*(0:10); Nrep=10000
MSE0=Results[1,]; MSE1=Results[2,]; MSE2=Results[3,]; MSE3=Results[4,]; MSE4=Results[5,]; MSE1c=Results[6,]; MSE2c=Results[7,]; MSE3c=Results[8,]; MSE4c=Results[9,]
Var0=Results[10,]; Var1=Results[11,]; Var2=Results[12,]; Var3=Results[13,]; Var4=Results[14,]; Var1c=Results[15,]; Var2c=Results[16,]; Var3c=Results[17,]; Var4c=Results[18,]

#n=200 MICE results, calculated using code below.
# MSEmice=c(0.1328424,0.1355257,0.1374677,0.1342965,0.1351721,0.1374025,0.1396817,0.1385765,0.1406878,0.1434849,0.1830758)/n
# Varmice=c(0.0001775325,0.0001786261,0.0001912602,0.0001841224,0.0001807514,0.0001874732,0.0001991309,0.0001847857,0.0002005651,0.0002101527,0.0003420312)/(n*Nrep)

#n=600 MICE results, calculated using code below.
# MSEmice=c(0.1317722,0.1337938,0.1331833,0.1353186,0.1387210,0.1389396,0.1389456,0.1412965,0.1446840,0.1472614,0.2010326)/n
# Varmice=c(0.00006074245,0.00005998430,0.00005881146,0.00006244006,0.00006318514,0.00006393699,0.00006381509,0.00006806720,0.00007030742,0.00006938662,0.00001358752)/(n*Nrep)


########### Setting (ii) ############

n=600; lambda1=10; lambda2=10
Mu=2*(0:10)
h=0.1
Nrep=10000

MSE0=rep(0,length(Mu)); Var0=rep(0,length(Mu))
MSE1=rep(0,length(Mu)); Var1=rep(0,length(Mu))
MSE2=rep(0,length(Mu)); Var2=rep(0,length(Mu))
MSE3=rep(0,length(Mu)); Var3=rep(0,length(Mu))
MSE4=rep(0,length(Mu)); Var4=rep(0,length(Mu))
MSE1c=rep(0,length(Mu)); Var1c=rep(0,length(Mu))
MSE2c=rep(0,length(Mu)); Var2c=rep(0,length(Mu))
MSE3c=rep(0,length(Mu)); Var3c=rep(0,length(Mu))
MSE4c=rep(0,length(Mu)); Var4c=rep(0,length(Mu))
for(i in 1:length(Mu)){
  print(i)
  mu=Mu[i]
  
  n0=1000000
  X=rMvdc(n0,mvdc(copula=claytonCopula(3),margins=c("beta","beta"),paramMargins = list(list(shape1=1,shape2=1),list(shape1=1,shape2=1))))
  Y=rbinom(n0,1,prob=exp(mu*(X[,1]+X[,2]-1))/(1+exp(mu*(X[,1]+X[,2]-1))))
  thetaTRUE=mean(Y)
  
  Theta0=rep(0,Nrep)
  Theta1=rep(0,Nrep)
  Theta2=rep(0,Nrep)
  Theta3=rep(0,Nrep)
  Theta4=rep(0,Nrep)
  Theta1c=rep(0,Nrep)
  Theta2c=rep(0,Nrep)
  Theta3c=rep(0,Nrep)
  Theta4c=rep(0,Nrep)
  for(nrep in 1:Nrep){
    X=rMvdc(n,mvdc(copula=claytonCopula(3),margins=c("beta","beta"),paramMargins = list(list(shape1=1,shape2=1),list(shape1=1,shape2=1))))
    Y=rbinom(n,1,prob=exp(mu*(X[,1]+X[,2]-1))/(1+exp(mu*(X[,1]+X[,2]-1))))
    
    X1=runif(n*lambda1); X2=runif(n*lambda2)
    
    Theta0[nrep]=mean(Y) 
    Theta1[nrep]=thetaDirect(X,Y,X1,X2,1,h)
    Theta2[nrep]=thetaDirect(X,Y,X1,X2,2,h)
    Theta3[nrep]=thetaDirect(X,Y,X1,X2,3,h)
    Theta4[nrep]=thetaDirect(X,Y,X1,X2,4,h)
    Theta1c[nrep]=thetac(X,Y,X1,X2,1,h)
    Theta2c[nrep]=thetac(X,Y,X1,X2,2,h)
    Theta3c[nrep]=thetac(X,Y,X1,X2,3,h)
    Theta4c[nrep]=thetac(X,Y,X1,X2,4,h)
    
  }
  MSE0[i]=mean((Theta0-thetaTRUE)^2); MSE1[i]=mean((Theta1-thetaTRUE)^2); MSE2[i]=mean((Theta2-thetaTRUE)^2); MSE3[i]=mean((Theta3-thetaTRUE)^2); MSE4[i]=mean((Theta4-thetaTRUE)^2)
  MSE1c[i]=mean((Theta1c-thetaTRUE)^2); MSE2c[i]=mean((Theta2c-thetaTRUE)^2); MSE3c[i]=mean((Theta3c-thetaTRUE)^2); MSE4c[i]=mean((Theta4c-thetaTRUE)^2)
  Var0[i]=(mean((Theta0-thetaTRUE)^4)-MSE0[i]^2)/Nrep; Var1[i]=(mean((Theta1-thetaTRUE)^4)-MSE1[i]^2)/Nrep; Var2[i]=(mean((Theta2-thetaTRUE)^4)-MSE2[i]^2)/Nrep
  Var3[i]=(mean((Theta3-thetaTRUE)^4)-MSE3[i]^2)/Nrep; Var4[i]=(mean((Theta4-thetaTRUE)^4)-MSE4[i]^2)/Nrep
  Var1c[i]=(mean((Theta1c-thetaTRUE)^4)-MSE1c[i]^2)/Nrep; Var2c[i]=(mean((Theta2c-thetaTRUE)^4)-MSE2c[i]^2)/Nrep
  Var3c[i]=(mean((Theta3c-thetaTRUE)^4)-MSE3c[i]^2)/Nrep; Var4c[i]=(mean((Theta4c-thetaTRUE)^4)-MSE4c[i]^2)/Nrep
}
plot(Mu,n*MSE0,type="l",ylim=c(0,0.4),ylab='n*MSE',xlab=TeX(r"($\mu$)"))
lines(Mu,n*MSE1,col=2)
lines(Mu,n*MSE2,col=3)
lines(Mu,n*MSE3,col=4)
lines(Mu,n*MSE4,col=5)
lines(Mu,n*MSE1c,col=2,lty=2)
lines(Mu,n*MSE2c,col=3,lty=2)
lines(Mu,n*MSE3c,col=4,lty=2)
lines(Mu,n*MSE4c,col=5,lty=2)
lines(Mu,n*MSEmice,lty=2)
for(i in 1:length(Mu)){
  points(Mu[i],n*MSE0[i],pch=1)
  segments(Mu[i],n*MSE0[i]-3*n*sqrt(Var0[i]),Mu[i],n*MSE0[i]+3*n*sqrt(Var0[i]))
  points(Mu[i],n*MSE1[i],pch=1,col=2)
  segments(Mu[i],n*MSE1[i]-3*n*sqrt(Var1[i]),Mu[i],n*MSE1[i]+3*n*sqrt(Var1[i]),col=2)
  points(Mu[i],n*MSE2[i],pch=1,col=3)
  segments(Mu[i],n*MSE2[i]-3*n*sqrt(Var2[i]),Mu[i],n*MSE2[i]+3*n*sqrt(Var2[i]),col=3)
  points(Mu[i],n*MSE3[i],pch=1,col=4)
  segments(Mu[i],n*MSE3[i]-3*n*sqrt(Var3[i]),Mu[i],n*MSE3[i]+3*n*sqrt(Var3[i]),col=4)
  points(Mu[i],n*MSE4[i],pch=1,col=5)
  segments(Mu[i],n*MSE4[i]-3*n*sqrt(Var4[i]),Mu[i],n*MSE4[i]+3*n*sqrt(Var4[i]),col=5)
  points(Mu[i],n*MSE1c[i],pch=2,col=2)
  segments(Mu[i],n*MSE1c[i]-3*n*sqrt(Var1c[i]),Mu[i],n*MSE1c[i]+3*n*sqrt(Var1c[i]),col=2)
  points(Mu[i],n*MSE2c[i],pch=2,col=3)
  segments(Mu[i],n*MSE2c[i]-3*n*sqrt(Var2c[i]),Mu[i],n*MSE2c[i]+3*n*sqrt(Var2c[i]),col=3)
  points(Mu[i],n*MSE3c[i],pch=2,col=4)
  segments(Mu[i],n*MSE3c[i]-3*n*sqrt(Var3c[i]),Mu[i],n*MSE3c[i]+3*n*sqrt(Var3c[i]),col=4)
  points(Mu[i],n*MSE4c[i],pch=2,col=5)
  segments(Mu[i],n*MSE4c[i]-3*n*sqrt(Var4c[i]),Mu[i],n*MSE4c[i]+3*n*sqrt(Var4c[i]),col=5)
  points(Mu[i],n*MSEmice[i],pch=2)
  segments(Mu[i],n*MSEmice[i]-3*n*sqrt(Varmice[i]),Mu[i],n*MSEmice[i]+3*n*sqrt(Varmice[i]))
}
legend(15,0.5,legend=c("CC",TeX(r"(Direct with $M = 1$)"),TeX(r"(Direct with $M = 2$)"),TeX(r"(Direct with $M = 3$)"),TeX(r"(Direct with $M = 4$)"),
                       TeX(r"(Cross fit with $M = 1$)"),TeX(r"(Cross fit with $M = 2$)")),
       col=c("black","red","green","blue","cyan","red","green"), lty=c(1,1,1,1,1,2,2))

Results=rbind(MSE0,MSE1,MSE2,MSE3,MSE4,MSE1c,MSE2c,MSE3c,MSE4c,
              Var0,Var1,Var2,Var3,Var4,Var1c,Var2c,Var3c,Var4c)
save(Results, file="Clayton600.RData")

n=600; Mu=2*(0:10); Nrep=10000
MSE0=Results[1,]; MSE1=Results[2,]; MSE2=Results[3,]; MSE3=Results[4,]; MSE4=Results[5,]; MSE1c=Results[6,]; MSE2c=Results[7,]; MSE3c=Results[8,]; MSE4c=Results[9,]
Var0=Results[10,]; Var1=Results[11,]; Var2=Results[12,]; Var3=Results[13,]; Var4=Results[14,]; Var1c=Results[15,]; Var2c=Results[16,]; Var3c=Results[17,]; Var4c=Results[18,]

#n=200 MICE results, calculated using code below.
# MSEmice=c(1.1031262,0.2705292,0.2011489,0.2044763,0.1992181,0.1844506,0.1847066,0.1850245,0.1777426,0.1723417,0.1639580)/n
# Varmice=c(0.0222337506,0.0007268991,0.0004141692,0.0004308469,0.0004158989,0.0003293372,0.0003518238,0.0003383010,0.0003103790,0.0002865324,0.0002566747)/(n*Nrep)

#n=600 MICE results, calculated using code below.
# MSEmice=c(2.4448079,0.2808012,0.2194846,0.2551857,0.2844701,0.2826332,0.2656920,0.2556021,0.2514752,0.2467162,0.2562607)/n
# Varmice=c(0.0462808168,0.0002496405,0.0001594218,0.0002188250,0.0002821502,0.0002566802,0.0002311650,0.0001980954,0.0001875369,0.0001799696,0.0001886138)/(n*Nrep)






########### Setting (iii) ############

n=600; lambda1=10; lambda2=10
sigY=0.3; Rho=0.1*(0:10)
h=0.1
Nrep=10000

MSE0=rep(0,length(Rho)); Var0=rep(0,length(Rho))
MSE1=rep(0,length(Rho)); Var1=rep(0,length(Rho))
MSE2=rep(0,length(Rho)); Var2=rep(0,length(Rho))
MSE3=rep(0,length(Rho)); Var3=rep(0,length(Rho))
MSE4=rep(0,length(Rho)); Var4=rep(0,length(Rho))
MSE1c=rep(0,length(Rho)); Var1c=rep(0,length(Rho))
MSE2c=rep(0,length(Rho)); Var2c=rep(0,length(Rho))
MSE3c=rep(0,length(Rho)); Var3c=rep(0,length(Rho))
MSE4c=rep(0,length(Rho)); Var4c=rep(0,length(Rho))
MSEmice=rep(0,length(Rho)); Varmice=rep(0,length(Rho))
for(i in 1:length(Rho)){
  print(i)
  
  n0=1000000
  X=pnorm(rmvnorm(n0,sigma=matrix(c(1,-Rho[i],-Rho[i],1),2,2)))
  Y=-5*X[,1]*X[,2]+rnorm(n0,sd=sigY)
  thetaTRUE=mean(Y)
  
  
  Theta0=rep(0,Nrep)
  Theta1=rep(0,Nrep)
  Theta2=rep(0,Nrep)
  Theta3=rep(0,Nrep)
  Theta4=rep(0,Nrep)
  Theta1c=rep(0,Nrep)
  Theta2c=rep(0,Nrep)
  Theta3c=rep(0,Nrep)
  Theta4c=rep(0,Nrep)
  Thetamice=rep(0,Nrep)
  for(nrep in 1:Nrep){
    X=pnorm(rmvnorm(n,sigma=matrix(c(1,-Rho[i],-Rho[i],1),2,2)))
    Y=-5*X[,1]*X[,2]+rnorm(n,sd=sigY)
    
    X1=runif(n*lambda1); X2=runif(n*lambda2)
    
    Theta0[nrep]=mean(Y) 
    Theta1[nrep]=thetaDirect(X,Y,X1,X2,1,h)
    Theta2[nrep]=thetaDirect(X,Y,X1,X2,2,h)
    Theta3[nrep]=thetaDirect(X,Y,X1,X2,3,h)
    Theta4[nrep]=thetaDirect(X,Y,X1,X2,4,h)
    Theta1c[nrep]=thetac(X,Y,X1,X2,1,h)
    Theta2c[nrep]=thetac(X,Y,X1,X2,2,h)
    Theta3c[nrep]=thetac(X,Y,X1,X2,3,h)
    Theta4c[nrep]=thetac(X,Y,X1,X2,4,h)
    
    Xcomp=cbind(X,Y)
    X1=cbind(X1,rep(NA,n*lambda1),rep(NA,n*lambda1))
    X2=cbind(rep(NA,n*lambda2),X2,rep(NA,n*lambda2))
    X=rbind(Xcomp,X1,X2)
    
    tempData=mice(X,printFlag=FALSE)
    estimates=with(tempData,mean(Y))$analyses
    Thetamice[nrep]=mean(unlist(estimates))
    
  }
  MSE0[i]=mean((Theta0-thetaTRUE)^2); MSE1[i]=mean((Theta1-thetaTRUE)^2); MSE2[i]=mean((Theta2-thetaTRUE)^2); MSE3[i]=mean((Theta3-thetaTRUE)^2); MSE4[i]=mean((Theta4-thetaTRUE)^2)
  MSE1c[i]=mean((Theta1c-thetaTRUE)^2); MSE2c[i]=mean((Theta2c-thetaTRUE)^2); MSE3c[i]=mean((Theta3c-thetaTRUE)^2); MSE4c[i]=mean((Theta4c-thetaTRUE)^2)
  MSEmice[i]=mean((Thetamice-thetaTRUE)^2)
  Var0[i]=(mean((Theta0-thetaTRUE)^4)-MSE0[i]^2)/Nrep; Var1[i]=(mean((Theta1-thetaTRUE)^4)-MSE1[i]^2)/Nrep; Var2[i]=(mean((Theta2-thetaTRUE)^4)-MSE2[i]^2)/Nrep
  Var3[i]=(mean((Theta3-thetaTRUE)^4)-MSE3[i]^2)/Nrep; Var4[i]=(mean((Theta4-thetaTRUE)^4)-MSE4[i]^2)/Nrep
  Var1c[i]=(mean((Theta1c-thetaTRUE)^4)-MSE1c[i]^2)/Nrep; Var2c[i]=(mean((Theta2c-thetaTRUE)^4)-MSE2c[i]^2)/Nrep
  Var3c[i]=(mean((Theta3c-thetaTRUE)^4)-MSE3c[i]^2)/Nrep; Var4c[i]=(mean((Theta4c-thetaTRUE)^4)-MSE4c[i]^2)/Nrep
  Varmice[i]=(mean((Thetamice-thetaTRUE)^4)-MSEmice[i]^2)/Nrep
}
plot(Rho,n*MSE0,type="l",ylim=c(0,1.5),ylab='n*MSE',xlab=TeX(r"($\rho$)"))
lines(Rho,n*MSE1,col=2)
lines(Rho,n*MSE2,col=3)
lines(Rho,n*MSE3,col=4)
lines(Rho,n*MSE4,col=5)
lines(Rho,n*MSE1c,col=2,lty=2)
lines(Rho,n*MSE2c,col=3,lty=2)
lines(Rho,n*MSE3c,col=4,lty=2)
lines(Rho,n*MSE4c,col=5,lty=2)
lines(Rho,n*MSEmice,col=1,lty=2)
for(i in 1:length(Rho)){
  points(Rho[i],n*MSE0[i],pch=1)
  segments(Rho[i],n*MSE0[i]-3*n*sqrt(Var0[i]),Rho[i],n*MSE0[i]+3*n*sqrt(Var0[i]))
  points(Rho[i],n*MSE1[i],pch=1,col=2)
  segments(Rho[i],n*MSE1[i]-3*n*sqrt(Var1[i]),Rho[i],n*MSE1[i]+3*n*sqrt(Var1[i]),col=2)
  points(Rho[i],n*MSE2[i],pch=1,col=3)
  segments(Rho[i],n*MSE2[i]-3*n*sqrt(Var2[i]),Rho[i],n*MSE2[i]+3*n*sqrt(Var2[i]),col=3)
  points(Rho[i],n*MSE3[i],pch=1,col=4)
  segments(Rho[i],n*MSE3[i]-3*n*sqrt(Var3[i]),Rho[i],n*MSE3[i]+3*n*sqrt(Var3[i]),col=4)
  points(Rho[i],n*MSE4[i],pch=1,col=5)
  segments(Rho[i],n*MSE4[i]-3*n*sqrt(Var4[i]),Rho[i],n*MSE4[i]+3*n*sqrt(Var4[i]),col=5)
  points(Rho[i],n*MSE1c[i],pch=2,col=2)
  segments(Rho[i],n*MSE1c[i]-3*n*sqrt(Var1c[i]),Rho[i],n*MSE1c[i]+3*n*sqrt(Var1c[i]),col=2)
  points(Rho[i],n*MSE2c[i],pch=2,col=3)
  segments(Rho[i],n*MSE2c[i]-3*n*sqrt(Var2c[i]),Rho[i],n*MSE2c[i]+3*n*sqrt(Var2c[i]),col=3)
  points(Rho[i],n*MSE3c[i],pch=2,col=4)
  segments(Rho[i],n*MSE3c[i]-3*n*sqrt(Var3c[i]),Rho[i],n*MSE3c[i]+3*n*sqrt(Var3c[i]),col=4)
  points(Rho[i],n*MSE4c[i],pch=2,col=5)
  segments(Rho[i],n*MSE4c[i]-3*n*sqrt(Var4c[i]),Rho[i],n*MSE4c[i]+3*n*sqrt(Var4c[i]),col=5)
  
  points(Rho[i],n*MSEmice[i],pch=2)
  segments(Rho[i],n*MSEmice[i]-3*n*sqrt(Varmice[i]),Rho[i],n*MSEmice[i]+3*n*sqrt(Varmice[i]))
}
legend(0,1,legend=c("CC",TeX(r"($M = 1$)"),TeX(r"($M = 2$)"),TeX(r"($M = 3$)"),TeX(r"($M = 4$)"),
                    TeX(r"(Direct)"),TeX(r"(Cross fit)")),
       col=c("black","red","green","blue","cyan","black","black"), lty=c(1,1,1,1,1,1,2))


Results=rbind(MSE0,MSE1,MSE2,MSE3,MSE4,MSE1c,MSE2c,MSE3c,MSE4c,MSEmice,
              Var0,Var1,Var2,Var3,Var4,Var1c,Var2c,Var3c,Var4c,Varmice)
save(Results, file="Quadratic600.RData")

n=600; Rho=0.1*(0:10); Nrep=10000
MSE0=Results[1,]; MSE1=Results[2,]; MSE2=Results[3,]; MSE3=Results[4,]; MSE4=Results[5,]; MSE1c=Results[6,]; MSE2c=Results[7,]; MSE3c=Results[8,]; MSE4c=Results[9,]; MSEmice=Results[10,]
Var0=Results[11,]; Var1=Results[12,]; Var2=Results[13,]; Var3=Results[14,]; Var4=Results[15,]; Var1c=Results[16,]; Var2c=Results[17,]; Var3c=Results[18,]; Var4c=Results[19,]; Varmice=Results[20,]







##################### MICE ############################


########### Setting (i) ############

n=200; lambda1=10; lambda2=10
sigY=0.3; Rho=0.1*(0:10)
Nrep=10000

bias=rep(0,length(Rho))
variance=rep(0,length(Rho))
MSE=rep(0,length(Rho))
Var=rep(0,length(Rho))
for(i in 1:length(Rho)){
  print(i)
  Theta=rep(0,Nrep)
  for(nrep in 1:Nrep){
    X=pnorm(rmvnorm(n,sigma=matrix(c(1,-Rho[i],-Rho[i],1),2,2)))
    Y=X[,1]-X[,2]+rnorm(n,sd=sigY)
    Xcomp=cbind(X,Y)
    X1=cbind(runif(n*lambda1),rep(NA,n*lambda1),rep(NA,n*lambda1))
    X2=cbind(rep(NA,n*lambda2),runif(n*lambda2),rep(NA,n*lambda2))
    X=rbind(Xcomp,X1,X2)
    
    tempData=mice(X,printFlag=FALSE)
    estimates=with(tempData,mean(Y))$analyses
    Theta[nrep]=mean(unlist(estimates))
    
    thetaTRUE=0
    
  }
  bias[i]=mean(Theta)-thetaTRUE
  variance[i]=var(Theta)
  MSE[i]=mean((Theta-thetaTRUE)^2)
  Var[i]=(mean((Theta-thetaTRUE)^4)-MSE[i]^2)/Nrep
}
plot(Rho,n*MSE,type="l",ylim=c(0,0.4),ylab='n*MSE',xlab=TeX(r"($\rho$)"))
lines(Rho,n*bias^2,col=2)
lines(Rho,n*variance,col=3)
for(i in 1:length(Rho)){
  points(Rho[i],n*MSE[i],pch=1)
  segments(Rho[i],n*MSE[i]-3*n*sqrt(Var[i]),Rho[i],n*MSE[i]+3*n*sqrt(Var[i]))
}

#n=200 results
# MSE=c(0.1328424,0.1355257,0.1374677,0.1342965,0.1351721,0.1374025,0.1396817,0.1385765,0.1406878,0.1434849,0.1830758)/n
# Var=c(0.0001775325,0.0001786261,0.0001912602,0.0001841224,0.0001807514,0.0001874732,0.0001991309,0.0001847857,0.0002005651,0.0002101527,0.0003420312)/(n*Nrep)


#n=600 results
# n*MSE=c(0.1317722,0.1337938,0.1331833,0.1353186,0.1387210,0.1389396,0.1389456,0.1412965,0.1446840,0.1472614,0.2010326)
# n*Nrep*Var=c(0.00006074245,0.00005998430,0.00005881146,0.00006244006,0.00006318514,0.00006393699,0.00006381509,0.00006806720,0.00007030742,0.00006938662,0.00001358752)




########### Setting (ii) ############

n=200; lambda1=10; lambda2=10
Mu=2*(0:10)
Nrep=10000

bias=rep(0,length(Mu))
variance=rep(0,length(Mu))
MSE=rep(0,length(Mu))
Var=rep(0,length(Mu))
for(i in 1:length(Mu)){
  print(i)
  mu=Mu[i]
  
  n0=1000000
  X=rMvdc(n0,mvdc(copula=claytonCopula(3),margins=c("beta","beta"),paramMargins = list(list(shape1=1,shape2=1),list(shape1=1,shape2=1))))
  Y=rbinom(n0,1,prob=exp(mu*(X[,1]+X[,2]-1))/(1+exp(mu*(X[,1]+X[,2]-1))))
  thetaTRUE=mean(Y)
  
  Theta=rep(0,Nrep)
  for(nrep in 1:Nrep){
    X=rMvdc(n,mvdc(copula=claytonCopula(3),margins=c("beta","beta"),paramMargins = list(list(shape1=1,shape2=1),list(shape1=1,shape2=1))))
    Y=rbinom(n,1,prob=exp(mu*(X[,1]+X[,2]-1))/(1+exp(mu*(X[,1]+X[,2]-1))))
    
    Xcomp=cbind(X,Y)
    X1=cbind(runif(n*lambda1),rep(NA,n*lambda1),rep(NA,n*lambda1))
    X2=cbind(rep(NA,n*lambda2),runif(n*lambda2),rep(NA,n*lambda2))
    X=rbind(Xcomp,X1,X2)
    
    tempData=mice(X,printFlag=FALSE)
    estimates=with(tempData,mean(Y))$analyses
    Theta[nrep]=mean(unlist(estimates))
  }
  bias[i]=mean(Theta)-thetaTRUE
  variance[i]=var(Theta)
  MSE[i]=mean((Theta-thetaTRUE)^2)
  Var[i]=(mean((Theta-thetaTRUE)^4)-MSE[i]^2)/Nrep
}
plot(Mu,n*MSE,type="l",ylim=c(0,0.4),ylab='n*MSE',xlab=TeX(r"($\mu$)"))
lines(Mu,n*bias^2,col=2)
lines(Mu,n*variance,col=3)
for(i in 1:length(Mu)){
  points(Mu[i],n*MSE[i],pch=1)
  segments(Mu[i],n*MSE[i]-3*n*sqrt(Var[i]),Mu[i],n*MSE[i]+3*n*sqrt(Var[i]))
}


#n=200 results
# MSE=c(1.1031262,0.2705292,0.2011489,0.2044763,0.1992181,0.1844506,0.1847066,0.1850245,0.1777426,0.1723417,0.1639580)/n
# Var=c(0.0222337506,0.0007268991,0.0004141692,0.0004308469,0.0004158989,0.0003293372,0.0003518238,0.0003383010,0.0003103790,0.0002865324,0.0002566747)/(n*Nrep)

#n=600 results
# MSE=c(2.4448079,0.2808012,0.2194846,0.2551857,0.2844701,0.2826332,0.2656920,0.2556021,0.2514752,0.2467162,0.2562607)/n
# Var=c(0.0462808168,0.0002496405,0.0001594218,0.0002188250,0.0002821502,0.0002566802,0.0002311650,0.0001980954,0.0001875369,0.0001799696,0.0001886138)/(n*Nrep)



########### Setting (iii) ############

n=200; lambda1=10; lambda2=10
sigY=0.3; Rho=0.1*(0:10)
Nrep=100

bias=rep(0,length(Rho))
variance=rep(0,length(Rho))
MSE=rep(0,length(Rho))
Var=rep(0,length(Rho))
for(i in 1:length(Rho)){
  print(i)
  
  n0=1000000
  X=pnorm(rmvnorm(n0,sigma=matrix(c(1,-Rho[i],-Rho[i],1),2,2)))
  Y=-5*X[,1]*X[,2]+rnorm(n0,sd=sigY)
  thetaTRUE=mean(Y)
  
  Theta=rep(0,Nrep)
  for(nrep in 1:Nrep){
    X=pnorm(rmvnorm(n,sigma=matrix(c(1,-Rho[i],-Rho[i],1),2,2)))
    Y=-5*X[,1]*X[,2]+rnorm(n,sd=sigY)
    
    Xcomp=cbind(X,Y)
    X1=cbind(runif(n*lambda1),rep(NA,n*lambda1),rep(NA,n*lambda1))
    X2=cbind(rep(NA,n*lambda2),runif(n*lambda2),rep(NA,n*lambda2))
    X=rbind(Xcomp,X1,X2)
    
    tempData=mice(X,printFlag=FALSE)
    estimates=with(tempData,mean(Y))$analyses
    Theta[nrep]=mean(unlist(estimates))
    
    
  }
  bias[i]=mean(Theta)-thetaTRUE
  variance[i]=var(Theta)
  MSE[i]=mean((Theta-thetaTRUE)^2)
  Var[i]=(mean((Theta-thetaTRUE)^4)-MSE[i]^2)/Nrep
}
Rho=Rho[-11]; MSE=MSE[-11]; Var=Var[-11]; bias=bias[-11]; variance=variance[-11]
plot(Rho,n*MSE,type="l",ylim=c(0,10),ylab='n*MSE',xlab=TeX(r"($\rho$)"))
lines(Rho,n*bias^2,col=2)
lines(Rho,n*variance,col=3)
for(i in 1:length(Rho)){
  points(Rho[i],n*MSE[i],pch=1)
  segments(Rho[i],n*MSE[i]-3*n*sqrt(Var[i]),Rho[i],n*MSE[i]+3*n*sqrt(Var[i]))
}
