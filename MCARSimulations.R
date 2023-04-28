library(MCARtest)
library(PKLMtest)

# Function to convert p-way contingency tables with sum n to nxp data matrices
tableinv <- function(x) {  
  stopifnot(is.table(x))
  obs <- as.data.frame(x)[rep(1:prod(dim(x)),c(round(x))),-length(dim(x))-1]
  rownames(obs) <- NULL; obs
}

# Parameters defining the 3-dimensional setting
bS=matrix(c(1,1,0, 1,0,1, 0,1,1),byrow=TRUE,ncol=3) # Observing all pairs of three variables
r=2; s=2; t=2; M=c(r,s,t) # The alphabet sizes of our three variables

alpha=0.05 # Desired significance level (Type I error rate)
Pd21=0.005*(50+0:25) # Range of values of $p_{\bullet 2 1}$
n0=200 # Number of complete cases 
n=200 # Number of observations of each bivariate distribution
Nrep=50 # Number of repetitions of the experiment
B=99 # Number of bootstrap samples for our test

Sim3d=function(bS,M,alpha,Pd21,n0,n,B,Nrep,CompatTest=TRUE,Fuchs=FALSE,PKLM=FALSE){
  if(CompatTest==TRUE) Power=rep(0,length(Pd21))
  if(Fuchs==TRUE) PowerFuchs=rep(0,length(Pd21))
  if(PKLM==TRUE) PowerPKLM=rep(0,length(Pd21))
  Rpop=rep(0,length(Pd21))
  for(index in 1:length(Pd21)){
    pd21=Pd21[index]
    r=M[1]; s=M[2]; t=M[3]
    pS=rep(0,r*s+s*t+r*t)
    
    pS[row_index(bS,M,3,c(2,1))]=pd21
    pS[row_index(bS,M,3,c(1,1))]=1/2-pd21
    pS[row_index(bS,M,3,c(2,2))]=1/2-pd21
    pS[row_index(bS,M,3,c(1,2))]=pd21
    
    for(i in 1:r){
      pS[row_index(bS,M,1,c(i,1))]=(1+(-1)^i)/(2*r)
      pS[row_index(bS,M,1,c(i,2))]=(1-(-1)^i)/(2*r)
      
      pS[row_index(bS,M,2,c(i,1))]=1/(2*r)
      pS[row_index(bS,M,2,c(i,2))]=1/(2*r)
    }
    
    P12=pS[1:(2*r)]; P13=pS[(2*r+1):(4*r)]; P23=pS[-(1:(4*r))]
    Rpop[index]=Rindex(pS,bS,M)
    
    p0=array(0.25/r,dim=c(r,2,2))-(0.25/r)*(( (-1)^(1:r) %o% (-1)^(1:2)) %o% rep(1,2))
    
    X12=t(rmultinom(Nrep,size=n,prob=P12)/n); X13=t(rmultinom(Nrep,size=n,prob=P13)/n); X23=t(rmultinom(Nrep,size=n,prob=P23)/n)
    pSh=cbind(X12,X13,X23)
    
    if(CompatTest==TRUE) Rej=rep(0,Nrep)
    if(Fuchs==TRUE)  RejFuchs=rep(0,Nrep)
    if(PKLM==TRUE) RejPKLM=rep(0,Nrep)
    for(nrep in 1:Nrep){
      
      if(CompatTest==TRUE){
        pval=ProjectionTest(pSh[nrep,],rep(n,3),bS,M,B)[[1]]
        Rej[nrep]=(pval<=alpha)
      }
      
      if(Fuchs==TRUE){
        p0h=array(rmultinom(1,n0,p0),dim=M)/n0
        RejFuchs[nrep]=(FuchsTest(p0h,n0,pSh[nrep,],rep(n,3),bS,M,100)<=alpha)
      }
      
      if(PKLM==TRUE){
        x0=array(rmultinom(1,n0,p0),dim=M)
        x12=n*matrix(X12[nrep,],nrow=M[1])
        x23=n*matrix(X23[nrep,],nrow=M[2])
        x13=n*matrix(X13[nrep,],nrow=M[1])
        
        data0=as.matrix(tableinv(as.table(x0)))
        data12=as.matrix(tableinv(as.table(x12))); data12=cbind(data12,rep(NA,n))
        data23=as.matrix(tableinv(as.table(x23))); data23=cbind(rep(NA,n),data23)
        data13=as.matrix(tableinv(as.table(x13))); data13=cbind(data13[,1],rep(NA,n),data13[,2])
        data=rbind(data0,data12,data23,data13)
        for(i in 1:r){
          data[data==LETTERS[i]]=i
        }
        data=matrix(as.numeric(data),ncol=3)
        
        pval=PKLMtest(data)
        RejPKLM[nrep]=(pval<=alpha)
      }
      
    }
    
    
    if(CompatTest==TRUE) Power[index]=sum(Rej)/Nrep
    if(Fuchs==TRUE) PowerFuchs[index]=sum(RejFuchs)/Nrep
    if(PKLM==TRUE) PowerPKLM[index]=sum(RejPKLM)/Nrep
  }
  out=list(Rpop)
  if(CompatTest==TRUE) out=c(out,list(Power))
  if(Fuchs==TRUE) out=c(out,list(PowerFuchs))
  if(PKLM==TRUE) out=c(out,list(PowerPKLM))
  return(out)
}

plot(out[[1]],out[[2]],type="l")
lines(out[[1]],out[[3]],col=2)
lines(out[[1]],out[[4]],col=3)




# Parameters defining the 5-dimensional setting

bS=matrix(c(1,1,1,1,0, 1,1,1,0,1, 1,1,0,1,1, 1,0,1,1,1, 0,1,1,1,1),byrow=TRUE,ncol=5)
# Observing all size-four subsets of our five variables
r=2; M=rep(r,5) # Each variable has alphabet size r

alpha=0.05 # Desired significance level (Type I error rate)
Epsilon=0.01*(20:35) # Range of values of $\epsilon$
N0=c(25,50,100,200) # Numbers of complete cases for Fuchs' test 
n=500 # Number of observations of each 4-dimensional distribution
Nrep=50 # Number of repetitions of the experiment
B=99 # Number of bootstrap samples for our test

perturb=(-1)^(1:r) %o% (-1)^(1:r) %o% (-1)^(1:r) %o% (-1)^(1:r)
Power=rep(0,length(Epsilon))
PowerFuchs1=rep(0,length(Epsilon))
PowerFuchs2=rep(0,length(Epsilon))
PowerFuchs3=rep(0,length(Epsilon))
PowerFuchs4=rep(0,length(Epsilon))
Rpop=rep(0,length(Epsilon))
for(index in 1:length(Epsilon)){
  epsilon=Epsilon[index]
  pS=c(rep(1+epsilon*as.vector(perturb),4)/r^4, (1-epsilon*as.vector(perturb))/r^4)
  p0=array(RindexDual(pS,bS,M)[[2]],dim=M); p0=p0/sum(p0)
  
  PosData=t(rmultinom(4*Nrep,size=n,prob=(1+epsilon*as.vector(perturb))/r^4)/n)
  NegData=t(rmultinom(Nrep,size=n,prob=(1-epsilon*as.vector(perturb))/r^4)/n)
  pSh=cbind(PosData[1:Nrep,],PosData[(Nrep+1):(2*Nrep),],PosData[(2*Nrep+1):(3*Nrep),],PosData[(3*Nrep+1):(4*Nrep),],NegData)
  
  
  Rej=rep(0,Nrep)
  RejFuchs1=rep(0,Nrep)
  RejFuchs2=rep(0,Nrep)
  RejFuchs3=rep(0,Nrep)
  RejFuchs4=rep(0,Nrep)
  for(nrep in 1:Nrep){
    pval=ProjectionTest(pSh[nrep,],rep(n,5),bS,M,B)[[1]]
    Rej[nrep]=(pval<=0.05)
    
    p0h=array(rmultinom(1,N0[1],prob=p0),dim=M)/N0[1]
    pval=FuchsTest(p0h,N0[1],pSh[nrep,],rep(n,5),bS,M,100)
    RejFuchs1[nrep]=(pval<=0.05)
    
    p0h=array(rmultinom(1,N0[2],prob=p0),dim=M)/N0[2]
    pval=FuchsTest(p0h,N0[2],pSh[nrep,],rep(n,5),bS,M,100)
    RejFuchs2[nrep]=(pval<=0.05)
    
    p0h=array(rmultinom(1,N0[3],prob=p0),dim=M)/N0[3]
    pval=FuchsTest(p0h,N0[3],pSh[nrep,],rep(n,5),bS,M,100)
    RejFuchs3[nrep]=(pval<=0.05)
    
    p0h=array(rmultinom(1,N0[4],prob=p0),dim=M)/N0[4]
    pval=FuchsTest(p0h,N0[4],pSh[nrep,],rep(n,5),bS,M,100)
    RejFuchs4[nrep]=(pval<=0.05)
  }
  Power[index]=sum(Rej)/Nrep
  PowerFuchs1[index]=sum(RejFuchs1)/Nrep
  PowerFuchs2[index]=sum(RejFuchs2)/Nrep
  PowerFuchs3[index]=sum(RejFuchs3)/Nrep
  PowerFuchs4[index]=sum(RejFuchs4)/Nrep
  print(c(epsilon,Power[index],PowerFuchs1[index],PowerFuchs2[index],PowerFuchs3[index],PowerFuchs4[index]))
}
plot(Epsilon,Power,type="l",ylim=c(0,1))
lines(Epsilon,PowerFuchs1,col=2)
lines(Epsilon,PowerFuchs2,col=3)
lines(Epsilon,PowerFuchs3,col=4)
lines(Epsilon,PowerFuchs4,col=5)
abline(h=0.05,col=3)

