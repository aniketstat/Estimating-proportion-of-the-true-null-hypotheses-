## Generating simulated dataset for fixed choice of m=no. of hypotheses=1000, b=no. of blocks
## in the dispersion matrix "\Sigma"=100, r=order of blocks dispersion matrix with AR structure
## =10, \nu=df of chi-square variate from which \sigma_i^2's are generated. Flexible argument: 
## pi0=proportion of true null hypotheses, n=sample size available corresponding to each test, u.l=upper limit 
## of the uniform distribution from which non-null (more than zero) \mu values are generated, 
## l.l= lower limit of the uniform distribution from which non-null (more than zero) \mu values 
## are generated, right.alloc= proportion of non-null \mu-values greater than zero. rho=parameter 
## in the auto-regressive structure. 

## Attaching packages
library(Matrix) 
library(MASS)
library(cp4p)
## Code for the function
simul.details=function(pi0,n,u.l,u.u,right.alloc,rho)
{
## Fixing m,b,r:
m=1000
b=100
r=m/b

## Generating the Dispersion matrix:
sig.rho=matrix(0,ncol=b,nrow=b)
for(k1 in 1:b)
{
for(k2 in 1:b)
{
sig.rho[k1,k2]=rho^(abs(k1-k2))
}
}
sigsq=rchisq(r,10)/10

mt.sigma=bdiag(sigsq[1]*sig.rho,sigsq[2]*sig.rho,sigsq[3]*sig.rho
,sigsq[4]*sig.rho,sigsq[5]*sig.rho,sigsq[6]*sig.rho,sigsq[7]*sig.rho,
sigsq[8]*sig.rho,sigsq[9]*sig.rho,sigsq[10]*sig.rho)

## Creating index for true and false null(>0 and <0 seperately)
m0=floor(m*pi0)#=0
m1=m-m0#!=0
m1.right=floor(m1*right.alloc)#>0
m1.left=m1-m1.right#<0

index0=1:m
index1=sample(index0,m0,replace=FALSE)#index set for \mu=0
index2=sample(setdiff(index0,index1),m1.left,replace=FALSE)#index set for \mu<0
index3=setdiff(setdiff(index0,index1),index2)#index set for \mu>0

## Feeding the null and non-null \mu-values(>0 and <0 seperately)
mu.left=runif(m1.left,u.l,u.u)
mu.right=runif(m1.right,-u.u,-u.l)
mu0=rep(0,m0)

mu=0
mu[index1]=mu0
mu[index2]=mu.left
mu[index3]=mu.right

## Generating the data matrix
normal.sample=mvrnorm(n,mu,mt.sigma)

## Getting output: Dataset, Index set of non-null \mu-values, Array of non-null \mu-values, Array of \mu-values
gen.data=normal.sample
index.nonnull=c(index2,index3)
mu.nonnull=c(mu.left,mu.right)
mu=c(mu0,mu.nonnull)
sig=mt.sigma[col(mt.sigma)==row(mt.sigma)]
return(list("gen.data"=gen.data,"index.nonnull"=index.nonnull,"mu.nonnull"=mu.nonnull,"mu"=mu,"sig"=sig))
}
####################################
#p=simul.details(0.5,10,0.5,1.5,0.5,0)
#data=p$gen.data
#n=10
#mu=p$mu
#sig=p$sig
#mu.nonnull=p$mu.nonnull
#index.nonnull=p$index.nonnull
####################################
##############################################################################################################################################
## Helping function:
expect.alt=function(samplesize,delvalue)
{
integrand=function(x)
{
pt(x,df=samplesize-1,ncp=delvalue)*dt(x,df=samplesize-1)
}
I1=integrate(integrand,lower=0,upper=Inf)$value
I2=integrate(integrand,lower=-Inf,upper=0)$value
return(2*(I1-I2))
}
## Main function returning e,\hat{e},\tilde{e} for a fixed dataset and supplied \pi_0-estimate:
e.details=function(n,data,pi0.init,mu,sig,mu.nonnull,index.nonnull)
# n is known from simul.data arg. 
# mu supplied from generated data function. sig is array of diagonal elements from mt.sigma matrix
# in the simulated data generation function. index.nonnull is also supplied. 
{
m=1000# As m is fixed in the numerical studies

## Calculating e:
e.true=0
for(i in 1:length(mu.nonnull))
{
del=mu.nonnull[i]/sqrt(sig[index.nonnull[i]])
e.true[i]=expect.alt(n,sqrt(n)*del)
}
e1=mean(e.true)

## Initialize some arrays for storing values under a for loop:
x.bar=0
s=0
delta.hat=0
for(i in 1:m)
{
  x.bar[i]=mean(data[,i])
  s[i]=sum((data[,i]-x.bar[i])^2)/(n-1)
  delta.hat[i]=sqrt(2/(n-1))*(factorial((n-1)/2-1)/factorial((n-2)/2-1))*(x.bar[i]/sqrt(s[i]))
}

## Calculating \hat{e}:
e.hat=0
for(i in 1:length(index.nonnull))
{
e.hat[i]=expect.alt(n,sqrt(n)*delta.hat[index.nonnull[i]])
}
e2=mean(e.hat)

## Calculating \tilde{e}:
m0.init=floor(m*pi0.init)
d=m-m0.init
e.hat=0
for(i in 1:m)
{
e.hat[i]=expect.alt(n,sqrt(n)*delta.hat[i])
}   
e=0
for(i in 1:d)
{
e[i]=sort(e.hat)[i]
}
e3=mean(e)


## Return:
true.e=e1
hat.e=e2
tilde.e=e3
return(list("true.e"=true.e,"hat.e"=hat.e,"tilde.e"=tilde.e)) 
}
############################
#p=e.details(n,data,0.5,mu,sig,mu.nonnull,index.nonnull)
#p$hat.e
#warnings()
############################
############################################################################################################

## Function for calculating e, \hat{e}, \tilde{e}(with bootstrap and Langaas as 
## initial estimators) with large repetition and fixed value of \pi_0:
e.pi0=function(pi0,n,u.l,u.u,right.alloc,rho)
{
m=1000
proc1=simul.details(pi0,n,u.l,u.u,right.alloc,rho)
data=proc1$gen.data
mu=proc1$mu
sig=proc1$sig
mu.nonnull=proc1$mu.nonnull
index.nonnull=proc1$index.nonnull

## Constructing array of p-values by m no. of single sample two-sided tests:
p.val=0
for(i in 1:m)
{
  test=t.test(data[,i],alternative="two.sided")
  p.val[i]=test$p.value
}
init.bootstrap=estim.pi0(p.val,pi0.method="st.boot")$pi0
init.langaas=estim.pi0(p.val,pi0.method="langaas")$pi0

proc2=e.details(n,data,init.bootstrap,mu,sig,mu.nonnull,index.nonnull)
proc3=e.details(n,data,init.langaas,mu,sig,mu.nonnull,index.nonnull)

## Estimates and Return:
true.e=proc2$true.e
hat.e=proc2$hat.e
tilde.e.boot=proc2$tilde.e
tilde.e.lang=proc3$tilde.e

pi0.est.boot.1=(mean(p.val)-tilde.e.boot)/(0.5-tilde.e.boot)
pi0.est.lang.1=(mean(p.val)-tilde.e.lang)/(0.5-tilde.e.lang)

tilde.e.boot.1=e.details(n,data,pi0.est.boot.1,mu,sig,mu.nonnull,index.nonnull)$tilde.e
tilde.e.lang.1=e.details(n,data,pi0.est.lang.1,mu,sig,mu.nonnull,index.nonnull)$tilde.e

pi0.est.boot.2=(mean(p.val)-tilde.e.boot.1)/(0.5-tilde.e.boot.1)
pi0.est.lang.2=(mean(p.val)-tilde.e.lang.1)/(0.5-tilde.e.lang.1)

return(list("true.e"=true.e,"hat.e"=hat.e,"tilde.e.boot"=tilde.e.boot,"tilde.e.lang"=tilde.e.lang,
"tilde.e.boot.1"=tilde.e.boot.1,"tilde.e.lang.1"=tilde.e.lang.1,"pi0.est.boot.1"=pi0.est.boot.1,
"pi0.est.lang.1"=pi0.est.lang.1,"pi0.est.boot.2"=pi0.est.boot.2,"pi0.est.lang.2"=pi0.est.lang.2,"data.fix"=data,
"p.val"=p.val))
}
####################
#e.pi0(0.5,10,0.5,1.5,0.5,0.2)
####################
############################################################################################################
W=function(x,cut)
{
sum(ifelse(x>cut,1,0))
}
cheng=function(n,data)
{
m=1000
p.val=0
x.bar=0
s=0
delta.hat=0
for(i in 1:m)
{
  test=t.test(data[,i],alternative="two.sided")
  p.val[i]=test$p.value
  x.bar[i]=mean(data[,i])
  s[i]=sum((data[,i]-x.bar[i])^2)/(n-1)
  delta.hat[i]=sqrt(2/(n-1))*(factorial((n-1)/2-1)/factorial((n-2)/2-1))*(x.bar[i]/sqrt(s[i]))
}
init.bootstrap=estim.pi0(p.val,pi0.method="st.boot")$pi0
m0.init=floor(m*init.bootstrap)
d=m-m0.init

lambda.index=seq(from=0.20,to=0.50,by=0.05)  ## Later need to take an index set viz. index.lambda.
pi0.cheng.lambda=0
Qhat=0
for(j in 1:length(lambda.index))
{
##calculating upper tail probability of 
Q.hat=0
for(i in 1:m)
{
Q.hat[i]=pt(qt(1-lambda.index[j]/2,df=n-1),df=n-1,ncp=sqrt(n)*delta.hat[i])- pt(qt(lambda.index[j]/2,df=n-1),df=n-1,ncp=sqrt(n)*delta.hat[i])
} 
Q=0
for(i in 1:d)
{
  Q[i]=sort(Q.hat)[i]
}  
Qhat=mean(Q)
pi0.cheng.lambda[j]=min(1,max(((W(p.val,lambda.index[j])-m*Qhat)/((1-lambda.index[j])-Qhat))/m,0))
}
pi0.cheng=mean(pi0.cheng.lambda)

return(pi0.cheng)
}

################
#cheng(n,data)
################
############################################################################################################
simul.e.pi0=function(sim.size,pi0,n,u.l,u.u,right.alloc,rho)
{
arr.true.e=0
arr.hat.e=0
arr.tilde.e.boot=0
arr.tilde.e.lang=0
arr.tilde.e.boot.1=0
arr.tilde.e.lang.1=0

arr.pi0.est.boot.1=0
arr.pi0.est.lang.1=0
arr.pi0.est.boot.2=0
arr.pi0.est.lang.2=0
arr.cheng=0


for(i in 1:sim.size)
{
proc=e.pi0(pi0,n,u.l,u.u,right.alloc,rho)
arr.true.e[i]=proc$true.e
arr.hat.e[i]=proc$hat.e
arr.tilde.e.boot[i]=proc$tilde.e.boot
arr.tilde.e.lang[i]=proc$tilde.e.lang
arr.tilde.e.boot.1[i]=proc$tilde.e.boot.1
arr.tilde.e.lang.1[i]=proc$tilde.e.lang.1

arr.pi0.est.boot.1[i]=proc$pi0.est.boot.1
arr.pi0.est.lang.1[i]=proc$pi0.est.lang.1
arr.pi0.est.boot.2[i]=proc$pi0.est.boot.2
arr.pi0.est.lang.2[i]=proc$pi0.est.lang.2

}
e.arr=c(mean(arr.true.e),mean(arr.hat.e),mean(arr.tilde.e.boot),mean(arr.tilde.e.lang),mean(arr.tilde.e.boot.1),mean(arr.tilde.e.lang.1))
mse.arr=c(mean((arr.pi0.est.boot.1-pi0)^2),mean((arr.pi0.est.lang.1-pi0)^2),mean((arr.pi0.est.boot.2-pi0)^2),mean((arr.pi0.est.lang.2-pi0)^2),mean((arr.cheng-pi0)^2))
bias.arr=c(mean(arr.pi0.est.boot.1-pi0),mean(arr.pi0.est.lang.1-pi0),mean(arr.pi0.est.boot.2-pi0),mean(arr.pi0.est.lang.2-pi0),mean(arr.cheng-pi0))

## Return

return(list("e.arr"=e.arr,"mse.arr"=mse.arr,"bias.arr"=bias.arr,
"arr.pi0.est.boot.1"=arr.pi0.est.boot.1,"arr.pi0.est.lang.1"=arr.pi0.est.lang.1,
"arr.pi0.est.boot.2"=arr.pi0.est.boot.2,"arr.pi0.est.lang.2"=arr.pi0.est.lang.2,"arr.cheng"=arr.cheng))
}
######################################################################################################################################################
