## Function for computing the sum based bias corrected estimator (two-sided t-tests)
## Attach libraries:
library(cp4p)
library(MASS)
library(Matrix)
## Helping function:
expect.alt=function(samplesize1,samplesize2,delvalue)
{
integrand=function(x)
{
pt(x,df=samplesize1+samplesize2-2,ncp=delvalue)*dt(x,df=samplesize1+samplesize2-2)
}
I1=integrate(integrand,lower=0,upper=Inf)$value
I2=integrate(integrand,lower=-Inf,upper=0)$value
return(2*(I1-I2))
}
## Main function:
hat.pi0.E=function(data1,data2,pi0.init)#data1: Group 1 and data2: group 2. no. of rows 
# are same but, no. of columns may be different. row: individual, col: gene.
# pi0.init: numerical value in (0,1)- initial estimator of the parameter \pi_0, usually 
# bootstrap estimator of the same. 
{
m=dim(data1)[1]
n1=dim(data1)[2]
n2=dim(data2)[2]
n.star=(n1*n2)/(n1+n2)
p.val=0
x1.bar=0
x2.bar=0
s1=0
s2=0
s=0
delta.hat=0
for(i in 1:m)
{
  test=t.test(data1[i,],data2[i,],alternative="two.sided")
  p.val[i]=test$p.value
  x1.bar[i]=mean(data1[i,])
  x2.bar[i]=mean(data2[i,])
  s1[i]=sum((data1[i,]-x1.bar[i])^2)/(n1-1)
  s2[i]=sum((data2[i,]-x2.bar[i])^2)/(n2-1)
  s[i]=((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  delta.hat[i]=sqrt(2/(n1+n2-2))*(factorial((n1+n2)/2-1)/factorial((n1+n2)/2-1.5))*((x1.bar[i]-x2.bar[i])/sqrt(s[i]))
}
m0.init=floor(pi0.init*m)
d=m-m0.init
e.hat=0
for(i in 1:m)
{
e.hat[i]=expect.alt(n1,n2,sqrt(n.star)*delta.hat[i])
}
e=0
for(i in 1:d)
{
e[i]=sort(e.hat)[i]
}
ehat=mean(e)
pi0.new=min(1,max((mean(p.val)-ehat)/(0.5-ehat),0))
## Return:
return(pi0.new)
}