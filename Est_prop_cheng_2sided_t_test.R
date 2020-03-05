## Code for computing Cheng's estimator for the proportion of true null hypotheses:
## Helping function:
W=function(x,cut)
{
sum(ifelse(x>cut,1,0))
}
## Main function
hat.pi0.U=function(data1,data2,pi0.init)#data1: Group 1 and data2: group 2. no. of rows 
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
lambda.index=seq(from=0.20,to=0.50,by=0.05)  ## Later need to take an index set viz. index.lambda.
m0.cheng.lambda=0
Qhat=0
for(j in 1:length(lambda.index))
{
##calculating upper tail probability of 
Q.hat=0
for(i in 1:m)
{
Q.hat[i]=pt(qt(1-lambda.index[j]/2,df=n1+n2-2),df=n1+n2-2,ncp=sqrt(n.star)*delta.hat[i])- pt(qt(lambda.index[j]/2,df=n1+n2-2),df=n1+n2-2,ncp=sqrt(n.star)*delta.hat[i])
} 
Q=0
for(i in 1:d)
{
  Q[i]=sort(Q.hat)[i]
} 
Qhat=mean(Q) 
m0.cheng.lambda[j]=min(1,max(((W(p.val,lambda.index[j])-m*Qhat)/((1-lambda.index[j])-Qhat))/m,0))
}
pi0.cheng=mean(m0.cheng.lambda)
# Return:
return(pi0.cheng)
}
