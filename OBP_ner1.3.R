obp_ner1.3<-function(y,x,area,ni,Ni,Xibar,xibar)
{
require(MASS)
ar<-sort(unique(area))
p<-ncol(x)
m<-length(ar)
fi=ni/Ni

yi<- function(i){ matrix(c(y[area==ar[i]]),ni[i],1)}
xi<- function(i){x[area==ar[i],]}
yibar<- rep(0,m)
yibar<-as.numeric(tapply(y,area,mean))
#parameter estimation 

mu.i.2.hat<-rep(0,m)
for(i in 1:m)mu.i.2.hat[i]=mean(yi(i)^2) - ((Ni[i]-1)/(Ni[i]*(ni[i]-1)))* sum((yi(i)-yibar[i])^2)

Q.gamma=function(gamma){ 
  tmp.Gammai = ((ni*gamma)/(1+ni*gamma)) 
  tmp.G=diag((fi+(1-fi)*tmp.Gammai))
  tmp.G.2=diag((fi+(1-fi)*tmp.Gammai)^2)
  
  tmp.H= t(Xibar-tmp.G%*%xibar)%*%(Xibar-tmp.G%*%xibar)
  tmp.L= (diag(1,m)-2*tmp.G)%*%Xibar +(tmp.G.2%*%xibar)

  Q<-t(yibar)%*%(tmp.G.2-tmp.L%*%ginv(tmp.H)%*%t(tmp.L))%*%yibar+t(rep(1,m))%*%(diag(1,m)-2*tmp.G)%*%mu.i.2.hat
  Q
}
fgamma=optimise(Q.gamma,interval=c(0.001,1000))
est.gamma=fgamma$minimum
Gammai = ((ni*est.gamma)/(1+ni*est.gamma)) 
G=diag((fi+(1-fi)*Gammai))
G.2=diag((fi+(1-fi)*Gammai)^2)
H= t(Xibar-G%*%xibar)%*%(Xibar-G%*%xibar)
L=(G.2%*%xibar)+(diag(1,m)-2*G)%*%Xibar
beta=c(ginv(H)%*%t(L)%*%yibar)
  
mu.i<-Xibar%*%beta+(fi+(1-fi)*(ni*est.gamma/(1+ni*est.gamma)))*(yibar-xibar%*%beta)
list(OBP=mu.i,beta=beta,est.gamma=est.gamma)
}