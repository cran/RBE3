ML.BE3 <-
function(data, tau=0.5, link.mu="logit", seeds.max=100, init.intercept=FALSE, seeds.valid=6)
{
dist.link <- switch(link.mu, logit="logis",
probit="norm", loglog="gumbel", cloglog="gumbel2")
g1.inv <- get(paste("q", dist.link, sep = ""), mode = "function")
      r1<-ncol(data$Z1);r2<-ncol(data$Z2);r3<-ncol(data$Z3)
if(init.intercept[1]==FALSE)
{
init.c<-rep(0,r1+r2+r3);desv<-5
}
if(init.intercept[1]!=FALSE)
{
init.c<-init.intercept
desv<-ifelse(init.c==0,5,1)
if(sum(init.c==0)==1)
{
ind=which(init.c==0)
seq.min=-10;seq.max=10;dif=1;init.aux=init.c
while(dif>=0.01)
{
bb=seq(seq.min,seq.max,by=dif)
ll=c()
for(j in 1:length(bb))
{init.aux[ind]=bb[j];ll[j]=-llike.BE3(init.aux, data, tau=tau, link.mu=link.mu)}
seq.min=bb[which.max(ll)]-dif;seq.max=bb[which.max(ll)]+dif
dif=dif/10
}
init.aux[ind]=bb[which.max(ll)]
init.c=init.aux
IF<-suppressWarnings(try(sqrt(diag(solve(hessian(llike.BE3, x0=init.c, data=data, tau=tau, link.mu=link.mu),tol=1e-10))),silent=TRUE))
if(!grepl("Error",IF)[1])
{
if(all(!is.na(IF)))
{
if(max(IF)<100 & max(IF)>0.01)
{
desv<-apply(rbind(IF,0.5),2,max)
}
}
}
if(!grepl("Error",IF)[1])
{
desv[ind]=ifelse(abs(bb[which.max(ll)])>5,5,1)
}
}
}
init<-init.c
IF.max=1;i<-1;flag=0; j<-1; llike.max=-llike.BE3(init, data, tau=tau, link.mu=link.mu)
auxi<-c();auxi$par=init.c;auxi$value=-llike.BE3(init.c, data, tau=tau, link.mu=link.mu)
while(i<=seeds.max & j<=seeds.valid)
{
method.s=ifelse(j<=ceiling(seeds.valid/2), "Nelder-Mead","BFGS")
#method.s="BFGS"
aux.ml<-try(optim(init,
fn=llike.BE3, data=data, 
tau=tau, link.mu=link.mu, method=method.s,
control=list(maxit=100000)),silent=TRUE)
if(!grepl("Error",aux.ml)[1])
{
if(aux.ml$convergence==0)
{
IF<-suppressWarnings(try(sqrt(diag(solve(hessian(llike.BE3, x0=aux.ml$par, data=data, tau=tau, link.mu=link.mu),tol=1e-10))),silent=TRUE))
if(!grepl("Error",IF)[1])
{
if(all(!is.na(IF)))
{
j<-j+1
if(max(IF)>0.2 & max(abs(aux.ml$par))<10)
{
if(llike.max < -aux.ml$value)
{llike.max<--aux.ml$value;auxi=aux.ml;IF.max=IF}
init.c<-aux.ml$par
desv2<-desv*0.5^(j-1)
desv<-apply(rbind(desv2,0.1),2,max)
}
}
}
}
}
cc<-paste(c("it",j-1,i,"cova",r1,r2,r3,"log-like",round(llike.max,2),round(min(desv),2),round(max(desv),2)),sep="")
#print(cc)
i<-i+1
init<-init.c+rnorm(length(init.c),sd=desv)
#init<-init.c+runif(r1+r2+r3,min=-desv/2,max=desv/2)
}
aux.ml=auxi;IF=IF.max
tt<-cbind(aux.ml$par, IF)
colnames(tt)<-c("Estimate","S.e")
rownames(tt)<-c(paste("beta1",1:ncol(data$Z1),sep=""),
paste("beta2",1:ncol(data$Z2),sep=""),
paste("beta3",1:ncol(data$Z3),sep=""))
list(estimate=tt, iter=i-1, loglike=-aux.ml$value)
}
