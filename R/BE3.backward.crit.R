BE3.backward.crit <-
function(data, tau=0.5, link.mu="logit", seeds.max=100) 
{
cov.allcombinations<-function(r1,r2,r3)
{
cov.allcombinations.aux<-function(r)
{
 write.ind<-function(x){ss<-length(x);sum(sort(x,decreasing=TRUE)*10^(0:(ss-1)))}
 comb<-1
 if(r>1)
 {
  for(a in 1:(r-1))
  {
  aux<-apply(combinations(r-1,a)+1,1,write.ind)
  comb<-c(comb,aux+10^a)
  }
 }
 comb
}
 cov1<-cov.allcombinations.aux(r1)
 cov2<-cov.allcombinations.aux(r2)
 cov3<-cov.allcombinations.aux(r3)
 cov<-c()
 for(k in 1:length(cov3))
 {
 for(j in 1:length(cov2))
 {
 for(i in 1:length(cov1))
 {
cov<-rbind(cov,c(cov1[i],cov2[j],cov3[k]))
 }
 }
 }
 cov
}
ext.ind<-function(x)
{
vec<-c()
ss=floor(log(x,base=10))+1
for(i in 1:(ss))
{
vec[i]<-(x-10^i*floor(x/10^i)-sum(vec[1:(i-1)]))
}
sort(vec/10^(1:ss-1))
}
  y<-data$y;Z1<-data$Z1;Z2<-data$Z2;Z3<-data$Z3
  Z<-cbind(Z1,Z2,Z3)
  r1<-ncol(Z1);r2<-ncol(Z2);r3<-ncol(Z3)
  aux<-cov.allcombinations(r1,r2,r3)
  resumen<-c()
  for(i in 1:nrow(aux))
  {
     ind1<-ext.ind(aux[i,1])
     ind2<-ext.ind(aux[i,2])
     ind3<-ext.ind(aux[i,3])
     data.t=list(y=y, Z1=Z1[,ind1,drop=FALSE],
         Z2=Z2[,ind2,drop=FALSE],
         Z3=Z3[,ind3,drop=FALSE])
      fit<-ML.BE3(data.t, tau=tau, link.mu=link.mu, seeds.max=seeds.max)
      resumen<-rbind(resumen, c(aux[i,],fit$loglike, length(c(ind1,ind2,ind3))))
  }
  AIC<--2*resumen[,4]+2*resumen[,5]
  BIC<--2*resumen[,4]+log(length(y))*resumen[,5]
  resumen<-cbind(resumen,AIC,BIC)
  colnames(resumen)<-c("cov1","cov2","cov3","loglike","parameters","AIC","BIC")
  resumen
}
