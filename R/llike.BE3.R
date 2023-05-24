llike.BE3 <-
function(theta, data, tau=0.5, link.mu="logit")
{
y<-data$y;Z1<-data$Z1;Z2<-data$Z2;Z3<-data$Z3
r1<-ncol(Z1);r2<-ncol(Z2);r3<-ncol(Z3)
beta1<-theta[1:r1]
beta2<-theta[1:r2+r1]
beta3<-theta[1:r3+r1+r2]

dist.link <- switch(link.mu, logit="logis",
probit="norm", loglog="gumbel", cloglog="gumbel2")
g1 <- get(paste("p", dist.link, sep = ""), mode = "function")

mu<-g1(Z1%*%beta1)
alpha<-exp(Z2%*%beta2)
beta<-exp(Z3%*%beta3)

q<-suppressWarnings(qbeta(tau, shape1=alpha, shape2=beta))

ll<-alpha*(log1p(q-1+1e-15)+log1p(-mu)-log1p(-q+1e-15)-log1p(mu-1+1e-15))-lbeta(alpha,beta)+(alpha-1)*log(y)+(beta-1)*log1p(-y)-(alpha+beta)*log1p(-(1-q*(1-mu)/((1-q)*mu))*y+1e-15)
min(-sum(ll),1e308)
#-sum(ll)
}
