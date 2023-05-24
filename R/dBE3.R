dBE3 <-
function(x, mu = 0.5, alpha = 1, beta = 1, tau=0.5, log=FALSE)
{
 if (is.null(x)) 
        stop("x must be specified")
    if (any(is.null(c(mu, alpha, beta, tau)))) 
        stop("mu, alpha, beta and tau must be specified")
    if (any(tau<=0) | any(tau >= 1)) 
        stop("tau must be between 0 and 1")
    if (any(mu <= 0)) 
        stop("mu must be positive")
    if (any(alpha <= 0)) 
        stop("alpha must be positive")
    if (any(beta <= 0)) 
        stop("beta must be positive")
z=qbeta(tau, shape1=alpha,shape2=beta)
lambda<-(1-mu)*z/(mu*(1-z))
log.f<-alpha*log(lambda)-(alpha+beta)*log1p(-(1-lambda)*x)+dbeta(x, shape1=alpha, shape2=beta, log=TRUE)
if(!log) log.f=exp(log.f)
log.f
}
