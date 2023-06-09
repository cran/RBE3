rBE3 <-
function(n, mu = 0.5, alpha = 1, beta = 1, tau=0.5)
{
 if (is.null(n)) 
        stop("sample size must be specified")
    if (any(is.null(c(mu, alpha, beta, tau)))) 
        stop("mu, alpha, beta and tau must be specified")
    if (round(n)!=n | n <= 0) 
        stop("sample size must be a positive integer")
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
X1<-rgamma(n, shape=alpha, rate=lambda)
X2<-rgamma(n, shape=beta, rate=1)
X1/(X1+X2)
}
