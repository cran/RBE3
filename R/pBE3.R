pBE3 <-
function(q, mu = 0.5, alpha = 1, beta = 1, tau=0.5,lower.tail = TRUE, log.p = FALSE)
{
if (is.null(q)) 
        stop("q must be specified")
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
lf=pbeta(lambda*q/(1+lambda*q-q),shape1=alpha, shape2=beta, log.p=TRUE)
    if (!log.p) 
        lf = exp(lf)
 if (!lower.tail) 
        lf = 1 - lf
    lf
}
