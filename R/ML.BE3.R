ML.BE3=function (data, tau = 0.5, link.mu = "logit")
{
    dist.link <- switch(link.mu, logit = "logis", probit = "norm",
        loglog = "gumbel", cloglog = "gumbel2")
    g1.inv <- get(paste("q", dist.link, sep = ""), mode = "function")
    r1 <- ncol(data$Z1)
    r2 <- ncol(data$Z2)
    r3 <- ncol(data$Z3)
    if(r1>1) aa=coef(lm(qlogis(y)~data$Z1[,-1,drop=FALSE]))
    if(r1==1) aa=coef(lm(qlogis(y)~1))
    init <- c(aa,rep(0, r2 + r3))
    llike.max = -llike.BE3(init, data, tau = tau, link.mu = link.mu)
    aux.ml <- suppressWarnings(try(optim(init, fn = llike.BE3, data = data,
         tau = tau, link.mu = link.mu, method = "BFGS",
         control = list(maxit = 1e+05),hessian=TRUE), silent = TRUE))
    IF <- suppressWarnings(try(sqrt(diag(solve(aux.ml$hessian))), silent = TRUE))
    llike.max <- -aux.ml$value
    tt <- cbind(aux.ml$par, IF)
    colnames(tt) <- c("Estimate", "S.e")
    rownames(tt) <- c(paste("beta1", 1:ncol(data$Z1), sep = ""),
        paste("beta2", 1:ncol(data$Z2), sep = ""), paste("beta3",
            1:ncol(data$Z3), sep = ""))
    list(estimate = tt, logLik = -aux.ml$value)
}
