\name{BE3}
\alias{dBE3}
\alias{pBE3}
\alias{qBE3}
\alias{rBE3}
\title{
The Generalized Beta Distribution
}
\description{
Density, distribution function, quantile function and random generation for the generalized beta distribution.
}
\usage{
dBE3(x, mu = 0.5, alpha = 1, beta = 1, tau = 0.5, log = FALSE)
pBE3(q, mu = 0.5, alpha = 1, beta = 1, tau = 0.5, lower.tail = TRUE, log.p = FALSE)
qBE3(p, mu = 0.5, alpha = 1, beta = 1, tau = 0.5, log.p = FALSE)
rBE3(n, mu = 0.5, alpha = 1, beta = 1, tau = 0.5)
}
\arguments{
  \item{x, q}{
vector of quantiles.
}
  \item{p}{
vector of probabilities.
}
  \item{n}{
number of observations.
}
  \item{mu}{
vector of \eqn{\tau}-quantiles of the distribution.
}
  \item{alpha, beta}{
shape parameters of the distribution
}
  \item{tau}{
corresponding quantile of the distribution (\eqn{0<\tau<1})}
  \item{log, log.p}{
logical; if TRUE, probabilities \eqn{p} are given as \eqn{\log p}.
}
  \item{lower.tail}{
logical; if TRUE (default), probabilities are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
}
}
\details{
The probability density function for the generalized beta distribution is
\deqn{
f(y;\lambda,\alpha,\beta)=\frac{\lambda^\alpha y^{\alpha-1}(1-y)^{\beta-1}}{B(\alpha, \beta)[1-(1-\lambda)y]^{\alpha+\beta}}, \quad 0<y<1, 
}
where \eqn{\alpha, \beta>0} and \eqn{\lambda>0}.
We consider the reparameterization in terms of the \eqn{\tau}-quantile of the distribution, say \eqn{0<\mu<1}, taking 
\deqn{
\lambda=\frac{(1-\mu)}{\mu}\frac{z_{\alpha,\beta}(\tau)}{[1-z_{\alpha,\beta}(\tau)]},
}
with \eqn{z_{\alpha,\beta}(\tau)} denoting the \eqn{\tau}-quantile of the usual beta distribution with shape parameters \eqn{\alpha} and \eqn{\beta}. 
The cumulative distribution function is given by
\deqn{
F(y;\lambda,\alpha,\beta)=I_{\lambda x/(1+\lambda x -x)}(\alpha, \beta), \quad 0<y<1,
}
where \eqn{I_x(\alpha,\beta)=B_x(\alpha,\beta)/B(\alpha,\beta)} is the incomplete beta funcion ratio, \eqn{B_x(\alpha,\beta)=\int_0^x w^{\alpha-1}(1-w)^{\beta-1}dw} 
is the incomplete beta function and \eqn{B(\alpha,\beta)=\Gamma(\alpha)\Gamma(\beta)/\Gamma(\alpha+\beta)} is the ordinary beta function.
The quantile of the distribution can be represented as
\deqn{
q(\tau;\lambda,\alpha,\beta)=\frac{z_{\alpha,\beta}(\tau)}{\lambda[1-z_{\alpha,\beta}(\tau)]+z_{\alpha,\beta}(\tau)}, \quad 0<\tau<1.
}
Random generation can be performed using the stochastic representation of the model. If \eqn{X_1 \sim \mbox{Gamma}(\alpha,\theta_1)} and 
\eqn{X_2 \sim \mbox{Gamma}(\beta,\theta_2)}, then
\deqn{
\frac{X_1}{X_1+X_2}\sim GB3(\alpha,\beta,\lambda),
}
where \eqn{\lambda=\theta_1/\theta_2.}
}
\value{
dBE3 gives the density, pBE3 gives the distribution function, qBE3 gives the quantile function, and rBE3 generates random deviates.

The length of the result is determined by n for rBE3, and is the maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.

}
\references{
Libby, D. L. and Novick, M. R. (1982). Multivariate generalized beta-distributions with applications
to utility assessment. Journal of Educational Statistics, 7.
}
\author{
Diego Gallardo and Marcelo Bourguignon
}
\examples{
rBE3(20, mu=0.5, alpha=2, beta=1)
dBE3(c(0.4,0.7), mu=0.5, alpha=2, beta=1)
pBE3(c(0.4,0.7), mu=0.5, alpha=2, beta=1)
}
