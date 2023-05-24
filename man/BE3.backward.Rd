\name{BE3.backward}
\alias{BE3.backward.crit}
%\alias{BE3.backward.sign}
\title{
backward stepwise regression for RBE3 model based on the AIC criterion or significance.
}
\description{
BE3.backward.crit implements the covariates selection based on backward and the Akaike's information criteria (AIC).
%BE3.backward.sign implements the covariates selection based on backward and significance of the covariates.
}
\usage{
BE3.backward.crit(data, tau = 0.5, link.mu = "logit", seeds.max = 100)
%BE3.backward.sign(data, tau = 0.5, link.mu = "logit", seeds.max = 100,
%      significance = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list containing the response vector (\eqn{y}), and the matrices to model \eqn{\mu}, the \eqn{\tau}-quantile of distribution,
and the shape parameters \eqn{\alpha} and \eqn{\beta}, labeled as \eqn{Z_1, Z_2} and \eqn{Z_3}, respectively.
}
  \item{tau}{
the quantile of the distribution to be modelled (\eqn{0<\tau<1}).
}
  \item{link.mu}{
link function to be used for \eqn{\mu}: logit (default), probit, loglog or cloglog. 
}
  \item{seeds.max}{
maximum of seeds to be reached. 
}
%\item{significance}{
%level of significance used. Only for BE3.backward.sign function.
%}
}
\value{
A list containg the covariates to be included for modelling \eqn{\mu}, \eqn{\alpha} and \eqn{\beta}, respectively.
}
\author{
Diego Gallardo and Marcelo Bourguignon.
}
\examples{
\donttest{
##Simulating two covariates
set.seed(2100)
x1<-rnorm(200); x2<-rbinom(200, size=1, prob=0.5)
##Desing matrices: Z1 includes x1 and x2, 
##Z2 includes only x1 and Z3 includes only x2
Z1=model.matrix(~x1+x2);Z2=model.matrix(~x1);Z3=model.matrix(~x2)
##Fixing parameters
theta=c(1, 0.2, -0.5); nu=c(0.5,-0.2); eta=c(-0.5, 0.3); tau=0.4
mu=plogis(Z1\%*\%theta); alpha=exp(Z2\%*\%nu); beta=exp(Z3\%*\%eta)
y=rBE3(200, mu, alpha, beta, tau=tau)
data=list(y=y, Z1=Z1, Z2=Z2, Z3=Z3)
BE3.backward.crit(data, tau = tau)
}
}

