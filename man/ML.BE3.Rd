\name{ML.BE3}
\alias{ML.BE3}
\title{
Perform the parameter estimation for the Generalized beta distribution
}
\description{
ML.BE3 computes the maximum likelihood estimates based on the maximum likelihood method. 
}
\usage{
ML.BE3(data, tau = 0.5, link.mu = "logit", seeds.max = 100, init.intercept = FALSE, 
       seeds.valid = 6)
}
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
  \item{init.intercept}{
initial values for the intercept terms. If no values have been provided, a vector of zeros will be assumed.
}
  \item{seeds.valid}{
maximum of estimation procedures with convergence to be reached. 
}
}
\details{
Covariates are included as \eqn{g_1(\mu_i(\tau))=\mathbf{Z}_{1i}^\top {\bm \theta}(\tau)}, \eqn{g_2(\alpha_i(\tau))=\mathbf{Z}_{2i}^\top {\bm \nu}(\tau)}
and \eqn{g_3(\beta_i(\tau))=\mathbf{Z}_{3i}^\top {\bm \eta}(\tau)}, where \eqn{{\bm \theta}(\tau)=(\theta_1(\tau),\ldots,\theta_{r_1}(\tau))},
\eqn{{\bm \nu}(\tau)=(\nu_1(\tau),\ldots,\nu_{r_2}(\tau))} and \eqn{{\bm \eta}(\tau)=(\eta_1(\tau),\ldots,\eta_{r_3}(\tau))}, where \eqn{r_1, r_2} and \eqn{r_3}
are the dimensions of \eqn{Z_1, Z_2} and \eqn{Z_3}, respectively.
Since in previous attempts to maximize the log-likelihood function of the BE3 model we have noticed that the result can strongly depend on the seed used, we have decided to implement an algorithm that uses the Nelder-Mead option and the BFGS option of the optimal function, each of them used half the times defined in seed.valid (6 by default). 
}
\value{
a list containing the following elements
\item{estimate}{A matrix with the estimates }

}
\references{
Bourguignon, M., Gallardo, D.I., Saulo, H. (2023) A parametric quantile beta regression for modeling case fatality rates of COVID-19. Submitted.
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
ML.BE3(data, tau=tau)
}
}


