#' Placeholder for development of autocorrelation function in neuroimaging
#' This will be necessary to compute corrected general linear models that are comparable to current implementations.
#' @reference
#' Kiebel et al (2003) A heuristic for the degrees of freedom of statistics based on multiple variance parameters
#'
#' rft.ACF.r
rft.ACF <-function(){

# 1. estimate beta using ordinary least-squares (OLS)
# 2. Compute F ratio

# y=X*beta+error
# In= the identity matrix of rank 'n'
# Xo denotes the subspace of the full model X
# M is a projector onto the subspace X that is not spanned by Xo
## M can be defined by direct specificationof Xo or through contrast weights 'c'

R=In-X*Xinv
Ro=In-XoXoinv
M=Ro-R

# classical statistics formed from log-likelihood ratios
p(y|X) = ((2*pi*sd)^(-n/2)) * exp(-(1/(2*sd))*yT*RT*Ry)
p(y|Xo) = ((2*pi*sd)^(-n/2)) * exp(-(1/(2*sd))*yT*RoT*Roy)
l = log((p(y|X))/(p(y|Xo)))
= (1/(2*sd))*(yT*RoT*Roy-yT*RT*Ry)
=(1/(2*sd))*(yT*My)

# scaling the log-likelihood ratio (l) so its expectation under the null hypothesis is unity
l/lo = (yT*My)/(sd*tr(M))

# Replacing the unknown error variance variance parameter with this estimator gives the classical F ratio
f = (yT*My)/(sd*tr(M))=((yT*My)/(tr(M)))/(yT*Ry/tr(R))

# ReML estimator of standard deviation
# 'tr' = trace operator

Rml = In - X*((XT*(C^-1)*X)^-1)*XT*(C^-1)
sd=(yT*RTml*(Q^-1)*Rml*y)/tr(Rml)

# nonspherical variance component
## for general linear model with single variance parameter and a nonspherical variance component Q F statistic and df can be estimated
f = ((yT*My)/(tr(MQ))/((yT*Ry)/tr(RQ))
vo = (tr(MQ)^2)/tr(MQMQ)
v = (tr(RQ)^2)/tr(RQRQ)
# variance parameter estimate
sd = (yT*Ry)/tr(RQ)
# <> denotes expectation and Var() the variance of a random variable (<sd>=expected sd)
Var(<sd>) = (2*sd*tr(RQRQ))/(tr(RQ)^2)

# Normal approximation
## Taylor series expansion for log-likelihood log (p(y|lambdaEst))

p(y|lambda) = N(lambdaEst, I(lambdaEst)^-1) 
# I(lambdaEst) is the observed information matrix with:
I(lambdaEst) = - (partial derivative lambda)*log(p(y|lambdaEst))


}
