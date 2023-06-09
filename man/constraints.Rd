% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constraints.R
\name{constraints}
\alias{constraints}
\title{Create New Constraints}
\usage{
constraints(linear_part, offset, covariance = NULL, mean = NULL, rank = NULL)
}
\arguments{
\item{linear_part}{The linear part, $A$ of the affine constraint}

\item{offset}{The offset part, $b$ of the affine constraint}

\item{covariance}{Covariance matrix of Gaussian distribution to be truncated
(Defaults to NULL)}

\item{mean}{Mean vector of Gaussian distribution to be truncated
(Defaults to NULL)}

\item{rank}{If not NULL, this should specify the rank of the covariance matrix
(Defaults to NULL)}
}
\value{
An affine constraint object, which is a list of:
\item{linear_part}{The linear part, $A$ of the affine constraint}
\item{offset}{The offset part, $b$ of the affine constraint}
\item{dim}{The dimension of the covariance matrix}
\item{rank}{The rank of the covariance matrix}
\item{covariance}{Covariance matrix of Gaussian distribution to be truncated}
\item{mean}{Mean vector of Gaussian distribution to be truncated}
}
\description{
This is the core object for affine selection procedures.
It is meant to describe sets of the form $C$, where C = {z: Az <= b}.
Its main purpose is to consider slices through C and the conditional
distribution of a Gaussian N(mu,Sigma) restricted to such slices. In this
case, we have a truncated normal distribution condition on the pre-test.
In initialization, we create a new inequality.
}
\examples{
# Fit the optimization model first
Z = matrix(rnorm(10*3), nrow = 10, ncol = 3) +
    matrix(replicate(3, matrix(0, nrow = 10, ncol = 1)),
    nrow = 10)
errorTerm = MASS::mvrnorm(n=10, mu=rep(0,2),
            Sigma=rbind(c(1, 0.8),c(0.8, 1)))
D = Z \%*\% rep(1, 3) + errorTerm[1:10, 2]
Y = D * 1 + errorTerm[1:10, 1]
gl <- group_lasso_iv(Y, D, Z)
model <- fit_tsls(gl)

# Construct a box constraint
A_scaling = -diag(1)
b_scaling = rep(0, 1)
affine_con <- conditionalInference::constraints(A_scaling,
                                                b_scaling,
                                                model$cond_mean,
                                                model$cond_cov)

# Construct a white box constraint
cov_sol = conditionalInference::covariance_factors(affine_con)
white_A = affine_con$linear_part \%*\% cov_sol$sqrt_cov
den = sqrt(rowSums(white_A^2))
white_r = affine_con$linear_part \%*\% affine_con$mean
white_b = as.matrix(affine_con$offset,
                    dim=c(nrow(white_r),
                    ncol(white_r))) - white_r
conditionalInference::constraints(white_A / as.matrix(den,ncol=1),
                                  white_b / den)
}
\references{
Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
Plausibly Valid Instruments with Application to Mendelian Randomization.

Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
Effects After Testing Instrument Strength in Linear Models.
}
\seealso{
\link{sample_from_constraints} for how to generate random samples from
this type of box constraint.
}
