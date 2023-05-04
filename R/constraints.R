#' Create New Constraints
#'
#' This is the core object for affine selection procedures.
#' It is meant to describe sets of the form $C$, where C = {z: Az <= b}.
#' Its main purpose is to consider slices through C and the conditional
#' distribution of a Gaussian N(mu,Sigma) restricted to such slices. In this
#' case, we have a truncated normal distribution condition on the pre-test.
#' In initialization, we create a new inequality.
#'
#' @param linear_part The linear part, $A$ of the affine constraint
#' @param offset The offset part, $b$ of the affine constraint
#' @param covariance Covariance matrix of Gaussian distribution to be truncated
#'                   (Defaults to NULL)
#' @param mean Mean vector of Gaussian distribution to be truncated
#'             (Defaults to NULL)
#' @param rank If not NULL, this should specify the rank of the covariance matrix
#'             (Defaults to NULL)
#'
#' @return An affine constraint object, which is a list of:
#' \item{linear_part}{The linear part, $A$ of the affine constraint}
#' \item{offset}{The offset part, $b$ of the affine constraint}
#' \item{dim}{The dimension of the covariance matrix}
#' \item{rank}{The rank of the covariance matrix}
#' \item{covariance}{Covariance matrix of Gaussian distribution to be truncated}
#' \item{mean}{Mean vector of Gaussian distribution to be truncated}
#' @export
#'
#' @references
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
#' Plausibly Valid Instruments with Application to Mendelian Randomization.
#'
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
#' Effects After Testing Instrument Strength in Linear Models.
#'
#' @seealso
#' \link{sample_from_constraints} for how to generate random samples from
#' this type of box constraint.
#'
#' @examples
#' # Fit the optimization model first
#' Z = matrix(rnorm(10*3), nrow = 10, ncol = 3) +
#'     matrix(replicate(3, matrix(0, nrow = 10, ncol = 1)),
#'     nrow = 10)
#' errorTerm = MASS::mvrnorm(n=10, mu=rep(0,2),
#'             Sigma=rbind(c(1, 0.8),c(0.8, 1)))
#' D = Z %*% rep(1, 3) + errorTerm[1:10, 2]
#' Y = D * 1 + errorTerm[1:10, 1]
#' gl <- group_lasso_iv(Y, D, Z)
#' model <- fit_tsls(gl)
#'
#' # Construct a box constraint
#' A_scaling = -diag(1)
#' b_scaling = rep(0, 1)
#' affine_con <- conditionalInference::constraints(A_scaling,
#'                                                 b_scaling,
#'                                                 model$cond_mean,
#'                                                 model$cond_cov)
#'
#' # Construct a white box constraint
#' cov_sol = conditionalInference::covariance_factors(affine_con)
#' white_A = affine_con$linear_part %*% cov_sol$sqrt_cov
#' den = sqrt(rowSums(white_A^2))
#' white_r = affine_con$linear_part %*% affine_con$mean
#' white_b = as.matrix(affine_con$offset,
#'                     dim=c(nrow(white_r),
#'                     ncol(white_r))) - white_r
#' conditionalInference::constraints(white_A / as.matrix(den,ncol=1),
#'                                   white_b / den)
constraints <- function(linear_part,
                        offset,
                        covariance=NULL,
                        mean=NULL,
                        rank=NULL) {

  linear_part = linear_part
  offset = array(offset)

  if ((length(dim(linear_part))) == 2) {
    dim = dim(linear_part)[2]
  } else {
    dim = dim(linear_part)[1]
  }

  if (is.null(rank)) {
    rank = dim
  }

  if (is.null(covariance)) {
    covariance = diag(dim)
  }

  if (is.null(mean)) {
    mean = rep(0, dim)
  }

  return_list <- list("linear_part" = linear_part,
                      "offset" = offset,
                      "dim" = dim,
                      "rank" = rank,
                      "covariance" = covariance,
                      "mean" = mean)
  return(return_list)
}


#' Whiten
#'
#' Return a whitened version of constraints in a different basis, and
#' a change of basis matrix. If con$covariance is rank deficient, the change-of
#' basis matrix will not be square.
#'
#' @param con The box constraint
#'
#' @return A whiten object, which is a list of:
#' \item{sqrt_cov}{Squared root of the covariance matrix}
#' \item{sqrt_inv}{Squared root of the inverse covariance matrix}
#' \item{rowspace}{Row space matrix}
#' \item{mu}{Mean vector of Gaussian distribution to be truncated}
#' \item{linear_part}{The linear part, $A$ of the affine constraint}
#' \item{offset}{The offset part, $b$ of the affine constraint}
#' \item{dim}{The dimension of the covariance matrix}
#' \item{rank}{The rank of the covariance matrix}
#' \item{covariance}{Covariance matrix of Gaussian distribution to be truncated}
#' \item{mean}{Mean vector of Gaussian distribution to be truncated}
#' @export
#'
#' @references
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
#' Plausibly Valid Instruments with Application to Mendelian Randomization.
#'
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
#' Effects After Testing Instrument Strength in Linear Models.
#'
#' @seealso
#' \link{covariance_factors} that finds covariance factor, which is a possibly
#' non-square square-root.
#'
#' \link{constraints} for constructing a box constraint for the optimization
#' problem.
#'
#' @examples
#' # Fit the group lasso optimization model
#' Z = matrix(rnorm(10*3), nrow = 10, ncol = 3) +
#'     matrix(replicate(3, matrix(0, nrow = 10, ncol = 1)),
#'     nrow = 10)
#' errorTerm = MASS::mvrnorm(n=10, mu=rep(0,2),
#'             Sigma=rbind(c(1, 0.8),c(0.8, 1)))
#' D = Z %*% rep(1, 3) + errorTerm[1:10, 2]
#' Y = D * 1 + errorTerm[1:10, 1]
#' gl <- group_lasso_iv(Y, D, Z)
#' model <- fit_tsls(gl)
#'
#' # Whiten the box constraint
#' whiten(model$sampler$affine_con)
whiten <- function(con) {

  cov_sol = covariance_factors(con)
  sqrt_cov = cov_sol$sqrt_cov
  sqrt_inv = cov_sol$sqrt_inv
  rowspace = cov_sol$rowspace

  new_A = con$linear_part %*% sqrt_cov
  den = sqrt(rowSums(new_A^2))
  new_r = con$linear_part %*% con$mean
  new_b = as.matrix(con$offset, dim=c(nrow(new_r),ncol(new_r))) - new_r
  new_con = conditionalInference::constraints(new_A / as.matrix(den,ncol=1),
                                              new_b / den)
  mu = con$mean

  return_list <- list("sqrt_cov" = sqrt_cov,
                      "sqrt_inv" = sqrt_inv,
                      "rowspace" = rowspace,
                      "mu" = mu,
                      "linear_part" = new_con$linear_part,
                      "offset" = new_con$offset,
                      "dim" = new_con$dim,
                      "rank" = new_con$rank,
                      "covariance" = new_con$covariance,
                      "mean" = new_con$mean)
  return(return_list)
}


#' Inverse Map
#'
#' A callable inverse_map.
#'
#' @param Z The p instruments
#' @param white_con Whiten constraints
#'
#' @return The inverse map of white_con.
#' @export
#'
#' @references
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
#' Plausibly Valid Instruments with Application to Mendelian Randomization.
#'
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
#' Effects After Testing Instrument Strength in Linear Models.
#'
#' @examples
#' # Fit the group lasso optimization model
#' Z = matrix(rnorm(10*3), nrow = 10, ncol = 3) +
#'     matrix(replicate(3, matrix(0, nrow = 10, ncol = 1)),
#'     nrow = 10)
#' errorTerm = MASS::mvrnorm(n=10, mu=rep(0,2),
#'             Sigma=rbind(c(1, 0.8),c(0.8, 1)))
#' D = Z %*% rep(1, 3) + errorTerm[1:10, 2]
#' Y = D * 1 + errorTerm[1:10, 1]
#' gl <- group_lasso_iv(Y, D, Z)
#' model <- fit_tsls(gl)
#'
#' # Whiten the box constraint
#' white_con = whiten(model$sampler$affine_con)
#' white_Y = forward_map(model$sampler$initial_point, white_con)
#' Y = as.matrix(model$sampler$initial_point)
#' direction_of_interest = rnorm(length(Y))
#' white_direction_of_interest =
#' forward_map(model$sampler$affine_con$covariance %*% direction_of_interest,
#' white_con)
#'
#' # Generate white samples from truncated normal distribution
#' white_samples = sample_truncnorm_white(white_con$linear_part,
#' white_con$offset, white_Y, white_direction_of_interest,
#' how_often=2000, ndraw=1000, burnin=1000, sigma=1,
#' use_constraint_directions=TRUE, use_random_directions=TRUE)
#'
#' # Construct inverse map of white samples
#' inverse_map(white_samples, white_con)
inverse_map <- function(Z, white_con) {

  sqrt_cov = white_con$sqrt_cov
  mu = white_con$mu

  if (length(dim(Z)) == 2) {
    mu = as.matrix(mu, ncol=1)
    result = apply(t(sqrt_cov %*% t(Z)), c(1,2), function(x) x + mu)
    return(result)
  } else {
    result = apply(sqrt_cov %*% Z, c(1,2), function(x) x + mu)
    return(result)
  }
}


#' Forward Map
#'
#' A callable forward_map.
#'
#' @param W The input matrix
#' @param white_con Whiten constraints
#'
#' @return The forward map of white_con.
#' @export
#'
#' @references
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
#' Plausibly Valid Instruments with Application to Mendelian Randomization.
#'
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
#' Effects After Testing Instrument Strength in Linear Models.
#'
#' @examples
#' # Fit the group lasso optimization model
#' Z = matrix(rnorm(10*3), nrow = 10, ncol = 3) +
#'     matrix(replicate(3, matrix(0, nrow = 10, ncol = 1)),
#'     nrow = 10)
#' errorTerm = MASS::mvrnorm(n=10, mu=rep(0,2),
#'             Sigma=rbind(c(1, 0.8),c(0.8, 1)))
#' D = Z %*% rep(1, 3) + errorTerm[1:10, 2]
#' Y = D * 1 + errorTerm[1:10, 1]
#' gl <- group_lasso_iv(Y, D, Z)
#' model <- fit_tsls(gl)
#'
#' # Whiten constraint and create forward_map
#' white_con = whiten(model$sampler$affine_con)
#' forward_map(model$sampler$initial_point, white_con)
forward_map <- function(W, white_con) {

  sqrt_inv = white_con$sqrt_inv
  mu = white_con$mu

  result = sqrt_inv %*% as.matrix(W - mu)
  return(result)
}


#' Covariance Factors
#'
#' Finding a possibly non-square square-root.
#'
#' @param con Constraints
#' @param force If True, force a recomputation of the covariance.
#'              If not, assumes that covariance has not changed.
#'
#' @return A covariance_factors object, which is a list of
#' \item{sqrt_cov}{Squared root of the covariance matrix}
#' \item{sqrt_inv}{Squared root of the inverse covariance matrix}
#' \item{rowspace}{Row space matrix}
#' @export
#'
#' @references
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
#' Plausibly Valid Instruments with Application to Mendelian Randomization.
#'
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
#' Effects After Testing Instrument Strength in Linear Models.
#'
#' @examples
#' # Fit the group lasso optimization model
#' Z = matrix(rnorm(10*3), nrow = 10, ncol = 3) +
#'     matrix(replicate(3, matrix(0, nrow = 10, ncol = 1)),
#'     nrow = 10)
#' errorTerm = MASS::mvrnorm(n=10, mu=rep(0,2),
#'             Sigma=rbind(c(1, 0.8),c(0.8, 1)))
#' D = Z %*% rep(1, 3) + errorTerm[1:10, 2]
#' Y = D * 1 + errorTerm[1:10, 1]
#' gl <- group_lasso_iv(Y, D, Z)
#' model <- fit_tsls(gl)
#'
#' # Covariance factors for the original box constraint
#' covariance_factors(model$sampler$affine_con)
#' # still compute sqrt_cov, sqrt_inv, rowspace as the constraint
#' # doesn't contain such attributes
#' covariance_factors(model$sampler$affine_con, force=FALSE)
#'
#' # Whiten the box constraint
#' white_con = whiten(model$sampler$affine_con)
#'
#' # Covariance factors for the white box constraint
#' covariance_factors(white_con)
#' covariance_factors(white_con, force=FALSE)
covariance_factors <- function(con, force=TRUE) {

  if (!("sqrt_cov" %in% con) || force) {
    sol = svd(con$covariance)
    U = sol$u
    D = sol$d
    D = sqrt(D)

    sqrt_cov = U * as.matrix(D, nrow=1)
    sqrt_inv = t(U / as.matrix(D, nrow=1))
    rowspace = U
  }

  return_list <- list("sqrt_cov" = sqrt_cov,
                      "sqrt_inv" = sqrt_inv,
                      "rowspace" = rowspace)
  return(return_list)
}
