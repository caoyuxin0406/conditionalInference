#' Affine Gaussian Sampler
#'
#' Sample from an affine truncated Gaussian distribution.
#'
#' @param affine_con The box constraint
#' @param initial_point Initial solution from optimization problem
#' @param observed_score_state The negated S matrix
#' @param log_density log_density function
#' @param logdens_linear Described how score enters log_density (linear part)
#' @param opt_offset Described how score enters log_density (offset part)
#' @param selection_info A list of (opt_sampler, opt_sample, target_cov, score_cov) objects
#'
#' @return An affine_gaussian sampler object, which is a list:
#' \item{affine_con}{The box constraint}
#' \item{initial_point}{Initial solution from optimization problem}
#' \item{observed_score_state}{The negated S matrix}
#' \item{selection_info}{A list of (opt_sampler, opt_sample, target_cov, score_cov) objects}
#' \item{log_density}{log_density function}
#' \item{logdens_linear}{Described how score enters log_density (linear part)}
#' \item{opt_offset}{Described how score enters log_density (offset part)}
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
#' # Build a new partial log_density function
#' new_log_density = pryr::partial(log_density,
#'                                 logdens_linear=model$logdens_linear,
#'                                 offset=model$opt_offset,
#'                                 cond_prec=model$cond_precision,
#'                                 active_directions=model$active_directions,
#'                                 lagrange=model$lagrange,
#'                                 Z=Z)
#'
#' # Construct the affine gaussian sampler
#' affine_gaussian_sampler(affine_con = model$affine_con,
#'                         initial_point = model$observed_opt_state,
#'                         observed_score_state = model$observed_score_state,
#'                         log_density = new_log_density,
#'                         logdens_linear = model$logdens_linear,
#'                         opt_offset = model$opt_offset)
affine_gaussian_sampler <- function(affine_con,
                                    initial_point,
                                    observed_score_state,
                                    log_density,
                                    logdens_linear,
                                    opt_offset,
                                    selection_info=NULL) {

  return_list <- list("affine_con" = affine_con,
                      "initial_point" = initial_point,
                      "observed_score_state" = observed_score_state,
                      "selection_info" = selection_info,
                      "log_density" = log_density,
                      "logdens_linear" = logdens_linear,
                      "opt_offset" = opt_offset)
  return(return_list)
}

#' Sample From Constraints
#'
#' Use Gibbs sampler to simulate from "con".
#'
#' @param con Constraints
#' @param Y Point satisfying the constraint
#' @param direction_of_interest Represents which projection is of the most interest
#' @param how_often Represents how often would the sampler make a move along "direction_of_interest"
#'                  If negative, defaults to ndraw+burnin (so it will never be used).
#' @param ndraw Defaults to 1000
#' @param burnin Defaults to 1000
#' @param white Binary indicator of whether con$covariance is equal to identity
#' @param use_constraint_directions Binary indicator of whether using the directions
#'                                  formed by the constraints as in the Gibbs scheme
#' @param use_random_directions Binary indicator of whether using additional
#'                              random directions in the Gibbs scheme
#' @param accept_reject_params If not NULL should be a tuple (num_trial, min_accept, num_draw).
#'                             In this case, we first try num_trial accept-reject samples,
#'                             if at least min_accept of them succeed, we just draw num_draw
#'                             accept_reject samples.
#'
#' @return Z matrix - sample from the Gaussian distribution conditioned
#'         on the constraints
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
#' \link{whiten} for whitening the box constraint.
#'
#' \link{forward_map} for constructing a forward_map.
#'
#' \link{accept_reject} for testing random samples.
#'
#' \link{sample_truncnorm_white} for sampling from a truncated normal
#' with white constraint.
#'
#' \link{inverse_map} for constructing an inverse_map.
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
#' # Sample from the box constraint of the optimization problem
#' sample_from_constraints(con=model$sampler$affine_con,
#'                         Y=model$sampler$initial_point)
sample_from_constraints <- function(con,
                                    Y,
                                    direction_of_interest=NULL,
                                    how_often=-1,
                                    ndraw=1000,
                                    burnin=1000,
                                    white=FALSE,
                                    use_constraint_directions=TRUE,
                                    use_random_directions=TRUE,
                                    accept_reject_params=NULL) {

  if (is.null(direction_of_interest)) {
    if (is.null(dim(Y))) {
      as.matrix(Y)
    }
    direction_of_interest = rnorm(length(Y))
  }
  if (how_often < 0) {
    how_often = ndraw + burnin
  }

  if (!white) {
    white_con = whiten(con)
    white_Y = conditionalInference::forward_map(Y, white_con)
    white_direction_of_interest =
      conditionalInference::forward_map(con$covariance %*% direction_of_interest,
                                        white_con)
  } else {
    white_con = con
    inverse_map = function(x, con) x
  }

  # try 100 draws of accept reject
  # if we get more than 50 good draws, then just return a smaller sample
  # of size (burnin+ndraw)/5
  if (!is.null(accept_reject_params)) {
    use_hit_and_run = FALSE
    num_trial = accept_reject_params[1]
    min_accept = accept_reject_params[2]
    num_draw = accept_reject_params[3]

    Z_sample = accept_reject(100,
                             white_con$linear_part,
                             white_con$offset)

    if (dim(Z_sample)[1] >= min_accept) {
      while (TRUE) {
        new_sample = accept_reject(num_draw / 5,
                                   white_con$linear_part,
                                   white_con$offset)
        Z_sample = rbind(Z_sample, new_sample)
        if (dim(Z_sample)[1] > num_draw) {
          break
        }
      }
      white_samples = Z_sample
    } else {
      use_hit_and_run = TRUE
    }
  } else {
    use_hit_and_run = TRUE
  }

  if (use_hit_and_run) {
    white_samples = sample_truncnorm_white(white_con$linear_part,
                                           white_con$offset,
                                           white_Y,
                                           white_direction_of_interest,
                                           how_often=how_often,
                                           ndraw=ndraw,
                                           burnin=burnin,
                                           sigma=1,
                                           use_constraint_directions=use_constraint_directions,
                                           use_random_directions=use_random_directions)
  }

  Z = conditionalInference::inverse_map(white_samples, white_con)
  return(Z)
}


#' Accept or Reject Samples
#'
#' Random samples from normal distribution that satisfy the box constraint.
#'
#' @param sample_size Number of samples generated
#' @param linear_part Linear part of the constraint
#' @param offset Offset part of the constraint
#'
#' @return Random samples that satisfy the box constraint.
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
#' # Testing samples
#' white_con = whiten(model$sampler$affine_con)
#' accept_reject(100,
#'               white_con$linear_part,
#'               white_con$offset)
accept_reject <- function(sample_size,
                          linear_part,
                          offset) {

  Z_sample = matrix(rnorm(sample_size), ncol = dim(linear_part)[2])
  temp = apply(Z_sample %*% t(linear_part), c(1,2), function(x) x - offset)
  constraint_satisfied = apply(temp, 1, max) < 0

  result = as.matrix(Z_sample[constraint_satisfied], ncol = 1)
  return(result)
}


#' Sample Whiten Truncated Normal
#'
#' Sample from a truncated normal with covariance equal to sigma^2 I.
#' Constraint is $Ax <= b$ where `A` has shape `(q,n)` with `q` the number
#' of constraints and `n` the number of random variables.
#'
#' @param A Linear part of affine constraints
#' @param b Offset part of affine constraints
#' @param initial Initial point for Gibbs draws, assumed to satisfy the constraints
#' @param bias_direction Represents the projection that is of most interest
#' @param how_often Indicates how often should the sampler make a move along "direction_of_interest"
#'                  If negative, defaults to ndraw+burnin (so it will never be used)
#' @param sigma Variance parameter
#' @param burnin Indicates how many iterations until we start recording samples
#' @param ndraw Inidicates how many samples should we return
#' @param use_constraint_directions Binary indicator of whether using the directions
#'                                  formed by the constraints as in the Gibbs scheme
#' @param use_random_directions Binary indicator of whether using additional
#'                              random directions in the Gibbs scheme
#' @param ignore_bound_violations Binary indicator of whether ignoring bound violations
#'
#' @return Truncated normal samples.
#' @export
#'
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom pracma dot
#'
#' @references
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
#' Plausibly Valid Instruments with Application to Mendelian Randomization.
#'
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
#' Effects After Testing Instrument Strength in Linear Models.
#'
#' @seealso
#' \link[stats]{pnorm} for the distribution function of normal distribution.
#'
#' \link[stats]{qnorm} for the quantile function of normal distribution.
#'
#' \link[stats]{rnorm} for generating random deviates of normal distribution.
#'
#' \link[stats]{runif} for generating random deviates of uniform distribution.
#'
#' \link{sample} for sampling process.
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
#' white_Y = forward_map(model$sampler$initial_point,
#'                       white_con)
#' Y = as.matrix(model$sampler$initial_point)
#' direction_of_interest = rnorm(length(Y))
#' white_direction_of_interest =
#' forward_map(model$sampler$affine_con$covariance %*%
#'             direction_of_interest,
#'             white_con)
#'
#' # Generate whiten samples from truncated normal distribution
#' sample_truncnorm_white(white_con$linear_part,
#'                        white_con$offset,
#'                        white_Y,
#'                        white_direction_of_interest,
#'                        how_often=2000,
#'                        ndraw=1000,
#'                        burnin=1000,
#'                        sigma=1,
#'                        use_constraint_directions=TRUE,
#'                        use_random_directions=TRUE)
sample_truncnorm_white <- function(A,
                                   b,
                                   initial,
                                   bias_direction, #eta
                                   how_often=1000,
                                   sigma=1,
                                   burnin=500,
                                   ndraw=1000,
                                   use_constraint_directions=TRUE,
                                   use_random_directions=FALSE,
                                   ignore_bound_violations=TRUE) {

  nvar = dim(A)[2]
  nconstraint = dim(A)[1]
  tol = 1.e-7
  trunc_sample = array(0, c(ndraw, nvar))

  state = initial
  U = pracma::dot(A, state) - b
  usample = runif(burnin + ndraw)

  # directions not parallel to coordinate axes
  dirs = list()
  if (use_constraint_directions) {
    dirs = append(dirs, A)
  }
  if (use_random_directions) {
    rand = rnorm(as.integer(nvar/5), nvar)
    dirs = append(dirs, as.matrix(rand))
  }

  b = array(bias_direction, dim = length(bias_direction))
  dirs = append(dirs, t(t(b)))

  directions = as.matrix(unlist(dirs))
  directions = directions / t(t(sqrt(rowSums(directions^2))))
  ndir = dim(directions)[1]
  alphas_dir = A %*% t(directions)
  alphas_coord = A

  colMax <- function(data) sapply(data, max)
  alphas_max_dir = colMax(abs(alphas_dir)) * tol
  alphas_max_coord = colMax(abs(alphas_coord)) * tol

  random_idx_dir = sample(x=1:ndir, size=(burnin+ndraw), replace=TRUE)
  random_idx_coord = sample(x=1:nvar, size=(burnin+ndraw), replace=TRUE)

  invperiod = 14
  docoord = 0
  iperiod = 1
  ibias = 1
  dobias = 0
  make_no_move = 0
  restart_idx = 1

  for (iter_count in 1:(ndraw + burnin)) {

    make_no_move = 0

    docoord = 1
    iperiod = iperiod + 1
    ibias = ibias + 1

    if (iperiod == invperiod) {
      docoord = 0
      iperiod = 0
      dobias = 0
    }

    if (ibias == how_often) {
      docoord = 0
      ibias = 0
      dobias = 1
    }

    if (docoord == 1) {
      idx = random_idx_coord[iter_count]
      V = state[idx]
    } else {
      if (!dobias) {
        idx = random_idx_dir[iter_count]
      } else {
        idx = dim(directions)[1]-1 # last row of directions is bias_direction
        V = 0
        for (ivar in 1:nvar) {
          V = V + directions[idx, ivar] * state[ivar]
        }
      }
    }

    lower_bound = -1e12
    upper_bound = 1e12
    for (irow in 1:nconstraint) {
      if (docoord == 1) {
        alpha = alphas_coord[irow,idx]
        val = -U[irow] / alpha + V
        if ((alpha > alphas_max_coord[idx]) & (val < upper_bound)) {
          upper_bound = val
        } else if ((alpha < -alphas_max_coord[idx]) & (val > lower_bound)) {
          lower_bound = val
        }
      } else {
        alpha = alphas_dir[irow,idx]
        val = -U[irow] / alpha + V
        if ((alpha > alphas_max_dir[idx]) & (val < upper_bound)) {
          upper_bound = val
        } else if ((alpha < -alphas_max_dir[idx]) & (val > lower_bound)) {
          lower_bound = val
        }
      }
    }

    if (lower_bound > V) {
      lower_bound = V - tol * sigma
    } else if (upper_bound < V) {
      upper_bound = V + tol * sigma
    }

    lower_bound = lower_bound / sigma
    upper_bound = upper_bound / sigma

    if (lower_bound > upper_bound) {
      if (!ignore_bound_violations) {
        print("BoundViolation")
      } else {
        make_no_move = 1
      }
      if (iter_count - burnin > 0) {
        restart_idx = iter_count - burnin / 2
        for (ivar in 1:nvar) {
          state[ivar] = trunc_sample[restart_idx, ivar]
        }
      } else {
        for (ivar in 1:nvar) {
          state[ivar] = initial[ivar]
        }
      }
    }

    if (upper_bound < -10) { # use Exp approximation
      # the approximation is that
      # Z | lower_bound < Z < upper_bound
      # is fabs(upper_bound) * (upper_bound - Z) = E approx Exp(1)
      # so Z = upper_bound - E / fabs(upper_bound)
      # and the truncation of the exponential is
      # E < fabs(upper_bound - lower_bound) * fabs(upper_bound) = D

      # this has distribution function (1 - exp(-x)) / (1 - exp(-D))
      # so to draw from this distribution
      # we set E = - log(1 - U * (1 - exp(-D))) where U is Unif(0,1)
      # and Z (= tnorm below) is as stated

      unif = usample[iter_count] * (1 - exp(-abs(
        (lower_bound - upper_bound) * upper_bound)))
      tnorm = (upper_bound + log(1 - unif) / abs(upper_bound)) * sigma
    } else if (lower_bound > 10) {
      # here Z = lower_bound + E / fabs(lower_bound) (though lower_bound is positive)
      # and D = fabs((upper_bound - lower_bound) * lower_bound)
      unif = usample[iter_count] * (1 - exp(-abs(
        (upper_bound - lower_bound) * lower_bound)))
      tnorm = (lower_bound - log(1 - unif) / lower_bound) * sigma
    } else if (lower_bound < 0) {
      cdfL = pnorm(lower_bound)
      cdfU = pnorm(upper_bound)
      unif = usample[iter_count] * (cdfU - cdfL) + cdfL
      if (unif < 0.5) {
        tnorm = qnorm(unif) * sigma
      } else {
        tnorm = -qnorm(1-unif) * sigma
      }
    } else {
      cdfL = pnorm(-lower_bound)
      cdfU = pnorm(-upper_bound)
      unif = usample[iter_count] * (cdfL - cdfU) + cdfU
      if (unif < 0.5) {
        tnorm = -qnorm(unif) * sigma
      } else {
        tnorm = qnorm(1-unif) * sigma
      }
    }

    if (docoord == 1) {
      state[idx] = tnorm
      tnorm = tnorm - V
      for (irow in 1:nconstraint) {
        U[irow] = U[irow] + tnorm * A[irow, idx]
      }
    } else {
      tnorm = tnorm - V
      for (ivar in 1:nvar) {
        state[ivar] = state[ivar] + tnorm * directions[idx,ivar]
        for (irow in 1:nconstraint) {
          U[irow] = (U[irow] + A[irow, ivar] *
                       tnorm * directions[idx,ivar])
        }
      }
    }
    if (iter_count >= burnin & !make_no_move) {
      for (ivar in 1:nvar) {
        trunc_sample[iter_count - burnin, ivar] = state[ivar]
      }
    }
  }

  return(trunc_sample)
}


#' Coefficient pvalues with IV for TSLS Test Statistic
#'
#' Construct selective p-values for each parameter of the target. The p-values
#' are computed for the optimization problem with instrumental variables.
#' We compute the p-values of the Two-Stage Least-Squares (TSLS) test statistic.
#'
#' @param observed_target A vector of parameters, representing coordinates of the target
#' @param target_cov Covariance of the target (e.g. TSLS)
#' @param score_cov Covariance of the target score (e.g. TSLS score)
#' @param opt_sampler A sampler generated from the optimization problem
#' @param parameter A vector of parameters at which to evaluate p-values. Defaults
#'                  to the zeros matrix
#' @param sample_args Default to NULL
#' @param sample Represents a sample of the target. Allows reuse of the same sample
#'               for construction of confidence intervals, hypothesis tests, etc.
#' @param alternatives A list of ['greater', 'less', 'twosided'],
#'                     which indicates what alternative to use
#'
#' @return Sampled p-values.
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
#' \link{coefficient_pvalues_iv_ar} that computes selective p-values for AR
#' statistic. The method to compute p-values for TSLS and AR statistics are similar.
#'
#' \link{optimization_intervals} to initialize an interval for the optimization
#' problem, which specify the observed target and its covariance as well as
#' optimization samples and sampler. We need to initialize an optimization_intervals
#' before computing p-values.
#'
#' \link{pivot_iv} for calculating the pivot quantities of the optimization
#' intervals before computing p-values.
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
#' # Estimate covariance and obtain summary statistics
#' cov = estimate_covariance(Y, D, Z)
#' s <- summary_tsls(gl, model$sampler,
#'                   Sigma_11 = cov[1,1], Sigma_12 = cov[1,2])
#'
#' # Compute p-values for TSLS test statistic
#' pivots = coefficient_pvalues_iv(observed_target=s$observed_target,
#'                                 target_cov=s$cov_target,
#'                                 score_cov=s$cov_target_score,
#'                                 opt_sampler=model$sampler,
#'                                 parameter=0,
#'                                 sample=s$opt_sample,
#'                                 alternatives=rep('twosided',1))
#' pivots
coefficient_pvalues_iv <- function(observed_target, # matrix
                                   target_cov,
                                   score_cov,
                                   opt_sampler,
                                   parameter=NULL,  # matrix
                                   sample_args=NULL,
                                   sample, # matrix
                                   alternatives=NULL) {

  if (is.null(alternatives)) {
    alternatives = rep('twosided', dim(observed_target)[1])
  }

  ndraw = dim(sample)[1]

  if (is.null(parameter)) {
    parameter = rep(0, dim(observed_target)[1])
  }

  opt_sampling_info = list(sets::tuple(opt_sampler, sample, target_cov, score_cov))
  intervals = conditionalInference::optimization_intervals(opt_sampling_info,
                                                           observed_target,
                                                           ndraw)

  #pvals = c()

  for (i in 1:nrow(observed_target)) {
    keep = matrix(0, nrow = nrow(observed_target), ncol = ncol(observed_target))
    keep[i] = 1
    curr_pval = conditionalInference::pivot_iv(opt_intervals=intervals,
                                               linear_func=keep,
                                               candidate=parameter[i],
                                               alternative=alternatives[i])
    #print(curr_pval)

    if (i == 1) {
      pvals = list(curr_pval)
    } else {
      pvals = append(pvals, curr_pval)
    }
  }

  return(array(pvals))
}


#' Coefficient pvalues with Instrumental Variables for AR Test Statistic
#'
#' Construct selective p-values for each parameter of the target. The p-values
#' are computed for the optimization problem with instrumental variables.
#' We compute the p-values of the Anderson-Rubin (AR) test statistic.
#'
#' @param observed_target A vector of parameters, representing coordinates of the target
#' @param target_cov Covariance of the target (e.g. AR)
#' @param score_cov Covariance of the target score (e.g. AR score)
#' @param linear_func Linear function
#' @param gl_ar The group lasso object for AR test statistic representing the
#'              penalized convex optimization equation (generated by
#'              group_lasso_iv function)
#' @param opt_sampler A sampler generated from the optimization problem
#' @param parameter A vector of parameters at which to evaluate p-values. Defaults
#'                  to the zeros matrix
#' @param sample_args Default to NULL
#' @param sample Represents a sample of the target. Allows reuse of the same sample
#'               for construction of confidence intervals, hypothesis tests, etc.
#' @param alternatives A list of ['greater', 'less', 'twosided'],
#'                     which indicates what alternative to use
#'
#' @return Sample p-values
#' @export
#'
#' @importFrom sets tuple
#'
#' @references
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
#' Plausibly Valid Instruments with Application to Mendelian Randomization.
#'
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
#' Effects After Testing Instrument Strength in Linear Models.
#'
#' @seealso
#' \link{coefficient_pvalues_iv} that computes selective p-values for TSLS
#' statistic. The method to compute p-values for TSLS and AR statistics are similar.
#'
#' \link{optimization_intervals} to initialize an interval for the optimization
#' problem, which specify the observed target and its covariance as well as
#' optimization samples and sampler. We need to initialize an optimization_intervals
#' before computing p-values.
#'
#' \link{pivot_iv_ar} for calculating the pivot quantities of the optimization
#' intervals before computing p-values.
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
#' gl_ar <- group_lasso_iv_ar(Y, D, Z)
#' model_ar <- fit_ar(gl_ar)
#'
#' # Estimate covariance and obtain summary statistics
#' cov = estimate_covariance(Y, D, Z)
#' s <- summary_ar(gl_ar, model_ar$sampler,
#' Sigma_11 = cov[1,1], Sigma_12 = cov[1,2])
#'
#' # Compute p-values for TSLS test statistic
#' coefficient_pvalues_iv_ar(observed_target=s$observed_target,
#'                           target_cov=s$cov_target,
#'                           score_cov=s$cov_target_score,
#'                           linear_func=s$K2,
#'                           gl_ar=gl_ar,
#'                           opt_sampler=model_ar$sampler,
#'                           parameter=0,
#'                           sample=s$opt_sample,
#'                           alternatives=rep('greater',1))
coefficient_pvalues_iv_ar <- function(observed_target,
                                      target_cov,
                                      score_cov,
                                      linear_func,
                                      gl_ar,
                                      opt_sampler,
                                      parameter=NULL,
                                      sample_args=NULL,
                                      sample,
                                      alternatives=NULL) {

  if (is.null(alternatives)) {
    alternatives = rep('greater', 1)
  }

  ndraw = dim(sample)[1]

  if (is.null(parameter)) {
    parameter = rep(0, 1)
  }

  opt_sampling_info = list(sets::tuple(opt_sampler, sample, target_cov, score_cov))
  intervals = conditionalInference::optimization_intervals(opt_sampling_info,
                                                           observed_target,
                                                           ndraw)

  curr_pval = conditionalInference::pivot_iv_ar(opt_intervals=intervals,
                                                linear_func=linear_func,
                                                parameter_reference=parameter,
                                                candidate=0,
                                                gl_ar=gl_ar,
                                                alternative='greater')
  pvals = list(curr_pval)

  return(array(pvals))
}


