#' Optimization Intervals
#'
#' Initialize the optimization_intervals object, which specify the observed
#' target, target covariance, optimization samples, and optimization sampler.
#'
#' @param opt_sampling_info A sequence of (opt_sampler, opt_sample, target_cov,
#' score_cov) objects
#' @param observed Observed samples
#' @param nsample How large a normal sample
#' @param target_cov Covariance of the target (e.g. TSLS)
#'
#' @return An optimization_intervals object, which is a list of:
#' \item{blahvals}{Repeat opt_samples if necessary}
#' \item{opt_sampling_info}{A sequence of (opt_sampler, opt_sample, target_cov,
#' score_cov) objects}
#' \item{logden}{Total log density of the optimization sampler}
#' \item{observed}{Observed samples}
#' \item{target_cov}{Covariance of the target}
#' \item{normal_sample}{Samples generated from the multivariate normal distribution}
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
#' \link{log_density} for calculating the log density of a given equation,
#' including the Jacobian term.
#'
#' \link[MASS]{mvrnorm} for sampling from multivariate normal distribution.
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
#' # Construct the optimization intervals before computing confidence intervals
#' # for TSLS test statistic
#' ndraw = nrow(s$opt_sample)
#' opt_sampling_info = list(sets::tuple(model$sampler,
#'                                      s$opt_sample,
#'                                      s$cov_target,
#'                                      s$cov_target_score))
#' optimization_intervals(opt_sampling_info,
#'                        s$observed_target,
#'                        ndraw)
optimization_intervals <- function(opt_sampling_info, # a list of
                                   # (opt_sampler,
                                   #  opt_sample,
                                   #  target_cov,
                                   #  score_cov) tuple objects
                                   #  in theory all target_cov
                                   #  should be about the same...
                                   observed,
                                   nsample, # how large a normal sample
                                   target_cov=NULL) {

  blahvals = c()
  # not all opt_samples will be of the same size as nsample
  # let's repeat them as necessary

  #tiled_sampling_info = list()
  idx = 1
  for (t in opt_sampling_info) {
    opt_sampler = t[[1]]
    opt_sample = t[[2]]
    t_cov = t[[3]]
    score_cov = t[[4]]

    if (!is.null(opt_sample)) {
      if (dim(opt_sample)[1] < nsample) {
        if (length(dim(opt_sample)) == 1) {
          tiled_opt_sample = rep(opt_sample, as.integer(ceiling(nsample / dim(opt_sample)[1])))[1:nsample]
        } else {
          tiled_opt_sample = as.matrix(rep(opt_sample, as.integer(ceiling(nsample / dim(opt_sample)[1]))), ncol=1)[1:nsample]
        }
      } else {
        tiled_opt_sample = as.matrix(opt_sample[1:nsample])
      }
    } else {
      tiled_sample = NULL
    }
    #tiled_sampling_info.append((opt_sampler, tiled_opt_sample, t_cov, score_cov))
    tiled_t = sets::tuple(opt_sampler, tiled_opt_sample, t_cov, score_cov)
    #tiled_t = t
    if (idx == 1) {
      tiled_sampling_info = list(tiled_t)
    } else {
      tiled_sampling_info = append(tiled_sampling_info, tiled_t)
    }
    idx = idx + 1
  }

  opt_sampling_info = tiled_sampling_info
  logden = rep(0, nsample)
  for (t in opt_sampling_info) {
    opt_sampler = t[[1]]
    opt_sample = t[[2]]
    currden = opt_sampler$log_density(opt_sampler$observed_score_state,
                                      opt_sample)
    logden = logden + currden
    if (dim(opt_sample)[1] < nsample) {
      logden = rep(logden, as.integer(ceiling(nsample / dim(opt_sample)[1])))[1:nsample]
    }
  }

  # average covariances in case they might be different

  if (is.null(target_cov)) {
    target_cov = 0
    for (t in opt_sampling_info) {
      curr_target_cov = t[[3]]
      target_cov = target_cov + curr_target_cov
    }
    target_cov = target_cov / length(opt_sampling_info)
  }

  normal_sample = MASS::mvrnorm(n = nsample,
                                mu = rep(0, dim(target_cov)[1]),
                                Sigma = target_cov)

  return_list <- list("blahvals" = blahvals,
                      "opt_sampling_info" = opt_sampling_info,
                      "logden" = logden,
                      "observed" = observed,
                      "target_cov" = target_cov,
                      "normal_sample" = normal_sample)
  return(return_list)
}


#' Pivot
#'
#' Calculate the pivot quantities of the optimization intervals.
#'
#' @param opt_intervals Optimization intervals
#' @param linear_func Linear function
#' @param candidate Candidate value
#' @param alternative A list of ['greater', 'less', 'twosided'],
#'                    which indicates what alternative to use
#'
#' @return Pivot quantities for TSLS test statistic.
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
#' \link{pivot_iv_ar} for calculating the pivot quantities/p-values of the
#' optimization intervals for AR test statistic. The process is similar to
#' that of TSLS test statistic.
#'
#' \link{weights_iv} that computes the weights of candidate pivot values
#' for TSLS test statistic.
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
#'                   Sigma_11=cov[1,1], Sigma_12=cov[1,2])
#'
#' # Construct the optimization intervals before computing confidence intervals
#' # for TSLS test statistic
#' ndraw = nrow(s$opt_sample)
#' opt_sampling_info = list(sets::tuple(model$sampler,
#'                                      s$opt_sample,
#'                                      s$cov_target,
#'                                      s$cov_target_score))
#' intervals <- optimization_intervals(opt_sampling_info,
#'                                     s$observed_target,
#'                                     ndraw)
#'
#' # Compute pivot value
#' keep = matrix(0, nrow(s$observed_target), ncol(s$observed_target))
#' keep[1] = 1
#' pivot_iv(opt_intervals=intervals,
#'          linear_func=keep,
#'          candidate=0,
#'          alternative='twosided')
#'
#' # Compute the pivot value for a sequence of candidate value
#' level = 0.95
#' sample_stat = intervals$normal_sample %*% keep
#' observed_stat = intervals$observed %*% keep
#' grid_min = -20 * sd(sample_stat)
#' grid_max = 20 * sd(sample_stat)
#' for (gamma in seq(grid_min, grid_max, length.out=40)) {
#'  candidate = gamma + observed_stat[1,1]
#'  lower <- pivot_iv(intervals, keep, candidate, 'less') - (1 - level) / 2
#'  upper <- pivot_iv(intervals, keep, candidate, 'less') - (1 + level) / 2
#'  print(c(lower, upper))
#' }
pivot_iv <- function(opt_intervals,
                     linear_func,
                     candidate,
                     alternative='twosided') {

  if (alternative != 'greater' && alternative != 'less' && alternative != 'twosided') {
    print("alternative should be one of ['greater', 'less', 'twosided']")
  }

  observed_stat = opt_intervals$observed %*% linear_func
  #ns = read.table("normal_sample.txt")
  #normal_sample = as.matrix(ns)
  #sample_stat = normal_sample %*% linear_func
  sample_stat = opt_intervals$normal_sample %*% linear_func

  target_cov = linear_func %*% (opt_intervals$target_cov %*% linear_func)

  idx = 1
  for (t in opt_intervals$opt_sampling_info) {
    opt_sampler = t[[1]]
    opt_sample = t[[2]]
    score_cov = t[[4]]

    # score_cov is a 1xp matrix
    cur_score_cov = linear_func %*% t(score_cov)

    # cur_nuisance is in the view's score coordinates
    offset = apply(cur_score_cov, c(1,2), function(x) x * observed_stat / target_cov)
    cur_nuisance = opt_sampler$observed_score_state - t(offset) # px1 matrix
    cur_translate_dirs = t(apply(cur_score_cov, c(1,2), function(x) x / target_cov))
    if (idx == 1) {
      nuisance = list(cur_nuisance)
      translate_dirs = list(cur_translate_dirs)
    } else {
      nuisance = append(nuisance, cur_nuisance)
      translate_dirs = append(translate_dirs, cur_translate_dirs)
    }
    idx = idx + 1
  }

  weights = conditionalInference::weights_iv(opt_intervals,
                                             sample_stat,  # normal sample
                                             candidate,    # candidate value
                                             nuisance,       # nuisance sufficient stats for each view
                                             translate_dirs) # points will be moved like sample * score_cov

  comp = apply(sample_stat + candidate, c(1,2), function(x) x <= observed_stat)
  pivot = mean(comp * weights) / mean(weights)

  if (alternative == 'twosided') {
    return(2 * min(pivot, 1 - pivot))
  } else if (alternative == 'less') {
    return(pivot)
  } else {
    return(1 - pivot)
  }
}


#' Weights
#'
#' Calculate the weights of candidate pivot values for TSLS test statisitc.
#'
#' @param opt_intervals Optimization intervals
#' @param sample_stat Normal sample
#' @param candidate Candidate value
#' @param nuisance Nuisance sufficient stats for each view
#' @param translate_dirs Points will be moved like sample * score_cov
#'
#' @return Weights of the sequence of pivot values
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
#' \link{weights_iv_ar} for calculating the weights of candidate pivot values
#' for AR test statistic. This process is similar to that for TSLS statistic.
#'
#' \link{log_density_ray} that computes the log density of normal sample.
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
#' s <- summary_tsls(gl, model$sampler, Sigma_11=cov[1,1],
#'                   Sigma_12=cov[1,2], ndraw=20, burnin=20)
#'
#' # Construct the optimization intervals before computing confidence intervals
#' # for TSLS test statistic
#' ndraw = nrow(s$opt_sample)
#' opt_sampling_info = list(sets::tuple(model$sampler,
#'                                      s$opt_sample,
#'                                      s$cov_target,
#'                                      s$cov_target_score))
#' intervals <- optimization_intervals(opt_sampling_info,
#'                                     s$observed_target,
#'                                     ndraw)
#'
#' # Setting up parameters for weights_iv method
#' keep = matrix(0, nrow(s$observed_target), ncol(s$observed_target))
#' keep[1] = 1
#' observed_stat = intervals$observed %*% keep
#' sample_stat = intervals$normal_sample %*% keep
#' cur_nuisance = matrix(c(0.6197473, 1.3126907,  0.8292342))
#' cur_translate_dirs = matrix(c(-3.048598, -6.457254, -4.079085))
#' nuisance = list(cur_nuisance)
#' translate_dirs = list(cur_translate_dirs)
#'
#' # Compute weights for candidate pivot values
#' weights <- weights_iv(opt_intervals=intervals,
#'                       sample_stat=sample_stat,
#'                       candidate=0,
#'                       nuisance=nuisance,
#'                       translate_dirs=translate_dirs)
#'
#' # Using the computed weights to calculate the pivot value
#' comp <- apply(sample_stat + 0, c(1,2), function(x) x <= observed_stat)
#' pivot <- mean(comp * weights) / mean(weights)
#' pivot
weights_iv <- function(opt_intervals, # self
                       sample_stat,
                       candidate,
                       nuisance, # list of matrices
                       translate_dirs) { # list of matrices

  score_sample = c()
  lognum = 0

  i = 1
  for (t in opt_intervals$opt_sampling_info) {
    opt_sampler = t[[1]]
    opt_sample = t[[2]]
    logden = conditionalInference::log_density_ray(candidate=candidate,
                                                   direction=translate_dirs[[i]],
                                                   nuisance=nuisance[[i]],
                                                   gaussian_sample=sample_stat,
                                                   opt_sampler=opt_sampler,
                                                   opt_sample=opt_sample)
    lognum = lognum + logden$lognum
    i = i + 1
  }

  logratio = lognum - opt_intervals$logden
  logratio = logratio - max(logratio)

  return(exp(logratio))
}


#' Log Density Ray
#'
#' Calculate the log density of normal sample.
#'
#' @param candidate Candidate value
#' @param direction Points will be moved like sample * score_cov
#' @param nuisance Nuisance sufficient stats for each view
#' @param gaussian_sample Gaussian/normal samples
#' @param opt_sampler Optimization sampler (e.g. affine gaussian sampler)
#' @param opt_sample Samples generated from the optimization problem
#'
#' @return A log_density_ray object, which is a list of:
#' \item{linear_term}{Linear term}
#' \item{quadratic_term}{Quadratic term}
#' \item{constant_term}{Constant term}
#' \item{lognum}{return value of log_density_ray}
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
#' \link{log_density} to calculate the log density of given equation,
#' including the Jacobian term.
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
#' s <- summary_tsls(gl, model$sampler, Sigma_11=cov[1,1],
#'                   Sigma_12=cov[1,2], ndraw=20, burnin=20)
#'
#' # Construct the optimization intervals before computing confidence intervals
#' # for TSLS test statistic
#' opt_sampling_info = list(sets::tuple(model$sampler,
#'                                      s$opt_sample,
#'                                      s$cov_target,
#'                                      s$cov_target_score))
#' intervals <- optimization_intervals(opt_sampling_info,
#'                                     s$observed_target,
#'                                     nrow(s$opt_sample))
#'
#' # Setting up parameters for weights_iv method
#' keep = matrix(0, nrow(s$observed_target), ncol(s$observed_target))
#' keep[1] = 1
#' sample_stat = intervals$normal_sample %*% keep
#' cur_nuisance = matrix(c(0.6197473, 1.3126907,  0.8292342))
#' cur_translate_dirs = matrix(c(-3.048598, -6.457254, -4.079085))
#' nuisance = list(cur_nuisance)
#' translate_dirs = list(cur_translate_dirs)
#'
#' # Calculate the log density of normal sample
#' log_density_ray(candidate=0,
#'                 direction=translate_dirs[[1]],
#'                 nuisance=nuisance[[1]],
#'                 gaussian_sample=sample_stat,
#'                 opt_sampler=model$sampler,
#'                 opt_sample=s$opt_sample)
log_density_ray <- function(candidate,
                            direction,
                            nuisance,
                            gaussian_sample,
                            opt_sampler,
                            opt_sample) {

  if (is.null(opt_sampler$direction) || (opt_sampler$direction != direction)) {

    logdens_lin = opt_sampler$logdens_linear
    logdens_offset = opt_sampler$opt_offset

    if (ncol(opt_sample) == 1) {

      prec = 1 / opt_sampler$affine_con$covariance
      quadratic_term = (logdens_lin %*% direction)**2 * prec
      temp_mat = apply(gaussian_sample, c(1,2), function(x) x * (logdens_lin %*% direction)) +
        opt_sample
      temp_off = logdens_lin %*% (nuisance + logdens_offset)
      arg = apply(temp_mat, c(1,2), function(x) x + temp_off)
      temp = apply(logdens_lin %*% direction, c(1, 2), function(x) x * prec)
      linear_term = apply(arg, c(1,2), function(x) x * temp)
      constant_term = arg**2 %*% prec

    } else {

      logdens_lin = opt_sampler$logdens_linear
      logdens_offset = opt_sampler$opt_offset

      cov = opt_sampler$affine_con$covariance
      prec = solve(cov)
      linear_part = logdens_lin %*% direction # A gamma

      if (nrow(opt_sample) == 1 || ncol(opt_sample) == 1) {
        # pass/do nothing
      }
      cov = opt_sampler$affine_con$covariance

      quadratic_term = (t(linear_part) %*% prec) %*% linear_part

      arg1 = t(opt_sample)
      arg2 = logdens_lin %*% (outer(direction, gaussian_sample) +
                               as.vector(nuisance + logdens_offset))
      arg = arg1 + arg2
      linear_term = t(linear_part) %*% prec %*% arg
      temp = apply(arg, c(1,2), function(x) (prec %*% arg) * x)
      constant_term = colSums(temp)
    }
  }

  temp_ans = candidate * linear_term + 0.5 * constant_term
  temp_off = -0.5 * candidate**2 * quadratic_term
  lognum = apply(temp_ans, c(1,2), function(x) temp_off - x)

  return_list <- list("linear_term" = linear_term,
                      "quadratic_term" = quadratic_term,
                      "constant_term" = constant_term,
                      "lognum" = lognum)
  return(return_list)
}


#' Confidence Intervals
#'
#' Construct selective confidence intervals for each observed target.
#' The confidence intervals are computed for the optimization problem with
#' instrumental variables. We compute the confidence intervals of the
#' Two-Stage Least-Squares (TSLS) test statistic.
#'
#' @param observed_target A vector of parameters, representing coordinates of the target
#' @param target_cov Covariance of the target (e.g. TSLS)
#' @param score_cov Covariance of the target score (e.g. TSLS score)
#' @param sample_args Default to NULL (optional)
#' @param opt_sampler Optimization sampler
#' @param sample Represents a sample of the target. Allows reuse of the same sample
#'               for construction of confidence intervals, hypothesis tests, etc.
#' @param level Specify the confidence level (optional)
#' @param initial_guess Initial guesses at upper and lower limits (optional)
#'
#' @return List of confidence intervals.
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
#' \link{optimization_intervals} to initialize an interval for the optimization
#' problem, which specify the observed target and its covariance as well as
#' optimization samples and sampler. We need to initialize an optimization_intervals
#' before computing confidence intervals.
#'
#' \link{confidence_interval_iv} for constructing a selective confidence
#' intervals for a specific parameter of the target.
#'
#' \link{confidence_intervals_iv_ar} for constructing selective confidence
#' intervals for AR test statistic.
#'
#' \link{confidence_interval_iv_ar} for constructing a selective
#' confidence interval for a specific parameter of the target for AR
#' test statistic.
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
#' # Compute confidence intervals for TSLS test statistic
#' confidence_intervals_iv(observed_target=s$observed_target,
#'                         target_cov=s$cov_target,
#'                         score_cov=s$cov_target_score,
#'                         opt_sampler=model$sampler,
#'                         sample=s$opt_sample,
#'                         level=0.95)
confidence_intervals_iv <- function(observed_target,
                                    target_cov,
                                    score_cov,
                                    sample_args=NULL,
                                    opt_sampler,
                                    sample,
                                    level=0.9,
                                    initial_guess=NULL) {

  ndraw = nrow(sample)

  opt_sampling_info = list(sets::tuple(opt_sampler, sample, target_cov, score_cov))
  intervals = conditionalInference::optimization_intervals(opt_sampling_info,
                                                           observed_target,
                                                           ndraw)

  for (i in 1:nrow(observed_target)) {
    keep = matrix(0, nrow(observed_target), ncol(observed_target))
    keep[i] = 1
    if (is.null(initial_guess)) {
      confint = confidence_interval_iv(intervals, keep, level=level)
      l = confint$lower[1,1]
      u = confint$upper[1,1]
    } else {
      confint = confidence_interval_iv(intervals, keep, level=level,
                                       guess=initial_guess[i])
      l = confint$lower[1,1]
      u = confint$upper[1,1]
    }

    if (i == 1) {
      limits = list(sets::tuple(l, u))
    } else {
      #limits = append(limits, sets::tuple(l, u))
      limits[[i]] = sets::tuple(l, u)
    }
  }

  return(limits)
}


#' Confidence Interval
#'
#' Construct a selective confidence interval for a specific observed target
#' statistic.
#'
#' @param intervals Optimization intervals
#' @param linear_func Linear function
#' @param level Specify the confidence level (optional)
#' @param how_many_sd Specify the number of standard deviations (optional)
#' @param guess Initial guess at upper and lower limit (optional)
#'
#' @return A confidence interval object, which is a list of:
#' \item{upper}{Upper limit of the confidence interval}
#' \item{lower}{Lower limit of the confidence interval}
#' @export
#'
#' @importFrom stats uniroot
#'
#' @references
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
#' Plausibly Valid Instruments with Application to Mendelian Randomization.
#'
#' Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
#' Effects After Testing Instrument Strength in Linear Models.
#'
#' @seealso
#' \link{optimization_intervals} to initialize an interval for the optimization
#' problem, which specify the observed target and its covariance as well as
#' optimization samples and sampler. We need to initialize an optimization_intervals
#' before computing confidence intervals.
#'
#' \link{confidence_intervals_iv} for constructing selective confidence
#' intervals for TSLS test statistic.
#'
#' \link{confidence_intervals_iv_ar} for constructing selective confidence
#' intervals for AR test statistic.
#'
#' \link{confidence_interval_iv_ar} for constructing a selective
#' confidence interval for a specific parameter of the target for AR
#' test statistic.
#'
#' \link{pivot_iv} for calculating the pivot quantities/p-values of the
#' optimization intervals for TSLS test statistic.
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
#' # Construct the optimization intervals before computing confidence intervals
#' # for TSLS test statistic
#' ndraw = nrow(s$opt_sample)
#' opt_sampling_info = list(sets::tuple(model$sampler,
#'                                      s$opt_sample,
#'                                      s$cov_target,
#'                                      s$cov_target_score))
#' intervals <- optimization_intervals(opt_sampling_info,
#'                                     s$observed_target,
#'                                     ndraw)
#'
#' # Compute the confidence interval for a specific observed target
#' keep = matrix(0, nrow(s$observed_target), ncol(s$observed_target))
#' keep[1] = 1
#' confidence_interval_iv(intervals, keep, level=0.95)
#' confidence_interval_iv(intervals, keep, level=0.95,
#'                        guess=c(0.5, 1.5))
confidence_interval_iv <- function(intervals, # self
                                   linear_func,
                                   level=0.90,
                                   how_many_sd=20,
                                   guess=NULL) {

  sample_stat = intervals$normal_sample %*% linear_func
  observed_stat = intervals$observed %*% linear_func

  rootU <- function(gamma) {
    # gamma = matrix(gamma, nrow = 1, ncol = length(gamma))
    # candidate = observed_stat + gamma
    # candidate = apply(gamma, c(1,2), function(x) x + observed_stat)
    candidate = gamma + observed_stat[1,1]
    #candidate=1.0165933993484308
    result = conditionalInference::pivot_iv(intervals,
                                            linear_func,
                                            candidate,
                                            alternative='less') - (1 - level) / 2
    return(result)
  }

  rootL <- function(gamma) {
    # gamma = matrix(gamma, nrow = 1, ncol = length(gamma))
    # candidate = observed_stat + gamma
    # candidate = apply(gamma, c(1,2), function(x) x + observed_stat)
    candidate = gamma + observed_stat[1,1]
    result = conditionalInference::pivot_iv(intervals,
                                            linear_func,
                                            candidate,
                                            alternative='less') - (1 + level) / 2
    return(result)
  }

  if (is.null(guess)) {
    grid_min = -how_many_sd * sd(sample_stat)
    grid_max = how_many_sd * sd(sample_stat)
    upper = stats::uniroot(rootU, lower=grid_min, upper=grid_max)$root
    lower = stats::uniroot(rootL, lower=grid_min, upper=upper)$root
  } else {
    delta = 0.5 * (guess[2] - guess[1])

    # find interval bracketing upper solution
    count = 0
    while(TRUE) {
      Lu = guess[2] - delta
      Uu = guess[2] + delta
      valU = rootU(Uu)
      valL = rootU(Lu)
      if (valU * valL < 0) {
        break
      }
      delta = delta * 2
      count = count + 1
    }
    upper = stats::uniroot(rootU, lower=Lu, upper=Uu)$root

    # find interval bracketing lower solution
    count = 0
    while (TRUE) {
      Ll = guess[1] - delta
      Ul = guess[1] + delta
      valU = rootL(Ul)
      valL = rootL(Ll)
      if (valU * valL < 0) {
        break
      }
      delta = delta * 2
      count = count + 1
    }
    lower = stats::uniroot(rootL, lower=Ll, upper=Ul)$root
  }

  upper = apply(observed_stat, c(1,2), function(x) x + upper)
  lower = apply(observed_stat, c(1,2), function(x) x + lower)

  return_list <- list("upper" = upper,
                      "lower" = lower)
  return(return_list)
}


