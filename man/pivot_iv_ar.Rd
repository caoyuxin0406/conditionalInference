% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intervals_ar.R
\name{pivot_iv_ar}
\alias{pivot_iv_ar}
\title{Pivot with Instrument Variables for AR Test Statistic}
\usage{
pivot_iv_ar(
  opt_intervals,
  linear_func,
  parameter_reference,
  candidate,
  gl_ar,
  alternative = "twosided"
)
}
\arguments{
\item{opt_intervals}{Optimization intervals}

\item{linear_func}{Linear function}

\item{parameter_reference}{Parameter for reference}

\item{candidate}{Candidate value}

\item{gl_ar}{The group lasso object for AR test statistic representing the
penalized convex optimization equation (generated by
group_lasso_iv function)}

\item{alternative}{A list of ['greater', 'less', 'twosided'],
which indicates what alternative to use}
}
\value{
A pivot quantity.
}
\description{
Calculate the pivot quantities of the optimization intervals.
}
\examples{
# Fit the group lasso optimization model for AR statisitc
Z = matrix(rnorm(10*3), nrow = 10, ncol = 3) +
    matrix(replicate(3, matrix(0, nrow = 10, ncol = 1)),
    nrow = 10)
errorTerm = MASS::mvrnorm(n=10, mu=rep(0,2),
            Sigma=rbind(c(1, 0.8),c(0.8, 1)))
D = Z \%*\% rep(1, 3) + errorTerm[1:10, 2]
Y = D * 1 + errorTerm[1:10, 1]
gl_ar <- group_lasso_iv_ar(Y, D, Z)
model_ar <- fit_ar(gl_ar)

# Estimate covariance and obtain summary statistics
cov = estimate_covariance(Y, D, Z)
s <- summary_ar(gl_ar, model_ar$sampler,
                Sigma_11=cov[1,1], Sigma_12=cov[1,2])

# Construct the optimization intervals before computing confidence intervals
# for AR test statistic
opt_sampling_info = list(sets::tuple(model_ar$sampler,
                                     s$opt_sample,
                                     s$cov_target,
                                     s$cov_target_score))
intervals <- optimization_intervals(opt_sampling_info,
                                    s$observed_target,
                                    nrow(s$opt_sample))

# Compute pivot value for AR test statistic
pivot_iv_ar(opt_intervals=intervals,
            linear_func=s$K2,
            parameter_reference=0,
            candidate=0,
            gl_ar=gl_ar,
            alternative='greater')

# Compute the pivot value for a sequence of candidate value
level = 0.95
grid_min = -4000 * 0.5946693  # -how_many_sd * sd(sample_stat)
grid_max = 4000 * 0.5946693 # how_many_sd * sd(sample_stat)

for (beta in seq(grid_min, grid_max, length.out=10)) {
  pivot_iv_ar(opt_intervals=intervals,
              linear_func=s$K2,
              parameter_reference=s$two_stage_ls[1,1],
              candidate=beta,
              gl_ar=gl_ar,
              alternative='greater') - (1-level) + (1-level)
}
}
\references{
Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2019). Inference After Selecting
Plausibly Valid Instruments with Application to Mendelian Randomization.

Bi, Nan & Kang, Hyunseung & Taylor, Jonathan. (2020). Inferring Treatment
Effects After Testing Instrument Strength in Linear Models.
}
\seealso{
\link{pivot_iv} for calculating the pivot quantities/p-values of the
optimization intervals for TSLS test statistic. This process is similar
to that of AR test statistic.

\link{weights_iv_ar} that computes the weights of candidate pivot values
for AR test statistic.

\link{test_stat} that returns the result of comparing parameter of
interest and observed target test statistic.
}
