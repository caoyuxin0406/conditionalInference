## Codes

library(conditionalInference)

### Section 2.3

set.seed(0)
# True beta value is 1
b <- bigaussian_instance(n=1000, p=10, gsnr=1, beta=1,
                         Sigma = rbind(c(1, 0.8), c(0.8, 1)))
# True beta value is 0
b1 <- bigaussian_instance(n=1000, p=10, gsnr=1, beta=0,
                          Sigma = rbind(c(1, 0.8), c(0.8, 1)))

### Section 3.1.2

set.seed(0)
Y <- b$Y
D <- b$D
Z <- b$Z
pre_test(Y, D, Z)
pre_test(Y, D, Z, C0=1500)

### Section 3.2.2

set.seed(0)
gl <- group_lasso_iv(Y, D, Z, C0=10)
gl1 <- group_lasso_iv(Y, D, Z, penalty=3,
                      randomizer_scale=0.8, C0=5)
Y1 <- b1$Y
D1 <- b1$D
Z1 <- b1$Z
gl2 <- group_lasso_iv(Y1, D1, Z1, C0=10)
gl$data_part
gl$penalty

### Section 3.3.2

set.seed(0)
model <- fit_tsls(gl)
model1 <- fit_tsls(gl1)
model$initial_soln
model$observed_opt_state
model$cond_mean
model$cond_cov

initial_omega2 = t(matrix(c(15.10398, 3.426183, 8.380046, 19.18674,
                            15.99021, -8.367544, 8.134745, -1.295934, -0.8837694, 3.515583)))
model2 <- fit_tsls(gl, perturb = initial_omega2)
model2$observed_opt_state
model2$cond_mean
model2$cond_cov

set.seed(0)
cov <- estimate_covariance(gl$Y, gl$D, gl$Z)
cov1 <- estimate_covariance(gl1$Y, gl1$D, gl1$Z)
cov

### Section 3.4.3

set.seed(0)
s <- summary_tsls(gl, model$sampler, Sigma_11 = cov[1,1], Sigma_12 = cov[1,2])
s$observed_target
matrix(s$opt_sample[1:10,1])
s$pivots[[1]]
s$pvalues[[1]]
s$intervals[[1]]

set.seed(0)
s1 <- summary_tsls(gl1, model1$sampler, Sigma_11 = cov1[1,1], Sigma_12 = cov1[1,2])
s1$observed_target
s1$pivots[[1]]
s1$pvalues[[1]]
s1$intervals[[1]]

set.seed(0)
s2 <- summary_tsls(gl, model$sampler, parameter = 1,
                   Sigma_11 = cov[1,1], Sigma_12 = cov[1,2])
s2$observed_target
s2$pivots[[1]]
s2$pvalues[[1]]
s2$intervals[[1]]

### Section 4.2

set.seed(0)
gl_ar <- group_lasso_iv_ar(Y, D, Z)
model_ar <- fit_ar(gl_ar)

### Section 4.3

set.seed(0)
s1_ar <- summary_ar(gl_ar, model_ar$sampler,
                    Sigma_11 = cov[1,1], Sigma_12 = cov[1,2])
s1_ar$observed_target
matrix(s1_ar$opt_sample[1:10,1])
s1_ar$pivots[[1]]
s1_ar$pvalues[[1]]

### Section 5.2

#### Card Data

library(ivmodel)

# One Instrument Anaylsis with Proximity to 4yr College as IV#
data(card.data)
Y=card.data[,"lwage"]
D=card.data[,"educ"]
Z=card.data[,"nearc4"]
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
        "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
        "reg668", "smsa66")
X=card.data[,Xname]

# One Instrument Anaylsis with Proximity to 4yr College as IV#
iv_model <- ivmodel(Y=Y,D=D,Z=Z,X=X)
iv_model

#### Data Processing

# Compare with ivmodel
Y=matrix(Y)
D=matrix(D)
Z=matrix(Z)

# demean
Y = Y - mean(Y)
D = D - mean(D)
Z  = Z - mean(Z)
col_mu_X = colMeans(X)
for (i in 1:nrow(X)) {
  X[i,1:14]  = X[i,1:14] - col_mu_X
}

# residual out X
X = as.data.frame(X)

glm_YX <- glm(Y ~ ., data = X)
rY <- glm_YX$residuals
glm_DX <- glm(D ~ ., data = X)
rD <- glm_DX$residuals
glm_ZX <- glm(Z ~ ., data = X)
rZ <- glm_ZX$residuals

rY = matrix(rY)
rD = matrix(rD)
rZ = matrix(rZ)

#### Naive Inference

result <- pre_test(rY, rD, rZ); result
cov_card <- estimate_covariance(rY, rD, rZ)
naive_inf <- naive_inference_tsls(rY, rD, rZ, pass_pre_test=result,
                                  Sigma_11=cov_card[1,1], compute_intervals=TRUE)

naive_inf$tsls
naive_inf$std
naive_inf$pval
naive_inf$interval

#### Conditional Inference

set.seed(0)
gl_card <- conditionalInference::group_lasso_iv(rY, rD, rZ, C0=10, randomizer_scale=0.0001, perturb = 0)
model_card <- conditionalInference::fit_tsls(gl_card)
cov_card <- conditionalInference::estimate_covariance(rY, rD, rZ)
s_card <- conditionalInference::summary_tsls(gl_card, model_card$sampler, Sigma_11 = cov_card[1,1], Sigma_12 = cov_card[1,2], ndraw=1000000, burnin=100000)

s_card$observed_target
s_card$cov_target
s_card$cov_target_score
s_card$pivots
s_card$pvalues
s_card$intervals


