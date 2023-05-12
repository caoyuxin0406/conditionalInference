# conditionalInference

conditionalInference is an R package for inferring the treatment effect with instrument variables after testing for instrument strength in linear models. The package includes a method to generate the dataset  for testing instrument strength and inferring treatment effect. It can also take the user-specified dataset to perform these tasks. The package contains methods to construct and fit the model, provide summary statistics, and estimate variance-covariance matrix. The package contains methods for TSLS statistic and the AR test. All methods in the software are outlined in Bi, Kang, and Taylor (2020).

To install this package in R from GitHub, run the following commands:

```
install.packages("devtools")
library(devtools) 
install_github("caoyuxin0406/conditionalInference")
```

## References
Bi, Nan, Hyunseung Kang, and Jonathan Taylor. "Inference after selecting plausibly valid instruments with
application to mendelian randomization." arXiv preprint arXiv:1911.03985 (2019).

Bi, Nan, Hyunseung Kang, and Jonathan Taylor. "Inferring treatment effects after testing instrument strength
in linear models." arXiv preprint arXiv:2003.06723 (2020).

Card, David. "Using geographic variation in college proximity to estimate the return to schooling." Working
Paper 4483, National Bureau of Economic Research. (1993).

Kang, Hyunseung, Yang Jiang, Qingyuan Zhao, and Dylan S. Small. "ivmodel: An R Package for Inference and Sensitivity Analysis of Instrumental Variables Models with One Endogenous Variable." Observational Studies 7, no. 2 (2021): 1-24. doi:10.1353/obs.2021.0029.
