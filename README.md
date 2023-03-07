## NIBAMM

### Overview

NIBAMM is a Bayesian framework integrating gene expression, DNA methylation and pathway information to identify significant biological features associated with disease outcomes. Our data includes gene expression, DNA methylation, a map from DNA methylation probes to genes, structure dependencies among genes and clinical covariates and response.

### User manual

The required R packages:

* Rcpp
* RcppArmadillo
* RcppGSL
* scales
* caret
* MASS

After constructing the list `data` and list `prior` as shown in *example.R*, you can call `NIBAMM` function by

```R
> source("R/NIBAMM.R")
> result <- NIBAMM(data, prior, max_iters = 10000)
```

The full MCMC chain is stored in `result`, which can be used in the following analysis and plots. For example, we can find PPI of variables by

```R
> n.burnin = 2000
> max_iters = 10000
> colMeans(result$gamma[n.burnin:max_iters, ])
```

If the response is patients' survival time, the log transformation with natural base is required firstly. Then the `data` can be constructed as shown in *example.R*. In this case, run MCMC by

```R
> source("R/survivalNIBAMM.R")
> result <- survivalNIBAMM(data, prior, max_iters = 10000)
```

If there is any bugs without being fixed timely, please contact author directly zhubencong@gmail.com