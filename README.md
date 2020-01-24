# ECLasso

This is the R package described in "ECLasso: Fitting equality-constrained, L1-penalized models with inexact ADMM to find gene pairs" by Lam Tran, Lan Luo, and Hui Jiang. Briefly, this package uses the inexact version (Wang et al. 2014) of the alternating direction method of multipliers introduced by Boyd et al. (2011) to fit a lasso-penalized model to provided data. Provided that the sum of the estimated parameters is constrained to equal 0, two parameters will remain for sufficiently large weight given to the lasso penalty, forming the "solution pair". As this lasso penalty weight decreases, more terms will enter, such that more than one "pair" is produced.

Although our intended application in the usage of ECLasso is to find pairs of genes associated with some outcome (i.e. decreased survival), ECLasso can be used in any applications where pairs of solutions, rather than single solutions, may be more useful. A more rigorous treatment of the math behind ECLasso as well as step-by-step descriptions of the algorithms in the package are included in the supplementary material of the paper (available at:). 

## Features
The novelty of ECLasso lies in its ability to, when the sum of the estimated parameters is constrained to equal 0, find a solution pair to supported types of regression models other than linear regression.  Some of the features of ECLasso are described below:
* Ability to fit models to multiple types of data response variables: continuous (linear regression), binary (logistic regression), count (Poisson regression), and survival (Cox regression)
* Option to perform inexact updates to the estimation step, greatly speeding up the linear case (a requirement for the logistic, Poisson, and Cox cases)  
* Automatically define a sequence of lasso penalty weights or allow a user-defined sequence of weights 
* Allow the user to specify the maximum iterations per weight, trading precision at the cost of run time
* Set constraints on the solution vector (or none at all), the stopping criteria, and whether to include an intercept term or a plot of the solution path

## Getting started
A simple example using ECLasso with survival data is described below:
1. Install the package from GitHub and load/attach it to R: 
````
install.packages("devtools")
library(devtools)
install_github("umich-biostatistics/ECLasso")
library(ECLasso)
````
2. Simulate some survival data:
`````````````
c0 = 0.1
beta1 = 1  ##true values of significant predictors
beta2 = -1
beta3 = 0.5
beta4 = -0.5
betazeroes = rep(0, 6) ##6 noise predictors 
betavector = c(beta1, beta2, beta3, beta4, t(betazeroes))
n = 50  ##sample size
xdata = matrix(rnorm(n * length(betavector)),ncol = length(betavector))  ##simulated values for covariates
times = rexp(n, c0 * exp(rowSums(t(t(xdata) * betavector))))  ##time to death
time.censor = rexp(n, c0 * exp(beta1 * xdata[, 1])) ##time to censoring
censorv = ifelse(times < time.censor, 1, 0) ##equal to 1 if time to death occurs before time to censoring 
time = ifelse(times < time.censor, times, time.censor) ##time is the lesser of death and censoring times
`````````````
3. Run the main ECLasso.fit function:
`
ECLasso.fit(x = xdata, y = time, family = "cox", intercept = FALSE, equality = TRUE, censor = censorv, inexact = TRUE)
`

## Understanding ECLasso's output
Upon completion, ECLasso outputs a list with 6 items by default: "time", "iter", "lambda", "intercept", "highindices", and "highparams".

* Time: Displays the time (in seconds) it took for ECLasso to complete. 

* Iter: At each lambda, displays the number of ADMM iterations it took to converge or the maximum iterations, whichever is smaller. Consider increasing the maximum number of iterations if by the end of the sequence, the default maximum of 2000 iterations is still used, in order to ensure proper convergence. 

* Lambda: The value of lambda along the sequence; the weight associated with the lasso penalty. 

* Intercept: If specified, returns the intercept of the fitted model.

* Highindicies: The indices of the predictors with the largest absolute parameter value upon ECLasso completion. Note that choosing predictors this way is only an approximation for selection based on the order predictors enter the regression model. 

* Highparams: The estimated values of the parameters associated with the predictors indexed by "highindicies".

ECLasso outputs a plot of the solution path of the data, proceeding from right to left. Of particular note is the point when the parameter estimates first become nonzero; assuming enough iteration steps were performed, we are guaranteed to have 2 predictors enter the model as a pair. For further details about interpretation, please consult the supplementary material. 
## Troubleshooting
Please direct questions to Lam Tran at lamtran@umich.edu


