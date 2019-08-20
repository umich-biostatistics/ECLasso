#' Solve Generalized Linear Models via Alternating Direction Method of Multipliers (ADMM).
#'
#' This function uses ADMM to approximate parameters in Gaussian, Binomial, Poisson, and Cox models.
#'
#' @param x The design matrix of the data.
#' @param y The response vector; must contain the same amount of values as rows of x.
#' @param lambda_seq User-specified sequence of lambdas the algorithm forced to iterate over, default NULL.
#' @param family The type of model the data comes from; gaussian, binomial, poisson, or Cox.
#' @param rho Weight given to equality constraint vs. loss function, default is 1 aka equal weighting.
#' @param intercept Determines whether to add an intercept to the design matrix.
#' @param equality Forces the parameters to sum to zero if TRUE.
#' @param censor A 0/1 vector for Cox family only, where 1 is exact and 0 is right-censored; default is NULL.
#' @param C Part of equality constraints, set to 1 if equality is TRUE, 0 otherwise.
#' @param d Part of equality constraints, always set to 0 unless user-defined.
#' @param inexact If true, solves a quadratic approximation of the objective function at each iteration; required to be TRUE for "cox" family.
#' @param tol Determines stopping criteria for the algorithm.
#' @param maxit Half the max amount of time each iteration allowed to run.
#' @param nlambda The number of lambdas the function iterates over, should be a multiple of 50.
#' @param lambda.min.ratio Defines what fraction of the maximum lambda the minimum is.
#' @param min.abs Any covariate effect with smaller absolute value than this is set to zero during iteration.
#' @param lower.limits Lower boundary of parameter, any iteration with term(s) below this is reset to this value.
#' @param upper.limits Upper boundary of parameter, any iteration with term(s) above this is reset to this value.
#' @param penalty.factor Penalty factor for shrinkage function, involved in z update.
#' @param toplot Determines if a plot of the solution path should be outputted
#' @param ParamCutoff Specify a cutoff above which predictors are returned at completion, if NULL by default, returns the 6 predictors with greatest magnitude
#'
#' @keywords ADMM
#'
#' @return Solution: parameter estimates at each iteration,
#'     time: total elapsed program time,
#'     Iter: Time taken for each iteration (consider increasing maxit if this is at maximum),
#'     Lambda: the lambda for each iteration,
#'     Intercept: intercept value if specified TRUE.
#'
#' @examples
#' c0=0.1
#' beta1=1
#' beta2=-1
#' beta3=0.5
#' beta4=-0.5
#' betazeroes<-rep(0,6)
#' betavector<-c(beta1,beta2,beta3,beta4,t(betazeroes))
#' n=200
#' xdata<-matrix(rnorm(n*length(betavector)),ncol=length(betavector))
#' times=rexp(n,c0*exp(rowSums(t(t(xdata)*betavector))))
#' time.censor=rexp(n,c0*exp(beta1*xdata[,1]))
#' censorv=ifelse(times<time.censor, 1, 0)
#' time <- ifelse(times<time.censor, times, time.censor)
#' ADMM(x=xdata,y=time,family="cox",intercept=FALSE,equality=TRUE,censor=censorv,inexact=TRUE)
#'
#' @export


ADMM <-
  function(x, y, lambda_seq = NULL, family = c("gaussian", "binomial", "poisson","cox"),
           rho = 1, intercept = TRUE, equality = FALSE, censor=NULL, C = NULL, d = NULL, inexact = FALSE,
           tol = 1e-4, maxit = 1000, nlambda = 100, lambda.min.ratio = ifelse(nobs < nvars, 0.01, 0.0001),
           min.abs = 0, lower.limits = -Inf, upper.limits = Inf,
           penalty.factor = 1, toplot = TRUE, ParamCutoff = NULL, exclude = rep(FALSE, nvars)) {
    #note that using the cox family requires intercept=FALSE and inexact=TRUE
    #start the clock
    ptm <- proc.time()

    family = match.arg(family)

    abs.tol = ifelse(inexact, tol * 1e-4, tol * 1e-2)
    rel.tol = ifelse(inexact, tol * 1e-2, tol)
    max.iter = ifelse(inexact, maxit * 2, maxit)

    A = as.matrix(x)
    b = y

    nobs = nrow(A);
    nvars = ncol(A);

    if (is.null(d)) d = 0;

    if (equality == TRUE) {
      if (is.null(C)) {
        #C is not specified, set it to a row of 1.
        C = matrix(1, nrow = 1, ncol = nvars)
      }
    } else {
      # no equality constraints
      C = matrix(0, nrow = 1, ncol = nvars)
      d = 0
    }

    if (length(min.abs) == 1) min.abs = rep(min.abs, nvars)
    if (length(lower.limits) == 1) lower.limits = rep(lower.limits, nvars)
    if (length(upper.limits) == 1) upper.limits = rep(upper.limits, nvars)
    if (length(penalty.factor) == 1) penalty.factor = rep(penalty.factor, nvars)
    if (!is.logical(exclude)) {
      exc = rep(F, nvars)
      exc[exclude] = T
      exclude = exc
    }

    if (intercept == TRUE) {
      A = cbind(1, A);
      C = cbind(0, C);
      min.abs = c(0, min.abs)
      lower.limits = c(-Inf, lower.limits)
      upper.limits = c(Inf, upper.limits)
      penalty.factor = c(0, penalty.factor)
      exclude = c(F, exclude)
    }

    p = ncol(A);
    #initial value of x and z
    x = matrix(0, nrow = p, ncol = 1);
    z = matrix(0, nrow = p, ncol = 1);
    u1 = 0;
    u2 = matrix(0, nrow = p, ncol = 1);

    #B=rbind(C,diag(p));

    if (family == "binomial") {
      F = -diag(b[, 1]) %*% A;
      #pre-calculate the -b_i%*%a_i in logistic loss
      if (is.null(lambda_seq)) {
        b_hat = matrix(0, nobs);
        ratio = sum(b == 1) / nobs;
        for (i in 1:nobs) {
          if (b[i] == 1) b_hat[i] = 1 - ratio else
            b_hat[i] = -ratio
        }
        lambda_max = norm(t(A) %*% b_hat, 'i');
      }
      if (inexact == TRUE) {
        J = t(F) %*% diag((exp(F %*% x) / (1 + exp(F %*% x)) ^ 2)[, 1]) %*% F;
      }
    } else if (family == "poisson") {
      if (is.null(lambda_seq)) {
        lambda_max = norm(t(A) %*% (b - exp(A %*% x)), 'i');
      }
      if (inexact == TRUE) {
        J = t(A) %*% diag(exp(A %*% x)[, 1]) %*% A;
      }
    } else if(family=="cox"){

      if (inexact == TRUE){
        #order the data by decreasing order of survival time
        ord = order(b, decreasing =T)
        b1 = b[ord]
        censor1 = censor[ord]
        A1 = as.matrix(A[ord,])
        #calculating terms for observed information, to determine lambda max start point
        theta = exp(A1 %*% x)
        numerator = apply(A1 * as.vector(theta), 2, cumsum)
        H = matrix(0, nvars, nvars)
        temp = matrix(0, nvars, nvars)
        theta1 = as.vector(cumsum(theta))
        for (i in 1:nobs) {
          temp = temp + A1[i,] %*% t(A1[i,] * theta1[i]) ##j
          if (censor1[i] == 1) {
            H = H - temp / theta1[i] - numerator[i, ] %*% t(numerator[i,])/theta1[i]^2
          }
        }
        J=-H
      }
      if (is.null(lambda_seq)) {
        w0 = revcumsum(((0:(nobs-1))/(1:nobs)^2))
        z0 = (censor1-revcumsum(1/(1:nobs)))/w0
        lambda_max = max(colSums(w0*A1*z0))
      }
    }else {
      Atb = t(A) %*% b;
      if (is.null(lambda_seq)) {
        lambda_max = norm(Atb, 'i');
      }
      #cache the factorization
      if (inexact == TRUE) {
        J = t(A) %*% A;
      } else {
        E = rbind(A, sqrt(rho) * C);
        LU = Choleski_factors(E, rho);
        L = LU[[1]];
        #m*m lower triangular
        U = LU[[2]];
      }
    }

    if (is.null(lambda_seq)) {
      lambda_min = lambda.min.ratio * lambda_max;
      lambda_seq = exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
    }

    if (inexact == TRUE) {
      eigenvalues = eigen(J)$values;
      h = max(eigenvalues);
      P = 1 / (h + rho) * (diag(p) - rho / (h + rho * (p + 1)) * t(C) %*% C);
    }

    iter_num = c()
    solution = c()

    if(family == "binomial"|family == "poisson"|family == "cox"){
      if(family == "binomial"){
        grad = t(F) %*% matrix(rep(0.5,nobs))
      } else if(family == "poisson"){
        grad = -(t(A) %*% b) + t(A) %*% matrix(rep(1,nobs))
      } else if(family == "cox"){
        numeratorIT = apply(as.matrix(A1), 2, cumsum)
        denominatorIT = cumsum(matrix(rep(1,nobs)))
        w = colSums((A1 - numeratorIT / denominatorIT) * censor1)
        grad = -w
      }

      u1 = 1e-2
      u2 = t(C) * (d - u1) - grad/rho

      iter_num = c(iter_num, 1)
      solution = cbind(solution, matrix(0, nrow = p, ncol = 1))
    }

    for (i in 2:length(lambda_seq)) {
      lambda = lambda_seq[i];

      for (k in 1:max.iter) {
        delta = rho * (z - u2 + t(C) %*% (d - u1))

        if (family == "gaussian") {
          if (inexact == TRUE) {
            if (p > 50 * nobs) {
              x = P %*% (delta + h * x - t(A) %*% (A %*% x) + Atb)
            } else {
              x = P %*% (delta + h * x - J %*% x + Atb)
            }
          } else {
            q = Atb + delta
            #nvars*1 vector
            if (nobs >= p) {
              x = backsolve(U, forwardsolve(L, q))
            } else {
              x = q / rho - (t(E) %*% (backsolve(U, forwardsolve(L, (E %*% q))))) / rho ^ 2
            }
          }

        } else if (family=="cox"){
          thetaIT = exp(A1 %*% x)
          numeratorIT = apply(as.matrix(A1) * as.vector(thetaIT), 2, cumsum)
          denominatorIT = cumsum(thetaIT)
          w = colSums((A1 - numeratorIT / denominatorIT) * censor1)
          x= P %*% (h * x+ delta +w)
        }else {
          if (inexact == TRUE) {
            if (family == "binomial") {
              x = P %*% (h * x + delta - t(F) %*% (exp(F %*% x) / (1 + exp(F %*% x))))
            } else {
              x = P %*% (h * x + delta + t(A) %*% b - t(A) %*% exp(A %*% x))
            }
          } else {
            x = NR_update_x(A, F, b, x, u1, u2, z, rho, C, d, family)
          }
        }

        x[abs(x) < min.abs] = 0
        x[x < lower.limits] = lower.limits
        x[x > upper.limits] = upper.limits
        x[exclude] = 0

        # z-update
        zold = z;
        z = x + u2;

        z = shrinkage(z, lambda * penalty.factor / rho)

        #u-update
        Cx = as.numeric(C %*% x);
        u1 = u1 + Cx - d;
        u2 = u2 + (x - z);

        #diagnostics, reporting, termination checks
        r_norm = vector.2.norm(rbind(Cx - d, x - z));
        #r_norm
        s_norm = vector.2.norm( - rho * (z - zold));
        #s_norm
        eps_pri = sqrt(p + 1) * abs.tol + rel.tol * max(vector.2.norm(rbind(Cx, x)), vector.2.norm( - z), abs(d));
        #eps_pri
        eps_dual = sqrt(p) * abs.tol + rel.tol * vector.2.norm(rho * (u1 + u2));
        #eps_dual

        #termination criterion
        if (r_norm < eps_pri && s_norm < eps_dual) {
          break
        }
      }
      iter_num = c(iter_num, k)

      solution = cbind(solution, t(t(x)))
    }
    if(toplot == TRUE){
      solnNorm<-(abs(solution[,nlambda])-min(abs(solution[,nlambda])))/(max(abs(solution[,nlambda]))-min(abs(solution[,nlambda])))
      plotlambdas<-plot(log(lambda_seq),solution[1,],type="l",xlim=c(log(lambda_seq[nlambda]),log(lambda_seq[1])),ylim=c(-1.2*max(abs(solution[,nlambda])),1.2*max(abs(solution[,nlambda]))),xlab="Log lambda",ylab="Param. Estimate",col=addalpha(rep("red4", 100), solnNorm[1]))
      for (i in 2:nvars){
        lines(log(lambda_seq),solution[i,],col=addalpha(rep("red4", 100), solnNorm[i]))
      }
      abline(v=c(log(lambda_seq[nlambda*0.1]),log(lambda_seq[nlambda*0.5]),log(lambda_seq[nlambda*0.9])),col="blue")
    }
    if(is.null(ParamCutoff) == TRUE){
      highindices<-head(order(abs(solution[,nlambda]),decreasing = TRUE))
    } else {
      highindices<-which(abs(solution[,nlambda]) > ParamCutoff)[order(abs(solution[,nlambda][which(abs(solution[,nlambda]) > ParamCutoff)]),decreasing = TRUE)]
    }
    highparams<-solution[highindices,nlambda]
    time = proc.time() - ptm
    return(list(solution = solution, time = time, iter = iter_num, lambda = lambda_seq,intercept = intercept,highindices=highindices,highparams=highparams))
  }
