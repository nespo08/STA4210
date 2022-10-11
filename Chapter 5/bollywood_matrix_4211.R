bolly1 <- read.csv(
    "http://www.stat.ufl.edu/~winner/data/bollywood_boxoffice2.csv")
attach(bolly1); names(bolly1)

Y1 <- Gross; Y2 <- log(Gross)
X1 <- Budget; X2 <- log(Budget)
n <- length(Y2)

## Fitting model based on logs of Y and X
## Regression Coefficients, Standard errors, t-tests/CIs

Y <- matrix(Y2, ncol=1)                ## Y is nx1 matrix of log(Gross)
X <- matrix(c(rep(1,n),X2), ncol=2)    ## X is nx2 matrix [1_n, X2]
XprimeX <- t(X) %*% X                  ## X'X is 2x2 matrix
XprimeY <- t(X) %*% Y                  ## X'Y is 2x1 vector
invXX <- solve(XprimeX)                ## (X'X)^(-1) is 2x2 inverse of X'X
b       <- invXX %*% XprimeY           ## b is 2x1 - OLS estimate of beta
YprimeY <- t(Y) %*% Y                  ## Y'Y = sum(Y_i^2)  
s2 <- (YprimeY - 
     t(b) %*% XprimeY) / (n-2)         ## s^2 is Unbiased estimate of sigma^2
V_b <- s2[1,1] * invXX                 ## Variance-Covariance matrix of b
SE_b <- sqrt(diag(V_b ))               ## SE{b0}, SE{b1}
t_b <- b / SE_b                        ## t-stats for beta0, beta1
p_t_b <- 2*(1-pt(abs(t_b),n-2))        ## P-values for t-stats
beta_LB <- b + qt(.025,n-2)*SE_b       ## CI Lower Bounds for beta
beta_UB <- b + qt(.975,n-2)*SE_b       ## CI Upper Bounds for beta

## Matrix output
matrix.out <- cbind(XprimeX, XprimeY, invXX, 
       b, V_b, SE_b)
colnames(matrix.out) <- c("X'X", "X'X", "X'Y", "inv(X'X)", "inv(X'X)",
         "b", "V{b}", "V{b}", "SE{b}")
round(matrix.out,3)

## Set up tables for output
## Regression coefficients/t-ests,CI's
coeff.out <- cbind(b, SE_b, t_b, p_t_b,
                   beta_LB, beta_UB)
colnames(coeff.out) <- c("Estimate", "Std Err", "t", "2P(>|t|)", "LB", "UB")
rownames(coeff.out) <- c("Intercept", "Log Budget")
round(coeff.out, 4)


## Analysis of Variance
H <- X %*% invXX %*% t(X)              ## H is nxn projection (hat) matrix
I <- diag(n)                           ## I is nxn identity matrix
J_1_n <- matrix(rep(1/n, n^2), ncol=n) ## J_1_n is nxn matrix of 1/n

Yhat <- H %*% Y                        ## nx1 vector of predicted values
e <- (I - H) %*% Y                     ## nx1 vector of residuals
Ybar <- J_1_n %*% Y                    ## nx1 vector of ybar

SSTO <- t(Y) %*% (I - J_1_n) %*% Y    ## Total (corrected) sum of squares
SSE <- t(Y) %*% (I - H) %*% Y         ## Error (Residual) sum of squares
SSR <- t(Y) %*% (H - J_1_n) %*% Y     ## Regression sum of squares
df_T <- n - 1                         ## Total (corrected) degrees of freedom
df_E <- n - ncol(X)                   ## Error degrees of freedom
df_R <- ncol(X) - 1                   ## Regression degrees of freedom
MSE <- SSE / df_E                     ## Error Mean Square
MSR <- SSR / df_R                     ## Regression Mean Square
F_obs <- MSR / MSE                    ## F-statistic
F_05 <- qf(.95, df_R, df_E)           ## Critical F value
p_F <- 1 - pf(F_obs, df_R, df_E)      ## P-value

## Set-up output of Anaysis of Variance Table
df.out <- rbind(df_R,df_E,df_T)
SS.out <- rbind(SSR,SSE,SSTO)
MS.out <- rbind(MSR,MSE,NA)
F.out <- rbind(F_obs,NA,NA)
F05.out <- rbind(F_05,NA,NA)
p.out <- rbind(p_F,NA,NA)
aovm.out <- cbind(df.out,SS.out,MS.out,F.out,F05.out,p.out)
colnames(aovm.out) <- c("df","SS","MS","F","F(.05)","P(>F)")
rownames(aovm.out) <- c("Regression", "Error", "Total")
round(aovm.out, 3)


## CI for Mean @ X*, PI for Individual

X_h <- matrix(c(1,log(100)), ncol=1)   ## 2x1 vector X_h with X=log(100)
Yhat_h <- t(X_h) %*% b                 ## Estimated mean/Predicted value
V_Yhat_h <- 
  t(X_h) %*% V_b %*% X_h               ## Variance of estimated mean
V_Ynew_h <- s2 + V_Yhat_h              ## Variance of prediction error
Yhat_h_mean_LB <-   Yhat_h + 
    qt(.025,df_E)*sqrt(V_Yhat_h)       ## Lower Bound for Mean
Yhat_h_mean_UB <-   Yhat_h + 
    qt(.975,df_E)*sqrt(V_Yhat_h)       ## Upper Bound for Mean
Yhat_h_indv_LB <-   Yhat_h + 
    qt(.025,df_E)*sqrt(V_Ynew_h)       ## Lower Bound for Individual
Yhat_h_indv_UB <-   Yhat_h + 
    qt(.975,df_E)*sqrt(V_Ynew_h)       ## Upper Bound for Individual


## Output Estimated Mean/Predictions of Individual Values
est_log.out <- cbind(X_h[2,], Yhat_h, sqrt(V_Yhat_h), sqrt(V_Ynew_h),
      Yhat_h_mean_LB, Yhat_h_mean_UB, Yhat_h_indv_LB, Yhat_h_indv_UB)
colnames(est_log.out) <- c("X_h", "Estimate", "SE mean", "SE indiv", 
                   "LB mean", "UB mean", "LB indiv", "UB indiv")
rownames(est_log.out) <- c("log(Gross)")
round(est_log.out, 3)

est.out <- cbind(exp(X_h[2,]), exp(Yhat_h), exp(Yhat_h_mean_LB),  
   exp(Yhat_h_mean_UB), exp(Yhat_h_indv_LB), exp(Yhat_h_indv_UB))
colnames(est.out) <- c("X_h", "Estimate",  
                     "LB mean", "UB mean", "LB indiv", "UB indiv")
rownames(est.out) <- c("Gross")
round(est.out, 3)


### Using lm function

bolly.mod1 <- lm(Y2 ~ X2)
summary(bolly.mod1)
anova(bolly.mod1)
confint(bolly.mod1)
(ci.mean <- predict(bolly.mod1, list(X2=log(100)), int="c"))
(pi.indiv <- predict(bolly.mod1, list(X2=log(100)), int="p"))
exp(ci.mean); exp(pi.indiv)





