# Homework 1


# Question 1 --------------------------------------------------------------

# Read in data

spdata <- read.csv("http://www.stat.ufl.edu/~winner/data/bioequiv_sulf.csv")
attach(spdata); names(spdata)
spdata


# P.1.A -------------------------------------------------------------------

## Create a variable that is AUC for the drug sulfadoxine and a variable for formulation

AUC.sulf <- y[measure==2 & drug==1]       # measure=AUC, drug=sulfadoxine
form21 <- form[measure==2 & drug==1]               # form=1 if Test, form=2 if Ref
form.AUC.sulf <- factor(form21, levels=1:2, labels=c("T","R"))
cbind(form.AUC.sulf, AUC.sulf)
plot(AUC.sulf ~ form.AUC.sulf, main="Sulfadoxine AUC by Formulation")

## Compute n, ybar, sd for Test and Reference

n.T <- length(AUC.sulf[form.AUC.sulf=="T"])
n.R <- length(AUC.sulf[form.AUC.sulf=="R"])
ybarl.T <- mean(AUC.sulf[form.AUC.sulf=="T"])
ybarl.R <- mean(AUC.sulf[form.AUC.sulf=="R"])
sdl.T <- sd(AUC.sulf[form.AUC.sulf=="T"])
sdl.R <- sd(AUC.sulf[form.AUC.sulf=="R"])

# Print results in tabular form
T.stats <- cbind(n.T, ybarl.T, sdl.T)
R.stats <- cbind(n.R, ybarl.R, sdl.R)
TR.stats <- rbind(T.stats, R.stats)
colnames(TR.stats) <- c("n", "mean", "SD")
rownames(TR.stats) <- c("Test", "Reference")
round(TR.stats, 5)

## Compute pooled SD, 95%CI for muT-muR and print results

s_p <- sqrt(((n.T-1)*sdl.T^2+(n.R-1)*sdl.R^2)/(n.T+n.R-2))   # Pooled SD
df <- n.T+n.R-2                                              # Degrees of Freedom

ybarl.diff <- ybarl.T - ybarl.R                              # Mean Difference
se.diff <- s_p * sqrt(1/n.T + 1/n.R)                         # Standard Error of mean diff
t.025 <- qt(.975, df)                                          # Critical t-value for 95% CI

mul.LB <- ybarl.diff - t.025 * se.diff                        # Lower Confidence Limit
mul.UB <- ybarl.diff + t.025 * se.diff                        # Upper Confidence Limit

# Print out summary of 95%CI for muT-muR
ci.out <- cbind(df, s_p, ybarl.diff, se.diff, t.025, mul.LB, mul.UB)
colnames(ci.out) <- c("df","pooled SD", "mean diff", "Std Err", "t(.975)", "Lower Bound", "Upper Bound")
round(ci.out, 5)

## Use t.test function for 95%CI for muT-muR

AUC.ttest <- t.test(AUC.sulf ~ form.AUC.sulf, var.equal=TRUE, conf.level=0.95)
AUC.ttest


# P.1.B -------------------------------------------------------------------

## Use var.test function for 95%CI for s1^2/s2^2

AUC.vartest <- var.test(AUC.sulf ~ form.AUC.sulf, var.equal=TRUE, conf.level=0.95)
AUC.vartest


# P.1.C -------------------------------------------------------------------

# Create var X to represent function for Test/Ref
X <- c(rep(1,23), rep(0,23))

lm_X <- lm(AUC.sulf ~ X)
lm_X

plot(AUC.sulf ~ X, main="AUC.sulf vs Test/Reference Formulation")
abline(lm_X)


# Question 3 --------------------------------------------------------------

# Read in data
spain_latlong <- read.csv(
  "https://users.stat.ufl.edu/~winner/data/spain_latlong.csv")
attach(spain_latlong); names(spain_latlong)
head(spain_latlong)
tail(spain_latlong)

# P.3.a. Obtain scatterplots of WGS ref - GIS - (Y) versus Ptolemy (X) for lat and long

# Define X and Y for lat and long
Y_long_wgs <- Long_wgs
X_long_ptol <- Long_ptol
Y_lat_wgs <- Lat_wgs
X_lat_ptol <- Lat_ptol

# Plots for Y versus X (lat and long)
par(mfrow=c(1,2))

plot(Y_long_wgs ~ X_long_ptol, main="WGS Longitude vs Ptolemy Longitude")

plot(Y_lat_wgs ~ X_lat_ptol, main="WGS Latitude vs Ptolemy Latitude")

# P.3.b. Fit simple linear regression models, relating GIS (Y) to Ptolemy (X)

lm_long <- lm(Y_long_wgs ~ X_long_ptol)
lm_lat <- lm(Y_lat_wgs ~ X_lat_ptol)

# Plot scatterplots again, now with linear models fitted
par(mfrow=c(1,2))

plot(Y_long_wgs ~ X_long_ptol, main="WGS Longitude vs Ptolemy Longitude")
abline(lm_long)

plot(Y_lat_wgs ~ X_lat_ptol, main="WGS Latitude vs Ptolemy Latitude")
abline(lm_lat)

# P.3.C. Test whether there is a positive association between GIS and Ptolemy

summary(lm_long)
# H0: B1 <= 0, Ha: B1 > 0. Since we have a p-value that is very small (<2e-16), we reject the null hypothesis
# Therefore, there is a positive association between WGS Longitude (GIS) and Ptolemy Longitude.

summary(lm_lat)
# H0: B1 <= 0, Ha: B1 > 0. Since we have a p-value that is very small (<2e-16), we reject the null hypothesis
# Therefore, there is a positive association between WGS Latitude (GIS) and Ptolemy Latitude.

# P.3.D. C95% for b1 > 0

confint(lm_long)
# The 95% confidence interval for B1 is (0.7213,0.8684). We are 95% confident the true value for B1 (the mean change) 
# falls between 0.7213 and 0.8684. Since the interval does not contain 0, we can say the mean change in GIS as Ptolemy 
# is “increased by 1 unit”, or that the mean change is positive.

confint(lm_lat)
# The 95% confidence interval for B1 is (0.7830,0.9498). We are 95% confident the true value for B1 (the mean change) 
# falls between 0.7830 and 0.9498. Since the interval does not contain 0, we can say the mean change in GIS as Ptolemy 
# is “increased by 1 unit”, or that the mean change is positive.

# P.3.E. ANOVAs, F-tests, coefficients of correlation and determination

# Linear model for longitude
anova(lm_long)
aov(Y_long_wgs ~ X_long_ptol)
cor(X_long_ptol, Y_long_wgs) # X, Y

# R-squared (coeff. of determination)
summary(lm_long)$r.squared

# F-test 
long.var.test <- var.test(X_long_ptol, Y_long_wgs, var.equal=TRUE, conf.level=0.95) # X, Y
long.var.test

# Linear model for latitude
anova(lm_lat)
aov(Y_lat_wgs ~ X_lat_ptol)
cor(X_lat_ptol, Y_lat_wgs) # X, Y

# R-squared (coeff. of determination)
summary(lm_lat)$r.squared
 
# F-test
lat.var.test <- var.test(X_lat_ptol, Y_lat_wgs, var.equal=TRUE, conf.level=0.95) # X, Y
lat.var.test