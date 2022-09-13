# Homework 1


# Question 1 --------------------------------------------------------------

# Read in data

spdata <- read.csv("http://www.stat.ufl.edu/~winner/data/bioequiv_sulf.csv")
attach(spdata); names(spdata)
spdata

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

AUC.ttest <- t.test(ref.AUC.sulf ~ form.AUC.sulf, var.equal=TRUE, conf.level=0.95)
AUC.ttest



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

# P.3.c. Test whether there is a positive association between GIS and Ptolemy