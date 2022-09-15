####################################################################
####### WNBA 2010-2019

wnba <- read.csv(
   "https://www.stat.ufl.edu/~winner/data/wnba_20102019_ATS_OU.csv")
attach(wnba); names(wnba)

tmPrd <- (OU - tmSprd)/2
summary(tmPrd)


##obtain number of games from 2010-2019
(N.games <- length(OU))

tm.reg1 <- lm(teamPts ~ tmPrd)
summary(tm.reg1)
anova(tm.reg1)

tmPrd.grid <- seq(60,110,.01)

  
win.graph(height=5.5, width=7.0)
plot(teamPts ~ tmPrd, pch=16, cex=.35, 
     main="Actual versus Predicted Team Points - WNBA 2010-2019")
abline(tm.reg1, col="blue", lwd=2)
lines(tmPrd.grid, predict(tm.reg1, list(tmPrd=tmPrd.grid), int="c")[,2], 
     col="red", lwd=2)
lines(tmPrd.grid, predict(tm.reg1, list(tmPrd=tmPrd.grid), int="c")[,3], 
     col="red", lwd=2)
lines(tmPrd.grid, predict(tm.reg1, list(tmPrd=tmPrd.grid), int="p")[,2], 
     col="purple", lwd=2)
lines(tmPrd.grid, predict(tm.reg1, list(tmPrd=tmPrd.grid), int="p")[,3], 
     col="purple", lwd=2)
legend("bottomright", c("Fitted", "CI for Mean", "PI for Game"), lty=1,
     col=c("blue", "red", "purple"))

tm.reg1 <- lm(teamPts ~ tmPrd)
summary(tm.reg1)
anova(tm.reg1)

tm.reg2 <- lm(teamPts ~ -1, offset = tmPrd)
anova(tm.reg2)

anova(tm.reg2, tm.reg1)

