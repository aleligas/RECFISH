library(mgcv)
library(fmsb)
library(Hmisc)

# setting the working directory and select the datasets
setwd("C:/Ale/PROGETTI IN CORSO/Tender EASME-EMFF/RECFISH/WP2 RECFISH/ARA GSA9")
data <- read.table("TATB9ARAsst.csv", sep=",",header = T)
str(data)
#     data=data[data$year>1994,]
grid <- read.table("grid_0.01_GSA9_SST.csv", sep=";",header = T)
str(grid)
# computation of mean haul depth
data$depth <- (data$shooting_depth+data$hauling_depth)/2

# computation of swept area
data$swept_area <- data$distance*data$wing_opening/10000000

# computation of aboundance index
data$kg_km2 <-data$ptot/data$swept_area/1000

# computation of density index
data$n_km2 <- data$nbtot/data$swept_area

# hour extraction
data$hour<-as.numeric(as.character(ifelse(nchar(data$shooting_time)==3, substr(data$shooting_time,1,1), substr(data$shooting_time,1,2))))
GSA  <- unique(data$area)

# selection of the variables for the analysis
data <- data.frame(kg_km2 = data$kg_km2, 
                   X = data$X, 
                   Y = data$Y, 
                   depth = data$depth, 
                   year = as.numeric(data$year), 
                   month = as.numeric(data$month), 
                   hour = data$hour, 
                   SST = data$SST,
                   vessel = as.numeric(data$vessel))

# Inspection of the plots by depths for the selection of the depth range:
for (i in 1:nrow(data)){
if (data$depth[i] >= 10 & data$depth[i] < 50){
data$stratum[i]=1
} else if (data$depth[i] >= 50 & data$depth[i] < 100) {
data$stratum[i]=2
} else if (data$depth[i] >= 100 & data$depth[i] < 200) {
data$stratum[i]=3
} else if (data$depth[i] >= 200 & data$depth[i] < 500){
data$stratum[i]=4
} else if (data$depth[i] >= 500 & data$depth[i] <= 800) {
data$stratum[i]=5}
}

plot(data$depth, data$kg_km2)
ggplot(data=data, aes(x=depth,y= kg_km2)) + geom_histogram(stat="identity",colour = "blue", fill = "blue", binwidth = 0.5) + facet_grid(stratum~ .)

# selection of depth range
#   data <-  data[data$depth>200 & data$depth<=500,]
data <-  data[data$depth>500,]
# VIF analysis
mod_lm <- lm(kg_km2 ~ X + Y + depth + year + month + hour + SST, data=data)
vif<-VIF(lm(mod_lm, data = data)); vif

# correlation matrix
rcorr(as.matrix(data)) 

#### PAIR PLOT FOR COLLINEARITY CHECK
panel.smooth2<-function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                         cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = 1, ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor ,...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  #r <- abs(cor(x, y))
  r <- (cor(x, y))
  #txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- format(c(r, 0.123456789), digits = 1)[1]
  txt <- paste(prefix, txt, sep = "")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * r)
  text(0.5, 0.5, txt, cex = 3*abs(r))
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
}

#tiff("pairplot.tiff", height = 20, width = 20, units = "cm", compression = "lzw", res = 600) 
#op<-par(mfrow=c(1,1), mar=c(4,4,2,2))
Z<-cbind(data$X, data$Y, data$depth, data$year, data$SST, data$month, data$hour, data$vessel)
colnames(Z)<- c("Lon", "Lat", "Depth", "Year", "SST", "Month", "Hour", "Vessel")
pairs(Z,las=1 ,lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist, cex=0.9, cex.labels=1.0)
# box()
# dev.off()


# transformation of dependent variable
# data$kg_sqrt <- sqrt(data$kg_km2)
data$kg_sqrt <- log(data$kg_km2+1)

# Stepwise forward inclusion for GAM modeling
mod1 <- gam(kg_sqrt ~ s(X)+0, family= gaussian (link="identity"),data=data, select=T) 
mod2 <- gam(kg_sqrt ~ s(Y)+0, family= gaussian (link="identity"),data=data, select=T)
mod3 <- gam(kg_sqrt ~ s(depth)+0, family= gaussian (link="identity"),data=data, select=T)
mod4 <- gam(kg_sqrt ~ s(year)+0, family= gaussian (link="identity"),data=data, select=T)
mod5 <- gam(kg_sqrt ~ factor(month)+0, family= gaussian (link="identity"),data=data, select=T)
mod6 <- gam(kg_sqrt ~ s(hour)+0, family= gaussian (link="identity"),data=data, select=T)
mod7 <- gam(kg_sqrt ~ s(SST)+0, family= gaussian (link="identity"),data=data, select=T)
mod8 <- gam(kg_sqrt ~ factor(vessel)+0, family= gaussian (link="identity"),data=data, select=T)
# from the summary the GCV and the explained deviance can be inspected
summary(mod1) 
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)
summary(mod8)

AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)

# once that the best basic model was detected, the variables should be included one at the time in the same way
#.....
mod1.1 <- gam(kg_sqrt ~ factor(vessel)+ s(X), family= gaussian (link="identity"),data=data, select=T) 
mod1.2 <- gam(kg_sqrt ~ factor(vessel)+ s(Y), family= gaussian (link="identity"),data=data, select=T)
mod1.3 <- gam(kg_sqrt ~ factor(vessel)+ s(depth), family= gaussian (link="identity"),data=data, select=T)
mod1.4 <- gam(kg_sqrt ~ factor(vessel)+ s(year), family= gaussian (link="identity"),data=data, select=T)
mod1.5 <- gam(kg_sqrt ~ factor(month)+ factor(vessel), family= gaussian (link="identity"),data=data, select=T)
mod1.6 <- gam(kg_sqrt ~ factor(vessel)+ s(hour), family= gaussian (link="identity"),data=data, select=T)
mod1.7 <- gam(kg_sqrt ~ factor(vessel)+ s(SST), family= gaussian (link="identity"),data=data, select=T)

summary(mod1.1) 
summary(mod1.2)
summary(mod1.3)
summary(mod1.4)
summary(mod1.5)
summary(mod1.6)
summary(mod1.7)

AIC(mod1.1,mod1.2,mod1.3,mod1.4,mod1.5,mod1.6,mod1.7)


mod2.1 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(X), family= gaussian (link="identity"),data=data, select=T)
mod2.2 <- gam(kg_sqrt ~ factor(vessel)+ s(X) + s(Y), family= gaussian (link="identity"),data=data, select=T)
mod2.3 <- gam(kg_sqrt ~ factor(vessel)+ s(X) + s(year), family= gaussian (link="identity"),data=data, select=T)
mod2.4 <- gam(kg_sqrt ~ factor(vessel)+ s(X) + s(SST), family= gaussian (link="identity"),data=data, select=T)
mod2.5 <- gam(kg_sqrt ~ factor(vessel)+ s(X) + s(hour), family= gaussian (link="identity"),data=data, select=T)
mod2.6 <- gam(kg_sqrt ~ factor(month)+ s(X) + factor(vessel), family= gaussian (link="identity"),data=data, select=T)

summary(mod2.1) 
summary(mod2.2)
summary(mod2.3)
summary(mod2.4)
summary(mod2.5)
summary(mod2.6)

AIC(mod2.1,mod2.2,mod2.3,mod2.4,mod2.5,mod2.6)

mod3.1 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(Y) + s(X), family= gaussian (link="identity"),data=data, select=T)
mod3.2 <- gam(kg_sqrt ~ factor(vessel)+ s(X) + s(Y) + s(SST), family= gaussian (link="identity"),data=data, select=T)
mod3.3 <- gam(kg_sqrt ~ factor(vessel)+ s(X) + s(Y) + s(year), family= gaussian (link="identity"),data=data, select=T)
mod3.4 <- gam(kg_sqrt ~ factor(vessel)+ s(X) + s(Y) + s(hour), family= gaussian (link="identity"),data=data, select=T)
mod3.5 <- gam(kg_sqrt ~ factor(month)+ s(X) + s(Y) + factor(vessel), family= gaussian (link="identity"),data=data, select=T)

summary(mod3.1) 
summary(mod3.2)
summary(mod3.3)
summary(mod3.4)
summary(mod3.5)

AIC(mod3.1,mod3.2,mod3.3,mod3.4,mod3.5)

mod4.1 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(Y) + s(X) + s(SST), family= gaussian (link="identity"),data=data, select=T)
mod4.2 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(Y) + s(X) + s(year), family= gaussian (link="identity"),data=data, select=T)
mod4.3 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(Y) + s(X) + s(hour), family= gaussian (link="identity"),data=data, select=T)
mod4.4 <- gam(kg_sqrt ~ factor(month)+ s(depth) + s(Y) + s(X) + factor(vessel), family= gaussian (link="identity"),data=data, select=T)

summary(mod4.1) 
summary(mod4.2)
summary(mod4.3)
summary(mod4.4)

AIC(mod4.1,mod4.2,mod4.3,mod4.4)

mod5.1 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(Y) + s(X) + s(year) + s(SST), family= gaussian (link="identity"),data=data, select=T)
mod5.2 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(Y) + s(X) + s(year) + s(hour), family= gaussian (link="identity"),data=data, select=T)
mod5.3 <- gam(kg_sqrt ~ factor(month)+ s(depth) + s(Y) + s(X) + s(year) + factor(vessel), family= gaussian (link="identity"),data=data, select=T)

summary(mod5.1) 
summary(mod5.2)
summary(mod5.3)

AIC(mod5.1,mod5.2,mod5.3)

mod6.1 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(Y) + s(X) + s(year) + s(hour) +s(SST), family= gaussian (link="identity"),data=data, select=T)
mod6.2 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(Y) + s(X) + s(year) + s(hour) + factor(month), family= gaussian (link="identity"),data=data, select=T)

summary(mod6.1) 
summary(mod6.2)

AIC(mod6.1,mod6.2)


mod2 <- gam(kg_sqrt ~ factor(vessel)+ s(depth) + s(Y) + s(X) + s(year) + s(hour) +s(SST), family= gaussian (link="identity"),data=data, select=T)
mod1 <- gam(kg_sqrt ~ s(depth) + s(Y) + s(X) + s(year) + s(hour) +s(SST), family= gaussian (link="identity"),data=data, select=T)
mod <- gam(kg_sqrt ~ s(depth) + s(Y) + s(X) + s(year) + s(SST) +factor(month), family= gaussian (link="identity"),data=data, select=T)

data$vessel
summary(mod)
summary(mod1)
summary(mod2)
AIC(mod, mod1,mod2)

# inspection of residuals and diagnostic plots
res <- residuals(mod)
plot(res, ylim=c(-5.5,5.5))
abline(0,0, col="red")
par(mfrow=c(2,2))
gam.check(mod)
par(mfrow=c(3,3))
plot(mod,select =1)
plot(mod, select =2)
plot(mod, select =3)
plot(mod, select =4)
plot(mod, select =5)
plot(mod, select =6)
#termplot(mod, se=T, ylim=c(1,4))
termplot(mod, se=T)


#residui pearson

E.m4<- resid(mod, type = "pearson")
Fit.m4 <- fitted(mod)
#e4<-resid(M3)
plot(x = Fit.m4, y = E.m4, xlab = "Fitted values",
     ylab = "Residuals")
#plot(x = Fit.m4, y = E.m4, xlab = "Fitted values",
#    ylab = "Residuals", main = "logKG.KM2")
abline(0,0)
tmp <- loess(E.m4 ~Fit.m4 ,span=0.75)
tmp2 <- predict(tmp,se=T)
I1 <- order(Fit.m4)
lines(Fit.m4[I1], tmp2$fit[I1], lty=1)
lines(Fit.m4[I1], tmp2$fit[I1] + 2*tmp2$se.fit[I1], lty = 2)
lines(Fit.m4[I1], tmp2$fit[I1] - 2*tmp2$se.fit[I1], lty = 2)

# ASSUNZIONE SU INDIPENDENZA, SI PLOTTANO I RESIDUI CONTRO LA VARIABILE INDIPENDENTE:
var1 <- data$year
plot(x = var1, y = E.m4, xlab = "Year",
     ylab = "Residuals")
abline(0,0)
tmp <- loess(E.m4 ~var1 ,span=0.8)
tmp2 <- predict(tmp,se=T)
I1 <- order(var1)
lines(var1[I1], tmp2$fit[I1], lty=1)
lines(var1[I1], tmp2$fit[I1] + 2*tmp2$se.fit[I1], lty = 2)
lines(var1[I1], tmp2$fit[I1] - 2*tmp2$se.fit[I1], lty = 2)

var1 <- data$Y
plot(x = var1, y = E.m4, xlab = "Latitude",
     ylab = "Residuals")
abline(0,0)
tmp <- loess(E.m4 ~var1 ,span=0.8)
tmp2 <- predict(tmp,se=T)
I1 <- order(var1)
lines(var1[I1], tmp2$fit[I1], lty=1)
lines(var1[I1], tmp2$fit[I1] + 2*tmp2$se.fit[I1], lty = 2)
lines(var1[I1], tmp2$fit[I1] - 2*tmp2$se.fit[I1], lty = 2)

var1 <- data$depth
plot(x = var1, y = E.m4, xlab = "Depth",
     ylab = "Residuals")
abline(0,0)
tmp <- loess(E.m4 ~var1 ,span=0.8)
tmp2 <- predict(tmp,se=T)
I1 <- order(var1)
lines(var1[I1], tmp2$fit[I1], lty=1)
lines(var1[I1], tmp2$fit[I1] + 2*tmp2$se.fit[I1], lty = 2)
lines(var1[I1], tmp2$fit[I1] - 2*tmp2$se.fit[I1], lty = 2)

par(mar=c(5,5,2,2))
boxplot(E.m4 ~ month, las=1, data = data)
abline(h = 0, lty = 2)
title(xlab = "Month",ylab = "Residuals")

par(mar=c(5,5,2,2))
boxplot(E.m4 ~ vessel, las=1, data = data)
abline(h = 0, lty = 2)
title(xlab = "Month",ylab = "Residuals")

#Normality
par(mar=c(5,5,2,2))
hist(E.m4, ylab = "Frequency", xlab = "Residuals", las=1,breaks=16, cex.lab=1.1, cex.axis=1.1)

#Or qq-plot
par(mar=c(5,5,2,2))
qqnorm(E.m4, main=NULL, lwd=1.5,cex.lab=1.1, las=1,cex.axis=1, bty="l", xlab="Theoretical quantiles", ylab="Std. Dev. quantiles")
abline(0.0,1, lty=1, lwd=1.5)

# cOOK'S DISTANCE
par(mfrow=c(1,1))
plot(cooks.distance(mod), ylab="Cook's distance", type = "h", ylim=c(0,1))
abline(h=1, col=1,lwd=2)

#analysis_stratum1 <- F
#analysis_stratum2 <- F
#analysis_stratum3 <- F
#analysis_stratum4 <- F # or T for 200-500
analysis_stratum5 <- T
#analysis_stratum6 <- T

#stratum_1 <- grid[grid$depth >= 10 & grid$depth < 50, ]
#stratum_2 <- grid[grid$depth >= 50 & grid$depth < 100, ]
#stratum_3 <- grid[grid$depth >= 100 & grid$depth < 200, ]
#stratum_4 <- grid[grid$depth >= 200 & grid$depth < 500, ]
stratum_5 <- grid[grid$depth >= 500 & grid$depth <= 800,]
#stratum_6 <- grid[grid$depth >= 500, ]

res_table <- data.frame(seq(1994,2017, 1)) # or 1995
#res_table$index_1 <- NA
#res_table$index_2 <- NA
#res_table$index_3 <- NA
#res_table$index_4 <- NA
res_table$index_5 <- NA
#res_table$index_6 <- NA

stratification <- read.table("Stratification_Scheme.csv", sep=";", header=T)
#area_s1 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==1,5])
#area_s2 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==2,5])
#area_s3 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==3,5])
#area_s4 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==4,5])
area_s5 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==5,5])
#area_s6 <- area_s5


env94<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_1994_SST, year=1994, month=6, hour=12, vessel=1)
env95<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_1995_SST, year=1995, month=6, hour=12, vessel=1)
env96<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_1996_SST, year=1996, month=6, hour=12, vessel=1)
env97<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_1997_SST, year=1997, month=6, hour=12, vessel=1)
env98<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_1998_SST, year=1998, month=6, hour=12, vessel=1)
env99<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_1999_SST, year=1999, month=6, hour=12, vessel=1)
env00<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2000_SST, year=2000, month=6, hour=12, vessel=1)
env01<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2001_SST, year=2001, month=6, hour=12, vessel=1)
env02<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2002_SST, year=2002, month=6, hour=12, vessel=1)
env03<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2003_SST, year=2003, month=6, hour=12, vessel=1)
env04<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2004_SST, year=2004, month=6, hour=12, vessel=1)
env05<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2005_SST, year=2005, month=6, hour=12, vessel=1)
env06<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2006_SST, year=2006, month=6, hour=12, vessel=1)
env07<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2007_SST, year=2007, month=6, hour=12, vessel=1)
env08<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2008_SST, year=2008, month=6, hour=12, vessel=1)
env09<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2009_SST, year=2009, month=6, hour=12, vessel=1)
env10<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2010_SST, year=2010, month=6, hour=12, vessel=1)
env11<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2011_SST, year=2011, month=6, hour=12, vessel=1)
env12<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2012_SST, year=2012, month=6, hour=12, vessel=1)
env13<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2013_SST, year=2013, month=6, hour=12, vessel=1)
env14<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2014_SST, year=2014, month=6, hour=12, vessel=1)
env15<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2015_SST, year=2015, month=6, hour=12, vessel=1)
env16<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2016_SST, year=2016, month=6, hour=12, vessel=1)
env17<-data.frame(X=grid$X, Y=grid$Y, depth=grid$depth, SST=grid$June_2017_SST, year=2017, month=6, hour=12, vessel=1)


pred94 <- predict(mod, newdata = data.frame(env94), type="response")
pred95 <- predict(mod, newdata = data.frame(env95), type="response")
pred96 <- predict(mod, newdata = data.frame(env96), type="response")
pred97 <- predict(mod, newdata = data.frame(env97), type="response")
pred98 <- predict(mod, newdata = data.frame(env98), type="response")
pred99 <- predict(mod, newdata = data.frame(env99), type="response")
pred00 <- predict(mod, newdata = data.frame(env00), type="response")
pred01 <- predict(mod, newdata = data.frame(env01), type="response")
pred02 <- predict(mod, newdata = data.frame(env02), type="response")
pred03 <- predict(mod, newdata = data.frame(env03), type="response")
pred04 <- predict(mod, newdata = data.frame(env04), type="response")
pred05 <- predict(mod, newdata = data.frame(env05), type="response")
pred06 <- predict(mod, newdata = data.frame(env06), type="response")
pred07 <- predict(mod, newdata = data.frame(env07), type="response")
pred08 <- predict(mod, newdata = data.frame(env08), type="response")
pred09 <- predict(mod, newdata = data.frame(env09), type="response")
pred10 <- predict(mod, newdata = data.frame(env10), type="response")
pred11 <- predict(mod, newdata = data.frame(env11), type="response")
pred12 <- predict(mod, newdata = data.frame(env12), type="response")
pred13 <- predict(mod, newdata = data.frame(env13), type="response")
pred14 <- predict(mod, newdata = data.frame(env14), type="response")
pred15 <- predict(mod, newdata = data.frame(env15), type="response")
pred16 <- predict(mod, newdata = data.frame(env16), type="response")
pred17 <- predict(mod, newdata = data.frame(env17), type="response")

p94 <- mean(exp(pred94)+1)
p95 <- mean(exp(pred95)+1)
p96 <- mean(exp(pred96)+1)
p97 <- mean(exp(pred97)+1)
p98 <- mean(exp(pred99)+1)
p99 <- mean(exp(pred99)+1)
p00 <- mean(exp(pred00)+1)
p01 <- mean(exp(pred01)+1)
p02 <- mean(exp(pred02)+1)
p03 <- mean(exp(pred03)+1)
p04 <- mean(exp(pred04)+1)
p05 <- mean(exp(pred05)+1)
p06 <- mean(exp(pred06)+1)
p07 <- mean(exp(pred07)+1)
p08 <- mean(exp(pred08)+1)
p09 <- mean(exp(pred09)+1)
p10 <- mean(exp(pred10)+1)
p11 <- mean(exp(pred11)+1)
p12 <- mean(exp(pred12)+1)
p13 <- mean(exp(pred13)+1)
p14 <- mean(exp(pred14)+1)
p15 <- mean(exp(pred15)+1)
p16 <- mean(exp(pred16)+1)
p17 <- mean(exp(pred17)+1)

predictions <- c(p94,p95,p96,p97,p98,p99,p00,p01,p02,p03,p04,p05, p06, p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17)

years <- c(1994,1995,1996, 1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017)

# 500-800m
time_series <- c(5.03,3.90,2.17,4.50,4.25,3.96,3.07,2.54,1.52,2.35,4.54,2.62,3.43,3.92,5.94,1.63,16.01,12.68,7.66,2.00,3.67,2.29,1.98,1.79) # calcolati dallo script indici medits nella cartella stecf/medits2015

par(mfrow=c(1,1))
plot(years, time_series, col="red", xlab="year", ylab= "kg/km^2", pch=16)
lines(years, time_series, col="red")
points(years, predictions, col="black", pch=16)
lines(years, predictions, col="black")
legend("topright", c("Time series","Prediction"), col=c("red", "black"), lwd=1, pch=16)


Index<-cbind(years, predictions)
colnames(Index)=c("Year","Index")
write.table(Index,"IndexARA500800.csv",sep=";",row.names=F)



############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################



analysis_stratum1 <- F
analysis_stratum2 <- F
analysis_stratum3 <- F
analysis_stratum4 <- F # or T for 200-500
analysis_stratum5 <- T
#analysis_stratum6 <- T

stratum_1 <- grid[grid$depth >= 10 & grid$depth < 50, ]
stratum_2 <- grid[grid$depth >= 50 & grid$depth < 100, ]
stratum_3 <- grid[grid$depth >= 100 & grid$depth < 200, ]
stratum_4 <- grid[grid$depth >= 200 & grid$depth < 500, ]
stratum_5 <- grid[grid$depth >= 500 & grid$depth <= 800,]
#stratum_6 <- grid[grid$depth >= 500, ]

res_table <- data.frame(seq(1994,2017, 1)) # or 1995
res_table$index_1 <- NA
res_table$index_2 <- NA
res_table$index_3 <- NA
res_table$index_4 <- NA
res_table$index_5 <- NA
#res_table$index_6 <- NA

stratification <- read.table("Stratification_Scheme.csv", sep=";", header=T)
area_s1 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==1,5])
area_s2 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==2,5])
area_s3 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==3,5])
area_s4 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==4,5])
area_s5 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==5,5])
#area_s6 <- area_s5

i=1
for(i in 1:length(res_table[,1])){
  if (analysis_stratum1 == T){
    stratum_1$year <- res_table[i,1]
    stratum_1$month <- 6
    stratum_1$hour <- 12
    stratum_1$vessel <- 1
    #stratum_1$SST <- grid$June_1994_SST
        stratum_1$pred <- exp(predict(mod, newdata = data.frame(stratum_1)))-1
    res_table[i,2] <- mean(stratum_5$pred) *area_s1} else {
      area_s1 = 0
    }
  
  if (analysis_stratum2 == T){
    stratum_2$year <- res_table[i,1]
    stratum_2$month <- 6
    stratum_2$hour <- 12
    stratum_2$vessel <- 1
    stratum_2$pred <- exp(predict(mod, newdata = data.frame(stratum_2)))-1
    res_table[i,3] <- mean(stratum_5$pred) *area_s2} else {
      area_s2 = 0
    }
  
  if (analysis_stratum3 == T){
    stratum_3$year <- res_table[i,1]
    stratum_3$month <- 6
    stratum_3$hour <- 12
    stratum_3$vessel <- 1
    stratum_3$pred <- exp(predict(mod, newdata = data.frame(stratum_3)))-1
    res_table[i,4] <- mean(stratum_5$pred) *area_s3} else {
      area_s3 = 0
    }
  
  if (analysis_stratum4 == T){
    stratum_4$year <- res_table[i,1]
    stratum_4$month <- 6
    stratum_4$hour <- 12
    stratum_4$vessel <- 1
    stratum_4$pred <- exp(predict(mod, newdata = data.frame(stratum_4)))-1
    res_table[i,5] <- mean(stratum_5$pred) *area_s4} else {
      area_s4 = 0
    }
  
  if (analysis_stratum5 == T){
   stratum_5$year <- res_table[i,1]
    stratum_5$month <- 6
    stratum_5$hour <- 12
    stratum_5$vessel <- 1
    stratum_5$pred <- exp(predict(mod, newdata = data.frame(stratum_5)))-1 
  res_table[i,6] <- mean(stratum_5$pred) *area_s5} else {
   area_s5 = 0
  }
  
  
 #sum_res <- c(res_table[i,2],res_table[i,3],res_table[i,4],res_table[i,5],res_table[i,6],res_table[i,7])
 sum_res <- c(res_table[i,2],res_table[i,3],res_table[i,4],res_table[i,5],res_table[i,6])
  #res_table[i, 8]<- sum(sum_res[!is.na(sum_res)])/sum(area_s1,area_s2,area_s3,area_s6)
  #res_table[i, 7]<- sum(sum_res[!is.na(sum_res)])/sum(area_s1,area_s2,area_s3,area_s4, area_s5)
 res_table[i, 7]<- sum(sum_res[!is.na(sum_res)])/area_s5
  }

colnames(res_table) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "Indices")
#colnames(res_table) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "Stratum 6", "Indices")
par(mfrow=c(1,1))
plot(res_table[,1], res_table[,7], type="b", xlab="year", ylab= "index")

# 500-800m
time_series <- c(5.03,3.90,2.17,4.50,4.25,3.96,3.07,2.54,1.52,2.35,4.54,2.62,3.43,3.92,5.94,1.63,16.01,12.68,7.66,2.00,3.67,2.29,1.98,1.79) # calcolati dallo script indici medits nella cartella stecf/medits2015

par(mfrow=c(1,1))
plot(res_table[,1], time_series, col="red", xlab="year", ylab= "kg/km^2", pch=16)
lines(res_table[,1], time_series, col="red")
points(res_table[,1], res_table[,7], col="black", pch=16)
lines(res_table[,1], res_table[,7], col="black")
legend("topright", c("Time series","Prediction"), col=c("red", "black"), lwd=1, pch=16)


Index=res_table[,c(1,8)]
colnames(Index)=c("Year","Index")
write.table(Index,"Index200800.csv",sep=";",row.names=F)