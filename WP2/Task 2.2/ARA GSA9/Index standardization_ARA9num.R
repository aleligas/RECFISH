library(mgcv)
library(fmsb)
library(Hmisc)

# setting the working directory and select the datasets
setwd("C:/Ale/PROGETTI IN CORSO/Tender EASME-EMFF/RECFISH/WP2 RECFISH/TASK 2.2/ARA GSA9")
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

# selection of the variables for the analysis (kg or n)
data <- data.frame(n_km2 = data$n_km2, 
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

plot(data$depth, data$n_km2)
ggplot(data=data, aes(x=depth,y= n_km2)) + geom_histogram(stat="identity",colour = "blue", fill = "blue", binwidth = 0.5) + facet_grid(stratum~ .)

# selection of depth range
#   data <-  data[data$depth>200 & data$depth<=500,]
data <-  data[data$depth>500,]
# VIF analysis
mod_lm <- lm(n_km2 ~ X + Y + depth + year + month + hour + SST, data=data)
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

############ CORVIF 
source("HighstatLibV10.R")
colnames(Z)<- c("Lon", "Lat", "Depth", "Year", "SST", "Month", "Hour", "Vessel")
corvif(Z)



# transformation of dependent variable
 data$kg_sqrt <- sqrt(data$n_km2)
####data$kg_sqrt <- log(data$n_km2+1) #### num instead of kg

# Stepwise forward inclusion for GAM modeling
mod1 <- gam(kg_sqrt ~ s(X, depth), family= gaussian (link="identity"),data=data, select=T)
mod2 <- gam(kg_sqrt ~ s(X) + s(depth), family= gaussian (link="identity"),data=data, select=T) 
summary(mod1) 
summary(mod2)

AIC(mod1,mod2)

mod1 <- gam(kg_sqrt ~ s(X, depth) + s(year), family= gaussian (link="identity"),data=data, select=T) 
mod2 <- gam(kg_sqrt ~ s(X, depth) + factor(month), family= gaussian (link="identity"),data=data, select=T)
mod3 <- gam(kg_sqrt ~ s(X, depth) + s(hour), family= gaussian (link="identity"),data=data, select=T)
mod4 <- gam(kg_sqrt ~ s(X, depth) + s(SST), family= gaussian (link="identity"),data=data, select=T)

# from the summary the GCV and the explained deviance can be inspected
summary(mod1) 
summary(mod2)
summary(mod3)
summary(mod4)
AIC(mod1,mod2,mod3,mod4)

mod1 <- gam(kg_sqrt ~ s(X, depth) + s(year) + factor(month), family= gaussian (link="identity"),data=data, select=T)
mod2 <- gam(kg_sqrt ~ s(X, depth) + s(year) + s(hour), family= gaussian (link="identity"),data=data, select=T) 
mod3 <- gam(kg_sqrt ~ s(X, depth) + s(year) + s(SST), family= gaussian (link="identity"),data=data, select=T) 

summary(mod1) 
summary(mod2)
summary(mod3)
AIC(mod1,mod2,mod3)

# THIS IS THE BEST MODEL IN TERMS OF BOTH RESIDUALS AND PREDICTIONS
data$kg_sqrt <- sqrt(data$n_km2)
mod <- gam(kg_sqrt ~ s(X,depth) +  s(year, fx=T) , family= gaussian (link="identity"),data=data, method="REML")


summary(mod)

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

var1 <- data$X
plot(x = var1, y = E.m4, xlab = "Longitude",
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
title(xlab = "Vessel",ylab = "Residuals")

#Normality
par(mar=c(5,5,2,2))
hist(E.m4, ylab = "Frequency", xlab = "Residuals", las=1,breaks=16, cex.lab=1.1, cex.axis=1.1, main=NULL)

#Or qq-plot
par(mar=c(5,5,2,2))
qqnorm(E.m4, main=NULL, lwd=1.5,cex.lab=1.1, las=1,cex.axis=1, bty="l", xlab="Theoretical quantiles", ylab="Std. Dev. quantiles")
abline(0.0,1.8, lty=1, lwd=1.5)

# cOOK'S DISTANCE
par(mfrow=c(1,1))
plot(cooks.distance(mod), ylab="Cook's distance", type = "h", ylim=c(0,1))
abline(h=1, col=1,lwd=2)

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

se_table <- data.frame(seq(1994,2017, 1))
se_table$index_1 <- NA
se_table$index_2 <- NA
se_table$index_3 <- NA
se_table$index_4 <- NA
se_table$index_5 <- NA
se_table$se <- NA


stratification <- read.table("Stratification_Scheme.csv", sep=";", header=T)
area_s1 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==1,5])
area_s2 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==2,5])
area_s3 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==3,5])
area_s4 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==4,5])
area_s5 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==5,5])
#area_s6 <- area_s5

std <- function(x) sd(x)/sqrt(length(x))


i=1
for(i in 1:length(res_table[,1])){
  if (analysis_stratum1 == T){
    stratum_1$year <- res_table[i,1]
    stratum_1$month <- 6
    stratum_1$hour <- 12
    stratum_1$vessel <- 1
    #stratum_1$SST <- grid$June_1994_SST
        stratum_1$pred <- (predict(mod, newdata = data.frame(stratum_1)))^2
    res_table[i,2] <- mean(stratum_5$pred) *area_s1
    se_strata1 <-  std(stratum_1$pred) 
    se_table[i,2] <- se_strata1 *area_s1 
    } else {
      area_s1 = 0
    }
  
  if (analysis_stratum2 == T){
    stratum_2$year <- res_table[i,1]
    stratum_2$month <- 6
    stratum_2$hour <- 12
    stratum_2$vessel <- 1
    stratum_2$pred <- (predict(mod, newdata = data.frame(stratum_2)))^2
    res_table[i,3] <- mean(stratum_5$pred) *area_s2
    se_strata2 <- std(stratum_2$pred) 
    se_table[i,3] <- se_strata2 *area_s2
    } else {
      area_s2 = 0
    }
  
  if (analysis_stratum3 == T){
    stratum_3$year <- res_table[i,1]
    stratum_3$month <- 6
    stratum_3$hour <- 12
    stratum_3$vessel <- 1
    stratum_3$pred <- (predict(mod, newdata = data.frame(stratum_3)))^2
    res_table[i,4] <- mean(stratum_5$pred) *area_s3
    se_strata3 <- std(stratum_3$pred) 
    se_table[i,4] <- se_strata3 *area_s3 
    } else {
      area_s3 = 0
    }
  
  if (analysis_stratum4 == T){
    stratum_4$year <- res_table[i,1]
    stratum_4$month <- 6
    stratum_4$hour <- 12
    stratum_4$vessel <- 1
    stratum_4$pred <- (predict(mod, newdata = data.frame(stratum_4)))^2
    res_table[i,5] <- mean(stratum_5$pred) *area_s4
    se_strata4 <- std(stratum_4$pred) 
    se_table[i,5] <- se_strata4 *area_s4
    } else {
      area_s4 = 0
    }
  
  if (analysis_stratum5 == T){
   stratum_5$year <- res_table[i,1]
    stratum_5$month <- 6
    stratum_5$hour <- 12
    stratum_5$vessel <- 1
    stratum_5$pred <- (predict(mod, newdata = data.frame(stratum_5)))^2 
  res_table[i,6] <- mean(stratum_5$pred) *area_s5
  se_strata5 <- std(stratum_5$pred) 
  se_table[i,6] <- se_strata5 *area_s5
  } else {
   area_s5 = 0
  }
  
  
  sum_res <- c(res_table[i,2],res_table[i,3],res_table[i,4],res_table[i,5],res_table[i,6])
  res_table[i, 7]<- sum(sum_res[!is.na(sum_res)])/sum(area_s1,area_s2,area_s3,area_s4,area_s5)
  
  sum_se <- c(se_table[i,2],se_table[i,3],se_table[i,4],se_table[i,5],se_table[i,6])
  se_table[i, 7]<- sum(sum_se[!is.na(sum_se)])/sum(area_s1,area_s2,area_s3,area_s4,area_s5)  
  
  }

colnames(res_table) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "Indices")
#colnames(res_table) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "Stratum 6", "Indices")
par(mfrow=c(1,1))
plot(res_table[,1], res_table[,7], type="b", xlab="year", ylab= "index")

# 500-800m
time_series <- c(186.20,162.779,83.577,149.943,149.664,173.949,111.04,72.48,55.544,81.90,145.88,79.68,95.31,113.71, 268.74,86.65,654.92,532.51,302.53,75.195,137.84,90.42,78.96,81.61) # calcolati dallo script indici medits nella cartella stecf/medits2015


par(mfrow=c(1,1))
plot(res_table[,1], time_series, col="red", xlab="year", ylab= "kg/km^2", pch=16)
lines(res_table[,1], time_series, col="red")
points(res_table[,1], res_table[,7], col="black", pch=16)
lines(res_table[,1], res_table[,7], col="black")
legend("topright", c("Time series","Prediction"), col=c("red", "black"), lwd=1, pch=16)

#######  PLOT CON ERRORI STANDARD
# 500-800m
time_series <- c(186.20,162.779,83.577,149.943,149.664,173.949,111.04,72.48,55.544,81.90,145.88,79.68,95.31,113.71, 268.74,86.65,654.92,532.51,302.53,75.195,137.84,90.42,78.96,81.61) 

se <- c(98.8,76.9,48.8,62.1,47.4,82.5,49.7,31.5,30.3,45.5,72.5,41.5,43.7,58.6,139.8,43.8,229.3,147.0,85.4,22.2,46.5,55.9,48.1,49.8)


#CI_up <- time_series + (se)
#CI_low <-  time_series - (se)

CI_up <- time_series + (1.96* se)
CI_low <-  time_series - (1.96* se)

par(mfrow=c(1,1))
plot(res_table[,1], time_series, col="red", xlab="year", ylab= "num/km^2", pch=16, ylim=c(0,1100))
lines(res_table[,1], time_series, col="red")
lines(res_table[,1], CI_up, col="red", lty=2)
lines(res_table[,1], CI_low, col="red", lty=2)

points(res_table[,1], res_table[,7], col="black", pch=16)
lines(res_table[,1], res_table[,7], col="black")
lines(res_table[,1], (res_table[,7] +(10.96*se_table[,7])), col="black", lty=2)
lines(res_table[,1], (res_table[,7] -(10.96*se_table[,7])), col="black", lty=2)
legend("topleft", c("Time series","Prediction"), col=c("red", "black"), lwd=1, pch=16)






Index=res_table[,c(1,8)]
colnames(Index)=c("Year","Index")
write.table(Index,"Index500800.csv",sep=";",row.names=F)