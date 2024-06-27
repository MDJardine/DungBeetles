## packages
library(ggplot2)
library(nlme)
library(lme4)
library(multcomp)
library(plotrix)
library(lsmeans)
library(aod)
library(Rmisc)


#### read in the data ####
setwd("C:/Users/user/Documents/GitHub/DungBeetles/data")
dung<-read.csv("dungdata.csv")
str(dung)
summary(dung)

## proportion transformation for canopy cover
dung$transcc <- asin(sqrt(dung$canopycover/100))

##########################################################
#### descriptive statistics ####


## getting means of removal per treatment, roung and plot
summarySE(dung, measurevar="vlost", groupvars=c("treatment"))
# Control = , fence, plate = 0.238

## how does this vary betweeen the rounds of the experiment?
summarySE(dung, measurevar="vlost", groupvars=c("treatment", "Round"))
# looks roughly consistant 

## how does this vary between the different locations?
# use tapply since not all conmbinations are used
Pmeans <- tapply(dung$vlost, dung$plot.lane, mean)
Pmeans

####################################################################
#### PART 1: comparing dung removal in all 3 treatments ####
## first box plot with all treatments separate

## wht does this look like?
boxplot(vlost ~ treatment, notch=T, data=dung)
## does seem to be a difference between the control and the two treatments

## does this vary between the rounds?
boxplot(vlost ~ Round, data=dung)
## yes but not toomuch or in any clear way

# ggplot
bp <- ggplot(dung, aes(x = treatment, y = vlost)) +
  theme_bw() + geom_boxplot() +
  scale_y_continuous(name="Volume of Dung Removed (L)") +
  theme(axis.title.y = element_text(face="bold", size=18), 
        axis.title.x = element_text(face= "bold", size=18), 
        axis.text.x=element_text(size=15), 
        axis.text.y=element_text(size=15)) + 
theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
      panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_x_discrete(name="Treatment", labels=c("Control", "Fence", "Plate"))
bp

## simple linear model as a test
# full model
testm1 <- lm(vlost ~ treatment + transcc + plot.lane + Round, data=dung)
summary(testm1)
anova (testm1, test="F")

#remove cc
testm2 <- lm(vlost ~ treatment + plot.lane + Round, data=dung)
summary(testm2)
anova (testm2, test="F")

# remove plot.lane
testm3 <- lm(vlost ~ treatment + Round, data=dung)
summary(testm3)
anova (testm3, test="F")

AIC(testm1, testm2, testm3)
anova(testm1, testm2, testm3)
##so it seems that the simplest model (testm1) is the best
## only treatment and round seem relavent

### mixed models for all 3 treatments
## want to do with nlme but can't work it out
# try now with lmer and try to fix later -> same results anyway just p value problem
model1 <- lmer(vlost ~ treatment*transcc + (1 | Round) + (1 | plot.lane), 
               data=dung)
summary(model1)
anova(model1)
model1

wald.test(b=fixef(model1), Sigma=vcov(model1), Terms = 3, df=70)
wald.test(b=fixef(model1), Sigma=vcov(model1), Terms = 4, df=70)


## interaction reallu not significant so remove
model2 <- lmer(vlost ~ treatment + transcc + (1 | Round) + (1 | plot.lane), 
               data=dung)
summary(model2)
anova(model2)

## canopy cover p-value
wald.test(b=fixef(model2), Sigma=vcov(model2), Terms = 2, df=72)

## doen't seen to work :(

## no effect of canopy cover so simplify again
model3 <- lmer(vlost ~ treatment + (1 | Round) + (1 | plot.lane), 
               data=dung)
summary(model3)
anova(model3)
## seems justified by the increasing F value but lets check with anova and AIC
anova(model1, model2, model3)
AIC(model1, model2, model3)
## better AIc, log-likelihood and justified Fvalue for model 3 = the best!!!

## confidence intervals
confint(model3, method='Wald')

##testing differences between the three treatments
# Tukey test
treatmentMCP <- mcp(treatment = 'Tukey')
tuk <- glht(model3, treatmentMCP)
confint(tuk)
tuk
par(mar=c(5.1, 6.2, 4.1, 2.1))
plot(tuk)
summary(tuk)
tuk
?glht

## finally lets chekc the model to ake sure its not rubbish
res<- resid(model3)
ran <- ranef(model3)
fit <- fitted(model3)
par(mfrow=c(1,4))
hist(ran$Round[,1])
hist(ran$plot.lane[,1])
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

### lsmeans from model 3
lsmeans(model3, "treatment")

wald.test(b=fixef(model3), Sigma=vcov(model3), Terms = 2, df=73)

#### PART 2: comparing dung removal in control and the combined treatments
## combining levels
#library(plyr)
#df2 <- data.frame(dung$Round, dung$plot, dung$treatment, dung$vlost)
#adder <- function(x){
#  x[x$dung.treatment == "fence",]$dung.vlost <- x[x$dung.treatment == "fence",]$dung.vlost + x[x$dung.treatment == "plate",]$dung.vlost
#  x <- x[!x$dung.treatment == "plate",]
#  levels(x$dung.treatment) <- droplevels(x$dung.treatment)
#  levels(x$dung.treatment) <- c("control", "combined")
#  return(data.frame(x$dung.treatment, x$dung.vlost))
#}
#
#df3 <- ddply(df2, "dung.plot", adder)
## problems with this - new data set that does the same thing
comb <- read.csv("F+P.csv")
summary(comb)
str(comb)

# clean data set to remove the outliers
ccomb = comb[which(comb$vlost < 0.7),] 

## mean F+P vlost
FPmeans <- tapply(ccomb$vlost, ccomb$treatment, mean)
FPmeans

#sd and se
FPsd <- tapply(ccomb$vlost, ccomb$treatment, sd)
FPsd
FPse <- 0.07716976/sqrt(29)
FPse

## second box plot with fence and plate merged
par(mfrow=c(1,1))
plot(vlost ~ treatment, notch=T, data=ccomb)
# ggplot stuffs
qplot(treatment, vlost, data=ccomb, notch=T, geom="boxplot")
bp2 <- ggplot(ccomb, aes(x = treatment, y = vlost)) +
  theme_bw() + geom_boxplot() +
  scale_y_continuous(name="Volume of Dung Removed (L)") +
  theme(axis.title.y = element_text(face="bold", size=18), 
        axis.title.x = element_text(face= "bold", size=18), 
        axis.text.x=element_text(size=15), 
        axis.text.y=element_text(size=15))
theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
      panel.background=element_blank(), axis.line=element_line(colour="black"))
bp2 
bp2 + scale_x_discrete(breaks=c("control", "F+P"), labels=c("control", "combined"))


### MODELS!!!
## already know that mixed models are need and to remove cc
## so straight on with the similar models

model4 <- lmer(vlost ~ treatment + (1|Round) + (1|plot.lane),
               data=ccomb)
summary(model4)
anova(model4)

## appears to be no effect of treatment so additive effect on vlost
## confidence intervals
confint(model4, method="Wald")

## lsmeans
lsmeans(model4, "treatment")

wald.test(b=fixef(model4), Sigma=vcov(model4), Terms = 2, df=39)

##testing differences between the three treatments
# Tukey test
combMCP <- mcp(treatment = 'Tukey')
combtuk <- glht(model4, combMCP)
confint(combtuk)
plot(combtuk)

## finally lets chekc the model to ake sure its not rubbish
combres<- resid(model4)
combran <- ranef(model4)
combfit <- fitted(model4)
par(mfrow=c(1,4))
hist(combran$Round[,1])
hist(combran$plot.lane[,1])
plot(combres ~ combfit)
qqnorm(combres, main='')
qqline(combres)
## although not great we have deleted 1/3 fo the data and it not too bad
### analysis done!!


#### fimally to create the combined plot for the paper
allplot <- read.csv("combined plot.csv")
str(allplot)

data = factor(allplot$treatment, c("fence", "plate", "combined", "control"))

cb <- ggplot(allplot, aes(x = data, y = vlost)) +
  theme_bw() + geom_boxplot() +
  scale_y_continuous(name="Volume of dung removed (Litres)") +
  theme(axis.title.y = element_text(face="bold", size=20), 
        axis.title.x = element_text(face="bold", size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_x_discrete(name="Treatment", labels=c("Fence", "Plate", "Combined", "Control")) + 
  coord_cartesian(ylim=c(0, 1))
cb
theme(axis.title.x = element_text(vjust=-0.5))



##################################################################
#### projections #####

## so we have some numbers about how much dung is removed in a set time frame

## can we therefore project what the ecoysstem service perfomed by these insects is per year?
## per area?  
## per elephant?

