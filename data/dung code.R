## packages
library(ggplot2)
library(lme4)
library(multcomp)
#library(plotrix)
#library(lsmeans)
#library(aod)
library(Rmisc)
library(RColorBrewer)
library(pbkrtest)

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
## yes but not too much or in any clear way

## histograms
hist(dung$vlost)
# its not overly normal, but we already know it has a bimodal distribution with the treatments
# kind of looks like two peaks?


## so we want to know if the treatment has an effect ion the volume of dung that is lost? 

## we used multiple rounds of the experiment
## and we had multiple sites
## both of these are random variables to control for
## can fix the random effects structure later


### mixed models for all 3 treatments
# try now with lmer and try to fix later -> same results anyway just p value problem
model1 <- lmer(vlost ~ treatment*transcc + (1 | Round) + (1 | plot.lane), data=dung)
## plots
res<- resid(model1)
ran <- ranef(model1)
fit <- fitted(model1)
par(mfrow=c(1,4))
hist(ran$Round[,1])
hist(ran$plot.lane[,1])
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
# result
summary(model1)
anova(model1)
model1
## so it looks like there is no interaction or anything to do with canopy cover 
# Random effects: round = plot.lane = 


#we'll get proper hypothesis testing from bootstrapping since we have random effects
### check with bootsrapping
treat3_1.1 <- lmer(vlost ~ 1 + (1 | Round) + (1 | plot.lane), data=dung)
treat3_2.1 <- lmer(vlost ~ treatment + (1 | Round) + (1 | plot.lane), data=dung)
treat3_3.1 <- lmer(vlost ~ transcc + (1 | Round) + (1 | plot.lane), data=dung)
treat3_4.1 <- lmer(vlost ~ treatment+transcc + (1 | Round) + (1 | plot.lane), data=dung)
treat3_5.1 <- lmer(vlost ~ treatment*transcc + (1 | Round) + (1 | plot.lane), data=dung)
#treatment
treat_test <- PBmodcomp(treat3_2.1, treat3_1.1)
treat_test
# p< 0.0001

# canopy cover
cc_test_3 <- PBmodcomp(treat3_3.1, treat3_1.1)
cc_test_3
# p = 0.05976

# interaction
intertest_3 <- PBmodcomp(treat3_5.1, treat3_4.1)
intertest_3
# p = 0.9944

# now we can plot this
ggplot(dung, aes(x = treatment, y = vlost, fill=treatment)) +
  theme_bw() + geom_boxplot() +
  scale_y_continuous(name="Volume of Dung Removed (L)") +
  theme(axis.title.y = element_text(face="bold", size=18), 
        axis.title.x = element_text(face= "bold", size=18), 
        axis.text.x=element_text(size=15), 
        axis.text.y=element_text(size=15)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black"))+
  scale_fill_brewer(palette = 'Dark2')
  scale_x_discrete(name="Treatment", labels=c("Control", "Fence", "Plate"))

## so its clear that the plate and fence treament have less volume lots

# our hypothesis would be that this is due to some guilds of dunf beetle beign blocked from acting on the dung

## but can we show the different betrween the levels statistically?

##testing differences between the three treatments
# Tukey test
treatmentMCP <- mcp(treatment = 'Tukey')
tuk <- glht(model1, treatmentMCP)
confint(tuk)
tuk
par(mar=c(5.1, 6.2, 4.1, 2.1))
plot(tuk)
summary(tuk)
## from this we can clearly see that the difference is between the controls nd either the plate or fence treatment
## not between those two treatments

### previously I did lsmeans to find the mean when account for rnadom effects
## as of 2024 this package no lobnger lseems to eb in effect - how to sort?
### lsmeans from model 3
lsmeans(model3, "treatment")



#### PART 2: comparing dung removal in control and the combined treatments ####
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

