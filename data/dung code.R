## packages
library(ggplot2)
library(lme4)
library(multcomp)
#library(plotrix)
#library(lsmeans)
#library(aod)
#library(plyr)
library(Rmisc)
library(RColorBrewer)
library(pbkrtest)
library(emmeans)

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
# Random effects account for: round = 20.96%, plot.lane = 14.26% of the residual variance


#we'll get proper hypothesis testing from bootstrapping since we have random effects
### check with bootstrapping
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
  scale_fill_brewer(palette = 'Dark2')+
  scale_x_discrete(name="Treatment", labels=c("Control", "Fence", "Plate"))

## so its clear that the plate and fence treatment have less volume lots

# our hypothesis would be that this is due to some guilds of dung beetle being blocked from acting on the dung

## but can we show the different between the levels statistically?

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

### previously I did lsmeans but this has been replaced by emmeans
em_mean_treat = emmeans(model1, ~ treatment | transcc)
em_mean_treat

# plot
plot(em_mean_treat, comparisons = T)


#### PART 2: comparing dung removal in control and the combined treatments ####
## combining levels
# since the control is much higher than the two treatments we want to see if the volume
# lost there is much different from when we add the two treatments together


## this is included in the dataet 'F+P.csv'
comb <- read.csv("F+P.csv")
summary(comb)
str(comb)

# have a look
hist(comb$vlost)

# some large outliers where we think things were disripted by guinea fowl
ccomb = comb[which(comb$vlost < 0.7),] 

## mean F+P vlost
summarySE(ccomb, measurevar="vlost", groupvars=c("treatment"))
## looks very similar now - possibly no difference

## build a linear mixed model similar to before where we look at the volume lost per treatment
# these won't include the canopy cover since it was impossible to include this when combing measures
# unless we took an average soe not sure this is best

comb_model <- lmer(vlost ~ treatment + (1|Round) + (1|plot.lane), data=ccomb)
#plot
res<- resid(comb_model)
ran <- ranef(comb_model)
fit <- fitted(comb_model)
par(mfrow=c(1,4))
hist(ran$Round[,1])
hist(ran$plot.lane[,1])
plot(res ~ combfit)
qqnorm(res, main='')
qqline(res)
# looks ok - now have a look
summary(comb_model)
anova(comb_model)
comb_model
# Random effects account for: round = 11.3%, plot.lane = 17.66% of the residual variance


# doesn't look like there is a difference between treatments this time
## check with bootstrapping
#we'll get proper hypothesis testing from bootstrapping since we have random effects
### check with bootstrapping
treat2_1.1 <- lmer(vlost ~ 1 + (1 | Round) + (1 | plot.lane), data=ccomb)
treat2_2.1 <- lmer(vlost ~ treatment + (1|Round) + (1|plot.lane), data=ccomb)
#treatment
treat_test_comb <- PBmodcomp(treat2_2.1, treat2_1.1)
treat_test_comb
# p = 0.6811

## Combining the plate and the fence results in a similar volume of dung lost by each pile than when neither are present
## indicates that the guilds could be salutatory and that dung removed by one isn't utilised by the other
## if this was the case then the volume of dung lost would have been higehr than in the control


## Now plot this to illustrate this result
# PLOT
ggplot(ccomb, aes(x = treatment, y = vlost, fill=treatment))+
  theme_bw() + geom_boxplot()+
  scale_y_continuous(name="Volume of Dung Removed (L)")+
  scale_x_discrete(name="Treatment", breaks=c("control", "F+P"), labels=c("Control", "Combined"))+
  theme(axis.title.y = element_text(face="bold", size=18), 
        axis.title.x = element_text(face= "bold", size=18), 
        axis.text.x=element_text(size=15), 
        axis.text.y=element_text(size=15))+
  scale_fill_manual(values=c("#1B9E77", "#E7298A"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black"))
  

#### finally to create the combined plot with all four included.

# again I've done this with excel - should write a proper method in plyr for this and above
allplot <- read.csv("combined plot.csv")
str(allplot)
allplot

# take out disrupted volumes
allplot = allplot[which(allplot$vlost < 0.7),] 

# colours need to be re-ordered to be consistent
ggplot(allplot, aes(x = treatment, y = vlost, fill=treatment))+
  theme_bw() + geom_boxplot() +
  scale_y_continuous(name="Volume of dung removed (L)")+
  theme(axis.title.y = element_text(face="bold", size=20), 
        axis.title.x = element_text(face="bold", size=20, vjust=-0.5), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black"))+ 
  scale_fill_manual(values=c("#1B9E77", "#E7298A", "#D95F02", "#7570B3"))+
   scale_x_discrete(name="Treatment", labels=c("Control", "Fence + Plate", "Fence", "Plate"))


### Done

## look at projections in another script