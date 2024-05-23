setwd("C:/Users/Michael/Documents/Tropical Biology")
dung<-read.csv("../Tropical Biology/Project/dungdata.csv")
str(dung)
summary(dung)
attach(dung)
m1<-lm(Volume.lost..cm.3. ~ treatment + canopy.cover + plot + site..LANE.ROUND.)
summary.lm(m1)
anova(m1, test="F")

m2<-lm(Volume.lost..cm.3. ~ treatment + canopy.cover + site..LANE.ROUND.)
plot(m2)
anova(m2, test="F")


m3<-lm(Volume.lost..cm.3. ~ treatment + site..LANE.ROUND.)
anova(m3, test="F")

plot(Volume.lost..cm.3. ~ treatment, notch=T)
summary(m3)

library(ggplot2)

qplot(treatment, Volume.lost..cm.3., geom="boxplot", notch=T)
