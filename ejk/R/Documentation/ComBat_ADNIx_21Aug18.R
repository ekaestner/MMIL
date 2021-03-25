library(splines)
library(boot)
library(MASS)
library(car)
library(rms)
library(sm)
library(readxl)
library(ggplot2)
setwd("C:/Users/Sean Hatton/Google Drive/Material and methods/Combat/ADNIx")
#setwd("~/Desktop/Combat/gm")

###source ComBat scripts
source("../scripts/utils.R")
source("../scripts/combat.R")
source("../scripts/plot.R")

datorig <- read_excel("MRI_all_ADNI.xlsx")
dat.v1_15 <- datorig[which(datorig$Visit == "V1_15"),]
dat.v2_15 <- datorig[which(datorig$Visit == "V2_15"),]
dat.v1_30 <- datorig[which(datorig$Visit == "V1_30"),]
dat.v2_30 <- datorig[which(datorig$Visit == "V2_30"),]
dat.bsc <- rbind(dat.v1_15,dat.v1_30)
#dat.bsc <- rbind(dat.v1_30,dat.v2_15)
#complete.cases(dat.bsc)
attach(dat.bsc)
#detach(dat.bsc)
Visit[Visit == "V1_15"] <-1
Visit[Visit == "V1_30"] <-2
Visit <- as.numeric(Visit)
dat.new<-rbind(`cort_thick-ctx-lh-bankssts`,`cort_thick-ctx-lh-caudalanteriorcingulate`,
               `cort_thick-ctx-lh-caudalmiddlefrontal`,`cort_thick-ctx-lh-cuneus`,
               `cort_thick-ctx-lh-entorhinal`,`cort_thick-ctx-lh-fusiform`,
               `cort_thick-ctx-lh-inferiorparietal`,`cort_thick-ctx-lh-inferiortemporal`,
               `cort_thick-ctx-lh-isthmuscingulate`,`cort_thick-ctx-lh-lateraloccipital`,
               `cort_thick-ctx-lh-lateralorbitofrontal`,`cort_thick-ctx-lh-lingual`,
               `cort_thick-ctx-lh-medialorbitofrontal`,`cort_thick-ctx-lh-middletemporal`,
               `cort_thick-ctx-lh-parahippocampal`,`cort_thick-ctx-lh-paracentral`,
               `cort_thick-ctx-lh-parsopercularis`,`cort_thick-ctx-lh-parsorbitalis`,
               `cort_thick-ctx-lh-parstriangularis`,`cort_thick-ctx-lh-pericalcarine`,
               `cort_thick-ctx-lh-postcentral`,`cort_thick-ctx-lh-posteriorcingulate`,
               `cort_thick-ctx-lh-precentral`,`cort_thick-ctx-lh-precuneus`,
               `cort_thick-ctx-lh-rostralanteriorcingulate`,
               `cort_thick-ctx-lh-rostralmiddlefrontal`,`cort_thick-ctx-lh-superiorfrontal`,
               `cort_thick-ctx-lh-superiorparietal`,`cort_thick-ctx-lh-superiortemporal`,
               `cort_thick-ctx-lh-supramarginal`,`cort_thick-ctx-lh-frontalpole`,
               `cort_thick-ctx-lh-temporalpole`,`cort_thick-ctx-lh-transversetemporal`,
               `cort_thick-ctx-lh-insula`,`cort_thick-ctx-rh-bankssts`,
               `cort_thick-ctx-rh-caudalanteriorcingulate`,
               `cort_thick-ctx-rh-caudalmiddlefrontal`,`cort_thick-ctx-rh-cuneus`,
               `cort_thick-ctx-rh-entorhinal`,`cort_thick-ctx-rh-fusiform`,
               `cort_thick-ctx-rh-inferiorparietal`,`cort_thick-ctx-rh-inferiortemporal`,
               `cort_thick-ctx-rh-isthmuscingulate`,`cort_thick-ctx-rh-lateraloccipital`,
               `cort_thick-ctx-rh-lateralorbitofrontal`,`cort_thick-ctx-rh-lingual`,
               `cort_thick-ctx-rh-medialorbitofrontal`,`cort_thick-ctx-rh-middletemporal`,
               `cort_thick-ctx-rh-parahippocampal`,`cort_thick-ctx-rh-paracentral`,
               `cort_thick-ctx-rh-parsopercularis`,`cort_thick-ctx-rh-parsorbitalis`,
               `cort_thick-ctx-rh-parstriangularis`,`cort_thick-ctx-rh-pericalcarine`,
               `cort_thick-ctx-rh-postcentral`,`cort_thick-ctx-rh-posteriorcingulate`,
               `cort_thick-ctx-rh-precentral`,`cort_thick-ctx-rh-precuneus`,
               `cort_thick-ctx-rh-rostralanteriorcingulate`,`cort_thick-ctx-rh-rostralmiddlefrontal`,
               `cort_thick-ctx-rh-superiorfrontal`,`cort_thick-ctx-rh-superiorparietal`,
               `cort_thick-ctx-rh-superiortemporal`,`cort_thick-ctx-rh-supramarginal`,
               `cort_thick-ctx-rh-frontalpole`,`cort_thick-ctx-rh-temporalpole`,
               `cort_thick-ctx-rh-transversetemporal`,`cort_thick-ctx-rh-insula`,
               `cort_thick-ctx-lh-mean`,`cort_thick-ctx-rh-mean`,`cort_thick-ctx-mean`) #create matrix, with subjects as columns

#dat.new <- rbind(cort_thickctxlhmean,cort_thickctxrhmean,cort_thickctxmean)

combat.test<-combat(dat=dat.new, batch=Visit)

export <- data.frame(SubjID = dat.bsc$SubjID,
                VisitID = dat.bsc$VisitID,
                Visit = dat.bsc$Visit,
                age = Days_from_baseline)
temp <- data.frame(t(combat.test$dat.combat))
export <- data.frame(export,temp)
write.csv(export, file = "ComBat_harmonized.csv")
rm(temp)

ggplot(export, aes(x=Visit, y=cort_thick.ctx.mean, fill=Visit)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_brewer(palette="Set1") +
  ylim(2.1, 2.8) +
  labs(title="Mean cortical thickness (ComBat)", x="Visit", y="Thickness (mm)")

ggplot(dat.bsc, aes(x=Visit, y=`cort_thick-ctx-mean`, fill=Visit)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_brewer(palette="Set1") +
  ylim(2.7, 3.6) +
  labs(title="Mean cortical thickness", x="Visit", y="Thickness (mm)")


#plot(age_v2,combat.test$dat.combat[1,])
#plot(age_v2,FA_AllFibers)
predictor<-age
dv<-combat.test$dat.combat[69,] # cort_thick.ctx.mean
dv2<-cort_thickctxmean
#dv<-combat.test$dat.combat[5,] # cort_thickctxlhentorhinal
#dv2<-cort_thickctxlhentorhinal
#dv<-combat.test$dat.combat[1,] # subcort_volWholeBrain
#dv2<-subcort_volWholeBrain


pred1<-lm(dv~bs(predictor,degree=2, df = 2), na.action = na.omit) 
pred1.lm<-lm(dv~predictor, na.action = na.omit) 
pred2<-lm(dv2~bs(predictor,degree=2, df = 2), na.action = na.omit) 
pred2.lm<-lm(dv2~predictor, na.action = na.omit) 


##make the figures
tiff("combat_meancortthick.tiff", units = 'in', width = 8, height = 6, res = 200)

myPredict <- predict(pred1, interval="confidence", level = 0.95, na.action = na.omit) #For each value of x, get the value of y estimated by the model, and the confidence interval around this value.
x<-pred1.lm$model[,2]
dv<-pred1$model[,1]
plot(x,dv, col = rgb(0, .5, .6, .3), xlab = "Age in Years", main = "Mean Cortical Thickness (Raw Data = Red; ComBat Corrected Data = Teal)", ylab = "Volume", pch=16, cex.lab = 1.5, cex.axis = 1.5)
ix <- sort(x,index.return=T, na.last = NA)$ix # sort x values
lines(x[ix], myPredict[ix, 1], col = rgb(0, .5, .6, 1), lwd=2 ) # place the curve on the plot 
polygon(c(rev(x[ix]), x[ix]), c(rev(myPredict[ ix,3]), myPredict[ ix,2]), col = rgb(0, .5, .6, 0.2) , border = NA) #create the confidence interval

myPredict <- predict(pred2, interval="confidence", level = 0.95, na.action = na.omit) #For each value of x, get the value of y estimated by the model, and the confidence interval around this value.
x<-pred2.lm$model[,2]
dv<-pred2$model[,1]
points(x,dv, col = rgb(.9, 0, 0, .3), pch=16)
ix <- sort(x,index.return=T, na.last = NA)$ix # sort x values
lines(x[ix], myPredict[ix , 1],  col = rgb(.9, 0, 0, 1), lwd=2 ) # place the curve on the plot 
polygon(c(rev(x[ix]), x[ix]), c(rev(myPredict[ ix,3]), myPredict[ ix,2]), col = rgb(.9, 0, 0,0.2) , border = NA) #create the confidence interval

dev.off()

tiff("site_diff_before.tiff", units = 'in', width = 8, height = 6, res = 200)


sm.density.compare(cort_thickctxlhbankssts, Scanner, xlab="Thickness (mm)")
title(main="Site Differences in Cortical Thickness Before Correction")
dev.off()

tiff("site_diff_after.tiff", units = 'in', width = 8, height = 6, res = 200)

sm.density.compare(combat.test$dat.combat[1,], Scanner, xlab="Thickness (mm)")
title(main="Site Differences in Cortical Thickness After Correction")
dev.off()

# Overall slope
ggplot(aes(x=age,y=cort_thickctxmean), data=dat.bsc) + 
  geom_point(aes(color=factor(Scanner))) +
  geom_smooth(colour="black") +
  xlab("Age (years)") + ylab("Mean cortical thickness") + ggtitle("Mean cortical thickness (uncorrected) across 10 years") + 
  scale_colour_manual(breaks = c("1", "2", "3", "4", "5"),
                      labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                                 "VETSA2 UCSD", "VETSA3 UCSD"),
                      values = c("lightblue", "blue", "red", 
                                 "orange", "darkgreen"))

ggplot(aes(x=age,y=cort_thickctxmean), data=export) + 
  geom_point(aes(color=factor(Scanner))) +
  geom_smooth(colour="black") +
  xlab("Age (years)") + ylab("Mean cortical thickness") + ggtitle("Mean cortical thickness (corrected) across 10 years") + 
  scale_colour_manual(breaks = c("1", "2", "3", "4", "5"),
                      labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                                 "VETSA2 UCSD", "VETSA3 UCSD"),
                      values = c("lightblue", "blue", "red", 
                                 "orange", "darkgreen"))

# Slope per scanner
ggplot(aes(x=age,y=cort_thickctxmean,colour=factor(Scanner)), data=dat.bsc) + 
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Age (years)") + ylab("Mean cortical thickness") + ggtitle("Mean cortical thickness (uncorrected) across 10 years") + 
  scale_colour_manual(breaks = c("1", "2", "3", "4", "5"),
                      labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                                 "VETSA2 UCSD", "VETSA3 UCSD"),
                      values = c("lightblue", "blue", "red", 
                                 "orange", "darkgreen"))

ggplot(aes(x=age,y=cort_thickctxmean,colour=factor(Scanner)), data=export) + 
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Age (years)") + ylab("Mean cortical thickness") + ggtitle("Mean cortical thickness (corrected) across 10 years") + 
  scale_colour_manual(breaks = c("1", "2", "3", "4", "5"),
                      labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                                 "VETSA2 UCSD", "VETSA3 UCSD"),
                      values = c("lightblue", "blue", "red", 
                                 "orange", "darkgreen"))
## Spaghetti plots
id <- factor(SubjID)
ggplot(aes(x=age,y=cort_thickctxmean,colour=factor(Scanner), group = id), data=dat.bsc) + 
  geom_line() +
  #geom_smooth(se = FALSE, method = "lm") +
  xlab("Age (years)") + ylab("Mean cortical thickness") + ggtitle("Mean cortical thickness (uncorrected) across 10 years") + 
  scale_colour_manual(breaks = c("1", "2", "3", "4", "5"),
                      labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                                 "VETSA2 UCSD", "VETSA3 UCSD"),
                      values = c("lightblue", "blue", "red", 
                                 "orange", "darkgreen"))

ggplot(aes(x=age,y=cort_thickctxmean,colour=factor(Scanner), group = id), data=export) + 
  geom_line() +
  #geom_smooth(se = FALSE, method = "lm") +
  xlab("Age (years)") + ylab("Mean cortical thickness") + ggtitle("Mean cortical thickness (corrected) across 10 years") + 
  scale_colour_manual(breaks = c("1", "2", "3", "4", "5"),
                      labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                                 "VETSA2 UCSD", "VETSA3 UCSD"),
                      values = c("lightblue", "blue", "red", 
                                 "orange", "darkgreen"))




ggplot(aes(x=age,y=cort_thickctxlhentorhinal), data=export) + 
  geom_point(aes(color=factor(Scanner))) +
  geom_smooth(colour="black") +
  xlab("Age (years)") + ylab("Mean cortical thickness") + ggtitle("Left entorhinal cortical thickness (corrected) across 10 years") + 
  scale_colour_manual(breaks = c("1", "2", "3", "4", "5"),
                      labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                                 "VETSA2 UCSD", "VETSA3 UCSD"),
                      values = c("lightblue", "blue", "red", 
                                 "orange", "darkgreen"))

dat.bsc$Scanner <- factor(dat.bsc$Scanner)
ggplot(dat.bsc, aes(x=Scanner, y=cort_thickctxmean, fill=Scanner)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +
  labs(title="Mean cortical thickness (uncorrected)", x="Scanner", y="Thickness (mm)") +
  ylim(1.9, 2.8) +
  scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                      labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                                 "VETSA2 UCSD", "VETSA3 UCSD"),
                      values = c("lightblue", "blue", "red", 
                                 "orange", "darkgreen")) +
  scale_x_discrete(labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                              "VETSA2 UCSD", "VETSA3 UCSD"))

export$Scanner <- factor(export$Scanner)
ggplot(export, aes(x=Scanner, y=cort_thickctxmean, fill=Scanner)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +
  labs(title="Mean cortical thickness (corrected)", x="Scanner", y="Thickness (mm)") +
  ylim(1.9, 2.8) +
  scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                    labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                               "VETSA2 UCSD", "VETSA3 UCSD"),
                    values = c("lightblue", "blue", "red", 
                               "orange", "darkgreen")) +
  scale_x_discrete(labels = c("VETSA1 BU", "VETSA1 UCSD", "VETSA2 BU",
                              "VETSA2 UCSD", "VETSA3 UCSD"))




