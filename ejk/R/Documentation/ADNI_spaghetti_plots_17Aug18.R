library(splines)
library(boot)
library(MASS)
library(car)
library(rms)
library(sm)
library(readxl)
library(ggplot2)

setwd("C:/Users/Sean Hatton/Google Drive/Material and methods/Combat/ADNI")
dat <- read_excel("consolidated_17Aug18.xlsx")
attach(dat)

## Spaghetti plots
id <- factor(SubjID_correction)
Visit <- as.factor(Visit)
ggplot(aes(x=Visit,y=mean_cort_thick,colour=factor(Correction), group = id), data=dat) + 
  geom_line(size=1) +
  expand_limits(x = 0) +
  #facet_grid(. ~ Manufacturer) +
  #geom_smooth(se = FALSE, method = "lm") +
  xlab("Visit") + ylab("Mean cortical thickness") + 
  ggtitle("Change in mean cortical thickness (intra/interscanner)") +
  scale_colour_manual(breaks = c("15_orig", "30_orig", "15_30_harm", "30_15_harm"),
                      labels = c("Same 1.5T scanner", "Same 3.0T Scanner",
                                 "1.5T to 3.0T Harm","3.0T to 1.5T Harm"),
                      values = c("red", "blue", "darkgreen","lightblue"))

