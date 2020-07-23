install.packages('tidyverse')
install.packages("car")
install.packages('rstatix')
install.packages('cowplot')
install.packages('ggpubr')
install.packages('datarium')

library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(datarium)

## Learning Data ############################################################################################
## Load Data ############################################################################################
mini_iris = iris[c(1, 51, 101), ]
gather(mini_iris, key = "flower_att", value = "measurement", Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
gather(mini_iris, key = "flower_att", value = "measurement", Sepal.Length, Sepal.Width, Petal.Length, Species)
gather(mini_iris, key = "flower_att", value = "measurement", Sepal.Length, Sepal.Width, Petal.Length )

anx_sub = sample_n_by( anxiety, group, size=1 )
gather( anx_sub, key="time", value="score", t1, t2, t3 )

anx_ren = gather( anxiety, key="time", value="score", t1, t2, t3 )
anx_ren %>%
  group_by( time, group ) %>%
  get_summary_stats( score, type='mean_sd')

bxp = ggboxplot( anx_ren, x="time", y="score", color="group", palette='jco')
bxp

anx_ren %>%
  group_by( time, group ) %>%
  identify_outliers( score )

anx_ren %>%
  group_by( time, group ) %>%
  identify_outliers( score )

anx_ren %>%
  group_by(time, group) %>%
  shapiro_test(score)

ggqqplot(anx_ren, "score", ggtheme = theme_bw()) +
  facet_grid(time ~ group)

anx_ren %>%
  group_by(time) %>%
  levene_test(score ~ group)

box_m(anx_ren[, "score", drop = FALSE], anx_ren$group)

res.aov <- anova_test( data = anx_ren, 
                       dv = score, wid = id,
                       between = group, within = time )
get_anova_table(res.aov)

# Effect of group at each time point
one.way <- anx_ren %>%
  group_by(time) %>%
  anova_test(dv = score, wid = id, between = group) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

# Pairwise comparisons between group levels
pwc <- anx_ren %>%
  group_by(time) %>%
  pairwise_t_test(score ~ group, p.adjust.method = "bonferroni")
pwc

one.way2 <- anx_ren %>%
  group_by(group) %>%
  anova_test(dv = score, wid = id, within = time) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

pwc2 <- anx_ren %>%
  group_by(group) %>%
  pairwise_t_test(
    score ~ time, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details
pwc2

anx_ren %>%
  pairwise_t_test( score ~ time, paired = TRUE, p.adjust.method = "bonferroni" )

anx_ren %>%
  pairwise_t_test( score ~ group, paired = FALSE, p.adjust.method = "bonferroni" )

pwc <- pwc %>% add_xy_position(x = "time")
pwc.filtered <- pwc %>% filter(time != "t1")
bxp + 
  stat_pvalue_manual(pwc.filtered, tip.length = 0, hide.ns = TRUE) +
  labs( subtitle = get_test_label(res.aov, detailed = TRUE),  caption = get_pwc_label(pwc) )

## PS Data ############################################################################################
## Load Data ############################################################################################
dta_tbl     = read.table('/home/ekaestner/gitrep/MMIL/CLASS/SCRATCH/data_id_no_missing.csv' , header=TRUE, sep=",")
dta_tbl_org = dta_tbl[1:63,]

## Format Data
dta_use_org = gather( dta_tbl_org, key="difficulty", value="lure_score", 
                      L1_1.old, L2_1.old, L3_1.old, L4_1.old, L5_1.old)

dta_use_new =  gather( dta_tbl, key="difficulty", value="lure_score", 
                       L1_1.old, L2_1.old, L3_1.old, L4_1.old, L5_1.old)

## Data summary
dta_use_org %>%
  group_by( GROUP, difficulty) %>%
  get_summary_stats(lure_score, type = "mean_sd")

## Boxplot visualization
bxp <- ggboxplot(
  dta_use_org, x = "difficulty", y = "lure_score",
  color = "GROUP", palette = "jco"
)
bxp

## Identify Outliers
dta_use_org %>%
  group_by( GROUP, difficulty) %>%
  identify_outliers(lure_score)

## Normality Assumption
dta_use_org %>%
  group_by(GROUP, difficulty) %>%
  shapiro_test(lure_score)

ggqqplot(dta_use_org, "lure_score", ggtheme = theme_bw()) +
  facet_grid( GROUP ~ difficulty )

## Homogeneity of variance
dta_use_org %>%
  group_by(difficulty) %>%
  levene_test(lure_score ~ GROUP)

dta_use_org %>%
  group_by(GROUP) %>%
  levene_test(lure_score ~ difficulty)

## Homogeneity of co-variance
box_m(dta_use_org[, "lure_score", drop = FALSE], dta_use_org$GROUP)

## Assumption of sphericity automatically checked

## ANOVA
res.aov <- anova_test(
  data = dta_use_org, dv = lure_score, wid = ID,
  between = GROUP, within = difficulty, type=3 )
get_anova_table(res.aov)
get_anova_table(res.aov, correction="GG")

# Effect of group at each time point
one.way <- dta_use_org %>%
  group_by(difficulty) %>%
  anova_test(dv = lure_score, wid = ID, between = GROUP) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

pwc <- dta_use_org %>%
  group_by(difficulty) %>%
  pairwise_t_test(lure_score ~ GROUP, p.adjust.method = "bonferroni")
pwc

# Visualization
pwc <- pwc %>% add_xy_position(x = "difficulty")
bxp <- ggboxplot(
  dta_use_org, x = "difficulty", y = "lure_score",
  color = "GROUP", palette = "jco"
)
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

# Effect of group at each group
one.way2 <- dta_use_org %>%
  group_by(GROUP) %>%
  anova_test(dv = lure_score, wid = ID, within = difficulty) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

pwc <- dta_use_org %>%
  group_by(GROUP) %>%
  pairwise_t_test(lure_score ~ difficulty, p.adjust.method = "bonferroni")
pwc

# Visualization
pwc <- pwc %>% add_xy_position(x = "GROUP")
bxp <- ggboxplot(
  dta_use_org, x = "GROUP", y = "lure_score",
  color = "difficulty", palette = "jco"
)
bxp + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

# Difficulty main effect
dta_use_org %>%
  pairwise_t_test(
    lure_score ~ difficulty, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )

# Group main effect
dta_use_org %>%
  pairwise_t_test(
    lure_score ~ GROUP, 
    p.adjust.method = "bonferroni"
  )

## Data New ############################################################################################
# Data summary
dta_use_new %>%
  group_by( GROUP, difficulty) %>%
  get_summary_stats(lure_score, type = "mean_sd")

# ANOVA
res.aov <- anova_test(
  data = dta_use_new, dv = lure_score, wid = ID,
  between = GROUP, within = difficulty, type=2 )
get_anova_table(res.aov)
get_anova_table(res.aov, correction="GG")

# Effect of group at each time point
one.way <- dta_use_new %>%
  group_by(difficulty) %>%
  anova_test(dv = lure_score, wid = ID, between = GROUP) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

pwc <- dta_use_new %>%
  group_by(difficulty) %>%
  pairwise_t_test(lure_score ~ GROUP, p.adjust.method = "bonferroni")
pwc

# Effect of timepoint at each group
one.way <- dta_use_new %>%
  group_by(GROUP) %>%
  anova_test(dv = lure_score, wid = ID, within = difficulty) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

# Group main effect
dta_use_new %>%
  pairwise_t_test(
    lure_score ~ GROUP, 
    p.adjust.method = "bonferroni"
  )



