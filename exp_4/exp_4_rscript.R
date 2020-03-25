########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 4 SCRIPT     #############
#############                              #############
########################################################
########################################################
########################################################
# load libraries:
library(lme4)
library(nlme)
library(boot)
library(car) 
library(reshape2)
library(ggplot2)
library(ez)
library(plyr)
library(ggsignif)
library(lsr)
library(sjmisc)
library(sjstats)
library(BayesFactor)
library(foreign)
library(dplyr)
library(lattice)
options(scipen=9999)


# load data
D = read.csv(file.choose(), header = TRUE)
D$X = NULL
D = D[c(1:32),]
D = as.data.frame(D[,c(1,2,3,7,5,9,4,8,6,10,11)])



# reshape the data
D_tall = reshape(D, varying = 3:10, v.names = "measure", 
                 timevar = "condition", idvar = "ID", 
                 new.row.names = 1:256, direction = "long")

# order data
D_tall = D_tall[order(D_tall$ID),]

# add q.type.cat column
D_tall$q.type.cat = as.factor(rep(c(1:2), each = 4, times = 32))


# set appropriate factor variables in "tall" data
D_tall$condition = as.factor(D_tall$condition)
D_tall$q.type = as.factor(D_tall$q.type)
D_tall$q.type.cat = as.factor(D_tall$q.type.cat)
D_tall$measure.2 = (100-D_tall$measure)


########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################

## NORMALITY CHECKS

par(mfrow=c(4,2)) 
for (ii in 1:8)  hist(D_tall$measure[D_tall$condition==ii], breaks=5)
par(mfrow=c(1,1)) 

# formal test of normality
shapiro.ps = rep(0,8)
for(i in 1:8) {
  shap.calc = shapiro.test(D_tall$measure[D_tall$condition==i])
  shapiro.ps[i] = shap.calc$p.value
}


## HOMOSKEDASTICITY CHECKS
# plot the boxplots
# perceptual
boxplot(D_tall$measure[D_tall$condition==c(1:4)]~D_tall$condition[D_tall$condition==c(1:4)])

# causal
boxplot(D_tall$measure[D_tall$condition==c(5:8)]~D_tall$condition[D_tall$condition==c(5:8)])

# formal test of equal variance
# perceptual
leveneTest(D_tall$measure[D_tall$condition==c(1:4)], as.factor(D_tall$condition[D_tall$condition==c(1:4)]), center=median) # used 'median' because it's a better measure of central tendency given the non-normality

# causal 
leveneTest(D_tall$measure[D_tall$condition==c(5:8)], as.factor(D_tall$condition[D_tall$condition==c(5:8)]), center=median) # used 'median' because it's a better measure of central tendency given the non-normality


## ASSUMPTION CHECK NOTES ##
# Given that there is evidence of non-normality, but NOT heteroskedasticity,
# non-parametric boostrapping and permuation
# testing will be used to estimate confidence intervals and for hypothesis checking


##############################
#### PRELIMINARY ANALYSIS ####
##############################
# analysis to determine effect of question type and 
lme.fit.prelim = lme(measure~q.type.cat+q.type+q.type.cat:q.type, random=~1|ID, data=D_tall)
anova.lme(lme.fit.prelim)


#########################
# FOLLOW UP COMPARISONS #
#########################

# PERCEPTUAL QUESTION FIRST: PERMUTATION
set.seed(2018)
b = rep(0,4000) 
for(i in 1:4000){
  y = sample(D_tall$measure.2, replace=TRUE)
  lm_1 = lmer(y[D_tall$q.type==0] ~ D_tall$q.type.cat[D_tall$q.type==0] + (1|D_tall$ID[D_tall$q.type==0]), data=D_tall) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lmer(D_tall$measure.2[D_tall$q.type==0] ~ D_tall$q.type.cat[D_tall$q.type==0] + (1|D_tall$ID[D_tall$q.type==0]), data=D_tall)
beta_actual = fixed.effects(lm.fit)[[2]]
c(beta_actual, sum(abs(b) > beta_actual)/4000, sum(abs(b) < beta_actual)/4000,
  sum(b > beta_actual)/4000, sum(b < beta_actual)/4000)
# -6.81250  1.00000  0.00000  0.95775  0.04200

    
    
# CAUSAL QUESTION FIRST: PERMUTATION
set.seed(2018)
b = rep(0,4000) 
for(i in 1:4000){
  y = sample(D_tall$measure.2, replace=TRUE)
  lm_1 = lmer(y[D_tall$q.type==1] ~ D_tall$q.type.cat[D_tall$q.type==1] + (1|D_tall$ID[D_tall$q.type==1]), data=D_tall) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lmer(D_tall$measure.2[D_tall$q.type==1] ~ D_tall$q.type.cat[D_tall$q.type==1] + (1|D_tall$ID[D_tall$q.type==1]), data=D_tall)
beta_actual = fixed.effects(lm.fit)[[2]]
c(beta_actual, sum(abs(b) > beta_actual)/4000, sum(abs(b) < beta_actual)/4000,
  sum(b > beta_actual)/4000, sum(b < beta_actual)/4000)



# CAUSAL AND PERCEPTUAL: BOOTSTRAPPED FUNCTION
boot_mean = as.data.frame(matrix(NA, nrow=4, ncol=3, byrow = FALSE)) 
k=1
for(i in 0:1){
  for(j in 1:2){
    set.seed(2018) 
    boot_func = function(data,b,formula, p){ 
      d= data[b,] 
      dif.1 =  mean(d$measure.2[d$q.type==i & d$q.type.cat==j], data=d) 
      return(dif.1)
    }
    GBGR_P_boot = boot(D_tall, boot_func, R=4000)
    boot_mean[k,] = c(GBGR_P_boot$t0, GBGR_P_boot$t0  + 1.96*-sd(GBGR_P_boot$t), 
                      GBGR_P_boot$t0  + 1.96*sd(GBGR_P_boot$t))
    k = k+1
  }
}

     V1       V2       V3
1 31.18750 26.60139 35.77361
2 24.37500 17.67544 31.07456
3 16.95312 11.22700 22.67925
4 20.78125 16.17914 25.38336  


#######################
#### MAIN ANALYSIS ####
#######################
# perceptual question
sub.percep = subset(D_tall, ! condition %in% c(5:8))
lme.fit.main.percep = lme(measure.2~as.factor(condition), random=~1|ID, data = sub.percep)
anova.lme(lme.fit.main.percep)
r2(lme.fit.main.percep)


# causal question
sub.causal = subset(D_tall, ! condition %in% c(1:4))
lme.fit.main.causal = lme(measure.2~condition, random=~1|ID, data = sub.causal)
anova.lme(lme.fit.main.causal)
r2(lme.fit.main.causal)



###################################
## FOLLOW-UP PLANNED COMPARISONS ##
###################################
#  permutation global function
perm_func = function(p,v){
  set.seed(2018)
  b = rep(0,4000) 
  for(i in 1:4000){
    x = factor(D_tall$condition, levels=c(p,v)) 
    y = sample(D_tall$measure.2, replace=TRUE)
    lm_1 = lmer(y ~ x + (1|ID), data=D_tall) 
    b[i] = fixed.effects(lm_1)[2]
  }
  
  lm.fit = lmer(D_tall$measure.2~factor(D_tall$condition, levels=c(p,v))+(1|ID), data=D_tall)
  beta_actual = fixed.effects(lm.fit)[[2]]
  c(beta_actual, sum(abs(b) > beta_actual)/4000, sum(abs(b) < beta_actual)/4000,
    sum(b > beta_actual)/4000, sum(b < beta_actual)/4000)
}



# bootstrap global function
boot_mean = as.data.frame(matrix(NA, nrow=8, ncol=3, byrow=TRUE))
for(i in 1:nrow(boot_mean)){ # want number of iterations to equal number of rows, especially because we're filling in by row
  set.seed(2018)
  boot_func = function(data,b,formula, p){ 
    d= data[b,] 
    x = d[d$condition==i,7]
    dif.1 =  mean(x, data=D_tall) 
    return(dif.1)
  }
  
  GBGR_P_boot = boot(D_tall, boot_func, R=4000) 
  boot_mean[i,] = c(GBGR_P_boot$t0, GBGR_P_boot$t0  + 1.96*-sd(GBGR_P_boot$t), 
                    GBGR_P_boot$t0  + 1.96*sd(GBGR_P_boot$t))
}




#####################################
## BOOT TESTS: PERCEPTUAL & CAUSAL ##
#####################################
V1        V2       V3
1 18.37500 10.144270 26.60573 # SFCF-P
2 27.65625 19.762143 35.55036 # SNCF-P
3 25.62500 18.196862 33.05314 # SFCN-P
4 24.62500 17.417065 31.83293 # SNCN-P
5 13.46875  6.681903 20.25560 # SFCF-C
6 23.81250 15.552217 32.07278 # SNCF-C
7 24.59375 16.485663 32.70184 # SFCN-C
8 28.43750 19.979172 36.89583 # SNCN-C


############################
## PERM TESTS: PERCEPTUAL ##
############################
# SFCF-P v. SNCF-P
perm_func(1,2)

# SFCF-P v. SFCN-P
perm_func(1,3)

# SFCF-P v. SNCN-P
perm_func(1,4)

# SNCF-P v. SFCN-P
perm_func(2,3)

# SNCF-P V. SNCN-P
perm_func(2,4)

# SFCN-P V. SNCN-P
perm_func(3,4)



############################
## PERM TESTS: CAUSAL ##
############################
# SFCF-C v. SNCF-C
perm_func(5,6)

# SFCF-C v. SFCN-C
perm_func(5,7)

# SFCF-C v. SNCN-C
perm_func(5,8)

# SNCF-C v. SFCN-C
perm_func(6,7)

# SNCF-C V. SNCN-C
perm_func(6,8)

# SFCN-C V. SNCN-C
perm_func(7,8)


################
# RELEVANT BFs for marginal differences in Perceptual main analysis #
################
# SFCF v SFCN
sub.percep.1 = subset(D_tall, ! condition %in% c(2,4:8))
# define the null and alternative models #
lm.null = lme(measure.2~1, random=~1|ID, data = sub.percep.1)
lm.alt = lme(measure.2~as.factor(condition), random=~1|ID, data = sub.percep.1)

#obtain BICs for the null and alternative models
null.bic = BIC(lm.null)
alt.bic = BIC(lm.alt)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01 = exp((alt.bic - null.bic)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10 = 1/BF01

# 5.162398


# SFCF v SNCN
sub.percep.2 = subset(D_tall, ! condition %in% c(2:3,5:8))
# define the null and alternative models #
lm.null = lme(measure.2~1, random=~1|ID, data = sub.percep.2)
lm.alt = lme(measure.2~as.factor(condition), random=~1|ID, data = sub.percep.2)

#obtain BICs for the null and alternative models
null.bic = BIC(lm.null)
alt.bic = BIC(lm.alt)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01 = exp((alt.bic - null.bic)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10 = 1/BF01

# 3.578657

#########################################
#### INDIVIDUAL DIFFERENCES ANALYSIS ####
#########################################
# INDIVIDUAL DIFFERENCE PLOTS FOR BOTH THE PERCEPTUAL AND CAUSAL QUESTIONS
condition_barplot = ggplot(F_tall, aes(condition, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 80)) +
  scale_fill_manual(values=c("#000000", "#999999")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# INDIVIDUAL DIFFERENCE PLOTS FOR PERCEPTUAL QUESTION
condition_barplot = ggplot(sub.percep, aes(condition, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 80)) +
  scale_fill_manual(values=c("#000000", "#999999")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials")


# INDIVIDUAL DIFFERENCE PLOTS FOR CAUSAL QUESTION
condition_barplot = ggplot(sub.causal, aes(condition, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 95)) +
  scale_fill_manual(values=c("#000000", "#999999")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials")


## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES ##
# install and load post-hoc chi-square test package (called follow up binomial tests)


# create contigency table with counts
chi.data = matrix(c(4,2,5,5,4,6,19,18,0,1),2)
dimnames(chi.data) = list(c("Perceptual", "Causal"), c("Color", "Shape", 
                                                       "All Familiar", "Other", "All Novel"))
## COUNTS
# Percep: C=4, S=5, AF=4, O=6,  L3N=13
# Causal: C=2, S=5, AF=6, O=7, L3N=12

# the chi data
              Color Shape All Familiar Other     L3N    Associative
Perceptual     4     5            4    6         13         17
Causal         2     5            6    7         12         18


# run chisq.test() on chi table
chi.test = chisq.test(chi.data, simulate.p.value = TRUE) 


## run post-hoc tests (where you test two particular cell means): binomial tests ##
# perceptual question

# Define a function to run post-hoc binomial tests
binom_func = function(x,y){
  bin_test = binom.test(x, x+y, p = 0.5, alternative = "two.sided")
  return(bin_test$p.value)
}


                    ##
## PERCEPTUAL CONDITION POST-HOC TESTS: ##
                    ##
# c v. s
binom_func(6,19)

# c v. af
binom_func(4,8)

# c v. o
binom_func(4,23)

#c v. an
binom_func(4,4)

                  ##
## PERCEPTUAL CONDITION POST-HOC TESTS: ##
                  ##
# c v. s
binom_func(2,7)

# c v. af
binom_func(4,8)

# c v. o
binom_func(4,23)

#c v. an
binom_func(4,4)






########################################################
########################################################
########################################################
#############                              #############
#############            Figures           #############
#############                              #############
########################################################
########################################################
########################################################

# Create 'F_tall' data frame to use for ggplot
F_tall = D_tall

# rename levels of 'condition' and 'q.type.cat' factors
F_tall$condition = revalue(x = as.factor(F_tall$condition), 
                           c("1" = "SFCF-P", "2"="SNCF-P", "3" = "SFCN-P", 
                             "4" = "SNCN-P", "5" = "SFCF-C", "6" = "SNCF-C",
                             "7" = "SFCN-C", "8" = "SNCN-C"))
F_tall$q.type.cat = revalue(x = as.factor(F_tall$q.type.cat), 
                            c("1" = "Perceptual Question", "2"="Causal Question"))


# OMNIBUS ANALYSIS FIGURE
condition_barplot = ggplot(F_tall, aes(condition, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  #facet_wrap(~q.type.cat, scales="free") + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  # PERCEPTUAL SIGNIFICANCE LINES
  geom_signif(comparisons = list(c("SFCF-P", "SNCF-P")), annotations=c("p < .05"), y_position = 37, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("SFCF-P", "SFCN-P")), annotations=c("p < .1"), y_position = 40, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("SFCF-P", "SNCN-P")), annotations=c("p < .15"), y_position = 43, tip_length = 0.00375) +
  # CAUSAL SIGNIFICANCE LINES
  geom_signif(comparisons = list(c("SFCF-C", "SNCF-C")), annotations=c("p < .05"), y_position = 36, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("SFCF-C", "SFCN-C")), annotations=c("p = .03"), y_position = 39, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("SFCF-C", "SNCN-C")), annotations=c("p < .01"), y_position = 42, tip_length = 0.00375) +
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 75)) +
  theme_classic() +
  scale_fill_manual(values=c("#000000", "#999999")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials")