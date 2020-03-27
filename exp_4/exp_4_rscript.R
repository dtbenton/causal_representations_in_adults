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

# convert all categorical variables to 'factors'
D_tall$q.type.cat = as.factor(rep(c(0,1), each = 4, times = 32))
D_tall$q.type.cat = revalue(x = as.factor(D_tall$q.type.cat),
                            c("0"="Perceptual", "1"="Causal"))

D_tall$q.type = revalue(x = as.factor(D_tall$q.type), 
                        c("0" = "Perceptual_First", 
                          "1"="Causal_First")) # this corresponds to whether participants
                                               # were asked the causal or perceptual question
                                               # first.

D_tall$exp = revalue(x = as.factor(D_tall$exp),
                     c("0" = "Kaily", "1" = "Marina"))


D_tall$test_trials = as.factor(rep(c(1:4), each = 1, 
                                   times = 64))
D_tall$test_trials = revalue(x = as.factor(D_tall$test_trials),
                             c("1" = "SfCf", "2" = "SnCf",
                               "3" = "SfCn", "4" = "SnCn"))

# names D_tall
# "ID"         "exp"        "q.type"     "condition"  "measure"    "q.type.cat"

# drop unnecessary columns and reorder 'D_tall'
D_tall$row.names = NULL
D_tall$condition = NULL


# reorder columns of 'D_tall'
D_tall = D_tall[,c(1,2,3,5,6,4)]
names(D_tall)



########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################

## NORMALITY CHECKS
par(mfrow=c(4,2)) 
for(i in c("Perceptual","Causal")){
  for(j in c("SfCf","SnCf","SfCn","SnCn")){
    hist(D_tall$measure[D_tall$test_trials==j], breaks=5)
  }
}
par(mfrow=c(1,1)) 


# formal test of normality
shapiro.ps = as.data.frame(matrix(NA, nrow=4, ncol=1, byrow = TRUE)) 
k = 1
for(i in c("Perceptual","Causal")){
  for(j in c("SfCf","SnCf","SfCn","SnCn")){
    shap.calc = shap.calc = shapiro.test(D_tall$measure[D_tall$test_trials==j & D_tall$q.type.cat==i])
    shapiro.ps[k,] = shap.calc$p.value
    k = k+1
  }
}
shapiro.ps


## HOMOSKEDASTICITY CHECKS
# plot the boxplots
init_box_plots = ggplot(D_tall, aes(x = test_trials, 
                                    y = measure)) + geom_boxplot() + facet_wrap(~q.type.cat)
init_box_plots

par(mfrow=c(4,2))
for (i in c("Perceptual","Causal")){
  for(j in c("GfCf","GfCn","GnCf","GnCn")){
    boxplot(D_tall$measure[D_tall$q.type.cat==i & D_tall$test_trials==j], breaks=21)
  }
} 
par(mfrow=c(1,1))

# formal test of equal variance
# perceptual
levene_test_vec = rep(0,2)
for (i in c("Perceptual","Causal")){
  calc = leveneTest(D_tall$measure[D_tall$q.type.cat==i], 
                    as.factor(D_tall$test_trials[D_tall$q.type.cat==i]),
                    center = median)$"Pr(>F)"[[1]] 
  levene_test_vec[i] = calc
}
levene_test_vec


## ASSUMPTION CHECK NOTES ##
# There is both evidence of nonnormality and heteroskedasticity.



####################
####################
# GLOBAL FUNCTIONS #
####################
####################
main_boot = function(a,x,y){
  set.seed(2020)
  trees_estimate = function(data, indices) {
    d = data[indices,]
    mean_height = mean(d[d[,x]==a,y])
    return (mean_height)
  }
  compute_boot = boot(D_tall,trees_estimate,R=4000)
  return(c(compute_boot$t0,compute_boot$t0 + 1.96*c(-sd(compute_boot$t), 
                                                    sd(compute_boot$t))))
}

# example: main_boot("Perceptual",5,7)


perm_func = function(data_1,data_2,x,y){
  c = rep(0,10000)
  for(i in 1:10000){
    d = sample(data_1, replace=TRUE)
    a = d[data_2==x]
    b = d[data_2==y]
    dif = a - b
    c[i] = mean(dif)
  }
  
  bb_dif = mean(data_1[data_2==x]-data_1[data_2==y])
  c(mean(c),bb_dif,sum(abs(c) > bb_dif)/length(c),sum(abs(c) < bb_dif)/length(c), 
    mean(data_1[data_2==x]),mean(data_1[data_2==y]))
}



################################################
################################################
# WHETHER TO INCLUDE FIXED- VS. RANDOM EFFECTS #
################################################
################################################

# intercept only model with random effects of subjects
model_1 = lme(measure ~ 1, random=~1|ID,
              data=D_tall)

model_2 = lme(measure ~ 1, random=~1|exp/ID,
              data=D_tall)

model_3 = lm(measure~1, data = D_tall)

anova(model_1,model_2)
anova(model_1,model_3) # include random effects
anova(model_2,model_3)


# This analysis indicates that a model in which each deviation of each subject's intercept
# around the group intercept is the best model. Here, this is model_2. Thus, in all subsequent
# analyses random-effect intercepts will be included for subjects.


########################################################
########################################################
########################################################
#############                              #############
#############            Models            #############
#############                              #############
########################################################
########################################################
########################################################

##############################
#### PRELIMINARY ANALYSIS ####
##############################
# analysis to determine effect of question type and 
lme.fit.prelim = lme(measure~q.type, random=~1|ID, data=D_tall)
anova.lme(lme.fit.prelim)


## PRELIMINARY ANALYSIS NOTES
# no effect of question type (perceptual first vs causal first) or location (red first or blue 
# first in training sequence)


#######################
#### MAIN ANALYSIS ####
#######################
lme.fit.main = lme(measure~(test_trials+q.type.cat)^2, random=~1|ID, 
                   data = D_tall)
anova.lme(lme.fit.main)



# Follow-up analyses for test trial (i.e., gBgR, gRgB, gBRg,gRBg) #
test_trial_means = rep(0,4)
for(i in c("SfCf", "SnCf", "SfCn", "SnCn")){
  calc = mean(D_tall$measure[D_tall$test_trials==i])
  test_trial_means[i] = calc
}
test_trial_means

# SfCf stats #
mean(D_tall$measure[D_tall$test_trials=="SfCf"])
main_boot("SfCf",5,6)

# SnCf stats #
mean(D_tall$measure[D_tall$test_trials=="SnCf"])
main_boot("SnCf",5,6)

# SfCn stats #
mean(D_tall$measure[D_tall$test_trials=="SfCn"])
main_boot("SfCn",5,6)

# SnCn stats #
mean(D_tall$measure[D_tall$test_trials=="SnCn"])
main_boot("SnCn",5,6)


# comparisons #

# SfCf v SnCf
perm_func(D_tall$measure,D_tall$test_trials,"SfCf","SnCf") # p = .01
BF_02 = ttestBF(D_tall$measure[D_tall$test_trials=="SfCf"],
                D_tall$measure[D_tall$test_trials=="SnCf"], paired=TRUE)
BF_02 # BF10 = 16


# SfCf v SfCn
perm_func(D_tall$measure,D_tall$test_trials,"SfCf","SfCn") # p = .02
BF_03 = ttestBF(D_tall$measure[D_tall$test_trials=="SfCf"],
                D_tall$measure[D_tall$test_trials=="SfCn"], paired=TRUE)
BF_03 # BF10 = 13


# SfCf v SnCn
perm_func(D_tall$measure,D_tall$test_trials,"SfCf","SnCn") # p = .01
BF_04 = ttestBF(D_tall$measure[D_tall$test_trials=="SfCf"],
                D_tall$measure[D_tall$test_trials=="SnCn"], paired=TRUE)
BF_04 # BF10 = 11.86


# SnCf v SfCn
perm_func(D_tall$measure,D_tall$test_trials,"SnCf","SfCn") # p = 1
BF_05 = ttestBF(D_tall$measure[D_tall$test_trials=="SnCf"],
                D_tall$measure[D_tall$test_trials=="SfCn"], paired=TRUE)
BF_05 # BF10 = .14


# SnCf v SnCn
perm_func(D_tall$measure,D_tall$test_trials,"SnCf","SnCn") # p = .84
BF_06 = ttestBF(D_tall$measure[D_tall$test_trials=="SnCf"],
                D_tall$measure[D_tall$test_trials=="SnCn"], paired=TRUE)
BF_06 # BF10 = .14


# SfCn v SnCn
perm_func(D_tall$measure,D_tall$test_trials,"SfCn","SnCn") # p = .73
BF_06 = ttestBF(D_tall$measure[D_tall$test_trials=="SfCn"],
                D_tall$measure[D_tall$test_trials=="SnCn"], paired=TRUE)
BF_06 # BF10 = .17



########################################################
########################################################
########################################################
#############                              #############
#############            Figures           #############
#############                              #############
########################################################
########################################################
########################################################
condition_barplot = ggplot(D_tall, aes(test_trials, measure, fill = test_trials)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = c("#999999", "#E69F00", "#56B4E9","black")) + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 110)) +
  scale_fill_manual(values=c("white", "white", "white", "white")) +
  theme_classic() +
  geom_signif(comparisons = list(c("SfCf", "SnCn")), 
              annotations=c("p = .01, BF = 11.86"), 
              y_position = 105, tip_length = 0.003) +
  geom_signif(comparisons = list(c("SfCf", "SfCn")), 
              annotations=c("p = .02, BF = 13"), 
              y_position = 100, tip_length = 0.003) +
  geom_signif(comparisons = list(c("SfCf", "SnCf")), 
              annotations=c("p = .01, BF = 16"), 
              y_position = 95, tip_length = 0.003) +
  theme(legend.position="none") +
  theme(legend.position="none") +
  theme(strip.text = element_text(colour = 'black', size = 12)) + # this changes the size and potentially weight of the facet labels
  labs(x = "Test trials")