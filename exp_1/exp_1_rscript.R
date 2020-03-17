########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 1 SCRIPT     #############
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
library(openxlsx)
options(scipen=9999)

# load data
D = read.csv(file.choose(), header = TRUE)

# reorder columns
D = as.data.frame(D[,c(1,2,3,9,5,7,4,10,6,8,11,12)])

# reshape the data
D_tall = reshape(D, varying = 3:10, v.names = "measure", timevar = "condition", 
                 idvar = "ID", 
                 new.row.names = 1:512, direction = "long")

# order data
D_tall = D_tall[order(D_tall$ID),]



# convert all categorical variables to 'factors'
D_tall$q.type.cat = as.factor(rep(c(0,1), each = 4, times = 64))
D_tall$q.type.cat = revalue(x = as.factor(D_tall$q.type.cat),
                            c("0"="Perceptual", "1"="Causal"))

D_tall$q.type = revalue(x = as.factor(D_tall$q.type), 
                                                c("0" = "Perceptual_First", 
                                                  "1"="Causal_First")) # this corresponds to whether participants
                                                                       # were asked the causal or perceptual question
                                                                       # first.

D_tall$exp = revalue(x = as.factor(D_tall$exp),
                       c("0" = "Kaily", "1" = "Marina"))


D_tall$group = revalue(x = as.factor(D_tall$group),
                     c("0" = "Red", "1" = "Blue"))


D_tall$test_trials = as.factor(rep(c(1:4), each = 1, 
                                   times = 128))
D_tall$test_trials = revalue(x = as.factor(D_tall$test_trials),
                       c("1" = "gBgR", "2" = "gRgB",
                         "3" = "gBRg", "4" = "gRBg"))



# drop unnecessary columns and reorder 'D_tall'
D_tall$row.names = NULL
D_tall$condition = NULL
D_tall = D_tall[,c(1,2,4,3,6,7,5)]
names(D_tall)
dim(D_tall)

# names of column in 'D_tall'
# "ID"         "exp"        "group"      "q.type"     "q.type.cat" "measure" 


########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################

## NORMALITY CHECKS
par(mfrow=c(4,2)) 
for(i in c("Perceptual","Causal")){
  for(j in c("gBgR","gRgB","gBRg","gRBg")){
    hist(D_tall$measure[D_tall$test_trials==j], breaks=5)
  }
}
par(mfrow=c(1,1)) 


# formal test of normality
shapiro.ps = rep(0,8)
for(i in c("Perceptual","Causal")){
  for(j in c("gBgR","gRgB","gBRg","gRBg")){
    shap.calc = shapiro.test(D_tall$measure[D_tall$test_trials==j])
    shapiro.ps[i] = shap.calc$p.value
  }
}


## HOMOSKEDASTICITY CHECKS
# plot the boxplots
init_box_plots = ggplot(D_tall, aes(x = test_trials, 
                                 y = measure)) + geom_boxplot() + facet_wrap(~q.type.cat)
init_box_plots

par(mfrow=c(4,2))
for (i in c("Perceptual","Causal")){
  for(j in c("gBgR","gRgB","gBRg","gRBg")){
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
# Despite the fact that there is no evidence of heteroskedasticity, because
# there is evidence of non-normality, non-parametric boostrapping and permuation
# testing will be used to estimate confidence intervals and for hypothesis checking


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
  for(i in 1:length(c)){
    set.seed(2020)
    sample_y = sample(data_1)
    vec_1 = sample_y[data_2==x]
    vec_2 = sample_y[data_2==y]
    dif = vec_1 - vec_2
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
lme.fit.prelim = lme(measure~(q.type+group)^2, random=~1|ID, data=D_tall)
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

# Follow-up analyses for rating type (i.e., perceptual v. causal rating) #
mean(D_tall$measure[D_tall$q.type.cat=="Perceptual"])
mean(D_tall$measure[D_tall$q.type.cat=="Causal"])
perm_func(D_tall$measure,D_tall$q.type.cat,"Causal","Perceptual")
BF_01 = ttestBF(D_tall$measure[D_tall$q.type.cat=="Perceptual"],
        D_tall$measure[D_tall$q.type.cat=="Causal"], paired=TRUE)
BF_01
  

# SD of causal question type
main_boot("Causal",5,7)

# SD of perceptual question type
main_boot("Perceptual",5,7)


# Follow-up analyses for test trial (i.e., gBgR, gRgB, gBRg,gRBg) #
test_trial_means = rep(0,4)
for(i in c("gBgR","gRgB","gBRg","gRBg")){
  calc = mean(D_tall$measure[D_tall$test_trials==i])
  test_trial_means[i] = calc
}
test_trial_means

# gBgR stats #
mean(D_tall$measure[D_tall$test_trials=="gBgR"])
main_boot("gBgR",6,7)

# gRgB stats #
mean(D_tall$measure[D_tall$test_trials=="gRgB"])
main_boot("gRgB",6,7)

# gBRg stats #
mean(D_tall$measure[D_tall$test_trials=="gBRg"])
main_boot("gBRg",6,7)

# gRBg stats #
mean(D_tall$measure[D_tall$test_trials=="gRBg"])
main_boot("gRBg",6,7)


# comparisons #

# gBgR v gBRg
perm_func(D_tall$measure,D_tall$test_trials,"gBgR","gBRg") # p < .0001
BF_02 = ttestBF(D_tall$measure[D_tall$test_trials=="gBgR"],
                D_tall$measure[D_tall$test_trials=="gBRg"], paired=TRUE)
BF_02 # BF10 = 31


# gBgR v gRBg
perm_func(D_tall$measure,D_tall$test_trials,"gBgR","gRBg") # p < .0001
BF_03 = ttestBF(D_tall$measure[D_tall$test_trials=="gBgR"],
                D_tall$measure[D_tall$test_trials=="gRBg"], paired=TRUE)
BF_03 # BF10 = 3



# gRgB v gBRg
perm_func(D_tall$measure,D_tall$test_trials,"gRgB","gBRg") # p < .0001
BF_04 = ttestBF(D_tall$measure[D_tall$test_trials=="gRgB"],
                D_tall$measure[D_tall$test_trials=="gBRg"], paired=TRUE)
BF_04 # BF10 = 291


# gRgB v gBRg
perm_func(D_tall$measure,D_tall$test_trials,"gRgB","gRBg") # p < .0001
BF_05 = ttestBF(D_tall$measure[D_tall$test_trials=="gRgB"],
                D_tall$measure[D_tall$test_trials=="gRBg"], paired=TRUE)
BF_05 # BF10 = 104


# gBgR v gRgB
perm_func(D_tall$measure,D_tall$test_trials,"gBgR","gRgB") # p < .0001
BF_06 = ttestBF(D_tall$measure[D_tall$test_trials=="gBgR"],
                D_tall$measure[D_tall$test_trials=="gRgB"], paired=TRUE)
BF_06 # BF10 = 31


# gBRg v gRBg
perm_func(D_tall$measure,D_tall$test_trials,"gBRg","gRBg") # p < .0001
BF_07 = ttestBF(D_tall$measure[D_tall$test_trials=="gBRg"],
                D_tall$measure[D_tall$test_trials=="gRBg"], paired=TRUE)
BF_07 # BF10 = 31


##################################
# INDIVIDUAL DIFFERENCE ANALYSIS #
##################################
# This analysis makes it possible to determine how many participants
# conformed with the independence relation predictions
B = D
B$ind_dif = rep(0,nrow(B))
for(i in 1:64){
  if(mean(c(B$gbgr.perceptual.rating[i],B$gbgr.causal.rating[i]))>mean(c(B$gbrg.perceptual.rating[i],B$gbrg.causal.rating[i]))&
     mean(c(B$gbgr.perceptual.rating[i],B$gbgr.causal.rating[i]))>mean(c(B$grbg.perceptual.rating[i],B$grbg.causal.rating[i]))&
     mean(c(B$grgb.perceptual.rating[i],B$grgb.causal.rating[i]))>mean(c(B$gbrg.perceptual.rating[i],B$gbrg.causal.rating[i]))&
     mean(c(B$grgb.perceptual.rating[i],B$grgb.causal.rating[i]))>mean(c(B$grbg.perceptual.rating[i],B$grbg.causal.rating[i]))){
    B$ind_dif[i] = 1
  }else{
    B$ind_dif[i] = 0
  }
}



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
  coord_cartesian(ylim=c(0, 100)) +
  scale_fill_manual(values=c("white", "white", "white", "white")) +
  theme_classic() +
  geom_signif(comparisons = list(c("gBgR", "gBRg")), 
              annotations=c("p < .0001"), y_position = 95, tip_length = 0.03) +
  geom_signif(comparisons = list(c("gBgR", "gRBg")), 
              annotations=c("p < .0001"), 
              y_position = 85, tip_length = 0.03)+
  theme(strip.text = element_text(colour = 'black', size = 12)) + # this changes the size and potentially weight of the facet labels
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials")