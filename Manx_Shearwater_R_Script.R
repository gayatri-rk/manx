##########################################################################################################
# Title - Examining vocal discrimination in Manx Shearwaters (Puffinus puffinus)

# Date 13/08/2026

# Note - final datasets for Lundy and Bardsey was used (birds with responses to  1 or 2 stimuli were excluded)


## clear working space
rm(list=ls())

# crt l to clear console 

# set working directory 
setwd()
getwd()

#load lundy data
manx_data <- read.csv("manx_data.csv")

library(dplyr)
library(tidyr)


## descriptive stats

## how many stimuli each bird saw
bird_count <- manx_data %>% group_by(bird_id) %>% summarise(count = n())


manx_data <- manx_data %>% distinct(bird_id, sex) %>% count(sex, name = "n_birds")

## mean latency for the 5 stimuli on Bardsey 
manx_stats_b <- manx_data %>% select(stimulus, latency, island) %>% filter(island == "Bardsey") %>% group_by(stimulus) %>% summarise(mean_latency =mean(latency, na.rm = TRUE)) 

## mean latency for the 5 stimuli on Lundy  
manx_stats_l <- manx_data %>% select(stimulus, latency, island) %>% filter(island == "Lundy") %>% group_by(stimulus) %>% summarise(mean_latency =mean(latency, na.rm = TRUE))


#range(manx_data$latency)

#range(manx_data$proportion)

#range(manx_data$duration)

## prepare data for statistical analysis 
manx_data <- manx_data %>% mutate( island = factor(island),
                                   sex = factor(sex),
                                   stimulus = factor(stimulus),
                                   order = as.numeric(order),
                                   bird_id = factor(bird_id),
                                   response_num = factor(response_num) )


### proportion - correction to adjust 0/ 1 values 

## manx_data <- manx_data %>%  mutate(prop_adj = (proportion * (n() - 1) + 0.5) /n())

## load required packages for analysis 

library("glmmTMB")
library("bbmle")
library("ggplot2")
library("ggeffects")
library("lme4")

######################## beta-GLMM-  interaction between island and stimulus type

glm_prop_island <- glmmTMB( prop_adj ~ island * stimulus + sex + order + (1| bird_id) ,
                            data = manx_data,
                            family = beta_family(link= "logit")
)


summary(glm_prop_island)

## model was not statistically significant 
# drop the interaction term and re-run the model

#########################################################################################

glm_prop <- glmmTMB( prop_adj ~ island + sex + stimulus + order + (1 | bird_id), 
                     data = manx_data,
                     family = beta_family(link ="logit"))

summary(glm_prop)

### model diagnosgtics for model fit and residual analysis 
library(DHARMa)
testDispersion(glm_prop_int)

sim_prop <- simulateResiduals(fittedModel = glm_prop_int, plot = FALSE)


testZeroInflation(sim_prop)
testUniformity(sim_prop)
testOutliers(sim_prop)


residuals(sim_prop)

residuals(sim_prop, quantileFunction = qnorm, outlierValues = c(-7, 7))


plot(sim_prop)

#### PLOT RESULTS

#get model predictions 
pred <- ggpredict(glm_prop, 
                  terms = c("stimulus" , "sex" ),  
                  type = "fixed" , bias_correction = TRUE)

# convert raw proportions to duration by multiplying by 100
pred$duration <- pred$predicted *100

pred$duration_low <- pred$conf.low *100

pred$duration_high <- pred$conf.high *100

## ensure the labels are correct 
pred$x <- factor(pred$x, 
                 levels = c("Control",
                            "Mate",
                            "Same-sex stranger",
                            "Self",
                            "Opposite sex stranger"),
                 labels = c("Control",
                            "Mate",
                            "Same-sex stranger" , 
                            "Self" ,
                            "Opposite sex stranger"))

## add significant stars to the plot for interpresation ease 
sig <- data.frame(
  stimulus = c("Control" , "Mate" , "Same-sex stranger"  , "Opposite sex stranger" , "Self"),
  label = c(" ", "** ", "***", "*", "***") , 
  y_pos = c(
    max(pred$duration_high[pred$x == "Control"]) + 5 ,
    max(pred$duration_high[pred$x == "Mate"]) +5 ,
    max(pred$duration_high[pred$x == "Same-sex stranger"]) +5 ,
    max(pred$duration_high[pred$x == "Opposite sex stranger"]) +5 ,
    max(pred$duration_high[pred$x == "Self"]) + 5 ))

## plot the results 

ggplot(pred, aes( x = x, y = duration, colour = group, group = group)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) + 
  geom_errorbar(
    aes(ymin = duration_low, ymax = duration_high),
    width = 0.2, 
    position = position_dodge(width = 0.4) 
  ) + 
  geom_text( data = sig,
             aes(x = stimulus, y= y_pos, label = label),
             colour = "black", 
             size = 6, 
             inherit.aes = FALSE) +
  labs(
    x = "Stimulus" ,
    y = "Response duration (%)" ,
    colour = "Sex"
  )  + theme_classic()  + theme(
    axis.title.x = element_text(size = 18, colour = "black"),
    axis.title.y  = element_text(size = 18, colour = "black"), 
    panel.background = element_rect(fill = "white" , colour = NA),
    plot.background = element_rect(fill = "white" , colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA , colour = "black", size = 1.5, linetype = "solid"),
    axis.text = element_text(colour = "black" , size = 15),
    legend.title = element_text(face = "bold" , size = 15, colour = "black") ,
    legend.text = element_text(size = 13, colour = "black")) +  theme(
      axis.title.x = element_text(vjust= -1.5) )



## save plot

ggsave("response_duration.png" , width = 25, height = 15, units = "cm" , dpi = 300)

## model interpretation
## intercept 
exp(-3.40056)
# # 0.03335459 log odds 

## back transforamtion
(exp(-3.40056))/ (1 + exp(-3.40056))
#  0.03227797 = 3.2% of the time birds responded to the female control

################# LATENCY


hist(manx_data$latency)
# data is not normally distributed - right skewed 

# log transform latency data 
hist(log(manx_data$latency))

## still not normally distributed 

## the data violates assumptions of a Gaussian GLMM

#### use GAMMA GLMM instead 

time_tmb <-glmmTMB(latency ~ island * stimulus + sex+ order + (1|bird_id),data=manx_data, family=Gamma(link = "log"))


summary(time_tmb)

## higher AIC (578.6) and interaction was n.s

######################################################################################

time_gamma <-glmmTMB(latency ~ island + stimulus + sex+ order + (1|bird_id),data=manx_data, family=Gamma(link = "log"))


summary(time_gamma)

## AIC = 573.1

### model diagnostics 
library(DHARMa)
testDispersion(time_gamma)

sim_time <- simulateResiduals(fittedModel = time_gamma, plot = FALSE)


testZeroInflation(sim_time)
testUniformity(sim_time)
testOutliers(sim_time)


residuals(sim_time)

residuals(sim_time, quantileFunction = qnorm, outlierValues = c(-7, 7))


plot(sim_time)


### no outliers or deviations from model assumptions 
## plot the results  

pred_3 <- ggpredict(time_gamma, 
                    terms = c("stimulus" , "sex"),  
                    type = "fixed" , bias_correction = FALSE)


pred_3$x <- factor(pred_3$x, 
                   levels = c("Control",
                              "Mate",
                              "Same-sex stranger",
                              "Self",
                              "Opposite sex stranger"),
                   labels = c("Control",
                              "Mate",
                              "Same-sex stranger" , 
                              "Self" ,
                              "Opposite sex stranger"))

## add ** to indicate significance level 

sig <- data.frame(
  stimulus = c("Control" , "Mate" , "Same-sex stranger"  , "Opposite sex stranger" , "Self"),
  label = c(" ", ". ", "**", ".", "***") , 
  y_pos = c(
    max(pred_3$conf.high[pred_3$x == "Control"]) + 5 ,
    max(pred_3$conf.high[pred_3$x == "Mate"]) + 6 ,
    max(pred_3$conf.high[pred_3$x == "Same-sex stranger"]) +5 ,
    max(pred_3$conf.high[pred_3$x == "Opposite sex stranger"]) + 6 ,
    max(pred_3$conf.high[pred_3$x == "Self"]) + 5 
  )
)


## plot the figure 
ggplot(pred_3, aes( x = x, y = predicted, colour = group, group = group)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) + 
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.3, 
    position = position_dodge(width = 0.4) 
  ) + 
  geom_text( data = sig,
             aes(x = stimulus, y= y_pos, label = label),
             colour = "black", 
             size = 8, 
             inherit.aes = FALSE) +
  labs(
    x = "Stimulus" ,
    y = "Response time (s)" ,
    colour = "Sex"
  )  + theme_classic()  + theme(
    axis.title.x = element_text(size = 13, colour = "black"),
    axis.title.y  = element_text(size = 13, colour = "black"), 
    panel.background = element_rect(fill = "white" , colour = NA),
    plot.background = element_rect(fill = "white" , colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA , colour = "black", size = 1.5, linetype = "solid"),
    axis.text = element_text(colour = "black" , size = 12),
    legend.title = element_text(face = "bold" , size = 12, colour = "black") ,
    legend.text = element_text(size = 11, colour = "black")) +  theme(
      axis.title.x = element_text(vjust= -1.5)
    )


## save plot

ggsave("latency.png" , width = 20, height = 13, units = "cm" , dpi = 300)

## model interpretation

## intercept 
exp(3.78)
## controls responded in 43.81604 s


################## descriptive stats for the results section

manx_sum <- manx_data %>% select(stimulus, sex, response_num, response)  %>% 
  group_by(stimulus, sex)  %>% summarise(count_0 = sum(response_num == 0),
                                         count_1 = sum(response_num == 1)) %>% group_by(stimulus, sex) 



ggplot(manx_sum, aes(fill=sex, y=count_1, x=stimulus)) + 
  geom_bar(position ='dodge', stat='identity') +
  
  labs(
    x = "Stimulus" ,
    y = "Response count" ,
    colour = "Sex"
  )  + theme_classic()  + theme(
    axis.title.x = element_text(size = 13, colour = "black"),
    axis.title.y  = element_text(size = 13, colour = "black"), 
    panel.background = element_rect(fill = "white" , colour = NA),
    plot.background = element_rect(fill = "white" , colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black" , size = 1),
    panel.border = element_rect(fill = NA , colour = "black", size = 1, linetype = "solid"),
    axis.text = element_text(colour = "black" , size = 12),
    legend.title = element_text(face = "bold" , size = 12, colour = "black") ,
    legend.text = element_text(size = 9, colour = "black")) 


## plotting the % of responses may be more useful

#load file 

manx_sum_stats <- read.csv("manx_summary.csv")

library(ggplot2)

ggplot(manx_sum_stats, aes(fill=sex, y=pct, x=stimulus)) + 
  geom_bar(position ='dodge', stat='identity') +
  
  labs(
    x = "Stimulus" ,
    y = "Response (%)" ,
    colour = "Sex"
  )  + theme_classic()  + theme(
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y  = element_text(size = 14, colour = "black"), 
    panel.background = element_rect(fill = "white" , colour = NA),
    plot.background = element_rect(fill = "white" , colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA , colour = "black", size = 1, linetype = "solid"),
    axis.text = element_text(colour = "black" , size = 12),
    legend.title = element_text(face = "bold" , size = 12, colour = "black") ,
    legend.text = element_text(size = 9, colour = "black"))  +  theme(
      axis.title.x = element_text(vjust= -1.5)
    ) + geom_text(aes(label =round(pct, 0) ,y=pct, x=stimulus),vjust=-.5,position=position_dodge(width =1))



## save plot

ggsave("manx_summary.png" , width = 25, height = 15, units = "cm" , dpi = 300)


## reference links used in this code

## https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf

## https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#binomial-data 

## https://www.statology.org/ggplot2-barplot-with-multiple-variables/

# https://www.geeksforgeeks.org/r-language/move-axis-labels-in-ggplot-in-r/

# https://luisdva.github.io/rstats/Labeled-barplots/
