##Full list of all potential models for E.coli analysis with Chisq comparisons of each model combination
##final_ecoli_analysis.r contains just the two models discussed in Figure 1 (Simple versus Final E.coli models)

d<-read.table("final_ecoli_binned_coding_table.txt",header=TRUE)

#Create dummy variable containing the R*logA value for each site
RlogA <- d$R * log10(d$abundance) 

sink("full_ecoli_analysis.txt", append = TRUE)
#Simple Model

glm_simple = glm(formula = E ~ R + RlogA + 0, family = poisson(link = identity), data = d)

summary(glm_simple)

###
#Condition Analysis, glm's containing condition as a variable for comparison
###
#Intercept by Condition, Single Slope
glm_cond_int = glm(formula = E ~ condition:R + RlogA + 0, family = poisson(link = identity), data = d)
summary(glm_cond_int)

#Intercept by Condition, Slope by Condition
full_cond = glm(formula = E ~ condition:R + condition:RlogA + 0, family = poisson(link = identity), data = d)
summary(full_cond)

###Chisq test comparisons
#Comparison of One Intercept/Slope with Different Intercepts by Condition, Single Slope
anova(glm_simple,glm_cond_int,test="Chisq")
#Comparison of One Intercept/Slope with Different Intercepts and Slopes by Condition
anova(glm_simple,full_cond,test="Chisq")
#Comparison of Separate Intercepts/One Slope with Different Intercepts and Slopes by Condition
anova(glm_cond_int,full_cond,test="Chisq")

sink()

###
#Substitution type Analysis, glm's containing subtype as a variable for comparison
###
#Some models unable to run without starting coefficients
#Start values for glms are simplified from the coefficients of a linear model of the corresponding formula, eg: lm(E ~ subtype:R + RlogA + 0)

sink("full_ecoli_analysis.txt", append = TRUE)

#Intercept by subtype, Single Slope

glm_sub_int = glm(formula = E ~ subtype:R + RlogA + 0, family = poisson(link = identity), data = d, start = c(-2.2e-06, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05))
summary(glm_sub_int)

#Intercept by subtype, Slope by subtype

full_sub = glm(formula = E ~ subtype:R + subtype:RlogA + 0, family = poisson(link = identity), data = d, start = c(1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, -2.2e-06, -2.2e-06, -2.2e-06, -2.2e-06, -2.2e-06, -2.2e-06, -2.2e-06, -2.2e-06, -2.2e-06, -2.2e-06, -2.2e-06))
summary(full_sub)

###Chisq test comparisons
#Comparison of One Intercept/Slope with Different Intercepts by Subtype, Single Slope
anova(glm_simple,glm_sub_int,test="Chisq")
#Comparison of One Intercept/Slope with Different Intercepts and Slopes by Subtype
anova(glm_simple,full_sub,test="Chisq")
#Comparison of Separate Intercepts/One Slope with Different Intercepts and Slopes by Subtype
anova(glm_sub_int,full_sub,test="Chisq")

sink()

###
#Combined Condition/Subtype Analysis
###
sink("full_ecoli_analysis.txt", append = TRUE)

#Condition/Subtype Intercept, Single Slope (Final Model Used for Analysis)

glm_set_slope = glm(formula = E ~ condition:R + subtype:R + RlogA + 0, family = poisson(link = identity), data = d,start = c(-2.2e-06,2e-05, 2e-05, 2e-05, 2e-05,1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05,1e-05, 1e-05))
summary(glm_set_slope)

###Combined intercept comparisons:
#Comparison of intercept by condition versus intercept by condition/subtype
anova(glm_cond_int,glm_set_slopetest="Chisq")
#Comparison of intercept by subtype versus intercept by condition/subtype
anova(glm_sub_int,glm_set_slope,test="Chisq")


##Multi-slope Analysis:

#Condition/Subtype Intercept, Slope by condition

glm_cond_slope = glm(formula = E ~ condition:R + subtype:R + condition:RlogA + 0, family = poisson(link = identity), data = d,start = c(2e-05, 2e-05, 2e-05, 2e-05,1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05,1e-05, 1e-05,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06))
summary(glm_cond_slope)

#Condition/Subtype Intercept, Slope by Subtype

glm_sub_slope = glm(formula = E ~ condition:R + subtype:R + subtype:RlogA + 0, family = poisson(link = identity), data = d,start = c(2e-05, 2e-05, 2e-05, 2e-05,1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05,1e-05, 1e-05,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06))
summary(glm_sub_slope)

#Full Possible Model, Condition and Subtype Intercept, Slope by Condtion and Subtype")

glm_full_slope = glm(formula = E ~ condition:R + subtype:R + condition:RlogA +subtype:RlogA + 0, family = poisson(link = identity), data = d,start = c(2e-05, 2e-05, 2e-05, 2e-05,1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05,1e-05, 1e-05,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06,-2.2e-06))
summary(glm_full_slope)

sink() 

###Multi-slope Comparisons:
sink("full_ecoli_analysis.txt", append = TRUE)

#Comparison of set slope(coniditon/subtype intercept) model to slope by condition
anova(glm_set_slope,glm_cond_slope,test="Chisq")
#Comparison of set slope(coniditon/subtype intercept) model to slope by subtype
anova(glm_set_slope,glm_sub_slope,test="Chisq")
#Comparison of set slope with Different Slopes by Condition and Subtype
anova(glm_set_slope,glm_full_slope,test="Chisq")

#Comparison of Slopes by Condition with Different Slopes by Condition and Subtype
anova(glm_cond_slope,glm_full_slope,test="Chisq")
#Comparison of Slopes by Subtype with Different Slopes by Condition and Subtype
anova(glm_sub_slope,glm_full_slope,test="Chisq")

sink()



