##Anlysis of the Simple and Final E.coli models with Chisq comparison
#full_ecoli_glm_analysis.r contains every potential model tested (all combinations of condition/subtype) and comparison tests
##
d<-read.table("final_ecoli_binned_coding_table.txt",header=TRUE)

#Create dummy variable containing the R*logA value for each site
RlogA <- d$R * log10(d$abundance) #Simple Model

sink("final_ecoli_analysis.txt", append = TRUE)

###Simple Model
glm_simple = glm(formula = E ~ R + RlogA + 0, family = poisson(link = identity), data = d)
summary(glm_simple)

###Final Model (Intercept by Condition/Substitution Type, Single Slope)
glm_set_slope = glm(formula = E ~ condition:R + subtype:R + RlogA + 0, family = poisson(link = identity), data = d,start = c(-2.2e-06,2e-05, 2e-05, 2e-05, 2e-05,1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05,1e-05, 1e-05))
summary(glm_set_slope)

anova(glm_simple,glm_set_slope,test="Chisq")

sink()