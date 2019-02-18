
d<-read.table('final_sacc_binned_table.txt',header=TRUE)

#Create dummy variable containing the R*logA value for each site
RlogA <- d$R * log10(d$abundance) 

sink("yeast_analysis_glm.txt", append = TRUE)

#Simple model of yeast data
simple_glm = glm(formula = E ~ R + RlogA + 0, family = poisson(link = identity), data = d)summary(simple_gro_glm)
summary(simple_glm)

#Subtype Model as a factor of intercept for yeast data
#No condition data is available for yeast so only substitution type is available for analysis
int_subtype_glm = glm(formula = E ~ subtype:R + RlogA + 0, family = poisson(link = identity), data = d)
summary(subtype_glm)
                                               

#Model with subtype as a factor of intercept and slope

#Start values derived from coefficients of linear model
#full_subtype_lm = lm(formula = E ~  subtype:R + subtype:RlogA + 0, data = d)
#summary(full_subtype_lm)
full_subtype_glm = glm(formula = E ~ subtype:R + subtype:RlogA + 0, family = poisson(link = identity), data = d, start = c(5e-07, 5e-07, 5e-07, 5e-07, 5e-07, 5e-07, 5e-07, 5e-07, 5e-07, 5e-07, 5e-07, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09))
summary(full_subtype_glm)
sink()   


#Chisq test comparison of potential models:

sink("yeast_analysis_glm.txt", append = TRUE)
#Comparison of Simple Model to Multiple Intercepts
anova(simple_glm, int_subtype_glm, test="Chisq")

#Comparison of Simple Model with Different Slopes by Subtype
anova(simple_glm, full_subtype_glm, test="Chisq")

#Comparison of Intercept by Subtype and Single Slope with Different Slopes by Subtype
anova(int_subtype_glm, full_subtype_glm, test="Chisq")
sink()
