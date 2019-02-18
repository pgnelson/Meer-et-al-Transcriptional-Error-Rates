#Generalized Linear Models used for the analysis and comparison of GroEL chaperone status on our exisiting model of mistranscription
d<-read.table('final_ecoli_binned_coding_table.txt',header=TRUE)

#Create dummy variable containing the R*logA value for each site
RlogA <- d$R * log10(d$abundance) 

sink("groEL_analysis_glm.txt", append = TRUE)

#Simple test of addition of groEL status to our "simple" model
simple_gro_glm = glm(formula = E ~ groEL:R + RlogA + 0, family = poisson(link = identity), data = d)
summary(simple_gro_glm)

#Start values derived from coefficients of linear model
#gro_lm = lm(formula = E ~ condition:R + subtype:R + groEL:R + RlogA + 0, data = d)
#summary(gro_lm)

#Addition of GroEL client status to our existing "Final Model"
final_gro_glm = glm(formula = E ~ condition:R + subtype:R + groEL:R + RlogA + 0, family = poisson(link = identity), data = d, start = c(-2e-06, 1e-04, 1e-04, 1e-04, 1e-04, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-04))
summary(final_gro_glm)

sink()


#Running our existing "Simple" and "Final" models for making comparisons to the GroEL models

sink("test_groEL_analysis_glm2.txt", append = TRUE)

#Simple Model
simple_glm = glm(formula = E ~ R + RlogA + 0, family = poisson(link = identity), data = d)
summary(simple_glm)

#Final Model
final_glm = glm(formula = E ~ condition:R + subtype:R + RlogA + 0, family = poisson(link = identity), data = d, start = c(-2.2e-06, 2e-05, 2e-05, 2e-05, 2e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05))
summary(final_glm)

#Chisq test to compare the model fit of the existing models and ones containing GroEL client status as a factor
anova(simple_glm, simple_gro_glm, test="Chisq")
anova(final_glm,final_gro_glm, test="Chisq")

sink()

