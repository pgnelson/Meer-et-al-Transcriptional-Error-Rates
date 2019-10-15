##Anlysis of the Simple and Final E.coli models with Chisq comparison
#full_ecoli_glm_analysis.r contains every potential model tested (all combinations of condition/subtype) and comparison tests
##
final_ecoli_binned_coding_table <- read.delim("final_ecoli_binned_coding_table.txt")

#Create dummy variable containing the R*logA value for each site
final_ecoli_binned_coding_table$RlogA<-final_ecoli_binned_coding_table$R*log(final_ecoli_binned_coding_table$abundance)

#First test whether G->A has a different slope with protein abundance than the other 10 error types. (C->U has already been excluded)
starting_vals = c(-10^-6, 0, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5,10^-5,10^-5)
##Create dummy variable that containing R*logA if an error is G->A and 0 otherwise
final_ecoli_binned_coding_table$type_GA<-ifelse(final_ecoli_binned_coding_table$subtype=='GA', final_ecoli_binned_coding_table$RlogA, 0)
ga_glm<-glm(E ~subtype:R+condition:R + RlogA+type_GA+0,start = starting_vals, family = poisson(link = identity), data = final_ecoli_binned_coding_table)
starting_vals = c(-10^-6, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5,10^-5,10^-5)
pooled_glm<-glm(E ~subtype:R+condition:R + RlogA+0,start = starting_vals, family = poisson(link = identity), data = final_ecoli_binned_coding_table)
print("p-value associate with giving G->A it's own slope")
pchisq(2*(logLik(ga_glm)[1] - logLik(pooled_glm)[1]), df=1, lower.tail=FALSE)

#Next test whether the model is improved by giving the remaining 10 error types their own slopes.
starting_vals_sep<- c(10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5 , 10^-5, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6 , -10^-6)
separate_slopes_glm<-glm(E ~subtype:R+condition:R + subtype:RlogA+0,start = starting_vals_sep, family = poisson(link = identity), data = final_ecoli_binned_coding_table)
summary(separate_slopes_glm)
print("p-value associated with including separate a separate slope for G->A errors")
pchisq(2*(logLik(separate_slopes_glm)[1]-logLik(ga_glm)[1]), df=9, lower.tail=FALSE)

#testing Synonymous vs. non-synonymous
starting_vals = c(-10^-6, 0,10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5,10^-5,10^-5)
summary(ga_syn_glm<-glm(E ~subtype:R+condition:R+SNS:R + RlogA+type_GA+0,start = starting_vals, family = poisson(link = identity), data = final_ecoli_binned_coding_table))
print("p-value associated with including separate intercepts for synonymous and non-synonymous errors")
pchisq(2*(logLik(ga_syn_glm)[1]-logLik(ga_glm)[1]), df=9, lower.tail=FALSE)

#testing groEL
starting_vals = c(-10^-6, 0,10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5,10^-5,10^-5)
summary(ga_groel_glm<-glm(E ~subtype:R+condition:R +groEL:R + RlogA+type_GA+0,start = starting_vals, family = poisson(link = identity), data = final_ecoli_binned_coding_table))
print("p-value associated with including separate intercepts for groEL client proteins")
pchisq(2*(logLik(ga_groel_glm)[1]-logLik(ga_glm)[1]), df=9, lower.tail=FALSE)
