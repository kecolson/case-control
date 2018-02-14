################################################################################################
# NAME: Old Code Dump
# AUTHORS: Ellie Matthay, Catherine Li, Chris Rowe
# DATE STARTED: 02/14/2018 
# PURPOSE: Place to store past code for case-control simulation project
#          
# UPDATES: 02/14/2018: CL added logistic outcome mechanism code
#                      CL added code from calculating true measures of association and dispersion
#
################################################################################################

## LOGISTIC OUTCOME MECHANISM CODE
# Generate Outcome as Function of Exposure and Covariates - Socio-Demographic Patterns Loosely Based on Incidence of Lung Cancer in U.S. (https://www.cdc.gov/cancer/lung/statistics/race.htm)
trueOR <- log(2)
baseline_Y <- log(0.005)
set.seed(3)
pop$Y <- rbinom(N, size = 1, prob = expit(baseline_Y + 
                                            trueOR*pop$A +
                                            (log(1.1))*pop$black + 
                                            (log(0.9))*pop$asian + 
                                            (log(1.1))*pop$hispanic +
                                            (log(1))*pop$otherrace + 
                                            (log(1.1))*pop$age_25_34 +
                                            (log(1.2))*pop$age_35_44 +
                                            (log(1.3))*pop$age_45_54 +
                                            (log(1.4))*pop$age_55_64 +
                                            (log(3))*pop$age_over64 +
                                            (log(3))*pop$male +
                                            (log(3))*pop$educ_ged +
                                            (log(1.1))*pop$educ_hs + 
                                            (log(1))*pop$educ_somecollege +
                                            (log(0.9))*pop$educ_associates +
                                            (log(0.8))*pop$educ_bachelors +
                                            (log(0.7))*pop$educ_advdegree))

## CALCULATING TRUE MEASURES OF ASSOCIATION AND DISPERSION
# Calculate True OR, CIR, and IDR for total population
# OR calculated using logistic regression model
trueOR_mod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='binomial')
trueOR <- as.numeric(exp(coef(trueOR_mod)["A"]))

# CIR calculated using log binomial model
trueCIR_lbmod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family= binomial(log))
trueCIR_lb <- as.numeric(exp(coef(trueCIR_lbmod)["A"]))

# CIR caluclated using poisson model
trueCIR_pmod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, family='poisson') # negative binomial model returned 1.894534
trueCIR_p <- as.numeric(exp(coef(trueCIR_pmod)["A"]))

# CIR calculated using negative binomial model
trueCIR_nbmod <- glm.nb(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, data=pop, control=glm.control(maxit=50))
trueCIR_nb <- as.numeric(exp(coef(trueCIR_nbmod)["A"]))

# IDR calculated with log binomial model
trueIDR_lbmod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(time), data=pop, family= binomial(log))
trueIDR_lb <- as.numeric(exp(coef(trueIDR_lbmod)["A"]))

# IDR calculated using poisson model with case occurrence offset
trueIDR_pmod <- glm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(time), data=pop, family='poisson')
trueIDR_p <- as.numeric(exp(coef(trueIDR_pmod)["A"]))

# IDR using modified poisson model
pop$id <- seq(1:nrow(pop))
trueIDR_mpmod <-geeglm(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree, offset = log(time), data=pop, id = id, family=poisson(link="log"), corstr="exchangeable")
trueIDR_mp <- as.numeric(exp(coef(trueIDR_mpmod)["A"]))

# IDR calculated using negative binomial model
trueIDR_nbmod <- glm.nb(Y ~ A + black + asian + hispanic + otherrace + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_over64 + male + educ_ged + educ_hs + educ_somecollege + educ_associates + educ_bachelors + educ_advdegree + offset(log(time)), data=pop, control=glm.control(maxit=100))
trueIDR_nb <- as.numeric(exp(coef(trueIDR_nbmod)["A"]))
# error with nb model: alternation limit reached/did not converge

# Looking at dispersion CIR
dispersiontest(trueCIR_pmod) 
pchisq(2 * (logLik(trueCIR_pmod) - logLik(trueCIR_nbmod)), df = 1, lower.tail = FALSE)

# Looking at dispersion IDR
dispersiontest(trueIDR_mod) 
pchisq(2 * (logLik(trueIDR_mod) - logLik(trueIDR_mod.nb)), df = 1, lower.tail = FALSE) #shows negative binomial model preferred when it converges?