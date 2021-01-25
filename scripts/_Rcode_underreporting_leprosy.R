########################################################################
#####   Correcting for Under-Reporting in Brazilian Leprosy data    ####
####                Paper: Oliveira et al. (2021)                   ####
####            Model specification and Implementation              ####
####    Based on supplementary materials of Stoner et al. (2019)    ####
########################################################################

## Directory ----
#getwd()
#setwd("??") #specify your working directory here

## Packages ----
#options(repos = c(CRAN = "http://cran.rstudio.com"))
#setRepositories(addURLs=c(CRANxtras = "https://cran-r.c3sl.ufpr.br/"))
if(!require(dplyr)){ install.packages("dplyr"); require(dplyr)}   
if(!require(stringr)){ install.packages("stringr"); require(stringr)} 
if(!require(devtools)){ install.packages("devtools"); require(devtools)} 
if(!require(coda)){ install.packages("coda"); require(coda)}   
if(!require(spdep)){ install.packages("spdep"); require(spdep)}   # For computing the neighbourhood and adjancency objects.
if(!require(nimble)){ install.packages("nimble"); require(nimble)} # For MCMC computation using NIMBLE.  


##
#---- Load in the data
load('_hansen_brasil_micro.RData') #variables for the model


##
#---- Part One: Model code for Nimble --

hansen_code=nimbleCode({ 
  for(i in 1:n_regions){
    epsilon[i] <- ilogit(alpha[1]+alpha[2]*W1[i]+gamma[i])
    theta[i] <- exp(log(POP[i])+beta[1]+beta[2]*X1[i]+beta[3]*X2[i]+
                      beta[4]*X3[i]+beta[5]*X4[i]+beta[6]*X5[i]+phi[i]+delta[i])
    mu[i] <- theta[i]*epsilon[i]
    Y[i] ~ dpois(mu[i])
    gamma[i]~dnorm(0,tau=sd_gamma)
    delta[i]~dnorm(0,tau=sd_delta)
  }
  phi[1:n_regions] ~ dcar_normal(adj=adj[1:l_adj], num=n_adj[1:n_regions], tau=nu, zero_mean=1)
  beta[1] ~ dnorm(-8,sd=1)
  for(j in 2:n_covariates_theta){
    beta[j] ~ dnorm(0, sd=sd_beta)
  }
  alpha[1] ~ dnorm(2.5, sd=0.3) #if precision, use tau=1/0.3^2
  for(k in 2:n_covariates_epsilon){
    alpha[k] ~ dnorm(0, sd=sd_alpha)
  }
  sd_gamma ~ dgamma(1,1)
  sd_delta ~ dgamma(1,1)
  nu ~ dgamma(1,1)
})

#
## Preditive analisys of prior elicited for \alpha_0:
if(!require(boot)){install.packages("boot"); require(boot)}   
alpha0_prior <- rnorm(1000, mean=2.5, sd=0.3)
hist(inv.logit(alpha0_prior))
summary(inv.logit(alpha0_prior))
quantile(inv.logit(alpha0_prior), prob=c(.05,.95))

#
##---- Part Two: Set up necessary data for the Nimble model--

# Set up necessary data for the ICAR prior.
adjacency=unlist(neighbourhood)
n_adj=card(neighbourhood)
l_adj=length(adjacency)
weights=rep(1,l_adj)

n_regions=length(n_adj) # Number of regions = 557.

# Centering covariates of reporting model according to region of reference
names(hansen_data_micro)
which(hansen_data_micro$NOME_MICRO=="RIOPRETODAEVA13") #microrregiao Rio Preto da Eva COD 13008 posição 15
pgrau2_centered=scale(hansen_data_micro$PROP_GRAU2,center=hansen_data_micro$PROP_GRAU2[15], scale=F)
pgrau2_centered = as.numeric(pgrau2_centered*100) 
print(hansen_data_micro$PROP_GRAU2*100)
summary(hansen_data_micro$PROP_GRAU2*100)

hansen_constants=list(n_regions = n_regions,
                      POP = hansen_data_micro$POP,
                      W1 = as.numeric(pgrau2_centered), 
                      X1 = as.numeric(scale(hansen_data_micro$CONT_EXAM*100, scale=F)),
                      X2 = as.numeric(scale(hansen_data_micro$COB_BF_POPALVO, scale=F)),
                      X3 = as.numeric(scale(hansen_data_micro$COB_PSF, scale=F)),
                      X4 = as.numeric(scale(hansen_data_micro$MED_MORA_DOM, scale=F)),
                      X5 = as.numeric(scale(hansen_data_micro$URBANISATION, scale=F)),
                      n_covariates_epsilon = 1 + 1, #"1" representa alpha0
                      n_covariates_theta = 1 + 5,  #"1" representa beta0
                      sd_beta = 10,
                      sd_alpha = 10,
                      adj = adjacency, 
                      l_adj=length(adjacency),
                      n_adj=n_adj,
                      w = weights
)


hansen_data=list(Y = hansen_data_micro$HANSEN)


#
##Regular Poisson  e Bonimial Negativo model just for checking
mPoisson <- glm(formula = hansen_data$Y ~ hansen_constants$X1 +  
                  hansen_constants$X2 + hansen_constants$X3 + 
                  hansen_constants$X4 + hansen_constants$X5 + 
                  offset(log(hansen_constants$POP)), 
                  family = poisson(link = "log"))
summary(mPoisson)
round(exp(mPoisson$coefficients),4)
round(exp(confint(mPoisson)),4)


if(!require(MASS)){ install.packages("MASS"); require(MASS)}   
mNegBinomial <- glm.nb(formula = hansen_data$Y ~ hansen_constants$X1 +  
                  hansen_constants$X2 + hansen_constants$X3 + 
                  hansen_constants$X4 + hansen_constants$X5 + 
                  offset(log(hansen_constants$POP)))
summary(mNegBinomial)
round(exp(mNegBinomial$coefficients),4)
round(exp(confint(mNegBinomial)),4)




#
##
# Set initial values.
#Chain 1
hansen_inits1=list(sd_delta=0.25, sd_gamma=0.25, nu=1,
                   alpha=c(2.9,rep(0.1, hansen_constants$n_covariates_epsilon-1)),
                   beta=c(-5,rep(0.1, hansen_constants$n_covariates_theta-1)),
                   phi = rnorm(n_regions,0,0.25),
                   delta=rnorm(n_regions,0,1),
                   gamma=rnorm(n_regions,0,.25)
)

#Chain 2
hansen_inits2=list(sd_delta=0.5, sd_gamma=0.75, nu=0.5,
                   alpha=c(2.1,rep(-0.1, hansen_constants$n_covariates_epsilon-1)),
                   beta=c(-9,rep(-0.1, hansen_constants$n_covariates_theta-1)),
                   phi = rnorm(n_regions,0,0.5),
                   delta=rnorm(n_regions,0,.5),
                   gamma=rnorm(n_regions,0,.5)
)

hansen_inits=list(chain1=hansen_inits1, chain2=hansen_inits2)


# Build the model.
hansen_model <- nimbleModel(hansen_code, hansen_constants, hansen_data, hansen_inits)
hansen_compiled_model <- compileNimble(hansen_model, resetFunctions = TRUE)


# Set up samplers.
hansen_mcmc_conf <- configureMCMC(hansen_model,monitors=c('beta','alpha','sd_delta',
                                                          'delta','gamma', 
                                                          'nu','sd_gamma',
                                                          'epsilon','theta','phi'), 
                                  useConjugacy = TRUE)
hansen_mcmc_conf$removeSamplers(c('beta[1]','alpha[1]','nu','sd_delta','sd_gamma'))
hansen_mcmc_conf$addSampler(target=c('beta[1]','alpha[1]','sd_gamma'),type='AF_slice')
hansen_mcmc_conf$addSampler(target=c('sd_delta','nu'),type='AF_slice')

hansen_mcmc <- buildMCMC(hansen_mcmc_conf)
hansen_compiled_mcmc <- compileNimble(hansen_mcmc, project = hansen_model, resetFunctions = TRUE)


#
##---- Part Three: Execution using nimble
## Run the model (a few hours ~ 4 hours).
##CAUTION!!! If you have the "output_leprosy" file, do not execute lines 183-188  
##Posterior samples are already saved in object "hansen_samples"
##Posterior summaries are already saved in file "output_leprosy"
#time <- proc.time()  #Iniciando tempo
#hansen_samples=runMCMC(hansen_compiled_mcmc, inits=hansen_inits,
#                       nchains = 2, nburnin=200000, niter = 600000, samplesAsCodaMCMC = TRUE, 
#                       thin=200, summary = FALSE, WAIC = FALSE, 
#                       setSeed=c(1,2)) 
#time = proc.time()-time; print(time)    #Finalizando tempo

# Check chains for convergence.
win.graph()
plot(hansen_samples[,c('alpha[1]','alpha[2]')])
savePlot(filename = "convergence_alphas" ,type="pdf")
dev.off()
win.graph()
plot(hansen_samples[,c('beta[1]','beta[2]','beta[3]')])
savePlot(filename = "convergence_betas1" ,type="pdf")
dev.off()
win.graph()
plot(hansen_samples[,c('beta[4]','beta[5]','beta[6]')])
savePlot(filename = "convergence_betas2" ,type="pdf")
dev.off()
win.graph()
plot(hansen_samples[,c('sd_delta','sd_gamma','nu')])
savePlot(filename = "convergence_precisions" ,type="pdf")
dev.off()


gelman.diag(hansen_samples[,c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]',
                              'alpha[1]','alpha[2]','sd_delta','sd_gamma','nu')])



#
##---- Part Four: Model Checking
# Combine MCMC chains.
hansen_mcmc=do.call('rbind',hansen_samples)

# Compute posterior quantities.
theta.names <- NULL; for(j in 1:n_regions){theta.names <-c(theta.names,paste("theta[",j,"]",sep=""))}
epsilon.names <- NULL; for(j in 1:n_regions){epsilon.names <-c(epsilon.names,paste("epsilon[",j,"]",sep=""))}
phi.names <- NULL; for(j in 1:n_regions){phi.names <-c(phi.names,paste("phi[",j,"]",sep=""))}
delta.names <- NULL; for(j in 1:n_regions){delta.names <-c(delta.names,paste("delta[",j,"]",sep=""))}
gamma.names <- NULL; for(j in 1:n_regions){gamma.names <-c(gamma.names,paste("gamma[",j,"]",sep=""))}
beta.names <- NULL; for(j in 1:hansen_constants$n_covariates_theta){beta.names <-c(beta.names,paste("beta[",j,"]",sep=""))}
alpha.names <- NULL; for(j in 1:hansen_constants$n_covariates_epsilon){alpha.names <-c(alpha.names,paste("alpha[",j,"]",sep=""))}

#selecting posterior sample of specific parameters:
posterior_beta=hansen_mcmc[,beta.names]
posterior_alpha=hansen_mcmc[,alpha.names]
#posterior_mu=hansen_mcmc[,mu.names]
posterior_theta=hansen_mcmc[,theta.names]
posterior_epsilon=hansen_mcmc[,epsilon.names]
#posterior_phi=hansen_mcmc[,phi.names]
#posterior_delta=hansen_mcmc[,delta.names]
#posterior_gamma=hansen_mcmc[,gamma.names]
#posterior_tau=hansen_mcmc[,"tau"]
#posterior_sd_delta=hansen_mcmc[,"sd_delta"]
#posterior_sd_gamma=hansen_mcmc[,"sd_gamma"]


# Simulate y.
posterior_y=t(apply(posterior_epsilon*posterior_theta,1,function(x)rpois(n_regions,x)))
posterior_lmse=apply(posterior_y,1,function(x) log(mean((x-hansen_data$Y)^2)))
mean(posterior_lmse)
mean(exp(posterior_lmse))

# Simulate t values.#predictive analysis of the total number of cases T_i
set.seed(12345)
posterior_t=t(apply(posterior_theta*(1-posterior_epsilon),1,function(x)rpois(n_regions,x)+hansen_data_micro$HANSEN))
posterior_total_t=(apply(posterior_t, 1, sum))
total_cases_corrected = round(t(apply(posterior_t, 2, mean)))
total_obs <- sum(hansen_data_micro$HANSEN)
posterior_total_extra_t=posterior_total_t-total_obs
hist(posterior_total_extra_t)
quantile(posterior_total_extra_t, seq(0,1,0.05))
mean(posterior_total_extra_t)
median(posterior_total_extra_t)
sd(posterior_total_extra_t)
HPDinterval(as.mcmc(posterior_total_extra_t),prob = .95)



#
##posterior odds ratio and irr
posterior_irr <- apply(posterior_beta,2,exp)
posterior_mean_irr <- apply(posterior_irr,2,mean)
HPDinterval(as.mcmc((posterior_irr)))
posterior_odds <- apply(posterior_alpha,2,exp)
posterior_mean_odds <- apply(posterior_odds,2,mean)
HPDinterval(as.mcmc((posterior_odds)))


#
##---- Part Five: Posterior results
# Create the function to find the posterior mode (if desired):
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
apply(posterior_alpha,2,getmode)
apply(posterior_beta,2,getmode)


#Following are the analysis based on the posterior mean:
post.mean <- round(apply(hansen_mcmc,2,mean),4) #media a posteriori de todos os parametros
post.mean[beta.names]
post.mean[alpha.names]


#saving important quantities:
posterior_total_rate=100000*as.numeric(post.mean[theta.names])/hansen_constants$POP
COD_UF <- str_sub(string = hansen_data_micro$ID_MICRO, start = 1, end = 2)
SIGLA_UF <- numeric(length(COD_UF))
SIGLA_UF[which(COD_UF=="11")] <- "RO"; SIGLA_UF[which(COD_UF=="12")] <- "AC"; SIGLA_UF[which(COD_UF=="13")] <- "AM"
SIGLA_UF[which(COD_UF=="14")] <- "RR"; SIGLA_UF[which(COD_UF=="15")] <- "PA"; SIGLA_UF[which(COD_UF=="16")] <- "AP"
SIGLA_UF[which(COD_UF=="17")] <- "TO"; SIGLA_UF[which(COD_UF=="21")] <- "MA"; SIGLA_UF[which(COD_UF=="22")] <- "PI"
SIGLA_UF[which(COD_UF=="23")] <- "CE"; SIGLA_UF[which(COD_UF=="24")] <- "RN"; SIGLA_UF[which(COD_UF=="25")] <- "PB"
SIGLA_UF[which(COD_UF=="26")] <- "PE"; SIGLA_UF[which(COD_UF=="27")] <- "AL"; SIGLA_UF[which(COD_UF=="28")] <- "SE"
SIGLA_UF[which(COD_UF=="29")] <- "BA"; SIGLA_UF[which(COD_UF=="31")] <- "MG"; SIGLA_UF[which(COD_UF=="32")] <- "ES"
SIGLA_UF[which(COD_UF=="33")] <- "RJ"; SIGLA_UF[which(COD_UF=="35")] <- "SP"; SIGLA_UF[which(COD_UF=="41")] <- "PR"
SIGLA_UF[which(COD_UF=="42")] <- "SC"; SIGLA_UF[which(COD_UF=="43")] <- "RS"; SIGLA_UF[which(COD_UF=="50")] <- "MS"
SIGLA_UF[which(COD_UF=="51")] <- "MT"; SIGLA_UF[which(COD_UF=="52")] <- "GO"; SIGLA_UF[which(COD_UF=="53")] <- "DF"


#saving data used in the model and posterior summaries:
output <- data.frame(ID_MICRO = hansen_data_micro$ID_MICRO,
                     NOME_MICRO = hansen_data_micro$NOME_MICRO, 
                     COD_UF = COD_UF,
                     SIGLA_UF = SIGLA_UF,
                     X1 = round(as.numeric(scale(hansen_data_micro$CONT_EXAM*100, scale=F)),5),
                     X2 = round(as.numeric(scale(hansen_data_micro$COB_BF_POPALVO, scale=F)),5),
                     X3 = round(as.numeric(scale(hansen_data_micro$COB_PSF, scale=F)),5),
                     X4 = round(as.numeric(scale(hansen_data_micro$MED_MORA_DOM, scale=F)),5),
                     X5 = round(as.numeric(scale(hansen_data_micro$URBANISATION, scale=F)),5),
                     W1 = round(as.numeric(pgrau2_centered),5),
                     POP = hansen_data_micro$POP,
                     total_cases_observed_y = as.numeric(hansen_data_micro$HANSEN),
                     incidence_observed = round(as.numeric(hansen_data_micro$INCIDENCE),5),
                     total_cases_corrected_t = as.numeric(total_cases_corrected),
                     incidence_corrected = round(as.numeric(posterior_total_rate),5),
                     epsilon_posterior_mean = round(as.numeric(post.mean[epsilon.names]),5),
                     theta_posterior_mean = round(as.numeric(post.mean[theta.names]),5),
                     phi_posterior_mean = round(as.numeric(post.mean[phi.names]),5),
                     delta_posterior_mean = round(as.numeric(post.mean[delta.names]),5),
                     gamma_posterior_mean = round(as.numeric(post.mean[gamma.names]),5),
                     PROP_GRAU2 = round(as.numeric(hansen_data_micro$PROP_GRAU2*100),5),
                     CONT_EXAM = round(as.numeric(hansen_data_micro$CONT_EXAM*100),5),
                     COB_BOLSAF = round(as.numeric(hansen_data_micro$COB_BF_POPALVO),5),
                     COB_PSF = round(as.numeric(hansen_data_micro$COB_PSF),5),
                     MED_MORA_DOM = round(as.numeric(hansen_data_micro$MED_MORA_DOM),5),
                     URBANISATION = round(as.numeric(hansen_data_micro$URBANISATION),5)
)
save(output, file="output_leprosy.RData")
write.csv2(output,  file="output_leprosy.csv", row.names=FALSE)


##
###obtaining some descriptive measures
head(output)
names(output)
statistics_by_state <- output %>%
  select(c(COD_UF, total_cases_observed_y, total_cases_corrected_t)) %>%
  group_by(COD_UF) %>%
  summarise_all(funs(sum))

statistics_by_brazil <- data.frame(COD_UF = "BR", 
                                   total_cases_observed_y = as.numeric(sum(statistics_by_state$total_cases_observed_y)), 
                                   total_cases_corrected_t = as.numeric(sum(statistics_by_state$total_cases_corrected_t)))

statistics_by_state <- rbind(statistics_by_state, statistics_by_brazil)

statistics_by_state <- statistics_by_state %>% 
  mutate(unobserved_cases_z = total_cases_corrected_t - total_cases_observed_y,
         leprosy_detection = total_cases_observed_y/total_cases_corrected_t,
         SIGLA_UF = c("RO", "AC", "AM", "RR", "PA", "AP", "TO", "MA", "PI", "CE", "RN",
                      "PB", "PE", "AL", "SE", "BA", "MG", "ES", "RJ", "SP", "PR", "SC",
                      "RS", "MS", "MT", "GO", "DF","BR"))
statistics_by_state
#

# rm(brasil_map, hansen_mcmc)
# rm(list = ls(pattern = "^posterior_"))
# rm(list = ls(pattern = ".names"))
# rm(brasil_micro, hansen_data_micro, neighbourhood)


