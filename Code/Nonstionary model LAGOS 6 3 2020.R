
# Code for running nonstationary model on LAGOS data including data processing,
# estimation of priors, running the nonstationary model, and checking output diagnostics.
# This code is for total lake phosphorus in the Northeast region. 
# Model parameters and priors for other response variables and regions can be
# found in supplementary information.

############## Data Processing ###################################

#only use use lakes with most recent data
D<-read.csv("LagosNonSecchi_16Jan2019.csv")

d_0<-na.omit(subset(D,D$state_name %in% c('New York')))[,-1] #removes *ALL* NAs


#take TP median per lake over time
library(reshape2)
melted<-melt(d_0,id=c(1:2,9:25),na.rm=TRUE)
d<-dcast(subset(melted,melted$variable=='tp'),
         nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
           hu8_MAP + hu8_MAT + hu8_runoff +
           ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
           maxdepth + lakeconnection + la_wa_ratio ~ variable,
         median,drop=TRUE)
summary(d)
d$log_tp<-log(d$tp)
rm(d_0,D,melted)

############ determine priors for stan model ########################
library(spaMM)
M_iso<-corrHLfit(log_tp ~ scale(hu8_totalndepo) + scale(ws_forest) + scale(maxdepth)+scale(la_wa_ratio) + 
                   Matern(1|nhd_long+nhd_lat),
                 ranFix=list(nu=0.5),
                 HLmethod="HL(1,1)",data=na.omit(d))
M_iso #1/rho here is phi in geoR likfit, so eff rng is 3/rho
AIC(M_iso) # priors: intercept(mean) = 2.85  / tot var = .2053+.4174= .62  / prior for range = 1.75


############ Run nonstationary models ########################

library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)

################# model 1 ################
#*stan analysis*#
y<-as.numeric(d$log_tp)
X<-model.matrix(~ -1 + scale(hu8_totalndepo) + scale(ws_forest)+ scale(maxdepth)+scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ -1 + scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #strech
X_T<-model.matrix(~ -1 + scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #angle of rotation
V_L<-ncol(X_L)
V_T<-ncol(X_T)
#**pass to stan**#
data<-list(y=y, h_xy=h_xy,  X=X,X_L=X_L,X_T=X_T,
           V_X=V_X,V_T=V_T,V_L=V_L, S=S,N=N) 
pars<-c("Int_mu","beta_mu","beta_L","beta_T",   
        "alpha_sq","nug_sq","rho","lambda_sig","theta_sig",
        "log_lik","bayes_R2")
str(data) #<-check 

#**run!**#
fit_tp_NS <- stan(
  file='NS_gaussian_model_H_2c_TP_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.93, max_treedepth=15),
  warmup = 300, iter = 1800, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tp_NS,"tp_ny_NS_samples_07102019_1.rds")
fit_tp_NS
check_hmc_diagnostics(fit_tp_NS)

################# model 2 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tp)
X<-model.matrix(~ -1 + scale(hu8_totalndepo) + scale(ws_forest)+ scale(maxdepth)+scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ -1 + scale(hu8_baseflow) + scale(ws_urban),data=d) #strech
X_T<-model.matrix(~ -1 + scale(hu8_baseflow) + scale(ws_urban),data=d) #angle of rotation
V_L<-ncol(X_L)
V_T<-ncol(X_T)
#**pass to stan**#
data<-list(y=y, h_xy=h_xy,  X=X,X_L=X_L,X_T=X_T,
           V_X=V_X,V_T=V_T,V_L=V_L, S=S,N=N) 
pars<-c("Int_mu","beta_mu","beta_L","beta_T",   
        "alpha_sq","nug_sq","rho","lambda_sig","theta_sig",
        "log_lik","bayes_R2")
str(data) #<-check 

#**run!**#
fit_tp_NS <- stan(
  file='NS_gaussian_model_H_2c_TP_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tp_NS,"tp_ny_NS_samples_07102019_2.rds")
fit_tp_NS
check_hmc_diagnostics(fit_tp_NS)

################# model 3 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tp)
X<-model.matrix(~ -1 + scale(hu8_totalndepo) + scale(ws_forest)+ scale(maxdepth)+scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ -1+scale(hu8_baseflow) + scale(ws_forest),data=d) #strech
X_T<-model.matrix(~ -1+scale(hu8_baseflow) + scale(ws_forest),data=d) #angle of rotation
V_L<-ncol(X_L)
V_T<-ncol(X_T)
#**pass to stan**#
data<-list(y=y, h_xy=h_xy,  X=X,X_L=X_L,X_T=X_T,
           V_X=V_X,V_T=V_T,V_L=V_L, S=S,N=N) 
pars<-c("Int_mu","beta_mu","beta_L","beta_T",   
        "alpha_sq","nug_sq","rho","lambda_sig","theta_sig",
        "log_lik","bayes_R2")
str(data) #<-check 

#**run!**#
fit_tp_NS <- stan(
  file='NNS_gaussian_model_H_2c_TP_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tp_NS,"tp_ny_NS_samples_01312020_3.rds")
fit_tp_NS
check_hmc_diagnostics(fit_tp_NS)


#################checking diagnostics ##########
t<-readRDS("tp_ny_NS_samples_07102019_1.rds")

launch_shinystan(t,rstudio = TRUE)

# post convergence analysis #
b_L<-extract(t,"beta_L")[[1]]
library(HDInterval)
hdi(b_L[,1],0.9)
hdi(b_L[,2],0.9)
hdi(b_L[,3],0.9)
hdi(b_L[,4],0.9)

b_T<-extract(t,"beta_T")[[1]]
library(HDInterval)
hdi(b_T[,1],0.9)
hdi(b_T[,2],0.9)
hdi(b_T[,3],0.9)
hdi(b_T[,4],0.9)

b_mu<-extract(t,"beta_mu")[[1]]
hdi(b_mu[,1],0.9)
hdi(b_mu[,2],0.9)
hdi(b_mu[,3],0.9)
hdi(b_mu[,4],0.9)
hdi(b_mu[,5],0.9)
hdi(b_mu[,6],0.9)
hdi(b_mu[,7],0.9)
hdi(b_mu[,8],0.9)

I_mu<-extract(t,"Int_mu")[[1]]
hdi(I_mu,0.9)

a<-extract(t,"alpha_sq")[[1]]
hdi(a,0.9)

a<-extract(t,"nug_sq")[[1]]
hdi(a,0.9)

a<-extract(t,"rho")[[1]]
hdi(a,0.9)

a<-extract(t,"bayes_R2")[[1]]
hdi(a,0.9)



