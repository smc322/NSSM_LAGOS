setwd("C://Users//charlotte//")

################################################################################ start with TP NY ###############################################
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


library(spaMM)
M_iso<-corrHLfit(log_tp ~ scale(hu8_totalndepo) + scale(ws_forest) + scale(maxdepth)+scale(la_wa_ratio) + 
                   Matern(1|nhd_long+nhd_lat),
                 ranFix=list(nu=0.5),
                 HLmethod="HL(1,1)",data=na.omit(d))
M_iso #1/rho here is phi in geoR likfit, so eff rng is 3/rho
AIC(M_iso) # priors: intercept(mean) = 2.85  / tot var = .2053+.4174= .62  / prior for range = 1.75


library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)

################# model 1 #################################################################
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
  file='NS_gaussian_model_H_2c_11052019_TP_NY.stan',
  #file='NS_gaussian_model_H_2c_TP_NY.stan',
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
  file='NS_gaussian_model_H_2c_11052019_TP_NY',
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
  file='NS_gaussian_model_H_2c_11052019_TP_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tp_NS,"tp_ny_NS_samples_01312020_3.rds")
fit_tp_NS
check_hmc_diagnostics(fit_tp_NS)

################################################################################ TN NY ###############################################
#only use use lakes with most recent data
D<-read.csv("LagosNonSecchi_16Jan2019.csv")

d_0<-na.omit(subset(D,D$state_name %in% c('New York')))[,-1] #removes *ALL* NAs


#take TN median per lake over time
library(reshape2)
melted<-melt(d_0,id=c(1:5,9:25),na.rm=TRUE)
d<-dcast(subset(melted,melted$variable=='tn_all'),
         nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
           hu8_MAP + hu8_MAT + hu8_runoff +
           ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
           maxdepth + lakeconnection + la_wa_ratio ~ variable,
         median,drop=TRUE)
summary(d)
d$log_tn<-log(d$tn_all)
rm(d_0,D,melted)

library(spaMM)
M_iso<-corrHLfit(log_tn ~ scale(hu8_totalndepo) + scale(ws_forest) + scale(maxdepth)+scale(la_wa_ratio) + 
                   Matern(1|nhd_long+nhd_lat),
                 ranFix=list(nu=0.5),
                 HLmethod="HL(1,1)",data=na.omit(d))
M_iso #1/rho here is phi in geoR likfit, so eff rng is 3/rho
AIC(M_iso) # priors: intercept(mean) = 6.07  / tot var = .1009+.09359 = 0.19449  / prior for range = 0.4163904


library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)

################# model 1  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tn)
X<-model.matrix(~ -1 + scale(hu8_totalndepo) + scale(ws_forest)+ scale(maxdepth)+scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #angle of rotation
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
fit_tn_NS <- stan(
  file='NS_gaussian_model_H_2c_TN_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=10),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tn_NS,"tn_ny_NS_samples_07012019_1.rds")
fit_tn_NS
check_hmc_diagnostics(fit_tn_NS)

################# model 2 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tn)
X<-model.matrix(~ -1 + scale(hu8_totalndepo) + scale(ws_forest)+ scale(maxdepth)+scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #angle of rotation
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
fit_tn_NS <- stan(
  file='NS_gaussian_model_H_2c_TN_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tn_NS,"tn_ny_NS_samples_07012019_2.rds")
fit_tn_NS
check_hmc_diagnostics(fit_tn_NS)

################# model 3 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tn)
X<-model.matrix(~ -1 + scale(hu8_totalndepo) + scale(ws_forest)+ scale(maxdepth)+scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #angle of rotation
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
fit_tn_NS <- stan(
  file='NS_gaussian_model_H_2c_TN_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tn_NS,"tn_ny_NS_samples_07012019_3.rds")
fit_tn_NS
check_hmc_diagnostics(fit_tn_NS)

################################################################################ TP MW ###############################################
#only use use lakes with most recent data
D<-read.csv("LagosNonSecchi_16Jan2019.csv")

d_0<-na.omit(subset(D,D$state_name %in% c('Iowa','Wisconsin','Illinois')))[,-1] #removes *ALL* NAs

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

####*** GAM exploratory analysis ***###
# interior knots
int_knots<-data.frame(fields::cover.design(cbind(d$nhd_long,d$nhd_lat),10,max.loop=50)$design)
names(int_knots)<-c("x","y")
library(mgcv)
sm_tp0<-gam(hu8_totalndepo~ s(nhd_long,nhd_lat,k=110,bs='tp'),
            data=na.omit(d),method="REML")
summary(sm_tp0)

library(spaMM)
M_iso<-corrHLfit(log_tp ~ scale(hu8_baseflow) + scale(ws_forest) + scale(maxdepth)+lakeconnection + 
                   Matern(1|nhd_long+nhd_lat),
                 ranFix=list(nu=0.5),
                 HLmethod="HL(1,1)",data=na.omit(d))
M_iso #1/rho here is phi in geoR likfit, so eff rng is 3/rho
AIC(M_iso) # priors: intercept(mean) = 4.15  / tot var = .3078+.2856 = .59  / prior for range = 2.34

install.packages("shinystan")
library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)

################# model 1 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tp)
X<-model.matrix(~ -1 + scale(hu8_baseflow) + scale(ws_forest) + scale(maxdepth)+lakeconnection , data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow,center=T) + scale(hu8_totalndepo,center=T),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow,center=T) + scale(hu8_totalndepo,center=T),data=d) #angle of rotation
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
  file='NS_gaussian_model_H_2c_TP_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=10),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tp_NS,"tp_mw_NS_samples_07012019_1.rds")
fit_tp_NS
check_hmc_diagnostics(fit_tp_NS)

################# model 2 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tp)
X<-model.matrix(~ -1 + scale(hu8_baseflow) + scale(ws_forest) + scale(maxdepth)+lakeconnection , data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #angle of rotation
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
  file='NS_gaussian_model_H_2c_TP_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tp_NS,"tp_mw_NS_samples_07082019_2.rds")
fit_tp_NS
check_hmc_diagnostics(fit_tp_NS)

################# model 3 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tp)
X<-model.matrix(~ -1 + scale(hu8_baseflow) + scale(ws_forest) + scale(maxdepth)+lakeconnection , data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #angle of rotation
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
  file='NS_gaussian_model_H_2c_TP_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tp_NS,"tp_mw_NS_samples_07042019_3.rds")
fit_tp_NS
check_hmc_diagnostics(fit_tp_NS)

################################################################################ TN MW ###############################################
#only use use lakes with most recent data
D<-read.csv("LagosNonSecchi_16Jan2019.csv")

d_0<-na.omit(subset(D,D$state_name %in% c('Iowa','Wisconsin','Illinois')))[,-1] #removes *ALL* NAs


#take TP median per lake over time
library(reshape2)
melted<-melt(d_0,id=c(1:5,9:25),na.rm=TRUE)
d<-dcast(subset(melted,melted$variable=='tn_all'),
         nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
           hu8_MAP + hu8_MAT + hu8_runoff +
           ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
           maxdepth + lakeconnection + la_wa_ratio ~ variable,
         median,drop=TRUE)
summary(d)
d$log_tn<-log(d$tn_all)
rm(d_0,D,melted)

library(spaMM)
M_iso<-corrHLfit(log_tn ~ scale(hu8_baseflow) + scale(hu8_totalndepo)+ scale(ws_forest) + scale(maxdepth)+lakeconnection + 
                   Matern(1|nhd_long+nhd_lat),
                 ranFix=list(nu=0.5),
                 HLmethod="HL(1,1)",data=na.omit(d))
M_iso #1/rho here is phi in geoR likfit, so eff rng is 3/rho
AIC(M_iso) # priors: intercept(mean) = 7.21/ tot var = .1977+.1734 = .59  / prior for range = 1.4206


library(rstan)
library(shinystan)
install.packages('tibble')
library(tibble)
install.packages('dplyr')
library(dplyr)
rstan_options(auto_write = TRUE)

################# model 1  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tn)
X<-model.matrix(~ -1 + scale(hu8_totalndepo) + scale(ws_forest)+ scale(maxdepth)+scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #angle of rotation
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
fit_tn_NS <- stan(
  file='NS_gaussian_model_H_2c_TN_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=10),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tn_NS,"tn_mw_NS_samples_07012019_1.rds")
fit_tn_NS
check_hmc_diagnostics(fit_tn_NS)

################# model 2 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tn)
X<-model.matrix(~ -1 + scale(hu8_totalndepo) + scale(ws_forest)+ scale(maxdepth)+scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #angle of rotation
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
fit_tn_NS <- stan(
  file='NS_gaussian_model_H_2c_TN_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 600, iter = 2600, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tn_NS,"tn_mw_NS_samples_09032019_2.rds")
fit_tn_NS
check_hmc_diagnostics(fit_tn_NS)

################# model 3 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_tn)
X<-model.matrix(~ -1 + scale(hu8_totalndepo) + scale(ws_forest)+ scale(maxdepth)+scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #angle of rotation
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
fit_tn_NS <- stan(
  file='NS_gaussian_model_H_2c_TN_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_tn_NS,"tn_mw_NS_samples_07012019_3.rds")
fit_tn_NS
check_hmc_diagnostics(fit_tn_NS)




###################################################### Chla NY #################################################
D<-read.csv("LagosNonSecchi_16Jan2019.csv")
table(D$state_name)

d_0<-na.omit(subset(D,D$state_name %in% c('New York')))[,-1] #removes *ALL* NAs

#take no2no3 median per lake over time
library(reshape2)
melted<-melt(d_0,id=c(1:2,4:6,9:25),na.rm=TRUE)
d<-dcast(subset(melted,melted$variable=='chla'),
         nhd_long + nhd_lat+ hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
           hu8_MAP + hu8_MAT + hu8_runoff +
           ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
           maxdepth + lakeconnection + la_wa_ratio ~ variable,
         median,drop=TRUE)

meltedp<-melt(d_0,id=c(1:2,9:25),na.rm=TRUE)
dp<-dcast(subset(meltedp,meltedp$variable=='tp'),
          nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
            hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio ~ variable,
          median,drop=TRUE)

dpa<-dp[order(dp$nhd_long),]
dpa<-dpa[order(dpa$nhd_lat),]

meltedn<-melt(d_0,id=c(1:5,9:25),na.rm=TRUE)
dn<-dcast(subset(meltedn,meltedn$variable=='tn_all'),
          nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
            hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio ~ variable,
          median,drop=TRUE)

summary(d)
d$log_chla<-log(d$chla)
d$tp<-dp$tp
d$tn_all<-dn$tn_all
rm(d_0,D,melted)

M0<-lm(log_chla ~ scale(tp)+scale(tn_all)+scale(hu8_baseflow)  + scale(hu8_totalndepo) +
         scale(ws_urban) + scale(ws_forest)  + 
         scale(lake_area) + scale(maxdepth) + lakeconnection + scale(la_wa_ratio),data=d)
summary(M0) 
vif(M0)

library(spaMM)
M_iso<-corrHLfit(log_chla ~ scale(tn_all)+scale(hu8_baseflow) + scale(la_wa_ratio)+ 
                   Matern(1|nhd_long+nhd_lat),
                 ranFix=list(nu=0.5),
                 HLmethod="HL(1,1)",data=na.omit(d))
M_iso #1/rho here is phi in geoR likfit, so eff rng is 3/rho
AIC(M_iso) # priors: intercept(mean) = 1.525199/ tot var =  0.3954+0.402 = 0.7974  / prior for range = 4.857961 
mean(d$log_chla)

################# model 1  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_chla)
X<-model.matrix(~ -1 + scale(tn_all)+scale(hu8_baseflow) + scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #angle of rotation
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
fit_chla_NS <- stan(
  file='NS_gaussian_model_H_2c_chla_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 600, iter = 2600, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_chla_NS,"chla_ny_NS_samples_09052019_1.rds")
fit_chla_NS
check_hmc_diagnostics(fit_chla_NS)

################# model 2 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_chla)
X<-model.matrix(~ -1 + scale(tn_all)+scale(hu8_baseflow) + scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #angle of rotation
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
fit_chla_NS <- stan(
  file='NS_gaussian_model_H_2c_chla_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=10),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_chla_NS,"chla_ny_NS_samples_07012019_2.rds")
fit_chla_NS
check_hmc_diagnostics(fit_chla_NS)

################# model 3 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_chla)
X<-model.matrix(~ -1 + scale(tn_all)+scale(hu8_baseflow) + scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #angle of rotation
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
fit_chla_NS <- stan(
  file='NS_gaussian_model_H_2c_chla_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=10),
  warmup = 600, iter = 2600, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_chla_NS,"chla_ny_NS_samples_8292019_3.rds")
fit_chla_NS
check_hmc_diagnostics(fit_chla_NS)


###################################################### Chla Mw ################################################
D<-read.csv("LagosNonSecchi_16Jan2019.csv")
table(D$state_name)

d_0<-na.omit(subset(D,D$state_name %in% c('Iowa','Wisconsin','Illinois')))[,-1] #removes *ALL* NAs

#take no2no3 median per lake over time
library(reshape2)
melted<-melt(d_0,id=c(1:2,4:6,9:25),na.rm=TRUE)
d<-dcast(subset(melted,melted$variable=='chla'),
         nhd_long + nhd_lat+ hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
           hu8_MAP + hu8_MAT + hu8_runoff +
           ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
           maxdepth + lakeconnection + la_wa_ratio ~ variable,
         median,drop=TRUE)

meltedp<-melt(d_0,id=c(1:2,9:25),na.rm=TRUE)
dp<-dcast(subset(meltedp,meltedp$variable=='tp'),
          nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
            hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio ~ variable,
          median,drop=TRUE)

dpa<-dp[order(dp$nhd_long),]
dpa<-dpa[order(dpa$nhd_lat),]

meltedn<-melt(d_0,id=c(1:5,9:25),na.rm=TRUE)
dn<-dcast(subset(meltedn,meltedn$variable=='tn_all'),
          nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
            hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio ~ variable,
          median,drop=TRUE)

summary(d)
d$log_chla<-log(d$chla)
d$tp<-dp$tp
d$tn_all<-dn$tn_all
rm(d_0,D,melted)

M0<-lm(log_chla ~ scale(tp)+scale(tn_all)+scale(hu8_baseflow)  + scale(hu8_totalndepo) +
         scale(ws_urban) + scale(ws_forest)  + 
         scale(lake_area) + scale(maxdepth) + lakeconnection + scale(la_wa_ratio),data=d)
summary(M0) 
vif(M0)

library(spaMM)
M_iso<-corrHLfit(log_chla ~ scale(tp)+scale(tn_all)+scale(hu8_baseflow) +
                   scale(ws_forest)  + scale(maxdepth) + lakeconnection + 
                   Matern(1|nhd_long+nhd_lat),
                 ranFix=list(nu=0.5),
                 HLmethod="HL(1,1)",data=na.omit(d))
M_iso #1/rho here is phi in geoR likfit, so eff rng is 3/rho
AIC(M_iso) # priors: intercept(mean) = 2.9898 / tot var =  0.5506 +0.5025 = 1.05  / prior for range = 0.849999 
mean(d$log_chla)

################# model 1  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_chla)
X<-model.matrix(~ -1 + scale(tp)+scale(tn_all)+scale(hu8_baseflow) +
                  scale(ws_forest)  + scale(maxdepth) + lakeconnection, data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #angle of rotation
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
fit_chla_NS <- stan(
  file='NS_gaussian_model_H_2c_chla_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_chla_NS,"chla_mw_NS_samples_07012019_1.rds")
fit_chla_NS
check_hmc_diagnostics(fit_chla_NS)

################# model 2 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_chla)
X<-model.matrix(~ -1 + scale(tp)+scale(tn_all)+scale(hu8_baseflow) +
                  scale(ws_forest)  + scale(maxdepth) + lakeconnection, data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #angle of rotation
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
fit_chla_NS <- stan(
  file='NS_gaussian_model_H_2c_chla_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_chla_NS,"chla_mw_NS_samples_07012019_2.rds")
fit_chla_NS
check_hmc_diagnostics(fit_chla_NS)

################# model 3 #################################################################
#*stan analysis*#
y<-as.numeric(d$log_chla)
X<-model.matrix(~ -1 + scale(tp)+scale(tn_all)+scale(hu8_baseflow) +
                  scale(ws_forest)  + scale(maxdepth) + lakeconnection, data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #angle of rotation
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
fit_chla_NS <- stan(
  file='NS_gaussian_model_H_2c_chla_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=10),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_chla_NS,"chla_mw_NS_samples_07012019_3.rds")
fit_chla_NS
check_hmc_diagnostics(fit_chla_NS)

###################################################################### N:P NY #############################################

D<-read.csv("LagosNonSecchi_16Jan2019.csv")
table(D$state_name)
D$N.P<-(D$tn_all/14.0067)/(D$tp/30.973762)


d_0<-na.omit(subset(D,D$state_name %in% c('New York')))[,-1]

#take N:P median per lake over time
library(reshape2)
melted<-melt(d_0,id=c(1:2,9:25),na.rm=TRUE)
d<-dcast(subset(melted,melted$variable=='N.P'),
         nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
           hu8_MAP + hu8_MAT + hu8_runoff +
           ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
           maxdepth + lakeconnection + la_wa_ratio ~ variable,
         median,drop=TRUE)
summary(d)
d$log_np<-log(d$N.P)
rm(d_0,D,melted)

M0<-lm(log_np ~ scale(hu8_baseflow)  + scale(hu8_totalndepo) +
         scale(ws_urban) + scale(ws_forest)  + 
         scale(lake_area) + scale(maxdepth) + lakeconnection + scale(la_wa_ratio),data=d)
summary(M0) 
vif(M0)

library(spaMM)
M_iso<-corrHLfit(log_np ~ scale(hu8_baseflow) +scale(hu8_totalndepo) +
                   scale(ws_forest) +scale(la_wa_ratio) +
                   Matern(1|nhd_long+nhd_lat),
                 ranFix=list(nu=0.5),
                 HLmethod="HL(1,1)",data=na.omit(d))
M_iso #1/rho here is phi in geoR likfit, so eff rng is 3/rho
AIC(M_iso) # priors: intercept(mean) = 4.004759 / tot var =  0.03426 + 0.3516 = .386  / prior for range = 5.160585 
mean(d$log_np)

################# model 1  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_np)
X<-model.matrix(~ -1 + scale(hu8_baseflow) +scale(hu8_totalndepo) +
                  scale(ws_forest) +scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #angle of rotation
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
fit_np_NS <- stan(
  file='NS_gaussian_model_H_2c_NP_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_np_NS,"np_ny_NS_samples_06152019_1.rds")
fit_np_NS
check_hmc_diagnostics(fit_np_NS)

################# model 2  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_np)
X<-model.matrix(~ -1 + scale(hu8_baseflow) +scale(hu8_totalndepo) +
                  scale(ws_forest) +scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #angle of rotation
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
fit_np_NS <- stan(
  file='NS_gaussian_model_H_2c_NP_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_np_NS,"np_ny_NS_samples_07222019_2.rds")
fit_np_NS
check_hmc_diagnostics(fit_np_NS)


################# model 3  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_np)
X<-model.matrix(~ -1 + scale(hu8_baseflow) +scale(hu8_totalndepo) +
                  scale(ws_forest) +scale(la_wa_ratio), data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #angle of rotation
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
  file='NS_gaussian_model_H_2c_NP_NY.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_np_NS,"np_ny_NS_samples_07232019_3.rds")
fit_np_NS
check_hmc_diagnostics(fit_np_NS)


################################################################ N:P MW ###################################

D<-read.csv("LagosNonSecchi_16Jan2019.csv")
table(D$state_name)
D$N.P<-(D$tn_all/14.0067)/(D$tp/30.973762)

d_0<-na.omit(subset(D,D$state_name %in% c('Iowa','Wisconsin','Illinois')))[,-1] #removes *ALL* NAs


#take N:P median per lake over time
library(reshape2)
melted<-melt(d_0,id=c(1:2,9:25),na.rm=TRUE)
d<-dcast(subset(melted,melted$variable=='N.P'),
         nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + 
           hu8_MAP + hu8_MAT + hu8_runoff +
           ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
           maxdepth + lakeconnection + la_wa_ratio ~ variable,
         median,drop=TRUE)
summary(d)
d$log_np<-log(d$N.P)
rm(d_0,D,melted)

M0<-lm(log_np ~ scale(hu8_baseflow)  + scale(hu8_totalndepo) +
         scale(ws_urban) + scale(ws_forest)  + 
         scale(lake_area) + scale(maxdepth) + lakeconnection + scale(la_wa_ratio),data=d)
summary(M0) 
vif(M0)

library(spaMM)
M_iso<-corrHLfit(log_np ~ scale(hu8_baseflow) +scale(maxdepth) +
                   lakeconnection +
                   Matern(1|nhd_long+nhd_lat),
                 ranFix=list(nu=0.5),
                 HLmethod="HL(1,1)",data=na.omit(d))
M_iso #1/rho here is phi in geoR likfit, so eff rng is 3/rho
AIC(M_iso) # priors: intercept(mean) = 3.851089 / tot var =   0.2461 +  0.3755 = .6216  / prior for range = 2.631478 
mean(d$log_np)

################# model 1  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_np)
X<-model.matrix(~ -1 + scale(hu8_baseflow)  +
                  scale(maxdepth) + lakeconnection, data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(hu8_totalndepo),data=d) #angle of rotation
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
  file='NS_gaussian_model_H_2c_NP_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=15),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_np_NS,"np_mw_NS_samples_06152019_1.rds")
fit_np_NS
check_hmc_diagnostics(fit_np_NS)

################# model 2  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_np)
X<-model.matrix(~ -1 + scale(hu8_baseflow)  +
                  scale(maxdepth) + lakeconnection, data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_urban),data=d) #angle of rotation
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
  file='NS_gaussian_model_H_2c_NP_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=10),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_np_NS,"np_mw_NS_samples_06152019_2.rds")
fit_np_NS
check_hmc_diagnostics(fit_np_NS)

################# model 3  #################################################################
#*stan analysis*#
y<-as.numeric(d$log_np)
X<-model.matrix(~ -1 + scale(hu8_baseflow)  +
                  scale(maxdepth) + lakeconnection, data=d)
V_X<-ncol(X)
N<-length(y)
h_xy<-as.matrix(d[,1:2]) 
S<-nrow(h_xy)
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #angle of rotation
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
fit_np_NS <- stan(
  file='NS_gaussian_model_H_2c_NP_MW.stan',
  model_name = 'Fully Hierarchical Gaussian N-S model 2c',
  data = data,pars=pars, seed=11191985,
  control=list(adapt_delta=0.9, max_treedepth=10),
  warmup = 300, iter = 1300, chains = 3,cores=3, thin=1, init_r=0.5)    
saveRDS(fit_np_NS,"np_mw_NS_samples_06152019_3.rds")
fit_np_NS
check_hmc_diagnostics(fit_np_NS)






##################################################checking diagnostics ##########
t<-readRDS("chla_mw_NS_samples_07142019_1.rds")

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

################################ plots ##############################################
################ TP ####################################
library(plotrix)
par(mfrow=c(1,2))
g<-read.csv("tp.csv")
par(mar=c(5,10,2,0))
plotCI(x=g$Estimates,y=c(4.2,4,3.2,3,2.2,2,1.2,1),li=g$HDIL,ui=g$HDIU, ylim=c(0.5,4.5),main="Total Phosphorus",  
       pch=c(1,16,1,16,1,16,1,16),xlim=c(-5,5), err="x", cex=2,cex.lab=2,xlab= "Estimate", ylab= " ", axes=FALSE)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.4,xpd = TRUE, at=c(4.5,4.1,3.1,2.1,1.1,.5), labels=c(" ","Baseflow", "Forest", "N deposition","Urban"," "), las=2)
abline(v=0, lty=3)
title(ylab="Covariate", line=8, cex.lab=2)
text(x=-5,y=4.5, "A", cex=2)


################ chla ####################################
library(plotrix)
g<-read.csv("chla.csv")
par(mar=c(5,10,2,0))
plotCI(x=g$Estimates,y=c(2.2,2,1.2,1),li=g$HDIL,ui=g$HDIU, ylim=c(.5,2.5),main="Chlorophyll A", pch=c(1,16,1,16),xlim=c(-5,5), err="x", cex=2,cex.lab=2,xlab= "Estimate", ylab= " ", axes=FALSE)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.4, at=c(2.48,2.1,1.1,.49), labels=c("", "Baseflow", "Forest", ""), las=2)
abline(v=0, lty=3)
title(ylab="Covariate", line=8, cex.lab=2)
text(x=-5,y=2.5, "B", cex=2)
legend(x=.2,y=2.5, pch=c(1,16), bty="n", c("Midwest", "New York") , cex=1.5, x.intersp=.3)


