
#Code for generating figures 3 and 4 showing the residual spatial process of 
#lake TP in the midwestern and northeastern regions.

library(rstan)
library(plot3D)
library(geoR)
library(car)
library(MBA)
library(plot3D)
library(maps)
library(matrixStats)
library(reshape2)

inv_logit<-function(x){exp(x)/(1+exp(x))}


#Example code for TP in the midwest region. Load different stan output for 
#other regions or response variables can be substituted and are included at the bottom

#load stan fit & coordinates
stan_mod<-readRDS("Data/ModelOutput/tp_mw_NS_samples_07042019_3.rds")
D<-read.csv("Data/LAGOS_Summaries/LagosNonSecchi_16Jan2019.csv")
d_0<-na.omit(subset(D,D$state_name %in% c('Illinois','Iowa', 'Wisconsin')))[,-1] #removes *ALL* NAs

#take TP median per lake over time
melted<-melt(d_0,id=c(1:2,9:25),na.rm=TRUE)
d<-dcast(subset(melted,melted$variable=='tp'), nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + 
           hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff + ws_urban + ws_rowcrop + ws_pasture +
           ws_forest + ws_wetland + ws_openh20 + lake_area + maxdepth + lakeconnection + 
           la_wa_ratio ~ variable, median,drop=TRUE)
summary(d)
d$log_tp<-log(d$tp)
rm(d_0,D,melted)

D_geo<-as.geodata(d, coords.col=1:2,data.col=20,covar.col = 3:19)
h_xy<-as.matrix(D_geo$coords)  
S<-nrow(h_xy)

M0<-lm(log_tp ~ scale(hu8_baseflow) + scale(ws_forest) + scale(maxdepth)+lakeconnection ,data=d)
summary(M0)

###### Recovering spatial random effect ##############
#*N-S spec*#
X_L<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #strech
X_T<-model.matrix(~ scale(hu8_baseflow) + scale(ws_forest),data=d) #angle of rotation

#Exctract necessary MCMC samples
b_mu<-extract(stan_mod,"beta_mu")[[1]]
b_L<-extract(stan_mod,"beta_L")[[1]]
b_T<-extract(stan_mod,"beta_T")[[1]]
alpha_sq<-extract(stan_mod,"alpha_sq")[[1]]
nug_sq<-extract(stan_mod,"nug_sq")[[1]]
rho<-extract(stan_mod,"rho")[[1]]

#compute S lambdas and thetas
lamb_mat<-matrix(NA,nrow=S,ncol=nrow(b_L)) #all matrices are S x MCMC samps 
ang_mat<-matrix(NA,nrow=S,ncol=nrow(b_T))
for(i in 1:nrow(b_L)){
  lamb_mat[,i]<-exp(X_L%*%b_L[i,])
  ang_mat[,i]<-inv_logit(X_T%*%b_T[i,])
}

#angle & lambda posterior per location
#set.seed(101119851)
par(mfrow=c(4,4))
loc<-sample(1:nrow(h_xy),size=16,replace = F) #select 25 locations at random
for(i in 1:length(loc)) { hist(lamb_mat[loc[i],],xlab='Lambda',xlim = c(0,40),
                               main=paste("Loc ",loc[i]))
  abline(v=c(mean(lamb_mat[loc[i],]),median(lamb_mat[loc[i],])),lty=1:2,lwd=2,col='red')}
for(i in 1:length(loc)) { hist(180*ang_mat[loc[i],],xlab='Angle',xlim=c(0,180),
                               main=paste("Loc ",loc[i]))
  abline(v=c(mean(180*ang_mat[loc[i],]),median(180*ang_mat[loc[i],])),lty=1:2,lwd=2,col='red')}

#pick one summary statistic - medians for lambda and modes for angles
lamb_S<-rowMedians(lamb_mat)
ang_S<-rowMedians(ang_mat)

#Compute N-S kernel function at each location
#N-S kernel function to rotate and scale at each location
kern_func_H2<-function(theta, lambda){
  G<-matrix(c(lambda*cos(theta),lambda*sin(theta),
              -sin(theta), cos(theta)), ncol=2,byrow=TRUE)
  H<- crossprod(G)
  return(H)
}

#Compute array of 2 x 2 covariances (ellipses) at each location
Sig_array<-array(NA,dim=c(2,2,S))
for(i in 1:S) { 
  if(lamb_S[i]<=1.0){      
    Sig_array[,,i]<-kern_func_H2((pi/2) + (pi*ang_S[i]), lamb_S[i]) 
  } else {
    Sig_array[,,i]<-kern_func_H2(pi*ang_S[i], lamb_S[i]) }} 

#plot ellipses
#using median angles: 101119851, 101119853 w 40 pts OK#### 11191991 for NY
set.seed(11191992)
n_pts<-50
example_pts<-fields::cover.design(h_xy,n_pts,max.loop=50)

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), heights=c(1.5,1.5,1,1))
par(oma=c(0,0,0,0), mar=c(2,2,1,1))
rm<-7
NY<-map("state",c("Illinois", "Iowa", "Wisconsin"),fill=TRUE,plot=F)
scatter2D(x=d$nhd_long,y=d$nhd_lat,colvar=resid(M0),pch=16,
          cex=1.3,cex.main=1.3,cex.lab=1.8,col=jet.col(32),
          ylim=c(38,48.6),xlim=c(-99,-84),
          main='a) OLS Residuals',xlab='',ylab='', bty='n',yaxt='n',
          ann=TRUE, xaxt='n', colkey = list(side=1))
mtext("Midwest", side=2, cex=1.8)
lines(NY,lwd=2)

#right panel
sp.d <- cbind(d$nhd_long,d$nhd_lat,d$ws_forest)
surf <- mba.surf(sp.d,75,75,extend = F)$xyz.est #extend=T extrapolates
scatter2D(x=d$nhd_long,y=d$nhd_lat,colvar=d$ws_forest,pch=16,cex=0,
          ylim=c(38,48.6),xlim=c(-99,-84),
          cex.main=1.3,col=topo.colors(32,0.3),
          main='b) % Forest and correlation kernels',xlab='',ylab='',
          bty='n', axes=FALSE, colkey = list(side=1))
lines(NY,lwd=2)
image2D(surf, xlab='',ylab='',add=TRUE,colkey=FALSE,
        alpha=0.4,col=topo.colors(32,0.3))
lines(NY,lwd=2)

j<-c(1:nrow(example_pts$design))
for(i in j){
  car::ellipse(center=as.numeric( h_xy[example_pts$best.id[i],] ),col='magenta',
               shape=Sig_array[,,as.numeric(example_pts$best.id[i])],lwd=3,
               radius=0.3,add=TRUE,draw=TRUE)
}


##STAN output for other response variables and regions

#stan_mod<-readRDS("Data/ModelOutput/tp_ny_NS_samples_07102019_2.rds")
#stan_mod<-readRDS("Data/ModelOutput/chla_mw_NS_samples_08072019_3.rds")
#stan_mod<-readRDS("Data/ModelOutput/chla_ny_NS_samples_8292019_3.rds")




