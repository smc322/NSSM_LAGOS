
########################################################## Starting with Secchi in multiple regions dilineated by states
D<-read.csv("LagosSecchi_02Jan2019.csv")
table(D$state_name)
#Connecticut      Illinois       Indiana          Iowa         Maine Massachusetts 
#       650           954          1262          1523         36520           481 
#Michigan     Minnesota      Missouri New Hampshire    New Jersey      New York 
#37404        284079          5543          1294           259         12519 
#Ohio  Pennsylvania  Rhode Island       Vermont     Wisconsin 
#1088           356          9686          9674         68869 


d_1<-subset(D,D$state_name %in% c('Indiana','Illinois', 'Ohio')) #based on Pavel's figures, midwest cities
d_2<-subset(D,D$state_name %in% c('Massachusetts','New Jersey', 'New York', 'Connecticut'))# east coast w big cities
d_3<-subset(D,D$state_name %in% c('Iowa', 'Wisconsin','Illinois'))# upper midwest (corn belt) lots of lakes in WI
d_4<-subset(D,D$state_name %in% c('Vermont', 'New Hampshire', 'Massachusetts','Connecticut'))# Northeast, smaller sample size states
d_5<-subset(D,D$state_name %in% c('Iowa', 'Missouri','Illinois'))# mid midwest
d_6<-subset(D,D$state_name %in% c('Pennsylvania','New York', 'Ohio')) # rust belt (for lack of a better term)

########################################################################################'Indiana','Illinois', 'Ohio'
library(reshape2)

melted<-melt(d_1,id=4:23,na.rm=TRUE) ##### 'Indiana','Illinois', 'Ohio'
#take secchi median per lake
d1<-dcast(subset(melted,melted$variable=='secchi'),
         lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
           ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
           maxdepth + lakeconnection + la_wa_ratio ~ variable,
         mean,drop=TRUE)
d1$log_s<-log(d1$secchi)

summary(d1)
names(d1)

#OLS#
M0<-lm(log_s ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d1) 
summary(M0)
car::qqPlot(M0) 
length(d1$secchi)# 733

#is there spatial corr?
library(geoR)
D_geo1<-as.geodata(na.omit(d1), coords.col=2:3,data.col=21,covar.col = 4:19)
plot(D_geo1) #max dist = 10.441298110 
summary(D_geo1)
#Number of data points: 571 

#variography
v_omni1<-variog(D_geo1,uvec=seq(0,5,0.2),maxdist=5,
               trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                 hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                 ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                 la_wa_ratio)
plot(v_omni1,type='b',main='Indiana,Illinois, Ohio(571 lakes)')
 
##########################################################################Massachusetts','New Jersey', 'New York', 'Connecticut
library(reshape2)

melted<-melt(d_2,id=4:23,na.rm=TRUE) ##### 'Massachusetts','New Jersey', 'New York', 'Connecticut
#take mean secchi  per lake
d2<-dcast(subset(melted,melted$variable=='secchi'),
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio ~ variable,
          mean,drop=TRUE)
d2$log_s<-log(d2$secchi)


#OLS#
M0<-lm(log_s ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d2) 
summary(M0)
car::qqPlot(M0) # not great
length(d2$secchi)#  956

#is there spatial corr?
library(geoR)
D_geo2<-as.geodata(na.omit(d2), coords.col=2:3,data.col=21,covar.col = 4:19)
plot(D_geo2) #diag dist =9.776063030
summary(D_geo2)
#Number of data points: 746 

#variography
v_omni2<-variog(D_geo2,uvec=seq(0,4,0.2),maxdist=4,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
pd2<-plot(v_omni2,type='b',main='Massachusetts,New Jersey, New York, Connecticut(746)')
##########################################################################'Iowa', 'Wisconsin','Illinois'
library(reshape2)

melted<-melt(d_3,id=4:23,na.rm=TRUE) ##### 'Iowa', 'Wisconsin','Illinois'
#take secchi median per lake
d3<-dcast(subset(melted,melted$variable=='secchi'),
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio ~ variable,
          mean,drop=TRUE)
d3$log_s<-log(d3$secchi)

#summary(d3)
#names(d3)

#OLS#
M0<-lm(d3$log1p_s ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d3) 
summary(M0) # doesn't work because log creates -inf, used log1P
d3$sqrt_s<-sqrt(d3$secchi)
d3$log1p_s<-log1p(d3$secchi)
car::qqPlot(M0) # still not great
length(d3$secchi)# 1923 

#is there spatial corr?
library(geoR)
D_geo3<-as.geodata(na.omit(d3), coords.col=2:3,data.col=23,covar.col = 4:19)# using d3$log1p secchi
plot(D_geo3) #max dist =9.8254262 
summary(D_geo3)
#Number of data points: 1833 

#variography
v_omni3<-variog(D_geo3,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni3,type='b',main='Iowa, Wisconsin,Illinois(1833 lakes)')

##########################################################################'Vermont', 'New Hampshire', 'Massachusetts','Connecticut'
library(reshape2)

melted<-melt(d_4,id=4:23,na.rm=TRUE) ##### 'Vermont', 'New Hampshire', 'Massachusetts','Connecticut'
#take secchi mean per lake
d4<-dcast(subset(melted,melted$variable=='secchi'),
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio ~ variable,
          mean,drop=TRUE)
d4$log_s<-log(d4$secchi)

#summary(d4)
#names(d4)

#OLS#
M0<-lm(d4$log_s ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d4) 
summary(M0) # 
car::qqPlot(M0) # not great, sqrt might be better
length(d4$secchi)# 784

#is there spatial corr?
library(geoR)
D_geo4<-as.geodata(na.omit(d4), coords.col=2:3,data.col=21,covar.col = 4:19)# using log
plot(D_geo4)
summary(D_geo4)# max dist = 4.825048247
#Number of data points: 674 

#variography
v_omni4<-variog(D_geo4,uvec=seq(0,2.5,0.2),maxdist=2.5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni4,type='b',main='Vermont, New Hampshire, Massachusetts,Connecticut(674 lakes)')

##########################################################################'Iowa', 'Missouri','Illinois'
library(reshape2)

melted<-melt(d_5,id=4:23,na.rm=TRUE) ##### 'Iowa', 'Missouri','Illinois'
#take secchi mean per lake
d5<-dcast(subset(melted,melted$variable=='secchi'),
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio ~ variable,
          mean,drop=TRUE)
d5$log_s<-log(d5$secchi)

#OLS#
M0<-lm(d5$log_s ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d5) 
summary(M0) # 
car::qqPlot(M0) # ok
length(d5$secchi)# 525

#is there spatial corr?
library(geoR)
D_geo5<-as.geodata(na.omit(d5), coords.col=2:3,data.col=21,covar.col = 4:19)# using log
plot(D_geo5)
summary(D_geo5)# max dist = 9.351287424 
#Number of data points: 456

#variography
v_omni5<-variog(D_geo5,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni5,type='b',main='Iowa, Missouri,Illinois(456 lakes)')

##########################################################################'Pennsylvania','New York', 'Ohio'
library(reshape2)

melted<-melt(d_6,id=4:23,na.rm=TRUE) ##### 'Pennsylvania','New York', 'Ohio'
#take secchi mean per lake
d6<-dcast(subset(melted,melted$variable=='secchi'),
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio ~ variable,
          mean,drop=TRUE)
d6$log_s<-log(d6$secchi)


#OLS#
M0<-lm(d6$log_s ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d6) 
summary(M0) # 
car::qqPlot(M0) # ok
length(d6$secchi)#  761

#is there spatial corr?
library(geoR)
D_geo6<-as.geodata(na.omit(d6), coords.col=2:3,data.col=21,covar.col = 4:19)# using log
plot(D_geo6)
summary(D_geo6)# max dist = 12.531754132 
#Number of data points: 442 

#variography
v_omni6<-variog(D_geo6,uvec=seq(0,6,0.2),maxdist=6,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni6,type='b',main='Pennsylvania,New York, Ohio(442 lakes)')

################################################################################# All Secchi together

par(mfrow=c(2,3))
plot(v_omni1,type='b',main='Sec-IN,IL, OH (571)')
plot(v_omni2,type='b',main='Sec-MA,NJ, NY, CT(746), abnormal')
plot(v_omni3,type='b',main='Sec-IA, WI,IL (1833), abnormal')
plot(v_omni4,type='b',main='Sec-VT, NH, MA,CT (674)')
plot(v_omni5,type='b',main='Sec-IA, MO,IL (456)')
plot(v_omni6,type='b',main='Sec-PA,NY, OH (442)')

########## working on other Responses################################################################# TP

setwd("C:/Users/Charlotte/OneDrive/Charlotte Work 2019/UWYO")

D<-read.csv("LagosNonSecchi_16Jan2019.csv")
table(D$state_name)
#Connecticut      Illinois       Indiana          Iowa         Maine Massachusetts 
#      895          1307          1258          1550          9091           530 
#Michigan     Minnesota      Missouri New Hampshire    New Jersey      New York 
#8305         47442          5499          5492           283         16372 
#Ohio  Pennsylvania  Rhode Island       Vermont     Wisconsin 
#1234           469          6863          7160         22395 



d_1<-subset(D,D$state_name %in% c('Indiana','Illinois', 'Ohio')) #based on Pavel's figures, midwest cities
d_2<-subset(D,D$state_name %in% c('Massachusetts','New Jersey', 'New York', 'Connecticut'))# east coast w big cities
d_3<-subset(D,D$state_name %in% c('Iowa', 'Wisconsin','Illinois'))# upper midwest, lots of lakes
d_4<-subset(D,D$state_name %in% c('Vermont', 'New Hampshire', 'Massachusetts','Connecticut'))# Northeast, smaller sample size states
d_5<-subset(D,D$state_name %in% c('Iowa', 'Missouri','Illinois'))# mid midwest
d_6<-subset(D,D$state_name %in% c('Pennsylvania','New York', 'Ohio')) # rust belt (for lack of a better term)

########################################################################################'Indiana','Illinois', 'Ohio'
library(reshape2)

melted<-melt(d_1,id=c(1:2,4,7:26)) #Need to adjust so response var, sampleyear and month are left out
levels(melted$variable)
#take mean TP per lake
d1<-dcast(subset(melted,melted$variable=='tp'), #<<<--pick outcome variable
         lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
           ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
           maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d1$log_p<-log(d1$tp)

summary(d1)
names(d1)

#OLS#
M0<-lm(log_p ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d1) ##### 
summary(M0)
car::qqPlot(M0) 
length(d1$log_p)#  745

#is there spatial corr?
library(geoR)
D_geo1<-as.geodata(na.omit(d1), coords.col=2:3,data.col=22,covar.col=4:19)
plot(D_geo1) 
summary(D_geo1)# Maxlength = 9.662869876 
#Number of data points: 406 

#variography
v_omni1<-variog(D_geo1,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni1,type='b',main='Indiana,Illinois, Ohio(406)')

##########################################################################Massachusetts','New Jersey', 'New York', 'Connecticut
library(reshape2)

melted<-melt(d_2,id=c(1:2,4,7:26)) ##### 'Massachusetts','New Jersey', 'New York', 'Connecticut
#take mean TP per lake
d2<-dcast(subset(melted,melted$variable=='tp'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d2$log_p<-log(d2$tp)


#OLS#
M0<-lm(log_p ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d2) 
summary(M0)
car::qqPlot(M0) # not perfect
length(d2$tp)#  1191

#is there spatial corr?
library(geoR)
D_geo2<-as.geodata(na.omit(d2), coords.col=2:3,data.col=22,covar.col=4:19)
plot(D_geo2) 
summary(D_geo2) # max dist =9.222398994 
#Number of data points: 330 

#variography
v_omni2<-variog(D_geo2,uvec=seq(0,4.5,0.2),maxdist=4.5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni2,type='b',main='Massachusetts,New Jersey, New York, Connecticut(330 lakes)')
##########################################################################'Iowa', 'Wisconsin','Illinois'
library(reshape2)

melted<-melt(d_3,id=c(1:2,4,7:26),na.rm=TRUE) ##### 'Iowa', 'Wisconsin','Illinois'
#take mean TP per lake
d3<-dcast(subset(melted,melted$variable=='tp'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d3$log_p<-log(d3$tp)


#OLS#
M0<-lm(d3$log_p ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d3) 
summary(M0) 
car::qqPlot(M0) # not great
length(d3$log_p)#  1633

#is there spatial corr?
library(geoR)
D_geo3<-as.geodata(na.omit(d3), coords.col=2:3,data.col=22,covar.col = 4:19)# 
plot(D_geo3)
summary(D_geo3) #max dist =9.889709408 
#Number of data points: 1548 

#variography
v_omni3<-variog(D_geo3,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni3,type='b',main='Iowa, Wisconsin,Illinois(1548 lakes)')

##########################################################################'Vermont', 'New Hampshire', 'Massachusetts','Connecticut'
library(reshape2)

melted<-melt(d_4,id=c(1:2,4,7:26),na.rm=TRUE) ##### 'Vermont', 'New Hampshire', 'Massachusetts','Connecticut'
#take mean TP per lake
d4<-dcast(subset(melted,melted$variable=='tp'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d4$log_p<-log(d4$tp)
d4$sqrt_p<-sqrt(d4$tp)


#OLS#
M0<-lm(d4$log_p ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d4) 
summary(M0) # 
car::qqPlot(M0) # REALLY BAD, trying sqrt, not helping :(

length(d4$tp)#  1011

#is there spatial corr?
library(geoR)
D_geo4<-as.geodata(na.omit(d4), coords.col=2:3,data.col=22,covar.col = 4:19)# using log
plot(D_geo4)
summary(D_geo4)# max dist = 4.825048247 
#Number of data points: 891 

#variography
v_omni4<-variog(D_geo4,uvec=seq(0,2.5,0.2),maxdist=2.5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni4,type='b',main='Vermont, New Hampshire, Massachusetts,Connecticut(891 lakes)')

##########################################################################'Iowa', 'Missouri','Illinois'
library(reshape2)

melted<-melt(d_5,id=c(1:2,4,7:26),na.rm=TRUE) ##### 'Iowa', 'Missouri','Illinois'
#take mean TP per lake
d5<-dcast(subset(melted,melted$variable=='tp'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d5$log_p<-log(d5$tp)


#OLS#
M0<-lm(d5$log_p ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d5) 
summary(M0) # 
car::qqPlot(M0) # ok
length(d5$tp)# 533

#is there spatial corr?
library(geoR)
D_geo5<-as.geodata(na.omit(d5), coords.col=2:3,data.col=22,covar.col = 4:19)# using log
plot(D_geo5)
summary(D_geo5)# max dist = 9.351287424  
#Number of data points: 462 

#variography
v_omni5<-variog(D_geo5,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni5,type='b',main='Iowa, Missouri,Illinois(462 lakes)')

##########################################################################'Pennsylvania','New York', 'Ohio'
library(reshape2)

melted<-melt(d_6,id=c(1:2,4,7:26),na.rm=TRUE) ##### 'Pennsylvania','New York', 'Ohio'
#take mean TP per lake
d6<-dcast(subset(melted,melted$variable=='tp'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d6$log_p<-log(d6$tp)
d6$log1p_p<-log1p(d6$tp)


#OLS#
M0<-lm(d6$log1p_p ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d6) 
summary(M0) # need to use log1p
car::qqPlot(M0) # not perfect 
length(d6$tp)#  671

#is there spatial corr?
library(geoR)
D_geo6<-as.geodata(na.omit(d6), coords.col=2:3,data.col=23,covar.col = 4:19)# using log
plot(D_geo6)
summary(D_geo6)# max dist = 12.531754132 
#Number of data points: 368 


#variography
v_omni6<-variog(D_geo6,uvec=seq(0,6,0.2),maxdist=6,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni6,type='b',main='Pennsylvania,New York, Ohio(368 lakes)')

################################################################################# All TP together

par(mfrow=c(2,3))
plot(v_omni1,type='b',main='TP-IN,IL, OH (406)')
plot(v_omni2,type='b',main='TP-MA,NJ, NY, CT (330)')
plot(v_omni3,type='b',main='TP-IA, WI,IL (1548), abnormal')
plot(v_omni4,type='b',main='TP-VT, NH, MA,CT (891),abnormal')
plot(v_omni5,type='b',main='TP-IA, MO,IL (462)')
plot(v_omni6,type='b',main='TP-PA,NY, OH (368), abnormal')

########## working on other Responses########################################################### TN_all


d_1<-subset(D,D$state_name %in% c('Indiana','Illinois', 'Ohio')) #based on Pavel's figures, midwest cities
d_2<-subset(D,D$state_name %in% c('Massachusetts','New Jersey', 'New York', 'Connecticut'))# east coast w big cities
d_3<-subset(D,D$state_name %in% c('Iowa', 'Wisconsin','Illinois'))# upper midwest, lots of lakes
d_4<-subset(D,D$state_name %in% c('Vermont', 'New Hampshire', 'Massachusetts','Connecticut'))# Northeast, smaller sample size states
d_5<-subset(D,D$state_name %in% c('Iowa', 'Missouri','Illinois'))# mid midwest
d_6<-subset(D,D$state_name %in% c('Pennsylvania','New York', 'Ohio')) # rust belt (for lack of a better term)

########################################################################################'Indiana','Illinois', 'Ohio'
library(reshape2)

melted<-melt(d_1,id=c(1:4,8:26)) #Need to adjust so response var, sampleyear and month are left out
levels(melted$variable)
#take mean tn_all per lake
d1<-dcast(subset(melted,melted$variable=='tn_all'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d1$log_n<-log(d1$tn_all)

summary(d1)
names(d1)

#OLS#
M0<-lm(log_n ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d1) ##### 
summary(M0)
car::qqPlot(M0) # ok
length(d1$log_n)#  745

#is there spatial corr?
library(geoR)
D_geo1<-as.geodata(na.omit(d1), coords.col=2:3,data.col=22,covar.col=4:19)
plot(D_geo1) 
summary(D_geo1)# Maxlength = 9.912227887 
#Number of data points: 367 

#variography
v_omni1<-variog(D_geo1,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni1,type='b',main='Indiana,Illinois, Ohio(367)')

##########################################################################Massachusetts','New Jersey', 'New York', 'Connecticut
library(reshape2)

melted<-melt(d_2,id=c(1:4,8:26)) ##### 'Massachusetts','New Jersey', 'New York', 'Connecticut
#take mean tn_all per lake
d2<-dcast(subset(melted,melted$variable=='tn_all'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d2$log_n<-log(d2$tn_all)


#OLS#
M0<-lm(log_n ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d2) 
summary(M0)
car::qqPlot(M0) # not perfect, but OK
length(d2$tn_all)#  1191

#is there spatial corr?
library(geoR)
D_geo2<-as.geodata(na.omit(d2), coords.col=2:3,data.col=22,covar.col=4:19)
plot(D_geo2) 
summary(D_geo2) # max dist =9.118595674 
#Number of data points: 263 

#variography
v_omni2<-variog(D_geo2,uvec=seq(0,4.5,0.2),maxdist=4.5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni2,type='b',main='Massachusetts,New Jersey, New York, Connecticut(263)')
##########################################################################'Iowa', 'Wisconsin','Illinois'
library(reshape2)

melted<-melt(d_3,id=c(1:4,8:26),na.rm=TRUE) ##### 'Iowa', 'Wisconsin','Illinois'
#take mean tn_all per lake
d3<-dcast(subset(melted,melted$variable=='tn_all'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d3$log_n<-log(d3$tn_all)
d3$sqrt_n<-sqrt(d3$tn_all)


#OLS#
M0<-lm(d3$sqrt_n~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d3) 
summary(M0) 
car::qqPlot(M0) # awful, better with sqrt
length(d3$log_n)#  1633

#is there spatial corr?
library(geoR)
D_geo3<-as.geodata(na.omit(d3), coords.col=2:3,data.col=24,covar.col = 4:19)# using sqrt_n
plot(D_geo3)
summary(D_geo3) #max dist =9.755039220  
#Number of data points: 809 


#variography
v_omni3<-variog(D_geo3,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni3,type='b',main='Iowa, Wisconsin,Illinois(809), abnormal')

##########################################################################'Vermont', 'New Hampshire', 'Massachusetts','Connecticut'
library(reshape2)

melted<-melt(d_4,id=c(1:4,8:26),na.rm=TRUE) ##### 'Vermont', 'New Hampshire', 'Massachusetts','Connecticut'
#take mean tn_all per lake
d4<-dcast(subset(melted,melted$variable=='tn_all'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d4$log_n<-log(d4$tn_all)
d4$sqrt_n<-sqrt(d4$tn_all)


#OLS#
M0<-lm(d4$log_n ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d4) 
summary(M0) # 
car::qqPlot(M0) # not great, but not awful

length(d4$tn_all)#  285

#is there spatial corr?
library(geoR)
D_geo4<-as.geodata(na.omit(d4), coords.col=2:3,data.col=22,covar.col = 4:19)# using log
plot(D_geo4)
summary(D_geo4)# max dist = 4.340333004
#Number of data points: 243 

#variography
v_omni4<-variog(D_geo4,uvec=seq(0,2.5,0.2),maxdist=2.5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni4,type='b',main='Vermont, New Hampshire, Massachusetts,Connecticut(243)')

##########################################################################'Iowa', 'Missouri','Illinois'
library(reshape2)

melted<-melt(d_5,id=c(1:4,8:26),na.rm=TRUE) ##### 'Iowa', 'Missouri','Illinois'
#take mean tn_all per lake
d5<-dcast(subset(melted,melted$variable=='tn_all'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d5$log_n<-log(d5$tn_all)

#OLS#
M0<-lm(d5$log_n ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d5) 
summary(M0) # 
car::qqPlot(M0) # a bit abnormal, but prob ok
length(d5$tn_all)# 489

#is there spatial corr?
library(geoR)
D_geo5<-as.geodata(na.omit(d5), coords.col=2:3,data.col=22,covar.col = 4:19)# using log
plot(D_geo5)
summary(D_geo5)# max dist = 9.351287424  
#Number of data points: 428 

#variography
v_omni5<-variog(D_geo5,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni5,type='b',main='Iowa, Missouri,Illinois(428 lakes)')

##########################################################################'Pennsylvania','New York', 'Ohio'
library(reshape2)

melted<-melt(d_6,id=c(1:4,8:26),na.rm=TRUE) ##### 'Pennsylvania','New York', 'Ohio'
#take mean tn_all per lake
d6<-dcast(subset(melted,melted$variable=='tn_all'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d6$log_n<-log(d6$tn_all)

#OLS#
M0<-lm(d6$log_n ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d6) 
summary(M0) # need to use log1p
car::qqPlot(M0) # not perfect, but prob ok
length(d6$tn_all)#  606

#is there spatial corr?
library(geoR)
D_geo6<-as.geodata(na.omit(d6), coords.col=2:3,data.col=22,covar.col = 4:19)# using log
plot(D_geo6)
summary(D_geo6)# max dist = 12.531754132 
#Number of data points: 351 


#variography
v_omni6<-variog(D_geo6,uvec=seq(0,6,0.2),maxdist=6,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni6,type='b',main='Pennsylvania,New York, Ohio(351)')

################################################################################# All TN together

par(mfrow=c(2,3))
plot(v_omni1,type='b',main='TN-IN,IL, OH (367)')
plot(v_omni2,type='b',main='TN-MA,NJ, NY, CT (263)')
plot(v_omni3,type='b',main='TN-IA, WI,IL (809), abnormal')
plot(v_omni4,type='b',main='TN-VT, NH, MA,CT (243)')
plot(v_omni5,type='b',main='TN-IA, MO,IL (428)')
plot(v_omni6,type='b',main='TN-PA,NY, OH (351)')

########## working on other Responses########################################################### chla


d_1<-subset(D,D$state_name %in% c('Indiana','Illinois', 'Ohio')) #based on Pavel's figures, midwest cities
d_2<-subset(D,D$state_name %in% c('Massachusetts','New Jersey', 'New York', 'Connecticut'))# east coast w big cities
d_3<-subset(D,D$state_name %in% c('Iowa', 'Wisconsin','Illinois'))# upper midwest, lots of lakes
d_4<-subset(D,D$state_name %in% c('Vermont', 'New Hampshire', 'Massachusetts','Connecticut'))# Northeast, smaller sample size states
d_5<-subset(D,D$state_name %in% c('Iowa', 'Missouri','Illinois'))# mid midwest
d_6<-subset(D,D$state_name %in% c('Pennsylvania','New York', 'Ohio')) # rust belt (for lack of a better term)

########################################################################################'Indiana','Illinois', 'Ohio'
library(reshape2)

melted<-melt(d_1,id=c(1,8:26)) #Need to adjust so response var, sampleyear and month are left out
levels(melted$variable)
#take mean chla per lake
d1<-dcast(subset(melted,melted$variable=='chla'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d1$log_chla<-log(d1$chla)

summary(d1)
names(d1)

#OLS#
M0<-lm(log_chla ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d1) ##### 
summary(M0)
car::qqPlot(M0) # ok
length(d1$chla)#  745

#is there spatial corr?
library(geoR)
D_geo1<-as.geodata(na.omit(d1), coords.col=2:3,data.col=22,covar.col=4:19)
plot(D_geo1) 
summary(D_geo1)# Maxlength = 9.657897826 
#Number of data points: 175 

#variography
v_omni1<-variog(D_geo1,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni1,type='b',main='Indiana,Illinois, Ohio(175)')

##########################################################################Massachusetts','New Jersey', 'New York', 'Connecticut
library(reshape2)

melted<-melt(d_2,id=c(1,8:26)) ##### 'Massachusetts','New Jersey', 'New York', 'Connecticut
#take mean chla per lake
d2<-dcast(subset(melted,melted$variable=='chla'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d2$log_chla<-log(d2$chla)

#summary(d2)
#names(d2)

#OLS#
M0<-lm(log_chla ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d2) 
summary(M0)
car::qqPlot(M0) # OK
length(d2$chla)#  1191

#is there spatial corr?
library(geoR)
D_geo2<-as.geodata(na.omit(d2), coords.col=2:3,data.col=22,covar.col=4:19)
plot(D_geo2) 
summary(D_geo2) # max dist =9.393272905
#Number of data points: 409 

#variography
v_omni2<-variog(D_geo2,uvec=seq(0,4.5,0.2),maxdist=4.5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni2,type='b',main='Massachusetts,New Jersey, New York, Connecticut(409)')
##########################################################################'Iowa', 'Wisconsin','Illinois'
library(reshape2)

melted<-melt(d_3,id=c(1,8:26),na.rm=TRUE) ##### 'Iowa', 'Wisconsin','Illinois'
#take mean chla per lake
d3<-dcast(subset(melted,melted$variable=='chla'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d3$log_chla<-log(d3$chla)
#d3$log1p_s<-log1p(d3$secchi)
#summary(d3)
#names(d3)

#OLS#
M0<-lm(d3$log_chla~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d3) 
summary(M0) 
car::qqPlot(M0) #not awful, but not perfect
length(d3$log_chla)#  1297

#is there spatial corr?
library(geoR)
D_geo3<-as.geodata(na.omit(d3), coords.col=2:3,data.col=22,covar.col = 4:19)# using sqrt_n
plot(D_geo3)
summary(D_geo3) #max dist =9.889709408 
#Number of data points: 1244 


#variography
v_omni3<-variog(D_geo3,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni3,type='b',main='Iowa, Wisconsin,Illinois(1244), a bit abnormal')

##########################################################################'Vermont', 'New Hampshire', 'Massachusetts','Connecticut'
library(reshape2)

melted<-melt(d_4,id=c(1,8:26),na.rm=TRUE) ##### 'Vermont', 'New Hampshire', 'Massachusetts','Connecticut'
#take mean chla per lake
d4<-dcast(subset(melted,melted$variable=='chla'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d4$log_chla<-log(d4$chla)
#d4$sqrt_n<-sqrt(d4$chla)
#d4$log1p_s<-log1p(d4$secchi)

#summary(d4)
#names(d4)

#OLS#
M0<-lm(d4$log_chla ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d4) 
summary(M0) # 
car::qqPlot(M0) # outliers?

length(d4$chla)#  346

#is there spatial corr?
library(geoR)
D_geo4<-as.geodata(na.omit(d4), coords.col=2:3,data.col=22,covar.col = 4:19)# using log
plot(D_geo4)
summary(D_geo4)# max dist = 4.691734775 
#Number of data points: 289 

#variography
v_omni4<-variog(D_geo4,uvec=seq(0,2.5,0.2),maxdist=2.5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni4,type='b',main='Vermont, New Hampshire, Massachusetts,Connecticut(289)')

##########################################################################'Iowa', 'Missouri','Illinois'
library(reshape2)

melted<-melt(d_5,id=c(1,8:26),na.rm=TRUE) ##### 'Iowa', 'Missouri','Illinois'
#take mean chla per lake
d5<-dcast(subset(melted,melted$variable=='chla'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d5$log_chla<-log(d5$chla)
#d5$sqrt_s<-sqrt(d5$secchi)
#d5$log1p_s<-log1p(d5$secchi)

#summary(d5)
#names(d5)

#OLS#
M0<-lm(d5$log_chla ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d5) 
summary(M0) # 
car::qqPlot(M0) # ok
length(d5$chla)# 488

#is there spatial corr?
library(geoR)
D_geo5<-as.geodata(na.omit(d5), coords.col=2:3,data.col=22,covar.col = 4:19)# using log
plot(D_geo5)
summary(D_geo5)# max dist = 9.028899459  
#Number of data points: 428

#variography
v_omni5<-variog(D_geo5,uvec=seq(0,5,0.2),maxdist=5,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni5,type='b',main='Iowa, Missouri,Illinois(428 lakes)')

##########################################################################'Pennsylvania','New York', 'Ohio'
library(reshape2)

melted<-melt(d_6,id=c(1,8:26),na.rm=TRUE) ##### 'Pennsylvania','New York', 'Ohio'
#take mean chla per lake
d6<-dcast(subset(melted,melted$variable=='chla'), #<<<--pick outcome variable
          lagoslakeid + nhd_long + nhd_lat + hu8_baseflow + hu8_no3depo + hu8_totalndepo + hu8_MAP + hu8_MAT + hu8_runoff +
            ws_urban + ws_rowcrop + ws_pasture + ws_forest + ws_wetland + ws_openh20 + lake_area +
            maxdepth + lakeconnection + la_wa_ratio + state_name ~ variable,
          mean,drop=TRUE)
d6$log_chla<-log(d6$chla)
#d6$sqrt_s<-sqrt(d6$secchi)
#d6$log1p_p<-log1p(d6$tp)

#summary(d6)
#names(d6)

#OLS#
M0<-lm(d6$log_chla ~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
         hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
         ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
         la_wa_ratio,data=d6) 
summary(M0) # need to use log1p
car::qqPlot(M0) # not perfect, but prob ok
length(d6$chla)#  758

#is there spatial corr?
library(geoR)
D_geo6<-as.geodata(na.omit(d6), coords.col=2:3,data.col=22,covar.col = 4:19)# using log
plot(D_geo6)
summary(D_geo6)# max dist = 12.53175413 
#Number of data points:  438  


#variography
v_omni6<-variog(D_geo6,uvec=seq(0,6,0.2),maxdist=6,
                trend=~ hu8_baseflow +hu8_no3depo+hu8_totalndepo+ hu8_MAP +
                  hu8_MAT+ hu8_runoff+ws_urban+ws_rowcrop+ws_pasture+ws_forest+
                  ws_wetland+ws_openh20+lake_area + maxdepth + lakeconnection+
                  la_wa_ratio)
plot(v_omni6,type='b',main='Pennsylvania,New York, Ohio(438)')

################################################################################# All chla together

par(mfrow=c(2,3))
plot(v_omni1,type='b',main='Chla-IN,IL, OH (175)')
plot(v_omni2,type='b',main='Chla-MA,NJ, NY, CT (409)')
plot(v_omni3,type='b',main='Chla-IA, WI,IL (1244),a bit abnormal')
plot(v_omni4,type='b',main='Chla-VT, NH, MA,CT (289) abnormal')
plot(v_omni5,type='b',main='Chla-IA, MO,IL (428)')
plot(v_omni6,type='b',main='Chla-PA,NY, OH (438)')

