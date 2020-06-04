
# Code for running gams on lake TP in Northeastern region and creating Figure 2. 

#
################################################################################ TP NY ###############################################
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


sm_X1<-gam(tp ~ s(nhd_long,nhd_lat,bs='tp'),
           
           data=na.omit(d),method="REML")

gam.check(sm_X1)

summary(sm_X1)
####*** GAM exploratory analysis ***###
# interior knots
int_knots<-data.frame(fields::cover.design(cbind(d$nhd_long,d$nhd_lat),10,max.loop=50)$design)
names(int_knots)<-c("x","y")
library(mgcv)
sm_tp0<-gam(ws_urban~ s(nhd_long,nhd_lat,k=110,bs='tp'),
            data=na.omit(d),method="REML")
summary(sm_tp0)


###################### PLOTS ############################
install.packages("plotrix")
library(plotrix)
g<-read.csv("plot_gams.csv")
par(mfrow=c(3,2), mar=c(4,5,2,1))
plotCI(x=g$Dev_exp,y=g$mdt, uiw=g$sedt,pch=c(0,1,2,15,16,17), cex=2,cex.lab=2, xlab= "", xlim=c(0,100), ylab= "Divergent Transitions", axes=FALSE)
axis(1, cex.axis=2, labels = FALSE)
axis(2, cex.axis=2)
text(x=4, y=240, "A", cex=2.2)

par(mar=c(4,3,2,3))
plotCI(x=g$edf,y=g$mdt,uiw=g$sedt, xlab= "", pch=c(0,1,2,15,16,17), cex=2,cex.lab=2, ylab= "", axes=FALSE)
axis(1, cex.axis=2, labels = FALSE)
axis(2, cex.axis=2, labels = FALSE)
text(x=11.5, y=240, "B", cex=2.2)
legend('topright', pch=c(0,1,2,15,16,17), bty="n", cex=1.5, c( "Urban - Midwest", "Urban - New York", "Forest - Midwest","Forest - New York", "N deposition - Midwest", 
                                                               "N deposition - New York"))

par(mar=c(4,5,2,1))
plotCI(x=g$Dev_exp,y=g$mneff,uiw=g$seneff, pch=c(0,1,2,15,16,17), cex=2,cex.lab=2,xlab= "",xlim=c(0,100), ylab= "Fewest Effective Samples", axes=FALSE)
axis(1, cex.axis=2, labels = FALSE)
axis(2, cex.axis=2)
text(x=4, y=290, "C", cex=2.2)

par(mar=c(4,3,2,3))
plotCI(x=g$edf,y=g$mneff,uiw=g$seneff, xlab= "", pch=c(0,1,2,15,16,17), cex=2,cex.lab=2,ylab= "", axes=FALSE)
axis(1, cex.axis=2, labels = FALSE)
axis(2, cex.axis=2, labels = FALSE)
text(x=11.5, y=290, "D", cex=2.2)

par(mar=c(5,5,2,1))
plotCI(x=g$Dev_exp,y=g$mrhat,uiw=g$serhat, pch=c(0,1,2,15,16,17), cex=2,cex.lab=2,xlab= "Deviance Explained (%)",xlim=c(0,100), ylab= "Highest r hat", axes=FALSE)
axis(1, cex.axis=2)
axis(2, cex.axis=2)
text(x=4, y=1.78, "E", cex=2.2)
par(mar=c(5,3,2,3))
plotCI(x=g$edf,y=g$mrhat,uiw=g$serhat, xlab= "EDF",pch=c(0,1,2,15,16,17), cex=2,cex.lab=2, ylab= "", axes=FALSE)
axis(1, cex.axis=2)
axis(2, cex.axis=2, labels = FALSE)
text(x=11.5, y=1.78, "F", cex=2.2)



