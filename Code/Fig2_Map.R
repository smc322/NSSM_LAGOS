library(lubridate)
library(tidyverse)
library(sf)
library(LAGOSNE)
library(maps)
library(vioplot)
library(ggplot2)

#save data files
dat<-na.omit(read.csv("Data/LAGOS_Summaries/LagosNonSecchi_16Jan2019.csv"))
dat.tp.midwest<-filter(dat, state_name== c("Wisconsin", "Iowa", "Illinois")) %>% group_by(lagoslakeid) %>% summarise(tpmed = median(tp, na.rm=T))

lakecords<-unique(dat[,c("lagoslakeid", "nhd_long", "nhd_lat")])

mapdat.mw<-merge(dat.tp.midwest, lakecords, by="lagoslakeid", all.x = T, all.y=F)

dat.tp.ny<-filter(dat, state_name== c("New York")) %>% group_by(lagoslakeid) %>% summarise(tpmed = median(tp, na.rm=T))

mapdat.ny<-merge(dat.tp.ny, lakecords, by="lagoslakeid", all.x = T, all.y = F)

mw<- filter(dat, state_name== c("Wisconsin", "Iowa", "Illinois")) 
ny<- filter(dat, state_name == c("New York"))

covars.mw<-mw[!duplicated(mw$lagoslakeid), ]
covars.mw$region<- "Midwest"
covars.ny<-ny[!duplicated(ny$lagoslakeid), ]
covars.ny$region<- "Northeast"

covars.all<-rbind(covars.mw,covars.ny)

covars.all.rel<- covars.all[,c("hu8_baseflow", "ws_urban", "ws_forest", "hu8_totalndepo", "region")]

keycol<- "covariate"
valuecol<- "value"
gathercols<-c("hu8_baseflow", "ws_urban", "ws_forest", "hu8_totalndepo")

covars.long<-gather_(covars.all.rel, keycol, valuecol, gathercols)

x=1/5
quants.mw<-quantile(mapdat.mw$tpmed, probs=c(x, 2*x, 3*x, 4*x))
quants.ny<-quantile(mapdat.ny$tpmed, probs=c(x, 2*x, 3*x, 4*x))
quants.all<-quantile(rbind(mapdat.mw,mapdat.ny)$tpmed, probs=c(x, 2*x, 3*x, 4*x))

alpha=255

colors<-c(rgb(5, 113, 176,max=255, alpha=alpha),
          rgb(146, 197, 222,max=255, alpha=alpha),
          rgb(247, 247, 247,max=255, alpha=alpha),
          rgb(244, 165, 130,max=255, alpha=alpha),
          rgb(202, 0, 32,max=255, alpha=alpha))

get.col.bins <- function(slopes, alpha=175) {
  z=slopes
  #x=1/7
  #quants<-quantile(z, probs = c(x, 2*x, 3*x, 4*x, 5*x, 6*x), na.rm=T)
  
  ii <- cut(z, breaks = c(-Inf, 13, 23, 60, 110, Inf), 
            include.lowest = TRUE)
  
  #ii <- cut(z, breaks = c(0,250, 500, 750, 1000, 1250, 1500,Inf), 
  #         include.lowest = TRUE)
  
  #purple blue green
  # levels(ii) <- c(rgb(246,239,247,max=255, alpha=alpha),
  #                 rgb(208,209,230,max=255, alpha=alpha),
  #                 rgb(166,189,219,max=255, alpha=alpha),
  #                 rgb(103,169,207,max=255, alpha=alpha),
  #                 rgb(54,144,192,max=255, alpha=alpha),
  #                 rgb(2,129,138,max=255, alpha=alpha),
  #                 rgb(1,100,80,max=255, alpha=alpha))
  #yellow green blue
  levels(ii) <- c(rgb(5, 113, 176,max=255, alpha=alpha),
                  rgb(146, 197, 222,max=255, alpha=alpha),
                  rgb(247, 247, 247,max=255, alpha=alpha),
                  rgb(244, 165, 130,max=255, alpha=alpha),
                  rgb(202, 0, 32,max=255, alpha=alpha))
  
  
  ii = as.character(ii)
  ii[is.na(ii)==TRUE] <- rgb(255,255,255,max=255)
  return(ii)
}


pdf("Figures/Fig2/Fig2_mapsviolins.pdf", width=6, height=6)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), heights=c(1.5,1))
par(oma=c(0,0,0,0), mar=c(0,0,0,0), xpd=NA)


map(database = "state", regions=c("Wisconsin", "Iowa", "Illinois"), fill = TRUE, col="white", fg="grey30", lwd=1)
points(mapdat.mw$nhd_long, mapdat.mw$nhd_lat, col= c("black"), bg=get.col.bins(mapdat.mw$tpmed), pch = c(21), cex=1.3)
text(-98, 47.3, "a) Midwestern Lakes", cex=1.5, pos=4)
text(-84, 47.3, "b) Northeastern Lakes", cex=1.5, pos=4)
text(-98, 36.5, "c) Covariates", cex=1.5, pos=4)

legend(-78, 35, legend=c("Midwest", "Northeast"), fill=c("grey80", "grey20"), bty='n')

map(database = "state", regions=c("New York"), fill = TRUE, col="white", fg="grey30", lwd=1)
points(mapdat.ny$nhd_long, mapdat.ny$nhd_lat, col= c("black"), bg=get.col.bins(mapdat.mw$tpmed), pch = c(21), cex=1.3)

points(x = seq(from = -81.5, to =-75.5, by = (1.5)), y= rep(40,5), pch = 22, cex = 5, bg = colors, col="grey30")
text(x = seq(from = -81.5, to = -75.5, by = (1.5)), y = rep(39.3,5), labels=c("13", "23", "60", "110", ">110"), cex = 1.2)
text(-80.2, 40.8, "TP (ug/L)", cex=1.3, pos=4)

## violin plot - combine midwest and NE


dev.off()

#make boxplot to add after. ugh.
baseflow=covars.long[covars.long$covariate=="hu8_baseflow",]
ndep=covars.long[covars.long$covariate=="hu8_totalndepo",]
urban=covars.long[covars.long$covariate=="ws_urban",]
forest=covars.long[covars.long$covariate=="ws_forest",]

pdf("Figures/Fig2/Fig2_boxplotinset.pdf", width=4, height=2.5)
par(mfrow=c(1,4),oma=c(0,0,0,0), mar=c(.5,2,2.5,.5), xpd=NA)

boxplot(value~region, data=baseflow, notch=TRUE,
        col=(c("gray80","gray20")), xaxt='n', xlab="", ylab="", cex.axis=1, frame.plot=FALSE,
        ylim=c(0,100))
segments(.4, -4, .4, 110, lwd=2)
segments(.4, 110, 12.3, 110, lwd=2)
segments(.4, -4, 12.3, -4, lwd=2)
segments(3.65, -4, 3.65, 110, lwd=2)
segments(6.85, -4, 6.85, 110, lwd=2)
segments(10.1, -4, 10.1, 110, lwd=2)
segments(12.3, -4, 12.3, 110, lwd=2)
text(2.025, 111, "Baseflow (%)", pos=3)
text(5.25, 111, "N Deposition (kg/ha)", pos=3)
text(8.475, 111, "% Urban", pos=3)
text(11.2, 111, "% Forest", pos=3)
boxplot(value~region, data=ndep, notch=TRUE,
        col=(c("gray80","gray20")), xaxt='n', xlab="", ylab="", cex.axis=1, frame.plot=FALSE,
        ylim=c(0,10))

boxplot(value~region, data=urban, notch=TRUE,
        col=(c("gray80","gray20")), xaxt='n', xlab="", ylab="", cex.axis=1, frame.plot=FALSE,
        ylim=c(0,100))
boxplot(value~region, data=forest, notch=TRUE,
        col=(c("gray80","gray20")), xaxt='n', xlab="", ylab="", cex.axis=1, frame.plot=FALSE,
        ylim=c(0,100))
#segments()

dev.off()




