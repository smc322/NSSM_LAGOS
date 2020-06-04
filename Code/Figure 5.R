
# PLots for Figure 5 showing the estimated effect of environmental covariates on 
# the residual spatial structure of lake TP and chla.

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