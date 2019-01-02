##want to make a dataset for Pavel & Marie to start testing out the model -- only spatial variation, not temporal. start here w/ Secchi b/c most data, and treat ts obs from 1990-present as replicates bc shouldn't change much. 

library(LAGOSNE)
library(dplyr)
lagos<-lagosne_load(version="1.087.1")

locus<-lagos$locus
sec<-lagos$secchi

data.sec<-sec[,c("lagoslakeid", "secchi", "sampleyear", "samplemonth")]

sec.90<-data.sec[data.sec$sampleyear>1989,]
sec.90.jas<-sec.90[sec.90$samplemonth == 7 | sec.90$samplemonth == 8 | sec.90$samplemonth == 9,]

#add in some predictors

hu8.chag<-lagos$hu8.chag
hu8.chag.rel<-hu8.chag[,c(1,4, 20, 68, 84, 92, 100)]

iws.lulc<-lagos$iws.lulc
iws.lulc.rel<-iws.lulc[,c(52, 54, 56, 58, 60, 64, 66, 68, 74, 76, 78, 80, 156)]
iws.lulc.rel$urban<-iws.lulc.rel$iws_nlcd2001_pct_21+iws.lulc.rel$iws_nlcd2001_pct_22+iws.lulc.rel$iws_nlcd2001_pct_23+iws.lulc.rel$iws_nlcd2001_pct_24
iws.lulc.rel$rowcrop<-iws.lulc.rel$iws_nlcd2001_pct_82
iws.lulc.rel$pasture<-iws.lulc.rel$iws_nlcd2001_pct_81
iws.lulc.rel$forest<-iws.lulc.rel$iws_nlcd2001_pct_41+iws.lulc.rel$iws_nlcd2001_pct_42+iws.lulc.rel$iws_nlcd2001_pct_43
iws.lulc.rel$wetland<-iws.lulc.rel$iws_nlcd2001_pct_90+iws.lulc.rel$iws_nlcd2001_pct_95
iws.lulc.rel$openh20<-iws.lulc.rel$iws_nlcd2001_pct_11

iws.lulc.rel.keep<-iws.lulc.rel[,c(13:19)]
iws.lulc.rel.keep.nona<-na.omit(iws.lulc.rel.keep)

area<-locus[,c(1,6)]
wsarea<-lagos$iws[,c(2,12)]
maxd<-lagos$lakes_limno[,c(1,6)]
connclass<-lagos$lakes.geo[,c(1,10, 31)]

lakewsa<-merge(area, wsarea, by="lagoslakeid", all.x=T, all.y=T)
areadepth<-merge(lakewsa, maxd, by="lagoslakeid", all.x=T, all.y=T)
areadepthconn<-merge(areadepth, connclass, by="lagoslakeid", all.x=T, all.y=T)
areadepthconn$la_wa_ratio<-areadepthconn$lake_area_ha/areadepthconn$iws_ha


#combine predictors to merge w n
hu.8id<-locus[,c(1,13)]
chag.id<-merge(hu8.chag.rel, hu.8id, by="hu8_zoneid", all.x=T, all.y=T)

chag.lulc<-merge(chag.id, iws.lulc.rel.keep.nona, by="lagoslakeid", all.x=T, all.y=T)

chag.lulc.lake<-merge(chag.lulc, areadepthconn, by="lagoslakeid", all.x=T, all.y=T)

#get rid of hu8id, lake area, watershed area
chag.lulc.lake$iws_ha=NULL
chag.lulc.lake$hu8_zoneid=NULL

names(chag.lulc.lake) <- c("lagoslakeid","hu8_baseflow",   "hu8_no3depo",    "hu8_totalndepo", "hu8_MAP", "hu8_MAT", "hu8_runoff", "ws_urban", "ws_rowcrop", "ws_pasture", "ws_forest", "ws_wetland", "ws_openh20", "lake_area", "maxdepth", "state_zoneid", "lakeconnection", "la_wa_ratio" )

data.preds<-merge(sec.90.jas, chag.lulc.lake, by="lagoslakeid", all.x=T, all.y=F)


#add latlong
ll.l<-lagos$locus[,c(1,4,5)]
data.preds.ll<-merge(data.preds, ll.l, by="lagoslakeid", all.x=T, all.y=F)

#add state name instead of statezoneid

states<-lagos$state[,c("state_name", "state_zoneid")]

data.sn<-merge(data.preds.ll, states, by="state_zoneid", all.x=T, all.y=F)
data.sn$state_zoneid=NULL

setwd("/Users/SarahiMac/Dropbox/Sarah_Work/Research/2018_A&S_SeedGrant_Pavel/Marie Postdoc Folder/Data/02Jan2019")
write.csv(data.sn, file="LagosSecchi_02Jan2019.csv")
