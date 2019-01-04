

library(ggplot2)
library(plotly)
library(data.table)
library(rgdal)
library(sf)

rootdir <- "C:/Mathis/ICSL/stormwater"
datadir <- file.path(rootdir, 'data')
setwd(file.path(rootdir,"results"))
OSMSeattle <- as.data.table(readOGR(dsn = 'PSOSM.gdb', layer = 'OSMSeattle_datajoin'))
OSMPierce <- as.data.table(readOGR(dsn = 'PSOSM.gdb', layer = 'OSMPierce_datajoin'))
OSMKing <- as.data.table(readOGR(dsn = 'PSOSM.gdb', layer = 'OSMKing_datajoin_sel'))
OSMSDOT <- as.data.table(st_read(dsn = 'PSOSM.gdb', layer ='OSM_WSDOT_joinstats'))

streets <- as.data.table(readOGR(file.path(datadir, 'CitySeattle_20180601/Seattle_Streets/Seattle_Streets.shp')))
streets[, `:=`(OBJECTID = as.integer(as.character(OBJECTID)),
               SLOPE_PCT = as.numeric(as.character(SLOPE_PCT)))]
elv19range <- as.data.table(st_read(dsn = 'PSOSM.gdb', layer ='Seattle_elv19range'))
elv19range_smooth <- as.data.table(st_read(dsn = 'PSOSM.gdb', layer ='Seattle_elv19range_smooth'))
elv13range <- as.data.table(st_read(dsn = 'PSOSM.gdb', layer ='Seattle_elv13range'))
elv13range_smooth <- as.data.table(st_read(dsn = 'PSOSM.gdb', layer ='Seattle_elv13range_smooth'))

##############################################################################################################
# #Relate Open Street Map functional class to traffic volume for Pierce County and Seattle traffic data
##############################################################################################################

##### ------------------------ Analyze Seattle's own classification ---------------------------------####
#Order factors
streets$ARTDESCRIP <- factor(streets$ARTDESCRIP, 
                             levels=c('Not Designated', 'Collector Arterial','Minor Arterial', 'Principal Arterial', 'State Route/Freeway','Interstate/Freeway'))

#Arterial classification only slightly discriminate among speed limits
ggplot(streets[!(streets$ARTDESCRIP %in% c('County Arterial', NA)),], aes(x=SPEEDLIMIT, fill=ARTDESCRIP)) + 
  geom_histogram()+
  facet_wrap(~ARTDESCRIP, scales='free_y', ncol=1) +
  theme_bw()

ggplot(streets[!(streets$ARTDESCRIP %in% c('County Arterial', NA)),], aes(x=ARTDESCRIP, y=SPEEDLIMIT, fill=ARTDESCRIP)) + 
  geom_boxplot() +
  theme_classic()

#Arterial classification moderately discriminates among AADT levels
setDT(streets)[,mean(AADT_inter), by=ARTDESCRIP]

ggplot(streets[!(streets$ARTDESCRIP %in% c('County Arterial', NA)),], aes(x=ARTDESCRIP, y=AADT_inter, fill=ARTDESCRIP)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  theme_classic()

##### ------------------------ Analyze OSM functional classification (fclass field) ---------------------------------####
setnames(OSMSeattle, c('AADT_avg', 'SPEEDLIMIT'), c('ADT', 'SpeedLimit'))
setnames(OSMSDOT, c('MEAN_AADT', 'FIRST_fclass'), c('ADT', 'fclass'))
setnames(OSMKing, 'SPEED_LIM', 'SpeedLimit')

OSMSeattle[, agency := 'City of Seattle']
OSMPierce[, agency := 'Pierce County']
OSMPierce <- OSMPierce[order(-OSMPierce$ADTYear),]
OSMPierce_nodupli <- OSMPierce[!duplicated(OSMPierce$osm_id),]
OSMSDOT[, agency := 'WSDOT']
OSMKing[, agency := 'King']

agfields <- c('osm_id', 'fclass','ADT', 'agency', 'maxspeed', 'SpeedLimit')
OSMbind <- rbind(OSMPierce_nodupli[intersper>0.80 & LENGTH_GEO>10, agfields, with=F],
                 OSMSeattle[intersper>0.80 & LENGTH_GEO>10, agfields, with=F],
                 cbind(OSMSDOT[,agfields[agfields %in% colnames(OSMSDOT)], with=F], data.frame(maxspeed=NA, SpeedLimit=NA)),
                 cbind(OSMKing[,agfields[agfields %in% colnames(OSMKing)], with=F], data.frame(ADT=NA)))

OSMbind[, `:=`(N=.N, medianADT=median(ADT, na.rm=T), meanADT=mean(ADT, na.rm=T)), by=.(fclass)]
OSMbind$fclass <- factor(OSMbind$fclass, levels = unique(OSMbind$fclass[order(OSMbind$medianADT)]))
ggplot(OSMbind, aes(x=fclass, y=ADT, fill=agency)) + 
  geom_boxplot() + 
  geom_text(aes(label=N, y = 5)) +
  geom_text(aes(label=paste0('median: ',medianADT), y=medianADT)) +
  #geom_text(aes(label=paste0('mean: ',round(meanADT)), y=meanADT)) +
  scale_y_log10() + 
  theme_classic()
OSMbind[order(OSMbind$medianADT), as.integer(median(ADT, na.rm=T)), by=.(fclass)]


OSMbind[,SpeedLimit:= as.numeric(as.character(SpeedLimit))]
OSMbind[, `:=`(N=.N, medianSPD=median(SpeedLimit, na.rm=T), meanSPD=mean(SpeedLimit, na.rm=T)), by=.(fclass)]
OSMbind$fclass <- factor(OSMbind$fclass, levels = unique(OSMbind$fclass[order(OSMbind$medianSPD)]))
ggplot(OSMbind, aes(x=fclass, y=SpeedLimit, fill=agency)) + 
  geom_boxplot() + 
  geom_text(aes(label=N, y = 5)) +
  geom_text(aes(label=paste0('median: ',medianSPD), y=medianSPD)) +
  #geom_text(aes(label=paste0('mean: ',round(meanADT)), y=meanADT)) +
  scale_y_log10() + 
  theme_classic()

ggplot(OSMbind, aes(x=maxspeed, y=SpeedLimit)) + 
  geom_point()

##############################################################################################################
# #Check accuracy of estimated road gradient using NED 1/3 and 1/9 arc seconds compared to Seattle gradient data
##############################################################################################################
#Compare Seattle estimated road gradient to actual one
streets_elv <- streets[elv19range, on='OBJECTID==Value']
setnames(streets_elv, 'RANGE', 'RANGE19')
streets_elv <- streets_elv[elv19range_smooth, on='OBJECTID==Value']
setnames(streets_elv, 'RANGE', 'RANGE19smooth')
streets_elv <- streets_elv[elv13range, on='OBJECTID==Value']
setnames(streets_elv, 'RANGE', 'RANGE13')
streets_elv <- streets_elv[elv13range_smooth, on='OBJECTID==Value']
setnames(streets_elv, 'RANGE', 'RANGE13smooth')

streets_elv[, `:=`(gradient19 = RANGE19/SHAPE_Leng,
                   gradient19_smooth = RANGE19smooth/SHAPE_Leng,
                   gradient13 = RANGE13/SHAPE_Leng,
                   gradient13_smooth = RANGE13smooth/SHAPE_Leng)]

grad19 <- ggplot(streets_elv, aes(x=100*gradient19_smooth, y=SLOPE_PCT, label=SHAPE_Leng)) + 
  geom_jitter()+
  geom_smooth(method='lm') + 
  geom_abline(slope=1)
grad19

grad13 <- ggplot(streets_elv, aes(x=100*gradient13_smooth, y=SLOPE_PCT, label=SHAPE_Leng)) + 
  geom_jitter()+
  geom_smooth(method='lm') + 
  geom_abline(slope=1)
grad13

ggplot(streets_elv, aes(x=gradient13_smooth, y=gradient19_smooth, label=SHAPE_Leng)) + 
  geom_jitter()

summary(lm(SLOPE_PCT~gradient19, data=streets_elv))
summary(lm(SLOPE_PCT~gradient19_smooth, data=streets_elv))
summary(lm(SLOPE_PCT~gradient13, data=streets_elv))
summary(lm(SLOPE_PCT~gradient13_smooth, data=streets_elv))
summary(lm(gradient19_smooth~gradient13_smooth, data=streets_elv))