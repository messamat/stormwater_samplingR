library(data.table)
library(ggplot2)

rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
datadir <- file.path(rootdir, "data")
inspectdir <- file.path(resdir, "data_inspection")
moddir <- file.path(resdir, "data_modeling")

#Import data
AQIsites <- as.data.table(readOGR(dsn = file.path(resdir, 'airdata'), layer = 'airsites'))
AQI2018<- fread(file.path(datadir, 'EPA_AirData_201902/annual_conc_by_monitor_2018/annual_conc_by_monitor_2018.csv'))
AQI2017<- fread(file.path(datadir, 'EPA_AirData_201902/annual_conc_by_monitor_2017/annual_conc_by_monitor_2017.csv'))
AQI2016<- fread(file.path(datadir, 'EPA_AirData_201902/annual_conc_by_monitor_2016/annual_conc_by_monitor_2016.csv'))
AQI2015<- fread(file.path(datadir, 'EPA_AirData_201902/annual_conc_by_monitor_2015/annual_conc_by_monitor_2015.csv'))

#Get scaling for pollution variables
scaling <- readRDS(file.path(moddir, 'fieldXRFmodels_scaling.rds'))

#Predict Zn for each site
AQIsites[, heatbing1902log300 := 10*hbglog300/scaling$heatbing1902log300]
AQIsites[, heatSPDlog300 := 100*hSPDlog300/scaling$heatSPDlog300]
qplot(pollutfieldclean_cast$heatSPDlog300)
qplot(AQIsites$heatSPDlog300)


AQIsites[, predZn := exp(predict(modlistlogZn[[42]], 
                                 data.frame(heatbing1902log300, heatSPDlog300)))]


#Summarize 2015-2018 air monitoring data
AQIdata <- do.call(rbind, list(AQI2018, AQI2017, AQI2016, AQI2015))
AQIdata[, UID := paste0(`State Code`, `County Code`, `Site Num`)]
AQIsummary <- AQIdata[, .(mean15_18 = sum(`Observation Count`*`Arithmetic Mean`)/sum(`Observation Count`),
                          totalcount = sum(`Observation Count`)),
                      by=c('UID', 'Parameter Name')]
AQIZn <- AQIsummary[grepl('.*Zinc.*', `Parameter Name`),] %>%
  .[AQIsites, on='UID', , nomatch=0]

unique(AQIZn$`Parameter Name`)

check <- AQIZn[`Parameter Name`=='Zinc PM2.5 LC',]
colnames(AQIZn)
qplot(AQIZn$totalcount) + scale_x_log10()
Zngam <- gam(mean15_18~ s(predZn, k=3), 
             data=AQIZn[`Parameter Name`=='Zinc PM2.5 LC' & totalcount>=350,])
summary(Zngam)
sefit <- predict(Zngam, se.fit = T)
MAE(AQIZn[`Parameter Name`=='Zinc PM2.5 LC' & totalcount>=350,mean15_18], fitted(Zngam))

AQIZn_pred <- ggplot(AQIZn[`Parameter Name`=='Zinc PM2.5 LC' & totalcount>=350,], 
       aes(x=predZn, y=mean15_18)) + 
  # geom_point(data = AQIZn[`Parameter Name`=='Zinc PM2.5 LC' & `Observation Count`>=48,],
  #            color='grey') +
  geom_point() +
  geom_line(aes(y=sefit$fit), size=1.3, color='red') +
  geom_ribbon(aes(ymin=sefit$fit-1.96*sefit$se.fit, ymax = sefit$fit+1.96*sefit$se.fit), fill='orange', alpha=1/4) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x='Predicted Zinc index', y='Zinc PM2.5 (ug/m3)') +
  #geom_quantile(quantiles = c(0.1,0.5,0.9)) +
  #geom_smooth(method='gam') +
  theme_classic() +
  theme(text=element_text(size=20))
png(file.path(resdir, 'airdata', 'Znmodel_scatterplot.png'), width=6, height=6, units='in', res=300) 
AQIZn_pred
dev.off()

SpatialPointsDataFrame(AQIZn, coordinates(AQIZn$Latitude, ))
writeOGR()





