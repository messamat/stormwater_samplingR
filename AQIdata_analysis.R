library(rprojroot)
library(data.table)
library(ggplot2)
library(rgdal)
library(magrittr)
library(mgcv)
library(DescTools)
library(MARSS)
library(tsibble)
library(tidyr)
library(tseries)
library(forecast)

rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
datadir <- file.path(rootdir, "data")
inspectdir <- file.path(resdir, "data_inspection")
moddir <- file.path(resdir, "data_modeling")
#load(file.path(rootdir, 'src/stormwater_samplingR/20190225dat.Rdata'))

# 2. Analyze annual summary data
############################################################################################################################################
#Import data
AQIsites <- as.data.table(readOGR(dsn = file.path(resdir, 'airdata'), layer = 'airsites20190315'))
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

length(unique(AQIZn$UID))
unique(AQIZn$`Parameter Name`)

check <- AQIZn[`Parameter Name`=='Zinc PM2.5 LC',]
check <- AQIZn[`Parameter Name`=='Zinc (TSP) STP',]

colnames(AQIZn)
qplot(AQIZn$totalcount) + scale_x_log10()
Zngam <- gam(mean15_18~ s(predZn, k=3), 
             data=AQIZn[`Parameter Name`=='Zinc (TSP) STP' & totalcount>=50,])
summary(Zngam)
sefit <- predict(Zngam, se.fit = T)
MAE(AQIZn[`Parameter Name`=='Zinc (TSP) STP',mean15_18], fitted(Zngam))

AQIZn_pred <- ggplot(AQIZn[`Parameter Name`=='Zinc (TSP) STP' & totalcount>=50,], 
       aes(x=predZn, y=mean15_18)) + 
  # geom_point(data = AQIZn[`Parameter Name`=='Zinc PM2.5 LC' & `Observation Count`>=48,],
  #            color='grey') +
  geom_point(aes(color=Location.S)) +
  #scale_color_distiller(palette='Spectral') +
  geom_line(aes(y=sefit$fit), size=1.3, color='red') +
  geom_ribbon(aes(ymin=sefit$fit-1.96*sefit$se.fit, ymax = sefit$fit+1.96*sefit$se.fit), fill='orange', alpha=1/4) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x='Predicted Zinc index', y='Zinc PM2.5 (ug/m3)') +
  #geom_quantile(quantiles = c(0.1,0.5,0.9)) +
  #geom_smooth(method='gam') +
  theme_classic() +
  theme(text=element_text(size=20))
#png(file.path(resdir, 'airdata', 'Znmodel_scatterplot.png'), width=6, height=6, units='in', res=300) 
AQIZn_pred
#dev.off()

AQIZn[`Parameter Name`=='Zinc PM2.5 LC' & totalcount>=350, 
      SAPE := 100*(sefit$fit-mean15_18)/((sefit$fit+mean15_18)/2)]
writeOGR(SpatialPointsDataFrame(coords=coordinates(AQIZn[,.( Longitude, Latitude)]),
                                data= AQIZn, proj4string=CRS('+init=epsg:4326')),
         dsn = file.path(resdir, 'airdata'),
         layer = 'predAQIZn',
         driver = 'ESRI Shapefile')


# 2. Analyze daily summary data
############################################################################################################################################
#Import data
AQIsites <- as.data.table(readOGR(dsn = file.path(resdir, 'airdata'), layer = 'airsites'))
AQIfiles <- list.files(file.path(datadir, 'EPA_AirData_201902'), full.names = T, recursive = T)

PM10 <- as.data.table(do.call(rbind, lapply(grep('.*10SPEC.*[.]csv', AQIfiles, value=T), fread)))
SPEC<- as.data.table(do.call(rbind, lapply(grep('.*_SPEC.*[.]csv', AQIfiles, value=T), fread)))
temp <- as.data.table(do.call(rbind, lapply(grep('.*_TEMP.*[.]csv', AQIfiles, value=T), fread)))
wind <- as.data.table(do.call(rbind, lapply(grep('.*_WIND.*[.]csv', AQIfiles, value=T), fread)))
rh <- as.data.table(do.call(rbind, lapply(grep('.*_RH_DP.*[.]csv', AQIfiles, value=T), fread)))
press <- as.data.table(do.call(rbind, lapply(grep('.*_PRESS.*[.]csv', AQIfiles, value=T), fread)))

#Format data
SPEC[, `:=`(UID = paste0(`State Code`, `County Code`, `Site Num`),
            UIDmeth = paste0(`State Code`, `County Code`, `Site Num`, `Method Code`),
            date = as.Date(`Date Local`),
            specmea = `Arithmetic Mean`)] #rename var for sake of simplificity
#table(SPEC$`Parameter Name`)

temp[, `:=`(UID = paste0(`State Code`, `County Code`, `Site Num`),
            date = as.Date(`Date Local`))][
              ,tempmea := mean(`Arithmetic Mean`), by=.(UID,date)]
#table(temp$`Parameter Name`)

wind[, `:=`(UID = paste0(`State Code`, `County Code`, `Site Num`),
            date = as.Date(`Date Local`),
            windmea = `Arithmetic Mean`)] 
#table(wind$`Parameter Name`)

rh[, `:=`(UID = paste0(`State Code`, `County Code`, `Site Num`),
            date = as.Date(`Date Local`),
            rhmea = `Arithmetic Mean`)] 
#table(rh$`Parameter Name`)

press[, `:=`(UID = paste0(`State Code`, `County Code`, `Site Num`),
            date = as.Date(`Date Local`),
            pressmea = `Arithmetic Mean`)] 
#table(press$`Parameter Name`)

#Inspect data
summary(SPEC)
summary(temp)

AQIZn <- SPEC[`Parameter Name`== 'Zinc PM2.5 LC',] #Isolate records of of Zinc PM2.5LC
AQIZn[, length(unique(UID))] #Number of stations
AQIZn[, yearN := .N, by=UIDmeth] #Compute # of records for each station
qplot(AQIZn$yearN) # Histogram of # of records

#Average across methods and remove UID-Date duplicates (because of multiple methods)
AQIZn[, methmean := mean(specmea), by=.(UID, date)][
  , methdiff := 100*(max(specmea)-min(specmea))/methmean, by=.(UID, date)]
     
ggplot(AQIZn, aes(methmean, methdiff)) +
         geom_point() +
         scale_y_sqrt() +
         scale_x_sqrt()

###############################################################################################
# Try fitting a model for one station
###############################################################################################
#------------- Format data -------------####
Znattri <- merge(AQIZn, rh[`Parameter Name`=='Relative Humidity ', .(UID, date, rhmea)], 
                 on=c('UID', 'date'), all.x=T, all.y=T) %>%
  merge(wind[`Parameter Name`=='Wind Speed - Resultant', .(UID, date, windmea)], 
        on=c('UID', 'date'), all.x=T, all.y=T) %>%
  merge(temp[!duplicated(temp[, .(UID,date)]),
             .(UID, date, tempmea)], on=c('UID', 'date'), all.x=T, all.y=T)

teststation <- Znattri[UID == 533330,]%>%
    complete(date = seq.Date(min(date), max(date), by='days'), 
           fill = list(value = NA)) %>%
  setDT

#Scale variables
meacols <- c("specmea", "tempmea", "windmea", "rhmea")
teststation[, (meacols) := lapply(.SD, scale), .SDcols=meacols]

#Assess % of missing data  
teststation[, lapply(.SD, function(x) round(100*sum(is.na(x))/.N, 2)), .SDcols=meacols]

#------------- Check weekly trend -------------####
weekday_levels = c('Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday')
teststation[, `:=`(weekday = factor(weekdays(date),levels= weekday_levels),
                   week = week(date))][
                     , week_standardized := (specmea-mean(specmea, na.rm=T))/mean(specmea, na.rm=T), by=week][
                       , weekday_mean := mean(week_standardized, na.rm=T), by=weekday]


ggplot(teststation, aes(x=weekday, y=week_standardized)) +
  geom_hline(yintercept=0, color='grey', size=1.3) +
  geom_point(alpha=1/2) +
  geom_line(aes(y=weekday_mean, group=1), color='red') +
  theme_classic()

#------------- Plot data -------------####
ggplot(teststation, aes(x=date, y=specmea, group=UID)) +
  #geom_point() +
  geom_line(aes(y=tempmea), color='red') +
  geom_line(aes(y=windmea), color = 'darkgreen') +
  geom_line(aes(y=rhmea), color='blue') +
  geom_point(size=1.3) +
  theme_classic()

ggplot(teststation, aes(x=date, y=specmea, group=UID)) +
  #geom_point() +
  geom_point(size=1.3) +
  geom_point(aes(y=tempmea), color='red') +
  geom_point(aes(y=windmea), color = 'darkgreen') +
  geom_point(aes(y=rhmea), color='blue') +
  theme_classic()
  
ggplot(teststation, aes(x=tempmea, y=specmea)) +
  geom_point()
ggplot(teststation, aes(x=windmea, y=specmea)) +
  geom_point()

#------------- Make union of time series -------------####
tsspec <- teststation[min(which(!is.na(specmea))):.N,ts(specmea, frequency=1, start=min(as.Date(date)))]
tsspec7 <- teststation[min(which(!is.na(specmea))):.N,ts(specmea, frequency=7, start=min(as.Date(date)))]
tstemp <- teststation[, ts(tempmea, frequency=1, start=min(as.Date(date)))]
tswind <- teststation[, ts(windmea, frequency=1, start=min(as.Date(date)))]
tsrh <- teststation[, ts(rhmea, frequency=1, start=min(as.Date(date)))]
tsformat <- ts.union(tsspec, tstemp, tswind, tsrh)
dim(tsformat)

#------------- Check auto.arima -------------####
acf(tsspec, na.action = na.pass)
auto.arima(tsspec) 
auto.arima(tsspec7) #Too many missing values

#------------- MARSS analysis without covariates -------------####
mod.list <- list(
  B=matrix('b'), U=matrix('u'), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("mu"), tinitx=0 )

spec_format <- as.matrix(dcast(teststation[min(which(!is.na(specmea))):.N, .(date, specmea)], 
                               .~date, value.var = 'specmea')[,-1])
fit <- MARSS(100*spec_format, model=mod.list)

ggplot(teststation[min(which(!is.na(specmea))):.N,], aes(x=date, y=specmea)) +
  geom_point() +
  geom_line(aes(y=fit$states[1,]), color='red', size=1, alpha=1/3)

par(mfrow=c(1,2))
resids <- residuals(fit)
plot(resids$model.residuals[1,], 
     ylab="model residual", xlab="", main="flat level")
abline(h=0)
plot(resids$state.residuals[1,], 
     ylab="state residual", xlab="", main="flat level")
abline(h=0)
acf(resids$model.residuals[1,], main="flat level v(t)")

#------------- Format covariates -------------####
covar_format <- t(teststation[min(which(!is.na(specmea))):.N, 
                              lapply(.SD, function(x) scale(na.interp(x))), 
                              .SDcols = c('tempmea', 'windmea')]) #'rhmea'
row.names(covar_format) <- c('tempmea', 'windmea')

#------------- MARSS model with covariates of observation error, process-error and observation-error -------------####
mod.list_covobs <- list(
  B=matrix('b'), U=matrix('u'), C=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), D='unconstrained', R=matrix('r'),
  d=covar_format)
fit_covobs <- MARSS(100*spec_format, model=mod.list_covobs)

ggplot(teststation[min(which(!is.na(specmea))):.N,], aes(x=date, y=specmea)) +
  geom_point() +
  geom_line(aes(y=fit_covobs$states[1,]), color='red', size=1, alpha=1/3) 

covobs_resids <- residuals(fit_covobs)
plot.ts(covobs_resids$model.residuals[1,], na.action=na.omit)
acf(covobs_resids$model.residuals[1,], na.action=na.omit)
plot.ts(covobs_resids$state.residuals[1,], na.action=na.omit)
acf(covobs_resids$state.residuals[1,], na.action=na.omit)

#------------- MARSS model with covariates of states, process-error only model -------------####
mod.list_covproc <- list(
  B=matrix('b'), U=matrix('u'), C='unconstrained', Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  c=covar_format)

fit_covproc <- MARSS(100*spec_format, model=mod.list_covproc)

covproc_resids <- residuals(fit_covproc)$model.residuals[1,]
plot.ts(covproc_resids)
acf(covproc_resids)

#Compare models
c(fit$AICc, fit_covobs$AICc, fit_covproc$AICc)

#Check residuals
auto.arima(ts(covproc_resids, frequency=7))
auto.arima(covobs_resids)
auto.arima(ts(covobs_resids, frequency=7))

#Check states without covariate effect
pred_nocovar <- fit_covproc$model$data[1,] -
  (fit_covproc$coef['C.(X.Y1,tempmea)']*covar_format[1,] +
     fit_covproc$coef['C.(X.Y1,windmea)']*covar_format[2,])
pred_nocovar <- data.frame(pred=pred_nocovar, teststation[min(which(!is.na(specmea))):.N, .(date, UID)])

ggplot(teststation[!is.na(specmea),], aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(group=1), color='black') +
  geom_line(data=pred_nocovar[!is.na(pred_nocovar$pred),],
            aes(x=date, y=pred, group=UID), color='red', size=1, alpha=1/2)

#------------- Add day of the week as fixed effect ---------------------####
#https://nwfsc-timeseries.github.io/atsa-labs/sec-msscov-season.html
# number of days of the week
period <- 7
TT <- dim(spec_format)[2]
# first week day (Sunday:1, Monday:2)
per.1st <- teststation[min(which(!is.na(specmea))), as.numeric(weekday)]
# create factors for days of the week
c.in <- diag(period)
for(i in 1:(ceiling(TT/period))) {c.in <- cbind(c.in,diag(period))}
# trim c.in to correct start & length
c.in <- c.in[,(1:TT)+(per.1st-1)]
# better row names
rownames(c.in) <- weekday_levels

#Add week days to covariates
covar_formatwd <- rbind(c.in, covar_format)
Cweekdays = c(weekday_levels, 'temp', 'wind') #'rh
  
  
#Model paramter list with week days
mod.list_covproc_wd <- list(
  B=matrix('b'), U=matrix(0), C=matrix(t(Cweekdays)), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), D=matrix(0), R=matrix(0),
  c=covar_formatwd)
fit_covproc_wd <- MARSS(spec_format, model=mod.list_covproc_wd)

ggplot(teststation[min(which(!is.na(specmea))):.N,], aes(x=date, y=specmea)) +
  geom_point() +
  geom_line(aes(y=fit_covproc_wd$states[1,]), color='red', size=1, alpha=1/3)

#Compare models
c(fit$AICc, fit_covobs$AICc, fit_covproc$AICc, fit_covproc_wd$AICc)


