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
              , tempmea := mean(`Arithmetic Mean`), by=.(UID,date)] %>%
  .[order(UID, date), 
    `:=`(tempmea_lag1 = shift(tempmea, n=1, type='lag', fill=NA),
         tempmea_lag2 = shift(tempmea, n=2, type='lag', fill=NA)), 
    by=.(UID, POC, `Parameter Name`)]
#table(temp$`Parameter Name`)

wind[, `:=`(UID = paste0(`State Code`, `County Code`, `Site Num`),
            date = as.Date(`Date Local`),
            windmea = `Arithmetic Mean`),] %>%
  .[order(UID, date), 
    `:=`(windmea_lag1 = shift(windmea, n=1, type='lag', fill=NA),
         windmea_lag2 = shift(windmea, n=2, type='lag', fill=NA)), 
    by=.(UID, POC, `Parameter Name`)]
#table(wind$`Parameter Name`)

rh[, `:=`(UID = paste0(`State Code`, `County Code`, `Site Num`),
          date = as.Date(`Date Local`),
          rhmea = `Arithmetic Mean`)] %>%
  .[order(UID, date), 
    `:=`(rhmea_lag1 = shift(rhmea, n=1, type='lag', fill=NA),
         rhmea_lag2 = shift(rhmea, n=2, type='lag', fill=NA)), 
    by=.(UID, POC, `Parameter Name`)]
#table(rh$`Parameter Name`)

press[, `:=`(UID = paste0(`State Code`, `County Code`, `Site Num`),
             date = as.Date(`Date Local`),
             pressmea = `Arithmetic Mean`)] %>%
  .[order(UID, date), 
    `:=`(pressmea_lag1 = shift(pressmea, n=1, type='lag', fill=NA),
         pressmea_lag2 = shift(pressmea, n=2, type='lag', fill=NA)), 
    by=.(UID, POC, `Parameter Name`)]
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
Znattri <- merge(AQIZn, rh[`Parameter Name`=='Relative Humidity ', 
                           .(UID, date, rhmea)], 
                 on=c('UID', 'date'), all.x=T, all.y=T) %>%
  merge(wind[`Parameter Name`=='Wind Speed - Resultant', 
             .(UID, date, windmea, windmea_lag1, windmea_lag2)], 
        on=c('UID', 'date'), all.x=T, all.y=T) %>%
  merge(temp[!duplicated(temp[, .(UID,date)]),
             .(UID, date, tempmea, tempmea_lag1, tempmea_lag2)],
        on=c('UID', 'date'), all.x=T, all.y=T)

teststation <- Znattri[UID == 533330,]%>% #Select station 
  complete(date = seq.Date(min(date), max(date), by='days'), 
           fill = list(value = NA)) %>% #Fill implicitly missing dates with explicit NAs
  setDT %>% #Set to data.table
  .[min(which(!is.na(specmea))):.N, ] #Only keep data after first recording of Zn

#Scale variables
meacols <- c("specmea", "tempmea", "tempmea_lag1", "tempmea_lag2",
             "windmea", "windmea_lag1", "windmea_lag2") #"rhmea", 'rhmea_lag1', "rhmea_lag2")
#teststation[, (meacols) := lapply(.SD, scale), .SDcols=meacols]

#Assess % of missing data  
teststation[, lapply(.SD, function(x) round(100*sum(is.na(x))/.N, 2)), .SDcols=meacols]

#------------- Check weekly trend -------------####
weekday_levels = c('Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday')
teststation[, `:=`(weekday = factor(weekdays(date),levels= weekday_levels),
                   week = week(date))][
                     !is.na(specmea),
                     `:=`(week_standardized = (specmea-mean(specmea, na.rm=T))/mean(specmea, na.rm=T),
                          weekN = .N), by=week][
                            weekN > 3
                            , weekday_mean := mean(week_standardized, na.rm=T), by=weekday]


ggplot(teststation[weekN>3,], aes(x=weekday, y=week_standardized)) +
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
tsspec <- teststation[,ts(specmea, frequency=1, start=min(as.Date(date)))]
tsspec7 <- teststation[,ts(specmea, frequency=7, start=min(as.Date(date)))]
tstemp <- teststation[, ts(tempmea, frequency=1, start=min(as.Date(date)))]
tswind <- teststation[, ts(windmea, frequency=1, start=min(as.Date(date)))]
tsrh <- teststation[, ts(rhmea, frequency=1, start=min(as.Date(date)))]
tsformat <- ts.union(tsspec, tstemp, tswind, tsrh)
tscovar <- ts.union(tstemp, tswind)
dim(tsformat)

#------------- Check auto.arima -------------####
acf(tsspec, na.action = na.pass)
auto.arima(tsspec) 
auto.arima(tsspec7) #Too many missing values

#------------- Try fitting ARIMAX ---------####
#Wind and temp
fit_arimaxwindtemp <- auto.arima(100*tsspec, xreg=teststation[, .(windmea, tempmea)]) 
fit_arimaxwindtemp
checkresiduals(fit_arimaxwind)

#Wind only
fit_arimaxwind <- auto.arima(100*tsspec, xreg=teststation[, .(windmea)]) 
fit_arimaxwind
checkresiduals(fit_arimaxwind)

#wind only + wind lag1
fit_arimaxwind <- auto.arima(100*tsspec, xreg=teststation[, .(windmea, windmea_lag1)]) 
fit_arimaxwind
checkresiduals(fit_arimaxwind)

fc <- forecast(fit_arimaxwind, xreg=teststation[, .(windmea, windmea_lag1)])$mean
ggplot(teststation, aes(x=date)) + 
  geom_point(aes(y= 100*specmea)) +
  geom_line(aes(date, y=fc, group=1), color='red')

#wind only + wind lag1 + wind lag2
fit_arimaxwind <- auto.arima(100*tsspec, xreg=teststation[, .(windmea, windmea_lag1, windmea_lag2)]) 
fit_arimaxwind
checkresiduals(fit_arimaxwind)

fc <- forecast(fit_arimaxwind, xreg=teststation[, .(windmea, windmea_lag1, windmea_lag2)])
ggplot(teststation, aes(x=date)) + 
  geom_point(aes(y= 100*specmea)) +
  geom_ribbon(aes(date, ymin=fc$lower[,'95%'], ymax=fc$upper[,'95%']), fill='orange', alpha=1/3) +
  geom_line(aes(date, y=fc$mean, group=1), color='red') 

#temp only
fit_arimaxtemp <- auto.arima(100*tsspec, xreg=teststation[, .(tempmea)]) 
fit_arimaxtemp

#temp only + temp lag1
fit_arimaxtemp <- auto.arima(100*tsspec, xreg=teststation[, .(tempmea, tempmea_lag1)]) 
fit_arimaxtemp

#------------- 1. MARSS model without covariates -------------####
mod.list <- list(
  B=matrix(1), U=matrix('u'), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("mu"), tinitx=0 )

spec_format <- as.matrix(dcast(teststation[, .(date, specmea)], 
                               .~date, value.var = 'specmea')[,-1])
fit <- MARSS(100*spec_format, model=mod.list)

ggplot(teststation[,], aes(x=date, y=100*specmea)) +
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
MARSScovar <- c('tempmea','tempmea_lag1', 'tempmea_lag2',
                'windmea', 'windmea_lag1', 'windmea_lag2')
covar_format <- t(teststation[, 
                              lapply(.SD, function(x) scale(na.interp(x))), 
                              .SDcols = MARSScovar]) #'rhmea'
row.names(covar_format) <- MARSScovar

#Format day of the week fourier series
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
covar_formatwd <- rbind(c.in, covar_format[c('tempmea', 'windmea'),])
Cweekdays = c(weekday_levels, 'tempmea', 'windmea') #'rh

#------------- 2a. MARSS ARMA model with covariates of observation error, process-error and observation-error -------------####
mod.list_covobs <- list(
  B=matrix('b'), U=matrix(0), C=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), D='unconstrained', R=matrix('r'),
  d=covar_format[c('tempmea', 'windmea'),])
fit_covobs <- MARSS(100*spec_format, model=mod.list_covobs, control=list(maxit=500))

#States (The expected value of x conditioned on the data)
ggplot(teststation, aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=fit_covobs$states[1,]), color='red', size=1, alpha=1/3) +
  geom_ribbon(aes(ymin=fit_covobs$states[1,] - 1.96*fit_covobs$states.se[1,],
                  ymax=fit_covobs$states[1,] + 1.96*fit_covobs$states.se[1,]), fill='orange', alpha=1/3)

#The expected value of y conditioned on the data. (just y for those y that are not missing)
covobs_p_ytt <- ggplot(teststation[,], aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=fit_covobs$ytT[1,]), color='red', size=1, alpha=1/3) 
covobs_p_ytt

#Simulated y (mean of 50 simulations)
sim50 <- MARSSsimulate(fit_covobs, nsim=50)$sim.data
sim50m <- apply(sim50[1,,], 1, mean)
covobs_p_sim <- ggplot(teststation, aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=sim50m), color='red', size=1, alpha=1/3) 
covobs_p_sim

#Kalman filter predictions of states
kf_covobs <- print(fit_covobs, what="kfs")
covobs_p_kf <- ggplot(teststation, aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=kf_covobs$xtt1[1,]), color='red', size=1) 
covobs_p_kf

#Check residuals
covobs_resids <- residuals(fit_covobs)
plot(covobs_resids$model.residuals[1,], na.action=na.omit)
acf(covobs_resids$model.residuals[1,], na.action=na.omit)
plot.ts(covobs_resids$state.residuals[1,], na.action=na.omit)
acf(covobs_resids$state.residuals[1,], na.action=na.omit)

#------------- 3a. MARSS ARMA model with covariates of states, process-error only model -------------####
mod.list_covproc <- list(
  B=matrix('b'), U=matrix(0), C='unconstrained', Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  c=covar_format[c('tempmea', 'windmea'),])
fit_covproc <- MARSS(100*spec_format, model=mod.list_covproc)

covproc_resids <- residuals(fit_covproc)$model.residuals[1,]
plot.ts(covproc_resids)
acf(covproc_resids)

#Compare models
c(fit$AICc, fit_covobs$AICc, fit_covproc$AICc)

#Check residuals
auto.arima(ts(covproc_resids, frequency=7))
auto.arima(covobs_resids$state.residuals[1,])
auto.arima(covobs_resids$model.residuals[1,])
auto.arima(ts(covobs_resids$state.residuals[1,], frequency=7))

#Kalman filter predictions of states
kf_covproc <- print(fit_covproc, what="kfs")
covproc_p_kf <- ggplot(teststation, aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=kf_covproc$xtt1[1,]), color='red', size=1) 
covproc_p_kf

#States
ggplot(teststation[,], aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=fit_covproc$states[1,]), color='red', size=1, alpha=1/3)

#States without covariate effect
pred_nocovar <- fit_covproc$model$data[1,] -
  (fit_covproc$coef['C.(X.Y1,tempmea)']*covar_format[1,] +
     fit_covproc$coef['C.(X.Y1,windmea)']*covar_format[2,])
pred_nocovar <- data.frame(pred=pred_nocovar, teststation[, .(date, UID)])

ggplot(teststation[!is.na(specmea),], aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(group=1), color='black') +
  geom_line(data=pred_nocovar[!is.na(pred_nocovar$pred),],
            aes(x=date, y=pred, group=UID), color='red', size=1, alpha=1/2)

#------------- Add day of the week as fixed effect to process but keep covariates as observation ---------------#
#Model paramter list with week days (will not converge with B=1 and U=matrix('u'))
mod.list_covproc_wd <- list(
  B=matrix('u'), U=matrix(0), C=matrix(t(Cweekdays)), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), D=matrix(0), R=matrix(0),
  c=covar_formatwd, d=covar_format[c('tempmea', 'windmea'),])
fit_covproc_wd <- MARSS(100*spec_format, model=mod.list_covproc_wd)

#States
ggplot(teststation[,], aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=fit_covproc_wd$states[1,]), color='red', size=1, alpha=1/3)

#Kalman filter predictions of states
kf_covproc_wd <- print(fit_covproc_wd, what="kfs")
covproc_wd_p_kf <- ggplot(teststation, aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=kf_covproc_wd$xtt1[1,]), color='red', size=1) 
covproc_wd_p_kf

#Check residuals
covproc_wd_resids <- residuals(fit_covproc_wd)$model.residuals[1,]
plot.ts(covproc_wd_resids)
acf(covproc_wd_resids)

#Compare models
c(fit$AICc, fit_covobs$AICc, fit_covproc$AICc, fit_covproc_wd$AICc)


#------------- 2b. MARSS model with observation covariates- add day of the week as fixed effect to process ---------------------####
#Model paramter list with week days (will not converge with B=1 and U=matrix('u'))
mod.list_covobs_wd <- list(
  B=matrix('u'), U=matrix(0), C=matrix(t(weekday_levels)), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), D='unconstrained', R=matrix('r'),
  c=c.in, d=covar_format[c('windmea', 'tempmea'),])
fit_covobs_wd <- MARSS(100*spec_format, model=mod.list_covobs_wd)

#States
ggplot(teststation[,], aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=fit_covobs_wd$states[1,]), color='red', size=1, alpha=1/3)

#Kalman filter predictions of states
kf_covobs_wd <- print(fit_covobs_wd, what="kfs")
covobs_wd_p_kf <- ggplot(teststation, aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=kf_covobs_wd$xtt1[1,]), color='red', size=1) 
covobs_wd_p_kf

#Check residuals
covobs_wd_resids <- residuals(fit_covobs_wd)$model.residuals[1,]
plot.ts(covobs_wd_resids)
acf(covobs_wd_resids)

#------------- 3b. MARSS model with state covariates - add day of the week as fixed effect to process ---------------------####
#Model paramter list with week days (will not converge with B=1 and U=matrix('u'))
mod.list_covproc_wd <- list(
  B=matrix('u'), U=matrix(0), C=matrix(t(Cweekdays)), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), D=matrix(0), R=matrix(0),
  c=covar_formatwd)
fit_covproc_wd <- MARSS(100*spec_format, model=mod.list_covproc_wd)

#States
ggplot(teststation[,], aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=fit_covproc_wd$states[1,]), color='red', size=1, alpha=1/3)

#Kalman filter predictions of states
kf_covproc_wd <- print(fit_covproc_wd, what="kfs")
covproc_wd_p_kf <- ggplot(teststation, aes(x=date, y=100*specmea)) +
  geom_point() +
  geom_line(aes(y=kf_covproc_wd$xtt1[1,]), color='red', size=1) 
covproc_wd_p_kf

#Check residuals
covproc_wd_resids <- residuals(fit_covproc_wd)$model.residuals[1,]
plot.ts(covproc_wd_resids)
acf(covproc_wd_resids)

#Compare models
c(fit$AICc, fit_covobs$AICc, fit_covproc$AICc, fit_covobs_wd$AICc, fit_covproc_wd$AICc)

#------------- Add lagged wind and temperature by one day  ---------------------####

#------------- Add # of antecedent days without winds  ---------------------####