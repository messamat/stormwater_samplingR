library(rprojroot)
library(data.table)
library(ggplot2)
library(GGally)
library(rgdal)
library(magrittr)
library(mgcv)
library(DescTools)
library(MARSS)
library(tsibble)
library(tidyr)
library(tseries)
library(zoo)
library(forecast)
library(lme4)
library(pls)
library(feather)

rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
datadir <- file.path(rootdir, "data")
AQIdir <- file.path(rootdir, 'results/airdata')
inspectdir <- file.path(resdir, "data_inspection")
moddir <- file.path(resdir, "data_modeling/Znspeciation_model")

sarlmmods <- readRDS(file.path(resdir,'data_modeling/fieldXRFmodels.rds'))
scaling <- readRDS(file.path(resdir, 'data_modeling/fieldXRFmodels_scaling.rds'))

if (!(dir.exists(moddir))) {
  dir.create(moddir)
}

#load(file.path(rootdir, 'src/stormwater_samplingR/20190225dat.Rdata'))

#### 1. Analyze daily summary data with NARR meteorological data ####
############################################################################################################################################
#------ Import and format NARR and Zn sensor data -----
airNARR <- fread(file.path(resdir, 'airdat_NARRjoin.csv'), colClasses = c('integer', 'character', 'character', 'character', 'character', 'character', 'numeric',
                                                                          'numeric', 'character', 'character', 'character', 'character', 'Date',
                                                                          'character', 'character', 'integer', 'numeric', 'numeric', 'numeric', 'numeric', 'character',
                                                                          'character', 'character', 'character', 'character', 'character', 'character', 'character', 
                                                                          'character', 'Date', 'character', 'Date', 'numeric', 'numeric', rep('numeric', 37)
))

colnames(airNARR) <- sapply(colnames(airNARR), function(x) gsub('\\s', '', x))
setnames(airNARR, c("DateLocal", "ArithmeticMean"), c("datelocal", "specmea")) 
airNARR[, datelocal := as.Date(datelocal, format="%Y-%m-%d")]


#Subset to only keep Zn and Al observations
NARRcols <- colnames(airNARR)[35:ncol(airNARR)] #NARR column names

castform <- paste0(paste(c('UID', 'datelocal', 'date', 'MethodCode', 'Latitude', 'Longitude', NARRcols), collapse=" + "),
                   ' ~ ParameterName')
airNARRZnAl <-  airNARR[!(duplicated(airNARR, by= c('UID', 'datelocal', 'ParameterName'))) &
                          (ParameterName == 'Zinc PM2.5 LC' | ParameterName == 'Aluminum PM2.5 LC'),] %>%
  dcast(formula = castform, value.var = 'specmea') %>%
  .[!is.na('Zinc PM2.5 LC'),]
setnames(airNARRZnAl, c('Zinc PM2.5 LC', 'Aluminum PM2.5 LC'), c("specmeaZn", "specmeaAl")) 

remove(airNARR)

#Fill implicitly missing dates with explicit NAs
# airNARRZnAlfill <- airNARRZn[, complete(.SD, datelocal = seq.Date(min(datelocal), max(datelocal), by='days'), 
#                                       fill = list(value = NA)), by=UID]
airNARRZnAlfill <- airNARRZnAl

#Floor values to just below minimum detectable quantity by method
table(airNARRZnAlfill$MethodCode)
airNARRZnAlfill[specmeaZn <= 0 & MethodCode =='800', specmeaZn := 0.0019]
airNARRZnAlfill[specmeaZn <= 0 & MethodCode =='811', specmeaZn := 0.0014]
airNARRZnAlfill[specmeaZn <= 0 & MethodCode =='846', specmeaZn := 0.0005]
airNARRZnAlfill[specmeaAl <= 0 & MethodCode =='800', specmeaAl := 0.0019]
airNARRZnAlfill[specmeaAl <= 0 & MethodCode =='811', specmeaAl := 0.01]
airNARRZnAlfill[specmeaAl <= 0 & MethodCode =='846', specmeaAl := 0.004]


#Compute weekdays and week-standardized Zn
weekday_levels = c('Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday')
airNARRZnAlfill[, `:=`(weekday = factor(weekdays(datelocal),levels= weekday_levels),
                       week = week(datelocal))][
                         !is.na(specmeaZn),
                         `:=`(week_standardizedZn = (specmeaZn-mean(specmeaZn, na.rm=T))/mean(specmeaZn, na.rm=T),
                              week_standardizedAl = (specmeaAl-mean(specmeaAl, na.rm=T))/mean(specmeaAl, na.rm=T),
                              weekN = .N), by=.(week, year(datelocal), UID)][
                                weekN > 3, 
                                `:=`(weekday_meanZn = mean(week_standardizedZn, na.rm=T),
                                     weekday_meanAl = mean(week_standardizedAl, na.rm=T)), by=.(weekday, UID)]

#Compute date number on record (1 corresponding to 2014/01/01 and 2191 to 2019/12/31)
datelist <- seq(as.Date("2014/01/01"), as.Date("2019/12/31"), by = "day")
datedt <- data.table(datelocal=datelist, datenum = seq_along(datelist))
airNARRZnAlfill <- airNARRZnAlfill[datedt, on='datelocal']

#qplot(airNARRZnAlfill[, sum(is.na(specmea))/.N, by=UID]$V1) #Check % missing 
qplot(airNARRZnAlfill[, mean(weekN, na.rm=T), by=UID]$V1) #Check average number of values per week



qplot(airNARRZnAlfill$specmeaZn) + scale_x_log10()
qplot(airNARRZnAlfill$specmeaAl) + scale_x_log10()
qplot(airNARRZnAlfill[,specmeaZn/specmeaAl])
qplot(airNARRZnAlfill[,specmeaZn/specmeaAl]) + scale_x_log10()

airNARRZnAlfill[!is.na(specmeaZn), .N, by=weekday]
airNARRZnAlfill[!is.na(specmeaZn/specmeaAl), .N, by=weekday]

mean(airNARRZnAlfill[!is.na(specmeaZn), .N, by=UID]$N)

#------ Import, format Zn predictor variables and join -----
predvars <- as.data.table(sf::st_read(file.path(AQIdir,'airsites.shp'))) %>%
  .[!nlcd_imp_A == 127,]

#Compute predZn
sarlmmod_Zn <- sarlmmods$logZnstandmod$coefficients
predvars[, aadtlog100frt := (aadtlog100/scaling$heatsubAADTlog100)^(1/4)][,
                                                                          predZn := exp(sarlmmod_Zn['(Intercept)'] + 
                                                                                          sarlmmod_Zn['heatsubAADTlog100frt']*aadtlog100frt +
                                                                                          sarlmmod_Zn['nlcd_imp_ps_mean']*nlcd_imp_A)]

#Join with NARR and speciation data
airNARRZnAljoin <- merge(airNARRZnAlfill, predvars, by='UID', all.y=F)

#Fill smoke indices with 0
airNARRZnAljoin[is.na(smokeindex_1d), smokeindex_1d := 0]
airNARRZnAljoin[is.na(smokeindex_3d), smokeindex_3d := 0]
airNARRZnAljoin[is.na(smokeindex_6d), smokeindex_6d := 0]

#Remove all records with NA values in meteorological variables
nrow(airNARRZnAljoin[!complete.cases(airNARRZnAljoin[, c('specmeaZn', 'specmeaAl', NARRcols), with=F]),])/nrow(airNARRZnAljoin) #Check % of records that will be removed
airNARRZnAljoin[,lapply(.SD, function(x) length(which(is.na(x)))), .SDcols = c('specmeaZn', 'specmeaAl', NARRcols)] #Check % of missing records by column
airNARRZnAl_nona <- airNARRZnAljoin[complete.cases(airNARRZnAljoin[,c('specmeaZn', 'specmeaAl', NARRcols), with=F]),]

#Compute Zn/Al ratio
airNARRZnAl_nona[, specmearatio := specmeaZn/specmeaAl]

#Select only stations with >200 records
Nobs <- airNARRZnAl_nona[!is.na(specmeaZn), .N, by=.(UID)]
selectedstations <- Nobs[N>200,UID]

#Scale NARR variables (z-standardize)
airNARRZnAl_nona[, (NARRcols) := lapply(.SD, scale), .SDcols=NARRcols]

#------ GLM models to control for week day and measurement method -----
#Compute for each station the proportion of each week day
stationweekdays <- airNARRZnAl_nona[,{
  stationN = .N
  .SD[,.(weekdayfrac=.N/stationN),by=weekday]
},by=UID]

ggplot(stationweekdays, aes(x=weekday, y=weekdayfrac)) + 
  geom_boxplot() +
  geom_abline(slope=0, intercept=1/7, color='red', alpha=1/2) + 
  theme_classic()

ggplot(airNARRZnAl_nona[UID %in% selectedstations & week_standardizedZn > 0.01,], aes(x=weekday, y=week_standardizedZn)) + 
  geom_boxplot() +
  scale_y_log10()

#Check relationship between method and mean recorded Zn
stationmethods <- airNARRZnAl_nona[,{
  stationN = .N
  .SD[,.(methodfrac=.N/stationN),by=MethodCode]
},by=UID]

ggplot(stationmethods, aes(x=MethodCode, y=methodfrac)) + 
  geom_boxplot() +
  theme_classic()

methodmode<- setkey(airNARRZnAl_nona[, list(freq = .N, meanZn = mean(specmeaZn), meanAl = mean(specmeaAl)), by=list(UID, MethodCode)], 
                    UID, freq)[J(unique(UID)), mult="last"]

ggplot(methodmode, aes(x=MethodCode, y=meanZn)) + 
  geom_boxplot() +
  scale_y_log10()

ggplot(methodmode, aes(x=MethodCode, y=meanAl)) + 
  geom_boxplot() +
  scale_y_log10()

#------ PLSR for daily data of Zn + NARR variables at each station with more than 200 obs -----
#station <- "4813916083899" #for debugging
for (station in selectedstations) { #will change to apply at some point
  print(station)
  #For each station
  stationdat <- airNARRZnAl_nona[UID==station,]
  
  #Develop plsr format to predict residuals from week days
  pls_format <- data.frame(
    weekday = stationdat$weekday,
    specmea = stationdat$specmeaZn,
    logspecmea = log(stationdat$specmeaZn),
    NARRvars = I(as.matrix(stationdat[,NARRcols, with=F]))
  )
  plsr1 <- plsr(logspecmea ~ NARRvars, ncomp = 5, data = pls_format, validation = "LOO")
  png(file.path(moddir, paste0('CVplot_', station, '.png')))
  compsel <- selectNcomp(plsr1, nperm=1000, alpha=0.1, method = "randomization", plot = TRUE)
  title(paste0('UID: ', station))
  dev.off()
  
  plsR2 <- R2(plsr1, ncomp=compsel)$val[compsel]
  predarray <- predict(plsr1, ncomp=compsel)
  
  if (compsel > 0){
    airNARRZnAl_nona[UID==station, plspred := exp(predarray)]
    png(file.path(moddir, paste0('residplot_', station, '.png')))
    plsplot <- ggplot(airNARRZnAl_nona[UID==station,], aes(x=plspred, y=specmeaZn)) + 
      geom_point() +
      geom_abline() + 
      annotate("text", x=-Inf, y=Inf, label =paste0('R2:', round(plsR2,2)), hjust = -0.25, vjust = 1) +
      labs(title=paste0('UID: ', station)) + 
      scale_y_log10() + 
      theme_classic() + 
      coord_fixed()
    print(plsplot)
    dev.off()
    
    #Create dataframe where all meteorological variables are 0 (US-wide average for prediction)
    predat <- data.frame(matrix(data=0, nrow=nrow(stationdat), ncol=length(NARRcols)))
    colnames(predat) <- NARRcols
    predat$NARRvars <- I(as.matrix(predat[,NARRcols]))
    #Generate predicted Zn level given mean meteorology + residuals
    airNARRZnAl_nona[UID==station, Znspecpred := exp(predict(plsr1, ncomp=compsel, newdata = predat))][
      UID==station, Znspecpred_nometeo := plspred + as.vector(residuals(plsr1, comp=compsel)[1:nrow(stationdat), 1, compsel])]
  } else {
    airNARRZnAl_nona[UID==station, Znspecpred_nometeo := specmeaZn]
  }
}

#Develop model based on weekday and long-term trend
glm_nometeo0 <- glm(Znspecpred_nometeo ~ 1, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo0)

glm_nometeo1 <- glm(Znspecpred_nometeo ~ weekday, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo1)

glm_nometeo2 <- glm(Znspecpred_nometeo ~ weekday + datenum, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo2)

glm_nometeo3 <- glm(Znspecpred_nometeo ~ 0 + UID + weekday + datenum, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo3)

glm_nometeo4 <- glm(Znspecpred_nometeo ~ 0 + predZn + weekday + datenum, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo4)

#
glm3coefs <- coefficients(glm_nometeo3)[grep('UID.*', names(coefficients(glm_nometeo3)), perl=T)]
dtformat <- data.table(UID = gsub('UID', '', names(glm3coefs)), Znspec_mean = glm3coefs)
predZnformat <- airNARRZnAl_nona[UID %in% selectedstations, .(predZn = max(predZn),
                                                              specmea_nometeomean = mean(Znspecpred_nometeo),
                                                              specmea_nometeomedian = median(Znspecpred_nometeo),
                                                              specmea_nometeo70th = quantile(Znspecpred_nometeo, 0.7),
                                                              specmea_nometeo90th = quantile(Znspecpred_nometeo, 0.9),
                                                              specmeaZn_mean = mean(specmeaZn),
                                                              specmeaZn_median = median(specmeaZn),
                                                              nlcd_imp_A = mean(nlcd_imp_A),
                                                              aadtlog100frt = mean(aadtlog100frt),
                                                              N = .N), by=UID]
NARRcolsmean <- airNARRZnAl_nona[UID %in% selectedstations, lapply(.SD, mean(x, na.rm=T)), 
                                 .SDcols = NARRcols, by=UID]

dtformatmerge <- merge(merge(dtformat, predZnformat, by='UID'),
                       NARRcolsmean, by='UID')

ggplot(dtformatmerge, aes(x=predZn, y=Znspec_mean)) + 
  geom_point() + 
  geom_smooth()

ggplot(dtformatmerge, aes(x=predZn, y=specmea_nometeomean)) + 
  geom_point() + 
  geom_smooth() +
  scale_y_log10() + 
  scale_x_log10()

lm_nometeo1 <- lm(log(specmea_nometeomean)~log(predZn), data=dtformatmerge)
summary(lm_nometeo1)

ggplot(dtformatmerge, aes(x=predZn, y=specmea_nometeomedian)) + 
  geom_point() + 
  geom_smooth()+
  scale_y_log10()

ggplot(dtformatmerge, aes(x=predZn, y=specmea_nometeo70th)) + 
  geom_point() + 
  geom_smooth()

#------ PLSR for long-term Zn data + NARR variables of all stations with more than 200 obs together -----
#Remove effects of weekday and method
glm0 <- glm(specmeaZn ~ 1, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm0)

glm_weekday <- glm(specmeaZn ~ weekday*UID, data=airNARRZnAl_nona[UID %in% selectedstations,], family= gaussian(link='log'))
summary(glm_weekday)
airNARRZnAl_nonaMonday <- airNARRZnAl_nona[UID %in% selectedstations,][, weekday := "Monday"]
airNARRZnAl_nonaMonday[, specmeaZn_monday := exp(predict(glm_weekday, airNARRZnAl_nonaMonday)) + (specmeaZn-exp(fitted(glm_weekday)))]

dtformatmerge <- merge(dtformatmerge, airNARRZnAl_nonaMonday[, .(specmeaZn_monday_mean = mean(specmeaZn_monday),
                                                                 specmeaZn_monday_median = median(specmeaZn_monday)),
                                                             by=.(UID)],
                       by = 'UID')


# glm_weekdaymethod <- glm(log(specmeaZn) ~ weekday*UID + MethodCode, data=airNARRZnAl_nona[UID %in% selectedstations,])
# summary(glm_weekdaymethod)
# coefficients(glm_weekdaymethod)

#Check long-term relationships with NARR variables
for (var in NARRcols) {
  print(var)
  png(file.path(moddir, paste0('NARRZnspecplot_', var, '.png')))
  narrplot <- ggplot(dtformatmerge, aes_string(x=var, y='specmeaZn_mean')) + 
    geom_point() + 
    geom_smooth() +
    theme_classic()
  print(narrplot)
  dev.off()
}

#Develop plsr format
plsall_format <- data.frame(
  specmeaZn = dtformatmerge$specmeaZn_mean,
  logspecmeaZn = log(dtformatmerge$specmeaZn_mean),
  specmeaZnMonday = dtformatmerge$specmeaZn_monday_mean,
  predZn = dtformatmerge$predZn,
  NARRvars = I(as.matrix(dtformatmerge[,NARRcols, with=F])),
  NARRvarsZn = I(as.matrix(dtformatmerge[,c(NARRcols, 'predZn'), with=F]))
)

#With just NARR variables
plsrall_narr <- plsr(logspecmeaZn ~ NARRvars, ncomp = 10, data = plsall_format, validation = "LOO")
compsel <- selectNcomp(plsrall_narr, nperm=1000, alpha=0.05, method = "randomization", plot = TRUE)
plsR2NARR <- R2(plsrall_narr, ncomp=compsel)$val[2]
predarrayNARR <- predict(plsrall_narr, ncomp=compsel)

dtformatmerge[, `:=`(plspredNARR = exp(predarrayNARR),
                     plsresidNARR = exp(predarrayNARR)-specmeaZn_mean)]


png(file.path(moddir, paste0('residplot_allnarr.png')))
plsplotNARR <- ggplot(dtformatmerge, aes(1000*plspredNARR, 1000*specmeaZn_mean)) + 
  geom_abline() + 
  geom_point(alpha = 0.6, size=1.5) +
  scale_y_log10(name = expression(Observed~average~Zn~PM[2.5]~(mg/m^3)),limits = c(0.7, 120)) +
  scale_x_log10(name = expression(Predicted~average~Zn~PM[2.5]~(mg/m^3)), limits = c(0.7, 120)) +
  annotate("text", x=10, y=100, vjust=-1, #hjust=0.80, 
           label =paste0('A.~Meteorology-only~~R^2 ==', round(plsR2NARR,2)), parse=T, size=4) + 
  coord_fixed(clip='off') +
  theme_classic() + 
  theme(text = element_text(size=12),
        #axis.title.x = element_blank(), 
        plot.margin = margin(0.6, 0.25, 0, 0.25, "cm"),
        rect = element_rect(fill = "transparent"))
print(plsplotNARR)
dev.off()

ggplot(dtformatmerge, aes(x=predZn, y=plsresidNARR)) + 
  geom_point() + 
  geom_smooth(method='lm') 
  
gam_plsrallresid<- (mgcv::gam(plsresid~s(predZn, k=4), data=dtformatmerge))
summary(gam_plsrallresid)
plot(gam_plsrallresid, rug = TRUE, residuals = TRUE, pch = 1, cex = 1)

#With NARR and Zn
plsrall_narrZn <- plsr(logspecmeaZn ~ NARRvarsZn, ncomp = 10, data = plsall_format, validation = "LOO")
compsel <- selectNcomp(plsrall_narrZn, nperm=1000, alpha=0.05, method = "randomization", plot = TRUE)
plsR2 <- R2(plsrall_narrZn)$val[compsel]
predarray <- predict(plsrall_narrZn, ncomp=compsel)

dtformatmerge[, `:=`(plspred = exp(predarray),
                     plsresid = exp(predarray)-specmeaZn_mean,
                     plsresidper = (exp(predarray)-specmeaZn_mean)/specmeaZn_mean)]

png(file.path(moddir, paste0('residplot_allnarrZn.png')))
plsplot <- ggplot(dtformatmerge, aes(1000*plspred, 1000*specmeaZn_mean)) + 
  geom_abline() +
  geom_point(alpha = 0.6, size=1.5) +
  scale_y_log10(name = expression(Observed~average~Zn~PM[2.5]~(mg/m^3)),limits = c(0.7, 120)) +
  scale_x_log10(name = expression(Predicted~average~Zn~PM[2.5]~(mg/m^3)), limits = c(0.7, 120)) +
  annotate("text", x=10, y=100, vjust=-1, #hjust=0.60, 
           label =paste0('B.~Meteorology+modeled~Zn~~R^2 ==', round(plsR2,2)), parse=T, size=4) + 
  coord_fixed(clip='off') +
  theme_classic() + 
  theme(text = element_text(size=12),
        plot.margin = margin(0.6, 0.25, 0, 0.25, "cm"),
        rect = element_rect(fill = "transparent"))
print(plsplot)
dev.off()
print(plsplot)
summary(plsrall_narrZn)
loadings(plsrall_narrZn)

png(file.path(moddir, 'AQImodels_plsrZn.png'), width=7, height=3.5, units='in', res=600)
grid.arrange(plsplotNARR, plsplot, ncol=2)
dev.off()

write.dbf(dtformatmerge[,.(UID, plspred, specmeaZn_mean, plsresid, plsresidper)], 
          file.path(moddir, 'AQIdata_analysis_dtformatmerge.dbf'))

#With NARR and Zn controlling for weekday
plsrall_narrZnMonday <- plsr(specmeaZnMonday ~ NARRvarsZn, ncomp = 10, data = plsall_format, validation = "LOO")
compsel <- selectNcomp(plsrall_narrZnMonday, nperm=1000, alpha=0.05, method = "randomization", plot = TRUE)
plsR2 <- R2(plsrall_narrZnMonday)$val[compsel]
predarray <- predict(plsrall_narrZnMonday, ncomp=compsel)

dtformatmerge[, `:=`(plspred = exp(predarray),
                     plsresid = specmeaZn_monday_mean-exp(predarray))]

png(file.path(moddir, paste0('residplot_allnarrZnMonday.png')))
plsplot <- qplot(dtformatmerge$plspred, dtformatmerge[,specmeaZn_monday_mean]) + 
  geom_abline() + 
  annotate("text", x=0.0015, y=0.1, label =paste0('R2:', round(plsR2,2))) + 
  coord_fixed() +
  theme_classic()
print(plsplot)
dev.off()
print(plsplot)
summary(plsrall_narrZnMonday)
loadings(plsrall_narrZnMonday)



############################################################################################################################################
############################################################################################################################################
###################################### NOT USED ############################################################################################
#------ PLSR for daily data of Zn/Al + NARR variables at each station with more than 200 obs -----
#station <- "4813916083899" #for debugging
for (station in selectedstations) { #will change to apply at some point
  print(station)
  #For each station
  stationdat <- airNARRZnAl_nona[UID==station,]
  
  #Develop plsr format to predict residuals from week days
  pls_format <- data.frame(
    weekday = stationdat$weekday,
    specmearatio = stationdat$specmearatio,
    logspecmearatio = log(stationdat$specmearatio),
    NARRvars = I(as.matrix(stationdat[,NARRcols, with=F]))
  )
  plsr1 <- plsr(logspecmearatio ~ NARRvars, ncomp = 10, data = pls_format, validation = "LOO")
  png(file.path(moddir, paste0('CVplotratio_', station, '.png')))
  compsel <- selectNcomp(plsr1, nperm=2000, alpha=0.1, method = "randomization", plot = TRUE)
  title(paste0('UID: ', station))
  dev.off()
  
  plsR2 <- R2(plsr1, ncomp=compsel)$val[compsel]
  predarray <- predict(plsr1, ncomp=compsel)
  
  if (compsel > 0){
    airNARRZnAl_nona[UID==station, plspred := exp(predarray)]
    png(file.path(moddir, paste0('residplotratio_', station, '.png')))
    plsplot <- ggplot(airNARRZnAl_nona[UID==station,], aes(x=plspred, y=specmearatio)) + 
      geom_point() +
      geom_abline() + 
      annotate("text", x=-Inf, y=Inf, label =paste0('R2:', round(plsR2,2)), hjust = -0.25, vjust = 1) +
      labs(title=paste0('UID: ', station)) + 
      scale_y_log10() + 
      scale_x_log10() + 
      theme_classic() 
    print(plsplot)
    dev.off()
    
    #Create dataframe where all meteorological variables are 0 (US-wide average for prediction)
    predat <- data.frame(matrix(data=0, nrow=nrow(stationdat), ncol=length(NARRcols)))
    colnames(predat) <- NARRcols
    predat$NARRvars <- I(as.matrix(predat[,NARRcols]))
    #Generate predicted Zn/Al level given mean meteorology + residuals
    airNARRZnAl_nona[UID==station, ratiospecpred := exp(predict(plsr1, ncomp=compsel, newdata = predat))][
      UID==station, `:=`(ratiospecpred_nometeo = plspred + as.vector(residuals(plsr1, comp=compsel)[1:nrow(stationdat), 1, compsel]),
                         meteocorrected = 1)]
  } else {
    airNARRZnAl_nona[UID==station, `:=`(ratiospecpred_nometeo  = specmearatio,
                                        meteocorrected = 0)]
  }
}

#Develop model based on weekday and long-term trend
glm_nometeo0 <- glm(ratiospecpred_nometeo ~ 1, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo0)

glm_nometeo1 <- glm(ratiospecpred_nometeo ~ weekday, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo1)

glm_nometeo2 <- glm(ratiospecpred_nometeo ~ weekday + datenum, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo2)

glm_nometeo3 <- glm(ratiospecpred_nometeo ~ 0 + UID + weekday + datenum, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo3)

glm_nometeo4 <- glm(ratiospecpred_nometeo ~ 0 + predZn + weekday + datenum, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_nometeo4)

#
glm3coefs <- coefficients(glm_nometeo3)[grep('UID.*', names(coefficients(glm_nometeo3)), perl=T)]
dtformat <- data.table(UID = gsub('UID', '', names(glm3coefs)), ratiospec_mean = glm3coefs)
predZnformat <- airNARRZnAl_nona[UID %in% selectedstations, .(predZn = max(predZn),
                                                              specmearatio_nometeomean = mean(ratiospecpred_nometeo),
                                                              specmearatio_nometeomedian = median(ratiospecpred_nometeo),
                                                              specmearatio_nometeo70th = quantile(ratiospecpred_nometeo, 0.7),
                                                              specmearatio_nometeo90th = quantile(ratiospecpred_nometeo, 0.9),
                                                              specmearatio_mean = mean(specmearatio),
                                                              specmearatio_median = median(specmearatio),
                                                              specmeaZn_mean = mean(specmeaZn),
                                                              specmeaZn_median = median(specmeaZn), 
                                                              nlcd_imp_A = mean(nlcd_imp_A),
                                                              aadtlog100frt = mean(aadtlog100frt),
                                                              N = .N), by=.(UID, meteocorrected)]
NARRcolsmean <- airNARRZnAl_nona[UID %in% selectedstations, lapply(.SD, mean(x, na.rm=T)), 
                                 .SDcols = NARRcols, by=UID]

dtformatmerge <- merge(merge(dtformat, predZnformat, by='UID'),
                       NARRcolsmean, by='UID')

ggplot(dtformatmerge, aes(x=predZn, y=ratiospec_mean, color=factor(meteocorrected))) + 
  geom_point() + 
  geom_smooth()

ggplot(dtformatmerge[meteocorrected == 1,], aes(x=predZn, y=ratiospec_mean)) + 
  geom_point() + 
  geom_smooth() +
  scale_y_log10() + 
  scale_x_log10()

lm_nometeonoweek <- lm(log(specmearatio_nometeomean)~log(predZn), data=dtformatmerge[meteocorrected == 1,])
summary(lm_nometeonoweek)

ggplot(dtformatmerge, aes(x=predZn, y=specmearatio_nometeomean, color=factor(meteocorrected))) + 
  geom_point() + 
  geom_smooth() +
  scale_y_log10() + 
  scale_x_log10()

lm_nometeo2 <- lm(log(specmearatio_nometeomean)~log(predZn), data=dtformatmerge[meteocorrected == 1,])
summary(lm_nometeo2)

ggplot(dtformatmerge, aes(x=predZn, y=specmearatio_nometeomedian, color=factor(meteocorrected))) + 
  geom_point() + 
  geom_smooth()+
  scale_y_log10() + 
  scale_x_log10()

ggplot(dtformatmerge[meteocorrected == 1,], aes(x=predZn, y=specmearatio_nometeo70th)) + 
  geom_point() + 
  geom_smooth() + 
  scale_y_log10() + 
  scale_x_log10()

#------ PLSR for long-term Zn/Al data + NARR variables of all stations with more than 200 obs together -----
#Remove effects of weekday and method
glm0 <- glm(specmearatio ~ 1, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm0)

glm_weekday <- glm(log(specmearatio) ~ weekday*UID, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_weekday)

glm_weekdaymethod <- glm(log(specmearatio) ~ weekday*UID + MethodCode, data=airNARRZnAl_nona[UID %in% selectedstations,])
summary(glm_weekdaymethod)
#coefficients(glm_weekdaymethod)

#Check long-term relationships with NARR variables
for (var in NARRcols) {
  print(var)
  png(file.path(moddir, paste0('NARRratiospecplot_', var, '.png')))
  narrplot <- ggplot(dtformatmerge, aes_string(x=var, y='specmearatio_mean', color='predZn')) + 
    geom_point() + 
    geom_smooth() +
    scale_color_distiller(palette='Spectral') +
    theme_classic()
  print(narrplot)
  dev.off()
}

#Develop plsr format to predict residuals from week days
plsall_format <- data.frame(
  specmearatio = dtformatmerge$specmearatio_mean,
  predZn = dtformatmerge$predZn,
  logspecmearatio = log(dtformatmerge$specmearatio_mean),
  NARRvars = I(as.matrix(dtformatmerge[,NARRcols, with=F])),
  NARRvarsZn = I(as.matrix(dtformatmerge[,c(NARRcols, 'predZn'), with=F]))
)

#With just NARR variables
plsrall_narr <- plsr(logspecmearatio ~ NARRvars, ncomp = 10, data = plsall_format, validation = "LOO")
compsel <- selectNcomp(plsrall_narr, nperm=5000, alpha=0.05, method = "randomization", plot = TRUE)
plsR2 <- R2(plsrall_narr, ncomp=compsel)$val[2]
predarray <- predict(plsrall_narr, ncomp=compsel)

dtformatmerge[, `:=`(plspred = exp(predarray),
                     plsresid = specmearatio_mean-exp(predarray))]

png(file.path(moddir, paste0('residplotratio_allnarr.png')))
plsplot <- qplot(dtformatmerge$plspred, dtformatmerge[,specmearatio_mean]) + 
  geom_abline() + 
  scale_y_log10() +
  scale_x_log10() +
  annotate("text", x=0.08, y=2, label =paste0('R2:', round(plsR2,2))) + 
  coord_fixed() +
  theme_classic()
print(plsplot)
dev.off()

ggplot(dtformatmerge, aes(x=predZn, y=plsresid)) + 
  geom_point() + 
  geom_smooth() + 
  scale_y_sqrt() + 
  scale_x_sqrt()
gam_plsrallresid<- (mgcv::gam(plsresid~s(predZn), data=dtformatmerge))
summary(gam_plsrallresid)
plot(gam_plsrallresid, rug = TRUE, residuals = TRUE, pch = 1, cex = 1)

#With NARR and Zn
plsrall_narrZn <- plsr(logspecmearatio ~ NARRvarsZn, ncomp = 10, data = plsall_format, validation = "LOO")
compsel <- selectNcomp(plsrall_narrZn, nperm=5000, alpha=0.05, method = "randomization", plot = TRUE)
plsR2 <- R2(plsrall_narrZn, ncomp=compsel)$val[2]
predarray <- predict(plsrall_narrZn, ncomp=compsel)

dtformatmerge[, `:=`(plspred = exp(predarray),
                     plsresid = specmearatio_mean-exp(predarray))]

png(file.path(moddir, paste0('residplotratio_allnarrZn.png')))
plsplot <- qplot(dtformatmerge$plspred, dtformatmerge[,specmearatio_mean]) + 
  geom_abline() + 
  scale_y_log10() +
  scale_x_log10() +
  annotate("text", x=0.08, y=2, label =paste0('R2:', round(plsR2,2))) + 
  coord_fixed() +
  theme_classic()
print(plsplot)
dev.off()
print(plsplot)
summary(plsrall_narrZn)
loadings(plsrall_narrZn)



dtformatmerge


#------ GLM/GAM for long-term Zn data without meteorology of all stations with more than 200 obs together -----
ggplot(dtformatmerge, aes(x=predZn, y=log(specmeaZn_mean))) + 
  geom_point() + 
  geom_smooth()
lm1 <- lm(log(specmeaZn_mean)~log(predZn), data=dtformatmerge)
summary(lm1)

glm1 <- glm(specmeaZn_mean~log(predZn), data=dtformatmerge, family=gaussian('log'))
summary(glm1)


gam1 <- mgcv::gam(log(specmeaZn_mean)~s(predZn, k=3), data=dtformatmerge)
summary(gam1)
plot(gam1, residuals=T)

ggplot(dtformatmerge, aes(x=predZn, y=log(specmeaZn_median))) + 
  geom_point() + 
  geom_smooth()
ggplot(dtformatmerge, aes(x=nlcd_imp_A, y=log(specmeaZn_median))) + 
  geom_point() + 
  geom_smooth()
ggplot(dtformatmerge, aes(x=aadtlog100frt, y=log(specmeaZn_median))) + 
  geom_point() + 
  geom_smooth()

lm1 <- lm(log(specmeaZn_median)~log(predZn), data=dtformatmerge)
summary(lm1)

lm2 <- lm(log(specmeaZn_median)~nlcd_imp_A + aadtlog100frt, data=dtformatmerge)
summary(lm2)

lm3 <- lm(log(specmeaZn_median)~nlcd_imp_A*aadtlog100frt, data=dtformatmerge)
summary(lm3)

gam1 <- mgcv::gam(log(specmeaZn_median)~s(predZn, k=4), data=dtformatmerge)
summary(gam1)
plot(gam1, residuals=T)

gam2 <- mgcv::gam(log(specmeaZn_median)~s(nlcd_imp_A, k=4), data=dtformatmerge)
summary(gam2)
plot(gam2, residuals=T)

gam3 <- mgcv::gam(log(specmeaZn_median)~s(nlcd_imp_A, k=4) + s(aadtlog100frt, k=3), data=dtformatmerge)
summary(gam3)
plot(gam3, residuals=T)


#------ GLM/GAM for long-term Zn/Al data without meteorology of all stations with more than 200 obs together -----
ggplot(dtformatmerge, aes(x=predZn, y=log(specmearatio_mean))) + 
  geom_point() + 
  geom_smooth()
lm1 <- lm(log(specmearatio_mean)~log(predZn), data=dtformatmerge)
summary(lm1)

glm1 <- glm(specmearatio_mean~log(predZn), data=dtformatmerge, family=gaussian('log'))
summary(glm1)


gam1 <- mgcv::gam(log(specmearatio_mean)~s(predZn, k=3), data=dtformatmerge)
summary(gam1)
plot(gam1, residuals=T)

ggplot(dtformatmerge, aes(x=predZn, y=log(specmearatio_median))) + 
  geom_point() + 
  geom_smooth()
ggplot(dtformatmerge, aes(x=nlcd_imp_A, y=log(specmearatio_median))) + 
  geom_point() + 
  geom_smooth()
ggplot(dtformatmerge, aes(x=aadtlog100frt, y=log(specmearatio_median))) + 
  geom_point() + 
  geom_smooth()

lm1 <- lm(log(specmearatio_median)~log(predZn), data=dtformatmerge)
summary(lm1)

lm2 <- lm(log(specmearatio_median)~nlcd_imp_A + aadtlog100frt, data=dtformatmerge)
summary(lm2)

lm3 <- lm(log(specmearatio_median)~nlcd_imp_A*aadtlog100frt, data=dtformatmerge)
summary(lm3)

gam1 <- mgcv::gam(log(specmearatio_median)~s(predZn, k=4), data=dtformatmerge)
summary(gam1)
plot(gam1, residuals=T)

gam2 <- mgcv::gam(log(specmearatio_median)~s(nlcd_imp_A, k=4), data=dtformatmerge)
summary(gam2)
plot(gam2, residuals=T)

gam3 <- mgcv::gam(log(specmearatio_median)~s(nlcd_imp_A, k=4) + s(aadtlog100frt, k=3), data=dtformatmerge)
summary(gam3)
plot(gam3, residuals=T)





#################################################
##### 1. Analyze annual summary data ####
############################################################################################################################################
#---- Import data ----
AQIsites <- as.data.table(readOGR(dsn = file.path(resdir, 'airdata'), layer = 'airsites20190315'))
AQI2018<- fread(file.path(datadir, 'EPA_AirData_201902/annual_conc_by_monitor_2018/annual_conc_by_monitor_2018.csv'))
AQI2017<- fread(file.path(datadir, 'EPA_AirData_201902/annual_conc_by_monitor_2017/annual_conc_by_monitor_2017.csv'))
AQI2016<- fread(file.path(datadir, 'EPA_AirData_201902/annual_conc_by_monitor_2016/annual_conc_by_monitor_2016.csv'))
AQI2015<- fread(file.path(datadir, 'EPA_AirData_201902/annual_conc_by_monitor_2015/annual_conc_by_monitor_2015.csv'))

#---- Get scaling for pollution variables ----
scaling <- readRDS(file.path(moddir, 'fieldXRFmodels_scaling.rds'))

#---- Predict Zn for each site ----
AQIsites[, heatbing1902log300 := 10*hbglog300/scaling$heatbing1902log300]
AQIsites[, heatSPDlog300 := 100*hSPDlog300/scaling$heatSPDlog300]
qplot(pollutfieldclean_cast$heatSPDlog300)
qplot(AQIsites$heatSPDlog300)


AQIsites[, predZn := exp(predict(modlistlogZn[[42]], 
                                 data.frame(heatbing1902log300, heatSPDlog300)))]


#---- Summarize 2015-2018 air monitoring data ----
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


##### 2. Analyze daily summary data with station-based meteorological data ####
############################################################################################################################################
#---- Import data ----
AQIsites <- as.data.table(readOGR(dsn = file.path(resdir, 'airdata'), layer = 'airsites'))
AQIfiles <- list.files(file.path(datadir, 'EPA_AirData_201902'), full.names = T, recursive = T)

PM10 <- as.data.table(do.call(rbind, lapply(grep('.*10SPEC.*[.]csv', AQIfiles, value=T), fread)))
SPEC<- as.data.table(do.call(rbind, lapply(grep('.*_SPEC.*[.]csv', AQIfiles, value=T), fread)))
temp <- as.data.table(do.call(rbind, lapply(grep('.*_TEMP.*[.]csv', AQIfiles, value=T), fread)))
wind <- as.data.table(do.call(rbind, lapply(grep('.*_WIND.*[.]csv', AQIfiles, value=T), fread)))
rh <- as.data.table(do.call(rbind, lapply(grep('.*_RH_DP.*[.]csv', AQIfiles, value=T), fread)))
press <- as.data.table(do.call(rbind, lapply(grep('.*_PRESS.*[.]csv', AQIfiles, value=T), fread)))

#---- Format data ----
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

#---- Inspect data ----
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


