source('00_packages.R')
source('00_functions.R')
source('00_dirstructure.R')

compute_elemcv <- function(dt, bycol=NULL)  {
  outdf <- as.data.frame(
    dt[!(dt$Elem %in% c('Rh','Pd','Ar')),
       paste0('m%CV: ',
              format(100*mean(cv, na.rm=T), digits=0),
              '%'),
       by=c('Elem', 'name', bycol)]
  )
  return(outdf)
}

# 1. Import data ------------------------------------
############################################################################################################################################
########### ---- A. Import and format field data ---- ####
fieldata <- read.csv(file.path(datadir,"field_data/field_data_raw_20190430_edit.csv"))
colnames(fieldata)[1] <- 'SiteID'
fieldata$SiteID <- as.character(fieldata$SiteID)
fieldata$Date <- as.Date(as.character(fieldata$Date), format='%m/%d/%Y')
fieldata_sel <- fieldata[!is.na(fieldata$XRFmin) & !is.na(fieldata$SiteID),] #Remove extraneous sites with no XRF data or just for TNC tour
fieldata_format <- data.frame()
#Create separate records for each XRF measurement (rather than one record with XRFmin and XRFmax)
for (row in seq(1,nrow(fieldata_sel))) {
  extract <- fieldata_sel[row,]
  for (xrf in seq(fieldata_sel[row,'XRFmin'], fieldata_sel[row,'XRFmax'])){
    #print(xrf)
    extract$XRFID <- xrf
    fieldata_format <- rbind(fieldata_format, extract)
  }
}
#Separate summer 2018 and spring 2019 records as XRF ids have been reset and are thus duplicates
fieldata_format[fieldata_format$Date<'2019/01/01', 'season'] <- 'summer2018'
fieldata_format[fieldata_format$Date>'2019/01/01', 'season'] <- 'spring2019'


########### ---- B. Import and format field XRF deconvolution results ---- ####
fieldXRF_summer2018 <- deconvolution_import(file.path(datadir, 'XRF20190501/PostBrukerCalibration/deconvolutionresults_XRF12_245_20181215'),
                                            idstart = 41, idend= 43)
fieldXRF_spring2019 <- deconvolution_import(file.path(datadir, 'XRF20190501/2019_MarchApril/deconvolutionresults_XRF39_161_20190513'),
                                            idstart = 41, idend= 43)

fieldXRF_summer2018[,'season'] <- 'summer2018'
fieldXRF_spring2019[,'season'] <- 'spring2019'

fieldXRF <- rbind(fieldXRF_summer2018, fieldXRF_spring2019)

########### ---- C. Import raw XRF lab data for inspection ---- ####
labxrf_list <- data.table(filename = grep('\\.csv$', list.files(file.path(datadir, 'XRF20190501/PelletMeasurements')) ,value=T))
for (i in labxrf_list$filename) {
  print(i)
  xrfrec <- read.csv(file.path(datadir, 'XRF20190501/PelletMeasurements', i))
  labxrf_list[filename == i, `:=`(cps = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Valid Count Last Packet',])),
                                  c_cum = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Valid Accumulated Counts',])),
                                  duration = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Duration Time',]))
  )]
}
setDT(labxrf_list)[, `:=`(c_cum_ps = c_cum/duration,
                          SiteID = regmatches(labxrf_list$filename, regexpr('[0-9A-Z]+[_][0-9]', labxrf_list$filename)))]
labxrf_list[, season := 'summer2018']

#List of initially problematic samples that had to be re-measured
# problem_recs <- setDT(labxrf_list)[SiteID %in% c('3A_1', '3A_2','6B_2','7A_1', '7A_2','19A_1', '19A_2', '20A_1',
#                                                  '20B_1', '20B_2','22A_1', '22A_2', '28B_2', '32A_1', '33B_1', 
#                                                  '33B_2','34A_1','36A_2', '37A_1', '37A_2','40A_1', '42A_1', 
#                                                  '42B_1','42B_2','44A_1', '44A_2', '45A_1', '45A_2','46A_1','46A_2',
#                                                  '48A_1','49A_1', '49A_2','53A_1', '53A_2', '54B_2','55B_1', '55B_2',
#                                                  '57A_1', '57A_2','59A_1', '59A_2', '61A_1', '61A_2'),]
labxrf_list[c_cum_ps > 90000,] #[!(labxrf_list[c_cum_ps > 90000, SiteID] %in% problem_recs$SiteID)] #Check whether any sites have anomalous photon counts

########### ---- D. Import and format lab XRF deconvolution results ---- ####
labXRF <- deconvolution_import(file.path(datadir, 'XRF20190501/PelletMeasurements/deconvolutionresults_labXRF_20181215'),
                               idstart = 29, idend= 33)
labXRF[, season := 'summer2018']
########### ---- E. Look at determinants of signalratio (noise) ----
fieldXRFsummary <- fieldXRF[,list(meannet = mean(Net),
                                  signalratio = mean(Net/Backgr.)), .(Element, Line, Energy.keV)] %>%
  .[,  `:=` (element_line = paste0(Element, '.', substr(Line, 1, 1)),
            type = 'In situ XRF')]

labXRFsummary <- labXRF[,list(meannet = mean(Net),
                              signalratio = mean(Net/Backgr.)), .(Element, Line, Energy.keV)] %>%
  .[, `:=` (element_line = paste0(Element, '.', substr(Line, 1, 1)),
            type = 'Laboratory XRF')]

XRFesummary <- rbind(fieldXRFsummary, labXRFsummary, use.names= T)

keVmaxnet <- unique(XRFesummary) %>% #Only keep Element-Line-Energy characteristics of each element
  setkey(meannet) #%>% #Classify from low to high energy 
  #.[,.SD[.N], by=Element] #Only keep line with the highest average net photon count

ggplot(keVmaxnet, aes(x=meannet, y=signalratio)) +
  geom_text(aes(label=element_line)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(span=1) +
  theme_bw() + 
  theme(text=element_text(size=14))

png(file.path(inspectdir, 'signalratio_energy.png'), 
    height = 10 , width = 8, units = 'in', res=600)
signalenergy_p <- ggplot(keVmaxnet[!(Element %in% c('Pd', 'Rh')) & type == 'In situ XRF',],
       aes(x=Energy.keV, y=signalratio)) +
  geom_smooth(span=1, method='gam', formula = y ~ s(x, bs = "cs", k=4), alpha = 1/4) + 
  geom_text(aes(label=element_line), alpha=0.75) +
  #labs(title ='In situ pXRF') + 
  scale_y_log10(name='In situ XRF - Net/Background photon count') +
  scale_x_continuous(limits= c(0,33.5), expand=c(0,0))  +
  theme_bw()+ 
  theme(text=element_text(size=14))

labfield_signal_p <- ggplot(dcast(keVmaxnet[!(Element %in% c('Pd', 'Rh')),],
             element_line~type, value.var = 'signalratio'),
       aes(x=`Laboratory XRF`, y=`In situ XRF`)) + 
  geom_text(aes(label=element_line), alpha=0.75) +
  scale_x_log10(name='Laboratory XRF - Net/Background photon count') + 
  scale_y_log10(name='In situ XRF - Net/Background photon count') + 
  geom_abline() +
  theme_bw()+ 
  theme(text=element_text(size=14))
patchwork <- signalenergy_p/labfield_signal_p
patchwork + plot_annotation(tag_levels = 'A')
dev.off()

########### ---- F. Import ICP-OES data ---- ####
ICPdat <- read_excel(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'),
                     sheet = 1)
colnames(ICPdat)[1] <- 'SAMPLE.SET'
ICPthresholds <- as.data.frame(
  read_excel(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'),
             sheet = 2))
rownames(ICPthresholds) <- gsub('\\W', '', ICPthresholds[,1]) #Clean sample IDs
ICPthresholds_format <- as.data.table(t(ICPthresholds[,-1])) %>%
  .[, Elem := colnames(ICPthresholds)[-1]]

########### ---- G. Import GIS data (including pollution variables) ---- ####
trees <- as.data.table(readOGR(dsn = file.path(resdir, 'pollution_variables.gdb'), layer = 'XRFsites_aea'))
#summary(trees)
heatcols <- colnames(trees)[grep('heat', colnames(trees))]
trees[, (heatcols) := lapply(.SD, function(x){x[is.na(x)] <- 0; x}), .SDcols = heatcols]
colnames(trees)[1] <- 'SiteID'
trees <- trees[!(SiteID %in% c('NA', NA)),]
trees[, NLCD_reclass_final_PS := as.factor(NLCD_reclass_final_PS)]

#Separate summer 2018 and spring 2019 records as XRF ids have been reset and are thus duplicates
trees[trees$Date<'2019/01/01', 'season'] <- 'summer2018'
trees[trees$Date>'2019/01/01', 'season'] <- 'spring2019'

########### ---- H. Define elements associated with car traffic ---- ####
#Elements present in brake linings and emitted brake dust (from Thorpe et al. 2008)
brakelining_elem <- c('Al', 'As', 'Ba', 'Ca', 'Cd', 'Co' ,'Cr', 'Cu', 'Fe', 'K', 'Li', 'Mg', 
                      'Mn', 'Mo', 'Na', 'Ni', 'Pb', 'Sb', 'Se', 'Sr', 'Zn') 
#Elements present in passenger car tyre tread (from Thorpe et al. 2008)
tire_elem <- c('Al', 'Ba', 'Ca', 'Cd', 'Co', 'Cr', 'Cu', 'Fe', 'K', 'Mg', 'Mn', 'Na', 'Ni',
               'Pb', 'Sb', 'Sr', 'Ti', 'Zn')

############################################################################################################################################


# 2. Format XRF and ICP data -------------------------------
############################################################################################################################################
########### ---- A. Format XRF data ---- ####
# ---- 1. Cast while summing net photon counts across electron transitions ----
fieldXRFcast <- dcast(setDT(fieldXRF), XRFID+season~Element, value.var='Net', fun.aggregate=sum) 
fieldXRFcast[, XRFID := as.numeric(gsub('[_]', '', XRFID))] #Format site number
fieldXRFcast[fieldXRFcast < 0] <- 0 #Floor negative net photon count to 0

labXRFcast <- dcast(setDT(labXRF), XRFID+season~Element, value.var='Net', fun.aggregate=sum)
labXRFcast[labXRFcast < 0] <- 0

# ---- 2. Normalize data by Rhodium photon count for field and lab results ----
fieldXRFcastnorm <- fieldXRFcast[, lapply(.SD, function(x) {x/Rh}), by = .(XRFID, season)]
labXRFcastnorm <- labXRFcast[, lapply(.SD, function(x) {x/Rh}), by = .(XRFID, season)]

# ---- 3. Merge datasets: lab XRF + field XRF + field variables ----
fieldt <- setDT(fieldata_format)[fieldXRFcastnorm, on=.(XRFID, season)]

labXRFcastnorm[, `:=`(SiteID = gsub('[A-B].*', '', XRFID),
                      Pair = gsub('[0-9_]+', '', XRFID),
                      XRFID = NULL)]
labdt <- setDT(fieldata_sel)[labXRFcastnorm, on =  .(SiteID, Pair)]

# ---- 4. Compute average, sd, and range for lab XRF results over multiple measurements for a given pellet ----
elemcols <- which(colnames(labdt) %in% periodicTable$symb)
lab_artaxstats <- labdt[, sapply(.SD, function(x) list(mean=mean(x, na.rm=T),
                                                       sd=sd(x, na.rm=T),
                                                       range=max(x,na.rm=TRUE)-min(x,na.rm=TRUE))), 
                        by=c('SiteID','Pair','Dry.start'), .SDcols = elemcols]
setnames(lab_artaxstats, c('SiteID', 'Pair','Dry.start', 
                           paste(rep(colnames(labdt)[elemcols],each=3),
                                 c('mean', 'sd', 'range'), sep='_')))

labXRF_format <- melt(lab_artaxstats, id.vars=c('SiteID','Pair','Dry.start'),
                      variable.name='Elem_stats') %>%
  .[, `:=`(Elem = sub('(.*)[_](.*)', '\\1', Elem_stats),
           stats = sub('(.*)[_](.*)', '\\2', Elem_stats))] %>% #Replaces the whole match with the first group in the regex pattern
  dcast(SiteID+Pair+Dry.start+Elem~stats, value.var = 'value') %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

# ---- 5. Compute average, sd, and range, then melt XRF field data ----
elemcols <- which(colnames(fieldt) %in% periodicTable$symb)
field_artaxstats <- fieldt[, sapply(.SD, function(x) list(mean=mean(x, na.rm=T),
                                                          sd=sd(x, na.rm=T),
                                                          range=max(x,na.rm=TRUE)-min(x,na.rm=TRUE))), 
                           by=c('SiteID','Pair','season', 'Date'), .SDcols = elemcols]
setnames(field_artaxstats, c('SiteID', 'Pair', 'season', 'Date', 
                             paste(rep(colnames(fieldt)[elemcols],each=3),
                                   c('mean', 'sd', 'range'), sep='_')))

fieldXRF_format <- melt(field_artaxstats,
                        id.vars=c('SiteID','Pair', 'season', 'Date'), 
                        variable.name='Elem_stats') %>%
  .[, `:=`(Elem = sub('(.*)[_](.*)', '\\1', Elem_stats),
           stats = sub('(.*)[_](.*)', '\\2', Elem_stats))] %>% #Replaces the whole match with the first group in the regex pattern
  dcast(SiteID+Pair+season+Elem+Date~stats, value.var = 'value') %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

# ---- 6. Write out mean XRF values (should write this directly to XRF sites proj) ----
field_artaxmean <- field_artaxstats[, .SD, .SDcols = c(1,2, grep('mean', colnames(field_artaxstats)))] 
colnames(field_artaxmean) <- gsub('_mean', '', colnames(field_artaxmean))
write.dbf(field_artaxmean, 'field_artaxmean_20190501.dbf')

lab_artaxmean <- lab_artaxstats[, .SD, .SDcols = c(1,2, grep('mean', colnames(lab_artaxstats)))] 
colnames(lab_artaxmean) <- gsub('_mean', '', colnames(lab_artaxmean))

# ---- 7. Check field XRF data distribution by element then transform and standardize ----
transcols <- c('transmean', 'tukeylambda')

#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
fieldXRF_format[Elem %in% fieldXRF_format[, length(unique(mean))>1, by=Elem][V1==T, Elem],
                (transcols):= transformTukey_lambda(mean, start = -2.5, end = 2.5,int = 0.025, 
                                                   rastertab=F,rep=100, verbose = FALSE, 
                                                   statistic = 1), by=Elem]

#Run a shapiro test on each element
field_transmean <- dcast(fieldXRF_format, SiteID+Pair~Elem, value.var = 'transmean') 
normtest <- shapiro_test_df(field_transmean[,sapply(field_transmean, class)=='numeric' &
                                              sapply(field_transmean, function(x) length(unique(x)))>1, with=F])
normsig <- data.frame(sig = normtest$significance)
normsig$Elem <- row.names(normsig)
#Plot histogram of transformed data color-coded by whether normally distributed or not
fieldnorm_join <- fieldXRF_format[normsig, on='Elem']
png(file.path(inspectdir, paste0('fieldXRF_TukeyTransform.png')), width = 20, height=12, units='in', res=300)
ggplot(fieldnorm_join, aes(x=transmean, fill=sig)) + 
  geom_histogram() + 
  facet_wrap(~Elem, scales='free') + 
  theme_classic()
dev.off()            
#Mo and Na only non-normal ones after transformation

#z-scale data by element
fieldXRF_format[, transmean:=scale(transmean), by=Elem]

# ---- 8. Check lab XRF data distribution by element then transform and standardize ----
#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
labXRF_format[Elem %in% labXRF_format[, length(unique(mean))>1, by=Elem][V1==T, Elem],
              (transcols) := transformTukey_lambda(mean, start = -2, end = 2,int = 0.025, 
                                                 rastertab=F,rep=100, verbose = FALSE, 
                                                 statistic = 1), by=Elem]

#Run a shapiro test on each element
lab_transmean <- dcast(labXRF_format, SiteID+Pair~Elem, value.var = 'transmean') 
normtest <- shapiro_test_df(lab_transmean[,sapply(lab_transmean, class)=='numeric' &
                                            sapply(lab_transmean, function(x) length(unique(x)))>1, with=F])
normsig <- data.frame(sig = normtest$significance)
normsig$Elem <- row.names(normsig)
#Plot histogram of transformed data color-coded by whether normally distributed or not
labnorm_join <- labXRF_format[normsig, on='Elem']
png(file.path(inspectdir, paste0('labXRF_TukeyTransform.png')), width = 20, height=12, units='in', res=300)
ggplot(labnorm_join, aes(x=transmean, fill=sig)) + 
  geom_histogram() + 
  facet_wrap(~Elem, scales='free') + 
  theme_classic()
dev.off()       

#As, Mo, Na, Se, and Si still really not normal
#Rb and Ar almost normal

#z-scale data by element
labXRF_format[, transmean:=scale(transmean), by=Elem]

# ---- 9. Re-cast data for multivariate analysis ----
field_transmean <- dcast(fieldXRF_format, SiteID+Pair~Elem, value.var = 'transmean') 
lab_transmean <- dcast(labXRF_format, SiteID+Pair~Elem, value.var = 'transmean') 

########### ---- B. Format ICP-OES data ---- ####
# ---- 1. Format data and check site overlap between field XRF and ICP----
ICPdat[ICPdat == 'TR'] <- '0'
ICPdat[ICPdat == 'ND'] <- '0'
ICPdat <- setDT(ICPdat)[!(SAMPLE.SET %in% c('BLK', 'QC-1', NA)),]
ICPdat <- ICPdat[, lapply(.SD, as.character), by=SAMPLE.SET]
ICPdat <- ICPdat[, lapply(.SD, as.numeric), by=SAMPLE.SET]

#Check what sites are in ICP but not in field XRF
ICPdat[,unique(SAMPLE.SET)][!(ICPdat[,unique(SAMPLE.SET)] %in% fieldXRF_format[,unique(paste0(SiteID, Pair))])] 
#Site 22: XRF analyzer stopped working
#Site 5A: moss too high to reach with gun
#Site 63B: typo by lab analyst, should be 63A
ICPdat[SAMPLE.SET=='63B', SAMPLE.SET:='63A']
ICPdat[SAMPLE.SET=='63B-DUP', SAMPLE.SET:='63A-DUP']

#Check what sites are in XRF but not in ICP
fieldXRF_format[,unique(paste0(SiteID, Pair))][!(fieldXRF_format[,unique(paste0(SiteID, Pair))] %in% ICPdat[,unique(SAMPLE.SET)])] 
table(ICPdat$SAMPLE.SET)
#1A and 1B were taken from the wrong moss
#61A and 53A: insufficient sample to process in ICP-OES
#36B: typo by lab analyst, put 36A twice
ICPdat[SAMPLE.SET=='36A' & Cu == 7.24874371859296, SAMPLE.SET:='36B']

#Set aside and compute mean over duplicate measurements
ICPdup <- ICPdat[substr(ICPdat$SAMPLE.SET, 1,3) %in% substr(grep('DUP', ICPdat$SAMPLE.SET, value=T), 1,3),]
ICPdat[, SAMPLE.SET := substr(SAMPLE.SET, 1,3)] 
ICPmean <- ICPdat[, lapply(.SD, mean), by= SAMPLE.SET]

# ---- 2. Melt dataset to merge with XRF data ----
ICPmelt <- melt(setDT(ICPmean), id.vars = 'SAMPLE.SET', variable.name = 'Elem', value.name = 'ICP')
ICPmelt[, `:=`(SiteID = gsub('[A-Z]', '', SAMPLE.SET),
               Pair = gsub('[0-9]', '', SAMPLE.SET))]

# ---- 3. Check ICP data distribution by element then transform and standardize ----
transcolsicp <- c('transICP', 'tukeylambda')
#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
ICPmelt[Elem %in% ICPmelt[, length(unique(ICP))>1, by=Elem][V1==T, Elem],
        (transcolsicp) := transformTukey_lambda(ICP, start = -2, end = 2,int = 0.025, 
                                                rastertab=F,rep=100, verbose = FALSE, 
                                                statistic = 1)[[1]], by=Elem]

#Run a shapiro test on each element
ICP_trans <- dcast(ICPmelt, SiteID+Pair~Elem, value.var = 'transICP') 
normtest <- shapiro_test_df(ICP_trans[,sapply(ICP_trans, class)=='numeric' &
                                        sapply(ICP_trans, function(x) length(unique(x)))>1, with=F])
normsig <- data.frame(sig = normtest$significance)
normsig$Elem <- row.names(normsig)
#Plot histogram of transformed data color-coded by whether normally distributed or not
ICPnorm_join <- ICPmelt[normsig, on='Elem']
png(file.path(inspectdir, paste0('ICP_TukeyTransform.png')), width = 20, height=12, units='in', res=300)
ggplot(ICPnorm_join, aes(x=transICP, fill=sig)) + 
  geom_histogram() + 
  facet_wrap(~Elem, scales='free') + 
  theme_classic()
dev.off()     

#As, Cd, Mo, and Se still very far from normality (mostly 0s with a few large values)
#Cr, Na, Ni. and Pb not great (still non-normal)

# z-standardize by element
ICPmelt[, transICP:=scale(transICP), by=Elem]

# ---- 4. Re-cast data for multivaiate analysis ----
ICP_trans <- dcast(ICPmelt, SiteID+Pair~Elem, value.var='transICP')

########### ---- C. Format pollutant data ---- ####
# ---- 1. Check pollutant data distribution by element then transform and standardize ----
#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
treesmelt <- melt(trees[!is.na(SiteID),], id.vars=grep('heat|NLCD_imp', colnames(trees),
                                                       ignore.case=T, value=T, invert=T), 
                  variable.name = 'pollutvar' , value.name = 'pollutvalue')
treesmelt[, transpollut := transformTukey_lambda(pollutvalue, start = -2, end = 2,int = 0.025, 
                                                 rastertab=F,rep=100, verbose = FALSE, 
                                                 statistic = 1)[[1]], by=pollutvar]
#Run a shapiro test on each element
trees_trans <- dcast(treesmelt, SiteID+Pair~pollutvar, value.var = 'transpollut') 

normtest <- shapiro_test_df(trees_trans[,sapply(trees_trans, class)=='numeric' &
                                          sapply(trees_trans, function(x) length(unique(x)))>1, with=F])
normsig <- data.frame(sig = normtest$significance)
normsig$pollutvar <- row.names(normsig)
#Plot histogram of transformed data color-coded by whether normally distributed or not
treesnorm_join <- treesmelt[normsig, on='pollutvar']
png(file.path(inspectdir, paste0('trees_TukeyTransform.png')), width = 20, height=12, units='in', res=300)
ggplot(treesnorm_join, aes(x=transpollut, fill=sig)) + 
  geom_histogram() + 
  facet_wrap(~pollutvar, scales='free') + 
  theme_classic()
dev.off()     

# z-standardize by element
treesmelt[, transpollut:=scale(transpollut), by= pollutvar]

# ---- 2. Re-cast data for multivariate analysis ----
trees_trans <- dcast(treesmelt, formula = SiteID+Pair~pollutvar,  value.var='transpollut')

########### ---- D. Merge datasets ---- ####
# ---- 1. Merge XRF field data with ICP data ----
ICPfieldmerge <- merge(ICPmelt, fieldXRF_format, by = c('SiteID', 'Pair', 'Elem'), all.x=F, all.y=F)
unique(ICPmelt$Elem)[!(unique(ICPmelt$Elem) %in% unique(fieldXRF_format$Elem))] #Check what elems are in vs out
unique(fieldXRF_format$Elem)[!(unique(fieldXRF_format$Elem) %in% unique(ICPmelt$Elem))] #Check what elems are in vs out

# ---- 2. Merge XRF field data with XRF lab data ----
joincols <- c('SiteID', 'Pair', 'Elem')
labfieldmerge <- labXRF_format[fieldXRF_format, on = joincols] %>%
  .[SiteID != 1 & Elem != 'Rh',]
setnames(labfieldmerge, 
         colnames(labXRF_format[,-joincols, with=F]),
         gsub('^(?!i)', 'lab_', colnames(labXRF_format[,-joincols, with=F]), perl=T))
setnames(labfieldmerge, 
         colnames(labfieldmerge[,-joincols, with=F]),
         gsub('^i[.]', 'field_', colnames(labfieldmerge[,-joincols, with=F]), perl=T))

# ---- 3. Merge XRF lab data with ICP data ----
ICPlabmerge <- merge(ICPmelt, labXRF_format, by = c('SiteID', 'Pair', 'Elem'), all.x=F, all.y=F)
unique(ICPmelt$Elem)[!(unique(ICPmelt$Elem) %in% unique(labXRF_format$Elem))] #Check what elems are in vs out
unique(labXRF_format$Elem)[!(unique(labXRF_format$Elem) %in% unique(ICPmelt$Elem))] #Check what elems are in vs out

# ---- 4. Join non-transformed/standardized pollutant data to sites ----
pollutfieldmerge <- fieldXRF_format[trees, on=c('SiteID','Pair')]
pollutlabmerge <- labXRF_format[trees, on=c('SiteID','Pair')]
pollutICPmerge <- ICPmelt[trees, on=c('SiteID','Pair')]

# ---- 5. Join transformed/standardized pollutant data to sites ----
trees_fieldxrf_trans <- trees_trans[field_transmean, on=c('SiteID','Pair')][
  !(is.na(Fe) | SiteID==1),]
trees_labxrf_trans <- trees_trans[lab_transmean, on=c('SiteID','Pair')][
  !(is.na(Fe) | SiteID==1),]
trees_ICP_trans <- trees_trans[ICP_trans, on=c('SiteID', 'Pair')][
  !(is.na(Fe) | SiteID==1),]

############################################################################################################################################
############################################################################################################################################


# 3. Inspect data and remove outliers ------------------
######################  ---- A. Field XRF ---- ###########
str(fieldXRF_format)
# ---- 1. Assess within-tree field variability ----
#Plot coefficient of variation distributions for every element
cvmean_labeltree <- compute_elemcv(fieldXRF_format) 

# #Plot relationship between elemental concentration and cv
# ggplotly(
#   ggplot(fieldXRF_format[!(fieldXRF_format$Elem %in% c('Rh','Pd','Ar')),],
#          aes(x=mean, y=cv, color=name, label=paste0(SiteID, Pair))) +
#     geom_point()+
#     geom_smooth() +
#     labs(x='Mean photon count (normalized)', y='Coefficient of variation') +
#     facet_wrap(~name, scales='free') +
#     theme_classic() +
#     theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
#           strip.text = element_text(size=14))
# )
# 
# #
for (elem in unique(fieldXRF_format[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
  print(elem)
  png(file.path(inspectdir, paste0('fieldXRF_withintree_',elem,'.png')), width = 20, height=12, units='in', res=300)
  print(
    ggplot(fieldt,
           aes(x=paste0(SiteID, Pair), y = get(elem), fill=Pair)) +
      geom_line(aes(group=paste0(SiteID, Pair)), color='black') +
      geom_point(size=5, colour='black', pch=21, alpha=0.75) +
      labs(x='Site', y='Mean net photon count') +
      theme_bw() +
      theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
            strip.text = element_text(size=14),
            panel.border = element_blank(),
            axis.text.x = element_text(angle=90),
            axis.line = element_line(color='black'))
  )
  dev.off()
}

# ---- 2. Assess within-site field variability ---- 
#TO DO: CREATE RANDOM PALETTE
field_artaxmeansite <- fieldt[,lapply(.SD, mean, na.rm=TRUE), by=c('SiteID'), .SDcols=31:59]
artaxsdsite <- fieldt[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), .SDcols=31:59]
fieldXRF_formatsite <- merge(melt(field_artaxmeansite, id.vars='SiteID', variable.name='Elem', value.name='mean'),
                             melt(artaxsdsite, id.vars='SiteID', variable.name='Elem', value.name='sd'),
                             by=c('SiteID','Elem')) %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')
cvmean_labelsite_field <- compute_elemcv(fieldXRF_formatsite) 

#Compare tree vs. site variability
duplisite <- unique(fieldXRF_format[duplicated(fieldXRF_format[, .(SiteID, name)]), SiteID])

treesitecv_compare <- merge(
  fieldXRF_formatsite[!(fieldXRF_formatsite$Elem %in% c('Rh','Pd','Ar')) &
                        SiteID %in% duplisite, 
                      list(elemcv_site = mean(cv, na.rm=T)), by=Elem],
  fieldXRF_format[!(fieldXRF_format$Elem %in% c('Rh','Pd','Ar')) &
                    SiteID %in% duplisite, 
                  list(elemcv_tree = mean(cv, na.rm=T)), by=Elem],
  by='Elem'
)


png(file.path(inspectdir, 'fieldXRF_CV.png'), 
    height = 7 , width = 6, units = 'in', res=300)
treecv_disp <- ggplot(fieldXRF_formatsite[!(fieldXRF_formatsite$Elem %in% c('Rh','Pd','Ar')),], 
                      aes(x=100*cv, fill=Elem)) + 
  geom_density()+
  geom_text(data=cvmean_labeltree, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~Elem, nrow = 6, ncol = 5) + 
  theme_bw() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        text = element_text(size=14),
        strip.text = element_text(size=10),
        panel.border = element_blank(),
        axis.text.x = element_text(angle=90),
        legend.position = "none",
        axis.line = element_line(color='black'))


treesite_cvp <- ggplot(treesitecv_compare, 
                       aes(x=100*elemcv_tree, y=100*elemcv_site)) + 
  # geom_text(aes(label = Elem), alpha=0.75) + 
  geom_abline(color='darkblue', alpha=1/2) +
  geom_text_repel(aes(label = Elem), alpha=0.75,
                  box.padding = 0, label.padding = 0, point.padding=0) +
  scale_x_log10(name='Tree - Mean % CV') + 
  scale_y_log10(name='Site - Mean % CV') + 
  theme_bw() + 
  theme(text = element_text(size=14))

grid.arrange(treecv_disp, treesite_cvp, layout_matrix = rbind(c(1), c(1), c(2)))
dev.off()
  
  
# for (elem in unique(fieldXRF_format[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
#   print(elem)
#   png(file.path(inspectdir, paste0('fieldXRF_withinsite_',elem,'.png')), width = 20, height=12, units='in', res=300)
#   print(
#     ggplot(fieldXRF_format[Elem == elem,], 
#            aes(x=SiteID, y = mean, fill=SiteID)) + 
#       geom_line(aes(group=SiteID), color='black') +
#       geom_point(size=5, colour='black', pch=21, alpha=0.75) +
#       geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color=SiteID)) +
#       labs(x='Element', y='Mean net photon count') + 
#       theme_bw() + 
#       theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
#             strip.text = element_text(size=14),
#             panel.border = element_blank(),
#             axis.line = element_line(color='black'))
#   )
#   dev.off()
# }

# ---- 3. Exclude columns that have high within-tree variability or that are unrelated to traffic pollution ----
fieldXRF_format[, `:=`(meanCV=mean(cv, na.rm=T),
                       sdCV = sd(cv, na.rm=T)), by=Elem]
excol <- unique(c(fieldXRF_format[meanCV>0.5,unique(Elem)], 
                  fieldXRF_format[!(fieldXRF_format$Elem %in% union(brakelining_elem, tire_elem)),Elem],
                  'Rh', 'Pd', 'Rb', 'Ag')) %>%
  setdiff('Zr') #Keep zirconium as appears to be correlated with traffic
excol2 <- c('As', 'Cd', 'Mo', 'Na', 'Se', 'Si') #Exclude columns that are very far from being normally distributed (mostly 0s in ICP data)

# ---- 4. Univariate flag based on within-tree CV for metals with CV < 0.5----
#NAcount <- function(x) length(which(!is.na(unlist(x))))
fieldXRF_format[(cv>(meanCV+2*sdCV)) & !(Elem %in% excol) & !is.na(cv),
                CVflag_count_field := .N, by=.(SiteID, Pair)][
                  , CVflag_count_field := ifelse(is.na(CVflag_count_field) | is.na(cv), 
                                                 as.integer(round(mean(CVflag_count_field, na.rm=T))),
                                                 CVflag_count_field),
                  by=.(SiteID, Pair)][
                    is.na(CVflag_count_field), CVflag_count_field := 0]

### ---- Multivariate outlier detection ---- 
"see https://stats.stackexchange.com/questions/213/what-is-the-best-way-to-identify-outliers-in-multivariate-data 
for reference discussion as well as https://rpubs.com/Treegonaut/301942
For other resources, see Research_resources/statistics/outliers"
# ---- 5. Compare classic vs robust Mahalanobis distance ----
outlierdat <- field_transmean[,-c('SiteID','Pair', excol), with=F]
par(mfrow=c(1,1))
distances <- dd.plot(outlierdat, quan=0.90, alpha=0.025)
outlierdat_dist <- cbind(field_transmean[,-excol, with=F], as.data.frame(distances))

ggplotly(
  ggplot(outlierdat_dist, aes(x=md.cla, y=md.rob, color=outliers, label=paste0(SiteID,Pair))) +
    #geom_point(size=3) +
    geom_text() + 
    theme_classic()
)

# ---- 6. Filzmoser et al. multivariate outlier detection ----
outliers <- aq.plot(outlierdat, delta=qchisq(0.975, df = ncol(outlierdat)), quan = 0.75, alpha = 0.05)
par(mfrow=c(1,1))
outlierdat_dist <- cbind(outlierdat_dist, aqplot_outliers = outliers$outliers)

outlierdat_distmelt <- melt(outlierdat_dist[, colnames(outlierdat_dist) %in% c('SiteID', 'Pair', 'aqplot_outliers', 'md.rob',periodicTable$symb), with=F], 
                            id.vars=c('SiteID','Pair','aqplot_outliers','md.rob'), variable.name='Elem')

#Check out univariate distribution of outliers detected by mvoutlier
ggplot(outlierdat_distmelt, aes(x='O', y=value, color=aqplot_outliers)) + 
  geom_jitter() +
  #scale_color_distiller(palette= 'Spectral') +.;iA
  facet_wrap(~Elem, scales='free') +
  theme_classic()

### ---- Outliers in relationship between field XRF and (lab XRF | ICP | pollutant drivers) ----
# ---- 7. Multivariate relationship to ICP-OES for the purpose of outlier detection ----
#RDA requirements: Y and X must be centered and Y must be standardized; collinearity among X variables should be reduced prior to RDA

#Only keep sites and elements that are in both ICP and field XRF
field_ICP_sites <- do.call(paste0,
                           c(intersect(ICP_trans[,c('SiteID','Pair'),with=F], field_transmean[,c('SiteID','Pair'),with=F])))
ICPrdaformat <- as.matrix(ICP_trans[paste0(SiteID,Pair) %in% field_ICP_sites,-c('SiteID','Pair', excol, excol2), with=F])
XRFrdaformat <- as.matrix(field_transmean[paste0(SiteID,Pair) %in% field_ICP_sites,-c('SiteID','Pair', excol, excol2), with=F])

#Run RDA
rda_fieldICP<- rda(ICPrdaformat ~ XRFrdaformat, scale=FALSE, na.action = na.fail)
summary(rda_fieldICP)
anova(rda_fieldICP) #Check RDA's significance through permutation
anova(rda_fieldICP, by='axis') #Check each RDA axis' significance through permutation

#RDA triplot based on different combinations of axes
plot(rda_fieldICP,choices=c(1,2),display=c('sp','bp'),scaling=1,
     xlim=c(-2, 10), ylim=c(-2, 3))
text(rda_fieldICP,choices=c(1,2),
     labels=ICP_trans[paste0(SiteID,Pair) %in% field_ICP_sites,paste0(SiteID, Pair)])

plot(rda_fieldICP,choices=c(1,3),display=c('sp','bp'),scaling=1,
     xlim=c(-2, 10), ylim=c(-2, 3))
text(rda_fieldICP,choices=c(1,3),
     labels=ICP_trans[paste0(SiteID,Pair) %in% field_ICP_sites,paste0(SiteID, Pair)])

plot(rda_fieldICP,choices=c(2,3),display=c('sp','bp'),scaling=1,
     xlim=c(-2, 10), ylim=c(-2, 3))
text(rda_fieldICP,choices=c(2,3),
     labels=ICP_trans[paste0(SiteID,Pair) %in% field_ICP_sites,paste0(SiteID, Pair)])

#Plot RDA studentized residuals (inspired from ordiresids)
rda_fitresid <- data.table(
  SiteIDPair = rep(ICP_trans[paste0(SiteID,Pair) %in% field_ICP_sites,paste0(SiteID, Pair)], ncol(ICPrdaformat)),
  elem = rep(colnames(ICPrdaformat), nrow(ICPrdaformat), 'each'),
  fitted = as.vector(sweep(fitted(rda_fieldICP, type = "working"), 2, sigma(rda_fieldICP), "/")),
  residuals = as.vector(rstudent(rda_fieldICP)))

ggplot(rda_fitresid, aes(x=fitted, y=residuals, label=paste0(SiteIDPair,'-', elem))) + 
  geom_text(position=position_jitter(width=1,height=0)) + 
  geom_abline(intercept=c(-2, 2), slope=0)

fieldICP_rda_resid <- rda_fitresid[, list(fieldICP_rda_resid = mean(residuals)), by=SiteIDPair]

# ---- 8. Multivariate relationship to lab XRF for the purpose of outlier detection ----
#Only keep sites and elements that are in both ICP and field XRF
field_lab_sites <- do.call(paste0,
                           c(intersect(lab_transmean[,c('SiteID','Pair'),with=F], field_transmean[,c('SiteID','Pair'),with=F])))
labrdaformat <- as.matrix(lab_transmean[paste0(SiteID,Pair) %in% field_lab_sites,-c('SiteID','Pair', excol, excol2), with=F])
fieldrdaformat <- as.matrix(field_transmean[paste0(SiteID,Pair) %in% field_lab_sites,-c('SiteID','Pair', excol, excol2), with=F])

#Run RDA
rda_fieldlab<- rda(labrdaformat ~ fieldrdaformat, scale=FALSE, na.action = na.fail)
summary(rda_fieldlab)
anova(rda_fieldlab) #Check RDA's significance through permutation
anova(rda_fieldlab, by='axis') #Check each RDA axis' significance through permutation

#RDA triplot based on different combinations of axes
plot(rda_fieldlab,choices=c(1,2),display=c('sp','bp'),scaling=1,
     xlim=c(-2, 10), ylim=c(-2, 3))
text(rda_fieldlab,choices=c(1,2),
     labels=lab_transmean[paste0(SiteID,Pair) %in% field_lab_sites,paste0(SiteID, Pair)])

plot(rda_fieldlab,choices=c(1,3),display=c('sp','bp'),scaling=1,
     xlim=c(-2, 10), ylim=c(-2, 3))
text(rda_fieldlab,choices=c(1,3),
     labels=lab_transmean[paste0(SiteID,Pair) %in% field_lab_sites,paste0(SiteID, Pair)])

plot(rda_fieldlab,choices=c(2,3),display=c('sp','bp'),scaling=1,
     xlim=c(-2, 10), ylim=c(-2, 3))
text(rda_fieldlab,choices=c(2,3),
     labels=lab_transmean[paste0(SiteID,Pair) %in% field_lab_sites,paste0(SiteID, Pair)])

#Plot RDA studentized residuals (inspired from ordiresids)
rda_fitresid <- data.table(
  SiteIDPair = rep(lab_transmean[paste0(SiteID,Pair) %in% field_lab_sites,paste0(SiteID, Pair)], ncol(labrdaformat)),
  elem = rep(colnames(labrdaformat), nrow(labrdaformat), 'each'),
  fitted = as.vector(sweep(fitted(rda_fieldlab, type = "working"), 2, sigma(rda_fieldlab), "/")),
  residuals = as.vector(rstudent(rda_fieldlab)))

ggplot(rda_fitresid, aes(x=fitted, y=residuals, label=paste0(SiteIDPair,'-', elem))) + 
  geom_text(position=position_jitter(width=1,height=0)) + 
  geom_abline(intercept=c(-2, 2), slope=0)


fieldlab_rda_resid <- rda_fitresid[, list(fieldlab_rda_resid = mean(residuals)), by=SiteIDPair]

# ---- 9. Multivariate relationship to pollution predictors for the purpose of outlier detection  ---- 
#Format data
trees_fieldxrf_trans_melt <- melt(trees_fieldxrf_trans[,-c(excol, excol2), with=F], 
                                  id.vars=colnames(trees_fieldxrf_trans)[!(colnames(trees_fieldxrf_trans) %in% periodicTable$symb)], 
                                  variable.name = 'Elem') %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

fieldrdaformat <- as.matrix(trees_fieldxrf_trans[, setdiff(
  colnames(trees_fieldxrf_trans)[colnames(trees_fieldxrf_trans) %in% periodicTable$symb], c(excol,excol2)), with=F])
pollutrdaformat <- as.matrix(trees_fieldxrf_trans[,grep('heat|NLCD', colnames(trees_fieldxrf_trans), value=T, ignore.case = T), 
                                                  with=F])

#Run RDA
rda_pollutfield<- rda(fieldrdaformat ~ pollutrdaformat, scale=FALSE, na.action = na.fail)
summary(rda_pollutfield)
anova(rda_pollutfield) #Check RDA's significance through permutation
anova(rda_pollutfield, by='axis') #Check each RDA axis' significance through permutation

#RDA triplot based on different combinations of axes
plot(rda_pollutfield,choices=c(1,2),display=c('sp','bp'),scaling=1,
     xlim=c(-2, 10), ylim=c(-2, 3))
text(rda_pollutfield,choices=c(1,2),
     labels=lab_transmean[paste0(SiteID,Pair) %in% field_lab_sites,paste0(SiteID, Pair)])

#Plot RDA studentized residuals (inspired from ordiresids)
rda_fitresid <- data.table(
  SiteIDPair = rep(trees_fieldxrf_trans[,paste0(SiteID, Pair)], ncol(fieldrdaformat)),
  elem = rep(colnames(fieldrdaformat), nrow(fieldrdaformat), 'each'),
  fitted = as.vector(sweep(fitted(rda_pollutfield, type = "working"), 2, sigma(rda_pollutfield), "/")),
  residuals = as.vector(rstudent(rda_pollutfield)))

ggplot(rda_fitresid, aes(x=fitted, y=residuals, label=paste0(SiteIDPair,'-', elem))) + 
  geom_text(position=position_jitter(width=1,height=0)) + 
  geom_abline(intercept=c(-2, 2), slope=0)

meanresid <- rda_fitresid[, mean(residuals), by=SiteIDPair]
# ---- 10. Univariate relationship to ICP-OES for the purpose of outlier detection ----
#Outlier diagnostic plots
ICPfieldmerge[, ICPfield_flags := 0]
for (chem in unique(ICPfieldmerge$Elem)) {
  print(chem)
  ICPfield_lm <- lm(ICP ~ mean, data = ICPfieldmerge[Elem == chem,])
  # ggsave(file.path(inspectdir, paste0('fieldXRFICP_regoutliers', chem, '.png')),
  #        chemregdiagnostic_custom(ICPfield_lm, ICPfieldmerge, chem,  
  #                                 flagcol = 'ICPfield_flags', tresids = TRUE),
  #        width = 20, height=12, units='in', dpi=300)
  ICPfieldmerge[Elem == chem, ICPfieldR2 := summary(ICPfield_lm)$adj.r.squared]
}

#Plot field XRF ~ ICP data
ICPfield_plot <- ggplot(ICPfieldmerge[!(is.na(ICP) |  
                                          Elem %in% c(excol, excol2)),], 
                        aes(x=mean, y=ICP, color=factor(ICPfield_flags), group=1)) + 
  geom_linerangeh(aes(xmin=mean-sd, xmax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
ICPfield_plot


#Investigate determinants of R2
ICPmeanR2 <- unique(ICPfieldmerge[, list(meanXRF = mean(mean, na.rm=T),
                                         meanICP = mean(ICP, na.rm=T),
                                         R2 = ICPfieldR2), by=Elem][
                                           keVmaxnet, on='Elem==Element'])
ggplot(ICPmeanR2, aes(x=meanXRF, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=meanICP, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=signalratio, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

summary(loess(R2~meanXRF, data=ICPmeanR2, span=1))
summary(loess(R2~meanXRF+signalratio, data=ICPmeanR2, span=1))

# ---- 11. Univariate relationship to lab XRF for the purpose of outlier detection ----
#Outlier diagnostic plots
labfieldmerge[, labfield_flags := 0]

for (chem in unique(labfieldmerge$Elem)) {
  print(chem)
  labfield_lm <- lm(lab_mean ~ field_mean, data = labfieldmerge[Elem == chem & !is.na(lab_mean),])
  # ggsave(file.path(inspectdir, paste0('fieldlab_regoutliers', chem, '.png')),
  #        chemregdiagnostic_custom(labfield_lm, labfieldmerge[!is.na(lab_mean),], chem,  flagcol = 'labfield_flags', tresids = TRUE),
  #        width = 20, height=12, units='in', dpi=300)
  labfieldmerge[Elem == chem & !is.na(lab_mean), labfieldR2 := summary(labfield_lm)$adj.r.squared]
}

#Plot field XRF ~ lab data
labfield_plot <- ggplot(labfieldmerge[!(is.na(lab_mean) |  
                                          Elem %in% c(excol, excol2)),], 
                        aes(x=field_mean, y=lab_mean, color=factor(labfield_flags), group=1)) + 
  geom_linerangeh(aes(xmin=field_mean-field_sd, xmax=field_mean+field_sd), alpha=1/2) + 
  geom_segment(aes(xend=field_mean, y=lab_mean-lab_sd, yend=lab_mean+lab_sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
labfield_plot

#Look at determinants of R2
labmeanR2 <- unique(labfieldmerge[, list(meanfield = mean(field_mean, na.rm=T),
                                         meanlab = mean(lab_mean, na.rm=T),
                                         R2 = labfieldR2), by=Elem])[
                                           keVmaxnet, on='Elem==Element']

ggplot(labmeanR2, aes(x=meanfield, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(labmeanR2, aes(x=meanlab, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(labmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(labmeanR2, aes(x=signalratio, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

# ---- 12. Univariate relationship to pollution predictors for the purpose of outlier detection ----
#Outlier diagnostic plots
pollutfieldmerge <- pollutfieldmerge[!is.na(Elem),]
pollutfieldmerge[, pollutfield_flags := 0]
for (chem in unique(pollutfieldmerge$Elem)) {
  print(chem)
  pollutfield_lm <- lm(mean ~ heatbing1902log300proj*heatsubAADTlog300*heatsubslopelog300*nlcd_imp_ps +
                         heatbustransitlog300 + heatbustransitlog300:heatsubslopelog300, 
                       data = pollutfieldmerge[Elem == chem ,])
  # ggsave(file.path(inspectdir, paste0('fieldpollution_regoutliers2', chem, '.png')),
  #        chemregdiagnostic_custom(pollutfield_lm, pollutfieldmerge, chem,  flagcol = 'pollutfield_flags', tresids = TRUE),
  #        width = 20, height=12, units='in', dpi=300)
  pollutfieldmerge[Elem == chem, `:=`(pollutpred = fitted(pollutfield_lm),
                                      pollutfieldR2 = summary(pollutfield_lm)$adj.r.squared)]
}

#Plot predicted field XRF ~ observed field XRF
pollutfield_plot <- ggplot(pollutfieldmerge[!(Elem %in% c(excol, excol2, NA)),], 
                           aes(x=pollutpred, y=mean, group=1, color=factor(pollutfield_flags))) + 
  geom_linerange(aes(ymin=mean-sd, ymax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  geom_smooth(method='lm') +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
pollutfield_plot


ggsave(file.path(inspectdir, 'fieldpollut_all.png'),
       pollutfield_plot,
       width = 20, height=12, units='in', dpi=600)

#Investigate determinants of R2
pollutmeanR2 <- unique(pollutfieldmerge[, list(meanXRF = mean(mean, na.rm=T),
                                               R2 = pollutfieldR2), by=Elem][
                                                 keVmaxnet, on='Elem==Element'])
ggplot(pollutmeanR2, aes(x=meanXRF, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=signalratio, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

### ---- 13. Compile outlier flags ----
#Merge all flag datasets
fieldXRF_formatflags <- merge(fieldXRF_format, ICPfieldmerge[, .(SiteID, Pair, Elem, ICPfield_flags, ICPfieldR2)], 
                              by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(labfieldmerge[, .(SiteID, Pair, Elem, labfield_flags, labfieldR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(pollutfieldmerge[, .(SiteID, Pair, Elem, pollutfield_flags, pollutfieldR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T)

ggplot(fieldXRF_formatflags, aes(x=pollutfieldR2, y=ICPfieldR2, label=Elem)) + 
  geom_text()
ggplot(fieldXRF_formatflags, aes(x=pollutfieldR2, y=labfieldR2, label=Elem)) + 
  geom_text()

#Compute total number of flags for each tree across all elements that have an R2 of at least 0.10 across models
#Don't want to count outliers from spurious models
flagelems <- fieldXRF_formatflags[, min(ICPfieldR2, labfieldR2, pollutfieldR2, na.rm=T)>0.10, by=Elem][
  V1==TRUE & Elem != 'Rh', Elem]
fieldXRF_formatflags[Elem %in% flagelems,
                     flagpartsum_field := sum(ICPfield_flags, labfield_flags, pollutfield_flags, na.rm=T), 
                     by=.(SiteID, Pair)][
                       , flagsum_field := sum(flagpartsum_field, CVflag_count_field, na.rm=T), by=.(SiteID, Pair)]

#Check total number of flags by each type of flag
fieldXRF_formatflags[Elem %in% flagelems,
                     `:=`(ICPfield_flags_sum = sum(ICPfield_flags, na.rm=T),
                          labfield_flags_sum = sum(labfield_flags, na.rm=T),
                          pollutfield_flags_sum = sum(pollutfield_flags, na.rm=T)),
                     by=.(SiteID, Pair)]

fieldXRF_formatflags_u <- unique(
  fieldXRF_formatflags[!is.na(flagpartsum_field),
                       .SD, 
                       .SDcols=c('SiteID', 'Pair', 'CVflag_count_field',
                                 grep('flag.*sum', colnames(fieldXRF_formatflags), value=T))])

ggplot(fieldXRF_formatflags_u, aes(x=ICPfield_flags_sum, y=labfield_flags_sum, 
                                   label=paste0(SiteID, Pair), color = CVflag_count_field)) + 
  scale_color_distiller(palette='YlGnBu') +
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5) 

ggplot(fieldXRF_formatflags_u, aes(x=ICPfield_flags_sum, y=pollutfield_flags_sum, 
                                   label=paste0(SiteID, Pair), color = CVflag_count_field)) + 
  scale_color_distiller(palette='YlGnBu') + 
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5)

#Records to inspect
flagelems
fieldXRF_inspect <- unique(fieldXRF_formatflags[(ICPfield_flags_sum>3) | 
                                                  (labfield_flags_sum > 3) | 
                                                  (pollutfield_flags_sum > 3), 
                                                .(SiteID, Pair)])
outliertrees <- trees[fieldXRF_inspect, on=c('SiteID', 'Pair')]
outlierlocs <- SpatialPointsDataFrame(coords = data.frame(outliertrees$POINT_X, outliertrees$POINT_Y),
                                      data= as.data.frame(outliertrees))
#View(outlierlocs@data)
# leaflet(data = outlierlocs) %>% addTiles() %>%
#   addMarkers(clusterOptions = markerClusterOptions(),
#              popup = ~paste0(SiteID, Pair))


###################### ---- B. Lab XRF ---- ####
# ---- 1. Assess within-pellet lab variability ----
#Plot coefficient of variation distributions for every element
cvmean_labelpellet_lab <- compute_elemcv(labXRF_format)

ggplot(labXRF_format[!(labXRF_format$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density() +
  facet_wrap(~name, scales='free') + 
  labs(x='Coefficient of variation', y='Count') + 
  geom_text(data=cvmean_labelpellet_lab, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  theme_classic() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14))

for (elem in unique(labXRF_format[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
  print(elem)
  png(file.path(inspectdir, paste0('labXRF_withintree_',elem,'.png')), width = 20, height=12, units='in', res=300)
  print(
    ggplot(labdt,
           aes(x=paste0(SiteID, Pair), y = get(elem), fill=Pair)) +
      geom_line(aes(group=paste0(SiteID, Pair)), color='black') +
      geom_point(size=5, colour='black', pch=21, alpha=0.75) +
      labs(x='Site', y='Mean net photon count') +
      theme_bw() +
      theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
            strip.text = element_text(size=14),
            panel.border = element_blank(),
            axis.line = element_line(color='black'))
  )
  dev.off()
}


# ---- 2. Assess within-site lab variability ---- 
lab_artaxmeansite <- labdt[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID'), .SDcols=30:52]
artaxsdsite <- labdt[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), .SDcols=30:52]
labXRF_formatsite <- merge(melt(lab_artaxmeansite, id.vars='SiteID', variable.name='Elem', value.name='mean'),
                           melt(artaxsdsite, id.vars='SiteID', variable.name='Elem', value.name='sd'),
                           by=c('SiteID','Elem')) %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')
cvmean_labelsite_lab <- compute_elemcv(labXRF_formatsite)

ggplot(labXRF_formatsite[!(labXRF_formatsite$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density()+
  geom_text(data=cvmean_labelsite_lab, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name, scales='free') + 
  theme_bw() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14),
        panel.border = element_blank(),
        axis.line = element_line(color='black'))

for (elem in unique(labXRF_format[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
  print(elem)
  png(file.path(inspectdir, paste0('labXRF_withinsite_',elem,'.png')), width = 20, height=12, units='in', res=300)
  print(
    ggplot(labXRF_format[Elem == elem,], 
           aes(x=SiteID, y = mean, fill=SiteID)) + 
      geom_line(aes(group=SiteID), color='black') +
      geom_point(size=5, colour='black', pch=21, alpha=0.75) +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color=SiteID)) +
      labs(x='Element', y='Mean net photon count') + 
      theme_bw() + 
      theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
            strip.text = element_text(size=14),
            panel.border = element_blank(),
            axis.line = element_line(color='black'))
  )
  dev.off()
}
# ---- 3. Univariate flag based on within-tree CV for metals with CV < 0.5----
labXRF_format[, `:=`(meanCV=mean(cv, na.rm=T),
                     sdCV = sd(cv, na.rm=T)), by=Elem]

#NAcount <- function(x) length(which(!is.na(unlist(x))))
labXRF_format[(cv>(meanCV+2*sdCV)) & !(Elem %in% excol) & !is.na(cv),
              CVflag_count_lab := .N, by=.(SiteID, Pair)][
                , CVflag_count_lab := ifelse(is.na(CVflag_count_lab) | is.na(cv), 
                                             as.integer(round(mean(CVflag_count_lab, na.rm=T))),
                                             CVflag_count_lab),
                by=.(SiteID, Pair)][
                  is.na(CVflag_count_lab), CVflag_count_lab := 0]


# ---- 4. Univariate relationship to ICP-OES for the purpose of outlier detection ----
#Outlier diagnostic plots
ICPlabmerge[, ICPlab_flags := 0]
for (chem in unique(ICPlabmerge$Elem)) {
  print(chem)
  ICPlab_lm <- lm(ICP ~ mean, data = ICPlabmerge[Elem == chem,])
  # ggsave(file.path(inspectdir, paste0('labXRFICP_regoutliers', chem, '.png')),
  #        chemregdiagnostic_custom(ICPlab_lm, ICPlabmerge, chem,  flagcol = 'ICPlab_flags', tresids = TRUE),
  #        width = 20, height=12, units='in', dpi=300)
  ICPlabmerge[Elem == chem, ICPlabR2 := summary(ICPlab_lm)$adj.r.squared]
}

#Plot lab XRF ~ ICP data
ICPlab_plot <- ggplot(ICPlabmerge[!(is.na(ICP) |  
                                      Elem %in% c(excol, excol2)),], 
                      aes(x=mean, y=ICP, color=factor(ICPlab_flags), group=1)) + 
  geom_linerangeh(aes(xmin=mean-sd, xmax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
ICPlab_plot

#Investigate determinants of R2
ICPmeanR2 <- unique(ICPlabmerge[, list(meanXRF = mean(mean, na.rm=T),
                                       meanICP = mean(ICP, na.rm=T),
                                       R2 = ICPlabR2), by=Elem][
                                         keVmaxnet, on='Elem==Element'])
ggplot(ICPmeanR2, aes(x=meanXRF, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=meanICP, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=signalratio, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

summary(loess(R2~meanXRF, data=ICPmeanR2, span=1))
summary(loess(R2~meanXRF+signalratio, data=ICPmeanR2, span=1))

# ---- 5. Univariate relationship to pollution predictors for the purpose of outlier detection ----
#Outlier diagnostic plots
pollutlabmerge[, pollutlab_flags := 0]
pollutlabmerge <- pollutlabmerge[!is.na(Elem),]
for (chem in unique(pollutlabmerge$Elem)) {
  print(chem)
  pollutlab_lm <- lm(mean ~  heatbing1902log300proj*heatsubAADTlog300*heatsubslopelog300*nlcd_imp_ps +
                       heatbustransitlog300 + heatbustransitlog300:heatsubslopelog300, 
                     data = pollutlabmerge[Elem == chem,])
  # ggsave(file.path(inspectdir, paste0('labpollution_regoutliers', chem, '.png')),
  #        chemregdiagnostic_custom(pollutlab_lm, pollutlabmerge, chem,  flagcol = 'pollutlab_flags', tresids = TRUE),
  #        width = 20, height=12, units='in', dpi=300)
  pollutlabmerge[Elem == chem, `:=`(pollutpred = fitted(pollutlab_lm),
                                    pollutlabR2 = summary(pollutlab_lm)$adj.r.squared)]
}

#Plot predicted lab XRF ~ observed lab XRF
pollutlab_plot <- ggplot(pollutlabmerge[!(Elem %in% c(excol, excol2, NA)),], 
                         aes(x=pollutpred, y=mean, group=1, color=factor(pollutlab_flags))) + 
  geom_linerange(aes(ymin=mean-sd, ymax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  geom_smooth(method='lm') +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()

ggsave(file.path(inspectdir, 'labpollution_all.png'),
       pollutlab_plot,
       width = 20, height=12, units='in', dpi=600)

#Investigate determinants of R2
pollutmeanR2 <- unique(pollutlabmerge[, list(meanXRF = mean(mean, na.rm=T),
                                             R2 = pollutlabR2), by=Elem][
                                               keVmaxnet, on='Elem==Element'])
ggplot(pollutmeanR2, aes(x=meanXRF, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=signalratio, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()


### ---- 6. Compile outlier flags ----
#Merge all flag datasets
labXRF_formatflags <- merge(labXRF_format, ICPlabmerge[, .(SiteID, Pair, Elem, ICPlab_flags, ICPlabR2)], 
                            by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(labfieldmerge[, .(SiteID, Pair, Elem, labfield_flags, labfieldR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(pollutlabmerge[, .(SiteID, Pair, Elem, pollutlab_flags, pollutlabR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T)

ggplot(labXRF_formatflags, aes(x=ICPlabR2, y=pollutlabR2, label=Elem)) + 
  geom_text()
ggplot(labXRF_formatflags, aes(x=labfieldR2, y=pollutlabR2, label=Elem)) + 
  geom_text()

#Compute total number of flags for each tree across all elements that have an R2 of at least 0.10 across models
#Don't want to count outliers from spurious models
flagelems <- labXRF_formatflags[, min(ICPlabR2, labfieldR2, pollutlabR2, na.rm=T)>0.20, by=Elem][
  V1==TRUE & Elem != 'Rh', Elem]
labXRF_formatflags[Elem %in% flagelems,
                   flagpartsum_lab := sum(ICPlab_flags, labfield_flags, pollutlab_flags, na.rm=T), 
                   by=.(SiteID, Pair)][
                     , flagsum_lab := sum(flagpartsum_lab, CVflag_count_lab, na.rm=T), by=.(SiteID, Pair)]

#Check total number of flags by each type of flag
labXRF_formatflags[Elem %in% flagelems,
                   `:=`(ICPlab_flags_sum = sum(ICPlab_flags, na.rm=T),
                        labfield_flags_sum = sum(labfield_flags, na.rm=T),
                        pollutlab_flags_sum = sum(pollutlab_flags, na.rm=T)),
                   by=.(SiteID, Pair)]

labXRF_formatflags_u <- unique(
  labXRF_formatflags[!is.na(flagpartsum_lab),
                     .SD, 
                     .SDcols=c('SiteID', 'Pair', 'CVflag_count_lab',
                               grep('flag.*sum', colnames(labXRF_formatflags), value=T))])

ggplot(labXRF_formatflags_u, aes(x=ICPlab_flags_sum, y=labfield_flags_sum, 
                                 label=paste0(SiteID, Pair), color = CVflag_count_lab)) + 
  scale_color_distiller(palette='YlGnBu') +
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5) 

ggplot(labXRF_formatflags_u, aes(x=ICPlab_flags_sum, y=pollutlab_flags_sum, 
                                 label=paste0(SiteID, Pair), color = CVflag_count_lab)) + 
  scale_color_distiller(palette='YlGnBu') + 
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5)

#Records to inspect
flagelems
labXRF_inspect <- unique(labXRF_formatflags[(ICPlab_flags_sum>3) | 
                                              (labfield_flags_sum > 3) | 
                                              (pollutlab_flags_sum > 3), 
                                            .(SiteID, Pair)])
outliertrees <- trees[labXRF_inspect, on=c('SiteID', 'Pair')]
outlierlocs <- SpatialPointsDataFrame(coords = data.frame(outliertrees$POINT_X, outliertrees$POINT_Y),
                                      data= as.data.frame(outliertrees))

###################### ---- C. ICP-OES ---- ####
# -----1. Assess duplicate variability (but n=2 only) ------
#For two duplicates compute CV
ICPCV <- melt(ICPdat[SAMPLE.SET %in% c('63A', '15A'),], variable.name = 'Elem') %>%
  .[, list(cv = sd(value)/mean(value)), by=.(SAMPLE.SET, Elem)] %>%
  .[, list(replicateICPcv = round(100*mean(cv))), by=Elem] %>%
  setkey(Elem)
  
# ---- 1. Assess within-site variability ----
ICPmean[, `:=`(SiteID = gsub('[A-Z]', '', SAMPLE.SET),
               Pair = gsub('[0-9]', '', SAMPLE.SET))]
ICPmeansite <- ICPmean[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID'), 
                       .SDcols=which(colnames(ICPmean) %in% periodicTable$symb)]
ICPsdsite <- ICPmean[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), 
                     .SDcols=which(colnames(ICPmean) %in% periodicTable$symb)]
ICP_formatsite <- merge(melt(ICPmeansite, id.vars='SiteID', variable.name='Elem', value.name='mean'),
                        melt(ICPsdsite, id.vars='SiteID', variable.name='Elem', value.name='sd'),
                        by=c('SiteID','Elem')) %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')
cvmean_labelsite <- as.data.frame(ICP_formatsite[!(ICP_formatsite$Elem %in% c('Rh','Pd','Ar')),
                                                 paste0('m%CV: ',
                                                        format(100*mean(cv, na.rm=T),digits=2)),
                                                 by=name])

ggplot(ICP_formatsite[!(ICP_formatsite$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density()+
  geom_text(data=cvmean_labelsite, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name, scales='free') + 
  theme_bw() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14),
        panel.border = element_blank(),
        axis.line = element_line(color='black'))

for (elem in unique(ICPmelt[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
  print(elem)
  png(file.path(inspectdir, paste0('ICP_withinsite_',elem,'.png')), width = 20, height=12, units='in', res=300)
  print(
    ggplot(ICPmelt[Elem == elem,], 
           aes(x=SiteID, y = ICP, fill=SiteID)) + 
      geom_line(aes(group=SiteID), color='black') +
      geom_point(size=5, colour='black', pch=21, alpha=0.75) +
      labs(x='Element', y='Concentration') + 
      theme_bw() + 
      theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
            strip.text = element_text(size=14),
            panel.border = element_blank(),
            axis.line = element_line(color='black'))
  )
  dev.off()
}

# ---- 2. Univariate relationship to pollution predictors for the purpose of outlier detection ----
#Outlier diagnostic plots
pollutICPmerge[, pollutICP_flags := 0]
pollutICPmerge <- pollutICPmerge[!is.na(Elem),]
for (chem in unique(pollutICPmerge[!is.na(Elem), Elem])) {
  print(chem)
  pollutICP_lm <- lm(ICP ~  heatbing1902log300proj*heatsubAADTlog300*heatsubslopelog300*nlcd_imp_ps +
                       heatbustransitlog300 + heatbustransitlog300:heatsubslopelog300, 
                     data = pollutICPmerge[Elem == chem,])
  # ggsave(file.path(inspectdir, paste0('ICPpollution_regoutliers', chem, '.png')),
  #        chemregdiagnostic_custom(pollutICP_lm, pollutICPmerge, chem,  flagcol = 'pollutICP_flags', tresids = TRUE),
  #        width = 20, height=12, units='in', dpi=300)
  pollutICPmerge[Elem == chem, `:=`(pollutpred = fitted(pollutICP_lm),
                                    pollutICPR2 = summary(pollutICP_lm)$adj.r.squared)]
}

#Plot predicted ICP ~ observed ICP 
pollutICP_plot <- ggplot(pollutICPmerge[!(Elem %in% c(excol, excol2, NA)),], 
                         aes(x=pollutpred, y=ICP, group=1, color=factor(pollutICP_flags))) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  geom_smooth(method='lm') +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
pollutICP_plot

ggsave(file.path(inspectdir, 'ICPpollution_all.png'),
       pollutICP_plot,
       width = 20, height=12, units='in', dpi=600)

#Investigate determinants of R2
pollutmeanR2 <- unique(pollutICPmerge[, list(meanICP = mean(ICP, na.rm=T),
                                             R2 = pollutICPR2), by=Elem])[
                                               ICPthresholds_format, on='Elem'
                                             ]
ggplot(pollutmeanR2, aes(x=meanICP, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=meanICP/QUANTLIM, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

### ---- 3. Compile outlier flags ----
#Merge all flag datasets
ICP_formatflags <- merge(fieldXRF_format, ICPfieldmerge[, .(SiteID, Pair, Elem, ICPfield_flags, ICPfieldR2)], 
                         by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(ICPlabmerge[, .(SiteID, Pair, Elem, ICPlab_flags, ICPlabR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(pollutICPmerge[, .(SiteID, Pair, Elem, pollutICP_flags, pollutICPR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T)

ggplot(ICP_formatflags, aes(x=ICPfieldR2, y=pollutICPR2, label=Elem)) + 
  geom_text()
ggplot(ICP_formatflags, aes(x=ICPlabR2, y=pollutICPR2, label=Elem)) + 
  geom_text()

#Compute total number of flags for each tree across all elements that have an R2 of at least 0.10 across models
#Don't want to count outliers from spurious models
flagelems <- ICP_formatflags[, min(ICPfieldR2, ICPlabR2, pollutICPR2, na.rm=T)>0.20, by=Elem][
  V1==TRUE & Elem != 'Rh', Elem]
ICP_formatflags[Elem %in% flagelems,
                flagsum_ICP := sum(ICPfield_flags, ICPlab_flags, pollutICP_flags, na.rm=T), 
                by=.(SiteID, Pair)]

#Check total number of flags by each type of flag
ICP_formatflags[Elem %in% flagelems,
                `:=`(ICPfield_flags_sum = sum(ICPfield_flags, na.rm=T),
                     ICPlab_flags_sum = sum(ICPlab_flags, na.rm=T),
                     pollutICP_flags_sum = sum(pollutICP_flags, na.rm=T)),
                by=.(SiteID, Pair)]

ICP_formatflags_u <- unique(
  ICP_formatflags[!is.na(flagsum_ICP),
                  .SD, 
                  .SDcols=c('SiteID', 'Pair',
                            grep('flag.*sum', colnames(ICP_formatflags), value=T))])

ggplot(ICP_formatflags_u, aes(x=ICPfield_flags_sum, y=ICPlab_flags_sum, 
                              label=paste0(SiteID, Pair))) + 
  scale_color_distiller(palette='YlGnBu') +
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5) 

ggplot(ICP_formatflags_u, aes(x=ICPfield_flags_sum, y=pollutICP_flags_sum, 
                              label=paste0(SiteID, Pair))) + 
  scale_color_distiller(palette='YlGnBu') + 
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5)

#Records to inspect
flagelems
ICP_inspect <- unique(ICP_formatflags[(ICPfield_flags_sum>3) | 
                                        (ICPlab_flags_sum > 3) | 
                                        (pollutICP_flags_sum > 3), 
                                      .(SiteID, Pair)])
outliertrees <- trees[ICP_inspect, on=c('SiteID', 'Pair')]
outlierlocs <- SpatialPointsDataFrame(coords = data.frame(outliertrees$POINT_X, outliertrees$POINT_Y),
                                      data= as.data.frame(outliertrees))
#View(outlierlocs@data)
# leaflet(data = outlierlocs) %>% addTiles() %>%
#   addMarkers(clusterOptions = markerClusterOptions(),
#              popup = ~paste0(SiteID, Pair))

###################### ---- D. Compile all flags and make a flag matrix/heatmap then inspect data and decide on their fate ---- ####
# ---- 1. Compile flags and make heatmap ----
allflags <- ICP_formatflags_u[labXRF_formatflags_u, on=.(SiteID, Pair)][
  fieldXRF_formatflags_u, on=.(SiteID, Pair)] 
allflags[, (grep('i[.]', colnames(allflags), value=T)) := NULL][
  , flagpartsum := sum(flagsum_ICP, flagpartsum_lab, flagpartsum_field, na.rm=T), by=.(SiteID, Pair)][
    , SiteIDPair := factor(paste0(SiteID, Pair), levels = unique(paste0(SiteID, Pair)[order(-flagpartsum)]))]

colnames(allflags)
#levels(allflags_melt$variable)
allflags_melt <- melt(allflags, id.vars = c('SiteIDPair', 'SiteID', 'Pair')) %>%
  .[, variable := factor(gsub('[_]|flag*|sum', '', variable), 
                         levels=gsub('[_]|flag|sum', '',
                                     c("flagsum_lab", "flagsum_field", "flagpartsum",
                                       "flagsum_ICP" ,"flagpartsum_lab", "flagpartsum_field",  
                                       "pollutICP_flags_sum", "ICPlab_flags_sum","pollutlab_flags_sum",
                                       "CVflag_count_lab", "labfield_flags_sum", "ICPfield_flags_sum",
                                       "pollutfield_flags_sum", "CVflag_count_field")))] 

ggplot(data = allflags_melt[!(variable %in% c('field', 'lab', 'part'))],  #Re-add part if needed
       aes(x=SiteIDPair, y=variable, fill=value)) + 
  geom_tile() +
  scale_fill_distiller(palette='Spectral', trans = "sqrt") + 
  theme(legend.position = c(0.9, 0.2),
        legend.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90, vjust=-0.07))

# ---- 2. Analysis of the main outliers ----
#112A: South Beacon Hill, forested area across from Boeing Field/I5 (closest to road)
#     5 measurements on this tree. High within tree variability but no clear outlier point. 
#     pollut-field outlier for Cr, Fe, Ti, Zn, Zr (model overestimates pollution; could be due to protection from surrounding trees and topo.)

#111A: Interbay, Magnolia Bridge exit
#     High within tree variability with one apparent XRF sample with much higher value than the two others
#     pollut-field outlier for Cr, Fe, Ti, Zn (model underestimates pollution)

#114A and 114B: 15th Ave NE & 43rd St, UW campus by bus stop
#     Not abnormal within-tree variability for outlying elements, similar between 114A and 114B in both range and mean
#     pollut-field outlier for Cu, Fe, Zn (model underestimates pollution, particularly for Cu)

#117A and 117B: 1st Ave S, Southern end of bridge over Union Pacific cargo yard in Industrial District
#     Normal within-tree variability
#     pollut-field outlier for Cr, Fe, Zn, Zr (model underestimates pollution)

#49A: on slope by intersection SW of Garfield School in Central District
#     outlier for all lab things but not for field stuff. Must indicate that something was wrong with lab.
#     no collection or analysis comments, moss appearance NTR, overestimated for both ICP and lab XRF.
#     results completely different from 49B which is not an outlier, suggesting erroneous aspect.
#     Must be either the subsample of moss taken or the time between collecting and processing?
#     -> Remove from ICP and lab results

#23A: Harvard I5 exit towards 520 east
#     not really a problem for either field, lab, or ICP related to pollution drivers.
#     no collection or analysis comments, moss appearance NTR
#     ICP-pollution: on the higher pollution level so leverage but not truly outside of the cloud. Surely mostly due to slight heteroscedaticity
#     field-lab outlier for Br, Ca, Cu, Fe, Mn, Ti, Zr (higher for lab)
#     field-ICP outlier for Cu, Fe, Mn, Ni, Pb, Ti 
#     lab-ICP - not an outlier
#     normal within tree variability for lab and field
#     -> NTR keep it 

#51A: Downhill side of street on 15th Ave NE UW
#     multiple outlier flags for pollut-ICP and pollut-lab but nothing with fields. 
#     Must indicate that something was wrong with lab.
#     no collection or analysis comments, moss appearance NTR
#     pollut-field outlier for Ca, K, and Zn (but only because high value and variability at these levels)
#     pollut-lab outlier for Fe, Mn, Ti, Zn 
#     pollut-ICP outlier for Fe, Ni, Zn, but nothing striking
#     lab-ICP - not an outlier
#     field-lab and field-ICP - not an outlier
#     normal within tree variability for lab and field 
#     -> NTR keep it

#54A: Monroe Highway 2 and Main St intersection
#     not outlier for ICPlab and pollut-field
#     low outlier: pollut-ICP and pollut-lab, higher outlier: field-lab and field-ICP
#     suggests something was wrong with lab. moss was very wet when collected, could have rotten.
#     no collection or analysis comments, moss appearance NTR
#     pollut-lab outlier: Fe, Ti (quite far)
#     pollut-ICP outlier: Fe, Ti(but not completely outside)
#     field-lab and -ICP outlier: Fe, Mn, Zn
#     field within-tree variability: one high Fe and Zr value; otherwise nothing much
#     -> NTR keep it

#53A: Monroe - Lewis St and Highway 2 intersection
#     Fine with ICP, outlier with pollut-lab and pollut-field
#     no collection or analysis comments, moss appearance NTR
#     field-pollution outlier: Ca, K, Ni, Pb, 
#     lab-pollution outlier: Fe, Pb, Ti
#     ICP-pollution: None 
#     high within-tree and  for field and lab (for K, Pb, Ti, Si) 
#     high within-site variability: Ca, Co, Cr, Cu, K, Mn, Ni, Ti, Zn
#     -> could just keep 53B instead?

#19A: along I5 in Eastlake
#     high within-tree variability for field (Cu, Fe, Mg, Ti, Zn, Zr)
#     no collection or analysis comments, moss appearance NTR
#     lab XRF data are a bit off (must mean that a clump was taken that does not reflec the full tree)
#     It seems that one of the three XRF samples overshoots for multiple of the flagelems
#     ICP-lab outlier for Cu, Fe, Sr, Zn (not by huge amount)
#     -> keep it for now?

#62A: Sultan under bridge by boat launch over Skykomish
#     Site with a lot of Al, Cr, Fe, Nickel and Zirkonium - must be because was under bridge
#     no collection or analysis comments, moss appearance NTR
#     ICP-lab and ICP-fields quite off: Cr seems just that there might be a non-linear pattern, Fe, Ti, and Zn ICP also higher, 
#     pollutfields and labfields a little off but not much
#     It might be because it was under the bridge and receiving more heavy dust?
#     -> NTR keep it

#53B: Monroe - Lewis St and Highway 2 intersection
#     A little off for everything aside from lab-fields
#     no collection or analysis comments, moss appearance NTR
#     Not really outlier for pollut-fields, a little off for Zn and Fe pollut-ICP, a bit of an outlier for Zn pollut-lab
#     not much else
#     -> NTR keep it

#7A: Stone Ave North of Pacific Ave intersection downhill side
#     field measurement is off pollut-field (most elements)
#     lab-field and  ICP-field (Ca, Cu, Mn mostly) are off
#     no collection or analysis comments, moss appearance NTR
#     very different from 7B from Ba, Cr, Cu, Mn, Mo, Pb, Ti, Zr
#     normal within tree variability
#     -> remove from field analysis

#20B: horizontal tree in Eastlake below I5
#     SMall outlier for ICPfields, labfields, ICPlab, pollutlabs - seems an XRF issue
#     no collection or analysis comments, moss appearance NTR
#     High amounts of Pb, Ti, and Zr.. must be like 19A, low to the ground new I5
#     outlier for Pb for pollut-ICP, pollut-lab, field-lab (also for Sr), field-ICP, etc. just whacky
#     lab-ICP perfect for Pb, off for Cu
#     -> Keep it, high Pb (and some others) must lead to higher variability

#33B: Cap Hill by volunteer park
#     little outlier for pollut-lab and pollut-ICP, larger outlier for pollut-field
#     no collection or analysis comments, moss appearance NTR
#     large outlier across all pollution predictions for Cu
#     small pollut-field outlier for Zn
#     -> very similar value to 33A for Cu, must mean that it's pollution from other source? Keep 

#16A: By viewpoint in blackberry bush in Discovery Park
#     no collection or analysis comments, moss appearance NTR
#     off for pollut-field: higher amounts of Cr, Fe, Ni, Ti, and Zr than predicted (pretty within-tree variability in those elements as well). 
#     underpredicted amounts in lab and ICP for a few elements as well. Contamination could be sea-borne?
#     -> NTR: keep

#20A: same as 20B. 
#     no collection or analysis comments, moss appearance NTR
#     -> NTR

#44B: by Denny Way near Hotel 
#     Very dirty-sooty moss though no collection or analysis comments, moss appearance NTR
#     pollutICP (Fe, Mn), ICPfields and labfields (Fe, Zn - very unpredicted by labs/overpredicted by field)
#     -> NTR

#52B: Monroe intersection of 522 and highway 2
#     Very wet moss
#     a little within-tree variability for both XRF 
#     ICP detected high levels of Se, Cd, As, and Ni hence outlier for lab-ICP and field-ICP
#     not really outlier for pollut-field
#     -> keep

#52A: Monroe intersection of 522 and highway 2
#     very wet moss
#     medium outlier for pollut-field: Fe, Zn, Zr (under-predicted compared to observed/higher iron than predicted)
#     large difference between 52A and 52B, 
#     within-tree over- vs. underestimate vary by element in Fe (under), Zn (over), Zr (over)
#     -> Remove, keep only 52B for pollut-field

#43A: downtown by Amazon building
#     high outlier for pollut-field: Ca, Cu, K, Mn (over-predicted compared to observed)
#     no collection or analysis comments, moss appearance NTR
#     -> NTR

#1B:  calm neighborhood by lake city - first site
#     wrong moss species
#     medium outlier for pollut-field: both 1A and 1B are outliers for Ca, K, and Sr - must be because different species
#     but these elements don't have strong relationship anyways
#     -> NTR

#3B: Lake City Way
#     no comments
#     medium outlier: underpredicted Fe and Zr (but matches ICP and lab well for these elements)

#25A:
#     high within-tree CV for Cr - nothing to report

#34A: Capitol Hill 
#     extreme variability for Pb and outlier for pollutfield, not sure why.
#     remember there might have been some construction there?

#44A: near Denny way downtown
#     pollutICPS off - see 44B- underpredicted levels of Cr, Fe, Ti, Cu, and Zn


#Inspect within tree variability for field XRF data:
#19A: it seems that one of the three XRF samples overshoots for multiple of the flagelems
#20A: NTR
#23A: some extra variability for Zn, Fe, and Cu
#35A: NTR
#44B; NTR
#49A: NTR
#51A: NTR
#53A: one underestimates Zn, NTR otherwise
#54A: Fe one much higher, NTR otherwise
#62A: One a bit higher Fe, NTR otherwise
#7A: one lower Ca, Cu, one higher Fe 
#34A has one point with extra variability for Pb

# ---- 3. Removal of most egregious outliers and those with a paired tree ----
#Remove outliers for model exploration
pollutfieldclean <- pollutfieldmerge[!(paste0(SiteID, Pair) %in% 
                                         c('7A','52A', '114A', '114B', paste0(64:69, 'A'))),] #Remove 114A and 114B as well because they are just duplicate measurements of 51A at another time
pollutlabclean <- pollutlabmerge[!(paste0(SiteID, Pair) %in% c('53A', '49A')),]
pollutICPclean <-  pollutICPmerge[!(paste0(SiteID, Pair) %in% c('49A')),]


############################################################################################################################################


# 4. Check relationship among and reduce variables  -----------------
############################################################################################################################################
extraoutliers <- c('51A') #Remove anomalous data points by 15th Ave NE and Luwam's (latter do not appear to be precise enough

#---- A. Field XRF - all elements against each other ----
castformula <- as.formula(paste0(
  paste0(colnames(pollutfieldclean)[!(colnames(pollutfieldclean) %in%
                                        c('Elem', 'name', 'mean','range', 'sd', 'cv', 'transmean',
                                          'pollutfield_flags', 'pollutfieldR2', 'pollutpred'))],
         collapse='+'),
  '~Elem'))
pollutfieldclean_cast <- dcast(pollutfieldclean[!is.na(mean)], 
                               formula = castformula,
                               value.var= 'mean')
pollutfieldclean_cast[, SiteIDPair := paste0(SiteID, Pair)]

#Create vector of columns - 1.G. and  3.A.3. for column selection
"***NOTE: consider removing all elems with net/background photon count < 0.1"
elemcols <- colnames(pollutfieldclean_cast)[colnames(pollutfieldclean_cast) %in% periodicTable$symb] %>%
  setdiff('Rh')
elemcols_sub <- elemcols %>%
  setdiff(c(excol, excol2))

#Scatterplot matrix 
outcatmat <- file.path(inspectdir, 'FieldElem_FieldElem_matrix.png')
if (!file.exists(outcatmat)) {
  png(outcatmat , width = 24, height=24, units='in', res=300)
  ggscatmat(as.data.frame(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), elemcols, with=F]),
            alpha=0.7)
  dev.off()
}


#Correlation heatmap
outfieldheatmat <- file.path(inspectdir, 'corheatmap_FieldElem_FieldElem.png')
if (!file.exists(outfieldheatmat)) {
  png(outfieldheatmat, width = 8, height=8, units='in', res=300)
  corr_heatmap(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), elemcols, with=F])
  dev.off()
}

#PCA with rrcov (see 4.2 p24 of Todorov and Filzmoser 2009)
"Interestingly, PcaGrid, a robust PCA, gives completely different results than classic PCA with PC1 and PC2 
being clearly correlated"
pca <- PcaClassic(~., data=pollutfieldclean_cast[,elemcols_sub,with=F], scale=TRUE) #z-standardize
summary(pca)
screeplot(pca, main="Screeplot: classic PCA", bstick=TRUE) #First PC significant compared to broken stick
ordi.monte(pollutfieldclean_cast[,elemcols_sub,with=F],ord='pca',dim=5) #2PCs significant with Monte-carlo test of component significance
#Plot
#Plot
#Plot
#Plot
#Plot
biplot(pca, main="Robust biplot", col=c("gray55", "red"))  
plot(pca)

pollutfieldclean_pca<- cbind(pollutfieldclean_cast, pca$scores)
setnames(pollutfieldclean_pca, colnames(pca$scores), paste0(colnames(pca$scores), '_scores'))

#Create a PCA biplot matrix where all components are graphed against each other
pcabiplot_grid(pollutfieldclean_cast, nPCs = 5, cols = elemcols_sub,
               idcols = c('SiteID', 'Pair'),  scload = 3)

#Inspect components to decide which ones to predict
loadings(pca)
#There are no very distinct patterns - everything is pretty correlated:
#First component: Cr, Fe, Ti, Zn, Zr (0.34-0.36);
#                 Ca, Cu a little less (~0.30)
#                 Ni, and Pb a little less (0.20-25);
#                 Sr, Mn, K, Ba not much (0.13-0.16)
#Second component:K loads the strongest (0.53)
#                 Sr, Ni, Ca, Ba load moderately (0.30-0.37)
#                 Pb and Mn load positively but not much (0.16)
#                 Cu is ~ 0
#                 Cr and Ti load moderately and negatively (~-0.15)
#                 Zn, Zr, Fe load negatively < -20

#Select individual elements to predict separately: Zn, Cu, Pb

#Variable clustering
elemtree <- hclustvar(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), elemcols, with=F])
png(file.path(resdir, 'data_inspection/varclustree.png'), width = 8, height=4, units='in', res=300)
plot(elemtree)
dev.off()

# stab <- stability(elemtree,B=40)
# plot(stab, main="Stability of the partitions")
# stab$matCR
# boxplot(stab$matCR, main="Dispersion of the ajusted Rand index")
# P14<-cutreevar(elemtree, 14, matsim=TRUE)
# print(P14)
# P14$var

#---- B. Create a synthetic sum-based pollution index (averaging centered and standardized elements)   ----
#Check that all can be transformed using the same transformation then transform
selcols <- c('Cu', 'Pb', 'Zn') #In three main different groups in PCA
selcols_trans <- data.trans(as.data.frame(pollutfieldclean_cast[, elemcols[elemcols %in% selcols], with=F]), 
                            method='power',exp=1/3, plot=F)
cbcols <- paste0(elemcols[elemcols %in% selcols], 'cubrt')
pollutfieldclean_cast[, (cbcols) := selcols_trans]
#standardize by max value (use of transformed index actually leads to heteroscedacity)
pollutfieldclean_cast[, cbpollution_index := 100*Reduce('+', lapply(.SD, function(x) (x/max(x))))/length(cbcols), 
                      .SDcols = cbcols][
                        ,`:=`(pollution_index = 100*Reduce('+', lapply(.SD, function(x) (x/max(x))))/length(cbcols),
                              Znstand = 100*Zn/max(Zn),
                              Custand = 100*Cu/max(Cu),
                              Pbstand = 100*Pb/max(Pb)),
                        .SDcols = elemcols[elemcols %in% selcols]]

#---- C. All pollution drivers against each other ----
#Define columns to analyze
pollutcols <- grep('heat|NLCD', colnames(pollutfieldclean), value=T, ignore.case = T)

#Rescale all pollutcols
pollutfieldclean_cast[, (heatcols) := lapply(.SD, function(x) 100*x/max(x)), .SDcols = heatcols]

#Correlation heatmap
outpollutheatmat <- file.path(inspectdir, 'corheatmap_PollutionDrivers_PollutionDrivers.png')
if (!file.exists(outpollutheatmat)) {
  png(outpollutheatmat, width = 30, height=30, units='in', res=300)
  corr_heatmap(pollutfieldclean_cast[, pollutcols[-which(pollutcols=='NLCD_reclass_final_PS')], with=F]) + 
    scale_x_discrete(labels=heatlab) + 
    scale_y_discrete(labels=heatlab)
  dev.off()
}

#---- D. All pollution drivers against field XRF all elems ----
#Excluding 51A, 114A and 114B
outfieldpollutheatmat <- file.path(inspectdir, 'corheatmap_PollutionDrivers_FieldElem.png')
if (!file.exists(outfieldpollutheatmat)) {
  png(outfieldpollutheatmat,width = 40, height=35, units='in', res=300)
  corr_heatmap(xmat=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), pollutcols[-which(pollutcols=='NLCD_reclass_final_PS')], with=F],
               ymat=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), c(elemcols, 'pollution_index', 'cbpollution_index'), with=F],
               clus = FALSE) +
    scale_y_discrete(labels=heatlab) + 
    theme(text=element_text(size=22))
  dev.off()
}

xmat <- pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), pollutcols[-which(pollutcols=='NLCD_reclass_final_PS')], with=F]
ymat <- pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), c(elemcols, 'pollution_index', 'cbpollution_index'), with=F]
cordf <- round(cor(x=xmat, y=ymat),2) %>%
  t() %>%
  as.data.frame()
cordf$group <- rownames(cordf)
cordf <- cordf[, c(ncol(cordf),
                   1:(ncol(cordf)-1))]

spidermetals <- c('Fe', 'Cu', 'Zn','Pb')
spidercols <-  colnames(cordf)[!(grepl('(^heatAADT.*)|(log[0-9]{2}$)|.*mean.*', colnames(cordf))) &
                                 grepl('(.*log.*)|.*nlcd.*|group', colnames(cordf))]
spiderlabels <- gsub('(heat_?)|log', '',
                     gsub('bing[1-9]*', 'Congestion ',
                          gsub('subslope', 'Gradient ',
                               gsub('subAADT', 'Volume ',
                                    gsub('SPD', 'Speed ',
                                         gsub('bustransit', 'Transit ',
                                              gsub('nlcd_imp_ps.*', 'Imperv.',
                                                   spidercols)))))))
funcCircleCoords <- function(center = c(0,0), r = 1, npoints = 100){
  #Adapted from Joran's response to http://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
grid.mid = 0.5
grid.min = 0
grid.max = 1
gridline75 <- funcCircleCoords(
  c(0,0),
  0.75 +abs(grid.min - ((1/9)*(grid.max-grid.min))),
  npoints = 360)
gridline25 <- funcCircleCoords(
  c(0,0),
  0.25 +abs(grid.min - ((1/9)*(grid.max-grid.min))),
  npoints = 360)
plot.data <- cordf[cordf$group %in% spidermetals,
                   spidercols[order(spidercols)]]
plot.data <- as.data.frame(plot.data)
plot.data[,1] <- as.factor(as.character(plot.data[,1]))
names(plot.data)[1] <- "group"
var.names <- colnames(plot.data)[-1]  #'Short version of variable names
plot.data.offset <- plot.data
plot.data.offset[,2:ncol(plot.data)]<- plot.data[,2:ncol(plot.data)]+abs(centre.y)
#print(plot.data.offset)
# (b) convert into radial coords
CalculateGroupPath <- function(df) {
  #Converts variable values into a set of radial x-y coordinates
  #Code adapted from a solution posted by Tony M to
  #http://stackoverflow.com/questions/9614433/creating-radar-chart-a-k-a-star-plot-spider-plot-using-ggplot2-in-r
  #Args:
  #  df: Col 1 -  group ('unique' cluster / group ID of entity)
  #      Col 2-n:  v1.value to vn.value - values (e.g. group/cluser mean or median) of variables v1 to v.n
  
  path <- df[,1]
  
  ##find increment
  angles = seq(from=0, to=2*pi, by=(2*pi)/(ncol(df)-1))
  ##create graph data frame
  graphData= data.frame(seg="", x=0,y=0)
  graphData=graphData[-1,]
  
  for(i in levels(path)){
    pathData = subset(df, df[,1]==i)
    for(j in c(2:ncol(df))){
      #pathData[,j]= pathData[,j]
      
      
      graphData=rbind(graphData, data.frame(group=i,
                                            x=pathData[,j]*sin(angles[j-1]),
                                            y=pathData[,j]*cos(angles[j-1])))
    }
    ##complete the path by repeating first pair of coords in the path
    graphData=rbind(graphData, data.frame(group=i,
                                          x=pathData[,2]*sin(angles[1]),
                                          y=pathData[,2]*cos(angles[1])))
  }
  #Make sure that name of first column matches that of input data (in case !="group")
  colnames(graphData)[1] <- colnames(df)[1]
  graphData #data frame returned by function
}
group <-NULL
group$path <- CalculateGroupPath(plot.data.offset)

#Create palette
library(RColorBrewer)
interspec <- colorRampPalette(brewer.pal(11, 'Spectral'))
intdat<- data.frame(y=round(group$path$y,2))
datrange <- seq(min(intdat), max(intdat), 0.01)
coldf <- data.frame(y = datrange, ycol=interspec(length(datrange)))
colvalues <- merge(intdat, coldf, on='y', all.y=F)
colvalues$yfac <- factor(as.character(colvalues$y),
                         levels= unique(as.character(colvalues$y)[order(-colvalues$y)]))

ggplot(colvalues, aes(x=y, y=y, color=yfac)) +
  geom_point() +
  scale_color_manual(values=unique(as.character(colvalues$ycol)))

ggradarplot <- ggradar(plot.data=cordf[cordf$group %in% c('Zn', 'Cu', 'Pb'),
                                       spidercols[order(spidercols)]],
                       grid.min=grid.min, grid.max=grid.max, grid.mid = grid.mid,
                       gridline.mid.linetype = "solid",
                       gridline.mid.colour= 'darkgrey',
                       grid.label.size = 0,
                       group.point.size = 2,
                       group.line.width = 1.1,
                       axis.labels = str_wrap(spiderlabels[order(spidercols)][-1], 5),
                       axis.label.size = 8,
                       axis.label.offset = 1.1,
                       #group.colours = '#8CE071',
                       background.circle.transparency=0.15) +
  geom_path(data=gridline75, aes(x=x, y=y), color='darkgrey') +
  geom_path(data=gridline25, aes(x=x, y=y), color='darkgrey') +
  # scale_colour_manual(values=unique(as.character(colvalues$ycol))) +
  # facet_wrap(~group) +
  theme_minimal() +
  theme(text= element_text(size = 18),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position=c(0.15, 0.15),
        plot.margin=unit(c(-1.5, -1.5, -1.5, -1.5),'cm'))
ggradarplot


extrafont::loadfonts()
pdf(file.path(resdir, 'data_inspection/spiderplot_20191219.pdf'), width=9, height=9)
#png(file.path(moddir, 'spiderplotZn_20190514_1.png'), width=9, height=9, units='in', res=450)
ggradarplot
dev.off()


"Summary of best univariate predictors for each category
- Zn: 
  - IMPERVIOUSNESS > BING > SPDL > AADT > SLOPE  > TRANSIT
  - AADT: 100 > 200-300 log~pw 
  - SPDL: 50~100~200>300>500
  - slope: 500>300>200>100 and log~pow
  - bing: 300=200>100>500 log~pow
  - transit: 500=300=200 > 100 log~pow
  - imp_mean > imp
"

############################################# Write out data frame for modelling steps ###############
write.fst(pollutfieldclean_cast, 
          file.path(resdir, 
                    paste0('pollutfieldclean_cast', Sys.Date(), '.fst')))
##############################################################################


# 5. Validate field XRF and lab XRF against ICP-OES -------------------
###################### ---- E. Make summary table of variability ----------------
#Find sites that are available for all three methods
commonsites <- intersect(
  intersect(ICPmean[,unique(paste0(SiteID, Pair))], 
            fieldXRF_format[,unique(paste0(SiteID, Pair))]
            ),
  labXRF_format[,unique(paste0(SiteID, Pair))]
)

commonduplisites <- intersect(
  intersect(ICPmean[duplicated(SiteID), unique(SiteID)], 
            fieldXRF_format[duplicated(SiteID), unique(SiteID)]
  ),
  labXRF_format[duplicated(SiteID), unique(SiteID)]
)


treefielddt <- fieldXRF_format[(!(Elem %in% c('Rh','Pd','Ar'))) & 
                                 (paste0(SiteID, Pair) %in% commonsites), 
                               list(treeCV_field = round(100*mean(cv, na.rm=T))), by=Elem] %>%
  setkey(Elem)
sitefielddt <- fieldXRF_format[(!(Elem %in% c('Rh','Pd','Ar'))) & 
                  (SiteID %in% commonduplisites),
                list(cvsite = sd(mean, na.rm=T)/mean(mean, na.rm=T)), by=.(SiteID, Elem)] %>%
  .[, list(siteCV_field = round(100*mean(cvsite, na.rm=T))), by=.(Elem)] %>%
  setkey(Elem)

pelletlabdt <- labXRF_format[!(Elem %in% c('Rh','Pd','Ar')) & 
                (paste0(SiteID, Pair) %in% commonsites), 
              list(pelletCV_lab = round(100*mean(cv, na.rm=T))), by=Elem] %>%
  setkey(Elem)
sitelabdt <-labXRF_format[(!(Elem %in% c('Rh','Pd','Ar'))) & 
                  (SiteID %in% commonduplisites),
                list(cvsite = sd(mean, na.rm=T)/mean(mean, na.rm=T)), by=.(SiteID, Elem)] %>%
  .[, list(siteCV_lab = round(100*mean(cvsite, na.rm=T))), by=.(Elem)] %>%
  setkey(Elem)

siteicpdt <- melt(ICPmean, id.vars = c('SAMPLE.SET', 'SiteID', 'Pair'),
                  variable.name = 'Elem')[
  (!(Elem %in% c('Rh','Pd'))) & 
    (SiteID %in% commonduplisites),
  list(cvsite = sd(value, na.rm=T)/mean(value, na.rm=T))
  , by=.(SiteID, Elem)] %>%
  .[, list(siteCV_ICP = round(100*mean(cvsite, na.rm=T))
           ), by=.(Elem)] %>%
  setkey(Elem)


icpfracstat <- melt(ICPmean, id.vars = c('SAMPLE.SET', 'SiteID', 'Pair'),
    variable.name = 'Elem')[
      (!(Elem %in% c('Rh','Pd'))) & 
        (paste0(SiteID, Pair) %in% commonsites),
      list(
        icpstat = paste0(round(mean(value, na.rm=T)), ' (', 
                       round(min(value, na.rm=T)),'-',
                       round(max(value, na.rm=T)), ')')
        ), by=Elem]
      

cvsummary_dt <- merge(treefielddt, pelletlabdt, all.x=T, all.y=T) %>%
  merge(ICPCV, all.x=T, all.y=T) %>%
  merge(sitefielddt, all.x=T, all.y=T) %>%
  merge(sitelabdt, all.x=T, all.y=T) %>%
  merge(siteicpdt, all.x=T, all.y=T) %>%
  merge(icpfracstat, all.x=T, all.y=T) %>%
  .[!(Elem %in% c('Ag')),] %>%
  setorder(Elem)

kable(cvsummary_dt, 
      format = "html", escape = F) %>%
  kable_styling("striped", full_width = T) %>%
  save_kable(file.path(resdir, 'CVsummary_20201028.doc'), self_contained=T)


###################### ---- F. Check relationship between methods ---------------
plot_elemvalid <- function(dt, elem, x, y, color, method='gam', outliers=NULL) {
  if (!is.null(outliers)) {
    dtout <- dt[(SiteIDPair %in% outliers),]
    dt <- dt[!(SiteIDPair %in% outliers),]
  }
  
  p <- ggplot(dt[Elem==elem,], aes_string(x=x, y=y)) + 
    geom_smooth(method=method, color='black') +
    geom_point(color=color, size=3, alpha=1/2) +
    scale_color_distiller(palette="Spectral") + 
    theme_classic()
  
  if (!is.null(outliers)) {
    p <- p + 
      geom_point(data=dtout[Elem==elem,], color='black', size=3, alpha=1/2)
  }
  return(p)
}

# ---- Relationships between field XRF and lab XRF ---- 
labfieldsub <- labfieldmerge[!(Elem %in% c(excol, excol2)) & !is.na(lab_mean),] %>%
  .[, drydelay := as.numeric(
    difftime(as.Date(substr(lab_Dry.start, 1, 10), format="%m/%d/%Y"),
             Date)
    )] %>%
  .[, SiteIDPair := paste0(SiteID, Pair)] %>%
  .[, `:=`(loglab = log10(lab_mean),
           logfield = log10(field_mean))]


qplot(labfieldsub$drydelay)

#Get correlation with outliers
labfieldsub[, list(r_all=cor(x=field_mean, y=lab_mean)), by=Elem]
labfieldsub[, list(r_all=cor(x=field_transmean, y=lab_transmean)), by=Elem]

# labfieldsub<- labfieldsub[, list(lab_mean=mean(lab_mean), field_mean=mean(field_mean)), 
#                           by=.(SiteID, Elem)]
# plot_elemvalid(dt=labfieldsub, elem='Cu', x='field_mean', y='lab_mean',
#                color='#238b45') +
#   geom_text(aes(label=SiteID))

# -------- Cu -------------
#Check relationship
plot_elemvalid(dt=labfieldsub, elem='Cu', x='field_mean', y='lab_mean',
               color='#238b45') +
  geom_text(aes(label=SiteIDPair))


labfieldmodCu_prelim <- glm(lab_mean~field_mean,
                   data = labfieldsub[Elem=='Cu',],
                   family=Gamma('log'))
GAMrescheck(labfieldmodCu_prelim)
Cuprelimdiag <- regdiagnostic_customtab(mod=labfieldmodCu_prelim,kCV = TRUE, 
                                        k=10, cvreps=50,
                                        labelvec = labfieldsub[Elem=='Cu', SiteIDPair],
                                        remove_outliers = 'outliers')
Cuoutliers_fieldlab <- strsplit(gsub('\\\\', '', Cuprelimdiag['outliers']), ',')$outliers

#Check relationship again
plot_elemvalid(dt=labfieldsub, outliers=Cuoutliers_fieldlab, 
               elem='Cu', x='field_mean', y='lab_mean',
               color='#238b45', method='lm') +
  scale_x_log10() + 
  scale_y_log10()

#GLM
labfieldmodCu <- list()
labfieldmodCu[[1]] <- glm(lab_mean~field_mean,
                          data = labfieldsub[Elem=='Cu',], 
                          family = gaussian('log'))
GAMrescheck(labfieldmodCu[[1]])

labfieldmodCu[[1]] <- glm(lab_mean~field_mean,
                  data = labfieldsub[Elem=='Cu',], 
                  family = Gamma('log'))
GAMrescheck(labfieldmodCu[[1]])

labfieldmodCu[[2]] <- glm(lab_mean~field_mean + drydelay,
                  data = labfieldsub[Elem=='Cu',], 
                  family = Gamma('log'))
GAMrescheck(labfieldmodCu[[2]])

labfieldmodCu[[3]] <- glm(lab_mean~logfield,
                          data = labfieldsub[Elem=='Cu',], 
                          family = Gamma('log'))
GAMrescheck(labfieldmodCu[[3]])


labfieldmodCu[[4]] <- glm(lab_mean~field_mean,
                          data = labfieldsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_fieldlab),], 
                          family = Gamma('log'))
GAMrescheck(labfieldmodCu[[4]])

labfieldmodCu[[5]] <- glm(lab_mean~logfield,
                          data = labfieldsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_fieldlab),], 
                          family = Gamma('log'))
GAMrescheck(labfieldmodCu[[5]])


labfieldmodCu[[3]] <- mgcv::gam(lab_mean~s(field_mean, k=4),
                        data = labfieldsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_fieldlab),], 
                        family = gaussian('log'))
GAMrescheck(labfieldmodCu[[3]])
GAMmultiplot(labfieldmodCu[[3]])

# -------- Pb -------------
#Check relationship
plot_elemvalid(dt=labfieldsub, elem='Pb', x='field_mean', y='lab_mean',
               color='#238b45', method='lm') + 
  scale_x_log10() + 
  geom_text(aes(label=SiteIDPair)) +
  scale_y_log10()

labfieldmodPb_prelim <-  glm(lab_mean~logfield,
                             data = labfieldsub[Elem=='Pb',], 
                             family = Gamma('log'))

Pbprelimdiag <- regdiagnostic_customtab(mod=labfieldmodPb_prelim,
                                        kCV = TRUE, 
                                        k=10, cvreps=50,
                                        labelvec = labfieldsub[Elem=='Pb', SiteIDPair],
                                        remove_outliers = 'outliers & leverage')
Pboutliers_fieldlab <- strsplit(gsub('\\\\', '', Pbprelimdiag['outliers']), ',')$outliers

#Check relationship again
plot_elemvalid(dt=labfieldsub, outliers=Pboutliers_fieldlab, 
               elem='Pb', x='field_mean', y='lab_mean',
               color='#238b45', method='lm') +
  scale_x_log10() + 
  scale_y_log10()

#GLM
labfieldmodPb <- list()
labfieldmodPb[[1]] <- glm(lab_mean~field_mean,
                  data = labfieldsub[Elem=='Pb',], 
                  family = Gamma('log'))
GAMrescheck(labfieldmodPb[[1]])

labfieldmodPb[[2]] <- glm(lab_mean~field_mean,
                  data = labfieldsub[Elem=='Pb' & !(SiteIDPair %in% Pboutliers_fieldlab),], 
                  family = Gamma('log'))
GAMrescheck(labfieldmodPb[[2]])

labfieldmodPb[[3]] <- glm(lab_mean~logfield,
                  data = labfieldsub[Elem=='Pb',], 
                  family = Gamma('log'))
GAMrescheck(labfieldmodPb[[3]])

labfieldmodPb[[4]] <- glm(lab_mean~logfield,
                          data = labfieldsub[Elem=='Pb' & !(SiteIDPair %in% Pboutliers_fieldlab),], 
                          family = Gamma('log'))
GAMrescheck(labfieldmodPb[[4]])

labfieldmodPb[[5]] <- mgcv::gam(lab_mean~s(field_mean, k=3),
                        data = labfieldsub[Elem=='Pb' & !(SiteIDPair %in% Pboutliers_fieldlab),], 
                        family = gaussian('log'))
GAMrescheck(labfieldmodPb[[5]])
GAMmultiplot(labfieldmodPb[[5]])


# -------- Zn -------------
#Check relationship
plot_elemvalid(dt=labfieldsub, elem='Zn', x='field_mean', y='lab_mean',
               color='#238b45')

labfieldmodZn_prelim <- glm(lab_mean~logfield,
                    data = labfieldsub[Elem=='Zn',],
                    family=Gamma('log'))
GAMrescheck(labfieldmodZn_prelim)
Znprelimdiag <- regdiagnostic_customtab(mod=labfieldmodZn_prelim,
                                        kCV = TRUE, 
                                        k=10, cvreps=50,
                                        labelvec = labfieldsub[Elem=='Zn', SiteIDPair],
                                        remove_outliers = 'outliers')
Znoutliers_fieldlab <- strsplit(gsub('\\\\', '', Znprelimdiag['outliers']), ',')$outliers

#Check relationship again
plot_elemvalid(dt=labfieldsub, outliers=Znoutliers_fieldlab, 
               elem='Zn', x='field_mean', y='lab_mean',
               color='#238b45', method='lm') +
  scale_x_log10() + 
  scale_y_log10()

#GLM
labfieldmodZn <- list()
labfieldmodZn[[1]] <- glm(lab_mean~field_mean,
                  data = labfieldsub[Elem=='Zn',], 
                  family = Gamma('log'))
GAMrescheck(labfieldmodZn[[1]])

labfieldmodZn[[2]] <- glm(lab_mean~field_mean,
                  data = labfieldsub[Elem=='Zn' & !(SiteIDPair %in% Znoutliers_fieldlab),], 
                  family = gaussian('log'))
GAMrescheck(labfieldmodZn[[2]])

labfieldmodZn[[3]] <- glm(lab_mean~field_mean + drydelay,
                  data = labfieldsub[Elem=='Zn',], 
                  family = Gamma('log'))
GAMrescheck(labfieldmodZn[[3]])

labfieldmodZn[[4]] <- mgcv::gam(lab_mean~s(field_mean, k=3),
                        data = labfieldsub[Elem=='Zn' & !(SiteIDPair %in% Znoutliers_fieldlab),], 
                        family = Gamma('log'))
GAMrescheck(labfieldmodZn[[4]])
GAMmultiplot(labfieldmodZn[[4]])

labfieldmodZn[[5]] <- lm(loglab~logfield,
                  data = labfieldsub[Elem=='Zn' & !(SiteIDPair %in% Znoutliers_fieldlab),])
ols_regress(labfieldmodZn[[5]])
#ols_plot_diagnostics(labfieldmodZn[[5]])




# ---- Relationships between field XRF and ICP OES ------
ICPfieldsub <- ICPfieldmerge[!(Elem %in% c(excol, excol2)) & !is.na(ICP),]  %>%
  #.[Elem=='Pb' & ICP == 0, ICP := 0.013] %>%
  .[, `:=`(SiteIDPair = paste0(SiteID, Pair),
           field_mean = mean)] %>%
  .[, `:=`(logICP = log10(ICP),
           logfield = log10(field_mean))]

#Get correlation with outliers
ICPfieldsub[, list(r_all=cor(x=field_mean, y=ICP)), by=Elem]

# -------- Cu -------------
#Check relationship
plot_elemvalid(dt=ICPfieldsub, elem='Cu', x='field_mean', y='ICP', color='#238b45')

ICPfieldmodCu_prelim <- glm(ICP~logfield,
                    data = ICPfieldsub[Elem=='Cu',],
                    family=Gamma('log'))
GAMrescheck(ICPfieldmodCu_prelim)
Cuprelimdiag <- regdiagnostic_customtab(mod=ICPfieldmodCu_prelim,
                                        kCV = TRUE, 
                                        k=10, cvreps=50,
                                        labelvec = ICPfieldsub[Elem=='Cu', SiteIDPair],
                                        remove_outliers = 'outliers')
Cuoutliers_fieldICP <- strsplit(gsub('\\\\', '', Cuprelimdiag['outliers']), ',')$outliers

#Check relationship again
plot_elemvalid(dt=ICPfieldsub, outliers=Cuoutliers_fieldICP, 
               elem='Cu', x='field_mean', y='ICP',
               color='#238b45', method='lm') +
  scale_x_log10() + 
  scale_y_log10()

#GLM
ICPfieldmodCu <- list()
ICPfieldmodCu[[1]] <- glm(ICP~logfield,
                  data = ICPfieldsub[Elem=='Cu',], 
                  family = Gamma('log'))
GAMrescheck(ICPfieldmodCu[[1]])

ICPfieldmodCu[[2]] <- glm(ICP~logfield,
                  data = ICPfieldsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_fieldICP),], 
                  family = Gamma('log'))
GAMrescheck(ICPfieldmodCu[[2]])

ICPfieldmodCu[[3]] <- mgcv::gam(ICP~s(field_mean, k=3),
                        data = ICPfieldsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_fieldICP),], 
                        family = Gamma('log'))
GAMrescheck(ICPfieldmodCu[[3]])
GAMmultiplot(ICPfieldmodCu[[3]])

ICPfieldmodCu[[5]] <- lm(logICP~logfield,
                 data = ICPfieldsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_fieldICP),])
ols_regress(ICPfieldmodCu[[5]])
#ols_plot_diagnostics(ICPfieldmodCu[[5]])

#Choose mod2

plot_elemvalid(dt=ICPfieldsub, outliers=Cuoutliers_fieldICP, 
               elem='Cu', x='field_mean', y='ICP',
               color='Date', method='none') +
  geom_line(data=cbind(ICPfieldsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_fieldICP),],
                       predict(ICPfieldmodCu[[2]], type='response')), aes(y=V2)) +
  scale_x_log10() + 
  scale_y_log10()


# -------- Pb -------------
#Check relationship
plot_elemvalid(dt=ICPfieldsub, elem='Pb', x='field_mean', y='ICP', color='#238b45')

ICPfieldmodPb_prelim <- glm(ICP~logfield,
                            data = ICPfieldsub[Elem=='Pb' & ICP > 0,],
                            family=Gamma('log'))
GAMrescheck(ICPfieldmodPb_prelim)
Pbprelimdiag <- regdiagnostic_customtab(mod=ICPfieldmodPb_prelim,
                                        kCV = TRUE, 
                                        k=10, cvreps=50,
                                        labelvec = ICPfieldsub[Elem=='Pb', SiteIDPair],
                                        remove_outliers = 'outliers')
Pboutliers_fieldICP <- strsplit(gsub('\\\\', '', Pbprelimdiag['outliers']), ',')$outliers

#Check relationship again
plot_elemvalid(dt=ICPfieldsub, outliers=Pboutliers_fieldICP, 
               elem='Pb', x='field_mean', y='ICP',
               color='#238b45', method='lm') +
  scale_x_log10() + 
  scale_y_log10()

#GLM
ICPfieldmodPb <- list()
ICPfieldmodPb[[1]] <- glm(ICP~logfield,
                          data = ICPfieldsub[Elem=='Pb' & ICP > 0,], 
                          family = Gamma('log'))
GAMrescheck(ICPfieldmodPb[[1]])

ICPfieldmodPb[[2]] <- glm(ICP~logfield,
                          data = ICPfieldsub[Elem=='Pb' & 
                                               !(SiteIDPair %in% Pboutliers_fieldICP) &
                                               ICP >0,], 
                          family = Gamma('log'))
GAMrescheck(ICPfieldmodPb[[2]])

ICPfieldmodPb[[3]] <- mgcv::gam(ICP~s(field_mean, k=3),
                                data = ICPfieldsub[Elem=='Pb' & !(SiteIDPair %in% Pboutliers_fieldICP) & ICP >0,], 
                                family = Gamma('log'))
GAMrescheck(ICPfieldmodPb[[3]])
GAMmultiplot(ICPfieldmodPb[[3]])

ICPfieldmodPb[[5]] <- lm(logICP~logfield,
                         data = ICPfieldsub[Elem=='Pb' & !(SiteIDPair %in% Pboutliers_fieldICP) & ICP >0,])
ols_regress(ICPfieldmodPb[[5]])

#Choose ICPfieldmodPb[[2]]

# -------- Zn -------------
#Check relationship
plot_elemvalid(dt=ICPfieldsub, elem='Zn', x='field_mean', y='ICP', color='#238b45')

ICPfieldmodZn_prelim <- glm(ICP~logfield,
                            data = ICPfieldsub[Elem=='Zn',],
                            family=Gamma('log'))
GAMrescheck(ICPfieldmodZn_prelim)
Znprelimdiag <- regdiagnostic_customtab(mod=ICPfieldmodZn_prelim,
                                        kCV = TRUE, 
                                        k=10, cvreps=50,
                                        labelvec = ICPfieldsub[Elem=='Zn', SiteIDPair],
                                        remove_outliers = 'outliers')
Znoutliers_fieldICP <- strsplit(gsub('\\\\', '', Znprelimdiag['outliers']), ',')$outliers 

#Check relationship again
plot_elemvalid(dt=ICPfieldsub, outliers=Znoutliers_fieldICP, 
               elem='Zn', x='field_mean', y='ICP',
               color='#238b45', method='lm') +
  scale_x_log10() + 
  scale_y_log10()

#GLM
ICPfieldmodZn <- list()
ICPfieldmodZn[[1]] <- glm(ICP~logfield,
                          data = ICPfieldsub[Elem=='Zn',], 
                          family = Gamma('log'))
GAMrescheck(ICPfieldmodZn[[1]])

ICPfieldmodZn[[2]] <- glm(ICP~logfield,
                          data = ICPfieldsub[Elem=='Zn' & !(SiteIDPair %in% Znoutliers_fieldICP),], 
                          family = Gamma('log'))
GAMrescheck(ICPfieldmodZn[[2]])

ICPfieldmodZn[[3]] <- mgcv::gam(ICP~s(field_mean, k=3),
                                data = ICPfieldsub[Elem=='Zn' & !(SiteIDPair %in% Znoutliers_fieldICP),], 
                                family = Gamma('log'))
GAMrescheck(ICPfieldmodZn[[3]])
GAMmultiplot(ICPfieldmodZn[[3]])

ICPfieldmodZn[[5]] <- lm(logICP~logfield,
                         data = ICPfieldsub[Elem=='Zn' & !(SiteIDPair %in% Znoutliers_fieldICP),])
ols_regress(ICPfieldmodZn[[5]])

#Choose ICPfieldmodZn[[2]]



# ---- Relationships between lab XRF and ICP OES ------
ICPlabsub <- ICPlabmerge[!(Elem %in% c(excol, excol2)) & !is.na(ICP),]  %>%
  #.[Elem=='Pb' & ICP == 0, ICP := 0.013] %>%
  .[, `:=`(SiteIDPair = paste0(SiteID, Pair),
           lab_mean = mean)] %>%
  .[, `:=`(logICP = log10(ICP),
           loglab = log10(lab_mean))]

#Get correlation with outliers
ICPlabsub[, list(r_all=cor(x=lab_mean, y=ICP)), by=Elem]

# -------- Cu -------------
#Check relationship
plot_elemvalid(dt=ICPlabsub, elem='Cu', x='lab_mean', y='ICP', color='#d94801')

ICPlabmodCu_prelim <- glm(ICP~loglab,
                            data = ICPlabsub[Elem=='Cu',],
                            family=Gamma('log'))
GAMrescheck(ICPlabmodCu_prelim)
Cuprelimdiag <- regdiagnostic_customtab(mod=ICPlabmodCu_prelim,
                                        kCV = TRUE, 
                                        k=5, cvreps=50,
                                        labelvec = ICPlabsub[Elem=='Cu', SiteIDPair],
                                        remove_outliers = 'outliers & leverage')
Cuoutliers_ICPlab <- strsplit(gsub('\\\\', '', Cuprelimdiag['outliers']), ',')$outliers

#Check relationship again
plot_elemvalid(dt=ICPlabsub, outliers=Cuoutliers_ICPlab, 
               elem='Cu', x='lab_mean', y='ICP',
               color='#d94801', method='gam') +
  scale_x_log10() + 
  scale_y_log10()

#GLM
ICPlabmodCu <- list()
ICPlabmodCu[[1]] <- glm(ICP~lab_mean,
                          data = ICPlabsub[Elem=='Cu',], 
                          family = Gamma('log'))
GAMrescheck(ICPlabmodCu[[1]])

ICPlabmodCu[[2]] <- glm(ICP~loglab + I(loglab^2),
                          data = ICPlabsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_ICPlab),], 
                          family = Gamma('log'))
GAMrescheck(ICPlabmodCu[[2]])

glm(ICP~poly(loglab, 2),
    data = ICPlabsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_ICPlab),], 
    family = Gamma('log'))

ICPlabmodCu[[3]] <- mgcv::gam(ICP~s(loglab, k=2),
                                data = ICPlabsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_ICPlab),], 
                                family = Gamma('log'))
GAMrescheck(ICPlabmodCu[[3]])
GAMmultiplot(ICPlabmodCu[[3]])

ICPlabmodCu[[5]] <- lm(logICP~loglab,
                         data = ICPlabsub[Elem=='Cu' & !(SiteIDPair %in% Cuoutliers_ICPlab),])
ols_regress(ICPlabmodCu[[5]])
#ols_plot_diagnostics(ICPlabmodCu[[5]])

#Choose mod2

plot_elemvalid(dt=ICPlabsub, outliers=Cuoutliers_ICPlab, 
               elem='Cu', x='lab_mean', y='ICP',
               color='#d94801', method='none') +
  geom_line(data=cbind(ICPlabsub[Elem=='Cu',],
                       predict(ICPlabmodCu[[2]], type='response')), aes(y=V2)) +
  scale_x_log10() + 
  scale_y_log10()


# -------- Pb -------------
#Check relationship
plot_elemvalid(dt=ICPlabsub, elem='Pb', x='lab_mean', y='ICP', color='#d94801')

ICPlabmodPb_prelim <- glm(ICP~loglab,
                            data = ICPlabsub[Elem=='Pb' & ICP > 0,],
                            family=gaussian('log'))
GAMrescheck(ICPlabmodPb_prelim)
Pbprelimdiag <- regdiagnostic_customtab(mod=ICPlabmodPb_prelim,
                                        kCV = TRUE, 
                                        k=5, cvreps=50,
                                        labelvec = ICPlabsub[Elem=='Pb', SiteIDPair],
                                        remove_outliers = 'outliers & leverage')
Pboutliers_ICPlab <- strsplit(gsub('\\\\', '', Pbprelimdiag['outliers']), ',')$outliers 

#Check relationship again
plot_elemvalid(dt=ICPlabsub, outliers=Pboutliers_ICPlab, 
               elem='Pb', x='lab_mean', y='ICP',
               color='#d94801', method='lm') +
  scale_x_log10() + 
  scale_y_log10()

#GLM
ICPlabmodPb <- list()
ICPlabmodPb[[2]] <- glm(ICP~loglab,
                          data = ICPlabsub[Elem=='Pb' & ICP >0,], 
                          family = gaussian('log'))
GAMrescheck(ICPlabmodPb[[2]])

ICPlabmodPb[[3]] <- mgcv::gam(ICP~s(lab_mean, k=3),
                                data = ICPlabsub[Elem=='Pb' &
                                                   !(SiteIDPair %in% Pboutliers_ICPlab) &
                                                   ICP >0,], 
                                family = Gamma('log'))
GAMrescheck(ICPlabmodPb[[3]])
GAMmultiplot(ICPlabmodPb[[3]])

ICPlabmodPb[[5]] <- lm(logICP~loglab,
                         data = ICPlabsub[Elem=='Pb' & !(SiteIDPair %in% Pboutliers_ICPlab) & ICP >0,])
ols_regress(ICPlabmodPb[[5]])

#Choose ICPlabmodPb[[2]]

# -------- Zn -------------
#Check relationship
plot_elemvalid(dt=ICPlabsub, elem='Zn', x='lab_mean', y='ICP', color='#d94801')

ICPlabmodZn_prelim <- glm(ICP~loglab,
                            data = ICPlabsub[Elem=='Zn',],
                            family=Gamma('log'))
GAMrescheck(ICPlabmodZn_prelim)
Znprelimdiag <- regdiagnostic_customtab(mod=ICPlabmodZn_prelim,
                                        kCV = TRUE, 
                                        k=5, cvreps=50,
                                        labelvec = ICPlabsub[Elem=='Zn', SiteIDPair],
                                        remove_outliers = 'outliers & leverage')
Znoutliers_ICPlab <- strsplit(gsub('\\\\', '', Znprelimdiag['outliers']), ',')$outliers 

#Check relationship again
plot_elemvalid(dt=ICPlabsub, outliers=Znoutliers_ICPlab, 
               elem='Zn', x='lab_mean', y='ICP',
               color='#d94801', method='lm') +
  scale_x_log10() + 
  scale_y_log10()

#GLM
ICPlabmodZn <- list()
ICPlabmodZn[[2]] <- glm(ICP~loglab,
                          data = ICPlabsub[Elem=='Zn',], 
                          family = Gamma('log'))
GAMrescheck(ICPlabmodZn[[2]])

ICPlabmodZn[[3]] <- mgcv::gam(ICP~s(lab_mean, k=3),
                                data = ICPlabsub[Elem=='Zn' & !(SiteIDPair %in% Znoutliers_ICPlab),], 
                                family = Gamma('log'))
GAMrescheck(ICPlabmodZn[[3]])
GAMmultiplot(ICPlabmodZn[[3]])

ICPlabmodZn[[5]] <- lm(logICP~loglab,
                         data = ICPlabsub[Elem=='Zn' & !(SiteIDPair %in% Znoutliers_ICPlab),])
ols_regress(ICPlabmodZn[[5]])

#Choose ICPlabmodZn[[2]]

# ---- Compilation of models --------------------------------------------------
#--------- Make table ---------------
final_modlist <- list(
  ICPfieldmodCu[[1]],
  ICPfieldmodCu[[2]],
  ICPfieldmodPb[[1]],
  ICPfieldmodPb[[2]],
  ICPfieldmodZn[[1]],
  ICPfieldmodZn[[2]],
  ICPlabmodCu[[2]],
  ICPlabmodPb[[2]],
  ICPlabmodZn[[2]]
)

lapply(final_modlist, function(mod) GAMrescheck(mod))

final_modsummary <- as.data.table(ldply(final_modlist, function(mod) {
  regdiagnostic_customtab(mod=mod, kCV = TRUE, k=10, cvreps=50)
}))
final_modsummary[, `:=`(Elem=c(rep(c('Cu', 'Pb', 'Zn'), each=2), c('Cu', 'Pb', 'Zn')),
                       type = c(rep(c('With outliers', 'Without outliers'), each=3),
                                rep('With outliers', 3)),
                       method = c(rep('in situ XRF - ICP', 6), 
                                  rep('lab. XRF - ICP', 3)),
                       N =  unlist(lapply(final_modlist, function(mod) nrow(mod$data)))
                       )]
cat(latex_format(final_modsummary),
    file = file.path(moddir, 'modeltable_XRFvalidation_20201029.tex'))


#--------- Make scatterplots ---------------
plot_modelpreds <- function(mod, pred, elem, dt, stdbreaks, outliers=NULL,
                            point.color = '#238443', bottom.row=TRUE, 
                            ylims) {
  p <- jtools::effect_plot(mod, pred = !!pred,
                           interval = TRUE, plot.points = TRUE, 
                           point.size=3, point.alpha=1/2, point.color=point.color) + 
    scale_x_continuous(breaks=log10(stdbreaks),
                       labels=stdbreaks) + 
    scale_y_log10(name = 'Concentration (ppm) | ICP-OES', limits=ylims) + 
    labs(title=elem) +
    theme_bw() + 
    theme(text=element_text(size=14))
    
  if (!bottom.row) {
    p <- p + 
      theme(axis.title.x = element_blank())
  }
  
  if(!is.null(outliers)) {
    p <- p + geom_point(data=dt[
      (SiteIDPair %in% outliers) &  Elem == elem,],
      aes(x=get(pred), y=ICP), color='black', size=3, alpha=1/2)
  }
  
  return(p)
}

pCu_1 <- plot_modelpreds(
  mod = ICPfieldmodCu[[2]],
  pred = 'logfield',
  dt = ICPfieldsub,
  stdbreaks = c(0.05, 0.1, 0.25, 0.5),
  ylims = c(3, 175),
  outliers = Cuoutliers_fieldICP,
  elem = 'Cu',
  bottom.row = FALSE) +
  theme(axis.title.y = element_blank())

pPb_1 <- plot_modelpreds(
  mod = ICPfieldmodPb[[2]],
  pred = 'logfield',
  dt = ICPfieldsub,
  stdbreaks = c(0.01, 0.05, 0.1, 0.25, 0.5),
  outliers = c(Pboutliers_fieldICP, 
               ICPfieldsub[ICP==0 & Elem=='Pb', SiteIDPair]),
  ylims = c(1, 100),
  elem = 'Pb',
  bottom.row = FALSE)

pZn_1 <- plot_modelpreds(
  mod = ICPfieldmodZn[[2]],
  pred = 'logfield',
  dt = ICPfieldsub,
  stdbreaks = c(0.1, 0.25, 0.5, 1),
  outliers = Znoutliers_fieldICP,
  ylims = c(10, 350),
  elem = 'Zn',
  bottom.row = TRUE)+
  xlab('Normalized count | in-situ XRF') +
  theme(axis.title.y = element_blank(),
        plot.background = element_blank())

pCu_2 <- plot_modelpreds(
  mod = ICPlabmodCu[[2]],
  pred = 'loglab',
  dt = ICPlabsub,
  stdbreaks = c(0.05, 0.1, 0.25, 0.5, 1),
  ylims = c(3, 175),
  elem = 'Cu', 
  point.color = '#225ea8',
  bottom.row = FALSE)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        title = element_blank())

pPb_2 <- plot_modelpreds(
  mod = ICPlabmodPb[[2]],
  pred = 'loglab',
  dt = ICPlabsub,
  stdbreaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 1),
  elem = 'Pb',
  ylims = c(1, 100),
  outliers = ICPlabsub[ICP==0 & Elem=='Pb', SiteIDPair],
  point.color = '#225ea8',
  bottom.row = FALSE)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        title = element_blank())

pZn_2 <- plot_modelpreds(
  mod = ICPlabmodZn[[2]],
  pred = 'loglab',
  dt = ICPlabsub,
  stdbreaks = c(0.025, 0.05, 0.1, 0.25, 0.5, 1),
  ylims = c(10, 350),
  elem = 'Zn',
  point.color = '#225ea8',
  bottom.row = TRUE)+
  xlab('Normalized count | lab. XRF') +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        plot.background = element_blank())

png(file.path(moddir, 
              paste0('XRFvalidation_20201029.png')), 
    width = 6, height=9, units='in', res=600)
(pCu_1 | pCu_2) /
(pPb_1 | pPb_2) /
(pZn_1 | pZn_2)
dev.off()
