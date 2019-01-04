#options(warn=-1)
library(data.table)
library(magrittr)
library(xlsx)
library(plyr)
library(stringr)
library(foreign)
library(rgdal)
library(raster)
library(ggplot2)
library(ggstance)
library(ggpmisc)
library(PeriodicTable)
library(kableExtra)
library(dplyr)
library(hexbin)
library(colorspace)
library(plotly)
library(listviewer)
library(mvoutlier)
library(rrcov)
data(periodicTable)

rootdir <- "C:/Mathis/ICSL/stormwater" #UPDATE
resdir <- file.path(rootdir, "results")
datadir <- file.path(rootdir, "data")
inspectdir <- file.path(resdir, "data_inspection")

if (dir.exists(inspectdir)) {
  warning(paste(inspectdir, 'already exists'))
} else {
  warning(paste('Create new directory:',inspectdir))
  dir.create(inspectdir)
}

#Element colors: http://jmol.sourceforge.net/jscolors/

#Define functions
deconvolution_import <- function(artaxdir, idstart, idend) {
  #Collate all results from ARTAX deconvolution outputs
  tablist <- list.files(artaxdir)
  for (i in 1:length(tablist)) {
    res <- read.csv(file.path(artaxdir,tablist[i]))
    res$XRFID <- str_extract(substr(tablist[i], idstart, idend), "[0-9A-B_]+")
    if (i==1){
      deconvres <- res
    } else{
      deconvres <- rbind(deconvres, res)
    }
  }
  return(setDT(deconvres))
}

vispal <- function(dat) {
  #Create a color palette graph
  barplot(rep(1,length(dat)), col=dat)
}

############################################################################################################################################
# Import data 
############################################################################################################################################
########### ---- Import and format field data ---- ####
fieldata <- read.csv(file.path(datadir,"field_data/field_data_raw_20180808_edit.csv"))
colnames(fieldata)[1] <- 'SiteID'
fieldata$SiteID <- as.character(fieldata$SiteID)
fieldata_sel <- fieldata[!is.na(fieldata$XRFmin) & !is.na(fieldata$SiteID),] #Remove extraneous sites with no XRF data or just for TNC tour
fieldata_format <- data.frame()
for (row in seq(1,nrow(fieldata_sel))) {
  extract <- fieldata_sel[row,]
  for (xrf in seq(fieldata_sel[row,'XRFmin'], fieldata_sel[row,'XRFmax'])){
    #print(xrf)
    extract$XRFID <- xrf
    fieldata_format <- rbind(fieldata_format, extract)
  }
}

########### ---- Import and format field XRF deconvolution results ---- ####
fieldXRF <- deconvolution_import(file.path(datadir, 'XRF20181210/PostBrukerCalibration/deconvolutionresults_XRF12_245_20181215'),
                                 idstart = 41, idend= 43)

########### ---- Import raw XRF lab data for inspection ---- ####
labxrf_list <- data.table(filename = grep('\\.csv$', list.files(file.path(datadir, 'XRF20181210/PelletMeasurements')) ,value=T))
for (i in labxrf_list$filename) {
  print(i)
  xrfrec <- read.csv(file.path(datadir, 'XRF20181210/PelletMeasurements', i))
  labxrf_list[filename == i, `:=`(cps = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Valid Count Last Packet',])),
                                  c_cum = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Valid Accumulated Counts',])),
                                  duration = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Duration Time',]))
  )]
}
setDT(labxrf_list)[, `:=`(c_cum_ps = c_cum/duration,
                          SiteID = regmatches(labxrf_list$filename, regexpr('[0-9A-Z]+[_][0-9]', labxrf_list$filename)))]

#List of initially problematic samples that had to be re-measured
# problem_recs <- setDT(labxrf_list)[SiteID %in% c('3A_1', '3A_2','6B_2','7A_1', '7A_2','19A_1', '19A_2', '20A_1',
#                                                  '20B_1', '20B_2','22A_1', '22A_2', '28B_2', '32A_1', '33B_1', 
#                                                  '33B_2','34A_1','36A_2', '37A_1', '37A_2','40A_1', '42A_1', 
#                                                  '42B_1','42B_2','44A_1', '44A_2', '45A_1', '45A_2','46A_1','46A_2',
#                                                  '48A_1','49A_1', '49A_2','53A_1', '53A_2', '54B_2','55B_1', '55B_2',
#                                                  '57A_1', '57A_2','59A_1', '59A_2', '61A_1', '61A_2'),]
labxrf_list[c_cum_ps > 90000,][!(labxrf_list[c_cum_ps > 90000, SiteID] %in% problem_recs$SiteID)]

########### ---- Import and format lab XRF deconvolution results ---- ####
labXRF <- deconvolution_import(file.path(datadir, 'XRF20181210/PelletMeasurements/deconvolutionresults_labXRF_20181215'),
                               idstart = 29, idend= 33)
########### ---- Import ICP-OES data ---- ####
ICPdat <- read.xlsx(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'),
                    1, stringsAsFactors = F)
ICPthresholds <- read.xlsx(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'),
                           2, stringsAsFactors = F)

########### ---- Import GIS data (including pollution variables) ---- ####
trees <- as.data.table(readOGR(dsn = file.path(resdir, 'Seattle_sampling.gdb'), layer = 'XRFsites_proj')) 
#summary(trees)
heatcols <- colnames(trees)[grep('heat', colnames(trees))]
trees[, (heatcols) := lapply(.SD, function(x){x[is.na(x)] <- 0; x}), .SDcols = heatcols]
############################################################################################################################################
#Format XRF data
############################################################################################################################################
#Cast while summing net photon counts across electron transitions
fieldXRFcast <- dcast(setDT(fieldXRF), XRFID~Element, value.var='Net', fun.aggregate=sum) 
fieldXRFcast[, XRFID := as.numeric(gsub('[_]', '', XRFID))] #Format site number
fieldXRFcast[fieldXRFcast < 0] <- 0 #Floor negative net photon count to 0

labXRFcast <- dcast(setDT(labXRF), XRFID~Element, value.var='Net', fun.aggregate=sum)
labXRFcast[labXRFcast < 0] <- 0

#Average lab XRF results over multiple measurements for a given pellet
labXRFcast <- labXRFcast[, `:=`(SiteID = gsub('[A-B].*', '', XRFID),
                                Pair = gsub('[0-9_]+', '', XRFID),
                                XRFID = NULL)] %>%
  .[, lapply(.SD, mean), by = .(SiteID, Pair)] 

#Normalize data by Rhodium photon count for field and lab results
fieldXRFcastnorm <- fieldXRFcast[, lapply(.SD, function(x) {x/Rh}), by = XRFID]
labXRFcastnorm <- labXRFcast[, lapply(.SD, function(x) {x/Rh}), by = .(SiteID, Pair)]

#Merge datasets: lab XRF + field XRF + field variables
fieldt <- setDT(fieldata_format)[fieldXRFcastnorm, on='XRFID']
labdt <- setDT(fieldata_sel)[labXRFcastnorm, on =  .(SiteID, Pair)]

#Compute average, sd, and range, then melt XRF field data
field_artaxstats <- fieldt[, sapply(.SD, function(x) list(mean=mean(x, na.rm=T),
                                                          sd=sd(x, na.rm=T),
                                                          range=max(x,na.rm=TRUE)-min(x,na.rm=TRUE))), 
                           by=c('SiteID','Pair'), .SDcols = 29:57]
setnames(field_artaxstats, c('SiteID', 'Pair', 
                             paste(rep(colnames(fieldt)[29:57],each=3),
                                   c('mean', 'sd', 'range'), sep='_')))

fieldXRF_format <- melt(field_artaxstats, id.vars=c('SiteID','Pair'), variable.name='Elem_stats') %>%
  .[, `:=`(Elem = sub('(.*)[_](.*)', '\\1', Elem_stats),
           stats = sub('(.*)[_](.*)', '\\2', Elem_stats))] %>% #Replaces the whole match with the first group in the regex pattern
  dcast(SiteID+Pair+Elem~stats, value.var = 'value') %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

#Melt XRF lab data
labXRF_format <- melt(labdt[, colnames(labdt) %in% c('SiteID', 'Pair', periodicTable$symb), with=F], 
                      id.vars=c('SiteID','Pair'), variable.name='Elem') %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

#Write out mean field XRF values (should write this directly to XRF sites proj)
field_artaxmean <- field_artaxstats[, .SD, .SDcols = c(1,2, grep('mean', colnames(field_artaxstats)))] 
colnames(field_artaxmean) <- gsub('_mean', '', colnames(field_artaxmean))
write.dbf(field_artaxmean, 'field_artaxmean_20180827.dbf')

############################################################################################################################################
# Inspect data and remove outliers
############################################################################################################################################
########### ---- Field XRF ---- ###########
str(fieldXRF_format)
# ---- Check distribution by element ----
fXRF_distrib <- ggplot(fieldXRF_format, aes(x=mean)) + 
  geom_histogram() + 
  facet_wrap(~Elem, scales='free') + 
  theme_classic()
ggplotly(fXRF_distrib)

# ---- Assess within-tree variability ----
#Plot coefficient of variation distributions for every element
cvmean_labeltree <- as.data.frame(fieldXRF_format[!(fieldXRF_format$Elem %in% c('Rh','Pd','Ar')),
                                                  paste0('Mean CV: ',format(mean(cv, na.rm=T),digits=2)),
                                                  by=name])
ggplot(fieldXRF_format[!(fieldXRF_format$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density()+
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name) + 
  geom_text(data=cvmean_label, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  theme_classic() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14))

#Plot relationship between elemental concentration and cv
ggplotly(
  ggplot(fieldXRF_format[!(fieldXRF_format$Elem %in% c('Rh','Pd','Ar')),], 
         aes(x=mean, y=cv, color=name, label=paste0(SiteID, Pair))) + 
    geom_point()+
    geom_smooth() +
    labs(x='Mean photon count (normalized)', y='Coefficient of variation') + 
    facet_wrap(~name, scales='free') +
    theme_classic() + 
    theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
          strip.text = element_text(size=14))
)

#
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
            axis.line = element_line(color='black'))
  )
  dev.off()
}

# ---- Assess within-site variability ---- 
#TO DO: CREATE RANDOM PALETTE
field_artaxmeansite <- fieldt[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID'), .SDcols=29:51]
artaxsdsite <- fieldt[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), .SDcols=29:51]
fieldXRF_formatsite <- merge(melt(field_artaxmeansite, id.vars='SiteID', variable.name='Elem', value.name='mean'),
                             melt(artaxsdsite, id.vars='SiteID', variable.name='Elem', value.name='sd'),
                             by=c('SiteID','Elem')) %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')
cvmean_labelsite <- as.data.frame(fieldXRF_formatsite[!(fieldXRF_formatsite$Elem %in% c('Rh','Pd','Ar')),
                                                      paste0('Mean CV: ',format(mean(cv, na.rm=T),digits=2)),
                                                      by=name])

ggplot(fieldXRF_formatsite[!(fieldXRF_formatsite$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density()+
  geom_text(data=cvmean_label, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name) + 
  theme_bw() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14),
        panel.border = element_blank(),
        axis.line = element_line(color='black'))

for (elem in unique(fieldXRF_format[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
  print(elem)
  png(file.path(inspectdir, paste0('fieldXRF_withinsite_',elem,'.png')), width = 20, height=12, units='in', res=300)
  print(
    ggplot(fieldXRF_format[Elem == elem,], 
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

# ---- Univariate flag based on within-tree CV for different metals ----
fieldXRF_format[, `:=`(meanCV=mean(cv, na.rm=T),
                       sdCV = sd(cv, na.rm=T)), by=Elem]
CV_flag <- dcast(fieldXRF_format[cv>(meanCV+2*sdCV),], SiteID+Pair~Elem, value.var='cv') 
NAcount <- function(x) length(which(!is.na(unlist(x))))
CV_flag[, flag_count := NAcount(.SD)-2, by = seq_len(nrow(CV_flag))]
ggplot(fieldXRF_format, aes(x=mean, y=cv)) + geom_point() +scale_x_sqrt()

#---- Multivariate outlier detection ---- 
#(see https://stats.stackexchange.com/questions/213/what-is-the-best-way-to-identify-outliers-in-multivariate-data for reference discussion)
#https://rpubs.com/Treegonaut/301942
#See Research_resources/statistics/outliers
excol <- c(fieldXRF_format[meanCV>0.5,unique(Elem)], 'Rh', 'Pd', 'Rb') #columns to exclude

#Compare classic vs robust Mahalanobis distance
outlierdat <- field_artaxmean[,-c('SiteID','Pair', excol), with=F]
par(mfrow=c(1,1))
distances <- dd.plot(outlierdat, quan=1/2, alpha=0.025)
outlierdat_dist <- cbind(field_artaxmean[,-excol, with=F], as.data.frame(distances))

ggplotly(
  ggplot(outlierdat_dist, aes(x=md.cla, y=md.rob, color=outliers, label=paste0(SiteID,Pair))) +
    geom_point(size=3)
)

#Filzmoser et al. multivariate outlier detection
outliers <- aq.plot(outlierdat, delta=qchisq(0.975, df = ncol(outlierdat)), quan = 1/2, alpha = 0.05)
par(mfrow=c(1,1))
outlierdat_dist <- cbind(outlierdat_dist, aqplot_outliers = outliers$outliers)

outlierdat_distmelt <- melt(outlierdat_dist[, colnames(outlierdat_dist) %in% c('SiteID', 'Pair', 'aqplot_outliers', 'md.rob',periodicTable$symb), with=F], 
                            id.vars=c('SiteID','Pair','aqplot_outliers','md.rob'), variable.name='Elem')

#Check out univariate distribution of outliers detected by mvoutlier
ggplot(outlierdat_distmelt, aes(x='O', y=value, color=aqplot_outliers)) + 
  geom_jitter() +
  #scale_color_distiller(palette= 'Spectral') +
  facet_wrap(~Elem, scales='free') +
  theme_classic()

#Robust PCA with rrcov (see 4.2 p24 of Todorov and Filzmoser 2009)
cols <- colnames(outlierdat_dist)[colnames(outlierdat_dist) %in% periodicTable$symb]
outlierdat_dist[, (cols) := lapply(.SD, scale), .SDcols=cols]

pca <- PcaClassic(~., data=outlierdat_dist[,cols,with=F])
summary(pca)
biplot(pca, main="Classic biplot", col=c("gray55", "red"))
rpca <- PcaGrid(~., data=outlierdat_dist[,cols,with=F])
summary(rpca)
screeplot(rpca, type="lines", main="Screeplot: robust PCA")
biplot(rpca, main="Robust biplot", col=c("gray55", "red"))
plot(rpca)
plot(PcaHubert(outlierdat_dist[,cols,with=F], k=3))

#FastPCS (rrcov package)? 

#Relate to lab XRF and ICP-OES, then look at residuals: 
#For each element: Cook's distance, Lund's test, mvoutlier corr.plot 
#Outlier in terms of multivariate relationship between field XRF and lab XRF or field XRF and ICP-OES
#MVN: An R Package for Assessing Multivariate Normality
#Look at outliers spatially

##### ---- Lab XRF ---- ####

##### ---- ICP-OES ---- ####

#Inspect outliers: where were they collected, is one measurement causing the main issue, how does the deconvolution fit?
#Test influence of method
#Consider excluding all variables with within-tree CV > 0.5

############################################################################################################################################
# Check whether tree species mattered at all
############################################################################################################################################
############################################################################################################################################
#Relate XRF data to ICP results
############################################################################################################################################
#Format data 
ICPdat[ICPdat == 'TR'] <- '0'
ICPdat[ICPdat == 'ND'] <- '0'
ICPdat <- setDT(ICPdat)[!(SAMPLE.SET %in% c('BLK', 'QC-1', NA)),]
ICPdat <- ICPdat[, lapply(.SD, as.character), by=SAMPLE.SET]
ICPdat <- ICPdat[, lapply(.SD, as.numeric), by=SAMPLE.SET]

#Set aside and compute mean over duplicate measurements
ICPdup <- ICPdat[substr(ICPdat$SAMPLE.SET, 1,3) %in% substr(grep('DUP', ICPdat$SAMPLE.SET, value=T), 1,3),]
ICPdat[, SAMPLE.SET := substr(SAMPLE.SET, 1,3)] 
ICPmean <- ICPdat[, lapply(.SD, mean), by= SAMPLE.SET]

#Melt dataset to merge with XRF data
ICPmelt <- melt(setDT(ICPdat), id.vars = 'SAMPLE.SET', variable.name = 'Elem', value.name = 'ICP')
ICPmelt[, `:=`(SiteID = gsub('[A-Z]', '', SAMPLE.SET),
               Pair = gsub('[0-9]', '', SAMPLE.SET))]

#----
#Relate XRF field data to ICP results
#----
#Merge with XRF field data
ICPfieldmerge <- ICPmelt[fieldXRF_format, on = .(SiteID, Pair, Elem)]
unique(ICPmelt$Elem)[!(unique(ICPmelt$Elem) %in% unique(fieldXRF_format$Elem))] #Check what elems are in vs out
unique(fieldXRF_format$Elem)[!(unique(fieldXRF_format$Elem) %in% unique(ICPmelt$Elem))] #Check what elems are in vs out

#Plot comparison 
outlier_sites <- c(20, 23, 43, 49, 52, 53, 54, 62) 
ICPfield_plot <- ggplot(ICPfieldmerge[!(is.na(ICP) | 
                                          ICPfieldmerge$SiteID %in% outlier_sites | 
                                          Elem %in% c('Na', 'Si', 'P')
),], 
aes(x=mean, y=ICP)) + 
  geom_point(aes(label = SiteID, color = cv)) +
  scale_color_distiller(palette= 'Spectral') +
  geom_smooth() +
  geom_linerangeh(aes(xmin=mean-sd, xmax=mean+sd)) + 
  scale_y_continuous(expand=c(0,0), limits=c(0, NA)) +
  scale_x_continuous(expand=c(0,0), limits=c(0, NA)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
#plotly_json(ICPfield_plot)
ggplotly(ICPfield_plot)%>%
  style(hoverinfo = "none", traces = 20:77)
ICPfieldmerge[!is.na(ICP)& !(SiteID %in% outlier_sites),
              round(summary(lm(ICP~mean))$adj.r.squared, 2),
              by=Elem]
#Think about removing sites with very high sd

#Plot comparison averaging within sites
ICPfieldmerge_sitemean <- ICPfieldmerge[!is.na(ICP),
                                        list(ICPmean = mean(ICP, na.rm=T), XRFmean = mean(mean, na.rm=T)),
                                        .(SiteID, Elem)]
ICPfield_plot <- ggplot(ICPfieldmerge_sitemean[!(is.na(ICPmean) | Elem %in% c('Na', 'Si', 'P') |
                                                   SiteID %in% outlier_sites),], 
                        aes(x=XRFmean, y=ICPmean)) + 
  geom_point(aes(label = SiteID)) +
  geom_smooth(method='lm') +
  facet_wrap(~Elem, scales='free') +
  scale_y_continuous(expand=c(0,0), limits=c(0, NA)) +
  scale_x_continuous(expand=c(0,0), limits=c(0, NA)) +
  theme_classic()
ggplotly(ICPfield_plot)

ICPfieldmerge_sitemean[!is.na(ICPmean)& !(SiteID %in% outlier_sites),
                       round(summary(lm(ICPmean~XRFmean))$adj.r.squared, 2),
                       by=Elem]

#----
#Relate XRF lab data to ICP results
#----
#Merge with XRF lab data
ICPlabmerge <- ICPmelt[labXRF_format, on = .(SiteID, Pair, Elem)]
unique(ICPmelt$Elem)[!(unique(ICPmelt$Elem) %in% unique(labXRF_format$Elem))] #Check what elems are in vs out
unique(labXRF_format$Elem)[!(unique(labXRF_format$Elem) %in% unique(ICPmelt$Elem))] #Check what elems are in vs out

#Plot comparison
ICPlab_plot <- ggplot(ICPlabmerge[!(is.na(ICP) | SiteID %in% c('49', '62')),], aes(x=value, y=ICP, label=SiteID)) + 
  geom_point() + 
  geom_smooth(span=0.9) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()

#plotly_json(ICPlab_plot)
ggplotly(ICPlab_plot)%>%
  style(hoverinfo = "none", traces = 24:62)

ICPlabmerge[!(is.na(ICP) | SiteID %in% c('49', '62')),
            round(summary(lm(ICP~value))$adj.r.squared, 2), by=Elem]

#--------------------------------------------------------------------------------------
#Relate XRF field data to XRF lab data
fieldlabmerge <- fieldXRF_format[labXRF_format, on = .(SiteID, Pair, Elem)]
#Plot comparison
fieldlab_plot <- ggplot(fieldlabmerge, aes(x=mean, y=value, label=SiteID)) + 
  geom_point() + 
  geom_smooth(span=0.9) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
ggplotly(fieldlab_plot)

############################################################################################################################################
#Relate XRF data to pollution predictors 
############################################################################################################################################

#Import shapefile of sites
colnames(trees)[1] <- 'SiteID'

#Join field data to site shapefile
treesxrf <- merge(trees, field_artaxmean, by=c('SiteID','Pair'))
treesxrf <- treesxrf[!is.na(treesxrf$Fe) | treesxrf$SiteID==1,]

#Format data
metal_select <- c('Br','Co','Cr','Cu','Fe','Mn','Ni','Pb','Sr','Ti','V','Zn')
#metal_select <- c('Fe','Zn','Co','Cr') #for Phil 2018/10/02 'goodresults'
#metal_select <- c('Ca','K','Sr','P') #for Phil 2018/10/02 'randomresults'

treesxrf_melt <- melt(treesxrf, 
                      id.vars=colnames(treesxrf)[!(colnames(treesxrf) %in% metal_select)], 
                      variable.name = 'Elem')
treesxrf_melt <- merge(treesxrf_melt, periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')


#Plot Bing Congestion index vs. XRF data
bing_xrf <- ggplot(treesxrf_melt, aes(x=heat_binglog300_1, y=value, label = SiteID)) + 
  # geom_point(aes(color=colsmix), size=3) +
  # scale_color_identity() +
  geom_point(aes(color=heatAADTlog100), size=3) +
  scale_color_distiller(palette='Spectral') +
  labs(x='Traffic congestion (Bing index)', y='Concentration in moss (normalized photon count)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', se = T, color='black') + 
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 5)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        text = element_text(size=16),
        axis.title = element_blank())
ggplotly(bing_xrf)


png(file.path(resdir, 'prelim_badXRFresults.png'), width=8, height=8, units = 'in', res=400)
#pdf(file.path(resdir, 'prelim_badXRFresults.pdf'), width=8, height=8)
bing_xrf
dev.off()

#Plot traffic volume index vs. XRF data
ggplot(treesxrf_melt, aes(x=heatAADTlog100, y=value)) + 
  geom_point() + 
  labs(x='Traffic volume index', y='Photon count (normalized)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic()

#Plot speed limit index vs. XRF data
ggplot(treesxrf_melt, aes(x=spdlm, y=value)) + 
  geom_point() + 
  labs(x='Speed limit index', y='Photon count (normalized)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic()


#Multiple linear regressions
Znlm <- lm(Zn~heat_binglog300_1 + sqrt(heatOSMAADTlog300) + sqrt(heatbustransitlog200), 
           data=treesxrf) #[!(treesxrf$SiteID %in% c('19','23')),]
summary(Znlm)
ggplotly(qplot(x=predict(Znlm), y=treesxrf$Zn, color=treesxrf$heatbustransitlog200, label=as.character(treesxrf$SiteID)) + geom_smooth(method='lm') )
#+ geom_text(aes(label=treesxrf$SiteID))

Felm <- lm(Fe~heat_binglog300_1 + sqrt(heatOSMAADTlog300) + sqrt(heatbustransitlog200), data=treesxrf)
summary(Felm)
qplot(x=predict(Felm), y=treesxrf$Fe) + geom_smooth(method='lm')

Culm <- lm(Cu~sqrt(AADT), data=treesxrf)
summary(Culm)
qplot(x=predict(Culm), y=treesxrf$Cu) + geom_smooth(method='lm') 

Pblm <- lm(Pb~bing+sqrt(AADT)+spdlm, data=treesxrf)
summary(Pblm)
ggplotly(qplot(x=predict(Pblm), y=treesxrf$Pb, label=treesxrf$SiteID) + geom_smooth(method='lm'))

############################################################################################################################################
# Relate lab XRF data to pollution predictors 
############################################################################################################################################
colstodelete <- colnames(treesxrf_melt) %in% c(periodicTable$symb, 'value', 'valuestand', 'valuefact')
treeslab <- merge(setDT(treesxrf_melt)[, .SD, .SDcols = !colstodelete], labXRF_format, by = c('SiteID', 'Pair', 'Elem'))

binglab<- ggplot(treeslab[SiteID != 49,], aes(x=bing, y=value, label= SiteID)) + 
  geom_point(aes(color = AADT)) +
  scale_color_distiller(palette='Spectral') +
  labs(x='Traffic congestion (bing index)', y='Concentration (ug/g)') + 
  facet_wrap(~name.x, nrow=3, scales='free_y') +
  geom_smooth(method='lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic()
ggplotly(binglab)

Znlm <- lm(value~heat_binglog300_1 + sqrt(heatOSMAADTlog300) + sqrt(heatbustransitlog200), 
           data=treeslab[Elem=='Zn' & !(treeslab$SiteID %in% c('44','49','51','52','53','54')),]) #[treesxrf$SiteID != '51',]
summary(Znlm)
ggplotly(qplot(x=predict(Znlm),
               y=as.list(treeslab[Elem=='Zn' &  !(treeslab$SiteID %in% c('49','51','52','53','54','44')),])$value, 
               label=as.list(treeslab[Elem=='Zn' &  !(treeslab$SiteID %in% c('49','51','52','53','54','44')),])$SiteID) + geom_smooth(method='lm') )


#--------------------------------------------------------------------------------------
# Relate ICP data to pollution predictors 
treesICP <- merge(setDT(treesxrf_melt)[, .SD, .SDcols = !colstodelete], ICPmelt, by = c('SiteID', 'Pair', 'Elem'))

#
bingICP<- ggplot(treesICP, aes(x=bing, y=ICP, label = paste(SiteID, Pair))) + 
  geom_point(aes(color=AADT)) +
  scale_color_distiller(palette='Spectral') +
  labs(x='Traffic congestion (bing index)', y='Concentration (ug/g)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic()
ggplotly(bingICP)

#Speed limit
ggplot(treesICP, aes(x=spdlm, y=ICP)) + 
  geom_point() + 
  labs(x='Speed limit index', y='Concentration (ug/g)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic()

#Traffic volume
ggplot(treesICP, aes(x=AADT, y=ICP)) + 
  geom_point() + 
  labs(x='Traffic volume index', y='Photon count (normalized)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic()

Znlm <- lm(ICP~bing+sqrt(AADT), data=treesICP[Elem=='Zn' & SiteID != 49,]) #[treesxrf$SiteID != '51',]
summary(Znlm)
ggplotly(qplot(x=predict(Znlm),
               y=as.list(treesICP[Elem=='Zn' & SiteID != 49,])$ICP, 
               label=as.list(treesICP[Elem=='Zn' & SiteID != 49,])$SiteID) + geom_smooth(method='lm') )


#--------------------------------------------------------------------------------------
# Output overall variability table 
elemvar <- data.frame(Elem=colnames(fieldXRFcastnorm[,-1]),
                      mean=colMeans(fieldXRFcastnorm[,-1]),
                      min =t(setDT(fieldXRFcastnorm[,-1])[ , lapply(.SD, min)]),
                      max =t(setDT(fieldXRFcastnorm[,-1])[ , lapply(.SD, max)]),
                      sd= t(setDT(fieldXRFcastnorm[,-1])[ , lapply(.SD, sd)]))

#Generate table of statistics
dt$SiteIDPair <- with(dt, paste0(SiteID,Pair))
allelem <- colnames(dt) %in% periodicTable$symb #Get all column which correspond to an element in the periodic table
#Apply to each element:
elemwise_comps <- sapply(colnames(dt)[allelem], function(x) {
  #Analyze error within tree samples
  withintree <- summary(aov(as.formula(paste(x, "~SiteIDPair")), data=dt))  
  #Analyze error among tree samples within sites
  withinsite <- summary(aov(as.formula(paste(x, "~SiteID")), data=dt)) 
  #Regress pollutant predictors against XRF data for each element
  bing_lm <- summary(lm(as.formula(paste(x, "~bing")), data=treesxrf))
  aadt_lm <- summary(lm(as.formula(paste(x, "~sqrt(AADT)")), data=treesxrf))
  spdlm_lm <- summary(lm(as.formula(paste(x, "~sqrt(spdlm)")), data=treesxrf))
  
  return(c(tree_error=withintree[[1]]$`Sum Sq`[2]/sum(withintree[[1]]$`Sum Sq`),
           site_error=(withinsite[[1]]$`Sum Sq`[2]-withintree[[1]]$`Sum Sq`[2])/sum(withinsite[[1]]$`Sum Sq`),
           total_error = withinsite[[1]]$`Sum Sq`[2]/sum(withinsite[[1]]$`Sum Sq`),
           bing_r2 = bing_lm$adj.r.squared, 
           bing_sig = bing_lm$coefficients[8],
           aadt_r2 = aadt_lm$adj.r.squared,
           aadt_sig = aadt_lm$coefficients[8],
           spdlm_r2 = spdlm_lm$adj.r.squared,
           spdlm_sig = spdlm_lm$coefficients[8]
  ))
})
elemvar <- cbind(elemvar, t(elemwise_comps))
#Get rid of useless elements (signals from XRF Tracer detector and collimator, not in sample)
elemvar <- elemvar[!(elemvar$Elem %in% c('Ar','Rh','Pd')),]
#Get full element names
elemvar <- merge(elemvar, periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

#Format and output table
options(knitr.table.format = "html") 
elemvar[order(-elemvar$bing_r2),] %>%
  mutate(
    congestion_r2 = cell_spec(format(bing_r2, digits=3), bold = ifelse(bing_sig < 0.05, TRUE, FALSE)),
    volume_r2 = cell_spec(format(aadt_r2, digits=3), bold = ifelse(aadt_sig < 0.05, TRUE, FALSE)),
    speedlimit_r2 = cell_spec(format(spdlm_r2, digits=3), bold = ifelse(spdlm_sig < 0.05, TRUE, FALSE))
  ) %>%
  select(-c(9:14)) %>%
  kable(format = "html", escape = F) %>%
  kable_styling("striped", full_width = F) 

#%>%
#  save_kable(file.path(resdir, 'XRFrestab_20180828.doc'), self_contained=T)


