#Notes:
# think about using GLMs rather than Multiple Linear Regressions to avoid transforming the response variable?
# check residuals for homogeneity of variance
# think how to deal with spatial autocollinearity? (2_d smoother of coordinates?)
# Univariate: Cook's distance, Lund's test, mvoutlier corr.plot 
# Look at potential outliers spatially
# Check collinearity among pollution variables and among elements
# Could also use OutliersO3 package (https://blog.revolutionanalytics.com/2018/03/outliers.html)

#---- Import libraries ----
#options(warn=-1)
library(rprojroot)
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
library(gginnards)
library(grid)
library(gridExtra)
library(PeriodicTable)
library(kableExtra)
library(dplyr)
library(hexbin)
library(colorspace)
library(plotly)
library(listviewer)
library(mvoutlier)
library(rrcov)
library(vegan)
library(olsrr) #https://cran.r-project.org/web/packages/olsrr/vignettes/intro.html
data(periodicTable)

#---- Define directory structure ----
rootdir <- find_root(has_dir("src"))
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

#---- Define functions----
#Source some internal functions
source(file.path(rootdir, 'src/stormwater_samplingR/internal_functions.R'))

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

shapiro_test_df <- function(df, bonf= F, alpha= 0.05) {
  #Perform a Shapiro-Wilks test on every variable in df from:
  #https://stackoverflow.com/questions/33489330/how-to-test-the-normality-of-many-variables-in-r-at-the-same-time
  l <- lapply(df, shapiro.test)
  s <- do.call("c", lapply(l, "[[", 1))
  p <- do.call("c", lapply(l, "[[", 2))
  if (bonf == TRUE) {
    sig <- ifelse(p > alpha / length(l), "H0", "Ha")
  } else {
    sig <- ifelse(p > alpha, "H0", "Ha")
  }
  return(list(statistic= s,
              p.value= p,
              significance= sig,
              method= ifelse(bonf == TRUE, "Shapiro-Wilks test with Bonferroni Correction",
                             "Shapiro-Wilks test without Bonferroni Correction")))
}

ols_plot_dfbetas_custom <- function(model) {
  #Custom adaptation of ols_plot_dfbetas: just returns ggplots instead of grobs
  obs <- NULL
  txt <- NULL
  dfb <- dfbetas(model)
  n <- nrow(dfb)
  np <- ncol(dfb)
  threshold <- 2/sqrt(n)
  myplots <- list()
  outliers <- list()
  for (i in seq_len(np)) {
    dbetas <- dfb[, i]
    df_data <- tibble(obs = seq_len(n), dbetas = dbetas)
    d <- ols_prep_dfbeta_data(df_data, threshold)
    f <- ols_prep_dfbeta_outliers(d)
    p <- eval(substitute(ggplot(d, aes(x = obs, y = dbetas, 
                                       label = txt, ymin = 0, ymax = dbetas)) + geom_linerange(colour = "blue") + 
                           geom_hline(yintercept = c(0, threshold, -threshold), 
                                      colour = "red") + geom_point(colour = "blue", 
                                                                   shape = 1) + xlab("Observation") + ylab("DFBETAS") + 
                           ggtitle(paste("Influence Diagnostics for", colnames(dfb)[i])) + 
                           geom_text(hjust = -0.2, nudge_x = 0.15, size = 2, 
                                     family = "serif", fontface = "italic", colour = "darkred", 
                                     na.rm = TRUE) + annotate("text", x = Inf, y = Inf, 
                                                              hjust = 1.5, vjust = 2, family = "serif", fontface = "italic", 
                                                              colour = "darkred", label = paste("Threshold:", round(threshold, 
                                                                                                                    2))), list(i = i)))
    myplots[[i]] <- p
    outliers[[i]] <- f
  }
  return(list(plots=myplots, outliers=outliers))
}

regdiagnostic_custom <- function(model, df, chem, flagcol) {
  "Set of diagnostic tests and plots to detect influential and outlying points.
  Note: modifies df in place by adding flags based on test results"
  df <- ICPfieldmerge
  k <- ols_prep_cdplot_data(model)
  d <- ols_prep_outlier_obs(k)
  df[which(df$Elem == chem)[as.data.table(d)[color=='outlier', obs]],
     (flagcol) := get(flagcol) + 1]
  cooksdchart <- ols_plot_cooksd_chart(model)$plot %>%
    delete_layers("GeomText") +
    geom_text(aes(label=df[Elem == chem, paste0(SiteID, Pair)]), size=4) + 
    theme_classic()
  
  #DFBETA (from ols_plot_dfbetas)
  dfbeta_mean <- ols_plot_dfbetas_custom(model)$plots[[1]] %>%
    delete_layers("GeomText") +
    geom_text(aes(label=df[Elem == chem, paste0(SiteID, Pair)]), size=4) + 
    theme_classic()
  dfbeta_int <- ols_plot_dfbetas_custom(model)$plots[[2]] %>%
    delete_layers("GeomText") +
    geom_text(aes(label=df[Elem == chem, paste0(SiteID, Pair)]), size=4) + 
    theme_classic()
  
  outliernums <- unique(bind_rows(ols_plot_dfbetas_custom(model)$outliers)$obs)
  df[which(df$Elem == chem)[outliernums], (flagcol) := get(flagcol) +1]
  
  #Externaly deleted studentized residuals vs. leverage
  residlevplot <- ols_plot_resid_lev(model)$plot %>%
    delete_layers("GeomText") +
    geom_text(aes(label=df[Elem == chem, paste0(SiteID, Pair)]), size=4) + 
    theme_classic()
  outliernums <- setDT(ols_plot_resid_lev(model)$plot$data)[color %in% c('outlier', 'outlier & leverage'), obs]
  df[which(df$Elem == chem)[outliernums], (flagcol) := get(flagcol) +1]
  
  #Straight resid
  residplot <- ols_plot_resid_fit(model) %>%
    delete_layers("GeomText") +
    geom_text(aes(label=df[Elem == chem, paste0(SiteID, Pair)]), size=4) + 
    theme_classic()
  
  #Regression with and without outliers plot
  regoutlier <- ggplot(df[Elem == chem,],
                       aes(x=mean, y=ICP, label=paste0(SiteID, Pair))) + 
    geom_text(aes(color=factor(get(flagcol))), size=5) +  
    geom_smooth(method='lm') + 
    labs(x = paste0('field XRF', chem), y= paste0('ICP', chem)) + 
    theme_classic()
  return(arrangeGrob(regoutlier, residplot, residlevplot, cooksdchart, dfbeta_mean, dfbeta_int))
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
labxrf_list[c_cum_ps > 90000,] #[!(labxrf_list[c_cum_ps > 90000, SiteID] %in% problem_recs$SiteID)]

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
colnames(trees)[1] <- 'SiteID'
trees <- trees[!(SiteID %in% c('NA', NA)),]

########### ---- Define elements associated with car traffic ---- ####
#Elements present in brake linings and emitted brake dust (from Thorpe et al. 2008)
brakelining_elem <- c('Al', 'As', 'Ba', 'Ca', 'Cd', 'Co' ,'Cr', 'Cu', 'Fe', 'K', 'Li', 'Mg', 
                      'Mn', 'Mo', 'Na', 'Ni', 'Pb', 'Sb', 'Se', 'Sr', 'Zn') 
#Elements present in passenger car tyre tread (from Thorpe et al. 2008)
tire_elem <- c('Al', 'Ba', 'Ca', 'Cd', 'Co', 'Cr', 'Cu', 'Fe', 'K', 'Mg', 'Mn', 'Na', 'Ni',
               'Pb', 'Sb', 'Sr', 'Ti', 'Zn')

############################################################################################################################################
#Format XRF and ICP data
############################################################################################################################################
########### ---- Format XRF data ---- ####
# ---- Cast while summing net photon counts across electron transitions ----
fieldXRFcast <- dcast(setDT(fieldXRF), XRFID~Element, value.var='Net', fun.aggregate=sum) 
fieldXRFcast[, XRFID := as.numeric(gsub('[_]', '', XRFID))] #Format site number
fieldXRFcast[fieldXRFcast < 0] <- 0 #Floor negative net photon count to 0

labXRFcast <- dcast(setDT(labXRF), XRFID~Element, value.var='Net', fun.aggregate=sum)
labXRFcast[labXRFcast < 0] <- 0

# ---- Normalize data by Rhodium photon count for field and lab results ----
fieldXRFcastnorm <- fieldXRFcast[, lapply(.SD, function(x) {x/Rh}), by = XRFID]
labXRFcastnorm <- labXRFcast[, lapply(.SD, function(x) {x/Rh}), by = XRFID]

# ---- Merge datasets: lab XRF + field XRF + field variables ----
fieldt <- setDT(fieldata_format)[fieldXRFcastnorm, on='XRFID']

labXRFcastnorm[, `:=`(SiteID = gsub('[A-B].*', '', XRFID),
                  Pair = gsub('[0-9_]+', '', XRFID),
                  XRFID = NULL)]
labdt <- setDT(fieldata_sel)[labXRFcastnorm, on =  .(SiteID, Pair)]

# ---- Compute average, sd, and range for lab XRF results over multiple measurements for a given pellet ----
elemcols <- which(colnames(labdt) %in% periodicTable$symb)
lab_artaxstats <- labdt[, sapply(.SD, function(x) list(mean=mean(x, na.rm=T),
                                                       sd=sd(x, na.rm=T),
                                                       range=max(x,na.rm=TRUE)-min(x,na.rm=TRUE))), 
                        by=c('SiteID','Pair'), .SDcols = elemcols]
setnames(lab_artaxstats, c('SiteID', 'Pair', 
                             paste(rep(colnames(labdt)[elemcols],each=3),
                                   c('mean', 'sd', 'range'), sep='_')))

labXRF_format <- melt(lab_artaxstats, id.vars=c('SiteID','Pair'), variable.name='Elem_stats') %>%
  .[, `:=`(Elem = sub('(.*)[_](.*)', '\\1', Elem_stats),
           stats = sub('(.*)[_](.*)', '\\2', Elem_stats))] %>% #Replaces the whole match with the first group in the regex pattern
  dcast(SiteID+Pair+Elem~stats, value.var = 'value') %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

# ---- Compute average, sd, and range, then melt XRF field data ----
elemcols <- which(colnames(fieldt) %in% periodicTable$symb)
field_artaxstats <- fieldt[, sapply(.SD, function(x) list(mean=mean(x, na.rm=T),
                                                          sd=sd(x, na.rm=T),
                                                          range=max(x,na.rm=TRUE)-min(x,na.rm=TRUE))), 
                           by=c('SiteID','Pair'), .SDcols = elemcols]
setnames(field_artaxstats, c('SiteID', 'Pair', 
                             paste(rep(colnames(fieldt)[elemcols],each=3),
                                   c('mean', 'sd', 'range'), sep='_')))

fieldXRF_format <- melt(field_artaxstats, id.vars=c('SiteID','Pair'), variable.name='Elem_stats') %>%
  .[, `:=`(Elem = sub('(.*)[_](.*)', '\\1', Elem_stats),
           stats = sub('(.*)[_](.*)', '\\2', Elem_stats))] %>% #Replaces the whole match with the first group in the regex pattern
  dcast(SiteID+Pair+Elem~stats, value.var = 'value') %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

# ---- Write out mean XRF values (should write this directly to XRF sites proj) ----
field_artaxmean <- field_artaxstats[, .SD, .SDcols = c(1,2, grep('mean', colnames(field_artaxstats)))] 
colnames(field_artaxmean) <- gsub('_mean', '', colnames(field_artaxmean))
write.dbf(field_artaxmean, 'field_artaxmean_20180827.dbf')

lab_artaxmean <- lab_artaxstats[, .SD, .SDcols = c(1,2, grep('mean', colnames(lab_artaxstats)))] 
colnames(lab_artaxmean) <- gsub('_mean', '', colnames(lab_artaxmean))

# ---- Check field XRF data distribution by element then transform and standardize ----
#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
fieldXRF_format[Elem %in% fieldXRF_format[, length(unique(mean))>1, by=Elem][V1==T, Elem],
                transmean := transformTukey_lambda(mean, start = -2, end = 2,int = 0.025, 
                                                   rastertab=F,rep=100, verbose = FALSE, 
                                                   statistic = 1)[[1]], by=Elem]
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

# ---- Check lab XRF data distribution by element then transform and standardize ----
#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
labXRF_format[Elem %in% labXRF_format[, length(unique(mean))>1, by=Elem][V1==T, Elem],
                transmean := transformTukey_lambda(mean, start = -2, end = 2,int = 0.025, 
                                                   rastertab=F,rep=100, verbose = FALSE, 
                                                   statistic = 1)[[1]], by=Elem]

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

# ---- Re-cast data for multivariate analysis ----
field_transmean <- dcast(fieldXRF_format, SiteID+Pair~Elem, value.var = 'transmean') 
lab_transmean <- dcast(labXRF_format, SiteID+Pair~Elem, value.var = 'transmean') 

########### ---- Format ICP-OES data ---- ####
# ---- Format data and check site overlap between field XRF and ICP----
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

# ---- Melt dataset to merge with XRF data ----
ICPmelt <- melt(setDT(ICPmean), id.vars = 'SAMPLE.SET', variable.name = 'Elem', value.name = 'ICP')
ICPmelt[, `:=`(SiteID = gsub('[A-Z]', '', SAMPLE.SET),
               Pair = gsub('[0-9]', '', SAMPLE.SET))]

# ---- Check ICP data distribution by element then transform and standardize ----
#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
ICPmelt[Elem %in% ICPmelt[, length(unique(ICP))>1, by=Elem][V1==T, Elem],
                transICP := transformTukey_lambda(ICP, start = -2, end = 2,int = 0.025, 
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

# ---- Re-cast data for multivariate analysis ----
ICP_trans <- dcast(ICPmelt, SiteID+Pair~Elem, value.var='transICP')

########### ---- Format pollutant data ---- ####
# ---- Check ICP data distribution by element then transform and standardize ----
#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
treesmelt <- melt(trees[!is.na(SiteID),], id.vars=grep('heat', colnames(trees), value=T, invert=T), 
                   variable.name = 'heatvar' , value.name = 'heatvalue')
treesmelt[, transheat := transformTukey_lambda(heatvalue, start = -2, end = 2,int = 0.025, 
                                          rastertab=F,rep=100, verbose = FALSE, 
                                          statistic = 1)[[1]], by=heatvar]
#Run a shapiro test on each element
trees_trans <- dcast(treesmelt, SiteID+Pair~heatvar, value.var = 'transheat') 

normtest <- shapiro_test_df(trees_trans[,sapply(trees_trans, class)=='numeric' &
                                        sapply(trees_trans, function(x) length(unique(x)))>1, with=F])
normsig <- data.frame(sig = normtest$significance)
normsig$heatvar <- row.names(normsig)
#Plot histogram of transformed data color-coded by whether normally distributed or not
treesnorm_join <- treesmelt[normsig, on='heatvar']
png(file.path(inspectdir, paste0('trees_TukeyTransform.png')), width = 20, height=12, units='in', res=300)
ggplot(treesnorm_join, aes(x=transheat, fill=sig)) + 
  geom_histogram() + 
  facet_wrap(~heatvar, scales='free') + 
  theme_classic()
dev.off()     

# z-standardize by element
treesmelt[, transheat:=scale(transheat), by= heatvar]

# ---- Re-cast data for multivariate analysis ----
trees_trans <- dcast(treesmelt, formula = SiteID+Pair~heatvar,  value.var='transheat')

########### ---- Merge datasets ---- ####
#Merge XRF field data with ICP data
ICPfieldmerge <- merge(ICPmelt, fieldXRF_format, by = c('SiteID', 'Pair', 'Elem'), all.x=F, all.y=F)
unique(ICPmelt$Elem)[!(unique(ICPmelt$Elem) %in% unique(fieldXRF_format$Elem))] #Check what elems are in vs out
unique(fieldXRF_format$Elem)[!(unique(fieldXRF_format$Elem) %in% unique(ICPmelt$Elem))] #Check what elems are in vs out

#Merge XRF field data with XRF lab data
joincols <- c('SiteID', 'Pair', 'Elem')
labfieldmerge <- labXRF_format[fieldXRF_format, on = joincols] %>%
  .[SiteID != 1,]
setnames(labfieldmerge, 
         colnames(labfieldmerge[,-joincols, with=F]),
         gsub('^(?!i)', 'lab_', colnames(labfieldmerge[,-joincols, with=F]), perl=T))
setnames(labfieldmerge, 
         colnames(labfieldmerge[,-joincols, with=F]),
         gsub('^i[.]', 'field_', colnames(labfieldmerge[,-joincols, with=F]), perl=T))

#Join pollutant data to sites
trees_fieldxrf <- trees_trans[field_transmean, on=c('SiteID','Pair')][
  !(is.na(Fe) | SiteID==1),]
trees_labxrf <- trees_trans[lab_transmean, on=c('SiteID','Pair')][
  !(is.na(Fe) | SiteID==1),]
trees_ICP <- trees_trans[ICP_trans, on=c('SiteID', 'Pair')][
  !(is.na(Fe) | SiteID==1),]

############################################################################################################################################
# Inspect data and remove outliers
############################################################################################################################################
######################  ---- Field XRF ---- ###########
str(fieldXRF_format)
# ---- Assess within-tree variability ----
#Plot coefficient of variation distributions for every element
cvmean_labeltree <- as.data.frame(fieldXRF_format[!(fieldXRF_format$Elem %in% c('Rh','Pd','Ar')),
                                                  paste0('Mean CV: ',format(mean(cv, na.rm=T),digits=2)),
                                                  by=name])
ggplot(fieldXRF_format[!(fieldXRF_format$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density()+
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name) + 
  geom_text(data=cvmean_labeltree, 
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
  geom_text(data=cvmean_labelsite, 
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
           aes(x=SiteID, y = logmean, fill=SiteID)) + 
      geom_line(aes(group=SiteID), color='black') +
      geom_point(size=5, colour='black', pch=21, alpha=0.75) +
      geom_errorbar(aes(ymin=logmean-sd, ymax=logmean+sd, color=SiteID)) +
      labs(x='Element', y='Mean net photon count') + 
      theme_bw() + 
      theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
            strip.text = element_text(size=14),
            panel.border = element_blank(),
            axis.line = element_line(color='black'))
  )
  dev.off()
}

# ---- Exclude columns that have high within-tree variability or that are unrelated to traffic pollution ----
fieldXRF_format[, `:=`(meanCV=mean(cv, na.rm=T),
                       sdCV = sd(cv, na.rm=T)), by=Elem]
excol <- unique(c(fieldXRF_format[meanCV>0.5,unique(Elem)], 
                  fieldXRF_format[!(Elem %in% union(brakelining_elem, tire_elem)),unique(Elem)],
                  'Rh', 'Pd', 'Rb', 'Ag'))
excol2 <- c('As', 'Cd', 'Mo', 'Na', 'Se', 'Si') #Exclude columns that are very far from being normally distributed (mostly 0s in ICP data)

# ---- Univariate flag based on within-tree CV for metals with CV < 0.5----
CV_flag <- dcast(fieldXRF_format[cv>(meanCV+2*sdCV),], SiteID+Pair~Elem, value.var='cv') 
ggplot(fieldXRF_format, aes(x=mean, y=cv)) + geom_point() +scale_x_sqrt()
NAcount <- function(x) length(which(!is.na(unlist(x))))
CV_flag[, flag_count := NAcount(.SD)-2,
        by = seq_len(nrow(CV_flag)), 
        .SDcols = !(colnames(CV_flag) %in% excol)]

### ---- Multivariate outlier detection ---- 
"see https://stats.stackexchange.com/questions/213/what-is-the-best-way-to-identify-outliers-in-multivariate-data 
for reference discussion as well as https://rpubs.com/Treegonaut/301942
For other resources, see Research_resources/statistics/outliers"
# ---- Compare classic vs robust Mahalanobis distance ----
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

# ---- Filzmoser et al. multivariate outlier detection ----
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

# ---- Check position of multivariate outliers in robust PCA with rrcov (see 4.2 p24 of Todorov and Filzmoser 2009) ---- 
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


outlier_rpca<- cbind(outlierdat_dist, rpca$scores)
setnames(outlier_rpca, colnames(rpca$scores), paste0(colnames(rpca$scores), '_scores'))
rpca_load <- as.data.table(rpca$loadings) 
rpca_load <- rpca_load[,(.SD)*rpca$eigenvalues, by = seq_len(nrow(rpca_load))]
setnames(rpca_load, colnames(rpca_load), paste0(colnames(rpca_load), '_loadings'))

pcs <- c('PC1', 'PC2')
ggplot(data=outlier_rpca, aes_string(x=paste0(pcs[1],'_scores'), y=paste0(pcs[2],'_scores'))) + 
  geom_text(aes(color=aqplot_outliers, label=paste0(SiteID, Pair))) +
  geom_segment(data=rpca_load, x=0, y=0, 
               aes_string(xend=paste0(pcs[1],'_loadings'), yend=paste0(pcs[2],'_loadings')),
               arrow = arrow(length = unit(0.1,"cm")))+
  geom_text(data=rpca_load, aes_string(x=paste0(pcs[1],'_loadings'), y=paste0(pcs[2],'_loadings')),label=rownames(rpca$loadings)) +
  theme_classic()


### ---- Outliers in relationship between field XRF and (lab XRF | ICP) ----
# ---- Multivariate relationship to ICP-OES for the purpose of outlier detection ----
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

# ---- Multivariate relationship to lab XRF for the purpose of outlier detection ----
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

# ---- Univariate relationship to ICP-OES for the purpose of outlier detection ----
#Outlier diagnostic plots
ICPfieldmerge[, ICPfield_flags := 0]
for (chem in unique(ICPfieldmerge$Elem)) {
  print(chem)
  ICPXRF_Zn <- lm(ICP ~ mean, data = ICPfieldmerge[Elem == chem,])
  ggsave(file.path(inspectdir, paste0('fieldXRFICP_regoutliers', chem, '.png')),
         regdiagnostic_custom(ICPXRF_Zn, ICPfieldmerge, chem,  flagcol = 'ICPfield_flags'),
         width = 20, height=12, units='in', dpi=300)
}

#Plot field XRF ~ ICP data
ICPfield_plot <- ggplot(ICPfieldmerge[!(is.na(ICP) |  
                                          Elem %in% c(excol, excol2)),], 
                        aes(x=mean, y=ICP, color=factor(ICPfield_flags))) + 
  geom_linerangeh(aes(xmin=mean-sd, xmax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  geom_smooth(aes(group=1),method='lm') +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
ICPfield_plot

# ---- Univariate relationship to lab XRF for the purpose of outlier detection ----
#Outlier diagnostic plots
ICPfieldmerge[, ICPfield_flags := 0]
for (chem in unique(ICPfieldmerge$Elem)) {
  print(chem)
  ICPXRF_Zn <- lm(ICP ~ mean, data = ICPfieldmerge[Elem == chem,])
  ggsave(file.path(inspectdir, paste0('fieldXRFICP_regoutliers', chem, '.png')),
         regdiagnostic_custom(ICPXRF_Zn, ICPfieldmerge, chem,  flagcol = 'ICPfield_flags'),
         width = 20, height=12, units='in', dpi=300)
}

#Plot field XRF ~ ICP data
ICPfield_plot <- ggplot(ICPfieldmerge[!(is.na(ICP) |  
                                          Elem %in% c(excol, excol2)),], 
                        aes(x=mean, y=ICP, color=factor(ICPfield_flags))) + 
  geom_linerangeh(aes(xmin=mean-sd, xmax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  geom_smooth(aes(group=1),method='lm') +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
ICPfield_plot




### ---- Field XRF - Pollution predictors outliers  ---- 
#Format data
trees_fieldxrf_melt <- melt(trees_fieldxrf[,-c(excol, excol2), with=F], 
                      id.vars=colnames(trees_fieldxrf)[!(colnames(trees_fieldxrf) %in% periodicTable$symb)], 
                      variable.name = 'Elem') %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

fieldrdaformat <- as.matrix(trees_fieldxrf[, setdiff(
  colnames(trees_fieldxrf)[colnames(trees_fieldxrf) %in% periodicTable$symb], c(excol,excol2)), with=F])
pollutrdaformat <- as.matrix(trees_fieldxrf[,grep('heat', colnames(trees_fieldxrf), value=T), with=F])

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
  SiteIDPair = rep(trees_fieldxrf[,paste0(SiteID, Pair)], ncol(fieldrdaformat)),
  elem = rep(colnames(fieldrdaformat), nrow(fieldrdaformat), 'each'),
  fitted = as.vector(sweep(fitted(rda_pollutfield, type = "working"), 2, sigma(rda_pollutfield), "/")),
  residuals = as.vector(rstudent(rda_pollutfield)))

ggplot(rda_fitresid, aes(x=fitted, y=residuals, label=paste0(SiteIDPair,'-', elem))) + 
  geom_text(position=position_jitter(width=1,height=0)) + 
  geom_abline(intercept=c(-2, 2), slope=0)

meanresid <- rda_fitresid[, mean(residuals), by=SiteIDPair]

## ---- Compile flags ----
CV_flag <- dcast(fieldXRF_format[cv>(meanCV+2*sdCV),], SiteID+Pair~Elem, value.var='cv') 
ggplot(fieldXRF_format, aes(x=mean, y=cv)) + geom_point() +scale_x_sqrt()


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


#----
#Relate XRF field data to ICP results
#----

#Plot comparison 
outlier_sites <- c() 
meanresid

ICPfield_plot <- ggplot(ICPfieldmerge[!(is.na(ICP) | 
                                          ICPfieldmerge$SiteID %in% outlier_sites | 
                                          Elem %in% c('Na', 'Si', 'P')),], 
                        aes(x=logmean, y=ICP)) + 
  geom_point(aes(label = paste0(SiteID,Pair), color = cv), size=3, alpha=1/2) +
  scale_color_distiller(palette= 'Spectral') +
  geom_smooth() +
  #geom_linerangeh(aes(xmin=mean-sd, xmax=mean+sd)) + 
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
#plotly_json(ICPfield_plot)

ICPfield_plot
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



#Plot speed limit index vs. XRF data
ggplot(trees_labxrf_melt, aes(x=spdlm, y=value)) + 
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
           data=trees_fieldxrf) #[!(trees_fieldxrf$SiteID %in% c('19','23')),]
summary(Znlm)
ggplotly(qplot(x=predict(Znlm), y=trees_fieldxrf$Zn, color=trees_fieldxrf$heatbustransitlog200, label=as.character(trees_fieldxrf$SiteID)) + geom_smooth(method='lm') )
#+ geom_text(aes(label=trees_fieldxrf$SiteID))

Felm <- lm(Fe~heat_binglog300_1 + sqrt(heatOSMAADTlog300) + sqrt(heatbustransitlog200), data=trees_fieldxrf)
summary(Felm)
qplot(x=predict(Felm), y=trees_fieldxrf$Fe) + geom_smooth(method='lm')

Culm <- lm(Cu~sqrt(AADT), data=trees_fieldxrf)
summary(Culm)
qplot(x=predict(Culm), y=trees_fieldxrf$Cu) + geom_smooth(method='lm') 

Pblm <- lm(Pb~bing+sqrt(AADT)+spdlm, data=trees_fieldxrf)
summary(Pblm)
ggplotly(qplot(x=predict(Pblm), y=trees_fieldxrf$Pb, label=trees_fieldxrf$SiteID) + geom_smooth(method='lm'))

############################################################################################################################################
# Relate lab XRF data to pollution predictors 
############################################################################################################################################
#Plot Bing Congestion index vs. XRF data
bing_xrf <- ggplot(trees_fieldxrf_melt, aes(x=heat_binglog500, y=value, label = SiteID)) + 
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

###########################################################################################
colstodelete <- colnames(trees_labxrf_melt) %in% c(periodicTable$symb, 'value', 'valuestand', 'valuefact')
treeslab <- merge(setDT(trees_labxrf_melt)[, .SD, .SDcols = !colstodelete], labXRF_format, by = c('SiteID', 'Pair', 'Elem'))

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
           data=treeslab[Elem=='Zn' & !(treeslab$SiteID %in% c('44','49','51','52','53','54')),]) #[trees_fieldxrf$SiteID != '51',]
summary(Znlm)
ggplotly(qplot(x=predict(Znlm),
               y=as.list(treeslab[Elem=='Zn' &  !(treeslab$SiteID %in% c('49','51','52','53','54','44')),])$value, 
               label=as.list(treeslab[Elem=='Zn' &  !(treeslab$SiteID %in% c('49','51','52','53','54','44')),])$SiteID) + geom_smooth(method='lm') )


#--------------------------------------------------------------------------------------
# Relate ICP data to pollution predictors 
treesICP <- merge(setDT(trees_labxrf_melt)[, .SD, .SDcols = !colstodelete], ICPmelt, by = c('SiteID', 'Pair', 'Elem'))

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

Znlm <- lm(ICP~bing+sqrt(AADT), data=treesICP[Elem=='Zn' & SiteID != 49,]) #[trees_fieldxrf$SiteID != '51',]
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
  bing_lm <- summary(lm(as.formula(paste(x, "~bing")), data=trees_fieldxrf))
  aadt_lm <- summary(lm(as.formula(paste(x, "~sqrt(AADT)")), data=trees_fieldxrf))
  spdlm_lm <- summary(lm(as.formula(paste(x, "~sqrt(spdlm)")), data=trees_fieldxrf))
  
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


