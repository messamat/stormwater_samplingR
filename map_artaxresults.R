library(data.table)
library(plyr)
library(stringr)
library(foreign)

rootdir= "C:/Mathis/ICSL/stormwater" #UPDATE
setwd(file.path(rootdir, "results")) 
datadir = file.path(rootdir, "data")

#Import and format field data TO DO: make faster
fielddata <- read.csv(file.path(datadir,"field_data/field_data_raw_20180808_edit.csv"))
colnames(fielddata)[1] <- 'SiteID'
fielddata$SiteID <- as.character(fielddata$SiteID)
sel <- fielddata[!is.na(fielddata$XRFmin) & !is.na(fielddata$SiteID),]
fieldata_format <- data.frame()
for (row in seq(1,nrow(sel))) {
  extract <- sel[row,]
  for (xrf in seq(sel[row,'XRFmin'], sel[row,'XRFmax'])){
    #print(xrf)
    extract$XRFID <- xrf
    fieldata_format <- rbind(fieldata_format, extract)
  }
}

#Import and format ARTAX deconvolution results
artaxdir <- file.path(datadir, 'XRF20180808/PostBrukerCalibration/deconvolutionresults_XRF12_245_20180827')
tablist <- list.files(artaxdir)
for (i in 1:length(tablist)) {
  res <- read.csv(file.path(artaxdir,tablist[i]))
  res$XRFID <- str_extract(substr(tablist[i],41,43), "\\d+")
  if (i==1){
    deconvres <- res
  } else{
    deconvres <- rbind(deconvres, res)
  }
}
deconvrescast <- dcast(deconvres, XRFID~Element, value.var='Net', fun.aggregate=sum)

#Merge datasets
df <- merge(fieldata_format, deconvrescast, by='XRFID')
#Compute average
colnames(df)
artaxmean <- setDT(df)[,lapply(.SD,mean,na.rm=TRUE), by='SiteID', .SDcols=29:51]
#Merge with original field data
fielddata_artaxmean <- merge(fielddata, artaxmean, by='SiteID')
#Write out
write.dbf(fielddata_artaxmean, 'fielddata_artaxmean_20180827.dbf')

