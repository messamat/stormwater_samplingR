library(data.table)
library(openxlsx)
library(plyr)

rootdir= "C:/Mathis/ICSL/stormwater" #UPDATE
setwd(file.path(rootdir, "results")) 
datadir = file.path(rootdir, "data")

#Import and format field data
fielddata <- read.xlsx(file.path(datadir,"field_data/field_data_raw_20180719.xlsx"), sheet=1,startRow=1, detectDates=T)
sel <- fielddata[!is.na(fielddata$XRFmin),]
fieldata_format <- data.frame()
for (row in seq(1,nrow(sel))) {
  extract <- sel[row,]
  for (xrf in seq(sel[row,'XRFmin'], sel[row,'XRFmax'])){
    #print(xrf)
    extract$XRFID <- xrf
    fieldata_format <- rbind(fieldata_format, extract)
  }
}

#Import ARTAX deconvolution results
prelimdata <- read.xlsx(file.path(datadir,'XRF20180717/Recalibrated/XRF12_41_results_20180719.xlsx'), sheet=2,startRow=1, detectDates=T)
TNCtourdata <- read.xlsx(file.path(datadir,'XRF20180717/Recalibrated/boardtour_data/Boardtour_results.xlsx'), sheet=2,startRow=1, detectDates=T)
artaxdata <- rbind(prelimdata, TNCtourdata)
#Normalize by Rh count
artaxdata[,-1] <- artaxdata[,-1]/ artaxdata[["Rh.K12"]]
artaxdata$XRFID <- as.numeric(substr(artaxdata$X1, 13, 14))
#Merge datasets
df <- merge(fieldata_format, artaxdata, by='XRFID')
#Compute average
colnames(df)
artaxmean <- setDT(df)[,lapply(.SD,mean,na.rm=TRUE), by='#', .SDcols=28:50]
#Merge with original field data
fielddata_artaxmean <- merge(fielddata, artaxmean, by='#')
#Write out
write.csv(fielddata_artaxmean, 'fielddata_artaxmean_20180719.csv', row.names=F)
