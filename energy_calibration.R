
#Only works for Tracer III SD data (PDZ v24 file format)
library(plyr)
library(ggplot2)
library(gridExtra)

#To do: 
#Install packrat and snapshot package configuration

##----------------------------------------------
#Define directory structure
rootdir <- "C:/Mathis/ICSL/stormwater"
datadir <- file.path(rootdir,'data')
resdir <- file.path(rootdir,'results')
bindir <- file.path(rootdir,'bin')
cloudcaldir <- file.path(bindir,'CloudCal/')

##----------------------------------------------
#Get data and functions from lee drake's CloudCal
globalR <- file.path(cloudcaldir, 'global.R')
if (!file.exists(globalR)) {
  warning("Download and source CloudCal",immediate.=T)
  gitrep <- "https://github.com/leedrake5/CloudCal"
  setwd(bindir)
  system(paste("git clone", gitrep)) #Clone github repository
  setwd(cloudcaldir)
  system("git status")
  source(file.path(cloudcaldir,'global.R'))
  setwd(resdir)
  options(digits=8)
  options(warn=0)
} else {
  warning("CloudCal already exists",immediate.=T)
}

#Run cloudcal locally
#shiny::runGitHub("leedrake5/CloudCal")

##----------------------------------------------
#Read PDZ files

#Redefine CloudCal's readPDZ24Data to correct for integer scaling issue
readPDZ24Data<- function(filepath, filename, evch=0.02){
  
  filename <- gsub(".pdz", "", filename)
  filename.vector <- rep(filename, 2048)
  
  nbrOfRecords <- 2048
  integers <- readPDZ24(filepath, start=357, size=nbrOfRecords)
  try <- data.frame(integers)
  sequence <- seq(1, length(integers), 1)
  
  #time.est <- integers[21] commented out from LD version
  
  channels <- sequence
  energy <- sequence*evch
  counts <- integers/min(integers[integers>0]) #Used to be integers/(integers[21]/10)
  
  data.frame(Channels=channels, Energy=energy, CPS=counts, Spectrum=filename.vector) #Added channels
}

##----------------------------------------------
#Recalibrate
#Find peak in empirical spectrum for a given element and electron transition. Range defines the area around the peak that is inspected
findpeak <- function(spec, elem, siegb,range) {
  lines <- fluorescence.lines[fluorescence.lines$Symbol==elem,siegb]
  subindx <- which.max(spec[spec$Energy > lines-range/2 &
                                  spec$Energy < lines+range/2,'CPS'])
  spec[spec$Energy > lines-range/2 &
                     spec$Energy < lines+range/2,'Channels'][subindx]
}

#Given a spectrum and data frame of elements, find peaks and output data frame
#with a row with the theoretical emitted energy and corresponding channel number for each element transition 
energycal <- function(spec, df) {
  #Get transition energy and channel number corresponding to each transition in spectrum
  adply(df, 1, function(x) {
    c(chpeak = findpeak(spec,x[,'elem'], x[,'siegb'], x[,'rsearch']),
      evpeak = fluorescence.lines[fluorescence.lines$Symbol==x[,'elem'],x[,'siegb']])
    })
}


#Compute evch and plot calibration curve and spectra
modplotcal <- function(spec, df, filepath, filename) {
  #Run linear regression
  evch_cor <- lm(evpeak~chpeak, data=caldf)
  evch_mod <- list(a = as.character(format(coef(evch_cor)[1], digits = 2)),
                   b = as.character(format(coef(evch_cor)[2], digits = 3)), 
                   r2 = as.character(format(summary(evch_cor)$r.squared, digits = 4)))
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,evch_mod)
  eq_exp <- as.character(as.expression(eq))
  
  speccor <- readPDZ24Data(filepath, filename, evch=evch_cor$coefficients[2])
  
  #Plot calibration regression
  calibrationplot <- ggplotGrob(ggplot(caldf, aes(x=chpeak, y=evpeak)) +
    geom_point(size=2)+
    geom_smooth(method='lm')+
    geom_text(x = min(caldf$chpeak)+diff(range(caldf$chpeak))/3, y =max(caldf$evpeak*0.8), label = eq_exp, parse = TRUE) +
    geom_text(aes(label=elem),hjust=0, vjust=0)+
    labs(x='Channel number', y="Observed peak's energy")+
    theme_classic())
  #Plot recalibrated spectrum
  spectraplot <- ggplot(speccor, aes(x=Energy, y=CPS))+
    geom_line(size=1.4) +
    geom_line(data=spec, color='grey') + 
    geom_vline(xintercept=caldf$evpeak, color='red') +
    annotate('text',label=caldf$elem, x=caldf$evpeak+0.1, y=0.8*max(speccor$CPS), color='darkred')+ 
    annotation_custom(grob = calibrationplot, xmin = 30, xmax = max(spec$Energy), ymin = 0.5*max(spec$CPS), ymax = max(spec$CPS)) +
    theme_bw()
  
  print(spectraplot)
  
  return(evch_cor)
}

#Write spectrum to ARTAX-compatible .txt with updated evch
writecal <- function(spec, evch_cor, filepath) {
  txtheader <- c("BeginHeader", "Elin=20 Eabs=0", "Fano=0.11 FWHM=150", "EndHeader")
  txtheader[2] <- paste0("Elin=",evch_cor$coefficients[2],' Eabs=0')
  outpath <- paste0(substr(filepath,1,nchar(filepath)-4),'_edit.txt')
  warning(paste('Writing calibrated spectrum to', outpath))
  writeLines(c(txtheader, spec$CPS), outpath)
}


##----------------------------------------------
#Run
#DF of elements and energy transitions to use in calibration
caldf <- data.frame(elem = as.character(c('Fe','Fe','Cu','Zn','Sr','Rh','Pd')),
                    siegb = as.character(c('Ka1','Kb1','Ka1','Ka1','Ka1','Ka1','Ka1')),
                    rsearch=c(0.4,0.3,0.3,0.3,1,0.6,1),
                    stringsAsFactors=FALSE)

num=12
pdz_filepath <- file.path(datadir,paste0('XRF20180808/Original/ANALYZE_EMP-',num,'.pdz'))
pdz_filename <- paste0('ANALYZE_EMP-',num,'.pdz')
spec <- readPDZ24Data(pdz_filepath,pdz_filename)

caldf <- energycal(spec, caldf)
evch_mod <- modplotcal(spec, caldf, pdz_filepath, pdz_filename)
writecal(spec, evch_cor, pdz_filepath)


