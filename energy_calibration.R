
#Only works for Tracer III SD data (PDZ v24 file format)
library(plyr)
library(ggplot2)

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

#Redefine readPDZ24Data to correct for integer scaling issue
readPDZ24Data<- function(filepath, filename, evch=0.02){
  
  filename <- gsub(".pdz", "", filename)
  filename.vector <- rep(filename, 2048)
  
  nbrOfRecords <- 2048
  integers <- readPDZ24(filepath, start=357, size=nbrOfRecords)
  try <- data.frame(integers)
  sequence <- seq(1, length(integers), 1)
  
  time.est <- integers[11] #Used to be 21
  
  channels <- sequence
  energy <- sequence*evch
  counts <- integers/(integers[11]) #Used to be integers[21]/10
  
  data.frame(Channels=channels, Energy=energy, CPS=counts, Spectrum=filename.vector) #Added channels
  
}
filepath=file.path(datadir,'XRF20180808/Original/ANALYZE_EMP-18.pdz')
filename <- 'ANALYZE_EMP-18.pdz'
spec18 <- readPDZ24Data(file.path(datadir,'XRF20180808/Original/ANALYZE_EMP-18.pdz'),'ANALYZE_EMP-18.pdz')


ggplot(spec18, aes(x=Energy, y=CPS))+
  geom_line() +
  geom_vline(xintercept=c(Fe.lines[['Ka1']],Fe.lines[['Kb1']], Zn.lines[['Ka1']],Pd.lines[['Ka1']],Rh.lines[['Ka1']]), color='red') +
  theme_bw()

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

#DF of elements and energy transitions to use in calibration
caldf <- data.frame(elem = as.character(c('Fe','Fe','Cu','Zn','Pd','Rh')), 
                    siegb = as.character(c('Ka1','Kb1','Ka1','Ka1','Ka1','Ka1')), 
                    rsearch=c(0.4,0.3,0.3,0.3,1,0.6),
                    stringsAsFactors=FALSE)
#Get theoretical transition energies and empirical energy levels
caldf <- energycal(spec18, caldf)
#Run linear regression
evch_cor <- lm(evpeak~chpeak, data=caldf)
eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(evch_cor)[1], digits = 2),
                      b = format(coef(evch_cor)[2], digits = 4), 
                      r2 = format(summary(evch_cor)$r.squared, digits = 4)))
eq_exp <- as.character(as.expression(eq))

spec18cor <- readPDZ24Data(file.path(datadir,'XRF20180808/Original/ANALYZE_EMP-18.pdz'),'ANALYZE_EMP-18.pdz', evch=evch_cor$coefficients[2])

#Plot calibration regression
calibrationplot <- ggplot(caldf, aes(x=chpeak, y=evpeak)) +
  geom_point(size=2)+
  geom_smooth(method='lm')+
  geom_text(x = min(caldf$chpeak)+diff(range(caldf$chpeak))/4, y =max(caldf$evpeak*0.8), label = eq_exp, parse = TRUE) +
  theme_classic()
#Plot recalibrated spectrum
spectraplot <- ggplot(spec18cor, aes(x=Energy, y=CPS))+
  geom_line() +
  geom_line(data=spec18, color='lightgrey') + 
  geom_vline(xintercept=caldf$evpeak, color='red') +
  annotate('text',label=caldf$elem, x=caldf$evpeak+0.1, y=1500, color='darkred')+ 
  theme_bw()

png()
