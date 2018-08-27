
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
  warning("CloudCal package is already installed.",immediate.=T)
}

#Run cloudcal locally
#shiny::runGitHub("leedrake5/CloudCal")

##----------------------------------------------
#Read PDZ files

#Redefine CloudCal's readPDZ24Data to correct for integer scaling issue
readPDZ24Data<- function(filepath, filename, evch=0.02){
  
  filename <- gsub(".pdz", "", filename)
  nbrOfRecords <- 2048 #Number of channels 
  filename.vector <- rep(filename, nbrOfRecords)
  
  integers <- readPDZ24(filepath, start=357, size=nbrOfRecords)
  sequence <- seq(0, length(integers)-1, 1) #Used to be seq(1, length(integers), 1)
  
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
energycaldat <- function(spec, df) {
  #Get transition energy and channel number corresponding to each transition in spectrum
  adply(df, 1, function(x) {
    c(chpeak = findpeak(spec,x[,'elem'], x[,'siegb'], x[,'rsearch']),
      evpeak = fluorescence.lines[fluorescence.lines$Symbol==x[,'elem'],x[,'siegb']])
    })
}

#Compute evch and plot calibration curve and spectra
modplotcal <- function(spec, df, filepath, filename, outdir=NULL) {
  #Run linear regression
  evch_cor <- lm(evpeak~chpeak, data=df)
  evch_mod <- list(a = as.character(format(coef(evch_cor)[1], digits = 4)),
                   b = as.character(format(coef(evch_cor)[2], digits = 3)), 
                   r2 = as.character(format(summary(evch_cor)$r.squared, digits = 4)))
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,evch_mod)
  eq_exp <- as.character(as.expression(eq))
  
  speccor <- readPDZ24Data(filepath, filename, evch=evch_cor$coefficients[2])
  
  #Plot calibration regression
  calibrationplot <- ggplotGrob(ggplot(df, aes(x=chpeak, y=evpeak)) +
    geom_point(size=2)+
    geom_smooth(method='lm')+
    geom_text(x = min(df$chpeak)+diff(range(df$chpeak))/3, y =max(df$evpeak*0.8), label = eq_exp, parse = TRUE) +
    geom_text(aes(label=elem),hjust=0.1, vjust=0)+
    labs(x='Channel number', y="Observed peak's energy")+
    theme_classic())
  
  #Plot recalibrated spectrum
  spectraplot <- ggplot(speccor, aes(x=Energy-evch_cor$coefficients[1], y=CPS))+
    geom_line(size=1.3) +
    geom_line(data=spec, color='grey') + 
    geom_vline(xintercept=df$evpeak, color='red') +
    scale_x_continuous(name='Energy (Kev)', expand=c(0,0))+
    scale_y_continuous(name='Pulses', expand=c(0,0))+
    annotate('text',label=df$elem, x=df$evpeak+0.1, y=0.8*max(speccor$CPS), color='darkred')+ 
    annotation_custom(grob = calibrationplot, xmin = 30, xmax = max(spec$Energy), ymin = 0.5*max(spec$CPS), ymax = max(spec$CPS)) +
    theme_bw()
  
  if (is.null(outdir)){
    print(spectraplot)
  } else{
    png(file.path(outdir,paste0(substr(filename,1,nchar(filepath)-4),'_plot.png')),width=16, height=8, units='in',res=300)
    print(spectraplot)
    dev.off()
  }
  
  return(evch_cor)
}

#Write spectrum to ARTAX-compatible .txt with updated evch
writecal <- function(spec, evch_cor, filename, outdir) {
  txtheader <- c("BeginHeader", "Elin=20 Eabs=0", "Fano=0.11 FWHM=150", "EndHeader")
  txtheader[2] <- paste0("Elin=",evch_cor$coefficients[2],' Eabs=',evch_cor$coefficients[1])
  outpath <- file.path(outdir, paste0(substr(filename,1,nchar(filename)-4),'_edit.txt'))
  warning(paste('Writing calibrated spectrum to', outpath),immediate.=T)
  writeLines(c(txtheader, spec$CPS), outpath)
}

#Run calibration over all .pdz files in a directory
batchcal <- function(dirpath, df) {
  
  #Get all .pdz files in directory
  speclist  <- list.files(dirpath, pattern='\\.pdz$')
  
  #Directory to write energy calibrated spectra to
  txt_outdir <- file.path(dirpath, paste0('ecalibrated_',format(Sys.time(),'%Y%m%d')))
  if (!file.exists(txt_outdir)){
    warning(paste('Create directory:',txt_outdir),immediate.=T)
    dir.create(txt_outdir)
  }
  
  #Run workflow for each spectrum
  for (i in 1:length(speclist)) {
    pdz_filepath <- file.path(dirpath, speclist[i])
    pdz_filename <- speclist[i]
    spec <- readPDZ24Data(pdz_filepath,pdz_filename) #Read
    caldfext <- energycaldat(spec, df) #Get data for calibration
    evch_mod <- modplotcal(spec, caldfext, pdz_filepath, pdz_filename, txt_outdir) 
    writecal(spec, evch_cor, pdz_filename, txt_outdir)
  }
}

##----------------------------------------------
#Run  for all spectra

#DF of elements and energy transitions to use in calibration
caldf <- data.frame(elem = as.character(c('Fe','Fe','Cu','Zn','Sr','Pd')),
                    siegb = as.character(c('Ka1','Kb1','Ka1','Ka1','Ka1','Ka1')),
                    rsearch=c(0.4,0.3,0.3,0.3,1,1),
                    stringsAsFactors=FALSE)
#Directory to get spectra from
batchcal(file.path(datadir,paste0('XRF20180808/Original_20180718')), caldf)
batchcal(file.path(datadir,paste0('XRF20180808/PostBrukerCalibration')), caldf)


#-------------------------------------------------
#Test method on duplex
duplexcaldf <- data.frame(elem = as.character(c('Cr','Fe','Fe','Mo')),
                    siegb = as.character(c('Ka1','Ka1','Ka1','Ka1')),
                    rsearch=c(0.5,0.5,0.5,0.5),
                    stringsAsFactors=FALSE)

batchcal(file.path(datadir,paste0('XRF20180808/Duplextests')), duplexcaldf)

#Import and format results
duplexspecdir <- file.path(datadir, 'XRF20180808\\Duplextests\\csv')
duplexresdir <- file.path(datadir, 'XRF20180808\\Duplextests\\results')

resfilepathvec <- file.path(duplexresdir,
                            c('duplex_test_20180816_120sec_result.csv',
                              'duplex_test_20180821_2_120sec_result.csv',
                              'duplex_test_20180821_120sec_result.csv',
                              'duplex_test_20180820_120sec_result.csv',
                              'duplex_test_20180822_132sec_result.csv',
                              'duplex_test_20180820_120sec_edit_result.csv',
                              'duplex_test_20180822_132sec_edit_result.csv'))
specfilepathvec <- file.path(duplexspecdir,
                             c('duplex_test_20180816_120sec.csv',
                               'duplex_test_20180821_2_120sec.csv',
                               'duplex_test_20180821_120sec.csv',
                               'duplex_test_20180820_120sec.csv',
                               'duplex_test_20180822_132sec.csv',
                               'duplex_test_20180820_120sec.csv',
                               'duplex_test_20180822_132sec.csv'))
for (i in 1:length(resfilepathvec)){
  spec <- read.csv(specfilepathvec[i])
  livetime <- as.numeric(as.character(spec[rownames(spec)=='Live Time',1]))
  res <- read.csv(resfilepathvec[i])
  res$ID <- substr(resfilepathvec[i], nchar(duplexresdir)+2, nchar(resfilepathvec[i])-4)
  res[,c('Net','Backgr.')] <- res[,c('Net','Backgr.')]/livetime
  if (i==1) {
    resdf <- res
  } else{
    resdf <- rbind(resdf,res)
  }
}

resdf[1:24,'type'] <- 'reference'
resdf[25:40,'type'] <- 'erroneous'
resdf[41:56,'type'] <- 'recalibrated'

fun_mean <- function(x){
  return(data.frame(y=mean(x),label=format(10^mean(x,na.rm=T), digits=4)))}

ggplot(resdf, aes(x=Element,y=Net)) + 
  geom_boxplot(aes(color=type)) + 
  scale_y_log10()+
  stat_summary(fun.data = fun_mean,aes(group=type),geom="text", vjust=-0.7, position=position_dodge(.9)) + 
  theme_bw()+
  theme(text=element_text(size=16))

sdref <- setDT(resdf[resdf$type=='reference'])[,sd(Net)/mean(Net),by=.(Element)]












