rootdir <- "C:/Mathis/ICSL/stormwater"
datadir <- file.path(rootdir,'data')
resdir <- file.path(rootdir,'results')
bindir <- file.path(rootdir,'bin')
cloudcaldir <- file.path(bindir,'CloudCal/')

source(file.path(cloudcaldir,'global.R'))

#From CloudCal
k.lines <- read.csv(file=file.path(cloudcaldir, "data/K Line-Table 1.csv"), sep=",")
l.lines <- read.csv(file=file.path(cloudcaldir,"data/L Line-Table 1.csv"), sep=",")
fluorescence.lines <- read.csv(file.path(cloudcaldir,"data/FluorescenceLines.csv"))

#Read PDZ files
readPDZ24Data<- function(filepath, filename){
  
  filename <- gsub(".pdz", "", filename)
  filename.vector <- rep(filename, 2020)
  
  nbrOfRecords <- 2020
  integers <- readPDZ24(filepath, start=357, size=nbrOfRecords)
  sequence <- seq(1, length(integers), 1)
  
  time.est <- integers[21]
  
  channels <- sequence
  energy <- sequence*.02
  counts <- integers/(integers[21]/10)
  
  data.frame(Energy=energy, CPS=counts, Spectrum=filename.vector)
  
}
filepath <- file.path(datadir,'XRF20180808/Postcalibration/ANALYZE_EMP-47.pdz')
filename <- 'ANALYZE_EMP-47.pdz'



