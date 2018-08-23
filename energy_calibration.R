rootdir <- "C:/Mathis/ICSL/stormwater"
datadir <- file.path(rootdir,'data')
resdir <- file.path(rootdir,'results')
bindir <- file.path(rootdir,'bin')
cloudcaldir <- file.path(bindir,'CloudCal/')

gitrep <- "https://github.com/leedrake5/CloudCal"
system(paste("cd", bindir))
system(paste("git clone", gitrep))
system("cd CloudCal")
system("git status")

#First download github repository to bin folder
setwd(cloudcaldir)
source(file.path(cloudcaldir,'global.R'))
setwd(resdir)


#Read PDZ files
filepath <- file.path(datadir,'XRF20180808/Original/ANALYZE_EMP-12.pdz')
filename <- 'ANALYZE_EMP-12.pdz'

#Redefine readPDZ24Data to correct for integer
readPDZ24Data<- function(filepath, filename){
  
  filename <- gsub(".pdz", "", filename)
  filename.vector <- rep(filename, 2048)
  
  nbrOfRecords <- 2048
  integers <- readPDZ24(filepath, start=357, size=nbrOfRecords)
  try <- data.frame(integers)
  sequence <- seq(1, length(integers), 1)
  
  time.est <- integers[18] #Used to be 21
  
  channels <- sequence
  energy <- sequence*.02
  counts <- integers/(integers[18]) #Used to be integers[21]/10
  
  data.frame(Energy=energy, CPS=counts, Spectrum=filename.vector)
  
}


spec12 <- readPDZData(filepath,filename)


