library(ggplot2)
library(data.table)

#Import raw spectra
spec42 <-read.csv('C:/Mathis/ICSL/stormwater/data/XRF20180717/Recalibrated/boardtour_data/ANALYZE_EMP-42@180718_152613.txt')
spec45 <-read.csv('C:/Mathis/ICSL/stormwater/data/XRF20180717/Recalibrated/boardtour_data/ANALYZE_EMP-45@180718_152613.txt')
specs <- merge(spec42, spec45)
#Format dfs
colnames(specs) <- c('Energy', 'XRF42', 'XRF45')
specs_melt <- melt(setDT(specs), id.vars='Energy')

#Import deconvolution model
spec42_deconv <-read.csv('C:/Mathis/ICSL/stormwater/data/XRF20180717/Recalibrated/boardtour_data/ANALYZE_EMP-42@180718_152613_result.csv')
spec42_deconv$val <- with(spec42_deconv, Net+Backgr.)
sel <- c('Fe','Cu', 'Zn', 'Pb', 'Sr')
spec42_deconv_sel <- spec42_deconv[spec42_deconv$Element %in% sel,]
labels <- data.frame(Element=sel, name=c('Iron','Copper','Zinc','Lead','Strontium'))
spec42_deconv_sel <- merge(spec42_deconv_sel, labels,by='Element')

exspec <- ggplot(specs_melt, aes(x=Energy, y=value, color=variable)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,30), name='Energy (keV)') +
  scale_y_continuous(expand=c(0,0), limits=c(0,3500), name='Photon counts') +
  scale_color_manual(name='Site',values=c('#bf812d','#35978f'), labels=c('15th Ave NW | NW 80th St','Carkeek Park'))+
  geom_line(size=0.8, alpha=0.7) + 
  annotate('text',label=as.character(spec42_deconv_sel$name), x=as.numeric(spec42_deconv_sel$Energy.keV), y=3300, angle=45)+ 
  geom_vline(xintercept = spec42_deconv_sel$Energy.keV, alpha=0.2) + 
  theme_classic() + 
  theme(legend.position=c(0.8,0.8),
        text=element_text(size=13))

pdf('C:/Mathis/ICSL/TNC_boardtour/example_spectra.pdf', width=9, height=6.5)
exspec
dev.off()
