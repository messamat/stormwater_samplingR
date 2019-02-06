#Author: Mathis Messager
#Contact information: messamat@uw.edu
#Creation date: July 2018
#Purpose: import, format, merge, and inspect data for pollution modelling project

############################################################################################################################################
#Planning Notes/To do:
# Test influence of method on residuals from relationships
# add spd limit 200 and 300 + gradient 300
# Think about using a percentile for synthetic pollution index to be less dependent on a single value

#Get places with air pollution monitoring or water pollution monitoring. Model for these areas
#Automatize workflow to be able to plug in shapefile and model for these polygons or in the basin of those areas
#Get nationwide AADT and speed limit averages for roads
#Create best model for each variable combination
#Start downloading for Cheasepeake Bay
#Parallelize API download to go much faster so that we can run multiple areas at once?

############################################################################################################################################
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
library(ggpubr)
library(GGally)
library(grid)
library(gridExtra)
library(PeriodicTable)
library(xtable)
library(tools)
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
library(leaflet)
library(devtools)
library(car)
library(AICcmodavg)
library(DescTools)
library(sjPlot)
library(caret)
library(gstat)
library(spdep)
library(ncf)

source_url("https://raw.githubusercontent.com/messamat/biostats_McGarigal/master/biostats.R")
data(periodicTable)

#---- Define directory structure ----
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
datadir <- file.path(rootdir, "data")
inspectdir <- file.path(resdir, "data_inspection")
moddir <- file.path(resdir, "data_modeling")

if (dir.exists(inspectdir)) {
  warning(paste(inspectdir, 'already exists'))
} else {
  warning(paste('Create new directory:',inspectdir))
  dir.create(inspectdir)
}


if (dir.exists(moddir)) {
  warning(paste(moddir, 'already exists'))
} else {
  warning(paste('Create new directory:',moddir))
  dir.create(moddir)
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

chemregdiagnostic_custom <- function(model, df, chem, flagcol, 
                                     CooksOut = FALSE, DFBETAOut = FALSE, tresids = FALSE) {
  "Set of diagnostic tests and plots to detect influential and outlying points.
  Note: modifies df in place by adding flags based on test results
  See https://cran.r-project.org/web/packages/olsrr/vignettes/influence_measures.html for plot examples"
  k <- ols_prep_cdplot_data(model)
  d <- ols_prep_outlier_obs(k)
  if (CooksOut) {
    df[which(df$Elem == chem)[as.data.table(d)[color=='outlier', obs]],
       (flagcol) := get(flagcol) + 1]
  }  
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
  if (DFBETAOut) {
    df[which(df$Elem == chem)[outliernums], (flagcol) := get(flagcol) +1]
  }
  
  #Externaly deleted studentized residuals vs. leverage
  residlevplot <- ols_plot_resid_lev(model)$plot %>%
    delete_layers("GeomText") +
    geom_text(aes(label=df[Elem == chem, paste0(SiteID, Pair)]), size=4) + 
    theme_classic()
  outliernums <- setDT(ols_plot_resid_lev(model)$plot$data)[color %in% c('outlier', 'outlier & leverage'), obs]
  if (tresids) {
    df[which(df$Elem == chem)[outliernums], (flagcol) := get(flagcol) +1]
  }
  
  
  #Straight resid
  residplot <- ols_plot_resid_fit(model) %>%
    delete_layers("GeomText") +
    geom_text(aes(label=df[Elem == chem, paste0(SiteID, Pair)]), size=4) + 
    theme_classic()
  
  #Regression with and without outliers plot
  
  
  if (length(attr(terms(model), 'term.labels')) == 1) { #If the model is a simple linear regression
    regoutlier <- ggplot(df[Elem == chem,],
                         aes_string(x=names(model$model)[2], y=names(model$model)[1])) + 
      geom_text(aes(label=paste0(SiteID, Pair), color=factor(get(flagcol))), size=5) +  
      geom_smooth(method='lm') + 
      labs(x = paste0(names(model$model)[2], chem), y= paste0(names(model$model)[1], chem)) + 
      theme_classic()
  } else {
    regoutlier <- ggplot(df[Elem == chem,],
                         aes(x=fitted(model), y=get(names(model$model)[1]))) + 
      geom_text(aes(label=paste0(SiteID, Pair), color=factor(get(flagcol))), size=5) +  
      geom_smooth(method='lm') + 
      labs(x = paste0('fitted', chem), y= paste0(names(model$model)[1], chem)) + 
      theme_classic()
  }
  return(arrangeGrob(regoutlier, residplot, residlevplot, cooksdchart, dfbeta_mean, dfbeta_int))
}

ols_noplot_resid_lev <- function(mod){
  #Run ols_plot_resid_lev without displaying a plot
  ff <- tempfile()
  png(filename=ff)
  res <- ols_plot_resid_lev(mod)
  dev.off()
  unlink(ff)
  res
}

#mod <- modlistA[[1]]
# for (mod in modlistlogZn) {
#   print(mod$call)
#   regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=1)
# }

regdiagnostic_customtab <- function(mod, maxpar = 20, remove_outliers = 'none', labelvec = NULL,
                                    kCV = TRUE, k=10, cvreps=20) {
  "Function to create the content of a latex formatted table comparing models including.
  - mod (required): lm model or glm model with family(gaussian, link='log') (where formula is fully spelled out with variable names)
  - maxpar: maximum number of parameters in lm model set
  - remove_outliers: whether to run models without outliers (or optionally just without outliers with high leverage)
  - labelvec: vector of IDs for each sample used for displaying which samples were removed as outliers
  - kCV: TRUE/FALSE whether to perform multiple k-fold validation
  - k: number of folds
  - cvreps: # of times to repeat k-fold validation
  "
  
  formcall <- mod$call$formula
  
  "Options for remove_outliers: none, outliers, outliers & leverage"
  
  if (remove_outliers != 'none') {
    if (remove_outliers == 'outliers') {
      outrows <- setDT(ols_prep_rstudlev_data(mod)$`levrstud`)[color %in% c('outlier','outlier & leverage'), obs]
    } 
    if (remove_outliers == 'outliers & leverage') {
      outrows <- setDT(ols_prep_rstudlev_data(mod)$`levrstud`)[color %in% c('outlier & leverage'), obs]
    }
    datasub <- if (length(outrows) > 0) mod$model[-outrows, , drop = FALSE] else mod$model
    
    if (all(class(mod)=='lm')) {
      mod <- lm(formula=as.formula(formcall), data=datasub)
    } 
    if (class(mod)[1] =='glm') {
      mod <- glm(formula=as.formula(formcall), data=datasub)
    }
    outlierlist <- if (!is.null(labelvec)) paste(labelvec[outrows], collapse='\\,') else NULL
  } else {
    outlierlist <- NULL
  }
  
  if (kCV) {
    if (length(mod$coefficients) > 1) {
      ctrl <- trainControl(method = "repeatedcv",
                           number = k,
                           repeats = cvreps)
      
      modcv <- train(form = as.formula(formcall), 
                     data = mod$model,
                     method = class(mod)[1],
                     family = if (class(mod)[1] == 'glm') gaussian(link='log') else NULL,
                     trControl = ctrl)
      
      kCVresults <- setDT(modcv$resample)[, sapply(.SD, function(x) {
        paste0(round(mean(x),2), '\\big(', 
               round(quantile(x, .20),2), '\\textendash', 
               round(quantile(x, .80), 2), '\\big)')}),
        .SDcols= c('Rsquared', 'RMSE', 'MAE')]
      names(kCVresults) <- c('R2cv', 'RMSEcv', 'MAEcv')
    } else {
      kCVresults <- c('R2cv'= '', 'RMSEcv'= '', 'MAEcv'= '')
    }
  } else {
    kCVresults <- NULL
  }
  
  sumod <- summary(mod)
  formula_out <- paste(signif(mod$coefficients, 2), 
                       gsub('heat_*|OSM|log|\\(Intercept\\)', '', names(mod$coefficients)))
  formula_out <- gsub('_','\\\\_', formula_out)
  formula_out <- gsub('\\^','\\\\^', formula_out)
  formula_out <- c(formula_out, rep('', maxpar-length(formula_out)))
  sigcoefs <- which(sumod$coefficients[,'Pr(>|t|)']<0.05)
  formula_out[sigcoefs] <- paste0('\\textbf{', formula_out[sigcoefs], '}')
  names(formula_out) <- paste0('V', 1:length(formula_out))
  
  #Prepare VIF columns
  if (length(mod$coefficients) > 2) {
    vifcols <- c(round(VIF(mod)), rep('', maxpar-(length(mod$coefficients)-1)))
  } else {
    vifcols <- rep('', maxpar)
  }
  names(vifcols) <- paste0('VIF', 2:(maxpar+1))
  

  #Prepare Breusch pagan test
  
  
  c(formula_out,
    nvars = length(unique(unlist(
      strsplit(gsub('\\s', '', names(mod$coefficients)[names(mod$coefficients)!='(Intercept)']),'[*+:]')
    ))), #Number of unique variables used (to check amount of pre-computing required)
    nsig = paste(length(sigcoefs),nrow(sumod$coefficients),sep='/'), #Proportion of significant coefficients
    ShapiroWilkp = round(shapiro.test(residuals(mod))$p.value, 2),
    'BreuschPagan\\_fitp' = if (all(class(mod)=='lm')) round(ols_test_breusch_pagan(mod)$p, 2) else NA, #Test of homoscedasticity of residuals vs. response, assumes indep and normally distributed - H0: variance is constant
    'Score\\_fitp' = if (all(class(mod)=='lm')) round(ols_test_score(mod)$p,2) else NA, #Test of homoscedasticity of residuals vs. fitted, assume iid - Ho: variance is homogenous
    AICc = round(AICcmodavg::AICc(mod)), 
    R2 =  if (all(class(mod)=='lm')) round(sumod$r.squared, 2) else NA, 
    R2adj = if (all(class(mod)=='lm')) round(sumod$adj.r.squared, 2) else NA,
    R2pred = if (all(class(mod)=='lm')) round(ols_pred_rsq(mod), 2) else NA,
    RMSE = round(DescTools::RMSE(mod),2), 
    MAE=round(DescTools::MAE(mod),2),
    kCVresults,
    vifcols,
    outliers = outlierlist)
}

heatlab <- function(x) {
  "quick function to remove keyword from axes"
  gsub('heat_?', '', x)
}

corr_heatmap <- function(xmat=matrix, ymat=NULL, clus=TRUE) {
  "bivariate clustered correlation heatmap from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization"
  cormat <- round(cor(x=xmat, y=ymat),2)
  
  melted_cormat <- melt(cormat)
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    if (is.null(ymat)) {
      dd <- as.dist((1-cormat)/2)
      hc <- hclust(dd)
      cormat <-cormat[hc$order, hc$order]
    } else {
      dd <- dist(cormat)
      hc <- hclust(dd)
      cormat <-cormat[hc$order,]
    }
  }
  # Reorder the correlation matrix
  if (clus) {cormat <- reorder_cormat(cormat)}
  
  if (is.null(ymat)) {
    upper_tri <- get_upper_tri(cormat)
    # Melt the correlation matrix
    melted_cormat <- melt(upper_tri, na.rm = TRUE)
  } else {
    melted_cormat <- melt(cormat, na.rm = TRUE)
  }
  
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed() + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  if(is.null(ymat)) {
    ggheatmap <- ggheatmap + theme(legend.position = c(0.6, 0.7))
  }
  print(ggheatmap)
}

pcabiplot_grid <- function(dt, cols, idcols, nPCs=NULL, scload=1) {
  "Create a plot matrix of pca principal components against each other"
  
  #Run PCA
  pca <- PcaClassic(dt[, cols, with=F], scale=T)
  
  #Add site ID labels
  pca_scores<- cbind(dt[, idcols, with=F], pca$scores)
  pca_scores[, ID := do.call(paste0,.SD), .SDcols=idcols]
  
  #Add chem element label
  pca_load <- as.data.table(pca$loadings) %>%
    .[, Elem := rownames(pca$loadings)]
  
  #Compute segment length based on component's eigenvalue
  pca_load <- pca_load[,(.SD)*pca$eigenvalues, by = Elem]
  
  #If number of PCs to graph was not entered as a parameter, default to all significant PCs
  if (is.null(nPCs)) {
    nPCs <- ncol(pca$scores)
  }
  
  #Iterate over every pair of PCs
  biplots_l <- sapply(seq_len(nPCs-1), function(i) {
    sapply(seq(2, nPCs), function(j) {
      xpc <- paste0('PC', i)
      ypc <- paste0('PC', j)
      
      if (i != j) {
        ggplotGrob(
          ggplot(data=pca_scores, aes_string(x = xpc, y= ypc)) + 
            geom_text(aes(label=ID), alpha=0.5) +
            geom_segment(data=pca_load, x=0, y=0, 
                         aes(xend=scload*get(xpc), 
                             yend=scload*get(ypc)),
                         arrow = arrow(length = unit(0.1,"cm")), color='red')+
            geom_text(data=pca_load, aes(x=scload*get(xpc), y=scload*get(ypc)),
                      label=rownames(pca$loadings), color='red', size=5, alpha=0.8) +
            theme_classic()
        )
      } else {
        #Write proportion of variance explained by PC
        text_grob(round(pca$eigenvalues[i]/sum(pca$eigenvalues), 2))
      }
    })
  })
  #Plot, only keeping unique combinations by grabbing lower triangle of plot matrix
  do.call("grid.arrange", list(grobs=biplots_l[lower.tri(biplots_l, diag=T)], ncol=nPCs-1)) 
}

latex_format <- function(summary_table) {
  "Function to format a DF with embedded latex formatting into a 
  full latex document that can directly turned into a PDF. Can adjust paperheight and paperwodth, margin, etc."
  ltable <- print(xtable(summary_table), type="latex", 
                  booktabs=TRUE, #format table
                  sanitize.text.function = identity, #Disable sanitizing to allow for latex formatting to be kept intact
                  comment = FALSE, #remove comment at beggining of table about xtable production
                  file='')
  #Embed table within a latex document, adjusting width and height, etc.
  ltable_edit <- gsub('\\begin{table}', 
                      paste('\\documentclass{article}',
                            '\\usepackage[paperheight=8.5in,paperwidth=28in,margin=0.1in,headheight=0.0in,footskip=0.5in,includehead,includefoot]{geometry}',
                            '\\usepackage{booktabs}',
                            '\\begin{document}',
                            '\\setlength{\\tabcolsep}{2pt}',
                            '\\begin{table}', sep='\n'), 
                      ltable, fixed=TRUE)
  ltable_edit <- gsub('\\end{table}', 
                      paste('\\end{table}','\\end{document}', sep='\n'),
                      ltable_edit, fixed=TRUE)
  return(ltable_edit)
}

weightmat_IDW <- function(coords, knb, mindist = 10) {
  "Function to build an inverse distance weight matrix, if knb is NULL returns an IDW matrix
  of all points; if knb is defined e.g. 10, returns an IDW matrix of knb neighbors.
  - coords must be a data.frame, data.table, or other two-column set of coordinates
  - mindist is a constant added to the distance matrix for distinct but overlapping points"
  
  spp <- SpatialPoints(coords)
  "See session 1 of Intro to mapping and spatial modelling in R by Richard Harris, creator of package spdep"
  distall_inv <- 1/(mindist+as.matrix(dist(coords)))^0.5 #Make inverse distance weight matrix
  diag(distall_inv) <- 0
  outmat <- mat2listw(distall_inv)
  
  if (!is.null(knb)) {
    knear <- knearneigh(spp, k=knb, RANN=F) #Find 10 nearest neighbors for each point
    dweights <- lapply(1:knear$np, function(i) { #For each point, only keep IDW of 10 nearest neighbors
      distall_inv[i, knear$nn[i,]]
    })
    knear_nb <- knn2nb(knear) #Convert kkn object to a neighbours list of class nb
    #plot(knear_nb, coords)
    outmat <- nb2listw(knear_nb, glist = dweights, style='C') #Create spatial weights matrix (globally standardized i.e. not by neighborhood)
  }
  return(outmat)
}
###########################################################################################################################################

# 1. Import data 
############################################################################################################################################
########### ---- A. Import and format field data ---- ####
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

########### ---- B. Import and format field XRF deconvolution results ---- ####
fieldXRF <- deconvolution_import(file.path(datadir, 'XRF20181210/PostBrukerCalibration/deconvolutionresults_XRF12_245_20181215'),
                                 idstart = 41, idend= 43)

########### ---- C. Import raw XRF lab data for inspection ---- ####
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

########### ---- D. Import and format lab XRF deconvolution results ---- ####
labXRF <- deconvolution_import(file.path(datadir, 'XRF20181210/PelletMeasurements/deconvolutionresults_labXRF_20181215'),
                               idstart = 29, idend= 33)
# ---- Summarize labXRF mean net for each element and Net/Background ratio ----
#Look at determinants of noise
fieldXRFsummary <- fieldXRF[,list(meannet = mean(Net),
                                  noise = mean(Net/Backgr.)), .(Element, Line, Energy.keV)]
keVmaxnet <- unique(fieldXRFsummary) %>% #Only keep Element-Line-Energy characteristics of each element
  setkey(meannet) %>% #Classify from low to high energy 
  .[,.SD[.N], by=Element] #Only keep line with the highest average net photon count
keVmineV <- unique(fieldXRFsummary) %>% #Only keep Element-Line-Energy characteristics of each element
  setkey(Energy.keV)%>% #Classify from low to high energy 
  .[,.SD[1], by=Element] #Only keep line with the lowest average net photon count

ggplot(keVmax, aes(x=log(meannet+1), y=noise)) +
  geom_text(aes(label=Element)) +
  scale_y_log10() +
  geom_smooth(span=1)

ggplot(keVmax, aes(x=Energy.keV, y=noise)) +
  geom_text(aes(label=Element)) +
  scale_y_log10() +
  geom_smooth(span=1)

ggplot(keVmin, aes(x=Energy.keV, y=noise)) +
  geom_text(aes(label=Element)) +
  scale_y_log10() +
  geom_smooth(span=1)

########### ---- E. Import ICP-OES data ---- ####
ICPdat <- read.xlsx(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'),
                    1, stringsAsFactors = F)
ICPthresholds <- read.xlsx(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'),
                           2, stringsAsFactors = F)
rownames(ICPthresholds) <- gsub('\\W', '', ICPthresholds[,1])
ICPthresholds_format <- as.data.table(t(ICPthresholds[,-1])) %>%
  .[, Elem := colnames(ICPthresholds)[-1]]

########### ---- F. Import GIS data (including pollution variables) ---- ####
trees <- as.data.table(readOGR(dsn = file.path(resdir, 'Seattle_sampling.gdb'), layer = 'XRFsites_proj'))
#summary(trees)
heatcols <- colnames(trees)[grep('heat', colnames(trees))]
trees[, (heatcols) := lapply(.SD, function(x){x[is.na(x)] <- 0; x}), .SDcols = heatcols]
colnames(trees)[1] <- 'SiteID'
trees <- trees[!(SiteID %in% c('NA', NA)),]

########### ---- G. Define elements associated with car traffic ---- ####
#Elements present in brake linings and emitted brake dust (from Thorpe et al. 2008)
brakelining_elem <- c('Al', 'As', 'Ba', 'Ca', 'Cd', 'Co' ,'Cr', 'Cu', 'Fe', 'K', 'Li', 'Mg', 
                      'Mn', 'Mo', 'Na', 'Ni', 'Pb', 'Sb', 'Se', 'Sr', 'Zn') 
#Elements present in passenger car tyre tread (from Thorpe et al. 2008)
tire_elem <- c('Al', 'Ba', 'Ca', 'Cd', 'Co', 'Cr', 'Cu', 'Fe', 'K', 'Mg', 'Mn', 'Na', 'Ni',
               'Pb', 'Sb', 'Sr', 'Ti', 'Zn')

############################################################################################################################################

# 2. Format XRF and ICP data
############################################################################################################################################
########### ---- A. Format XRF data ---- ####
# ---- 1. Cast while summing net photon counts across electron transitions ----
fieldXRFcast <- dcast(setDT(fieldXRF), XRFID~Element, value.var='Net', fun.aggregate=sum) 
fieldXRFcast[, XRFID := as.numeric(gsub('[_]', '', XRFID))] #Format site number
fieldXRFcast[fieldXRFcast < 0] <- 0 #Floor negative net photon count to 0

labXRFcast <- dcast(setDT(labXRF), XRFID~Element, value.var='Net', fun.aggregate=sum)
labXRFcast[labXRFcast < 0] <- 0

# ---- 2. Normalize data by Rhodium photon count for field and lab results ----
fieldXRFcastnorm <- fieldXRFcast[, lapply(.SD, function(x) {x/Rh}), by = XRFID]
labXRFcastnorm <- labXRFcast[, lapply(.SD, function(x) {x/Rh}), by = XRFID]

# ---- 3. Merge datasets: lab XRF + field XRF + field variables ----
fieldt <- setDT(fieldata_format)[fieldXRFcastnorm, on='XRFID']

labXRFcastnorm[, `:=`(SiteID = gsub('[A-B].*', '', XRFID),
                      Pair = gsub('[0-9_]+', '', XRFID),
                      XRFID = NULL)]
labdt <- setDT(fieldata_sel)[labXRFcastnorm, on =  .(SiteID, Pair)]

# ---- 4. Compute average, sd, and range for lab XRF results over multiple measurements for a given pellet ----
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

# ---- 5. Compute average, sd, and range, then melt XRF field data ----
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

# ---- 6. Write out mean XRF values (should write this directly to XRF sites proj) ----
field_artaxmean <- field_artaxstats[, .SD, .SDcols = c(1,2, grep('mean', colnames(field_artaxstats)))] 
colnames(field_artaxmean) <- gsub('_mean', '', colnames(field_artaxmean))
write.dbf(field_artaxmean, 'field_artaxmean_20180827.dbf')

lab_artaxmean <- lab_artaxstats[, .SD, .SDcols = c(1,2, grep('mean', colnames(lab_artaxstats)))] 
colnames(lab_artaxmean) <- gsub('_mean', '', colnames(lab_artaxmean))

# ---- 7. Check field XRF data distribution by element then transform and standardize ----
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

# ---- 8. Check lab XRF data distribution by element then transform and standardize ----
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

# ---- 9. Re-cast data for multivariate analysis ----
field_transmean <- dcast(fieldXRF_format, SiteID+Pair~Elem, value.var = 'transmean') 
lab_transmean <- dcast(labXRF_format, SiteID+Pair~Elem, value.var = 'transmean') 

########### ---- B. Format ICP-OES data ---- ####
# ---- 1. Format data and check site overlap between field XRF and ICP----
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

# ---- 2. Melt dataset to merge with XRF data ----
ICPmelt <- melt(setDT(ICPmean), id.vars = 'SAMPLE.SET', variable.name = 'Elem', value.name = 'ICP')
ICPmelt[, `:=`(SiteID = gsub('[A-Z]', '', SAMPLE.SET),
               Pair = gsub('[0-9]', '', SAMPLE.SET))]

# ---- 3. Check ICP data distribution by element then transform and standardize ----
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

# ---- 4. Re-cast data for multivariate analysis ----
ICP_trans <- dcast(ICPmelt, SiteID+Pair~Elem, value.var='transICP')

########### ---- C. Format pollutant data ---- ####
# ---- 1. Check pollutant data distribution by element then transform and standardize ----
#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
treesmelt <- melt(trees[!is.na(SiteID),], id.vars=grep('heat|NLCD', colnames(trees),
                                                       ignore.case=T, value=T, invert=T), 
                  variable.name = 'pollutvar' , value.name = 'pollutvalue')
treesmelt[, transpollut := transformTukey_lambda(pollutvalue, start = -2, end = 2,int = 0.025, 
                                                 rastertab=F,rep=100, verbose = FALSE, 
                                                 statistic = 1)[[1]], by=pollutvar]
#Run a shapiro test on each element
trees_trans <- dcast(treesmelt, SiteID+Pair~pollutvar, value.var = 'transpollut') 

normtest <- shapiro_test_df(trees_trans[,sapply(trees_trans, class)=='numeric' &
                                          sapply(trees_trans, function(x) length(unique(x)))>1, with=F])
normsig <- data.frame(sig = normtest$significance)
normsig$pollutvar <- row.names(normsig)
#Plot histogram of transformed data color-coded by whether normally distributed or not
treesnorm_join <- treesmelt[normsig, on='pollutvar']
png(file.path(inspectdir, paste0('trees_TukeyTransform.png')), width = 20, height=12, units='in', res=300)
ggplot(treesnorm_join, aes(x=transpollut, fill=sig)) + 
  geom_histogram() + 
  facet_wrap(~pollutvar, scales='free') + 
  theme_classic()
dev.off()     

# z-standardize by element
treesmelt[, transpollut:=scale(transpollut), by= pollutvar]

# ---- 2. Re-cast data for multivariate analysis ----
trees_trans <- dcast(treesmelt, formula = SiteID+Pair~pollutvar,  value.var='transpollut')

########### ---- D. Merge datasets ---- ####
# ---- 1. Merge XRF field data with ICP data ----
ICPfieldmerge <- merge(ICPmelt, fieldXRF_format, by = c('SiteID', 'Pair', 'Elem'), all.x=F, all.y=F)
unique(ICPmelt$Elem)[!(unique(ICPmelt$Elem) %in% unique(fieldXRF_format$Elem))] #Check what elems are in vs out
unique(fieldXRF_format$Elem)[!(unique(fieldXRF_format$Elem) %in% unique(ICPmelt$Elem))] #Check what elems are in vs out

# ---- 2. Merge XRF field data with XRF lab data ----
joincols <- c('SiteID', 'Pair', 'Elem')
labfieldmerge <- labXRF_format[fieldXRF_format, on = joincols] %>%
  .[SiteID != 1 & Elem != 'Rh',]
setnames(labfieldmerge, 
         colnames(labfieldmerge[,-joincols, with=F]),
         gsub('^(?!i)', 'lab_', colnames(labfieldmerge[,-joincols, with=F]), perl=T))
setnames(labfieldmerge, 
         colnames(labfieldmerge[,-joincols, with=F]),
         gsub('^i[.]', 'field_', colnames(labfieldmerge[,-joincols, with=F]), perl=T))

# ---- 3. Merge XRF lab data with ICP data ----
ICPlabmerge <- merge(ICPmelt, labXRF_format, by = c('SiteID', 'Pair', 'Elem'), all.x=F, all.y=F)
unique(ICPmelt$Elem)[!(unique(ICPmelt$Elem) %in% unique(labXRF_format$Elem))] #Check what elems are in vs out
unique(labXRF_format$Elem)[!(unique(labXRF_format$Elem) %in% unique(ICPmelt$Elem))] #Check what elems are in vs out

# ---- 4. Join non-transformed/standardized pollutant data to sites ----
pollutfieldmerge <- fieldXRF_format[trees, on=c('SiteID','Pair')]
pollutlabmerge <- labXRF_format[trees, on=c('SiteID','Pair')]
pollutICPmerge <- ICPmelt[trees, on=c('SiteID','Pair')]

# ---- 5. Join transformed/standardized pollutant data to sites ----
trees_fieldxrf_trans <- trees_trans[field_transmean, on=c('SiteID','Pair')][
  !(is.na(Fe) | SiteID==1),]
trees_labxrf_trans <- trees_trans[lab_transmean, on=c('SiteID','Pair')][
  !(is.na(Fe) | SiteID==1),]
trees_ICP_trans <- trees_trans[ICP_trans, on=c('SiteID', 'Pair')][
  !(is.na(Fe) | SiteID==1),]

############################################################################################################################################


# 3. Inspect data and remove outliers
############################################################################################################################################
######################  ---- A. Field XRF ---- ###########
str(fieldXRF_format)
# ---- 1. Assess within-tree variability ----
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
# ggplotly(
#   ggplot(fieldXRF_format[!(fieldXRF_format$Elem %in% c('Rh','Pd','Ar')),], 
#          aes(x=mean, y=cv, color=name, label=paste0(SiteID, Pair))) + 
#     geom_point()+
#     geom_smooth() +
#     labs(x='Mean photon count (normalized)', y='Coefficient of variation') + 
#     facet_wrap(~name, scales='free') +
#     theme_classic() + 
#     theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
#           strip.text = element_text(size=14))
# )
# 
# #
# for (elem in unique(fieldXRF_format[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
#   print(elem)
#   png(file.path(inspectdir, paste0('fieldXRF_withintree_',elem,'.png')), width = 20, height=12, units='in', res=300)
#   print(
#     ggplot(fieldt, 
#            aes(x=paste0(SiteID, Pair), y = get(elem), fill=Pair)) + 
#       geom_line(aes(group=paste0(SiteID, Pair)), color='black') +
#       geom_point(size=5, colour='black', pch=21, alpha=0.75) +
#       labs(x='Site', y='Mean net photon count') + 
#       theme_bw() + 
#       theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
#             strip.text = element_text(size=14),
#             panel.border = element_blank(),
#             axis.line = element_line(color='black'))
#   )
#   dev.off()
# }

# ---- 2. Assess within-site variability ---- 
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

# ---- 3. Exclude columns that have high within-tree variability or that are unrelated to traffic pollution ----
fieldXRF_format[, `:=`(meanCV=mean(cv, na.rm=T),
                       sdCV = sd(cv, na.rm=T)), by=Elem]
excol <- unique(c(fieldXRF_format[meanCV>0.5,unique(Elem)], 
                  fieldXRF_format[!(Elem %in% union(brakelining_elem, tire_elem)),Elem],
                  'Rh', 'Pd', 'Rb', 'Ag')) %>%
  setdiff('Zr') #Keep zirconium as appears to be correlated with traffic
excol2 <- c('As', 'Cd', 'Mo', 'Na', 'Se', 'Si') #Exclude columns that are very far from being normally distributed (mostly 0s in ICP data)

# ---- 4. Univariate flag based on within-tree CV for metals with CV < 0.5----
#NAcount <- function(x) length(which(!is.na(unlist(x))))
fieldXRF_format[(cv>(meanCV+2*sdCV)) & !(Elem %in% excol) & !is.na(cv),
                CVflag_count_field := .N, by=.(SiteID, Pair)][
                  , CVflag_count_field := ifelse(is.na(CVflag_count_field) | is.na(cv), 
                                                 as.integer(round(mean(CVflag_count_field, na.rm=T))),
                                                 CVflag_count_field),
                  by=.(SiteID, Pair)][
                    is.na(CVflag_count_field), CVflag_count_field := 0]

### ---- Multivariate outlier detection ---- 
"see https://stats.stackexchange.com/questions/213/what-is-the-best-way-to-identify-outliers-in-multivariate-data 
for reference discussion as well as https://rpubs.com/Treegonaut/301942
For other resources, see Research_resources/statistics/outliers"
# ---- 5. Compare classic vs robust Mahalanobis distance ----
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

# ---- 6. Filzmoser et al. multivariate outlier detection ----
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

### ---- Outliers in relationship between field XRF and (lab XRF | ICP | pollutant drivers) ----
# ---- 7. Multivariate relationship to ICP-OES for the purpose of outlier detection ----
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

# ---- 8. Multivariate relationship to lab XRF for the purpose of outlier detection ----
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

# ---- 9. Multivariate relationship to pollution predictors for the purpose of outlier detection  ---- 
#Format data
trees_fieldxrf_trans_melt <- melt(trees_fieldxrf_trans[,-c(excol, excol2), with=F], 
                                  id.vars=colnames(trees_fieldxrf_trans)[!(colnames(trees_fieldxrf_trans) %in% periodicTable$symb)], 
                                  variable.name = 'Elem') %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

fieldrdaformat <- as.matrix(trees_fieldxrf_trans[, setdiff(
  colnames(trees_fieldxrf_trans)[colnames(trees_fieldxrf_trans) %in% periodicTable$symb], c(excol,excol2)), with=F])
pollutrdaformat <- as.matrix(trees_fieldxrf_trans[,grep('heat|NLCD', colnames(trees_fieldxrf_trans), value=T, ignore.case = T), 
                                                  with=F])

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
  SiteIDPair = rep(trees_fieldxrf_trans[,paste0(SiteID, Pair)], ncol(fieldrdaformat)),
  elem = rep(colnames(fieldrdaformat), nrow(fieldrdaformat), 'each'),
  fitted = as.vector(sweep(fitted(rda_pollutfield, type = "working"), 2, sigma(rda_pollutfield), "/")),
  residuals = as.vector(rstudent(rda_pollutfield)))

ggplot(rda_fitresid, aes(x=fitted, y=residuals, label=paste0(SiteIDPair,'-', elem))) + 
  geom_text(position=position_jitter(width=1,height=0)) + 
  geom_abline(intercept=c(-2, 2), slope=0)

meanresid <- rda_fitresid[, mean(residuals), by=SiteIDPair]
# ---- 10. Univariate relationship to ICP-OES for the purpose of outlier detection ----
#Outlier diagnostic plots
ICPfieldmerge[, ICPfield_flags := 0]
for (chem in unique(ICPfieldmerge$Elem)) {
  print(chem)
  ICPfield_lm <- lm(ICP ~ mean, data = ICPfieldmerge[Elem == chem,])
  ggsave(file.path(inspectdir, paste0('fieldXRFICP_regoutliers', chem, '.png')),
         chemregdiagnostic_custom(ICPfield_lm, ICPfieldmerge, chem,  flagcol = 'ICPfield_flags', tresids = TRUE),
         width = 20, height=12, units='in', dpi=300)
  ICPfieldmerge[Elem == chem, ICPfieldR2 := summary(ICPfield_lm)$adj.r.squared]
}

#Plot field XRF ~ ICP data
ICPfield_plot <- ggplot(ICPfieldmerge[!(is.na(ICP) |  
                                          Elem %in% c(excol, excol2)),], 
                        aes(x=mean, y=ICP, color=factor(ICPfield_flags), group=1)) + 
  geom_linerangeh(aes(xmin=mean-sd, xmax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
ICPfield_plot


#Investigate determinants of R2
ICPmeanR2 <- unique(ICPfieldmerge[, list(meanXRF = mean(mean, na.rm=T),
                                         meanICP = mean(ICP, na.rm=T),
                                         R2 = ICPfieldR2), by=Elem][
                                           keVmaxnet, on='Elem==Element'])
ggplot(ICPmeanR2, aes(x=meanXRF, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=meanICP, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=noise, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

summary(loess(R2~meanXRF, data=ICPmeanR2, span=1))
summary(loess(R2~meanXRF+noise, data=ICPmeanR2, span=1))

# ---- 11. Univariate relationship to lab XRF for the purpose of outlier detection ----
#Outlier diagnostic plots
labfieldmerge[, labfield_flags := 0]
for (chem in unique(labfieldmerge$Elem)) {
  print(chem)
  labfield_lm <- lm(lab_mean ~ field_mean, data = labfieldmerge[Elem == chem,])
  ggsave(file.path(inspectdir, paste0('fieldlab_regoutliers', chem, '.png')),
         chemregdiagnostic_custom(labfield_lm, labfieldmerge, chem,  flagcol = 'labfield_flags', tresids = TRUE),
         width = 20, height=12, units='in', dpi=300)
  labfieldmerge[Elem == chem, labfieldR2 := summary(labfield_lm)$adj.r.squared]
}

#Plot field XRF ~ lab data
labfield_plot <- ggplot(labfieldmerge[!(is.na(lab_mean) |  
                                          Elem %in% c(excol, excol2)),], 
                        aes(x=field_mean, y=lab_mean, color=factor(labfield_flags), group=1)) + 
  geom_linerangeh(aes(xmin=field_mean-field_sd, xmax=field_mean+field_sd), alpha=1/2) + 
  geom_segment(aes(xend=field_mean, y=lab_mean-lab_sd, yend=lab_mean+lab_sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
labfield_plot

#Look at determinants of R2
labmeanR2 <- unique(labfieldmerge[, list(meanfield = mean(field_mean, na.rm=T),
                                         meanlab = mean(lab_mean, na.rm=T),
                                         R2 = labfieldR2), by=Elem])[
                                           keVmaxnet, on='Elem==Element']

ggplot(labmeanR2, aes(x=meanfield, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(labmeanR2, aes(x=meanlab, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(labmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(labmeanR2, aes(x=noise, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

# ---- 12. Univariate relationship to pollution predictors for the purpose of outlier detection ----
#Outlier diagnostic plots
pollutfieldmerge[, pollutfield_flags := 0]
for (chem in unique(pollutfieldmerge$Elem)) {
  print(chem)
  pollutfield_lm <- lm(mean ~ heat_binglog300*heatOSMAADTlog300*heatOSMgradientlog200*nlcd_imp_ps +
                         heatbustransitlog300 + heatbustransitlog300:heatOSMgradientlog200, 
                       data = pollutfieldmerge[Elem == chem,])
  ggsave(file.path(inspectdir, paste0('fieldpollution_regoutliers2', chem, '.png')),
         chemregdiagnostic_custom(pollutfield_lm, pollutfieldmerge, chem,  flagcol = 'pollutfield_flags', tresids = TRUE),
         width = 20, height=12, units='in', dpi=300)
  pollutfieldmerge[Elem == chem, `:=`(pollutpred = fitted(pollutfield_lm),
                                      pollutfieldR2 = summary(pollutfield_lm)$adj.r.squared)]
}

#Plot predicted field XRF ~ observed field XRF
pollutfield_plot <- ggplot(pollutfieldmerge[!(Elem %in% c(excol, excol2, NA)),], 
                           aes(x=pollutpred, y=mean, group=1, color=factor(pollutfield_flags))) + 
  geom_linerange(aes(ymin=mean-sd, ymax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  geom_smooth(method='lm') +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
pollutfield_plot


ggsave(file.path(inspectdir, 'fieldpollut_all.png'),
       pollutfield_plot,
       width = 20, height=12, units='in', dpi=600)

#Investigate determinants of R2
pollutmeanR2 <- unique(pollutfieldmerge[, list(meanXRF = mean(mean, na.rm=T),
                                               R2 = pollutfieldR2), by=Elem][
                                                 keVmaxnet, on='Elem==Element'])
ggplot(pollutmeanR2, aes(x=meanXRF, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=noise, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

### ---- 13. Compile outlier flags ----
#Merge all flag datasets
fieldXRF_formatflags <- merge(fieldXRF_format, ICPfieldmerge[, .(SiteID, Pair, Elem, ICPfield_flags, ICPfieldR2)], 
                              by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(labfieldmerge[, .(SiteID, Pair, Elem, labfield_flags, labfieldR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(pollutfieldmerge[, .(SiteID, Pair, Elem, pollutfield_flags, pollutfieldR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T)

ggplot(fieldXRF_formatflags, aes(x=pollutfieldR2, y=ICPfieldR2, label=Elem)) + 
  geom_text()
ggplot(fieldXRF_formatflags, aes(x=pollutfieldR2, y=labfieldR2, label=Elem)) + 
  geom_text()

#Compute total number of flags for each tree across all elements that have an R2 of at least 0.10 across models
#Don't want to count outliers from spurious models
flagelems <- fieldXRF_formatflags[, min(ICPfieldR2, labfieldR2, pollutfieldR2, na.rm=T)>0.20, by=Elem][
  V1==TRUE & Elem != 'Rh', Elem]
fieldXRF_formatflags[Elem %in% flagelems,
                     flagpartsum_field := sum(ICPfield_flags, labfield_flags, pollutfield_flags, na.rm=T), 
                     by=.(SiteID, Pair)][
                       , flagsum_field := sum(flagpartsum_field, CVflag_count_field, na.rm=T), by=.(SiteID, Pair)]

#Check total number of flags by each type of flag
fieldXRF_formatflags[Elem %in% flagelems,
                     `:=`(ICPfield_flags_sum = sum(ICPfield_flags, na.rm=T),
                          labfield_flags_sum = sum(labfield_flags, na.rm=T),
                          pollutfield_flags_sum = sum(pollutfield_flags, na.rm=T)),
                     by=.(SiteID, Pair)]

fieldXRF_formatflags_u <- unique(
  fieldXRF_formatflags[!is.na(flagpartsum_field),
                       .SD, 
                       .SDcols=c('SiteID', 'Pair', 'CVflag_count_field',
                                 grep('flag.*sum', colnames(fieldXRF_formatflags), value=T))])

ggplot(fieldXRF_formatflags_u, aes(x=ICPfield_flags_sum, y=labfield_flags_sum, 
                                   label=paste0(SiteID, Pair), color = CVflag_count_field)) + 
  scale_color_distiller(palette='YlGnBu') +
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5) 

ggplot(fieldXRF_formatflags_u, aes(x=ICPfield_flags_sum, y=pollutfield_flags_sum, 
                                   label=paste0(SiteID, Pair), color = CVflag_count_field)) + 
  scale_color_distiller(palette='YlGnBu') + 
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5)

#Records to inspect
flagelems
fieldXRF_inspect <- unique(fieldXRF_formatflags[(ICPfield_flags_sum>5) | 
                                                  (labfield_flags_sum > 5) | 
                                                  (pollutfield_flags_sum > 5), 
                                                .(SiteID, Pair)])
outliertrees <- trees[fieldXRF_inspect, on=c('SiteID', 'Pair')]
outlierlocs <- SpatialPointsDataFrame(coords = data.frame(outliertrees$POINT_X, outliertrees$POINT_Y),
                                      data= as.data.frame(outliertrees))
#View(outlierlocs@data)
leaflet(data = outlierlocs) %>% addTiles() %>%
  addMarkers(clusterOptions = markerClusterOptions(),
             popup = ~paste0(SiteID, Pair))


###################### ---- B. Lab XRF ---- ####
# ---- 1. Assess within-pellet variability ----
#Plot coefficient of variation distributions for every element
cvmean_labeltree <- as.data.frame(labXRF_format[!(labXRF_format$Elem %in% c('Rh','Pd','Ar')),
                                                paste0('Mean CV: ',format(mean(cv, na.rm=T),digits=2)),
                                                by=name])
ggplot(labXRF_format[!(labXRF_format$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density() +
  facet_wrap(~name, scales='free') + 
  labs(x='Coefficient of variation', y='Count') + 
  geom_text(data=cvmean_labeltree, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  theme_classic() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14))

for (elem in unique(labXRF_format[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
  print(elem)
  png(file.path(inspectdir, paste0('labXRF_withintree_',elem,'.png')), width = 20, height=12, units='in', res=300)
  print(
    ggplot(labdt,
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


# ---- 2. Assess within-site variability ---- 
lab_artaxmeansite <- labdt[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID'), .SDcols=29:51]
artaxsdsite <- labdt[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), .SDcols=29:51]
labXRF_formatsite <- merge(melt(lab_artaxmeansite, id.vars='SiteID', variable.name='Elem', value.name='mean'),
                           melt(artaxsdsite, id.vars='SiteID', variable.name='Elem', value.name='sd'),
                           by=c('SiteID','Elem')) %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')
cvmean_labelsite <- as.data.frame(labXRF_formatsite[!(labXRF_formatsite$Elem %in% c('Rh','Pd','Ar')),
                                                    paste0('Mean CV: ',format(mean(cv, na.rm=T),digits=2)),
                                                    by=name])

ggplot(labXRF_formatsite[!(labXRF_formatsite$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density()+
  geom_text(data=cvmean_labelsite, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name, scales='free') + 
  theme_bw() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14),
        panel.border = element_blank(),
        axis.line = element_line(color='black'))

for (elem in unique(labXRF_format[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
  print(elem)
  png(file.path(inspectdir, paste0('labXRF_withinsite_',elem,'.png')), width = 20, height=12, units='in', res=300)
  print(
    ggplot(labXRF_format[Elem == elem,], 
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
# ---- 3. Univariate flag based on within-tree CV for metals with CV < 0.5----
labXRF_format[, `:=`(meanCV=mean(cv, na.rm=T),
                     sdCV = sd(cv, na.rm=T)), by=Elem]

#NAcount <- function(x) length(which(!is.na(unlist(x))))
labXRF_format[(cv>(meanCV+2*sdCV)) & !(Elem %in% excol) & !is.na(cv),
              CVflag_count_lab := .N, by=.(SiteID, Pair)][
                , CVflag_count_lab := ifelse(is.na(CVflag_count_lab) | is.na(cv), 
                                             as.integer(round(mean(CVflag_count_lab, na.rm=T))),
                                             CVflag_count_lab),
                by=.(SiteID, Pair)][
                  is.na(CVflag_count_lab), CVflag_count_lab := 0]


# ---- 4. Univariate relationship to ICP-OES for the purpose of outlier detection ----
#Outlier diagnostic plots
ICPlabmerge[, ICPlab_flags := 0]
for (chem in unique(ICPlabmerge$Elem)) {
  print(chem)
  ICPlab_lm <- lm(ICP ~ mean, data = ICPlabmerge[Elem == chem,])
  ggsave(file.path(inspectdir, paste0('labXRFICP_regoutliers', chem, '.png')),
         chemregdiagnostic_custom(ICPlab_lm, ICPlabmerge, chem,  flagcol = 'ICPlab_flags', tresids = TRUE),
         width = 20, height=12, units='in', dpi=300)
  ICPlabmerge[Elem == chem, ICPlabR2 := summary(ICPlab_lm)$adj.r.squared]
}

#Plot lab XRF ~ ICP data
ICPlab_plot <- ggplot(ICPlabmerge[!(is.na(ICP) |  
                                      Elem %in% c(excol, excol2)),], 
                      aes(x=mean, y=ICP, color=factor(ICPlab_flags), group=1)) + 
  geom_linerangeh(aes(xmin=mean-sd, xmax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
ICPlab_plot

#Investigate determinants of R2
ICPmeanR2 <- unique(ICPlabmerge[, list(meanXRF = mean(mean, na.rm=T),
                                       meanICP = mean(ICP, na.rm=T),
                                       R2 = ICPlabR2), by=Elem][
                                         keVmaxnet, on='Elem==Element'])
ggplot(ICPmeanR2, aes(x=meanXRF, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=meanICP, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(ICPmeanR2, aes(x=noise, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

summary(loess(R2~meanXRF, data=ICPmeanR2, span=1))
summary(loess(R2~meanXRF+noise, data=ICPmeanR2, span=1))

# ---- 5. Univariate relationship to pollution predictors for the purpose of outlier detection ----
#Outlier diagnostic plots
pollutlabmerge[, pollutlab_flags := 0]
for (chem in unique(pollutlabmerge$Elem)) {
  print(chem)
  pollutlab_lm <- lm(mean ~  heat_binglog300*heatOSMAADTlog300*heatOSMgradientlog200*nlcd_imp_ps +
                       heatbustransitlog300 + heatbustransitlog300:heatOSMgradientlog200, 
                     data = pollutlabmerge[Elem == chem,])
  ggsave(file.path(inspectdir, paste0('labpollution_regoutliers', chem, '.png')),
         chemregdiagnostic_custom(pollutlab_lm, pollutlabmerge, chem,  flagcol = 'pollutlab_flags', tresids = TRUE),
         width = 20, height=12, units='in', dpi=300)
  pollutlabmerge[Elem == chem, `:=`(pollutpred = fitted(pollutlab_lm),
                                    pollutlabR2 = summary(pollutlab_lm)$adj.r.squared)]
}

#Plot predicted lab XRF ~ observed lab XRF
pollutlab_plot <- ggplot(pollutlabmerge[!(Elem %in% c(excol, excol2, NA)),], 
                         aes(x=pollutpred, y=mean, group=1, color=factor(pollutlab_flags))) + 
  geom_linerange(aes(ymin=mean-sd, ymax=mean+sd), alpha=1/2) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  geom_smooth(method='lm') +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()

ggsave(file.path(inspectdir, 'labpollution_all.png'),
       pollutlab_plot,
       width = 20, height=12, units='in', dpi=600)

#Investigate determinants of R2
pollutmeanR2 <- unique(pollutlabmerge[, list(meanXRF = mean(mean, na.rm=T),
                                             R2 = pollutlabR2), by=Elem][
                                               keVmaxnet, on='Elem==Element'])
ggplot(pollutmeanR2, aes(x=meanXRF, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=Energy.keV, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=noise, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()



### ---- 6. Compile outlier flags ----
#Merge all flag datasets
labXRF_formatflags <- merge(labXRF_format, ICPlabmerge[, .(SiteID, Pair, Elem, ICPlab_flags, ICPlabR2)], 
                            by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(labfieldmerge[, .(SiteID, Pair, Elem, labfield_flags, labfieldR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(pollutlabmerge[, .(SiteID, Pair, Elem, pollutlab_flags, pollutlabR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T)

ggplot(labXRF_formatflags, aes(x=ICPlabR2, y=pollutlabR2, label=Elem)) + 
  geom_text()
ggplot(labXRF_formatflags, aes(x=labfieldR2, y=pollutlabR2, label=Elem)) + 
  geom_text()

#Compute total number of flags for each tree across all elements that have an R2 of at least 0.10 across models
#Don't want to count outliers from spurious models
flagelems <- labXRF_formatflags[, min(ICPlabR2, labfieldR2, pollutlabR2, na.rm=T)>0.20, by=Elem][
  V1==TRUE & Elem != 'Rh', Elem]
labXRF_formatflags[Elem %in% flagelems,
                   flagpartsum_lab := sum(ICPlab_flags, labfield_flags, pollutlab_flags, na.rm=T), 
                   by=.(SiteID, Pair)][
                     , flagsum_lab := sum(flagpartsum_lab, CVflag_count_lab, na.rm=T), by=.(SiteID, Pair)]

#Check total number of flags by each type of flag
labXRF_formatflags[Elem %in% flagelems,
                   `:=`(ICPlab_flags_sum = sum(ICPlab_flags, na.rm=T),
                        labfield_flags_sum = sum(labfield_flags, na.rm=T),
                        pollutlab_flags_sum = sum(pollutlab_flags, na.rm=T)),
                   by=.(SiteID, Pair)]

labXRF_formatflags_u <- unique(
  labXRF_formatflags[!is.na(flagpartsum_lab),
                     .SD, 
                     .SDcols=c('SiteID', 'Pair', 'CVflag_count_lab',
                               grep('flag.*sum', colnames(labXRF_formatflags), value=T))])

ggplot(labXRF_formatflags_u, aes(x=ICPlab_flags_sum, y=labfield_flags_sum, 
                                 label=paste0(SiteID, Pair), color = CVflag_count_lab)) + 
  scale_color_distiller(palette='YlGnBu') +
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5) 

ggplot(labXRF_formatflags_u, aes(x=ICPlab_flags_sum, y=pollutlab_flags_sum, 
                                 label=paste0(SiteID, Pair), color = CVflag_count_lab)) + 
  scale_color_distiller(palette='YlGnBu') + 
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5)

#Records to inspect
flagelems
labXRF_inspect <- unique(labXRF_formatflags[(ICPlab_flags_sum>5) | 
                                              (labfield_flags_sum > 5) | 
                                              (pollutlab_flags_sum > 5), 
                                            .(SiteID, Pair)])
outliertrees <- trees[labXRF_inspect, on=c('SiteID', 'Pair')]
outlierlocs <- SpatialPointsDataFrame(coords = data.frame(outliertrees$POINT_X, outliertrees$POINT_Y),
                                      data= as.data.frame(outliertrees))
#View(outlierlocs@data)
leaflet(data = outlierlocs) %>% addTiles() %>%
  addMarkers(clusterOptions = markerClusterOptions(),
             popup = ~paste0(SiteID, Pair))

###################### ---- C. ICP-OES ---- ####
# ---- 1. Assess within-site variability ----
ICPmean[, `:=`(SiteID = gsub('[A-Z]', '', SAMPLE.SET),
               Pair = gsub('[0-9]', '', SAMPLE.SET))]
ICPmeansite <- ICPmean[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID'), 
                       .SDcols=which(colnames(ICPmean) %in% periodicTable$symb)]
ICPsdsite <- ICPmean[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), 
                     .SDcols=which(colnames(ICPmean) %in% periodicTable$symb)]
ICP_formatsite <- merge(melt(ICPmeansite, id.vars='SiteID', variable.name='Elem', value.name='mean'),
                        melt(ICPsdsite, id.vars='SiteID', variable.name='Elem', value.name='sd'),
                        by=c('SiteID','Elem')) %>%
  .[, cv := sd/mean] %>%
  merge(periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')
cvmean_labelsite <- as.data.frame(ICP_formatsite[!(ICP_formatsite$Elem %in% c('Rh','Pd','Ar')),
                                                 paste0('Mean CV: ',format(mean(cv, na.rm=T),digits=2)),
                                                 by=name])

ggplot(ICP_formatsite[!(ICP_formatsite$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density()+
  geom_text(data=cvmean_labelsite, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name, scales='free') + 
  theme_bw() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14),
        panel.border = element_blank(),
        axis.line = element_line(color='black'))

for (elem in unique(ICPmelt[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
  print(elem)
  png(file.path(inspectdir, paste0('ICP_withinsite_',elem,'.png')), width = 20, height=12, units='in', res=300)
  print(
    ggplot(ICPmelt[Elem == elem,], 
           aes(x=SiteID, y = ICP, fill=SiteID)) + 
      geom_line(aes(group=SiteID), color='black') +
      geom_point(size=5, colour='black', pch=21, alpha=0.75) +
      labs(x='Element', y='Concentration') + 
      theme_bw() + 
      theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
            strip.text = element_text(size=14),
            panel.border = element_blank(),
            axis.line = element_line(color='black'))
  )
  dev.off()
}

# ---- 2. Univariate relationship to pollution predictors for the purpose of outlier detection ----
#Outlier diagnostic plots
pollutICPmerge[, pollutICP_flags := 0]
for (chem in unique(pollutICPmerge[!is.na(Elem), Elem])) {
  print(chem)
  pollutICP_lm <- lm(ICP ~  heat_binglog300*heatOSMAADTlog300*heatOSMgradientlog200*nlcd_imp_ps +
                       heatbustransitlog300 + heatbustransitlog300:heatOSMgradientlog200, 
                     data = pollutICPmerge[Elem == chem,])
  ggsave(file.path(inspectdir, paste0('ICPpollution_regoutliers', chem, '.png')),
         chemregdiagnostic_custom(pollutICP_lm, pollutICPmerge, chem,  flagcol = 'pollutICP_flags', tresids = TRUE),
         width = 20, height=12, units='in', dpi=300)
  pollutICPmerge[Elem == chem, `:=`(pollutpred = fitted(pollutICP_lm),
                                    pollutICPR2 = summary(pollutICP_lm)$adj.r.squared)]
}

#Plot predicted ICP ~ observed ICP 
pollutICP_plot <- ggplot(pollutICPmerge[!(Elem %in% c(excol, excol2, NA)),], 
                         aes(x=pollutpred, y=ICP, group=1, color=factor(pollutICP_flags))) + 
  geom_text(aes(label = paste0(SiteID,Pair)), size=4, vjust="inward",hjust="inward") +
  geom_smooth(method='lm') +
  #scale_color_distiller(palette= 'Spectral') +
  scale_color_brewer(palette= 'Set1') +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()
pollutICP_plot

ggsave(file.path(inspectdir, 'ICPpollution_all.png'),
       pollutICP_plot,
       width = 20, height=12, units='in', dpi=600)

#Investigate determinants of R2
pollutmeanR2 <- unique(pollutICPmerge[, list(meanICP = mean(ICP, na.rm=T),
                                             R2 = pollutICPR2), by=Elem])[
                                               ICPthresholds_format, on='Elem'
                                               ]
ggplot(pollutmeanR2, aes(x=meanICP, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

ggplot(pollutmeanR2, aes(x=meanICP/QUANTLIM, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

### ---- 3. Compile outlier flags ----
#Merge all flag datasets
ICP_formatflags <- merge(fieldXRF_format, ICPfieldmerge[, .(SiteID, Pair, Elem, ICPfield_flags, ICPfieldR2)], 
                         by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(ICPlabmerge[, .(SiteID, Pair, Elem, ICPlab_flags, ICPlabR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T) %>%
  merge(pollutICPmerge[, .(SiteID, Pair, Elem, pollutICP_flags, pollutICPR2)], 
        by=c('SiteID', 'Pair', 'Elem'), all.x=T)

ggplot(ICP_formatflags, aes(x=ICPfieldR2, y=pollutICPR2, label=Elem)) + 
  geom_text()
ggplot(ICP_formatflags, aes(x=ICPlabR2, y=pollutICPR2, label=Elem)) + 
  geom_text()

#Compute total number of flags for each tree across all elements that have an R2 of at least 0.10 across models
#Don't want to count outliers from spurious models
flagelems <- ICP_formatflags[, min(ICPfieldR2, ICPlabR2, pollutICPR2, na.rm=T)>0.20, by=Elem][
  V1==TRUE & Elem != 'Rh', Elem]
ICP_formatflags[Elem %in% flagelems,
                flagpartsum_ICP := sum(ICPfield_flags, ICPlab_flags, pollutICP_flags, na.rm=T), 
                by=.(SiteID, Pair)][
                  , flagsum_ICP := sum(flagpartsum_ICP, CVflag_count_lab, na.rm=T)]

#Check total number of flags by each type of flag
ICP_formatflags[Elem %in% flagelems,
                `:=`(ICPfield_flags_sum = sum(ICPfield_flags, na.rm=T),
                     ICPlab_flags_sum = sum(ICPlab_flags, na.rm=T),
                     pollutICP_flags_sum = sum(pollutICP_flags, na.rm=T)),
                by=.(SiteID, Pair)]

ICP_formatflags_u <- unique(
  ICP_formatflags[!is.na(flagpartsum_ICP),
                  .SD, 
                  .SDcols=c('SiteID', 'Pair',
                            grep('flag.*sum', colnames(ICP_formatflags), value=T))])

ggplot(ICP_formatflags_u, aes(x=ICPfield_flags_sum, y=ICPlab_flags_sum, 
                              label=paste0(SiteID, Pair))) + 
  scale_color_distiller(palette='YlGnBu') +
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5) 

ggplot(ICP_formatflags_u, aes(x=ICPfield_flags_sum, y=pollutICP_flags_sum, 
                              label=paste0(SiteID, Pair))) + 
  scale_color_distiller(palette='YlGnBu') + 
  geom_text(position = position_jitter(width=0.25, height=0.25), size=5)

#Records to inspect
flagelems
ICP_inspect <- unique(ICP_formatflags[(ICPfield_flags_sum>5) | 
                                        (ICPlab_flags_sum > 5) | 
                                        (pollutICP_flags_sum > 5), 
                                      .(SiteID, Pair)])
outliertrees <- trees[ICP_inspect, on=c('SiteID', 'Pair')]
outlierlocs <- SpatialPointsDataFrame(coords = data.frame(outliertrees$POINT_X, outliertrees$POINT_Y),
                                      data= as.data.frame(outliertrees))
#View(outlierlocs@data)
leaflet(data = outlierlocs) %>% addTiles() %>%
  addMarkers(clusterOptions = markerClusterOptions(),
             popup = ~paste0(SiteID, Pair))

###################### ---- D. Compile all flags and make a flag matrix/heatmap then inspect data and decide on their fate ---- ####
# ---- 1. Compile flags and make heatmap ----
allflags <- ICP_formatflags_u[labXRF_formatflags_u, on=.(SiteID, Pair)][
  fieldXRF_formatflags_u, on=.(SiteID, Pair)] 
allflags[, (grep('i[.]', colnames(allflags), value=T)) := NULL][
  , flagpartsum := sum(flagpartsum_ICP, flagpartsum_lab, flagpartsum_field, na.rm=T), by=.(SiteID, Pair)][
    , SiteIDPair := factor(paste0(SiteID, Pair), levels = unique(paste0(SiteID, Pair)[order(-flagpartsum)]))]

colnames(allflags)
#levels(allflags_melt$variable)
allflags_melt <- melt(allflags, id.vars = c('SiteIDPair', 'SiteID', 'Pair')) %>%
  .[, variable := factor(gsub('[_]|flag*|sum', '', variable), 
                         levels=gsub('[_]|flag|sum', '',
                                     c("flagsum_lab", "flagsum_field", "flagpartsum",
                                       "flagpartsum_ICP" ,"flagpartsum_lab", "flagpartsum_field",  
                                       "pollutICP_flags_sum", "ICPlab_flags_sum","pollutlab_flags_sum",
                                       "CVflag_count_lab", "labfield_flags_sum", "ICPfield_flags_sum",
                                       "pollutfield_flags_sum", "CVflag_count_field")))] 

ggplot(data = allflags_melt[!(variable %in% c('field', 'lab', 'part'))],  #Re-add part if needed
       aes(x=SiteIDPair, y=variable, fill=value)) + 
  geom_tile() +
  scale_fill_distiller(palette='Spectral', trans = "sqrt") + 
  theme(legend.position = c(0.9, 0.2),
        legend.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=45))

# ---- 2. Analysis of the main outliers ----
#49A: on slope by intersection SW of Garfield School in Central District
#     outlier for all lab things but not for field stuff. Must indicate that something was wrong with lab.
#     no collection or analysis comments, moss appearance NTR, overestimated for both ICP and lab XRF.
#     results completely different from 49B which is not an outlier, suggesting erroneous aspect.
#     Must be either the subsample of moss taken or the time between collecting and processing?
#     -> Remove from ICP and lab results

#23A: Harvard I5 exit towards 520 east
#     not really a problem for either field, lab, or ICP related to pollution drivers.
#     no collection or analysis comments, moss appearance NTR
#     ICP-pollution: on the higher pollution level so leverage but not truly outside of the cloud. Surely mostly due to slight heteroscedaticity
#     field-lab outlier for Br, Ca, Cu, Fe, Mn, Ti, Zr (higher for lab)
#     field-ICP outlier for Cu, Fe, Mn, Ni, Pb, Ti 
#     lab-ICP - not an outlier
#     normal within tree variability for lab and field
#     -> NTR keep it 

#51A: Downhill side of street on 15th Ave NE UW
#     multiple outlier flags for pollut-ICP and pollut-lab but nothing with fields. 
#     Must indicate that something was wrong with lab.
#     no collection or analysis comments, moss appearance NTR
#     pollut-field outlier for Ca, K, and Zn (but only because high value and variability at these levels)
#     pollut-lab outlier for Fe, Mn, Ti, Zn 
#     pollut-ICP outlier for Fe, Ni, Zn, but nothing striking
#     lab-ICP - not an outlier
#     field-lab and field-ICP - not an outlier
#     normal within tree variability for lab and field 
#     -> NTR keep it

#54A: Monroe Highway 2 and Main St intersection
#     not outlier for ICPlab and pollut-field
#     low outlier: pollut-ICP and pollut-lab, higher outlier: field-lab and field-ICP
#     suggests something was wrong with lab. moss was very wet when collected, could have rotten.
#     no collection or analysis comments, moss appearance NTR
#     pollut-lab outlier: Fe, Ti (quite far)
#     pollut-ICP outlier: Fe, Ti(but not completely outside)
#     field-lab and -ICP outlier: Fe, Mn, Zn
#     field within-tree variability: one high Fe and Zr value; otherwise nothing much
#     -> NTR keep it

#53A: Monroe - Lewis St and Highway 2 intersection
#     Fine with ICP, outlier with pollut-lab and pollut-field
#     no collection or analysis comments, moss appearance NTR
#     field-pollution outlier: Ca, K, Ni, Pb, 
#     lab-pollution outlier: Fe, Pb, Ti
#     ICP-pollution: None 
#     high within-tree and  for field and lab (for K, Pb, Ti, Si) 
#     high within-site variability: Ca, Co, Cr, Cu, K, Mn, Ni, Ti, Zn
#     -> could just keep 53B instead?

#19A: along I5 in Eastlake
#     high within-tree variability for field (Cu, Fe, Mg, Ti, Zn, Zr)
#     no collection or analysis comments, moss appearance NTR
#     lab XRF data are a bit off (must mean that a clump was taken that does not reflec the full tree)
#     It seems that one of the three XRF samples overshoots for multiple of the flagelems
#     ICP-lab outlier for Cu, Fe, Sr, Zn (not by huge amount)
#     -> keep it for now?

#62A: Sultan under bridge by boat launch over Skykomish
#     Site with a lot of Al, Cr, Fe, Nickel and Zirkonium - must be because was under bridge
#     no collection or analysis comments, moss appearance NTR
#     ICP-lab and ICP-fields quite off: Cr seems just that there might be a non-linear pattern, Fe, Ti, and Zn ICP also higher, 
#     pollutfields and labfields a little off but not much
#     It might be because it was under the bridge and receiving more heavy dust?
#     -> NTR keep it

#53B: Monroe - Lewis St and Highway 2 intersection
#     A little off for everything aside from lab-fields
#     no collection or analysis comments, moss appearance NTR
#     Not really outlier for pollut-fields, a little off for Zn and Fe pollut-ICP, a bit of an outlier for Zn pollut-lab
#     not much else
#     -> NTR keep it

#7A: Stone Ave North of Pacific Ave intersection downhill side
#     field measurement is off pollut-field (most elements)
#     lab-field and  ICP-field (Ca, Cu, Mn mostly) are off
#     no collection or analysis comments, moss appearance NTR
#     very different from 7B from Ba, Cr, Cu, Mn, Mo, Pb, Ti, Zr
#     normal within tree variability
#     -> remove from field analysis

#20B: horizontal tree in Eastlake below I5
#     SMall outlier for ICPfields, labfields, ICPlab, pollutlabs - seems an XRF issue
#     no collection or analysis comments, moss appearance NTR
#     High amounts of Pb, Ti, and Zr.. must be like 19A, low to the ground new I5
#     outlier for Pb for pollut-ICP, pollut-lab, field-lab (also for Sr), field-ICP, etc. just whacky
#     lab-ICP perfect for Pb, off for Cu
#     -> Keep it, high Pb (and some others) must lead to higher variability

#33B: Cap Hill by volunteer park
#     little outlier for pollut-lab and pollut-ICP, larger outlier for pollut-field
#     no collection or analysis comments, moss appearance NTR
#     large outlier across all pollution predictions for Cu
#     small pollut-field outlier for Zn
#     -> very similar value to 33A for Cu, must mean that it's pollution from other source? Keep 

#16A: By viewpoint in blackberry bush in Discovery Park
#     no collection or analysis comments, moss appearance NTR
#     off for pollut-field: higher amounts of Cr, Fe, Ni, Ti, and Zr than predicted (pretty within-tree variability in those elements as well). 
#     underpredicted amounts in lab and ICP for a few elements as well. Contamination could be sea-borne?
#     -> NTR: keep

#20A: same as 20B. 
#     no collection or analysis comments, moss appearance NTR
#     -> NTR

#44B: by Denny Way near Hotel 
#     Very dirty-sooty moss though no collection or analysis comments, moss appearance NTR
#     pollutICP (Fe, Mn), ICPfields and labfields (Fe, Zn - very unpredicted by labs/overpredicted by field)
#     -> NTR

#52B: Monroe intersection of 522 and highway 2
#     Very wet moss
#     a little within-tree variability for both XRF 
#     ICP detected high levels of Se, Cd, As, and Ni hence outlier for lab-ICP and field-ICP
#     not really outlier for pollut-field
#     -> keep

#52A: Monroe intersection of 522 and highway 2
#     very wet moss
#     medium outlier for pollut-field: Fe, Zn, Zr (under-predicted compared to observed/higher iron than predicted)
#     large difference between 52A and 52B, 
#     within-tree over- vs. underestimate vary by element in Fe (under), Zn (over), Zr (over)
#     -> Remove, keep only 52B for pollut-field

#43A: downtown by Amazon building
#     high outlier for pollut-field: Ca, Cu, K, Mn (over-predicted compared to observed)
#     no collection or analysis comments, moss appearance NTR
#     -> NTR

#1B:  calm neighborhood by lake city - first site
#     wrong moss species
#     medium outlier for pollut-field: both 1A and 1B are outliers for Ca, K, and Sr - must be because different species
#     but these elements don't have strong relationship anyways
#     -> NTR

#3B: Lake City Way
#     no comments
#     medium outlier: underpredicted Fe and Zr (but matches ICP and lab well for these elements)

#25A:
#     high within-tree CV for Cr - nothing to report

#34A: Capitol Hill 
#     extreme variability for Pb and outlier for pollutfield, not sure why.
#     remember there might have been some construction there?

#44A: near Denny way downtown
#     pollutICPS off - see 44B- underpredicted levels of Cr, Fe, Ti, Cu, and Zn


#Inspect within tree variability for field XRF data:
#19A: it seems that one of the three XRF samples overshoots for multiple of the flagelems
#20A: NTR
#23A: some extra variability for Zn, Fe, and Cu
#35A: NTR
#44B; NTR
#49A: NTR
#51A: NTR
#53A: one underestimates Zn, NTR otherwise
#54A: Fe one much higher, NTR otherwise
#62A: One a bit higher Fe, NTR otherwise
#7A: one lower Ca, Cu, one higher Fe 
#34A has one point with extra variability for Pb

# ---- 3. Removal of most egregious outliers and those with a paired tree ----
#Remove outliers for model exploration
pollutfieldclean <- pollutfieldmerge[!(paste0(SiteID, Pair) %in% c('7A','52A')),]
pollutlabclean <- pollutlabmerge[!(paste0(SiteID, Pair) %in% c('53A', '49A')),]
pollutICPclean <-  pollutICPmerge[!(paste0(SiteID, Pair) %in% c('49A')),]

############################################################################################################################################

# 4. Check relationship among and reduce variables 
############################################################################################################################################
#---- A. Field XRF - all elements against each other ----
castformula <- as.formula(paste0(
  paste0(colnames(pollutfieldclean)[!(colnames(pollutfieldclean) %in%
                                        c('Elem', 'name', 'mean','range', 'sd', 'cv', 'transmean',
                                          'pollutfield_flags', 'pollutfieldR2', 'pollutpred'))],
         collapse='+'),
  '~Elem'))
pollutfieldclean_cast <- dcast(pollutfieldclean[!is.na(mean)], 
                               formula = castformula,
                               value.var= 'mean')

#Create vector of columns - 1.G. and  3.A.3. for column selection
"***NOTE: consider removing all elems with net/background photon count < 0.1"
elemcols <- colnames(pollutfieldclean_cast)[colnames(pollutfieldclean_cast) %in% periodicTable$symb] %>%
  setdiff('Rh')
elemcols_sub <- elemcols %>%
  setdiff(c(excol, excol2))

#Scatterplot matrix 
png(file.path(inspectdir, 'FieldElem_FieldElem_matrix.png'), width = 24, height=24, units='in', res=300)
ggscatmat(as.data.frame(pollutfieldclean_cast[, elemcols, with=F]),
          alpha=0.7)
dev.off()

#Correlation heatmap
png(file.path(inspectdir, 'corheatmap_FieldElem_FieldElem.png'), width = 12, height=12, units='in', res=300)
corr_heatmap(pollutfieldclean_cast[, elemcols, with=F])
dev.off()

#PCA with rrcov (see 4.2 p24 of Todorov and Filzmoser 2009)
"Interestingly, PcaGrid, a robust PCA, gives completely different results than classic PCA with PC1 and PC2 
being clearly correlated"
pca <- PcaClassic(~., data=pollutfieldclean_cast[,elemcols_sub,with=F], scale=TRUE) #z-standardize
summary(pca)
screeplot(pca, main="Screeplot: classic PCA", bstick=TRUE) #First PC significant compared to broken stick
ordi.monte(pollutfieldclean_cast[,elemcols_sub,with=F],ord='pca',dim=5) #2PCs significant with Monte-carlo test of component significance
#Plot
#Plot
#Plot
#Plot
#Plot
biplot(pca, main="Robust biplot", col=c("gray55", "red"))  
plot(pca)

pollutfieldclean_pca<- cbind(pollutfieldclean_cast, pca$scores)
setnames(pollutfieldclean_pca, colnames(pca$scores), paste0(colnames(pca$scores), '_scores'))
setnames(pca_load, colnames(pca_load), paste0(colnames(pca_load), '_loadings'))

#Create a PCA biplot matrix where all components are graphed against each other
pcabiplot_grid(pollutfieldclean_cast, nPCs = 5, cols = elemcols_sub,
               idcols = c('SiteID', 'Pair'),  scload = 3)

#Inspect components to decide which ones to predict
loadings(pca)
#There are no very distinct patterns - everything is pretty correlated:
#First component: Zn, Zr, Ti, Ni, Fe, Cu, and Cr load equally (0.33-0.36);
#                 Ca, K, and Pb a little less (0.20-25);
#                 Sr, and Mn not much (0.11-0.12)
#Second component:Ca, K, and Sr load the strongest (0.47-0.49)
#                 Pb and Mn load positively but not much (0.16)
#                 Cu and Ni are ~ 0
#                 Zn, Zr, Ti, Fe and Cr load negatively -0.11 to -0.30

#Select individual elements to predict separately: Zn, Fe, Cu, Pb

#---- B. Create a synthetic sum-based pollution index (averaging centered and standardized elements)   ----
#Check that all can be transformed using the same transformation then transform
selcols <- c('Cu', 'Fe', 'Pb', 'Zn')
selcols_trans <- data.trans(as.data.frame(pollutfieldclean_cast[, elemcols[elemcols %in% selcols], with=F]), 
                            method='power',exp=1/3, plot=F)
cbcols <- paste0(elemcols[elemcols %in% selcols], 'cubrt')
pollutfieldclean_cast[, (cbcols) := selcols_trans]
#standardize by max value (use of transformed index actually leads to heteroscedacity)
pollutfieldclean_cast[, cbpollution_index := 100*Reduce('+', lapply(.SD, function(x) (x/max(x))))/length(cbcols), 
                      .SDcols = cbcols][
                        ,pollution_index := 100*Reduce('+', lapply(.SD, function(x) (x/max(x))))/length(cbcols),
                        .SDcols = elemcols[elemcols %in% selcols]]

#---- C. All pollution drivers against each other ----
#Define columns to analyze
pollutcols <- grep('heat|NLCD', colnames(pollutfieldclean), value=T, ignore.case = T)

#Correlation heatmap
png(file.path(inspectdir, 'corheatmap_PollutionDrivers_PollutionDrivers.png'),
    width = 20, height=20, units='in', res=300)
corr_heatmap(pollutfieldclean_cast[, pollutcols, with=F]) + 
  scale_x_discrete(labels=heatlab) + 
  scale_y_discrete(labels=heatlab)
dev.off()

"Overall, all predictors are correlated with increasing correlation with increasing kernel size
as it picks up the common effect of increased road density on the indices. 
Within categories, indices are correlated > 0.90 with indices of 1 for pow1,2 and log.
Across categories (for kernel > 50 m):
gradient-AADT: 0.3-0.6; gradient_bing: 0.2-0.3; gradient-LU: 0.15-0.20; gradient-transit: 0.25-0.4; gradient-SPD:0.4-0.5
AADT-bing: 0.45-0.75 depending on kernel correspondence; AADT-LU: 0.25-0.40; AADT-transit: 0.65-0.80; AADT-SPD:0.6-0.7 particularly for OSMAADT
bing-LU: 0.5, bing-transit: 0.45-0.60; bing-SPD: 0.5-0.60
transit-LU: 0.35-0.45"

#---- D. All pollution drivers against field XRF all elems ----
png(file.path(inspectdir, 'corheatmap_PollutionDrivers_FieldElem.png'),
    width = 22, height=20, units='in', res=300)
corr_heatmap(xmat=pollutfieldclean_cast[, pollutcols, with=F],
             ymat=pollutfieldclean_cast[, c(elemcols, 'pollution_index', 'cbpollution_index'), with=F],
             clus = FALSE) +
  scale_y_discrete(labels=heatlab) + 
  theme(text=element_text(size=24))
dev.off()

"Best univariate predictors for each category:
- heat_binglog300
- nlcd_imp_ps (not filtered)
- heatbustransitlog300 (comparable to two others though much better than 100 for zinc)
- heatOSMAADTlog300
- heatOSMgradientlog200 (should produce 300 version)
- heatOSMSPDlog100 (should probably get 200 when available)"

#---- E. Selected pollution drivers against selected field XRF all elems ----
selcols <- c(selcols, 'heat_binglog300', 'nlcd_imp_ps', 'heatbustransitlog300',
             'heatOSMAADTlog300', 'heatOSMgradientlog200', 'heatOSMSPDlog100')

png(file.path(inspectdir, 'cormatrix_SelFieldElem_SelPollutionDrivers.png'), 
    width = 12, height=12, units='in', res=300)
ggscatmat(as.data.frame(pollutfieldclean_cast[, selcols, with=F]),
          alpha=0.7)
dev.off()

#---- F. Create a synthetic PCA-based pollution driver index  ----
pollutpca <- PcaClassic(~., data=pollutfieldclean_cast[,selcols[selcols %in% pollutcols],with=F], scale=TRUE) #z-standardize
summary(pollutpca)
screeplot(pollutpca, main="Screeplot: classic PCA", bstick=TRUE) #First PC significant compared to broken stick
ordi.monte(pollutfieldclean_cast[,selcols[selcols %in% pollutcols],with=F],ord='pca',dim=5) #2PCs significant with Monte-carlo test of component significance
pcabiplot_grid(pollutfieldclean_cast, nPCs = 3, cols = selcols[selcols %in% pollutcols],
               idcols = c('SiteID', 'Pair'),  scload = 3)
loadings(pollutpca)

############################################################################################################################################

# 5. Model selection
############################################################################################################################################
#--------------- A. Synthetic index ~ separate predictors ----
modlistA <- list() #List to hold models
modlistA[[1]] <- lm(pollution_index ~ 1, data = pollutfieldclean_cast) #Null/Intercept model

#------ 1. Single parameter models --------
modlistA[[2]] <- lm(pollution_index ~ heatOSMAADTlog300, data = pollutfieldclean_cast)
ols_regress(modlistA[[2]])
#ols_plot_diagnostics(modlistA[[2]])

modlistA[[3]] <- lm(pollution_index ~ heat_binglog300, data = pollutfieldclean_cast)
ols_regress(modlistA[[3]])
#ols_plot_diagnostics(modlistA[[3]])

modlistA[[4]] <- lm(pollution_index ~ nlcd_imp_ps, data = pollutfieldclean_cast)
ols_regress(modlistA[[4]])
#ols_plot_diagnostics(modlistA[[4]])

modlistA[[5]] <- lm(pollution_index ~ heatbustransitlog300, data = pollutfieldclean_cast)
ols_regress(modlistA[[5]])
#ols_plot_diagnostics(modlistA[[5]])

modlistA[[6]] <- lm(pollution_index ~ heatOSMSPDlog100, data = pollutfieldclean_cast)
ols_regress(modlistA[[6]])
#ols_plot_diagnostics(modlistA[[6]])

modlistA[[7]] <- lm(pollution_index ~ heatOSMSPDlog100 + I(heatOSMSPDlog100^2), data = pollutfieldclean_cast)
predf <- cbind(pollutfieldclean_cast, predict(modlistA[[7]], interval = 'confidence'))
ggplot(pollutfieldclean_cast, aes(x=heatOSMSPDlog100)) + 
  geom_ribbon(data = predf, aes(ymin=lwr, ymax=upr), fill='orange') +
  geom_point(aes(y=pollution_index)) +
  geom_line(data = predf,aes(y=fit)) +
  geom_smooth(aes(y=pollution_index), span=1) +
  theme_classic()
ols_regress(modlistA[[7]])
#ols_plot_diagnostics(modlistA[[7]])

modlistA[[8]] <- lm(pollution_index ~ heatOSMgradientlog200, data = pollutfieldclean_cast)
ols_regress(modlistA[[8]])
#ols_plot_diagnostics(modlistA[[8]])

#------ 2. Multiparameter models --------
modlistA[[9]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300, data = pollutfieldclean_cast)
ols_regress(modlistA[[9]])
#ols_plot_diagnostics(modlistA[[9]])
ols_coll_diag(modlistA[[9]])
ols_correlations(modlistA[[9]])

modlistA[[10]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300 + heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[10]])
#ols_plot_diagnostics(modlistA[[10]])
ols_coll_diag(modlistA[[10]])
ols_correlations(modlistA[[10]])

modlistA[[11]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
# subdat <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% c('23A', '34A', '33B')),]
# modlistA[[11]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100, 
#                      data = subdat)
ols_regress(modlistA[[11]])
#ols_plot_diagnostics(modlistA[[11]])
ols_coll_diag(modlistA[[11]])
ols_correlations(modlistA[[11]])
# theme_set(theme_sjplot())
# plot_model(modlistA[[11]], type='pred', terms = c('heat_binglog300', 'heatOSMSPDlog100'))


modlistA[[12]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMAADTlog300, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[12]])
#ols_plot_diagnostics(modlistA[[12]])
ols_coll_diag(modlistA[[12]])
ols_correlations(modlistA[[12]])

modlistA[[13]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300*heatOSMAADTlog300, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[13]])
#ols_plot_diagnostics(modlistA[[13]])
ols_coll_diag(modlistA[[13]])
ols_correlations(modlistA[[13]])
plot_model(modlistA[[13]], type='pred', terms = c('heat_binglog300', 'heatOSMAADTlog100'))


modlistA[[14]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMAADTlog300 + heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[14]])
#ols_plot_diagnostics(modlistA[[14]])
ols_coll_diag(modlistA[[14]])
ols_correlations(modlistA[[14]])

modlistA[[15]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300*heatOSMgradientlog200 + 
                       heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[15]])
#ols_plot_diagnostics(modlistA[[15]])
ols_coll_diag(modlistA[[15]])
ols_correlations(modlistA[[15]])

modlistA[[16]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMSPDlog100*heatOSMgradientlog200, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[16]])
#ols_plot_diagnostics(modlistA[[16]])
ols_coll_diag(modlistA[[16]])
ols_correlations(modlistA[[16]])

modlistA[[17]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMSPDlog100*heatOSMgradientlog200 + nlcd_imp_ps, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[17]])
#ols_plot_diagnostics(modlistA[[17]])
ols_coll_diag(modlistA[[17]])
ols_correlations(modlistA[[17]])

modlistA[[18]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMSPDlog100*heatOSMgradientlog200 + nlcd_imp_ps:heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[18]])
#ols_plot_diagnostics(modlistA[[18]])
ols_coll_diag(modlistA[18])
ols_correlations(modlistA[[18]])

modlistA[[19]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300*nlcd_imp_ps + 
                       heatOSMSPDlog100*heatOSMgradientlog200, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[19]])
#ols_plot_diagnostics(modlistA[[19]])
ols_coll_diag(modlistA[19])
ols_correlations(modlistA[[19]])

modlistA[[20]] <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMSPDlog100 + nlcd_imp_ps, 
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[20]])
#ols_plot_diagnostics(modlistA[[20]])
ols_coll_diag(modlistA[20])
ols_correlations(modlistA[[20]])

#Build models without bus transit
modlistA[[21]] <- lm(pollution_index ~ heatOSMAADTlog300 + heat_binglog300 + heatOSMSPDlog100,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[21]])
#ols_plot_diagnostics(modlistA[[21]])
ols_coll_diag(modlistA[[21]])
ols_correlations(modlistA[[21]])

modlistA[[22]] <- lm(pollution_index ~ heatOSMAADTlog300*heat_binglog300,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[22]])
#ols_plot_diagnostics(modlistA[[22]])
ols_coll_diag(modlistA[[22]])
ols_correlations(modlistA[[22]])

modlistA[[23]] <- lm(pollution_index ~ heatOSMSPDlog100*heat_binglog300,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[23]])
#ols_plot_diagnostics(modlistA[[23]])
ols_coll_diag(modlistA[[23]])
ols_correlations(modlistA[[23]])

modlistA[[24]] <- lm(pollution_index ~ heatOSMAADTlog300*heat_binglog300 + heatOSMSPDlog100,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[24]])
#ols_plot_diagnostics(modlistA[[24]])
ols_coll_diag(modlistA[[24]])
ols_correlations(modlistA[[24]])

modlistA[[25]] <- lm(pollution_index ~ heatOSMAADTlog300*heat_binglog300*heatOSMSPDlog100,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[25]])
#ols_plot_diagnostics(modlistA[[25]])
ols_coll_diag(modlistA[[25]])
ols_correlations(modlistA[[25]])


modlistA[[26]] <- lm(pollution_index ~ heatOSMAADTlog300 + heat_binglog300*heatOSMgradientlog200,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[26]])
#ols_plot_diagnostics(modlistA[[26]])
ols_coll_diag(modlistA[[26]])
ols_correlations(modlistA[[26]])

modlistA[[27]] <- lm(pollution_index ~ heatOSMAADTlog300 + heat_binglog300 + nlcd_imp_ps,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[27]])
#ols_plot_diagnostics(modlistA[[27]])
ols_coll_diag(modlistA[[27]])
ols_correlations(modlistA[[27]])

#Build models without congestion
modlistA[[28]] <- lm(pollution_index ~  heatbustransitlog300 + heatOSMAADTlog300,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[28]])
#ols_plot_diagnostics(modlistA[[28]])
ols_coll_diag(modlistA[[28]])
ols_correlations(modlistA[[28]])

modlistA[[29]] <- lm(pollution_index ~  heatbustransitlog300 + heatOSMSPDlog100,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[29]])
#ols_plot_diagnostics(modlistA[[29]])
ols_coll_diag(modlistA[[29]])
ols_correlations(modlistA[[29]])

modlistA[[30]] <- lm(pollution_index ~  heatbustransitlog300 + heatOSMSPDlog100*heatOSMgradientlog200,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[30]])
#ols_plot_diagnostics(modlistA[[30]])
ols_coll_diag(modlistA[[30]])
ols_correlations(modlistA[[30]])

modlistA[[31]] <- lm(pollution_index ~  heatbustransitlog300 + heatOSMSPDlog100*heatOSMgradientlog200 +
                       nlcd_imp_ps,
                     data = pollutfieldclean_cast)
ols_regress(modlistA[[31]])
#ols_plot_diagnostics(modlistA[[30]])
ols_coll_diag(modlistA[[31]])
ols_correlations(modlistA[[31]])

#------ 3. Make latex model summary table ----
vnum <- max(sapply(modlistA, function(mod) {length(mod$coefficients)}))
model_summary<- as.data.table(
  ldply(modlistA, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary, -R2pred, AICc)  
cat(latex_format(model_summary), file = file.path(moddir, 'pollutionindex_modeltable.tex'))
setwd(moddir)
texi2pdf('pollutionindex_modeltable.tex')

#------ 4. Make latex model summary table when excluding outliers ----
model_summary_nooutliers <- as.data.table(
  ldply(modlistA, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, 
                                                         remove_outliers = 'outliers',
                                                         labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                                         kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summary_nooutliers), file = file.path(moddir, 'pollutionindex_modeltable_nooutliers.tex'))
setwd(moddir)
texi2pdf('pollutionindex_modeltable_nooutliers.tex')

#------ 5. Compare final selected models ----
#Compare model 11 to model 10 for equal removal of outliers
mod10_nooutliers <- regdiagnostic_customtab(modlistA[[10]], maxpar=vnum, 
                                            remove_outliers = 'outliers',
                                            labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                            kCV = TRUE, k=10, cvreps=50)
subdat <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                                    strsplit(gsub('\\\\', '', mod10_nooutliers['outliers']), ',')$outliers),]
mod10_nooutliersub <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300 + heatOSMSPDlog100, 
                          data = subdat)
ols_regress(mod10_nooutlierssub)
ols_plot_diagnostics(mod10_nooutlierssub)
ols_coll_diag(mod10_nooutlierssub)
ols_correlations(mod10_nooutlierssub)
AICc(mod10_nooutlierssub)
qplot(subdat$heatbustransitlog300, mod10_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat$heat_binglog300, mod10_nooutliersub$residuals) +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat$heatOSMSPDlog100, mod10_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')

subdat <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                                  strsplit(gsub('\\\\', '', mod10_nooutliers['outliers']), ',')$outliers),]
mod11_nooutliersub <- lm(pollution_index ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100, 
                          data = subdat)
ols_regress(mod11_nooutlierssub)
ols_plot_diagnostics(mod11_nooutlierssub)
ols_coll_diag(mod11_nooutlierssub)
ols_correlations(mod11_nooutlierssub)
AICc(mod11_nooutlierssub)
qplot(subdat$heatbustransitlog300, mod11_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat$heat_binglog300, mod11_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat$heatOSMSPDlog100, mod11_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')

#For the sake of parsimony, go for model #10

#------ 6. Check spatial and temporal autocorrelation of residuals for full and 'robust' datasets -------
"Fron Anselin 2006: ignoring spatially correlated errors is mostly a problem of efficiency, in the
sense that the OLS coefficient standard error estimates are biased, but the
coefficient estimates themselves remain unbiased. However, to the extent that
the spatially correlated errors mask an omitted variable, the consequences of
ignoring this may be more serious."
#---- For model 10 with all data: pollution_index ~ heatbustransitlog300 + heat_binglog300 + heatOSMSPDlog100 -----
resnorm <- rstandard(modlistA[[10]]) #Get standardized residuals from model
#Make bubble map of residuals
bubbledat <- data.frame(resnorm, pollutfieldclean_cast$coords.x1, pollutfieldclean_cast$coords.x2)
coordinates(bubbledat) <- c("pollutfieldclean_cast.coords.x1","pollutfieldclean_cast.coords.x2")
bubble(bubbledat, "resnorm", col = c("blue","red"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")
#Check semi-variogram of residuals
plot(variogram(resnorm~1, bubbledat, cutoff=2000, width=50)) #isotropic
plot(variogram(resnorm~1, bubbledat, cutoff= 2000, width=50, alpha = c(0, 45, 90,135))) #anisotropic
#Check spline correlogram ()
plot(spline.correlog(x=coordinates(bubbledat)[,1], y=coordinates(bubbledat)[,2],
                     z=bubbledat$resnorm, resamp=500, quiet=TRUE, xmax = 5000))
#Compute a spatial weight matrix based on IDW
weightmat_k <- lapply(1:10, function(i) {
  weightmat_IDW(pollutfieldclean_cast[, .(coords.x1, coords.x2)], knb = i, mindist = 10)}) #Based on 10 nearest neighbors
weightmat_all <- weightmat_IDW(pollutfieldclean_cast[, .(coords.x1, coords.x2)], knb = NULL, mindist = 10) #Based on all points

#Moran plots
#lag_resnorm <- lag.listw(weightmat_all, resnorm) #Can be used to create customized Moran plot by plotting residuals against matrix
moran.plot(resnorm, weightmat_all, labels=pollutfieldclean_cast[,paste0(SiteID, Pair)], pch=19)
moran.plot(resnorm, weightmat_k[[2]], labels=pollutfieldclean_cast[,paste0(SiteID, Pair)], pch=19)

#Compute Moran's I
"Should always only use lm.morantest for residuals from regression, see http://r-sig-geo.2731867.n2.nabble.com/Differences-between-moran-test-and-lm-morantest-td7591336.html
for an explanation"
lm.morantest(modlistA[[10]], listw = listw2U(weightmat_k[[2]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistA[[10]], listw = listw2U(weightmat_all)) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)

#Test for need for spatial regression model using Lagrange Multiplier (LM) tests
lm.LMtests(modlistA[[10]], listw = listw2U(weightmat_k[[2]]), test=c("LMerr","RLMerr", "SARMA"))

#Spatial simultaneous autoregressive error model estimation with 2 nearest neighbors
sarlm_mod10 <- errorsarlm(modlistA[[10]]$call$formula, data = modlistA[[10]]$model, 
                          listw = listw2U(weightmat_k[[2]]))
summary(sarlm_mod10)
#Compare AIC
AIC(sarlm_mod10)
AIC(modlistA[[10]])
#Compare pseudo-R2
cor(modlistA[[10]]$model$pollution_index, fitted(sarlm_mod10))^2
cor(modlistA[[10]]$model$pollution_index, fitted(modlistA[[10]]))^2

#---- For model 10 without outliers: pollution_index ~ heatbustransitlog300 + heat_binglog300 + heatOSMSPDlog100 -----
resnorm_sub <- rstandard(mod10_nooutliersub) #Get standardized residuals from model
resnorm_subdf <- data.frame(resnorm_sub, subdat$coords.x1, subdat$coords.x2)
coordinates(resnorm_subdf) <- c("subdat.coords.x1","subdat.coords.x2")

#Check semi-variogram of residuals
plot(variogram(resnorm_sub~1, resnorm_subdf, cutoff=1000, width=50)) #isotropic
plot(variogram(resnorm_sub~1, resnorm_subdf, cutoff= 1000, width=50, alpha = c(0, 45, 90,135))) #anisotropic
#Check spline correlogram ()
plot(spline.correlog(x=coordinates(resnorm_subdf)[,1], y=coordinates(resnorm_subdf)[,2],
                     z=resnorm_subdf$resnorm_sub, resamp=500, quiet=TRUE, xmax = 5000))
#Compute a spatial weight matrix based on IDW
weightmat_ksub <- lapply(1:10, function(i) {
  weightmat_IDW(subdat[, .(coords.x1, coords.x2)], knb = i, mindist = 10)}) #Based on 10 nearest neighbors
weightmat_allsub <- weightmat_IDW(subdat[, .(coords.x1, coords.x2)], knb = NULL, mindist = 10) #Based on all points

#Moran plots
#lag_resnorm_sub <- lag.listw(weightmat_all, resnorm_sub) #Can be used to create customized Moran plot by plotting residuals against matrix
moran.plot(resnorm_sub, weightmat_allsub, labels=subdat[,paste0(SiteID, Pair)], pch=19)
moran.plot(resnorm_sub, weightmat_ksub[[1]], labels=subdat[,paste0(SiteID, Pair)], pch=19)

#Compute Moran's I
"Should always only use lm.morantest for residuals from regression, see http://r-sig-geo.2731867.n2.nabble.com/Differences-between-moran-test-and-lm-morantest-td7591336.html
for an explanation"
lm.morantest(mod10_nooutliersub, listw = listw2U(weightmat_ksub[[2]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod10_nooutliersub, listw = listw2U(weightmat_allsub)) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)

#Test for need for spatial regression model using Lagrange Multiplier (LM) tests
lm.LMtests(mod10_nooutliersub, listw = listw2U(weightmat_ksub[[2]]), test=c("LMerr","RLMerr", "SARMA"))

#Spatial simultaneous autoregressive error model estimation with 2 nearest neighbors
sarlm_mod10sub <- errorsarlm(mod10_nooutliersub$call$formula, data = mod10_nooutliersub$model, 
                          listw = listw2U(weightmat_ksub[[2]]))
summary(sarlm_mod10sub)
#Compare AIC
AIC(sarlm_mod10sub)
AIC(mod10_nooutliersub)
#Compare pseudo-R2
cor(mod10_nooutliersub$model$pollution_index, fitted(sarlm_mod10sub))^2
cor(mod10_nooutliersub$model$pollution_index, fitted(mod10_nooutliersub))^2
#Compare MAE
DescTools::MAE(mod10_nooutliersub$model$pollution_index, fitted(sarlm_mod10sub))
DescTools::MAE(mod10_nooutliersub$model$pollution_index, fitted(mod10_nooutliersub))

#Compare observed~predicted for full-no outlier model and for aspatial and spatial model
fullsub_comparisonplot <- ggplot(subdat, aes(x=fitted(sarlm_mod10sub), y=pollution_index)) + 
  geom_point(data=modlistA[[10]]$model, aes(x=fitted(modlistA[[10]])), size=2, color='black') +
  geom_point(aes(x=predict(modlistA[[10]], subdat)), size=2, color='grey') +
  geom_point(aes(x=fitted(mod10_nooutliersub)), size=2, alpha=1/2, color='orange') +
  geom_abline(size=1, slope=1, intercept=0, color='red') + 
  coord_fixed() +
  theme_classic()

spatial_comparisonplot <- ggplot(subdat, aes(x=fitted(sarlm_mod10sub), y=pollution_index)) + 
  geom_point(aes(x=fitted(mod10_nooutliersub)), size=2, alpha=1/2, color='orange') +
  geom_point(size=2, alpha=1/2, color='red') + 
  geom_abline(size=1, slope=1, intercept=0, color='red') + 
  coord_fixed() +
  theme_classic()

grid.arrange(fullsub_comparisonplot, spatial_comparisonplot)

#--------------- B. Zn ~ separate predictors ----
#Multiply Zn to make coefficients more readable
pollutfieldclean_cast[, Zn := 1000*Zn]
modlistZn <- list() #List to hold models
modlistZn[[1]] <- lm(Zn ~ 1, data = pollutfieldclean_cast) #Null/Intercept model
#------ 1. Single parameter models --------
modlistZn[[2]] <- lm(Zn ~ heatOSMAADTlog300, data = pollutfieldclean_cast)
ols_regress(modlistZn[[2]])
#ols_plot_diagnostics(modlistZn[[2]])

modlistZn[[3]] <- lm(Zn ~ heat_binglog300, data = pollutfieldclean_cast)
ols_regress(modlistZn[[3]])
#ols_plot_diagnostics(modlistZn[[3]])

modlistZn[[4]] <- lm(Zn ~ nlcd_imp_ps, data = pollutfieldclean_cast)
ols_regress(modlistZn[[4]])
#ols_plot_diagnostics(modlistZn[[4]])

modlistZn[[5]] <- lm(Zn ~ heatbustransitlog300, data = pollutfieldclean_cast)
ols_regress(modlistZn[[5]])
#ols_plot_diagnostics(modlistZn[[5]])

modlistZn[[6]] <- lm(Zn ~ heatOSMSPDlog100, data = pollutfieldclean_cast)
ols_regress(modlistZn[[6]])
#ols_plot_diagnostics(modlistZn[[6]])

modlistZn[[7]] <- lm(Zn ~ heatOSMSPDlog100 + I(heatOSMSPDlog100^2), data = pollutfieldclean_cast)
predf <- cbind(pollutfieldclean_cast, predict(modlistZn[[7]], interval = 'confidence'))
ggplot(pollutfieldclean_cast, aes(x=heatOSMSPDlog100)) + 
  geom_ribbon(data = predf, aes(ymin=lwr, ymax=upr), fill='orange') +
  geom_point(aes(y=Zn)) +
  geom_line(data = predf,aes(y=fit)) +
  geom_smooth(aes(y=Zn), span=1) +
  theme_classic()
ols_regress(modlistZn[[7]])
#ols_plot_diagnostics(modlistZn[[7]])

modlistZn[[8]] <- lm(Zn ~ heatOSMgradientlog200, data = pollutfieldclean_cast)
ols_regress(modlistZn[[8]])
#ols_plot_diagnostics(modlistZn[[8]])

#------ 2. Multiparameter models --------
modlistZn[[9]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300, data = pollutfieldclean_cast)
ols_regress(modlistZn[[9]])
#ols_plot_diagnostics(modlistZn[[9]])
ols_coll_diag(modlistZn[[9]])
ols_correlations(modlistZn[[9]])

modlistZn[[10]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300 + heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[10]])
#ols_plot_diagnostics(modlistZn[[10]])
ols_coll_diag(modlistZn[[10]])
ols_correlations(modlistZn[[10]])

modlistZn[[11]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
# subdat <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% c('23A', '34A', '33B')),]
# modlistZn[[11]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100, 
#                      data = subdat)
ols_regress(modlistZn[[11]])
ols_plot_diagnostics(modlistZn[[11]])
ols_coll_diag(modlistZn[[11]])
ols_correlations(modlistZn[[11]])
# theme_set(theme_sjplot())
# plot_model(modlistZn[[11]], type='pred', terms = c('heat_binglog300', 'heatOSMSPDlog100'))


modlistZn[[12]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMAADTlog300, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[12]])
#ols_plot_diagnostics(modlistZn[[12]])
ols_coll_diag(modlistZn[[12]])
ols_correlations(modlistZn[[12]])

modlistZn[[13]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300*heatOSMAADTlog300, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[13]])
#ols_plot_diagnostics(modlistZn[[13]])
ols_coll_diag(modlistZn[[13]])
ols_correlations(modlistZn[[13]])
plot_model(modlistZn[[13]], type='pred', terms = c('heat_binglog300', 'heatOSMAADTlog100'))


modlistZn[[14]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMAADTlog300 + heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[14]])
#ols_plot_diagnostics(modlistZn[[14]])
ols_coll_diag(modlistZn[[14]])
ols_correlations(modlistZn[[14]])

modlistZn[[15]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300*heatOSMgradientlog200 + 
                       heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[15]])
#ols_plot_diagnostics(modlistZn[[15]])
ols_coll_diag(modlistZn[[15]])
ols_correlations(modlistZn[[15]])

modlistZn[[16]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMSPDlog100*heatOSMgradientlog200, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[16]])
#ols_plot_diagnostics(modlistZn[[16]])
ols_coll_diag(modlistZn[[16]])
ols_correlations(modlistZn[[16]])

modlistZn[[17]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMSPDlog100*heatOSMgradientlog200 + nlcd_imp_ps, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[17]])
#ols_plot_diagnostics(modlistZn[[17]])
ols_coll_diag(modlistZn[[17]])
ols_correlations(modlistZn[[17]])

modlistZn[[18]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMSPDlog100*heatOSMgradientlog200 + nlcd_imp_ps:heatOSMSPDlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[18]])
ols_plot_diagnostics(modlistZn[[18]])
ols_coll_diag(modlistZn[18])
ols_correlations(modlistZn[[18]])

modlistZn[[19]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300*nlcd_imp_ps + 
                       heatOSMSPDlog100*heatOSMgradientlog200, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[19]])
#ols_plot_diagnostics(modlistZn[[19]])
ols_coll_diag(modlistZn[19])
ols_correlations(modlistZn[[19]])

modlistZn[[20]] <- lm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                       heatOSMSPDlog100 + nlcd_imp_ps, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[20]])
#ols_plot_diagnostics(modlistZn[[20]])
ols_coll_diag(modlistZn[20])
ols_correlations(modlistZn[[20]])

#Build models without bus transit
modlistZn[[21]] <- lm(Zn ~ heatOSMAADTlog300 + heat_binglog300 + heatOSMSPDlog100,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[21]])
#ols_plot_diagnostics(modlistZn[[21]])
ols_coll_diag(modlistZn[[21]])
ols_correlations(modlistZn[[21]])

modlistZn[[22]] <- lm(Zn ~ heatOSMAADTlog300*heat_binglog300,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[22]])
#ols_plot_diagnostics(modlistZn[[22]])
ols_coll_diag(modlistZn[[22]])
ols_correlations(modlistZn[[22]])

modlistZn[[23]] <- lm(Zn ~ heatOSMSPDlog100*heat_binglog300,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[23]])
#ols_plot_diagnostics(modlistZn[[23]])
ols_coll_diag(modlistZn[[23]])
ols_correlations(modlistZn[[23]])

modlistZn[[24]] <- lm(Zn ~ heatOSMAADTlog300*heat_binglog300 + heatOSMSPDlog100,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[24]])
#ols_plot_diagnostics(modlistZn[[24]])
ols_coll_diag(modlistZn[[24]])
ols_correlations(modlistZn[[24]])

modlistZn[[25]] <- lm(Zn ~ heatOSMAADTlog300*heat_binglog300*heatOSMSPDlog100,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[25]])
#ols_plot_diagnostics(modlistZn[[25]])
ols_coll_diag(modlistZn[[25]])
ols_correlations(modlistZn[[25]])


modlistZn[[26]] <- lm(Zn ~ heatOSMAADTlog300 + heat_binglog300*heatOSMgradientlog200,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[26]])
#ols_plot_diagnostics(modlistZn[[26]])
ols_coll_diag(modlistZn[[26]])
ols_correlations(modlistZn[[26]])

modlistZn[[27]] <- lm(Zn ~ heatOSMAADTlog300 + heat_binglog300 + nlcd_imp_ps,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[27]])
#ols_plot_diagnostics(modlistZn[[27]])
ols_coll_diag(modlistZn[[27]])
ols_correlations(modlistZn[[27]])

#Build models without congestion
modlistZn[[28]] <- lm(Zn ~  heatbustransitlog300 + heatOSMAADTlog300,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[28]])
#ols_plot_diagnostics(modlistZn[[28]])
ols_coll_diag(modlistZn[[28]])
ols_correlations(modlistZn[[28]])

modlistZn[[29]] <- lm(Zn ~  heatbustransitlog300 + heatOSMSPDlog100,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[29]])
#ols_plot_diagnostics(modlistZn[[29]])
ols_coll_diag(modlistZn[[29]])
ols_correlations(modlistZn[[29]])

modlistZn[[30]] <- lm(Zn ~  heatbustransitlog300 + heatOSMSPDlog100*heatOSMgradientlog200,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[30]])
#ols_plot_diagnostics(modlistZn[[30]])
ols_coll_diag(modlistZn[[30]])
ols_correlations(modlistZn[[30]])

modlistZn[[31]] <- lm(Zn ~  heatbustransitlog300 + heatOSMSPDlog100*heatOSMgradientlog200 +
                       nlcd_imp_ps,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[31]])
#ols_plot_diagnostics(modlistZn[[30]])
ols_coll_diag(modlistZn[[31]])
ols_correlations(modlistZn[[31]])

#------ 3. Make latex model summary table ----
vnum <- max(sapply(modlistZn, function(mod) {length(mod$coefficients)}))
model_summaryZn<- as.data.table(
  ldply(modlistZn, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summaryZn, -R2pred, AICc)  
cat(latex_format(model_summaryZn), file = file.path(moddir, 'modeltable_Zn.tex'))
setwd(moddir)
texi2pdf('modeltable_Zn.tex')
#------ 4. Transform Zn ----
pollutfieldclean_cast[, getlambda(Zn)]
grid.arrange(
  ggplot(pollutfieldclean_cast, aes(x=Zn)) + 
    geom_histogram(fill='blue', alpha=1/2),
  ggplot(pollutfieldclean_cast, aes(x=log(Zn))) + 
    geom_histogram(fill='red', alpha=1/2))

#lambda = 0.2, almost log. 
pollutfieldclean_cast[, logZn := log(Zn)]

modlistlogZn <- list() #List to hold models
modlistlogZn[[1]] <- lm(logZn ~ 1, data = pollutfieldclean_cast) #Null/Intercept model

#------ 5. logZn - Single parameter models --------
modlistlogZn[[2]] <- lm(logZn ~ heatOSMAADTlog300, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[2]])
#ols_plot_diagnostics(modlistlogZn[[2]])

modlistlogZn[[3]] <- lm(logZn ~ heat_binglog300, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[3]])
#ols_plot_diagnostics(modlistlogZn[[3]])

modlistlogZn[[4]] <- lm(logZn ~ nlcd_imp_ps, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[4]])
#ols_plot_diagnostics(modlistlogZn[[4]])

modlistlogZn[[5]] <- lm(logZn ~ heatbustransitlog300, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[5]])
#ols_plot_diagnostics(modlistlogZn[[5]])

modlistlogZn[[6]] <- lm(logZn ~ heatOSMSPDlog100, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[6]])
#ols_plot_diagnostics(modlistlogZn[[6]])

modlistlogZn[[7]] <- lm(logZn ~ heatOSMSPDlog100 + I(heatOSMSPDlog100^2), data = pollutfieldclean_cast)
predf <- cbind(pollutfieldclean_cast, predict(modlistlogZn[[7]], interval = 'confidence'))
ggplot(pollutfieldclean_cast, aes(x=heatOSMSPDlog100)) + 
  geom_ribbon(data = predf, aes(ymin=lwr, ymax=upr), fill='orange') +
  geom_point(aes(y=logZn)) +
  geom_line(data = predf,aes(y=fit)) +
  geom_smooth(aes(y=logZn), span=1) +
  theme_classic()
ols_regress(modlistlogZn[[7]])
#ols_plot_diagnostics(modlistlogZn[[7]])

modlistlogZn[[8]] <- lm(logZn ~ heatOSMgradientlog200, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[8]])
#ols_plot_diagnostics(modlistlogZn[[8]])

#------ 6. logZn - Multiparameter models --------
modlistlogZn[[9]] <- lm(logZn ~ heatbustransitlog300 + heat_binglog300, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[9]])
#ols_plot_diagnostics(modlistlogZn[[9]])
ols_coll_diag(modlistlogZn[[9]])
ols_correlations(modlistlogZn[[9]])
# ols_plot_added_variable(modlistlogZn[[9]])
# ols_plot_comp_plus_resid(modlistlogZn[[9]])

modlistlogZn[[10]] <- lm(logZn ~ heatbustransitlog300 + heat_binglog300 + heatOSMSPDlog100, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[10]])
ols_plot_diagnostics(modlistlogZn[[10]])
ols_coll_diag(modlistlogZn[[10]])
ols_correlations(modlistlogZn[[10]])

modlistlogZn[[11]] <- lm(logZn ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100, 
                      data = pollutfieldclean_cast)
# subdat <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% c('23A', '34A', '33B')),]
# modlistlogZn[[11]] <- lm(logZn ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100, 
#                      data = subdat)
summary(modlistlogZn[[11]])
ols_plot_diagnostics(modlistlogZn[[11]])
ols_coll_diag(modlistlogZn[[11]])
ols_correlations(modlistlogZn[[11]])
theme_set(theme_sjplot())
plot_model(modlistlogZn[[11]], type='pred', terms = c('heat_binglog300', 'heatOSMSPDlog100'))

modlistlogZn[[12]] <- lm(logZn ~ heatbustransitlog300 + heat_binglog300 + heatOSMSPDlog100 + I(heatOSMSPDlog100^2), 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[12]])
ols_plot_diagnostics(modlistlogZn[[12]])
ols_coll_diag(modlistlogZn[[12]])
ols_correlations(modlistlogZn[[12]])


modlistlogZn[[13]] <- lm(logZn ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100 + I(heatOSMSPDlog100^2), 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[13]])
ols_plot_diagnostics(modlistlogZn[[13]])
ols_coll_diag(modlistlogZn[[13]])
ols_correlations(modlistlogZn[[13]])


modlistlogZn[[14]] <- lm(logZn ~ heatbustransitlog300 + heat_binglog300 + 
                        heatOSMAADTlog300, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[14]])
#ols_plot_diagnostics(modlistlogZn[[12]])
ols_coll_diag(modlistlogZn[[14]])
ols_correlations(modlistlogZn[[14]])

modlistlogZn[[15]] <- lm(logZn ~ heatbustransitlog300 + heat_binglog300*heatOSMAADTlog300, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[15]])
#ols_plot_diagnostics(modlistlogZn[[13]])
ols_coll_diag(modlistlogZn[[15]])
ols_correlations(modlistlogZn[[15]])
#plot_model(modlistlogZn[[13]], type='pred', terms = c('heat_binglog300', 'heatOSMAADTlog100'))

modlistlogZn[[16]] <- lm(logZn ~ heatbustransitlog300 + heat_binglog300*heatOSMAADTlog300*heatOSMSPDlog100, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[16]])
ols_plot_diagnostics(modlistlogZn[[16]])
ols_coll_diag(modlistlogZn[[16]])
ols_correlations(modlistlogZn[[16]])

#------ 7. log(Zn) - Make latex model summary table ----
vnum <- max(sapply(modlistlogZn, function(mod) {length(mod$coefficients)}))
model_summarylogZn<- as.data.table(
  ldply(modlistlogZn, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summarylogZn, -R2pred, AICc)  
cat(latex_format(model_summarylogZn), file = file.path(moddir, 'modeltable_logZn.tex'))
setwd(moddir)
texi2pdf('modeltable_logZn.tex')


#Could add Synthetic index ~ PC or Individual pollutants ~ PC but could just add difficulty
#------ 8. glm(Zn, gaussian log link) - Single parameter models  --------
modlistglmZn <- list() #List to hold models
modlistglmZn[[1]] <- glm(Zn ~ 1, data = pollutfieldclean_cast, family=gaussian(link='log')) #Null/Intercept model
par(mfrow=c(2,2))

modlistglmZn[[2]] <- glm(Zn ~ heatOSMAADTlog300, data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[2]])
plot(modlistglmZn[[2]])

modlistglmZn[[3]] <- glm(Zn ~ heat_binglog300, data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[3]])
plot(modlistglmZn[[3]])

modlistglmZn[[4]] <- glm(Zn ~ nlcd_imp_ps, data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[4]])
plot(modlistglmZn[[4]])

modlistglmZn[[5]] <- glm(Zn ~ heatbustransitlog300, data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[5]])
plot(modlistglmZn[[5]])

modlistglmZn[[6]] <- glm(Zn ~ heatOSMSPDlog100, data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[6]])
plot(modlistglmZn[[6]])

modlistglmZn[[7]] <- glm(Zn ~ heatOSMSPDlog100 + I(heatOSMSPDlog100^2), data = pollutfieldclean_cast, family=gaussian(link='log'))
ggplot(pollutfieldclean_cast, aes(x=heatOSMSPDlog100)) + 
  geom_point(aes(y=Zn)) +
  geom_line(data = predf,aes(y=fitted(modlistglmZn[[7]])), color='red', size=1.5) +
  geom_smooth(aes(y=Zn), span=1) +
  theme_classic()
summary(modlistglmZn[[7]])
plot(modlistglmZn[[7]])

modlistglmZn[[8]] <- glm(Zn ~ heatOSMgradientlog200, data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[8]])
plot(modlistglmZn[[8]])

#------ 9. glm(Zn, gaussian log link)  - Multiparameter models --------
modlistglmZn[[9]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300, data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[9]])
plot(modlistglmZn[[9]])

modlistglmZn[[10]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300 + heatOSMSPDlog100, 
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[10]])
plot(modlistglmZn[[10]])

modlistglmZn[[11]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100, 
                         data = pollutfieldclean_cast, family=gaussian(link='log', variance='identity'))
summary(modlistglmZn[[11]])
plot(modlistglmZn[[11]])
theme_set(theme_sjplot())
plot_model(modlistglmZn[[11]], type='pred', terms = c('heat_binglog300', 'heatOSMSPDlog100'))


modlistglmZn[[12]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                           heatOSMAADTlog300, 
                         data = pollutfieldclean_cast, family=gaussian(link = "log"))
summary(modlistglmZn[[12]])
plot(modlistglmZn[[12]])

modlistglmZn[[13]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300*heatOSMAADTlog300, 
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[13]])
plot(modlistglmZn[[13]])

modlistglmZn[[14]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300*heatOSMSPDlog100 + 
                           heatOSMAADTlog300, 
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[14]])
plot(modlistglmZn[[14]])


modlistglmZn[[15]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300*heatOSMgradientlog200 + 
                           heatOSMSPDlog100, 
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[15]])
plot(modlistglmZn[[15]])


modlistglmZn[[16]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                           heatOSMSPDlog100*heatOSMgradientlog200, 
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[16]])
plot(modlistglmZn[[16]])


modlistglmZn[[17]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                           heatOSMSPDlog100*heatOSMgradientlog200 + nlcd_imp_ps, 
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[17]])
plot(modlistglmZn[[17]])

modlistglmZn[[18]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                           heatOSMSPDlog100*heatOSMgradientlog200 + nlcd_imp_ps:heatOSMSPDlog100, 
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[18]])
plot(modlistglmZn[[18]])


modlistglmZn[[19]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300*nlcd_imp_ps + 
                           heatOSMSPDlog100*heatOSMgradientlog200, 
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[19]])
plot(modlistglmZn[[19]])


modlistglmZn[[20]] <- glm(Zn ~ heatbustransitlog300 + heat_binglog300 + 
                           heatOSMSPDlog100 + nlcd_imp_ps, 
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[20]])
plot(modlistglmZn[[20]])


#Build models without bus transit
modlistglmZn[[21]] <- glm(Zn ~ heatOSMAADTlog300 + heat_binglog300 + heatOSMSPDlog100,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[21]])
plot(modlistglmZn[[21]])


modlistglmZn[[22]] <- glm(Zn ~ heatOSMAADTlog300*heat_binglog300,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[22]])
plot(modlistglmZn[[22]])


modlistglmZn[[23]] <- glm(Zn ~ heatOSMSPDlog100*heat_binglog300,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[23]])
plot(modlistglmZn[[23]])


modlistglmZn[[24]] <- glm(Zn ~ heatOSMAADTlog300*heat_binglog300 + heatOSMSPDlog100,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[24]])
plot(modlistglmZn[[24]])


modlistglmZn[[25]] <- glm(Zn ~ heatOSMAADTlog300*heat_binglog300*heatOSMSPDlog100,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[25]])
plot(modlistglmZn[[25]])

modlistglmZn[[26]] <- glm(Zn ~ heatOSMAADTlog300 + heat_binglog300*heatOSMgradientlog200,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[26]])
plot(modlistglmZn[[26]])

modlistglmZn[[27]] <- glm(Zn ~ heatOSMAADTlog300 + heat_binglog300 + nlcd_imp_ps,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[27]])
plot(modlistglmZn[[27]])


#Build models without congestion
modlistglmZn[[28]] <- glm(Zn ~  heatbustransitlog300 + heatOSMAADTlog300,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[28]])
plot(modlistglmZn[[28]])


modlistglmZn[[29]] <- glm(Zn ~  heatbustransitlog300 + heatOSMSPDlog100,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[29]])
plot(modlistglmZn[[29]])


modlistglmZn[[30]] <- glm(Zn ~  heatbustransitlog300 + heatOSMSPDlog100*heatOSMgradientlog200,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[30]])
plot(modlistglmZn[[30]])


modlistglmZn[[31]] <- glm(Zn ~  heatbustransitlog300 + heatOSMSPDlog100*heatOSMgradientlog200 +
                           nlcd_imp_ps,
                         data = pollutfieldclean_cast, family=gaussian(link='log'))
summary(modlistglmZn[[31]])
plot(modlistglmZn[[31]])


#------ 10. glm(Zn, gaussian log link)  - Make latex model summary table ----
vnum <- max(sapply(modlistglmZn, function(mod) {length(mod$coefficients)}))
model_summaryglmZn<- as.data.table(
  ldply(modlistglmZn, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
model_summaryglmZn <- model_summaryglmZn[order(as.numeric(AICc), as.numeric(MAE)),]
cat(latex_format(model_summaryglmZn), file = file.path(moddir, 'modeltable_glmZn.tex'))
setwd(moddir)
texi2pdf('modeltable_glmZn.tex')

#------ 10. glm(Zn, gaussian log link)  - Make latex model summary table when excluding outliers----
model_summaryglmZn_nooutliers <- as.data.table(
  ldply(modlistglmZn, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50,
                                                             remove_outliers = 'outliers',
                                                             labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)])}))

model_summaryglmZn_nooutliers <- model_summaryglmZn_nooutliers[order(as.numeric(AICc), as.numeric(MAE)),]
cat(latex_format(model_summaryglmZn_nooutliers), file = file.path(moddir, 'modeltable_glmZn_nooutliers.tex'))
setwd(moddir)
texi2pdf('modeltable_glmZn_nooutliers.tex')

#Could add Synthetic index ~ PC or Individual pollutants ~ PC but could just add difficulty




####################################### JUNK #############################################################


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


