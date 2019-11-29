#Author: Mathis Messager
#Contact information: messamat@uw.edu
#Creation date: July 2018
#Purpose: import, format, merge, inspect, and model development of field XRF data

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
library(mgcv)
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
library(ggradar) #devtools::install_github("rmadillo/ggradar")
#download : "https://github.com/ricardo-bion/ggtech/blob/master/Circular%20Air-Light%203.46.45%20PM.ttf" and copy to "C:/WINDOWS/Fonts"
extrafont::font_import(pattern = 'Circular', prompt=FALSE)
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
  
  #Studentized deleted residuals vs. leverage
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
  
  if (length(attr(terms(model), 'term.labels')) == 1) { #If the model is a simple linear regression
    regoutlier <- ggplot(df[Elem == chem,],
                         aes_string(x=names(model$model)[2], y=names(model$model)[1])) + 
      geom_text(aes(label=paste0(SiteID, Pair), color=factor(get(flagcol))), size=5) +  
      geom_smooth(method='lm') + 
      labs(x = paste0(names(model$model)[2], chem), y= paste0(names(model$model)[1], chem)) + 
      theme_classic()
  } else { #If model is multiple linear regression
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

# mod <- modlistCu[[15]]
# for (mod in modlistCu) {
#   print(mod$call)
#   regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=1, remove_outliers = 'outliers & leverage')
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
  sigcoefs <- which(sumod$coefficients[,'Pr(>|t|)']<0.10)
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
    coord_fixed(ratio=1) + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
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
                            '\\usepackage[paperheight=10in,paperwidth=35in,margin=0.1in,headheight=0.0in,footskip=0.5in,includehead,includefoot]{geometry}',
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

GAMrescheck <- function(model) {
  print(summary(model))
  print(paste0("AICc:", AICc(model)))
  par(mfrow=c(2,2))
  gam.check(model)
  par(mfrow=c(1,1))
}

GAMmultiplot <- function(model) {
  pdim <- ceiling(sqrt(length(attr(model$terms, "term.labels"))))
  par(mfrow=c(pdim, pdim))
  plot(model,residuals=TRUE,shade=T, cex=6)
  par(mfrow=c(1,1))
}

GAMinfluence <- function(mgcvmodel, data) {
  gamite <- ldply(1:nrow(data), function(index) {
    subdat <- data[-index,]
    submod <- mgcv::gam(formula = mgcvmodel$formula, family=mgcvmodel$family, data=subdat)
    influence <- sum((fitted(mgcvmodel)[-index] - fitted(submod))^2)
    return(cbind(data[index, 'SiteIDPair'], influence))
  })
  return(gamite)
}
###########################################################################################################################################

# 1. Import data 
############################################################################################################################################
########### ---- A. Import and format field data ---- ####
fieldata <- read.csv(file.path(datadir,"field_data/field_data_raw_20190430_edit.csv"))
colnames(fieldata)[1] <- 'SiteID'
fieldata$SiteID <- as.character(fieldata$SiteID)
fieldata$Date <- as.Date(as.character(fieldata$Date), format='%m/%d/%Y')
fieldata_sel <- fieldata[!is.na(fieldata$XRFmin) & !is.na(fieldata$SiteID),] #Remove extraneous sites with no XRF data or just for TNC tour
fieldata_format <- data.frame()
#Create separate records for each XRF measurement (rather than one record with XRFmin and XRFmax)
for (row in seq(1,nrow(fieldata_sel))) {
  extract <- fieldata_sel[row,]
  for (xrf in seq(fieldata_sel[row,'XRFmin'], fieldata_sel[row,'XRFmax'])){
    #print(xrf)
    extract$XRFID <- xrf
    fieldata_format <- rbind(fieldata_format, extract)
  }
}
#Separate summer 2018 and spring 2019 records as XRF ids have been reset and are thus duplicates
fieldata_format[fieldata_format$Date<'2019/01/01', 'season'] <- 'summer2018'
fieldata_format[fieldata_format$Date>'2019/01/01', 'season'] <- 'spring2019'


########### ---- B. Import and format field XRF deconvolution results ---- ####
fieldXRF_summer2018 <- deconvolution_import(file.path(datadir, 'XRF20190501/PostBrukerCalibration/deconvolutionresults_XRF12_245_20181215'),
                                 idstart = 41, idend= 43)
fieldXRF_spring2019 <- deconvolution_import(file.path(datadir, 'XRF20190501/2019_MarchApril/deconvolutionresults_XRF39_161_20190513'),
                                 idstart = 41, idend= 43)

fieldXRF_summer2018[,'season'] <- 'summer2018'
fieldXRF_spring2019[,'season'] <- 'spring2019'

fieldXRF <- rbind(fieldXRF_summer2018, fieldXRF_spring2019)

########### ---- C. Import raw XRF lab data for inspection ---- ####
labxrf_list <- data.table(filename = grep('\\.csv$', list.files(file.path(datadir, 'XRF20190501/PelletMeasurements')) ,value=T))
for (i in labxrf_list$filename) {
  print(i)
  xrfrec <- read.csv(file.path(datadir, 'XRF20190501/PelletMeasurements', i))
  labxrf_list[filename == i, `:=`(cps = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Valid Count Last Packet',])),
                                  c_cum = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Valid Accumulated Counts',])),
                                  duration = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Duration Time',]))
  )]
}
setDT(labxrf_list)[, `:=`(c_cum_ps = c_cum/duration,
                          SiteID = regmatches(labxrf_list$filename, regexpr('[0-9A-Z]+[_][0-9]', labxrf_list$filename)))]
labxrf_list[, season := 'summer2018']

#List of initially problematic samples that had to be re-measured
# problem_recs <- setDT(labxrf_list)[SiteID %in% c('3A_1', '3A_2','6B_2','7A_1', '7A_2','19A_1', '19A_2', '20A_1',
#                                                  '20B_1', '20B_2','22A_1', '22A_2', '28B_2', '32A_1', '33B_1', 
#                                                  '33B_2','34A_1','36A_2', '37A_1', '37A_2','40A_1', '42A_1', 
#                                                  '42B_1','42B_2','44A_1', '44A_2', '45A_1', '45A_2','46A_1','46A_2',
#                                                  '48A_1','49A_1', '49A_2','53A_1', '53A_2', '54B_2','55B_1', '55B_2',
#                                                  '57A_1', '57A_2','59A_1', '59A_2', '61A_1', '61A_2'),]
labxrf_list[c_cum_ps > 90000,] #[!(labxrf_list[c_cum_ps > 90000, SiteID] %in% problem_recs$SiteID)]

########### ---- D. Import and format lab XRF deconvolution results ---- ####
labXRF <- deconvolution_import(file.path(datadir, 'XRF20190501/PelletMeasurements/deconvolutionresults_labXRF_20181215'),
                               idstart = 29, idend= 33)
labXRF[, season := 'summer2018']
########### ---- E. Look at determinants of signalratio (noise) ----
fieldXRFsummary <- fieldXRF[,list(meannet = mean(Net),
                                  signalratio = mean(Net/Backgr.)), .(Element, Line, Energy.keV)]
keVmaxnet <- unique(fieldXRFsummary) %>% #Only keep Element-Line-Energy characteristics of each element
  setkey(meannet) %>% #Classify from low to high energy 
  .[,.SD[.N], by=Element] #Only keep line with the highest average net photon count
keVmineV <- unique(fieldXRFsummary) %>% #Only keep Element-Line-Energy characteristics of each element
  setkey(Energy.keV)%>% #Classify from low to high energy 
  .[,.SD[1], by=Element] #Only keep line with the lowest average net photon count


ggplot(keVmaxnet, aes(x=log(meannet+1), y=signalratio)) +
  geom_text(aes(label=Element)) +
  scale_y_log10() +
  geom_smooth(span=1)

ggplot(keVmaxnet, aes(x=Energy.keV, y=signalratio)) +
  geom_text(aes(label=Element)) +
  scale_y_log10() +
  geom_smooth(span=1)

ggplot(keVmineV, aes(x=Energy.keV, y=signalratio)) +
  geom_text(aes(label=Element)) +
  scale_y_log10() +
  geom_smooth(span=1)

########### ---- F. Import ICP-OES data ---- ####
ICPdat <- read.xlsx(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'),
                    1, stringsAsFactors = F)
ICPthresholds <- read.xlsx(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'),
                           2, stringsAsFactors = F)
rownames(ICPthresholds) <- gsub('\\W', '', ICPthresholds[,1]) #Clean sample IDs
ICPthresholds_format <- as.data.table(t(ICPthresholds[,-1])) %>%
  .[, Elem := colnames(ICPthresholds)[-1]]

########### ---- G. Import GIS data (including pollution variables) ---- ####
trees <- as.data.table(readOGR(dsn = file.path(resdir, 'pollution_variables.gdb'), layer = 'XRFsites_aea'))
#summary(trees)
heatcols <- colnames(trees)[grep('heat', colnames(trees))]
trees[, (heatcols) := lapply(.SD, function(x){x[is.na(x)] <- 0; x}), .SDcols = heatcols]
colnames(trees)[1] <- 'SiteID'
trees <- trees[!(SiteID %in% c('NA', NA)),]
trees[, NLCD_reclass_final_PS := as.factor(NLCD_reclass_final_PS)]

########### ---- H. Define elements associated with car traffic ---- ####
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
fieldXRFcast <- dcast(setDT(fieldXRF), XRFID+season~Element, value.var='Net', fun.aggregate=sum) 
fieldXRFcast[, XRFID := as.numeric(gsub('[_]', '', XRFID))] #Format site number
fieldXRFcast[fieldXRFcast < 0] <- 0 #Floor negative net photon count to 0

labXRFcast <- dcast(setDT(labXRF), XRFID+season~Element, value.var='Net', fun.aggregate=sum)
labXRFcast[labXRFcast < 0] <- 0

# ---- 2. Normalize data by Rhodium photon count for field and lab results ----
fieldXRFcastnorm <- fieldXRFcast[, lapply(.SD, function(x) {x/Rh}), by = .(XRFID, season)]
labXRFcastnorm <- labXRFcast[, lapply(.SD, function(x) {x/Rh}), by = .(XRFID, season)]

# ---- 3. Merge datasets: lab XRF + field XRF + field variables ----
fieldt <- setDT(fieldata_format)[fieldXRFcastnorm, on=.(XRFID, season)]

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
write.dbf(field_artaxmean, 'field_artaxmean_20190501.dbf')

lab_artaxmean <- lab_artaxstats[, .SD, .SDcols = c(1,2, grep('mean', colnames(lab_artaxstats)))] 
colnames(lab_artaxmean) <- gsub('_mean', '', colnames(lab_artaxmean))

# ---- 7. Check field XRF data distribution by element then transform and standardize ----
#Use Tukey ladder computation to determine each variable's transformation (only on numeric columns with more than 1 unique observation)
fieldXRF_format[Elem %in% fieldXRF_format[, length(unique(mean))>1, by=Elem][V1==T, Elem],
                transmean := transformTukey_lambda(mean, start = -2.5, end = 2.5,int = 0.025, 
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
treesmelt <- melt(trees[!is.na(SiteID),], id.vars=grep('heat|NLCD_imp', colnames(trees),
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
############################################################################################################################################


# 3. Inspect data and remove outliers
######################  ---- A. Field XRF ---- ###########
str(fieldXRF_format)
# ---- 1. Assess within-tree field variability ----
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

# #Plot relationship between elemental concentration and cv
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
            axis.text.x = element_text(angle=90),
            axis.line = element_line(color='black'))
  )
  dev.off()
}

# ---- 2. Assess within-site field variability ---- 
#TO DO: CREATE RANDOM PALETTE
field_artaxmeansite <- fieldt[,lapply(.SD, mean, na.rm=TRUE), by=c('SiteID'), .SDcols=31:59]
artaxsdsite <- fieldt[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), .SDcols=31:59]
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

# for (elem in unique(fieldXRF_format[!(Elem %in% c('Rh','Pd','Ar')), Elem])) {
#   print(elem)
#   png(file.path(inspectdir, paste0('fieldXRF_withinsite_',elem,'.png')), width = 20, height=12, units='in', res=300)
#   print(
#     ggplot(fieldXRF_format[Elem == elem,], 
#            aes(x=SiteID, y = mean, fill=SiteID)) + 
#       geom_line(aes(group=SiteID), color='black') +
#       geom_point(size=5, colour='black', pch=21, alpha=0.75) +
#       geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color=SiteID)) +
#       labs(x='Element', y='Mean net photon count') + 
#       theme_bw() + 
#       theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
#             strip.text = element_text(size=14),
#             panel.border = element_blank(),
#             axis.line = element_line(color='black'))
#   )
#   dev.off()
# }

# ---- 3. Exclude columns that have high within-tree variability or that are unrelated to traffic pollution ----
fieldXRF_format[, `:=`(meanCV=mean(cv, na.rm=T),
                       sdCV = sd(cv, na.rm=T)), by=Elem]
excol <- unique(c(fieldXRF_format[meanCV>0.5,unique(Elem)], 
                  fieldXRF_format[!(fieldXRF_format$Elem %in% union(brakelining_elem, tire_elem)),Elem],
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

ggplot(ICPmeanR2, aes(x=signalratio, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

summary(loess(R2~meanXRF, data=ICPmeanR2, span=1))
summary(loess(R2~meanXRF+signalratio, data=ICPmeanR2, span=1))

# ---- 11. Univariate relationship to lab XRF for the purpose of outlier detection ----
#Outlier diagnostic plots
labfieldmerge[, labfield_flags := 0]

for (chem in unique(labfieldmerge$Elem)) {
  print(chem)
  labfield_lm <- lm(lab_mean ~ field_mean, data = labfieldmerge[Elem == chem & !is.na(lab_mean),])
  ggsave(file.path(inspectdir, paste0('fieldlab_regoutliers', chem, '.png')),
         chemregdiagnostic_custom(labfield_lm, labfieldmerge[!is.na(lab_mean),], chem,  flagcol = 'labfield_flags', tresids = TRUE),
         width = 20, height=12, units='in', dpi=300)
  labfieldmerge[Elem == chem & !is.na(lab_mean), labfieldR2 := summary(labfield_lm)$adj.r.squared]
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

ggplot(labmeanR2, aes(x=signalratio, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

# ---- 12. Univariate relationship to pollution predictors for the purpose of outlier detection ----
#Outlier diagnostic plots
pollutfieldmerge <- pollutfieldmerge[!is.na(Elem),]
pollutfieldmerge[, pollutfield_flags := 0]
for (chem in unique(pollutfieldmerge$Elem)) {
  print(chem)
  pollutfield_lm <- lm(mean ~ heatbing1902log300proj*heatsubAADTlog300*heatsubslopelog300*nlcd_imp_ps +
                         heatbustransitlog300 + heatbustransitlog300:heatsubslopelog300, 
                       data = pollutfieldmerge[Elem == chem ,])
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

ggplot(pollutmeanR2, aes(x=signalratio, y=R2)) + 
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
flagelems <- fieldXRF_formatflags[, min(ICPfieldR2, labfieldR2, pollutfieldR2, na.rm=T)>0.10, by=Elem][
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
fieldXRF_inspect <- unique(fieldXRF_formatflags[(ICPfield_flags_sum>3) | 
                                                  (labfield_flags_sum > 3) | 
                                                  (pollutfield_flags_sum > 3), 
                                                .(SiteID, Pair)])
outliertrees <- trees[fieldXRF_inspect, on=c('SiteID', 'Pair')]
outlierlocs <- SpatialPointsDataFrame(coords = data.frame(outliertrees$POINT_X, outliertrees$POINT_Y),
                                      data= as.data.frame(outliertrees))
#View(outlierlocs@data)
# leaflet(data = outlierlocs) %>% addTiles() %>%
#   addMarkers(clusterOptions = markerClusterOptions(),
#              popup = ~paste0(SiteID, Pair))


###################### ---- B. Lab XRF ---- ####
# ---- 1. Assess within-pellet lab variability ----
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


# ---- 2. Assess within-site lab variability ---- 
lab_artaxmeansite <- labdt[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID'), .SDcols=30:52]
artaxsdsite <- labdt[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), .SDcols=30:52]
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

ggplot(ICPmeanR2, aes(x=signalratio, y=R2)) + 
  geom_text(aes(label=Elem)) +
  scale_x_log10() + 
  geom_smooth(span=1) + 
  theme_classic()

summary(loess(R2~meanXRF, data=ICPmeanR2, span=1))
summary(loess(R2~meanXRF+signalratio, data=ICPmeanR2, span=1))

# ---- 5. Univariate relationship to pollution predictors for the purpose of outlier detection ----
#Outlier diagnostic plots
pollutlabmerge[, pollutlab_flags := 0]
pollutlabmerge <- pollutlabmerge[!is.na(Elem),]
for (chem in unique(pollutlabmerge$Elem)) {
  print(chem)
  pollutlab_lm <- lm(mean ~  heatbing1902log300proj*heatsubAADTlog300*heatsubslopelog300*nlcd_imp_ps +
                       heatbustransitlog300 + heatbustransitlog300:heatsubslopelog300, 
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

ggplot(pollutmeanR2, aes(x=signalratio, y=R2)) + 
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
labXRF_inspect <- unique(labXRF_formatflags[(ICPlab_flags_sum>3) | 
                                              (labfield_flags_sum > 3) | 
                                              (pollutlab_flags_sum > 3), 
                                            .(SiteID, Pair)])
outliertrees <- trees[labXRF_inspect, on=c('SiteID', 'Pair')]
outlierlocs <- SpatialPointsDataFrame(coords = data.frame(outliertrees$POINT_X, outliertrees$POINT_Y),
                                      data= as.data.frame(outliertrees))

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
pollutICPmerge <- pollutICPmerge[!is.na(Elem),]
for (chem in unique(pollutICPmerge[!is.na(Elem), Elem])) {
  print(chem)
  pollutICP_lm <- lm(ICP ~  heatbing1902log300proj*heatsubAADTlog300*heatsubslopelog300*nlcd_imp_ps +
                       heatbustransitlog300 + heatbustransitlog300:heatsubslopelog300, 
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
                flagsum_ICP := sum(ICPfield_flags, ICPlab_flags, pollutICP_flags, na.rm=T), 
                by=.(SiteID, Pair)]

#Check total number of flags by each type of flag
ICP_formatflags[Elem %in% flagelems,
                `:=`(ICPfield_flags_sum = sum(ICPfield_flags, na.rm=T),
                     ICPlab_flags_sum = sum(ICPlab_flags, na.rm=T),
                     pollutICP_flags_sum = sum(pollutICP_flags, na.rm=T)),
                by=.(SiteID, Pair)]

ICP_formatflags_u <- unique(
  ICP_formatflags[!is.na(flagsum_ICP),
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
ICP_inspect <- unique(ICP_formatflags[(ICPfield_flags_sum>3) | 
                                        (ICPlab_flags_sum > 3) | 
                                        (pollutICP_flags_sum > 3), 
                                      .(SiteID, Pair)])
outliertrees <- trees[ICP_inspect, on=c('SiteID', 'Pair')]
outlierlocs <- SpatialPointsDataFrame(coords = data.frame(outliertrees$POINT_X, outliertrees$POINT_Y),
                                      data= as.data.frame(outliertrees))
#View(outlierlocs@data)
# leaflet(data = outlierlocs) %>% addTiles() %>%
#   addMarkers(clusterOptions = markerClusterOptions(),
#              popup = ~paste0(SiteID, Pair))

###################### ---- D. Compile all flags and make a flag matrix/heatmap then inspect data and decide on their fate ---- ####
# ---- 1. Compile flags and make heatmap ----
allflags <- ICP_formatflags_u[labXRF_formatflags_u, on=.(SiteID, Pair)][
  fieldXRF_formatflags_u, on=.(SiteID, Pair)] 
allflags[, (grep('i[.]', colnames(allflags), value=T)) := NULL][
  , flagpartsum := sum(flagsum_ICP, flagpartsum_lab, flagpartsum_field, na.rm=T), by=.(SiteID, Pair)][
    , SiteIDPair := factor(paste0(SiteID, Pair), levels = unique(paste0(SiteID, Pair)[order(-flagpartsum)]))]

colnames(allflags)
#levels(allflags_melt$variable)
allflags_melt <- melt(allflags, id.vars = c('SiteIDPair', 'SiteID', 'Pair')) %>%
  .[, variable := factor(gsub('[_]|flag*|sum', '', variable), 
                         levels=gsub('[_]|flag|sum', '',
                                     c("flagsum_lab", "flagsum_field", "flagpartsum",
                                       "flagsum_ICP" ,"flagpartsum_lab", "flagpartsum_field",  
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
        axis.text.x = element_text(angle=90, vjust=-0.07))

# ---- 2. Analysis of the main outliers ----
#112A: South Beacon Hill, forested area across from Boeing Field/I5 (closest to road)
#     5 measurements on this tree. High within tree variability but no clear outlier point. 
#     pollut-field outlier for Cr, Fe, Ti, Zn, Zr (model overestimates pollution; could be due to protection from surrounding trees and topo.)

#111A: Interbay, Magnolia Bridge exit
#     High within tree variability with one apparent XRF sample with much higher value than the two others
#     pollut-field outlier for Cr, Fe, Ti, Zn (model underestimates pollution)

#114A and 114B: 15th Ave NE & 43rd St, UW campus by bus stop
#     Not abnormal within-tree variability for outlying elements, similar between 114A and 114B in both range and mean
#     pollut-field outlier for Cu, Fe, Zn (model underestimates pollution, particularly for Cu)

#117A and 117B: 1st Ave S, Southern end of bridge over Union Pacific cargo yard in Industrial District
#     Normal within-tree variability
#     pollut-field outlier for Cr, Fe, Zn, Zr (model underestimates pollution)

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
pollutfieldclean_cast[, SiteIDPair := paste0(SiteID, Pair)]

#Create vector of columns - 1.G. and  3.A.3. for column selection
"***NOTE: consider removing all elems with net/background photon count < 0.1"
elemcols <- colnames(pollutfieldclean_cast)[colnames(pollutfieldclean_cast) %in% periodicTable$symb] %>%
  setdiff('Rh')
elemcols_sub <- elemcols %>%
  setdiff(c(excol, excol2))

#Scatterplot matrix 
outcatmat <- file.path(inspectdir, 'FieldElem_FieldElem_matrix.png')
if (!file.exists(outcatmat)) {
  png(outcatmat , width = 24, height=24, units='in', res=300)
  ggscatmat(as.data.frame(pollutfieldclean_cast[, elemcols, with=F]),
            alpha=0.7)
  dev.off()
}


#Correlation heatmap
outfieldheatmat <- file.path(inspectdir, 'corheatmap_FieldElem_FieldElem.png')
if (!file.exists(outfieldheatmat)) {
  png(outfieldheatmat, width = 12, height=12, units='in', res=300)
  corr_heatmap(pollutfieldclean_cast[, elemcols, with=F])
  dev.off()
}

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

#Create a PCA biplot matrix where all components are graphed against each other
pcabiplot_grid(pollutfieldclean_cast, nPCs = 5, cols = elemcols_sub,
               idcols = c('SiteID', 'Pair'),  scload = 3)

#Inspect components to decide which ones to predict
loadings(pca)
#There are no very distinct patterns - everything is pretty correlated:
#First component: Cr, Fe, Ti, Zn, Zr (0.34-0.36);
#                 Ca, Cu a little less (~0.30)
#                 Ni, and Pb a little less (0.20-25);
#                 Sr, Mn, K, Ba not much (0.13-0.16)
#Second component:K loads the strongest (0.53)
#                 Sr, Ni, Ca, Ba load moderately (0.30-0.37)
#                 Pb and Mn load positively but not much (0.16)
#                 Cu is ~ 0
#                 Cr and Ti load moderately and negatively (~-0.15)
#                 Zn, Zr, Fe load negatively < -20

#Select individual elements to predict separately: Zn, Cu, Pb

#---- B. Create a synthetic sum-based pollution index (averaging centered and standardized elements)   ----
#Check that all can be transformed using the same transformation then transform
selcols <- c('Cu', 'Pb', 'Zn') #In three main different groups in PCA
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

#Rescale all pollutcols
pollutfieldclean_cast[, (heatcols) := lapply(.SD, function(x) 100*x/max(x)), .SDcols = heatcols]

#Correlation heatmap
outpollutheatmat <- file.path(inspectdir, 'corheatmap_PollutionDrivers_PollutionDrivers.png')
if (!file.exists(outpollutheatmat)) {
  png(outpollutheatmat, width = 30, height=30, units='in', res=300)
  corr_heatmap(pollutfieldclean_cast[, pollutcols[-which(pollutcols=='NLCD_reclass_final_PS')], with=F]) + 
    scale_x_discrete(labels=heatlab) + 
    scale_y_discrete(labels=heatlab)
  dev.off()
}

#---- D. All pollution drivers against field XRF all elems ----
outfieldpollutheatmat <- file.path(inspectdir, 'corheatmap_PollutionDrivers_FieldElem.png')
if (!file.exists(outfieldpollutheatmat)) {
  png(outfieldpollutheatmat,width = 40, height=35, units='in', res=300)
  corr_heatmap(xmat=pollutfieldclean_cast[, pollutcols[-which(pollutcols=='NLCD_reclass_final_PS')], with=F],
               ymat=pollutfieldclean_cast[, c(elemcols, 'pollution_index', 'cbpollution_index'), with=F],
               clus = FALSE) +
    scale_y_discrete(labels=heatlab) + 
    theme(text=element_text(size=22))
  dev.off()
}

xmat <- pollutfieldclean_cast[, pollutcols[-which(pollutcols=='NLCD_reclass_final_PS')], with=F]
ymat <- pollutfieldclean_cast[, c(elemcols, 'pollution_index', 'cbpollution_index'), with=F]
cordf <- round(cor(x=xmat, y=ymat),2) %>%
  t() %>%
  as.data.frame()
cordf$group <- rownames(cordf)
cordf <- cordf[, c(ncol(cordf),
                         1:(ncol(cordf)-1))]

spidermetals <- c('Fe',  'Cu', 'Zn','Pb')
spidercols <-  colnames(cordf)[!(grepl('(^heatAADT.*)|(log[0-9]{2}$)|.*mean.*', colnames(cordf))) &
                                 grepl('(.*log.*)|.*nlcd.*|group', colnames(cordf))]
spiderlabels <- gsub('(heat_?)|log', '',
                     gsub('bing[1-9]*', 'Congestion ',
                          gsub('subslope', 'Gradient ',
                               gsub('subAADT', 'Volume ',
                                    gsub('SPD', 'Speed ',
                                         gsub('bustransit', 'Transit ',
                                              gsub('nlcd_imp_ps.*', 'Imperv.',
                                                   spidercols)))))))
funcCircleCoords <- function(center = c(0,0), r = 1, npoints = 100){
  #Adapted from Joran's response to http://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
grid.mid = 0.5
grid.min = 0
grid.max = 1
gridline75 <- funcCircleCoords(
  c(0,0),
  0.75 +abs(grid.min - ((1/9)*(grid.max-grid.min))),
  npoints = 360)
gridline25 <- funcCircleCoords(
  c(0,0),
  0.25 +abs(grid.min - ((1/9)*(grid.max-grid.min))),
  npoints = 360)
# plot.data <- cordf[cordf$group %in% spidermetals,
#                    spidercols[order(spidercols)]]
# plot.data <- as.data.frame(plot.data)
# plot.data[,1] <- as.factor(as.character(plot.data[,1]))
# names(plot.data)[1] <- "group"
# var.names <- colnames(plot.data)[-1]  #'Short version of variable names
# plot.data.offset <- plot.data
# plot.data.offset[,2:ncol(plot.data)]<- plot.data[,2:ncol(plot.data)]+abs(centre.y)
# #print(plot.data.offset)
# # (b) convert into radial coords
# CalculateGroupPath <- function(df) {
#   #Converts variable values into a set of radial x-y coordinates
#   #Code adapted from a solution posted by Tony M to
#   #http://stackoverflow.com/questions/9614433/creating-radar-chart-a-k-a-star-plot-spider-plot-using-ggplot2-in-r
#   #Args:
#   #  df: Col 1 -  group ('unique' cluster / group ID of entity)
#   #      Col 2-n:  v1.value to vn.value - values (e.g. group/cluser mean or median) of variables v1 to v.n
#
#   path <- df[,1]
#
#   ##find increment
#   angles = seq(from=0, to=2*pi, by=(2*pi)/(ncol(df)-1))
#   ##create graph data frame
#   graphData= data.frame(seg="", x=0,y=0)
#   graphData=graphData[-1,]
#
#   for(i in levels(path)){
#     pathData = subset(df, df[,1]==i)
#     for(j in c(2:ncol(df))){
#       #pathData[,j]= pathData[,j]
#
#
#       graphData=rbind(graphData, data.frame(group=i,
#                                             x=pathData[,j]*sin(angles[j-1]),
#                                             y=pathData[,j]*cos(angles[j-1])))
#     }
#     ##complete the path by repeating first pair of coords in the path
#     graphData=rbind(graphData, data.frame(group=i,
#                                           x=pathData[,2]*sin(angles[1]),
#                                           y=pathData[,2]*cos(angles[1])))
#   }
#   #Make sure that name of first column matches that of input data (in case !="group")
#   colnames(graphData)[1] <- colnames(df)[1]
#   graphData #data frame returned by function
# }
# group <-NULL
# group$path <- CalculateGroupPath(plot.data.offset)
#
# #Create palette
# library(RColorBrewer)
# interspec <- colorRampPalette(brewer.pal(11, 'Spectral'))
# intdat<- data.frame(y=round(group$path$y,2))
# datrange <- seq(min(intdat), max(intdat), 0.01)
# coldf <- data.frame(y = datrange, ycol=interspec(length(datrange)))
# colvalues <- merge(intdat, coldf, on='y', all.y=F)
# colvalues$yfac <- factor(as.character(colvalues$y),
#                          levels= unique(as.character(colvalues$y)[order(-colvalues$y)]))
#
# ggplot(colvalues, aes(x=y, y=y, color=yfac)) +
#   geom_point() +
#   scale_color_manual(values=unique(as.character(colvalues$ycol)))

ggradarplot <- ggradar(plot.data=cordf[cordf$group == 'Zn',
                                       spidercols[order(spidercols)]],
                       grid.min=grid.min, grid.max=grid.max, grid.mid = grid.mid,
                       gridline.mid.linetype = "solid",
                       gridline.mid.colour= 'darkgrey',
                       grid.label.size = 0,
                       group.point.size = 2,
                       group.line.width = 1.1,
                       axis.labels = str_wrap(spiderlabels[order(spidercols)][-1], 5),
                       axis.label.size = 8,
                       axis.label.offset = 1.1,
                       group.colours = '#8CE071',
                       background.circle.transparency=0.15) +
  geom_path(data=gridline75, aes(x=x, y=y), color='darkgrey') +
  geom_path(data=gridline25, aes(x=x, y=y), color='darkgrey') +
  #scale_colour_manual(values=unique(as.character(colvalues$ycol))) +
  #facet_wrap(~group) +
  theme_minimal() +
  theme(text= element_text(size = 18),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position=c(0.15, 0.15),
        plot.margin=unit(c(-1.5, -1.5, -1.5, -1.5),'cm'))
ggradarplot

#pdf(file.path(moddir, 'spiderplot_20190223_1.pdf'), width=9, height=9)
png(file.path(moddir, 'spiderplotZn_20190514_1.png'), width=9, height=9, units='in', res=450)
ggradarplot
dev.off()


"Summary of best univariate predictors for each category
- Zn: 
  - IMPERVIOUSNESS > BING > SPDL > AADT > SLOPE  > TRANSIT
  - AADT: 100 > 200-300 log~pw 
  - SPDL: 50~100~200>300>500
  - slope: 500>300>200>100 and log~pow
  - bing: 300=200>100>500 log~pow
  - transit: 500=300=200 > 100 log~pow
  - imp_mean > imp
"
#---- E. Create a synthetic PCA-based pollution driver index  ----
pollutpca <- PcaClassic(~., data=pollutfieldclean_cast[,selcols[selcols %in% pollutcols],with=F], scale=TRUE) #z-standardize
summary(pollutpca)
screeplot(pollutpca, main="Screeplot: classic PCA", bstick=TRUE) #First PC significant compared to broken stick
ordi.monte(pollutfieldclean_cast[,selcols[selcols %in% pollutcols],with=F],ord='pca',dim=5) #2PCs significant with Monte-carlo test of component significance
#Plot
#Plot
#Plot
#Plot
#Plot
pcabiplot_grid(pollutfieldclean_cast, nPCs = 3, cols = selcols[selcols %in% pollutcols],
               idcols = c('SiteID', 'Pair'),  scload = 3)
loadings(pollutpca)

############################################################################################################################################


# 5. Model selection
############################################################################################################################################
############## Prepare data ###################
pollutfieldclean_cast[, `:=`(heatsubAADTlog100sqrt = sqrt(heatsubAADTlog100),
                             heatsubAADTlog200sqrt = sqrt(heatsubAADTlog200),
                             heatsubAADTlog300sqrt = sqrt(heatsubAADTlog300),
                             heatbing1902log300projsqrt = sqrt(heatbing1902log300proj),
                             heatbing1902log500projsqrt = sqrt(heatbing1902log500proj),
                             heatsubspdlpow100_3sqrt = sqrt(heatsubspdlpow100_3),
                             heatsubspdllog500sqrt = heatsubspdllog500^(1/2),
                             heatsubAADTlog100thd = heatsubAADTlog100^(1/3),
                             heatsubAADTlog200thd = heatsubAADTlog200^(1/3),
                             heatsubAADTlog100frt = heatsubAADTlog100^0.25,
                             heatsubAADTlog200frt = heatsubAADTlog200^0.25,
                             heatsubAADTlog300frt = heatsubAADTlog300^0.25,
                             heatsubAADTlog100fth = heatsubAADTlog100^0.20,
                             heatsubAADTpow100_1frt = heatsubAADTpow100_1^0.25,
                             heatsubAADTpow100_2frt = heatsubAADTpow100_2^0.25,
                             heatbustransitpow100_1sqrt = heatbustransitpow100_1**0.5,
                             heatbustransitpow200_1sqrt = heatbustransitpow200_1**0.5,
                             heatbustransitpow300_1sqrt = heatbustransitpow300_1**0.5,
                             heatbustransitlog200sqrt = heatbustransitlog200**0.5,
                             heatbustransitlog300sqrt = heatbustransitlog300**0.5,
                             heatbustransitpow100_1thd = heatbustransitlog200**(1/3),
                             heatbustransitpow200_1thd = heatbustransitlog300**(1/3),
                             heatbustransitpow100_2thd = heatbustransitlog200**(1/3),
                             heatbustransitpow200_2thd = heatbustransitlog300**(1/3),
                             logZn = log(Zn),
                             logPI = log(pollution_index)
                             )]

pollutfieldclean_cast[, NLCD_reclass_final_PS_edit := NLCD_reclass_final_PS]
pollutfieldclean_cast[!(NLCD_reclass_final_PS %in% c('96', '97', '98', '99')), NLCD_reclass_final_PS_edit := '1']

extraoutliers <- c('51A', '114A', '114B', paste0(64:69, 'A')) #Remove anomalous data points by 15th Ave NE and Luwam's (latter do not appear to be precise enough


#--------------- A. Zn ~ separate predictors ----
#Multiply Zn to make coefficients more readable
modlistZn <- list() #List to hold models
modlistZn[[1]] <- lm(Zn ~ 1, data = pollutfieldclean_cast) #Null/Intercept model
#------ 1. Zn - Single parameter models --------
modlistZn[[2]] <- lm(Zn ~ heatsubAADTlog100, data = pollutfieldclean_cast)
ols_regress(modlistZn[[2]])
#ols_plot_diagnostics(modlistZn[[2]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog100^(1/4), y=Zn, label=paste0(SiteID, Pair),
                                  color = heatbing1902log300proj)) + 
  geom_text() + 
  scale_color_distiller(palette='Spectral')

modlistZn[[3]] <- lm(Zn ~ heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistZn[[3]])
#ols_plot_diagnostics(modlistZn[[3]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902log300proj, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()
ggplot(pollutfieldclean_cast[as.numeric(SiteID) < 63,], aes(x=heatbing1902log300proj, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistZn[[4]] <- lm(Zn ~ nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistZn[[4]])
#ols_plot_diagnostics(modlistZn[[4]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps_mean, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistZn[[5]] <- lm(Zn ~ nlcd_imp_ps, data = pollutfieldclean_cast)
ols_regress(modlistZn[[5]])
#ols_plot_diagnostics(modlistZn[[5]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistZn[[6]] <- lm(Zn ~ heatbustransitlog200, data = pollutfieldclean_cast)
ols_regress(modlistZn[[6]])
#ols_plot_diagnostics(modlistZn[[6]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitlog200, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistZn[[7]] <- lm(Zn ~ heatsubspdlpow100_3, data = pollutfieldclean_cast)
ols_regress(modlistZn[[7]])
#ols_plot_diagnostics(modlistZn[[7]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow100_3, y=Zn, 
                                  label=paste0(SiteID, Pair), color=nlcd_imp_ps_mean)) + 
  geom_text() +
  geom_smooth(method='lm') +
  scale_color_distiller(palette='Spectral')

modlistZn[[8]] <- lm(Zn ~ heatsubslopepow500_1, data = pollutfieldclean_cast)
ols_regress(modlistZn[[8]])
#ols_plot_diagnostics(modlistZn[[8]])
ggplot(pollutfieldclean_cast, aes(x=heatsubslopepow500_1, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

#------ 2. Zn - Multiparameter models --------
modlistZn[[9]] <- lm(Zn ~ heatsubspdlpow100_3 + heatsubAADTlog100, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[9]])
#ols_plot_diagnostics(modlistZn[[9]])
ols_coll_diag(modlistZn[[9]])
ols_correlations(modlistZn[[9]])

modlistZn[[10]] <- lm(Zn ~ heatsubspdlpow100_3 + heatbing1902log300proj, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[10]])
#ols_plot_diagnostics(modlistZn[[10]])
ols_coll_diag(modlistZn[[10]])
ols_correlations(modlistZn[[10]])

modlistZn[[11]] <- lm(Zn ~  heatsubspdlpow100_3 + nlcd_imp_ps_mean, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[11]])
#ols_plot_diagnostics(modlistZn[[11]])
ols_coll_diag(modlistZn[[11]])
ols_correlations(modlistZn[[11]])

modlistZn[[12]] <- lm(Zn ~ heatsubspdlpow100_3 + heatbustransitlog200, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[12]])
#ols_plot_diagnostics(modlistZn[[12]])
ols_coll_diag(modlistZn[[12]])
ols_correlations(modlistZn[[12]])

modlistZn[[13]] <- lm(Zn ~ heatsubspdlpow100_3 + heatsubslopepow500_1, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[13]])
#ols_plot_diagnostics(modlistZn[[13]])
ols_coll_diag(modlistZn[[13]])
ols_correlations(modlistZn[[13]])

modlistZn[[14]] <- lm(Zn ~ heatsubspdlpow100_3*heatsubAADTlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[14]])
#ols_plot_diagnostics(modlistZn[[14]])
ols_coll_diag(modlistZn[[14]])
ols_correlations(modlistZn[[14]])

modlistZn[[15]] <- lm(Zn ~ heatsubspdlpow100_3*heatbing1902log300proj, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[15]])
#ols_plot_diagnostics(modlistZn[[15]])
ols_coll_diag(modlistZn[[15]])
ols_correlations(modlistZn[[15]])

modlistZn[[16]] <- lm(Zn ~ heatsubspdlpow100_3*heatbing1902log300proj + nlcd_imp_ps_mean , 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[16]])
#ols_plot_diagnostics(modlistZn[[16]])
ols_coll_diag(modlistZn[[16]])
ols_correlations(modlistZn[[16]])

modlistZn[[17]] <- lm(Zn ~ heatsubspdlpow100_3 + heatbing1902log300proj + nlcd_imp_ps_mean, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[17]])
#ols_plot_diagnostics(modlistZn[[17]])
ols_coll_diag(modlistZn[[17]])
ols_correlations(modlistZn[[17]])

modlistZn[[18]] <- lm(Zn ~ heatsubspdlpow100_3*nlcd_imp_ps_mean + heatbing1902log300proj,
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[18]])
#ols_plot_diagnostics(modlistZn[[18]])
ols_coll_diag(modlistZn[18])
ols_correlations(modlistZn[[18]])

modlistZn[[19]] <- lm(Zn ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[19]])
#ols_plot_diagnostics(modlistZn[[19]])
ols_correlations(modlistZn[[19]])

modlistZn[[20]] <- lm(Zn ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean + 
                        heatbustransitlog200, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[20]])
#ols_plot_diagnostics(modlistZn[[20]])
ols_coll_diag(modlistZn[20])
ols_correlations(modlistZn[[20]])


#------ 3. Zn - Transformed predictors --------
modlistZn[[21]] <- lm(Zn ~ heatsubAADTlog200sqrt, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZn[[21]])
#ols_plot_diagnostics(modlistZn[[21]])
ols_coll_diag(modlistZn[[21]])
ols_correlations(modlistZn[[21]])

modlistZn[[22]] <- lm(Zn ~ heatsubspdlpow100_3sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[22]])
#ols_plot_diagnostics(modlistZn[[22]])
ols_coll_diag(modlistZn[[22]])
ols_correlations(modlistZn[[22]])

modlistZn[[23]] <- lm(Zn ~ heatbing1902log300projsqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[23]])
#ols_plot_diagnostics(modlistZn[[23]])
ols_coll_diag(modlistZn[[23]])
ols_correlations(modlistZn[[23]])

modlistZn[[24]] <- lm(Zn ~ heatsubAADTlog300sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[24]])
#ols_plot_diagnostics(modlistZn[[24]])
ols_coll_diag(modlistZn[[24]])
ols_correlations(modlistZn[[24]])

modlistZn[[25]] <- lm(Zn ~ heatsubAADTlog100sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[25]])
#ols_plot_diagnostics(modlistZn[[25]])
ols_coll_diag(modlistZn[[25]])
ols_correlations(modlistZn[[25]])

modlistZn[[26]] <- lm(Zn ~ heatsubAADTlog100sqrt + heatsubspdlpow100_3sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[26]])
#ols_plot_diagnostics(modlistZn[[26]])
ols_coll_diag(modlistZn[[26]])
ols_correlations(modlistZn[[26]])

modlistZn[[27]] <- lm(Zn ~ heatsubAADTlog100sqrt*heatsubspdlpow100_3, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[27]])
#ols_plot_diagnostics(modlistZn[[27]])
ols_coll_diag(modlistZn[[27]])
ols_correlations(modlistZn[[27]])

modlistZn[[28]] <- lm(Zn ~ heatsubAADTlog100sqrt*heatsubspdlpow100_3 + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[28]])
#ols_plot_diagnostics(modlistZn[[28]])
ols_coll_diag(modlistZn[[28]])
ols_correlations(modlistZn[[28]])

modlistZn[[29]] <- lm(Zn ~ (heatsubAADTlog100sqrt)*heatsubspdlpow100_3 + heatbing1902log300proj, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[29]])
#ols_plot_diagnostics(modlistZn[[29]])
ols_coll_diag(modlistZn[[29]])
ols_correlations(modlistZn[[29]])

modlistZn[[30]] <- lm(Zn ~ (heatsubAADTlog100sqrt)*heatsubspdlpow100_3 + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubAADTlog100sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[30]])
#ols_plot_diagnostics(modlistZn[[30]])
ols_coll_diag(modlistZn[[30]])
ols_correlations(modlistZn[[30]])


modlistZn[[31]] <- lm(Zn ~ (heatsubAADTlog100sqrt) + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubAADTlog100sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[31]])
#ols_plot_diagnostics(modlistZn[[31]])
ols_coll_diag(modlistZn[[31]])
ols_correlations(modlistZn[[31]])

modlistZn[[32]] <- lm(Zn ~ (heatsubspdlpow100_1) + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubspdlpow100_3, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[32]])
#ols_plot_diagnostics(modlistZn[[32]])
ols_coll_diag(modlistZn[[32]])
ols_correlations(modlistZn[[32]])

modlistZn[[33]] <- lm(Zn ~ (heatsubAADTlog100sqrt) + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubAADTlog100sqrt + nlcd_imp_ps, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[33]])
#ols_plot_diagnostics(modlistZn[[33]])
ols_coll_diag(modlistZn[[33]])
ols_correlations(modlistZn[[33]])

modlistZn[[34]] <- lm(Zn ~ (heatsubAADTlog100sqrt) + heatbing1902log300proj + nlcd_imp_ps +
                        heatbustransitlog200, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[34]])
#ols_plot_diagnostics(modlistZn[[34]])
ols_coll_diag(modlistZn[[34]])
ols_correlations(modlistZn[[34]])

modlistZn[[35]] <- lm(Zn ~ (heatsubAADTlog100sqrt) + heatbing1902log300proj + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[35]])
#ols_plot_diagnostics(modlistZn[[35]])
ols_coll_diag(modlistZn[[35]])
ols_correlations(modlistZn[[35]])

modlistZn[[36]] <- lm(Zn ~ heatsubAADTlog100sqrt + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZn[[36]])
#ols_plot_diagnostics(modlistZn[[36]])
ols_coll_diag(modlistZn[[36]])
ols_correlations(modlistZn[[36]])

#------ 4. Zn - Make latex model summary table ----
vnum <- max(sapply(modlistZn, function(mod) {length(mod$coefficients)}))
model_summaryZn<- as.data.table(
  ldply(modlistZn, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summaryZn, -R2pred, AICc)  
cat(latex_format(model_summaryZn), file = file.path(moddir, 'modeltable_Zn_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_Zn_2019.tex')

#------ 5. Zn - Make latex model summary table when excluding outliers----
vnum <- max(sapply(modlistZn, function(mod) {length(mod$coefficients)}))
model_summaryZn_nooutliers<- as.data.table(
  ldply(modlistZn, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summaryZn_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summaryZn_nooutliers), file = file.path(moddir, 'modeltable_Zn_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_Zn_nooutliers_2019.tex')

#------ 6. logZn - Single parameter models --------
modlistlogZn <- list() #List to hold models
modlistlogZn[[1]] <- lm(logZn ~ 1, data = pollutfieldclean_cast) #Null/Intercept model

modlistlogZn[[2]] <- lm(logZn ~ heatsubAADTlog100, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[2]])
#ols_plot_diagnostics(modlistlogZn[[2]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog100^(1/4), y=Zn, label=paste0(SiteID, Pair),
                                  color = heatbing1902log300proj)) + 
  geom_text() + 
  scale_color_distiller(palette='Spectral')

modlistlogZn[[3]] <- lm(logZn ~ heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[3]])
#ols_plot_diagnostics(modlistlogZn[[3]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902log300proj, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()
ggplot(pollutfieldclean_cast[as.numeric(SiteID) < 63,], aes(x=heatbing1902log300proj, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistlogZn[[4]] <- lm(logZn ~ nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[4]])
#ols_plot_diagnostics(modlistlogZn[[4]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps_mean, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistlogZn[[5]] <- lm(logZn ~ nlcd_imp_ps, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[5]])
#ols_plot_diagnostics(modlistlogZn[[5]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistlogZn[[6]] <- lm(logZn ~ heatbustransitlog200, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[6]])
#ols_plot_diagnostics(modlistlogZn[[6]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitlog200, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistlogZn[[7]] <- lm(logZn ~ heatsubspdlpow100_3, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[7]])
#ols_plot_diagnostics(modlistlogZn[[7]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow100_3, y=Zn, 
                                  label=paste0(SiteID, Pair), color=nlcd_imp_ps_mean)) + 
  geom_text() +
  geom_smooth(method='lm') +
  scale_color_distiller(palette='Spectral')

modlistlogZn[[8]] <- lm(logZn ~ heatsubslopepow500_1, data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[8]])
#ols_plot_diagnostics(modlistlogZn[[8]])
ggplot(pollutfieldclean_cast, aes(x=heatsubslopepow500_1, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text()

#------ 7. logZn - Multiparameter models --------
modlistlogZn[[9]] <- lm(logZn ~ heatsubspdlpow100_3 + heatsubAADTlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[9]])
#ols_plot_diagnostics(modlistlogZn[[9]])
ols_coll_diag(modlistlogZn[[9]])
ols_correlations(modlistlogZn[[9]])

modlistlogZn[[10]] <- lm(logZn ~ heatsubspdlpow100_3 + heatbing1902log300proj, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[10]])
#ols_plot_diagnostics(modlistlogZn[[10]])
ols_coll_diag(modlistlogZn[[10]])
ols_correlations(modlistlogZn[[10]])

modlistlogZn[[11]] <- lm(logZn ~  heatsubspdlpow100_3 + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[11]])
#ols_plot_diagnostics(modlistlogZn[[11]])
ols_coll_diag(modlistlogZn[[11]])
ols_correlations(modlistlogZn[[11]])

modlistlogZn[[12]] <- lm(logZn ~ heatsubspdlpow100_3 + heatbustransitlog200, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[12]])
#ols_plot_diagnostics(modlistlogZn[[12]])
ols_coll_diag(modlistlogZn[[12]])
ols_correlations(modlistlogZn[[12]])

modlistlogZn[[13]] <- lm(logZn ~ heatsubspdlpow100_3 + heatsubslopepow500_1, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[13]])
#ols_plot_diagnostics(modlistlogZn[[13]])
ols_coll_diag(modlistlogZn[[13]])
ols_correlations(modlistlogZn[[13]])

modlistlogZn[[14]] <- lm(logZn ~ heatsubspdlpow100_3*heatsubAADTlog100, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[14]])
#ols_plot_diagnostics(modlistlogZn[[14]])
ols_coll_diag(modlistlogZn[[14]])
ols_correlations(modlistlogZn[[14]])

modlistlogZn[[15]] <- lm(logZn ~ heatsubspdlpow100_3*heatbing1902log300proj, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[15]])
#ols_plot_diagnostics(modlistlogZn[[15]])
ols_coll_diag(modlistlogZn[[15]])
ols_correlations(modlistlogZn[[15]])

modlistlogZn[[16]] <- lm(logZn ~ heatsubspdlpow100_3*heatbing1902log300proj + nlcd_imp_ps_mean , 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[16]])
#ols_plot_diagnostics(modlistlogZn[[16]])
ols_coll_diag(modlistlogZn[[16]])
ols_correlations(modlistlogZn[[16]])

modlistlogZn[[17]] <- lm(logZn ~ heatsubspdlpow100_3 + heatbing1902log300proj + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[17]])
#ols_plot_diagnostics(modlistlogZn[[17]])
ols_coll_diag(modlistlogZn[[17]])
ols_correlations(modlistlogZn[[17]])

modlistlogZn[[18]] <- lm(logZn ~ heatsubspdlpow100_3*nlcd_imp_ps_mean + heatbing1902log300proj,
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[18]])
#ols_plot_diagnostics(modlistlogZn[[18]])
ols_coll_diag(modlistlogZn[18])
ols_correlations(modlistlogZn[[18]])

modlistlogZn[[19]] <- lm(logZn ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[19]])
#ols_plot_diagnostics(modlistlogZn[[19]])
ols_correlations(modlistlogZn[[19]])

modlistlogZn[[20]] <- lm(logZn ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean + 
                        heatbustransitlog200, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[20]])
#ols_plot_diagnostics(modlistlogZn[[20]])
ols_coll_diag(modlistlogZn[20])
ols_correlations(modlistlogZn[[20]])


#------ 8. logZn - Transformed predictors --------
modlistlogZn[[21]] <- lm(logZn ~ heatsubAADTlog200frt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[21]])
#ols_plot_diagnostics(modlistlogZn[[21]])
ols_coll_diag(modlistlogZn[[21]])
ols_correlations(modlistlogZn[[21]])

modlistlogZn[[22]] <- lm(logZn ~ heatsubspdlpow100_3frt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[22]])
#ols_plot_diagnostics(modlistlogZn[[22]])
ols_coll_diag(modlistlogZn[[22]])
ols_correlations(modlistlogZn[[22]])

modlistlogZn[[23]] <- lm(logZn ~ heatbing1902log300projfrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[23]])
#ols_plot_diagnostics(modlistlogZn[[23]])
ols_coll_diag(modlistlogZn[[23]])
ols_correlations(modlistlogZn[[23]])

modlistlogZn[[24]] <- lm(logZn ~ heatsubAADTlog300frt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[24]])
#ols_plot_diagnostics(modlistlogZn[[24]])
ols_coll_diag(modlistlogZn[[24]])
ols_correlations(modlistlogZn[[24]])

modlistlogZn[[25]] <- lm(logZn ~ heatsubAADTlog100frt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[25]])
#ols_plot_diagnostics(modlistlogZn[[25]])
ols_coll_diag(modlistlogZn[[25]])
ols_correlations(modlistlogZn[[25]])

modlistlogZn[[26]] <- lm(logZn ~ heatsubAADTlog100frt + heatsubspdlpow100_3sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[26]])
#ols_plot_diagnostics(modlistlogZn[[26]])
ols_coll_diag(modlistlogZn[[26]])
ols_correlations(modlistlogZn[[26]])

modlistlogZn[[27]] <- lm(logZn ~ heatsubAADTlog100frt*heatsubspdlpow100_3, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[27]])
#ols_plot_diagnostics(modlistlogZn[[27]])
ols_coll_diag(modlistlogZn[[27]])
ols_correlations(modlistlogZn[[27]])

modlistlogZn[[28]] <- lm(logZn ~ heatsubAADTlog100frt*heatsubspdlpow100_3 + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[28]])
#ols_plot_diagnostics(modlistlogZn[[28]])
ols_coll_diag(modlistlogZn[[28]])
ols_correlations(modlistlogZn[[28]])

modlistlogZn[[29]] <- lm(logZn ~ (heatsubAADTlog100frt)*heatsubspdlpow100_3 + heatbing1902log300proj, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[29]])
#ols_plot_diagnostics(modlistlogZn[[29]])
ols_coll_diag(modlistlogZn[[29]])
ols_correlations(modlistlogZn[[29]])

modlistlogZn[[30]] <- lm(logZn ~ (heatsubAADTlog100frt)*heatsubspdlpow100_3 + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubAADTlog100frt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[30]])
#ols_plot_diagnostics(modlistlogZn[[30]])
ols_coll_diag(modlistlogZn[[30]])
ols_correlations(modlistlogZn[[30]])


modlistlogZn[[31]] <- lm(logZn ~ (heatsubAADTlog100frt) + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubAADTlog100frt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[31]])
#ols_plot_diagnostics(modlistlogZn[[31]])
ols_coll_diag(modlistlogZn[[31]])
ols_correlations(modlistlogZn[[31]])

modlistlogZn[[32]] <- lm(logZn ~ (heatsubspdlpow100_1) + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubspdlpow100_3, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[32]])
#ols_plot_diagnostics(modlistlogZn[[32]])
ols_coll_diag(modlistlogZn[[32]])
ols_correlations(modlistlogZn[[32]])

modlistlogZn[[33]] <- lm(logZn ~ (heatsubAADTlog100frt) + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubAADTlog100frt + nlcd_imp_ps, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[33]])
#ols_plot_diagnostics(modlistlogZn[[33]])
ols_coll_diag(modlistlogZn[[33]])
ols_correlations(modlistlogZn[[33]])

modlistlogZn[[34]] <- lm(logZn ~ (heatsubAADTlog100frt) + heatbing1902log300proj + nlcd_imp_ps +
                        heatbustransitlog200, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[34]])
#ols_plot_diagnostics(modlistlogZn[[34]])
ols_coll_diag(modlistlogZn[[34]])
ols_correlations(modlistlogZn[[34]])

modlistlogZn[[35]] <- lm(logZn ~ (heatsubAADTlog100frt) + heatbing1902log300proj + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[35]])
#ols_plot_diagnostics(modlistlogZn[[35]])
ols_coll_diag(modlistlogZn[[35]])
ols_correlations(modlistlogZn[[35]])

modlistlogZn[[36]] <- lm(logZn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[36]])
#ols_plot_diagnostics(modlistlogZn[[36]])
ols_coll_diag(modlistlogZn[[36]])
ols_correlations(modlistlogZn[[36]])

modlistlogZn[[37]] <- lm(logZn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + heatbing1902log500proj,
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[37]])
#ols_plot_diagnostics(modlistlogZn[[37]])
ols_coll_diag(modlistlogZn[[37]])
ols_correlations(modlistlogZn[[37]])

modlistlogZn[[38]] <- lm(logZn ~ heatsubAADTlog200frt + nlcd_imp_ps_mean + heatbustransitlog200,
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[38]])
#ols_plot_diagnostics(modlistlogZn[[38]])
ols_coll_diag(modlistlogZn[[38]])
ols_correlations(modlistlogZn[[38]])

modlistlogZn[[39]] <- lm(logZn ~ heatsubAADTlog200frt + nlcd_imp_ps_mean + heatsubspdlpow300_3,
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZn[[39]])
#ols_plot_diagnostics(modlistlogZn[[39]])
ols_coll_diag(modlistlogZn[[39]])
ols_correlations(modlistlogZn[[39]])

#------ 9. logZn - Make latex model summary table ----
vnum <- max(sapply(modlistlogZn, function(mod) {length(mod$coefficients)}))
model_summarylogZn<- as.data.table(
  ldply(modlistlogZn, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summarylogZn, -R2pred, AICc)  
cat(latex_format(model_summarylogZn), file = file.path(moddir, 'modeltable_logZn_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_logZn_2019.tex')

#------ 10. logZn - Make latex model summary table when excluding outliers----
vnum <- max(sapply(modlistlogZn, function(mod) {length(mod$coefficients)}))
model_summarylogZn_nooutliers<- as.data.table(
  ldply(modlistlogZn, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50,
                                                             remove_outliers = 'outliers', # & leverage',
                                                             labelvec = pollutfieldclean_cast[, SiteIDPair])}))
setorder(model_summarylogZn_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summarylogZn_nooutliers), file = file.path(moddir, 'modeltable_logZn_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_logZn_nooutliers_2019.tex')

#------ 11. logZn - test selected model ----
summary(modlistlogZn[[36]])
ols_plot_diagnostics(modlistlogZn[[36]])
regdiagnostic_customtab(mod=modlistlogZn[[36]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistlogZn[[36]]))
MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistlogZn[[36]]))
RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistlogZn[[36]]))

qplot(pollutfieldclean_cast[, heatbustransitlog200],modlistlogZn[[36]]$residuals) + 
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[, heatbing1902log500proj], modlistlogZn[[36]]$residuals) + 
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[, NLCD_reclass_final_PS], modlistlogZn[[36]]$residuals) + 
  geom_smooth(method='lm', color='red')


summary(modlistlogZn[[38]])
GAMrescheck(modlistlogZn[[38]])
regdiagnostic_customtab(mod=modlistlogZn[[38]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistlogZn[[38]]))
MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistlogZn[[38]]))
RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistlogZn[[38]]))

qplot(abs(fitted(modlistlogZn[[36]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn]), 
      abs(fitted(modlistlogZn[[38]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model 36 - aadt10frt + nlcd', 
       y='Absolute error for model 38 - aadt100frt + nlcd + transit200')

#------ 12. GLM Zn - run all models for table ----
modlistglmZn <- list()

modlistglmZn[[1]] <- glm(Zn ~ 1,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[2]] <- glm(Zn ~ heatsubAADTlog200,
                        data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))
modlistglmZn[[3]] <- glm(Zn ~ heatbing1902log300proj,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))
modlistglmZn[[4]] <- glm(Zn ~ nlcd_imp_ps_mean,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[5]] <- glm(Zn ~ nlcd_imp_ps,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))
modlistglmZn[[6]] <- glm(Zn ~ heatbustransitlog200,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[7]] <- glm(Zn ~ heatsubspdlpow100_3,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[8]] <- glm(Zn ~ heatsubslopepow500_3,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[9]] <- glm(Zn ~ heatsubspdlpow100_3,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[10]] <- glm(Zn ~ heatsubspdlpow100_3 + heatsubAADTlog200,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[11]] <- glm(Zn ~ heatsubspdlpow100_3+ heatbing1902log300proj,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[12]] <- glm(Zn ~ heatsubspdlpow100_3 + nlcd_imp_ps_mean,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[13]] <- glm(Zn ~ heatsubspdlpow100_3 + heatbustransitlog200,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[14]] <- glm(Zn ~ heatsubspdlpow100_3 + heatsubslopepow500_1,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[15]] <- glm(Zn ~ heatsubspdlpow100_3*heatsubAADTlog200,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))


modlistglmZn[[16]] <- glm(Zn ~ heatsubspdlpow100_3*heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[17]] <- glm(Zn ~ heatsubspdlpow100_3*heatbing1902log300proj + nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[18]] <- glm(Zn ~ heatsubspdlpow100_3*nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[19]] <- glm(Zn ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[20]] <- glm(Zn ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean + heatbustransitlog200,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[21]] <- glm(Zn ~ heatsubspdlpow100_3*heatbustransitlog200 + heatbing1902log300proj*nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))


modlistglmZn[[22]] <- glm(Zn ~ heatsubspdlpow100_3*heatbustransitlog200 + heatbing1902log300proj*nlcd_imp_ps,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[23]] <- glm(Zn ~ heatsubAADTlog200frt,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[24]] <- glm(Zn ~ heatsubAADTlog100frt,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[25]] <- glm(Zn ~ heatsubAADTlog300frt,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))


modlistglmZn[[26]] <- glm(Zn ~ heatsubAADTpow100_1frt,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[27]] <- glm(Zn ~ heatsubAADTlog100frt + heatsubspdlpow100_3,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[28]] <- glm(Zn ~ heatsubAADTlog100frt*heatsubspdlpow100_3,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[29]] <- glm(Zn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[29]] <- glm(Zn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[30]] <- glm(Zn ~ heatsubAADTlog100frt*nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[31]] <- glm(Zn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[32]] <- glm(Zn ~ heatsubAADTlog100frt*nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[33]] <- glm(Zn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean*heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[34]] <- glm(Zn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + heatbing1902log300proj + heatbustransitlog200,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[35]] <- glm(Zn ~ heatsubAADTlog100frt*heatsubslopepow500_1 + nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[36]] <- glm(Zn ~ heatsubAADTlog100frt*heatbustransitlog200 + nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[37]] <- glm(Zn ~ heatsubAADTlog100frt*heatbustransitlog200 + heatsubspdllog100 + 
                            nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[38]] <- glm(Zn ~ heatsubAADTlog100frt + heatbustransitlog200 + nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[39]] <- glm(Zn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + heatbing1902log500proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZn[[40]] <- glm(Zn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean*heatbing1902log500proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

# modlistglmZn[[41]] <- glm(Zn ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + NLCD_reclass_final_PS_edit,
#                           data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))
# 
# modlistglmZn[[42]] <- glm(Zn ~ heatsubAADTlog100frt + heatbing1902log500proj + NLCD_reclass_final_PS_edit,
#                           data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))


#------ 13. GLM Zn - Make latex model summary table ----
vnum <- max(sapply(modlistglmZn, function(mod) {length(mod$coefficients)}))

model_summaryglmZn<- as.data.table(
  ldply(modlistglmZn, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
model_summaryglmZn[, AICc := as.numeric(AICc)]
model_summaryglmZn[nvars==0, `:=`(
  MAEcv = as.character(round(
    MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistglmZn[[1]])),2)),
  RMSEcv = as.character(round(
    RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistglmZn[[1]])),2))
  )]
setorder(model_summaryglmZn, AICc)  
cat(latex_format(model_summaryglmZn[, -c("BreuschPagan\\_fitp", "Score\\_fitp", "R2", "R2adj", "R2pred", "RMSE", "MAE", "VIF8")]),
    file = file.path(moddir, 'modeltable_glmZn_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_glmZn_2019.tex')

#Format table by hand
#pow100\_1 to linear (and other kernel size)
#spdl to speed
#nlcd\_imp\_ps to imperviousness
#heatbing1902 to congestion
#bustransit to transit
#remove sub
#reduce age size to 25 in width and 8 in height

texi2pdf('modeltable_glmZn_2019edit.tex')

#------ 14. GLM Zn - test selected model ----
summary(modlistglmZn[[29]])
GAMrescheck(modlistglmZn[[29]])
regdiagnostic_customtab(mod=modlistglmZn[[29]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistglmZn[[29]]))
MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistglmZn[[29]]))
RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistglmZn[[29]]))

qplot(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), heatbustransitlog200],
      modlistglmZn[[29]]$residuals) + 
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), heatbing1902log500proj],
      modlistglmZn[[29]]$residuals) + 
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), NLCD_reclass_final_PS],
      modlistglmZn[[29]]$residuals) + 
  geom_smooth(method='lm', color='red')


summary(modlistglmZn[[38]])
GAMrescheck(modlistglmZn[[38]])
regdiagnostic_customtab(mod=modlistglmZn[[38]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistglmZn[[38]]))
MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistglmZn[[38]]))
RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn], fitted(modlistglmZn[[38]]))

qplot(abs(fitted(modlistglmZn[[29]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn]), 
      abs(fitted(modlistglmZn[[38]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model 29 - aadt10frt + nlcd', 
       y='Absolute error for model 38 - aadt100frt + nlcd + transit200')


GAMrescheck(modlistglmZn[[39]])
qplot(abs(fitted(modlistglmZn[[29]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn]), 
      abs(fitted(modlistglmZn[[39]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Zn])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model 29 - aadt100frt + nlcd', 
       y='Absolute error for model 39 - aadt100frt + nlcd + bing500')

# plot_model(modlistglmZn[[33]], type='int', mdrt.values = 'meansd') +
#   scale_x_continuous(name = 'Imperviousness (%)', expand = c(0,0)) +
#   scale_y_continuous(name = 'Zn index', expand=c(0,0), 
#                      breaks = c(0,0.25, 0.5, 0.75, 1.0), labels = c(0, 25, 50, 75, 100)) +
#   coord_cartesian(xlim=c(0,100), ylim=c(0,1)) +
#   theme_classic() +
#   theme(plot.title=element_blank(),
#         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
#         legend.position = c(0.8, 0.2), 
#         text =element_text(size=14))


#plot_model(modlistglmZn[[33]], type='pred', mdrt.values = 'meansd')

#------ 15. GAM Zn - Single parameter models -------
modlistGAMZn <- list() #List to hold models
modlistGAMZn[[1]] <- lm(logZn ~ 1, data = pollutfieldclean_cast) #Null/Intercept model

modlistGAMZn[[2]] <- mgcv::gam(Zn~s(heatsubAADTlog100, k=4), family=gaussian(link='log'), method='ML',
                               data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[2]])
plot(modlistGAMZn[[2]],residuals=TRUE,shade=T, cex=6)

pollutfieldclean_cast$predgamZn <- predict(modlistGAMZn[[2]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog200, y=log(Zn), label=paste0(SiteID, Pair))) +
  geom_text() + 
  geom_line(aes(y=predgamZn))

modlistGAMZn[[3]] <- mgcv::gam(Zn~s(heatsubAADTlog100, k=4), family=Gamma(link='log'), method='ML',
                              data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[3]])
plot(modlistGAMZn[[3]],residuals=TRUE,shade=T, cex=6)

pollutfieldclean_cast$predgamZn <- predict(modlistGAMZn[[3]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog200, y=log(Zn), label=paste0(SiteID, Pair))) +
  geom_text() + 
  geom_line(aes(y=predgamZn), color='red', size=1.2) +
  theme_classic()

modlistGAMZn[[4]] <- mgcv::gam(Zn~s(heatsubAADTlog200, k=4), family=gaussian(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[4]])
plot(modlistGAMZn[[4]],residuals=TRUE,shade=T, cex=6)

modlistGAMZn[[5]] <- mgcv::gam(Zn~s(heatsubAADTlog200, k=4), family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[5]])
plot(modlistGAMZn[[5]],residuals=TRUE,shade=T, cex=6)

modlistGAMZn[[6]] <- mgcv::gam(Zn~s(heatsubspdlpow100_3, k=4), family=gaussian(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[6]])
plot(modlistGAMZn[[6]],residuals=TRUE,shade=T, cex=6)

modlistGAMZn[[7]] <- mgcv::gam(Zn~s(heatbing1902log300proj, k=4), 
                               family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[7]])
plot(modlistGAMZn[[7]],residuals=TRUE,shade=T, cex=6)

modlistGAMZn[[8]] <- mgcv::gam(Zn~s(heatbustransitlog200, k=4), 
                               family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[8]])
plot(modlistGAMZn[[8]],residuals=TRUE,shade=T, cex=6)

modlistGAMZn[[9]] <- mgcv::gam(Zn~s(nlcd_imp_ps, k=3), 
                               family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[9]])
plot(modlistGAMZn[[9]],residuals=TRUE,shade=T, cex=6)

modlistGAMZn[[10]] <- mgcv::gam(Zn~s(nlcd_imp_ps_mean, k=3), 
                               family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[10]])
plot(modlistGAMZn[[10]],residuals=TRUE,shade=T, cex=6)

modlistGAMZn[[11]] <- mgcv::gam(Zn~s(heatsubslopepow500_1, k=4), 
                               family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[11]])
plot(modlistGAMZn[[11]],residuals=TRUE,shade=T, cex=6)

modlistGAMZn[[12]] <- mgcv::gam(Zn~s(heatbing1902log500proj, k=4), 
                               family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[12]])
plot(modlistGAMZn[[12]],residuals=TRUE,shade=T, cex=6)

#------ 16. GAM Zn - Multiparameter models -------
modlistGAMZn[[11]] <- mgcv::gam(Zn~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3),
                                family=Gamma(link='log'), 
                                data=pollutfieldclean_cast,
                                select=TRUE)
GAMrescheck(modlistGAMZn[[11]])
GAMmultiplot(modlistGAMZn[[11]])
concurvity(modlistGAMZn[[11]])

pollutfieldclean_cast[, predgamZn := fitted(modlistGAMZn[[11]], type='response')]
ggplot(pollutfieldclean_cast, aes(x=predgamZn, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text() +
  geom_abline(intercept=0, slope=1) +
  #scale_x_continuous(limits=c(0,2)) +
  coord_fixed()

modlistGAMZn[[12]] <- mgcv::gam(Zn~s(heatsubAADTlog100thd) + s(heatsubAADTlog100thd, k=4, by=nlcd_imp_ps_mean) + 
                                  s(nlcd_imp_ps, k=3), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[12]])
GAMmultiplot(modlistGAMZn[[12]])


modlistGAMZn[[13]] <- mgcv::gam(Zn~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbing1902log300proj, k=4), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[13]])
GAMmultiplot(modlistGAMZn[[13]])
concurvity(modlistGAMZn[[13]])

modlistGAMZn[[14]] <- mgcv::gam(Zn~s(heatsubAADTlog100,k=4) + s(heatsubAADTlog100, k=4, by=heatbing1902log500proj) + 
                                  s(nlcd_imp_ps, k=3), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[14]])
GAMmultiplot(modlistGAMZn[[14]])

modlistGAMZn[[15]] <- mgcv::gam(Zn~s(heatsubAADTlog100thd,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[15]])
GAMmultiplot(modlistGAMZn[[15]])
#Plot
#Plot

modlistGAMZn[[16]] <- mgcv::gam(Zn~s(heatsubAADTlog100, k=4) + s(heatsubspdlpow100_3, k=4)+ 
                                  s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[16]])
GAMmultiplot(modlistGAMZn[[16]])

modlistGAMZn[[17]] <- mgcv::gam(Zn~s(heatsubAADTlog100, k=4) + s(heatsubAADTlog100, k=4, by=heatsubspdlpow100_3)+ 
                                  s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[17]])
GAMmultiplot(modlistGAMZn[[17]])
concurvity(modlistGAMZn[[17]])

pollutfieldclean_cast[, predgamZn := fitted(modlistGAMZn[[17]], type='response')]
ggplot(pollutfieldclean_cast, aes(x=predgamZn, y=Zn, label=paste0(SiteID, Pair))) + 
  geom_text() +
  geom_abline(intercept=0, slope=1) +
  scale_x_continuous(limits=c(0,2)) +
  coord_fixed()

modlistGAMZn[[18]] <- mgcv::gam(Zn~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200) +
                                  s(heatsubAADTlog100, k=4, by=heatsubslopepow500_1), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[18]])
GAMmultiplot(modlistGAMZn[[18]])
concurvity(modlistGAMZn[[18]])

modlistGAMZn[[19]] <- mgcv::gam(Zn~s(heatsubAADTlog100thd,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200, k=4) +
                                  s(heatbing1902log500proj, k=3), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[19]])
GAMmultiplot(modlistGAMZn[[19]])
concurvity(modlistGAMZn[[19]])

modlistGAMZn[[20]] <- mgcv::gam(Zn~s(heatsubAADTlog100thd,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog300, k=4) +
                                  s(heatbing1902log500proj, k=3), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast[])
GAMrescheck(modlistGAMZn[[20]])
GAMmultiplot(modlistGAMZn[[20]])
concurvity(modlistGAMZn[[20]])

modlistGAMZn[[21]] <- mgcv::gam(Zn~s(heatsubspdllog100,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZn[[21]])
GAMmultiplot(modlistGAMZn[[21]])
concurvity(modlistGAMZn[[21]])

#------ 17. GAM Zn - Multiparameter models without Luwam data or 15th Ave outliers -------
GAMZn_sub11 <- mgcv::gam(Zn~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3), 
                         family=Gamma(link='log'), data=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),])
GAMrescheck(GAMZn_sub11)
GAMmultiplot(GAMZn_sub11)
concurvity(GAMZn_sub11)

GAMZn_sub13 <- mgcv::gam(Zn~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbing1902log300proj, k=4), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),])
GAMrescheck(GAMZn_sub13)
GAMmultiplot(GAMZn_sub13)
concurvity(GAMZn_sub13)

GAMZn_sub20 <- mgcv::gam(Zn~s(heatsubAADTlog100,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                         family=Gamma(link='log'),data=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),])
GAMrescheck(GAMZn_sub20)
GAMmultiplot(GAMZn_sub20)

#------ 18. GAM Zn - Multiparameter models without outliers based on influence-------
#Model 11
influencemod11 <- GAMinfluence(modlistGAMZn[[11]], pollutfieldclean_cast)

ggplot(pollutfieldclean_cast, aes(x=fitted(modlistGAMZn[[11]]), y=Zn, 
                                  color=influencemod11$influence, label=SiteIDPair)) + 
  geom_text() +
  scale_color_distiller(palette='Spectral')
ggplot(influencemod11, aes(x=influence, label=SiteIDPair)) + 
  geom_histogram()
pollutfieldclean_cast[!(SiteIDPair %in% 
                          influencemod11[influencemod11$influence>0.04, 'SiteIDPair'])]

GAMZn_sub11 <- mgcv::gam(Zn~s(heatsubAADTlog100,k=4) + s(nlcd_imp_ps_mean, k=3), 
                         family=Gamma(link='log'), 
                         data=pollutfieldclean_cast[!(SiteIDPair %in% 
                                                        influencemod11[influencemod11$influence>0.04, 'SiteIDPair']),])
GAMrescheck(GAMZn_sub11)
GAMmultiplot(GAMZn_sub11)
GAMrescheck(modlistGAMZn[[11]])

#-------19. Zn - Compare final selected models -------
mod29_nooutliers <- regdiagnostic_customtab(modlistglmZn[[29]], maxpar=vnum, 
                                            remove_outliers = 'outliers',
                                            labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                            kCV = TRUE, k=10, cvreps=50)
#Determined by deleted studentized residuals
subdat29 <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                                    strsplit(gsub('\\\\', '', mod29_nooutliers['outliers']), ',')$outliers),]
#Simply excluding 15th Ave bus stop and Luwam's data
subdat29_2 <- pollutfieldclean_cast[!(paste0(SiteIDPair) %in% extraoutliers),]
  
mod29_nooutliersub <- glm(Zn ~  heatsubAADTlog100frt + nlcd_imp_ps_mean, 
                         data = subdat29_2, family=Gamma('log'))
GAMrescheck(mod29_nooutliersub)
plot(mod29_nooutliersub)
#Plot
#Plot
#Plot
#Plot
MAPE(subdat29_2[, Zn], fitted(mod29_nooutliersub))
MAE(subdat29_2[, Zn], fitted(mod29_nooutliersub))
RMSE(subdat29_2[, Zn], fitted(mod29_nooutliersub))

AICc(mod29_nooutliersub)
qplot(subdat29_2$heatbustransitpow500_1, mod29_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat29_2$heatsubslopepow500_1, mod29_nooutliersub$residuals) +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat29_2$heatsubspdlpow100_1, mod29_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat29_2$heatsubAADTlog300, mod29_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')

#Plot of model29
mod29pred <- as.data.frame(predict(mod29_nooutliersub, 
                                   pollutfieldclean_cast[!(SiteIDPair %in% c('64A', '65A', '66A', '67A', '68A', '69A'))],
                                   type='response', se.fit=T))
Znmod29plot <- ggplot(data=pollutfieldclean_cast[!(SiteIDPair %in% c('64A', '65A', '66A', '67A', '68A', '69A')),],
                      aes(x=mod29pred$fit,  y=Zn)) +
  geom_point(alpha=0.75, size=2, color='red') + 
  geom_ribbon(aes(ymin=mod29pred$fit + 1.96*mod29pred$se.fit, ymax=mod29pred$fit - 1.96*mod29pred$se.fit), fill='orange', alpha=1/4) + 
  geom_point(data=subdat29_2, aes(x=fitted(mod29_nooutliersub), y=Zn), color='black', alpha=1/2, size=2) +
  geom_abline(intercept=0, slope=1, size=1.3, color='red') +
  #geom_smooth(method='lm', se=FALSE) +
  #geom_text(aes(label=paste0(SiteID,Pair))) + 
  scale_x_continuous(limits=c(0,2.1), expand=c(0,0), breaks=seq(0,2,0.5)) +
  scale_y_continuous(limits=c(0,2.1), expand=c(0,0), breaks=seq(0,2,0.5)) +
  coord_fixed() +
  labs(x='Predicted Zn index', y='Observed Zn index') +
  theme_classic() + 
  theme(text= element_text(size=20))
png(file.path(moddir, 'scatterplot_Zn_mod29.png'), width=9, height=9, units='in', res=300)
Znmod29plot
dev.off()

#Compare glm 29 with logZn 36
#Compare with predictions without transformation

qplot(abs(fitted(modlistglmZn[[29]])-pollutfieldclean_cast[!(paste0(SiteIDPair) %in% extraoutliers), Zn]), 
      abs(exp(predict(modlistlogZn[[36]], newdata=pollutfieldclean_cast[!(paste0(SiteIDPair) %in% extraoutliers)]))-
            pollutfieldclean_cast[!(paste0(SiteIDPair) %in% extraoutliers), Zn])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model glm 29 - glm  nlcd_imp_ps_mean + AADTlog100frt ', 
       y='Absolute error for model log 36 - nlcd_imp_ps_mean + AADTlog100frt')

#------ 20. Zn - Check spatial and temporal autocorrelation of residuals for full dataset -------
#Other good resource: https://eburchfield.github.io/files/Spatial_regression_LAB.html
resnorm <- rstandard(modlistlogZn[[36]]) #Get standardized residuals from model
#Make bubble map of residuals
bubbledat <- data.frame(resnorm, 
                        pollutfieldclean_cast$coords.x1, 
                        pollutfieldclean_cast$coords.x2)
coordinates(bubbledat) <- c("pollutfieldclean_cast.coords.x1","pollutfieldclean_cast.coords.x2")
bubble(bubbledat, "resnorm", col = c("blue","red"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")
#Check semi-variogram of residuals
plot(variogram(resnorm~1, bubbledat, cutoff=2000, width=100)) #isotropic
plot(variogram(resnorm~1, bubbledat, cutoff= 2000, width=100, alpha = c(0, 45, 90,135))) #anisotropic
#Check spline correlogram ()
plot(spline.correlog(x=coordinates(bubbledat)[,1], y=coordinates(bubbledat)[,2],
                     z=bubbledat$resnorm, resamp=500, quiet=TRUE, xmax = 5000))
#Compute a spatial weight matrix based on IDW
weightmat_k <- lapply(1:10, function(i) {
  weightmat_IDW(pollutfieldclean_cast[,.(coords.x1, coords.x2)], knb = i, mindist = 10)}) #Based on 10 nearest neighbors
weightmat_all <- weightmat_IDW(pollutfieldclean_cast[,.(coords.x1, coords.x2)], knb = NULL, mindist = 10) #Based on all points

#Moran plots
#lag_resnorm <- lag.listw(weightmat_all, resnorm) #Can be used to create customized Moran plot by plotting residuals against matrix
moran.plot(resnorm, weightmat_all, labels=pollutfieldclean_cast[, SiteIDPair], pch=19)
moran.plot(resnorm, weightmat_k[[2]], labels=pollutfieldclean_cast[, SiteIDPair], pch=19)

#Compute Moran's I
"Should always only use lm.morantest for residuals from regression, see http://r-sig-geo.2731867.n2.nabble.com/Differences-between-moran-test-and-lm-morantest-td7591336.html
for an explanation"
lm.morantest(modlistlogZn[[36]], listw = listw2U(weightmat_k[[1]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZn[[36]], listw = listw2U(weightmat_k[[2]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZn[[36]], listw = listw2U(weightmat_k[[3]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZn[[36]], listw = listw2U(weightmat_k[[4]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZn[[36]], listw = listw2U(weightmat_k[[5]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZn[[36]], listw = listw2U(weightmat_all)) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)

#Test for need for spatial regression model using Lagrange Multiplier (LM) tests
lm.LMtests(modlistlogZn[[36]], listw = listw2U(weightmat_k[[1]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(modlistlogZn[[36]], listw = listw2U(weightmat_k[[2]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(modlistlogZn[[36]], listw = listw2U(weightmat_k[[3]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))

#Spatial simultaneous autoregressive error model estimation with 1 nearest neighbors
sarlm_modlogZn <- errorsarlm(modlistlogZn[[36]]$call$formula, data = modlistlogZn[[36]]$model, 
                          listw = listw2U(weightmat_k[[1]]))
summary(sarlm_modlogZn)
bptest.sarlm(sarlm_modlogZn)

#Compare pseudo-R2
cor(modlistlogZn[[36]]$model$logZn, fitted(sarlm_modlogZn))^2
cor(modlistlogZn[[36]]$model$logZn, fitted(modlistlogZn[[36]]))^2
#Compare MAE
DescTools::MAE(modlistlogZn[[36]]$model$logZn, fitted(sarlm_modlogZn))
DescTools::MAE(modlistlogZn[[36]]$model$logZn, fitted(modlistlogZn[[36]]))

#Compare observed~predicted for full-no outlier model and for aspatial and spatial model
spatial_comparisonplot <- ggplot(pollutfieldclean_cast, aes(x=exp(fitted(sarlm_modlogZn)), y=exp(logZn))) + 
  geom_point(aes(x=exp(fitted(modlistlogZn[[36]]))), size=2, alpha=1/2, color='orange') +
  geom_point(size=2, alpha=1/2, color='red') + 
  geom_point(aes(x=exp(-1.8256154+0.4308591*heatsubAADTlog100frt+0.0112312*nlcd_imp_ps_mean)),
             size=2,color='darkgreen', alpha=1/2) +
  geom_abline(size=1, slope=1, intercept=0, color='red') + 
  #geom_text(aes(label=paste0(SiteID, Pair))) +
  coord_fixed() +
  theme_classic()
spatial_comparisonplot

resnorm_postsarlm <- residuals(sarlm_modlogZn) #Get standardized residuals from model
#Make bubble map of residuals
bubbledat_postsarlm <- data.frame(resnorm_postsarlm, pollutfieldclean_cast$coords.x1, pollutfieldclean_cast$coords.x2)
coordinates(bubbledat_postsarlm) <- c("pollutfieldclean_cast.coords.x1","pollutfieldclean_cast.coords.x2")
bubble(bubbledat_postsarlm, "resnorm_postsarlm", col = c("blue","red"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")

### FINAL VERDICT: GO WITH SARLM COEFFICIENTS FOR MODLOG36-----

#--------------- B. Cu ~ separate predictors ----
#Exclude 33A and 33B as way too much of an outlier site
modlistCu <- list() #List to hold models
subdatCu <- pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),]
modlistCu[[1]] <- lm(Cu ~ 1, data = subdatCu) #Null/Intercept model
#------ 1. Cu - Single parameter models --------
modlistCu[[2]] <- lm(Cu ~ heatsubAADTlog200, data = subdatCu)
ols_regress(modlistCu[[2]])
#ols_plot_diagnostics(modlistCu[[2]])
ggplot(subdatCu, aes(x=heatsubAADTlog200, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[3]] <- lm(Cu ~ heatbing1902log300proj, data = subdatCu)
ols_regress(modlistCu[[3]])
#ols_plot_diagnostics(modlistCu[[3]])
ggplot(subdatCu, aes(x=heatbing1902log300proj, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[4]] <- lm(Cu ~ heatbing1902log500proj, data = subdatCu)
ols_regress(modlistCu[[4]])
#ols_plot_diagnostics(modlistCu[[4]])
ggplot(subdatCu, aes(x=heatbing1902log500proj, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[5]] <- lm(Cu ~ nlcd_imp_ps, data = subdatCu)
ols_regress(modlistCu[[5]])
#ols_plot_diagnostics(modlistCu[[5]])
ggplot(subdatCu, aes(x=nlcd_imp_ps, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[6]] <- lm(Cu ~ nlcd_imp_ps_mean, data = subdatCu)
ols_regress(modlistCu[[6]])
#ols_plot_diagnostics(modlistCu[[6]])
ggplot(subdatCu, aes(x=nlcd_imp_ps_mean, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[7]] <- lm(Cu ~ heatbustransitlog200, data = subdatCu)
ols_regress(modlistCu[[7]])
#ols_plot_diagnostics(modlistCu[[7]])
ggplot(subdatCu, aes(x=heatbustransitlog200, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[8]] <- lm(Cu ~ heatbustransitlog200sqrt, data = subdatCu)
ols_regress(modlistCu[[8]])
#ols_plot_diagnostics(modlistCu[[8]])
ggplot(subdatCu, aes(x=heatbustransitlog200sqrt, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[9]] <- lm(Cu ~ heatsubspdlpow50_1, data = subdatCu)
ols_regress(modlistCu[[9]])
#ols_plot_diagnostics(modlistCu[[9]])
ggplot(subdatCu, aes(x=heatsubspdlpow50_1, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[10]] <- lm(Cu ~ heatsubslopepow500_1, data = subdatCu)
ols_regress(modlistCu[[10]])
#ols_plot_diagnostics(modlistCu[[10]])
ggplot(subdatCu, aes(x=heatsubspdlpow500_1, y=Cu, label=SiteIDPair)) + 
  geom_text()


modlistCu[[11]] <- lm(Cu ~ heatsubslopepow500_1, data = subdatCu)
ols_regress(modlistCu[[11]])
#ols_plot_diagnostics(modlistCu[[11]])
ggplot(subdatCu, aes(x=heatsubspdlpow500_1, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[12]] <- lm(Cu ~ heatsubAADTlog100, data = subdatCu)
ols_regress(modlistCu[[12]])
#ols_plot_diagnostics(modlistCu[[12]])
ggplot(subdatCu, aes(x=heatsubAADTlog100, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[13]] <- lm(Cu ~ heatbustransitlog300sqrt, data = subdatCu)
ols_regress(modlistCu[[13]])
#ols_plot_diagnostics(modlistCu[[13]])
ggplot(subdatCu, aes(x=heatbustransitlog300sqrt, y=Cu, label=SiteIDPair)) + 
  geom_text()

modlistCu[[14]] <- lm(Cu ~ heatsubslopepow300_1, data = subdatCu)
ols_regress(modlistCu[[14]])
#ols_plot_diagnostics(modlistCu[[14]])
ggplot(subdatCu, aes(x=heatsubspdlpow300_1, y=Cu, label=SiteIDPair)) + 
  geom_text()


modlistCu[[15]] <- lm(Cu ~ heatsubslopepow300_1, data = subdatCu)
ols_regress(modlistCu[[15]])
#ols_plot_diagnostics(modlistCu[[15]])
ggplot(subdatCu, aes(x=heatsubspdlpow300_1, y=Cu, label=SiteIDPair)) + 
  geom_text()

#------ 2. Cu - Multiparameter models --------
modlistCu[[16]] <- lm(Cu ~ heatbustransitlog300sqrt + heatbing1902log500proj, data = subdatCu)
ols_regress(modlistCu[[16]])
#ols_plot_diagnostics(modlistCu[[16]])
ols_coll_diag(modlistCu[[16]])
ols_correlations(modlistCu[[16]])

modlistCu[[17]] <- lm(Cu ~ heatbustransitlog200sqrt + heatbing1902log500proj, data = subdatCu)
ols_regress(modlistCu[[17]])
#ols_plot_diagnostics(modlistCu[[17]])
ols_coll_diag(modlistCu[[17]])
ols_correlations(modlistCu[[17]])

modlistCu[[18]] <- lm(Cu ~ heatbustransitlog200sqrt + heatbing1902log300proj, data = subdatCu)
ols_regress(modlistCu[[18]])
#ols_plot_diagnostics(modlistCu[[18]])
ols_coll_diag(modlistCu[[18]])
ols_correlations(modlistCu[[18]])

modlistCu[[19]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean, data = subdatCu)
ols_regress(modlistCu[[19]])
#ols_plot_diagnostics(modlistCu[[19]])
ols_coll_diag(modlistCu[[19]])
ols_correlations(modlistCu[[19]])

modlistCu[[20]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatbing1902log300proj, data = subdatCu)
ols_regress(modlistCu[[20]])
#ols_plot_diagnostics(modlistCu[[20]])
ols_coll_diag(modlistCu[[20]])
ols_correlations(modlistCu[[20]])

modlistCu[[21]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubslopepow300_1, data = subdatCu)
ols_regress(modlistCu[[21]])
#ols_plot_diagnostics(modlistCu[[21]])
ols_coll_diag(modlistCu[[21]])
ols_correlations(modlistCu[[21]])

modlistCu[[22]] <- lm(Cu ~ heatsubspdlpow50_1 + nlcd_imp_ps_mean, data = subdatCu)
ols_regress(modlistCu[[22]])
#ols_plot_diagnostics(modlistCu[[22]])
ols_coll_diag(modlistCu[[22]])
ols_correlations(modlistCu[[22]])

modlistCu[[23]] <- lm(Cu ~ heatsubspdlpow50_1 + heatbing1902log300proj, data = subdatCu)
ols_regress(modlistCu[[23]])
#ols_plot_diagnostics(modlistCu[[23]])
ols_coll_diag(modlistCu[[23]])
ols_correlations(modlistCu[[23]])

modlistCu[[24]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubspdlpow50_1, data = subdatCu)
ols_regress(modlistCu[[24]])
#ols_plot_diagnostics(modlistCu[[24]])
ols_coll_diag(modlistCu[[24]])
ols_correlations(modlistCu[[24]])

modlistCu[[25]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubspdlpow50_1, data = subdatCu)
ols_regress(modlistCu[[25]])
#ols_plot_diagnostics(modlistCu[[25]])
ols_coll_diag(modlistCu[[25]])
ols_correlations(modlistCu[[25]])

modlistCu[[26]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatbing1902log500proj, data = subdatCu)
ols_regress(modlistCu[[26]])
#ols_plot_diagnostics(modlistCu[[26]])
ols_coll_diag(modlistCu[[26]])
ols_correlations(modlistCu[[26]])

modlistCu[[27]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean*heatbing1902log300proj, data = subdatCu)
ols_regress(modlistCu[[27]])
#ols_plot_diagnostics(modlistCu[[27]])
ols_coll_diag(modlistCu[[27]])
ols_correlations(modlistCu[[27]])

modlistCu[[28]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog100, data = subdatCu)
ols_regress(modlistCu[[28]])
#ols_plot_diagnostics(modlistCu[[28]])
ols_coll_diag(modlistCu[[28]])
ols_correlations(modlistCu[[28]])

modlistCu[[29]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog200, data = subdatCu)
ols_regress(modlistCu[[29]])
#ols_plot_diagnostics(modlistCu[[29]])
ols_coll_diag(modlistCu[[29]])
ols_correlations(modlistCu[[29]])

modlistCu[[30]] <- lm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog500, data = subdatCu)
ols_regress(modlistCu[[30]])
#ols_plot_diagnostics(modlistCu[[30]])
ols_coll_diag(modlistCu[[30]])
ols_correlations(modlistCu[[30]])

#------ 3. Cu - Make latex model summary table ----
vnum <- max(sapply(modlistCu, function(mod) {length(mod$coefficients)}))
model_summaryCu<- as.data.table(
  ldply(modlistCu, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50,
                                                          labelvec = subdatCu[, SiteIDPair])}))
model_summaryCu[, AICc := as.numeric(AICc)]
setorder(model_summaryCu, AICc, -R2pred)  
cat(latex_format(model_summaryCu), file = file.path(moddir, 'modeltable_Cu.tex'))
setwd(moddir)
texi2pdf('modeltable_Cu.tex')

#------ 4. Cu - Make latex model summary table when excluding outliers ----
vnum <- max(sapply(modlistCu, function(mod) {length(mod$coefficients)}))
model_summaryCu_nooutliers<- as.data.table(
  ldply(modlistCu, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50,
                                                             remove_outliers = 'outliers', # & leverage',
                                                             labelvec = pollutfieldclean_cast[, SiteIDPair])}))
setorder(model_summaryCu_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summaryCu_nooutliers), file = file.path(moddir, 'modeltable_Cu_nooutliers.tex'))
setwd(moddir)
texi2pdf('modeltable_Cu_nooutliers.tex')

#-------5. Cu - Check final selected model -------
mod19_nooutliers <- regdiagnostic_customtab(modlistCu[[19]], maxpar=vnum, 
                                           remove_outliers = 'outliers',
                                           labelvec = subdatCu[, paste0(SiteID, Pair)],
                                           kCV = TRUE, k=10, cvreps=50)
subdat <- subdatCu[!(paste0(SiteID, Pair) %in% 
                       strsplit(gsub('\\\\', '', mod19_nooutliers['outliers']), ',')$outliers),]
mod19_nooutliersub <- lm(Cu ~ heatbustransitlog200sqrt +  nlcd_imp_ps_mean, 
                        data = subdat)
summary(mod19_nooutliersub)
ols_plot_diagnostics(mod19_nooutliersub)
ols_regress(mod19_nooutliersub)
par(mfrow=c(2,2))
plot(mod19_nooutliersub)
ols_coll_diag(mod19_nooutliersub)
ols_correlations(mod19_nooutliersub)
AICc(mod19_nooutliersub)

qplot(subdat$heatbing1902log300proj, mod19_nooutliersub$residuals) +
  geom_smooth(method='lm', color='red')
qplot(subdat$heatsubspdllog300, mod19_nooutliersub$residuals) +
  geom_smooth(method='lm', color='red')
qplot(subdat$heatsubAADTlog300, mod19_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat$heatsubslopelog500, mod19_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')

ggplot(subdat, aes(x=NLCD_reclass_final_PS, y=mod19_nooutliersub$residuals)) +
  geom_boxplot() +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')

ggplot(subdatCu, aes(y=predict(modlistCu[[19]], subdatCu), x=Cu)) +
  geom_point(alpha=1/4, size=2) + 
  #geom_point(data=subdat, aes(y = predict(modlistCu[[19]], subdatCu)), alpha=1/2, size=2) + 
  geom_point(data=subdat, aes(x=Cu, y = predict(mod19_nooutliersub, subdat)), color='red', alpha=1/2, size=2) + 
  geom_abline(intercept=0, slope=1) +
  theme_classic()


#------ 6. GLM Cu - run all models for table ----
modlistglmCu <- list()

modlistglmCu[[1]] <- glm(Cu ~ heatsubAADTlog200, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[2]] <- glm(Cu ~ heatsubAADTlog200, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[2]] <- glm(Cu ~ heatsubAADTlog200, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[3]] <- glm(Cu ~ heatbing1902log300proj, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[4]] <- glm(Cu ~ heatbing1902log500proj, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[5]] <- glm(Cu ~ nlcd_imp_ps, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[6]] <- glm(Cu ~ nlcd_imp_ps_mean, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[7]] <- glm(Cu ~ heatbustransitlog200, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[8]] <- glm(Cu ~ heatbustransitlog200sqrt, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[9]] <- glm(Cu ~ heatsubspdlpow50_1, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[10]] <- glm(Cu ~ heatsubslopepow500_1, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[11]] <- glm(Cu ~ heatsubslopepow500_1, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[12]] <- glm(Cu ~ heatsubAADTlog100, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[13]] <- glm(Cu ~ heatbustransitlog300sqrt, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[14]] <- glm(Cu ~ heatsubslopepow300_1, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[15]] <- glm(Cu ~ heatsubslopepow300_1, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[16]] <- glm(Cu ~ heatbustransitlog300sqrt + heatbing1902log500proj, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[17]] <- glm(Cu ~ heatbustransitlog200sqrt + heatbing1902log500proj, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[18]] <- glm(Cu ~ heatbustransitlog200sqrt + heatbing1902log300proj, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[19]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[20]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatbing1902log300proj, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[21]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubslopepow300_1, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[22]] <- glm(Cu ~ heatsubspdlpow50_1 + nlcd_imp_ps_mean, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[23]] <- glm(Cu ~ heatsubspdlpow50_1 + heatbing1902log300proj, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[24]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubspdlpow50_1, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[25]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubspdlpow50_1, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[26]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatbing1902log500proj, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[27]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean*heatbing1902log300proj, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[28]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog100, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[29]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog200, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[30]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog300, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[31]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog500, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[32]] <- glm(Cu ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean*heatsubAADTlog500, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[33]] <- glm(Cu ~ heatbustransitlog200sqrt + heatsubAADTlog300, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[34]] <- glm(Cu ~ heatbustransitlog200sqrt + heatsubAADTlog500, data=subdatCu, family=Gamma(link='log'))

modlistglmCu[[35]] <- glm(Cu ~ nlcd_imp_ps_mean + heatsubAADTlog300, data=subdatCu, family=Gamma(link='log'))

#------ 7. GLM Cu - Make latex model summary table ----
vnum <- max(sapply(modlistglmCu, function(mod) {length(mod$coefficients)}))

model_summaryglmCu<- as.data.table(
  ldply(modlistglmCu, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
model_summaryglmCu[, AICc := as.numeric(AICc)]
model_summaryglmCu[nvars==0, `:=`(
  MAEcv = as.character(round(
    MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Cu], fitted(modlistglmCu[[1]])),2)),
  RMSEcv = as.character(round(
    RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Cu], fitted(modlistglmCu[[1]])),2))
)]
setorder(model_summaryglmCu, AICc)  
cat(latex_format(model_summaryglmCu[, -c("BreuschPagan\\_fitp", "Score\\_fitp", "R2", "R2adj", "R2pred", "RMSE", "MAE", "VIF8")]),
    file = file.path(moddir, 'modeltable_glmCu_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_glmCu_2019.tex')

#Format table by hand
#pow100\_1 to linear (and other kernel size)
#spdl to speed
#nlcd\_imp\_ps to imperviousness
#heatbing1902 to congestion
#bustransit to transit
#remove sub
#reduce age size to 25 in width and 8 in height

texi2pdf('modeltable_glmCu_2019edit.tex')

#------ 8. GLM Cu - Make latex model summary table when excluding outliers ----
vnum <- max(sapply(modlistglmCu, function(mod) {length(mod$coefficients)}))
model_summaryglmCu_nooutliers<- as.data.table(
  ldply(modlistglmCu, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50,
                                                          remove_outliers = 'outliers', # & leverage',
                                                          labelvec = subdatCu[, SiteIDPair])}))
setorder(model_summaryglmCu_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summaryglmCu_nooutliers), file = file.path(moddir, 'modeltable_glmCu_nooutliers.tex'))
setwd(moddir)
texi2pdf('modeltable_glmCu_nooutliers.tex')

#------ 9. GLM Cu - test selected model ----
summary(modlistglmCu[[30]])
GAMrescheck(modlistglmCu[[30]])
regdiagnostic_customtab(mod=modlistglmCu[[30]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(subdatCu[, Cu], fitted(modlistglmCu[[30]]))
MAE(subdatCu[, Cu], fitted(modlistglmCu[[30]]))
RMSE(subdatCu[, Cu], fitted(modlistglmCu[[30]]))

qplot(subdatCu[, heatbing1902log500proj],
      modlistglmCu[[30]]$residuals) +
  geom_smooth(method='lm', color='red')

qplot(subdatCu[, heatsubspdlpow50_1],
      modlistglmCu[[30]]$residuals) +
  geom_smooth(method='lm', color='red')

qplot(subdatCu[, heatsubslopepow300_1],
      modlistglmCu[[30]]$residuals) +
  geom_smooth(method='lm', color='red')

qplot(subdatCu[, NLCD_reclass_final_PS],
      modlistglmCu[[30]]$residuals) +
  geom_smooth(method='lm', color='red')

#Compare with predictions without transformation
qplot(abs(fitted(modlistglmCu[[30]])-subdatCu[, Zn]), 
      abs(fitted(modlistCu[[19]])-subdatCu[, Zn])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model glm 29 - glm transit200 + nlcd + AADT300 ', 
       y='Absolute error for model 19 - transit200 + nlcd')

#------ 10. Cu - Check spatial and temporal autocorrelation of residuals for dataset without outliers -------
par(mfrow=c(1,1))
resnorm <- rstandard(mod19_nooutliersub) #Get standardized residuals from model
#Make bubble map of residuals
bubbledat <- data.frame(resnorm, subdat$coords.x1, subdat$coords.x2)
coordinates(bubbledat) <- c("subdat.coords.x1","subdat.coords.x2")
bubble(bubbledat, "resnorm", col = c("blue","red"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")
#Check semi-variogram of residuals
plot(variogram(resnorm~1, bubbledat, cutoff=2000, width=100)) #isotropic
plot(variogram(resnorm~1, bubbledat, cutoff= 2000, width=100, alpha = c(0, 45, 90,135))) #anisotropic
#Check spline correlogram ()
plot(spline.correlog(x=coordinates(bubbledat)[,1], y=coordinates(bubbledat)[,2],
                     z=bubbledat$resnorm, resamp=500, quiet=TRUE, xmax = 5000))
#Compute a spatial weight matrix based on IDW
weightmat_k <- lapply(1:10, function(i) {
  weightmat_IDW(subdat[, .(coords.x1, coords.x2)], knb = i, mindist = 10)}) #Based on 1-10 nearest neighbors
weightmat_all <- weightmat_IDW(subdat[, .(coords.x1, coords.x2)], knb = NULL, mindist = 10) #Based on all points

#Moran plots
#lag_resnorm <- lag.listw(weightmat_all, resnorm) #Can be used to create customized Moran plot by plotting residuals against matrix
moran.plot(resnorm, weightmat_all, labels=subdat[,paste0(SiteID, Pair)], pch=19)
moran.plot(resnorm, weightmat_k[[2]], labels=subdat[,paste0(SiteID, Pair)], pch=19)

#Compute Moran's I
"Should always only use lm.morantest for residuals from regression, see http://r-sig-geo.2731867.n2.nabble.com/Differences-between-moran-test-and-lm-morantest-td7591336.html
for an explanation"
lm.morantest(mod19_nooutliersub, listw = listw2U(weightmat_k[[1]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod19_nooutliersub, listw = listw2U(weightmat_k[[2]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod19_nooutliersub, listw = listw2U(weightmat_k[[3]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod19_nooutliersub, listw = listw2U(weightmat_k[[4]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod19_nooutliersub, listw = listw2U(weightmat_k[[5]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod19_nooutliersub, listw = listw2U(weightmat_all)) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)

#Test for need for spatial regression model using Lagrange Multiplier (LM) tests
lm.LMtests(mod19_nooutliersub, listw = listw2U(weightmat_k[[1]]), test=c("LMerr","RLMerr", "SARMA"))
lm.LMtests(mod19_nooutliersub, listw = listw2U(weightmat_k[[2]]), test=c("LMerr","RLMerr", "SARMA"))
lm.LMtests(mod19_nooutliersub, listw = listw2U(weightmat_k[[3]]), test=c("LMerr","RLMerr", "SARMA"))
lm.LMtests(mod19_nooutliersub, listw = listw2U(weightmat_k[[4]]), test=c("LMerr","RLMerr", "SARMA"))
lm.LMtests(mod19_nooutliersub, listw = listw2U(weightmat_k[[5]]), test=c("LMerr","RLMerr", "SARMA"))
lm.LMtests(mod19_nooutliersub, listw = listw2U(weightmat_all), test=c("LMerr","RLMerr", "SARMA"))

#Spatial simultaneous autoregressive error model estimation with 1 nearest neighbors
sarlm_modCu <- errorsarlm(mod19_nooutliersub$call$formula, data = mod19_nooutliersub$model, 
                          listw = listw2U(weightmat_k[[1]]))
summary(sarlm_modCu)
bptest.sarlm(sarlm_modCu)
qplot(fitted(sarlm_modCu), sarlm_modCu$residuals)

#Compare pseudo-R2
cor(mod19_nooutliersub$model$Cu, fitted(sarlm_modCu))^2
cor(mod19_nooutliersub$model$Cu, fitted(mod19_nooutliersub))^2
#Compare MAE
DescTools::MAE(mod19_nooutliersub$model$Cu, fitted(sarlm_modCu))
DescTools::MAE(mod19_nooutliersub$model$Cu, fitted(mod19_nooutliersub))

#Compare observed~predicted for full-no outlier model and for aspatial and spatial model
spatial_comparisonplot <- ggplot(subdat, aes(x=fitted(sarlm_modCu), y=Cu)) + 
  geom_point(aes(x=fitted(mod19_nooutliersub)), size=2, alpha=1/2, color='orange') +
  geom_point(size=2, alpha=1/2, color='red') + 
  geom_abline(size=1, slope=1, intercept=0, color='red') + 
  #geom_text(aes(label=paste0(SiteID, Pair))) +
  coord_fixed() +
  theme_classic()
spatial_comparisonplot

resnorm_postsarlm <- residuals(sarlm_modCu) #Get standardized residuals from model
#Make bubble map of residuals
bubbledat_postsarlm <- data.frame(resnorm_postsarlm, subdat$coords.x1, subdat$coords.x2)
coordinates(bubbledat_postsarlm) <- c("subdat.coords.x1","subdat.coords.x2")
bubble(bubbledat_postsarlm, "resnorm_postsarlm", col = c("blue","red"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")



### FINAL VERDICT: GO WITH MOD 19 TRAINED WITHOUT OUTLIERS-----

#--------------- C. Synthetic pollution index (PI) ~ separate predictors ----
subdatPI <- pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers | SiteIDPair == '107A'),]

#------ 1. PI - Single parameter models --------
modlistPI <- list() #List to hold models
modlistPI[[1]] <- lm(pollution_index ~ 1, data = subdatPI) #Null/Intercept model

modlistPI[[2]] <- lm(pollution_index ~ nlcd_imp_ps_mean, data = subdatPI)
ols_regress(modlistPI[[2]])
#ols_plot_diagnostics(modlistPI[[2]])
ggplot(subdatPI, aes(x=nlcd_imp_ps_mean, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[3]] <- lm(pollution_index ~ nlcd_imp_ps, data = subdatPI)
ols_regress(modlistPI[[3]])
#ols_plot_diagnostics(modlistPI[[3]])
ggplot(subdatPI, aes(x=nlcd_imp_ps, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[4]] <- lm(pollution_index ~ heatbing1902log300proj, data = subdatPI)
ols_regress(modlistPI[[4]])
#ols_plot_diagnostics(modlistPI[[4]])
ggplot(subdatPI, aes(x=heatbing1902log300proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[5]] <- lm(pollution_index ~ heatbing1902log200proj, data = subdatPI)
ols_regress(modlistPI[[5]])
#ols_plot_diagnostics(modlistPI[[5]])
ggplot(subdatPI, aes(x=heatbing1902log200proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[6]] <- lm(pollution_index ~ heatbing1902log500proj, data = subdatPI)
ols_regress(modlistPI[[6]])
#ols_plot_diagnostics(modlistPI[[6]])
ggplot(subdatPI, aes(x=heatbing1902log500proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[7]] <- lm(pollution_index ~ heatbing1902pow200_1proj, data = subdatPI)
ols_regress(modlistPI[[7]])
#ols_plot_diagnostics(modlistPI[[7]])
ggplot(subdatPI, aes(x=heatbing1902pow200_1proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[8]] <- lm(pollution_index ~ heatbing1902pow300_1proj, data = subdatPI)
ols_regress(modlistPI[[8]])
#ols_plot_diagnostics(modlistPI[[8]])
ggplot(subdatPI, aes(x=heatbing1902pow200_1proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[9]] <- lm(pollution_index ~ heatbustransitpow200_1, data = subdatPI)
ols_regress(modlistPI[[9]])
#ols_plot_diagnostics(modlistPI[[9]])
ggplot(subdatPI, aes(x=heatbustransitpow200_1, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[10]] <- lm(pollution_index ~ heatbustransitpow200_2, data = subdatPI)
ols_regress(modlistPI[[10]])
#ols_plot_diagnostics(modlistPI[[10]])
ggplot(subdatPI, aes(x=heatbustransitpow200_2, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[11]] <- lm(pollution_index ~ heatsubslopelog500, data = subdatPI)
ols_regress(modlistPI[[11]])
#ols_plot_diagnostics(modlistPI[[11]])
ggplot(subdatPI, aes(x=heatsubslopelog500, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[12]] <- lm(pollution_index ~ heatsubspdlpow500_2, data = subdatPI)
ols_regress(modlistPI[[12]])
#ols_plot_diagnostics(modlistPI[[12]])
ggplot(subdatPI, aes(x=heatsubspdlpow500_2, y=pollution_index, label=SiteIDPair)) + 
  geom_text()


modlistPI[[13]] <- lm(pollution_index ~ heatsubspdllog500, data = subdatPI)
ols_regress(modlistPI[[13]])
#ols_plot_diagnostics(modlistPI[[13]])
ggplot(subdatPI, aes(x=heatsubspdllog500, y=pollution_index, label=SiteIDPair)) + 
  geom_text()


modlistPI[[14]] <- lm(pollution_index ~ heatsubAADTpow500_2, data = subdatPI)
ols_regress(modlistPI[[14]])
#ols_plot_diagnostics(modlistPI[[14]])
ggplot(subdatPI, aes(x=heatsubAADTlog200^(1/3), y=pollution_index, label=SiteIDPair)) + 
  geom_text()


modlistPI[[15]] <- lm(pollution_index ~ heatsubAADTlog500, data = subdatPI)
ols_regress(modlistPI[[15]])
#ols_plot_diagnostics(modlistPI[[15]])
ggplot(subdatPI, aes(x=heatsubAADTlog500, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[16]] <- lm(pollution_index ~ heatsubAADTlog300, data = subdatPI)
ols_regress(modlistPI[[16]])
#ols_plot_diagnostics(modlistPI[[16]])
ggplot(subdatPI, aes(x=heatsubAADTlog300, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[17]] <- lm(pollution_index ~ heatbustransitpow200_1thd, data = subdatPI)
ols_regress(modlistPI[[17]])
#ols_plot_diagnostics(modlistPI[[17]])
ggplot(subdatPI, aes(x=heatbustransitpow200_1^(1/3), y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[18]] <- lm(pollution_index ~ heatbustransitpow100_1thd, data = subdatPI)
ols_regress(modlistPI[[18]])
#ols_plot_diagnostics(modlistPI[[18]])
ggplot(subdatPI, aes(x=heatbustransitpow100_1^(1/3), y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[19]] <- lm(pollution_index ~ heatbustransitpow100_2thd, data = subdatPI)
ols_regress(modlistPI[[19]])
#ols_plot_diagnostics(modlistPI[[19]])
ggplot(subdatPI, aes(x=heatbustransitpow100_2^(1/3), y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[20]] <- lm(pollution_index ~ heatsubAADTlog200thd, data = subdatPI)
ols_regress(modlistPI[[20]])
#ols_plot_diagnostics(modlistPI[[20]])
ggplot(subdatPI, aes(x=heatsubAADTlog200thd, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[21]] <- lm(pollution_index ~ heatsubAADTlog100thd, data = subdatPI)
ols_regress(modlistPI[[21]])
#ols_plot_diagnostics(modlistPI[[21]])
ggplot(subdatPI, aes(x=heatsubAADTlog100thd, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

#------ 2. PI - Multiparameter models --------
modlistPI[[22]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1thd, data = subdatPI)
ols_regress(modlistPI[[22]])
#ols_plot_diagnostics(modlistPI[[22]])
ols_coll_diag(modlistPI[[22]])
ols_correlations(modlistPI[[22]])

modlistPI[[23]] <- lm(pollution_index ~ heatsubAADTlog200thd + heatbustransitpow100_1thd, data = subdatPI)
ols_regress(modlistPI[[23]])
#ols_plot_diagnostics(modlistPI[[23]])
ols_coll_diag(modlistPI[[23]])
ols_correlations(modlistPI[[23]])

modlistPI[[24]] <- lm(pollution_index ~ heatsubAADTlog100thd + nlcd_imp_ps_mean, data = subdatPI)
ols_regress(modlistPI[[24]])
#ols_plot_diagnostics(modlistPI[[24]])
ols_coll_diag(modlistPI[[24]])
ols_correlations(modlistPI[[24]])

modlistPI[[25]] <- lm(pollution_index ~ heatsubAADTlog200thd + nlcd_imp_ps_mean, data = subdatPI)
ols_regress(modlistPI[[25]])
#ols_plot_diagnostics(modlistPI[[25]])
ols_coll_diag(modlistPI[[25]])
ols_correlations(modlistPI[[25]])

modlistPI[[26]] <- lm(pollution_index ~ heatsubAADTlog100thd*nlcd_imp_ps_mean, data = subdatPI)
ols_regress(modlistPI[[26]])
#ols_plot_diagnostics(modlistPI[[26]])
ols_coll_diag(modlistPI[[26]])
ols_correlations(modlistPI[[26]])

modlistPI[[27]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbing1902log300proj, data = subdatPI)
ols_regress(modlistPI[[27]])
#ols_plot_diagnostics(modlistPI[[27]])
ols_coll_diag(modlistPI[[27]])
ols_correlations(modlistPI[[27]])


modlistPI[[28]] <- lm(pollution_index ~ heatsubAADTlog100thd*heatbing1902log300proj, data = subdatPI)
ols_regress(modlistPI[[28]])
#ols_plot_diagnostics(modlistPI[[28]])
ols_coll_diag(modlistPI[[28]])
ols_correlations(modlistPI[[28]])

modlistPI[[29]] <- lm(pollution_index ~ heatbustransitpow100_1thd + nlcd_imp_ps_mean, data = subdatPI)
ols_regress(modlistPI[[29]])
#ols_plot_diagnostics(modlistPI[[29]])
ols_coll_diag(modlistPI[[29]])
ols_correlations(modlistPI[[29]])

modlistPI[[30]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean, data = subdatPI)
ols_regress(modlistPI[[30]])
#ols_plot_diagnostics(modlistPI[[30]])
ols_coll_diag(modlistPI[[30]])
ols_correlations(modlistPI[[30]])

modlistPI[[31]] <- lm(pollution_index ~ heatsubAADTlog100 + heatbustransitpow100_1thd + nlcd_imp_ps_mean, data = subdatPI)
ols_regress(modlistPI[[31]])
#ols_plot_diagnostics(modlistPI[[30]])
ols_coll_diag(modlistPI[[30]])
ols_correlations(modlistPI[[30]])

modlistPI[[32]] <- lm(pollution_index ~ heatsubAADTlog100*heatbustransitpow100_1 + nlcd_imp_ps_mean, data = subdatPI)
ols_regress(modlistPI[[32]])
#ols_plot_diagnostics(modlistPI[[32]])
ols_coll_diag(modlistPI[[32]])
ols_correlations(modlistPI[[32]])

modlistPI[[33]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean + heatbing1902log300proj, data = subdatPI)
ols_regress(modlistPI[[33]])
#ols_plot_diagnostics(modlistPI[[33]])
ols_coll_diag(modlistPI[[33]])
ols_correlations(modlistPI[[33]])

modlistPI[[34]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean + heatsubspdllog500, data = subdatPI)
ols_regress(modlistPI[[34]])
#ols_plot_diagnostics(modlistPI[[34]])
ols_coll_diag(modlistPI[[34]])
ols_correlations(modlistPI[[34]])

modlistPI[[35]] <- lm(pollution_index ~ heatsubAADTlog100thd*heatsubspdllog500 + heatbustransitpow100_1 + nlcd_imp_ps_mean, data = subdatPI)
ols_regress(modlistPI[[35]])
#ols_plot_diagnostics(modlistPI[[35]])
ols_coll_diag(modlistPI[[35]])
ols_correlations(modlistPI[[35]])


#------ 3. PI - Make latex model summary table ----
vnum <- max(sapply(modlistPI, function(mod) {length(mod$coefficients)}))
model_summary<- as.data.table(
  ldply(modlistPI, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary, -R2pred, AICc)  
cat(latex_format(model_summary), file = file.path(moddir, 'modeltable_pollutionindex_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_pollutionindex_2019.tex')

#------ 4. PI - Make latex model summary table when excluding outliers ----
model_summary_nooutliers <- as.data.table(
  ldply(modlistPI, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, 
                                                         remove_outliers = 'outliers & leverage',
                                                         labelvec = subdatPI[, paste0(SiteID, Pair)],
                                                         kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summary_nooutliers), file = file.path(moddir, 'modeltable_pollutionindex_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_pollutionindex_nooutliers_2019.tex')

#------ 5. log PI - Run models for table -----
modlistlogPI <- list() #List to hold models
modlistlogPI[[1]] <- lm(logPI ~ 1, data = subdatPI) #Null/Intercept model

modlistlogPI[[2]] <- lm(logPI ~ nlcd_imp_ps_mean, data = subdatPI)

modlistlogPI[[3]] <- lm(logPI ~ nlcd_imp_ps, data = subdatPI)

modlistlogPI[[4]] <- lm(logPI ~ heatbing1902log300proj, data = subdatPI)

modlistlogPI[[5]] <- lm(logPI ~ heatbing1902log200proj, data = subdatPI)

modlistlogPI[[6]] <- lm(logPI ~ heatbing1902log500proj, data = subdatPI)

modlistlogPI[[7]] <- lm(logPI ~ heatbing1902pow200_1proj, data = subdatPI)

modlistlogPI[[8]] <- lm(logPI ~ heatbing1902pow300_1proj, data = subdatPI)

modlistlogPI[[9]] <- lm(logPI ~ heatbustransitpow200_1, data = subdatPI)

modlistlogPI[[10]] <- lm(logPI ~ heatbustransitpow200_2, data = subdatPI)

modlistlogPI[[11]] <- lm(logPI ~ heatsubslopelog500, data = subdatPI)

modlistlogPI[[12]] <- lm(logPI ~ heatsubspdlpow500_2, data = subdatPI)

modlistlogPI[[13]] <- lm(logPI ~ heatsubspdllog500, data = subdatPI)

modlistlogPI[[14]] <- lm(logPI ~ heatsubAADTpow500_2, data = subdatPI)

modlistlogPI[[15]] <- lm(logPI ~ heatsubAADTlog500, data = subdatPI)

modlistlogPI[[16]] <- lm(logPI ~ heatsubAADTlog300, data = subdatPI)

modlistlogPI[[17]] <- lm(logPI ~ heatbustransitpow200_1thd, data = subdatPI)

modlistlogPI[[18]] <- lm(logPI ~ heatbustransitpow100_1thd, data = subdatPI)

modlistlogPI[[19]] <- lm(logPI ~ heatbustransitpow100_2thd, data = subdatPI)

modlistlogPI[[20]] <- lm(logPI ~ heatsubAADTlog200thd, data = subdatPI)

modlistlogPI[[21]] <- lm(logPI ~ heatsubAADTlog100thd, data = subdatPI)

modlistlogPI[[22]] <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1thd, data = subdatPI)

modlistlogPI[[23]] <- lm(logPI ~ heatsubAADTlog200thd + heatbustransitpow100_1thd, data = subdatPI)

modlistlogPI[[24]] <- lm(logPI ~ heatsubAADTlog100thd + nlcd_imp_ps_mean, data = subdatPI)

modlistlogPI[[25]] <- lm(logPI ~ heatsubAADTlog200thd + nlcd_imp_ps_mean, data = subdatPI)

modlistlogPI[[26]] <- lm(logPI ~ heatsubAADTlog100thd*nlcd_imp_ps_mean, data = subdatPI)

modlistlogPI[[27]] <- lm(logPI ~ heatsubAADTlog100thd + heatbing1902log300proj, data = subdatPI)

modlistlogPI[[28]] <- lm(logPI ~ heatsubAADTlog100thd*heatbing1902log300proj, data = subdatPI)

modlistlogPI[[29]] <- lm(logPI ~ heatbustransitpow100_1thd + nlcd_imp_ps_mean, data = subdatPI)

modlistlogPI[[30]] <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean, data = subdatPI)

modlistlogPI[[31]] <- lm(logPI ~ heatsubAADTlog100 + heatbustransitpow100_1thd + nlcd_imp_ps_mean, data = subdatPI)

modlistlogPI[[32]] <- lm(logPI ~ heatsubAADTlog100*heatbustransitpow100_1 + nlcd_imp_ps_mean, data = subdatPI)

modlistlogPI[[33]] <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean + heatbing1902log300proj, data = subdatPI)

modlistlogPI[[34]] <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean + heatsubspdllog500, data = subdatPI)

modlistlogPI[[35]] <- lm(logPI ~ heatsubAADTlog100thd*heatsubspdllog500 + heatbustransitpow100_1 + nlcd_imp_ps_mean, data = subdatPI)

#------ 3. log PI - Make latex model summary table ----
model_summary_nooutliers <- as.data.table(
  ldply(modlistPI, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, 
                                                          remove_outliers = 'outliers & leverage',
                                                          labelvec = subdatPI[, paste0(SiteID, Pair)],
                                                          kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summary_nooutliers), file = file.path(moddir, 'modeltable_pollutionindex_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_pollutionindex_nooutliers_2019.tex')

vnum <- max(sapply(modlistlogPI, function(mod) {length(mod$coefficients)}))
model_summary<- as.data.table(
  ldply(modlistlogPI, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary, -R2pred, AICc)  
cat(latex_format(model_summary), file = file.path(moddir, 'modeltable_logpollutionindex_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_logpollutionindex_2019.tex')

#------ 4. log PI - Make latex model summary table when excluding outliers ----
model_summary_nooutliers <- as.data.table(
  ldply(modlistlogPI, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, 
                                                          remove_outliers = 'outliers & leverage',
                                                          labelvec = subdatPI[, paste0(SiteID, Pair)],
                                                          kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summary_nooutliers), file = file.path(moddir, 'modeltable_logpollutionindex_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_logpollutionindex_nooutliers_2019.tex')

#------ 9. Compare final selected models ----
logmod34_nooutliers <- regdiagnostic_customtab(modlistlogPI[[34]], maxpar=vnum, 
                                            remove_outliers = 'outliers',
                                            labelvec = subdatPI[, paste0(SiteID, Pair)],
                                            kCV = TRUE, k=10, cvreps=50)
subdat <- subdatPI[!(paste0(SiteID, Pair) %in% 
                                    strsplit(gsub('\\\\', '', logmod34_nooutliers['outliers']), ',')$outliers),]
logmod34_nooutliersub <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean + heatsubspdllog500, 
                          data = subdat)
summary(logmod34_nooutliersub)
logmod34_nooutliers['outliers']
ols_plot_diagnostics(logmod34_nooutliersub)
logmod34_nooutliers
AICc(logmod34_nooutliersub)
qplot(subdat$heatbing1902log500proj, logmod34_nooutliersub$residuals) +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat$heatsubslopelog500, logmod34_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')


logmod30_nooutliers <- regdiagnostic_customtab(modlistlogPI[[30]], maxpar=vnum, 
                                               remove_outliers = 'outliers',
                                               labelvec = subdatPI[, paste0(SiteID, Pair)],
                                               kCV = TRUE, k=10, cvreps=50)
subdat <- subdatPI[!(paste0(SiteID, Pair) %in% 
                       strsplit(gsub('\\\\', '', logmod30_nooutliers['outliers']), ',')$outliers),]
logmod30_nooutliersub <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean, 
                            data = subdat)
summary(logmod30_nooutliersub)
logmod30_nooutliers['outliers']
ols_plot_diagnostics(logmod30_nooutliersub)
logmod30_nooutliers
AICc(logmod30_nooutliersub)
qplot(subdat$heatbing1902log500proj, logmod30_nooutliersub$residuals) +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat$heatsubslopelog500, logmod30_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')


ggplot(pollutfieldclean_cast, aes(x=exp(predict(logmod34_nooutliersub, pollutfieldclean_cast, type='response')), 
                     y=pollution_index)) +
  geom_point(alpha=1/2, size=2) + 
  geom_point(data = subdatPI, aes(x = exp(predict(logmod30_nooutliersub, subdatPI))), color='red', alpha=1/2, size=2) + 
  geom_abline(intercept=0, slope=1) +
  theme_classic()

#Compare with predictions without transformation
qplot(abs(exp(fitted(modlistlogPI[[34]]))-subdatPI[, pollution_index]), 
      abs(exp(fitted(modlistlogPI[[30]]))-subdatPI[, pollution_index])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model glm 34 - glm  nlcd_imp_ps_mean + AADTlog100thd + transitpow200_1thd ', 
       y='Absolute error for model 30 - nlcd_imp_ps_mean + AADTlog100thd + transitpow100_1')

qplot(abs(fitted(modlistPI[[30]])-subdatPI[, pollution_index]), 
      abs(exp(predict(logmod30_nooutliersub, subdatPI))-subdatPI[, pollution_index])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model 30', 
       y='Absolute error for model log 30')

#------ 5. GLM PI - Single parameter models ------
modlistglmPI <- list() #List to hold models
modlistglmPI[[1]] <- glm(pollution_index ~ 1, data = subdatPI, family = Gamma('log')) #Null/Intercept model

modlistglmPI[[2]] <- glm(pollution_index ~ nlcd_imp_ps_mean, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[2]])

modlistglmPI[[3]] <- glm(pollution_index ~ nlcd_imp_ps, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[3]])

modlistglmPI[[4]] <- glm(pollution_index ~ heatbing1902log300proj, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[4]])

modlistglmPI[[5]] <- glm(pollution_index ~ heatbing1902log500proj, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[5]])

modlistglmPI[[6]] <- glm(pollution_index ~ heatbing1902log500projsqrt, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[6]])

modlistglmPI[[7]] <- glm(pollution_index ~ heatbing1902pow200_1proj, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[7]])

modlistglmPI[[8]] <- glm(pollution_index ~ heatbing1902pow500_1proj, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[8]])

modlistglmPI[[9]] <- glm(pollution_index ~ heatbustransitpow200_1, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[9]])

modlistglmPI[[10]] <- glm(pollution_index ~ heatbustransitpow200_2, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[10]])

modlistglmPI[[11]] <- glm(pollution_index ~ heatsubslopelog500, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[11]])

modlistglmPI[[12]] <- glm(pollution_index ~ heatsubspdllog500, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[12]])

modlistglmPI[[13]] <- glm(pollution_index ~ heatsubspdllog500sqrt, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[13]])

modlistglmPI[[14]] <- glm(pollution_index ~ heatsubAADTpow500_2, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[14]])

modlistglmPI[[15]] <- glm(pollution_index ~ heatsubAADTlog500, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[15]])

modlistglmPI[[16]] <- glm(pollution_index ~ heatsubAADTlog300, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[16]])

modlistglmPI[[17]] <- glm(pollution_index ~ heatbustransitpow200_1thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[17]])

modlistglmPI[[18]] <- glm(pollution_index ~ heatbustransitpow100_1thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[18]])

modlistglmPI[[19]] <- glm(pollution_index ~ heatbustransitpow100_2thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[19]])

modlistglmPI[[20]] <- glm(pollution_index ~ heatsubAADTlog200thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[20]])

modlistglmPI[[21]] <- glm(pollution_index ~ heatsubAADTlog100thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[21]])

#------ 6. GLM PI - Multiparameter models --------
modlistglmPI[[22]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[22]])

modlistglmPI[[23]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog200thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[23]])

modlistglmPI[[24]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbustransitpow100_1thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[24]])

modlistglmPI[[25]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbustransitpow200_1thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[25]])

modlistglmPI[[26]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubspdllog500sqrt, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[26]])

modlistglmPI[[27]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbing1902log500proj, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[27]])

modlistglmPI[[28]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbing1902log500projsqrt, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[28]])

modlistglmPI[[29]] <- glm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[29]])

modlistglmPI[[30]] <- glm(pollution_index ~ heatsubAADTlog200thd + heatbustransitpow100_1, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[30]])

modlistglmPI[[31]] <- glm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[31]])

modlistglmPI[[32]] <- glm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow200_1thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[32]])

modlistglmPI[[33]] <- glm(pollution_index ~ heatsubAADTlog100 + heatbustransitpow100_1thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[33]])

modlistglmPI[[34]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd + heatbustransitpow200_1thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[34]])

modlistglmPI[[35]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd + heatbing1902log500projsqrt, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[35]])

modlistglmPI[[36]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbustransitpow200_1thd + heatbing1902log500projsqrt, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[36]])

modlistglmPI[[37]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd + heatsubspdllog500sqrt, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[37]])

modlistglmPI[[38]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd + heatbustransitpow200_1thd + heatsubspdllog500sqrt, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[38]])

modlistglmPI[[39]] <- glm(pollution_index ~ nlcd_imp_ps_mean*heatsubAADTlog100thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[39]])

modlistglmPI[[40]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd* heatbing1902log500projsqrt, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[40]])

modlistglmPI[[41]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd* heatbustransitpow200_1thd, data = subdatPI, family = Gamma('log'))
GAMrescheck(modlistglmPI[[41]])



#------ 7. GLM PI - Make latex model summary table ----
vnum <- max(sapply(modlistglmPI, function(mod) {length(mod$coefficients)}))
model_summary<- as.data.table(
  ldply(modlistglmPI, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary, -R2pred, AICc)  
cat(latex_format(model_summary), file = file.path(moddir, 'modeltable_glmpollutionindex_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_glmpollutionindex_2019.tex')

#------ 8. GLM PI - Make latex model summary table when excluding outliers ----
model_summary_nooutliers <- as.data.table(
  ldply(modlistglmPI, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, 
                                                         remove_outliers = 'outliers',
                                                         labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                                         kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summary_nooutliers), file = file.path(moddir, 'modeltable_glmpollutionindex_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_glmpollutionindex_nooutliers_2019.tex')

#------ 9. Compare final selected models ----
mod34_nooutliers <- regdiagnostic_customtab(modlistglmPI[[34]], maxpar=vnum, 
                                            remove_outliers = 'outliers',
                                            labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                            kCV = TRUE, k=10, cvreps=50)
subdat <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                                    strsplit(gsub('\\\\', '', mod34_nooutliers['outliers']), ',')$outliers),]
mod34_nooutliersub <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd + heatbustransitpow200_1thd, 
                         data = subdatPI, family=Gamma(link='log'))
GAMrescheck(mod34_nooutliersub)
AICc(mod34_nooutliersub)
qplot(subdatPI$heatbing1902log500proj, mod34_nooutliersub$residuals) +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdatPI$heatsubspdllog500sqrt, mod34_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')


mod30_nooutliers <- regdiagnostic_customtab(modlistPI[[30]], maxpar=vnum, 
                                            remove_outliers = 'outliers',
                                            labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                            kCV = TRUE, k=10, cvreps=50)
subdat <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                                    strsplit(gsub('\\\\', '', mod30_nooutliers['outliers']), ',')$outliers),]
mod30_nooutliersub <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean, 
                         data = subdatPI)
ols_regress(mod30_nooutliersub)
ols_plot_diagnostics(mod30_nooutliersub)
ols_coll_diag(mod30_nooutliersub)
ols_correlations(mod30_nooutliersub)
AICc(mod30_nooutliersub)
qplot(subdatPI$heatbing1902log500proj, mod30_nooutliersub$residuals) +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdatPI$heatsubspdllog500sqrt, mod30_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')


ggplot(subdatPI, aes(x=predict(mod34_nooutliersub, subdatPI, type='response'), 
                                  y=pollution_index)) +
  geom_point(alpha=1/2, size=2) + 
  geom_point(aes(x = predict(mod30_nooutliersub, subdatPI)), color='red', alpha=1/2, size=2) + 
  geom_abline(intercept=0, slope=1) +
  theme_classic()

 #Compare with predictions without transformation
qplot(abs(fitted(modlistglmPI[[34]])-subdatPI[, pollution_index]), 
      abs(fitted(modlistPI[[30]])-subdatPI[, pollution_index])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model glm 34 - glm  nlcd_imp_ps_mean + AADTlog100thd + transitpow200_1thd ', 
       y='Absolute error for model 30 - nlcd_imp_ps_mean + AADTlog100thd + transitpow100_1')

#------ 10. Check spatial and temporal autocorrelation of residuals for full and robust datasets -------
"Fron Anselin 2006: ignoring spatially correlated errors is mostly a problem of efficiency, in the
sense that the OLS coefficient standard error estimates are biased, but the
coefficient estimates themselves remain unbiased. However, to the extent that
the spatially correlated errors mask an omitted variable, the consequences of
ignoring this may be more serious."
#Other good resource: https://eburchfield.github.io/files/Spatial_regression_LAB.html

par(mfrow=c(1,1))
resnorm <- rstandard(modlistPI[[30]]) #Get standardized residuals from model
#Make bubble map of residuals
bubbledat <- data.frame(resnorm, subdatPI$coords.x1, subdatPI$coords.x2)
coordinates(bubbledat) <- c("subdatPI.coords.x1","subdatPI.coords.x2")
bubble(bubbledat, "resnorm", col = c("blue","red"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")
#Check semi-variogram of residuals
plot(variogram(resnorm~1, bubbledat, cutoff=2000, width=100)) #isotropic
plot(variogram(resnorm~1, bubbledat, cutoff= 2000, width=100, alpha = c(0, 45, 90,135))) #anisotropic
#Check spline correlogram ()
plot(spline.correlog(x=coordinates(bubbledat)[,1], y=coordinates(bubbledat)[,2],
                     z=bubbledat$resnorm, resamp=500, quiet=TRUE, xmax = 5000))
#Compute a spatial weight matrix based on IDW
weightmat_k <- lapply(1:10, function(i) {
  weightmat_IDW(subdatPI[, .(coords.x1, coords.x2)], knb = i, mindist = 10)}) #Based on 1-10 nearest neighbors
weightmat_all <- weightmat_IDW(subdatPI[, .(coords.x1, coords.x2)], knb = NULL, mindist = 10) #Based on all points

#Moran plots
#lag_resnorm <- lag.listw(weightmat_all, resnorm) #Can be used to create customized Moran plot by plotting residuals against matrix
moran.plot(resnorm, weightmat_all, labels=subdatPI[,paste0(SiteID, Pair)], pch=19)
moran.plot(resnorm, weightmat_k[[2]], labels=subdatPI[,paste0(SiteID, Pair)], pch=19)

#Compute Moran's I
"Should always only use lm.morantest for residuals from regression, see http://r-sig-geo.2731867.n2.nabble.com/Differences-between-moran-test-and-lm-morantest-td7591336.html
for an explanation"
lm.morantest(modlistPI[[30]], listw = listw2U(weightmat_k[[1]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistPI[[30]], listw = listw2U(weightmat_k[[2]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistPI[[30]], listw = listw2U(weightmat_k[[3]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistPI[[30]], listw = listw2U(weightmat_k[[4]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistPI[[30]], listw = listw2U(weightmat_k[[5]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistPI[[30]], listw = listw2U(weightmat_all)) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)

#Test for need for spatial regression model using Lagrange Multiplier (LM) tests
lm.LMtests(modlistPI[[30]], listw = listw2U(weightmat_k[[1]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(modlistPI[[30]], listw = listw2U(weightmat_k[[2]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(modlistPI[[30]], listw = listw2U(weightmat_k[[3]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(modlistPI[[30]], listw = listw2U(weightmat_k[[4]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))                     
lm.LMtests(modlistPI[[30]], listw = listw2U(weightmat_k[[5]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))  
lm.LMtests(modlistPI[[30]], listw = listw2U(weightmat_all), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))

#Spatial simultaneous autoregressive error model estimation with 1 nearest neighbors
sarlm_modPI <- errorsarlm(modlistPI[[30]]$call$formula, data = modlistPI[[30]]$model, 
                          listw = listw2U(weightmat_k[[1]]))
summary(sarlm_modPI)
bptest.sarlm(sarlm_modPI)

#Compare pseudo-R2
cor(modlistPI[[30]]$model$pollution_index, fitted(sarlm_modPI))^2
cor(modlistPI[[30]]$model$pollution_index, fitted(modlistPI[[30]]))^2
#Compare MAE
DescTools::MAE(modlistPI[[30]]$model$pollution_index, fitted(sarlm_modPI))
DescTools::MAE(modlistPI[[30]]$model$pollution_index, fitted(modlistPI[[30]]))

#Compare observed~predicted for full-no outlier model and for aspatial and spatial model
spatial_comparisonplot <- ggplot(subdatPI, aes(x=fitted(sarlm_modPI), y=pollution_index)) + 
  geom_point(aes(x=fitted(modlistPI[[30]])), size=2, alpha=1/2, color='orange') +
  geom_point(size=2, alpha=1/2, color='red') + 
  geom_point(aes(x=6.992689+3.488446*heatsubAADTlog100thd+0.322152*heatbustransitpow100_1+0.182477*nlcd_imp_ps_mean), 
             size=2,color='darkgreen', alpha=1/2) +
  geom_abline(size=1, slope=1, intercept=0, color='red') + 
  #geom_text(aes(label=paste0(SiteID, Pair))) +
  coord_fixed() +
  theme_classic()
spatial_comparisonplot

resnorm_postsarlm <- residuals(sarlm_modPI) #Get standardized residuals from model
#Make bubble map of residuals
bubbledat_postsarlm <- data.frame(resnorm_postsarlm, subdatPI$coords.x1, subdatPI$coords.x2)
coordinates(bubbledat_postsarlm) <- c("subdatPI.coords.x1","subdatPI.coords.x2")
bubble(bubbledat_postsarlm, "resnorm_postsarlm", col = c("blue","red"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")

###----- FINAL VERDICT: GO WITH SARLM COEFFICIENTS FOR MOD 30-----

###################################################################################################################################
#------ 12. Check how well model predicts site rank -------
subdat[, pollution_index_mod48pred := 
         predict(sarlm_mod48sub, newdata = subdat, pred.type='trend')]
# subdat[, `:=`(pollution_index_mod48predmean = mean(pollution_index_mod48pred),
#               pollution_indexmean = mean(pollution_index)), by=SiteID]
#Plot prediction against observed
ggplot(subdat, aes(x=pollution_index_mod48pred, y=pollution_index)) + 
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  coord_fixed() + 
  theme_classic()

#Plot predicted rank against observed rank
ggplot(subdat, aes(x=rank(pollution_index_mod48pred), y=rank(pollution_index))) + 
  geom_point() +
  geom_abline(intercept=0, slope=1)

#Plot pollution index against rank error
ggplot(subdat, aes(x=pollution_index, y=rank(pollution_index_mod48pred)-rank(pollution_index))) + 
  geom_point() +
  geom_hline(yintercept=0)

#Plot pollution index against absolute rank error
ggplot(subdat, 
       aes(x=pollution_index, y=abs(rank(pollution_index_mod48pred)-rank(pollution_index)))) + 
  geom_point() +
  geom_smooth()
ggplot(subdat,
       aes(x=rank(pollution_index), y=abs(rank(pollution_index_mod48pred)-rank(pollution_index)))) + 
  geom_point() +
  geom_smooth()

#Compute mean rank error
subdat[, sum(abs(rank(pollution_index_mod48pred)-rank(pollution_index)))/.N]
subdat[, sum(abs(rank(pollution_index_mod48pred)-rank(pollution_index)))/(.N^2)] #Relative to total number of ranks


# 6. Export models and data for subsequent analysis
############################################################################################################################################
#--------------- A. Export models ----
exportlist <- list(logZnmod = sarlm_modlogZn, Cumod = sarlm_modCu, PImod = sarlm_modPI) 
                   
                   #logCumod = modlistlogCu[[23]])
saveRDS(exportlist, file = file.path(moddir, 'fieldXRFmodels.rds'))

#Export scaling variable for all pollution drivers
pollutmaxscale <- dcast(pollutfieldclean[!is.na(mean)], 
                        formula = castformula,
                        value.var= 'mean')[, lapply(.SD, max), .SDcols = heatcols]
predtab[SiteIDPair == '123A',]

saveRDS(pollutmaxscale, file = file.path(moddir, 'fieldXRFmodels_scaling.rds'))

#--------------- C. Check predictions  -----
predtab<- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='sitescheck')) 
predtab <- predtab[,SiteIDPair := factor(paste0(F_, Pair))][pollutfieldclean_cast, on='SiteIDPair']
predtab[, `:=`(sarlm_modPIaadt =  sarlm_modPI$coefficients[2]*heatsubAADTlog100thd,
               sarlm_modPIbus = sarlm_modPI$coefficients[3]*i.heatbustransitpow100_1,
               sarlm_modPInlcd = sarlm_modPI$coefficients[4]*nlcd_imp_ps_mean)]

predtab[, sarlmpredpi := sarlm_modPI$coefficients[1]+
          sarlm_modPI$coefficients[2]*heatsubAADTlog100thd+
          sarlm_modPI$coefficients[3]*i.heatbustransitpow100_1 +
          sarlm_modPI$coefficients[4]*nlcd_imp_ps_mean
          ]

ggplot(predtab, aes(x=predpiaadt, y=sarlm_modPIaadt)) + 
  geom_point() + 
  geom_abline()

ggplot(predtab, aes(x=predpinlcd, y=sarlm_modPInlcd)) + 
  geom_point() + 
  geom_abline()

ggplot(predtab, aes(x=predpitransit, y=sarlm_modPIbus)) + 
  geom_point() + 
  geom_abline()

ggplot(predtab, aes(x=predpi30, y=sarlmpredpi)) + 
  geom_point() + 
  geom_abline()


sarlm_modPI$coefficients[1]+sarlm_modPI$coefficients[2]

qplot(pollutfieldclean_cast[, 
                            as.integer(100*(exp(-1.8058766514 + 0.3384508984*heatsubAADTlog100frt + 
                                                  0.0125654964*nlcd_imp_ps_mean + 0.0133221944*heatbing1902log300proj +
                                    -0.0001217836*nlcd_imp_ps_mean*heatbing1902log300proj)+0.5))], pollutfieldclean_cast$Zn) 

check2 <- fread(file.path(resdir, 'testznint2.csv'))
check2[, SiteIDPair := paste0(F_, Pair)]
subdat33_2[,predzn := fitted(modlistglmZn[[33]])]
check2_join <- subdat33_2[check2, on='SiteIDPair']


qplot(check2_join[, predzn33/100], check2_join[,predzn]) + geom_abline()
qplot(check2_join[,i.heatbing1902log300proj], check2_join$heatbing1902log300proj)
qplot(check2_join[,i.nlcd_imp_ps_mean], check2_join$nlcd_imp_ps_mean) + geom_abline()
qplot(check2_join[,testbing], 0.0133221944*check2_join$heatbing1902log300proj) + geom_abline()

qplot(check2_join[,(-0.0001217836*100*i.heatbing1902log300proj/5183)*i.nlcd_imp_ps_mean],
      check2_join[,-0.0001217836*nlcd_imp_ps_mean*heatbing1902log300proj]) + geom_abline()

#TRY THIS - AFTER RERUNNING FROM L4324
qplot(check2_join[,testinterac],
      check2_join[,-0.0001217836*heatbing1902log300proj]) + geom_abline()

check3 <- fread(file.path(resdir, 'testznint3.csv'))
check3[, SiteIDPair := paste0(F_, Pair)]
check3_join <- subdat33_2[check3, on='SiteIDPair']
qplot(check3_join[, testaadt], check3_join[,heatsubAADTlog100frt]) + geom_abline()


qplot(check2[,(-0.0001345*nlcd_imp_ps_mean*heatbing1902log300proj)-checkinterac], check2$checkinterac)

#--------------- B. Export data for Luwam ----
write.csv(pollutfieldclean_cast, 
          file.path(resdir, 'map_forluwam/XRFsites_pollutiondata_forLuwam_20190530.csv'),
          row.names = F)

#--------------- C. Make graph of % area vs % total pollution -----
#This is a preliminary figure to show % area vs % total pollution for study area extent
predzntab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predzn36_tab'))
#Percentile # of cells
predzntab[order(-Value), `:=`(cumcount = cumsum(Count),
                                cumvalue = cumsum((Value)*Count))]
areacumpollution <- unlist(
  sapply(round(seq(0, sum(predzntab$Count), sum(predzntab$Count)/20)), function(x) { 
  predzntab[cumcount >= x, ][which.min(cumcount), cumvalue - (Value)*(cumcount-x)]/predzntab[, max(cumvalue)]
})
)
ggplot(data.table(cumpollution = 100*areacumpollution, percarea = seq(0,100, 5)), aes(x=percarea, y=cumpollution)) +
  geom_bar(stat = 'identity') + 
  scale_y_continuous(expand=c(0,0), name = 'Cumulative Zn pollution (%)') +
  scale_x_continuous(expand=c(0,0), name = 'Cumulative area, from most to least polluted (%)', limits = c(0,105)) +
  coord_fixed() +
  theme_classic() +
  theme(text = element_text(size=14),
        plot.margin= margin(0.25, 0.25, 0.25, 0.25, "cm"))


####################################### JUNK #############################################################
# sitelocs <- SpatialPointsDataFrame(coords = data.frame(pollutfieldclean_cast$POINT_X, pollutfieldclean_cast$POINT_Y),
#                                       data= as.data.frame(pollutfieldclean_cast))
# binpal <- colorBin("Reds", sitelocs$heat_binglog500, 10, pretty = FALSE)
# 
# #View(outlierlocs@data)
# leaflet(data = sitelocs) %>% addTiles() %>%
#   addCircleMarkers(color = ~binpal(heat_binglog500),
#                    #clusterOptions = markerClusterOptions(),
#                    popup = ~paste0(SiteID, Pair))

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


