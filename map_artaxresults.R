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
library(ClustOfVar)
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
                     family = if (class(mod)[1] == 'glm') Gamma(link='log') else NULL,
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
pollutfieldclean <- pollutfieldmerge[!(paste0(SiteID, Pair) %in% c('7A','52A', '114A', '114B', paste0(64:69, 'A'))),] #Remove 114A and 114B as well because they are just duplicate measurements of 51A at another time
pollutlabclean <- pollutlabmerge[!(paste0(SiteID, Pair) %in% c('53A', '49A')),]
pollutICPclean <-  pollutICPmerge[!(paste0(SiteID, Pair) %in% c('49A')),]

############################################################################################################################################


# 4. Check relationship among and reduce variables 
############################################################################################################################################
extraoutliers <- c('51A') #Remove anomalous data points by 15th Ave NE and Luwam's (latter do not appear to be precise enough

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
  ggscatmat(as.data.frame(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), elemcols, with=F]),
            alpha=0.7)
  dev.off()
}


#Correlation heatmap
outfieldheatmat <- file.path(inspectdir, 'corheatmap_FieldElem_FieldElem.png')
if (!file.exists(outfieldheatmat)) {
  png(outfieldheatmat, width = 8, height=8, units='in', res=300)
  corr_heatmap(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), elemcols, with=F])
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

#Variable clustering
elemtree <- hclustvar(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), elemcols, with=F])
png(file.path(resdir, 'data_inspection/varclustree.png'), width = 8, height=4, units='in', res=300)
plot(elemtree)
dev.off()

# stab <- stability(elemtree,B=40)
# plot(stab, main="Stability of the partitions")
# stab$matCR
# boxplot(stab$matCR, main="Dispersion of the ajusted Rand index")
# P14<-cutreevar(elemtree, 14, matsim=TRUE)
# print(P14)
# P14$var

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
                        ,`:=`(pollution_index = 100*Reduce('+', lapply(.SD, function(x) (x/max(x))))/length(cbcols),
                              Znstand = 100*Zn/max(Zn),
                              Custand = 100*Cu/max(Cu),
                              Pbstand = 100*Pb/max(Pb)),
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
#Excluding 51A, 114A and 114B
outfieldpollutheatmat <- file.path(inspectdir, 'corheatmap_PollutionDrivers_FieldElem.png')
if (!file.exists(outfieldpollutheatmat)) {
  png(outfieldpollutheatmat,width = 40, height=35, units='in', res=300)
  corr_heatmap(xmat=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), pollutcols[-which(pollutcols=='NLCD_reclass_final_PS')], with=F],
               ymat=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), c(elemcols, 'pollution_index', 'cbpollution_index'), with=F],
               clus = FALSE) +
    scale_y_discrete(labels=heatlab) + 
    theme(text=element_text(size=22))
  dev.off()
}

xmat <- pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), pollutcols[-which(pollutcols=='NLCD_reclass_final_PS')], with=F]
ymat <- pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), c(elemcols, 'pollution_index', 'cbpollution_index'), with=F]
cordf <- round(cor(x=xmat, y=ymat),2) %>%
  t() %>%
  as.data.frame()
cordf$group <- rownames(cordf)
cordf <- cordf[, c(ncol(cordf),
                   1:(ncol(cordf)-1))]

spidermetals <- c('Fe', 'Cu', 'Zn','Pb')
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
plot.data <- cordf[cordf$group %in% spidermetals,
                   spidercols[order(spidercols)]]
plot.data <- as.data.frame(plot.data)
plot.data[,1] <- as.factor(as.character(plot.data[,1]))
names(plot.data)[1] <- "group"
var.names <- colnames(plot.data)[-1]  #'Short version of variable names
plot.data.offset <- plot.data
plot.data.offset[,2:ncol(plot.data)]<- plot.data[,2:ncol(plot.data)]+abs(centre.y)
#print(plot.data.offset)
# (b) convert into radial coords
CalculateGroupPath <- function(df) {
  #Converts variable values into a set of radial x-y coordinates
  #Code adapted from a solution posted by Tony M to
  #http://stackoverflow.com/questions/9614433/creating-radar-chart-a-k-a-star-plot-spider-plot-using-ggplot2-in-r
  #Args:
  #  df: Col 1 -  group ('unique' cluster / group ID of entity)
  #      Col 2-n:  v1.value to vn.value - values (e.g. group/cluser mean or median) of variables v1 to v.n

  path <- df[,1]

  ##find increment
  angles = seq(from=0, to=2*pi, by=(2*pi)/(ncol(df)-1))
  ##create graph data frame
  graphData= data.frame(seg="", x=0,y=0)
  graphData=graphData[-1,]

  for(i in levels(path)){
    pathData = subset(df, df[,1]==i)
    for(j in c(2:ncol(df))){
      #pathData[,j]= pathData[,j]


      graphData=rbind(graphData, data.frame(group=i,
                                            x=pathData[,j]*sin(angles[j-1]),
                                            y=pathData[,j]*cos(angles[j-1])))
    }
    ##complete the path by repeating first pair of coords in the path
    graphData=rbind(graphData, data.frame(group=i,
                                          x=pathData[,2]*sin(angles[1]),
                                          y=pathData[,2]*cos(angles[1])))
  }
  #Make sure that name of first column matches that of input data (in case !="group")
  colnames(graphData)[1] <- colnames(df)[1]
  graphData #data frame returned by function
}
group <-NULL
group$path <- CalculateGroupPath(plot.data.offset)

#Create palette
library(RColorBrewer)
interspec <- colorRampPalette(brewer.pal(11, 'Spectral'))
intdat<- data.frame(y=round(group$path$y,2))
datrange <- seq(min(intdat), max(intdat), 0.01)
coldf <- data.frame(y = datrange, ycol=interspec(length(datrange)))
colvalues <- merge(intdat, coldf, on='y', all.y=F)
colvalues$yfac <- factor(as.character(colvalues$y),
                         levels= unique(as.character(colvalues$y)[order(-colvalues$y)]))

ggplot(colvalues, aes(x=y, y=y, color=yfac)) +
  geom_point() +
  scale_color_manual(values=unique(as.character(colvalues$ycol)))

ggradarplot <- ggradar(plot.data=cordf[cordf$group %in% c('Zn', 'Cu', 'Pb'),
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
                       #group.colours = '#8CE071',
                       background.circle.transparency=0.15) +
  geom_path(data=gridline75, aes(x=x, y=y), color='darkgrey') +
  geom_path(data=gridline25, aes(x=x, y=y), color='darkgrey') +
  # scale_colour_manual(values=unique(as.character(colvalues$ycol))) +
  # facet_wrap(~group) +
  theme_minimal() +
  theme(text= element_text(size = 18),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position=c(0.15, 0.15),
        plot.margin=unit(c(-1.5, -1.5, -1.5, -1.5),'cm'))
ggradarplot


extrafont::loadfonts()
pdf(file.path(resdir, 'data_inspection/spiderplot_20191219.pdf'), width=9, height=9)
#png(file.path(moddir, 'spiderplotZn_20190514_1.png'), width=9, height=9, units='in', res=450)
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
                             heatbing1902log300projfrt = (heatbing1902log300proj)^0.25,
                             heatbing1902log300projsqrt = sqrt(heatbing1902log300proj),
                             heatbing1902log500projsqrt = sqrt(heatbing1902log500proj),
                             heatsubspdlpow100_3sqrt = sqrt(heatsubspdlpow100_3),
                             heatsubspdlpow100_3frt = (heatsubspdlpow100_3)^0.25,
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
                             logZnstand = log(Znstand),
                             logPI = log(pollution_index)
)]

pollutfieldclean_cast[, NLCD_reclass_final_PS_edit := NLCD_reclass_final_PS]
pollutfieldclean_cast[!(NLCD_reclass_final_PS %in% c('96', '97', '98', '99')), NLCD_reclass_final_PS_edit := '1']

#--------------- A. Znstand ~ separate predictors ----
modlistZnstand <- list() #List to hold models
modlistZnstand[[1]] <- lm(Znstand ~ 1, data = pollutfieldclean_cast) #Null/Intercept model
#------ 1. Znstand - Single parameter models --------
modlistZnstand[[2]] <- lm(Znstand ~ heatsubAADTlog100, data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[2]])
#ols_plot_diagnostics(modlistZnstand[[2]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog100^(1/4), y=Znstand, label=paste0(SiteID, Pair),
                                  color = heatbing1902log300proj)) + 
  geom_text() + 
  scale_color_distiller(palette='Spectral')

modlistZnstand[[3]] <- lm(Znstand ~ heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[3]])
#ols_plot_diagnostics(modlistZnstand[[3]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902log300proj, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()
ggplot(pollutfieldclean_cast[as.numeric(SiteID) < 63,], aes(x=heatbing1902log300proj, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistZnstand[[4]] <- lm(Znstand ~ nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[4]])
#ols_plot_diagnostics(modlistZnstand[[4]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps_mean, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistZnstand[[5]] <- lm(Znstand ~ nlcd_imp_ps, data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[5]])
#ols_plot_diagnostics(modlistZnstand[[5]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistZnstand[[6]] <- lm(Znstand ~ heatbustransitlog200, data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[6]])
#ols_plot_diagnostics(modlistZnstand[[6]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitlog200, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistZnstand[[7]] <- lm(Znstand ~ heatsubspdlpow100_3, data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[7]])
#ols_plot_diagnostics(modlistZnstand[[7]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow100_3, y=Znstand, 
                                  label=paste0(SiteID, Pair), color=nlcd_imp_ps_mean)) + 
  geom_text() +
  geom_smooth(method='lm') +
  scale_color_distiller(palette='Spectral')

modlistZnstand[[8]] <- lm(Znstand ~ heatsubslopepow500_1, data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[8]])
#ols_plot_diagnostics(modlistZnstand[[8]])
ggplot(pollutfieldclean_cast, aes(x=heatsubslopepow500_1, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

#------ 2. Znstand - Multiparameter models --------
modlistZnstand[[9]] <- lm(Znstand ~ heatsubspdlpow100_3 + heatsubAADTlog100, 
                     data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[9]])
#ols_plot_diagnostics(modlistZnstand[[9]])
ols_coll_diag(modlistZnstand[[9]])
ols_correlations(modlistZnstand[[9]])

modlistZnstand[[10]] <- lm(Znstand ~ heatsubspdlpow100_3 + heatbing1902log300proj, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[10]])
#ols_plot_diagnostics(modlistZnstand[[10]])
ols_coll_diag(modlistZnstand[[10]])
ols_correlations(modlistZnstand[[10]])

modlistZnstand[[11]] <- lm(Znstand ~  heatsubspdlpow100_3 + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[11]])
#ols_plot_diagnostics(modlistZnstand[[11]])
ols_coll_diag(modlistZnstand[[11]])
ols_correlations(modlistZnstand[[11]])

modlistZnstand[[12]] <- lm(Znstand ~ heatsubspdlpow100_3 + heatbustransitlog200, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[12]])
#ols_plot_diagnostics(modlistZnstand[[12]])
ols_coll_diag(modlistZnstand[[12]])
ols_correlations(modlistZnstand[[12]])

modlistZnstand[[13]] <- lm(Znstand ~ heatsubspdlpow100_3 + heatsubslopepow500_1, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[13]])
#ols_plot_diagnostics(modlistZnstand[[13]])
ols_coll_diag(modlistZnstand[[13]])
ols_correlations(modlistZnstand[[13]])

modlistZnstand[[14]] <- lm(Znstand ~ heatsubspdlpow100_3*heatsubAADTlog100, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[14]])
#ols_plot_diagnostics(modlistZnstand[[14]])
ols_coll_diag(modlistZnstand[[14]])
ols_correlations(modlistZnstand[[14]])

modlistZnstand[[15]] <- lm(Znstand ~ heatsubspdlpow100_3*heatbing1902log300proj, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[15]])
#ols_plot_diagnostics(modlistZnstand[[15]])
ols_coll_diag(modlistZnstand[[15]])
ols_correlations(modlistZnstand[[15]])

modlistZnstand[[16]] <- lm(Znstand ~ heatsubspdlpow100_3*heatbing1902log300proj + nlcd_imp_ps_mean , 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[16]])
#ols_plot_diagnostics(modlistZnstand[[16]])
ols_coll_diag(modlistZnstand[[16]])
ols_correlations(modlistZnstand[[16]])

modlistZnstand[[17]] <- lm(Znstand ~ heatsubspdlpow100_3 + heatbing1902log300proj + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[17]])
#ols_plot_diagnostics(modlistZnstand[[17]])
ols_coll_diag(modlistZnstand[[17]])
ols_correlations(modlistZnstand[[17]])

modlistZnstand[[18]] <- lm(Znstand ~ heatsubspdlpow100_3*nlcd_imp_ps_mean + heatbing1902log300proj,
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[18]])
#ols_plot_diagnostics(modlistZnstand[[18]])
ols_coll_diag(modlistZnstand[18])
ols_correlations(modlistZnstand[[18]])

modlistZnstand[[19]] <- lm(Znstand ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[19]])
#ols_plot_diagnostics(modlistZnstand[[19]])
ols_correlations(modlistZnstand[[19]])

modlistZnstand[[20]] <- lm(Znstand ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean + 
                        heatbustransitlog200, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[20]])
#ols_plot_diagnostics(modlistZnstand[[20]])
ols_coll_diag(modlistZnstand[20])
ols_correlations(modlistZnstand[[20]])


#------ 3. Znstand - Transformed predictors --------
modlistZnstand[[21]] <- lm(Znstand ~ heatsubAADTlog200sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[21]])
#ols_plot_diagnostics(modlistZnstand[[21]])
ols_coll_diag(modlistZnstand[[21]])
ols_correlations(modlistZnstand[[21]])

modlistZnstand[[22]] <- lm(Znstand ~ heatsubspdlpow100_3sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[22]])
#ols_plot_diagnostics(modlistZnstand[[22]])
ols_coll_diag(modlistZnstand[[22]])
ols_correlations(modlistZnstand[[22]])

modlistZnstand[[23]] <- lm(Znstand ~ heatbing1902log300projsqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[23]])
#ols_plot_diagnostics(modlistZnstand[[23]])
ols_coll_diag(modlistZnstand[[23]])
ols_correlations(modlistZnstand[[23]])

modlistZnstand[[24]] <- lm(Znstand ~ heatsubAADTlog300sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[24]])
#ols_plot_diagnostics(modlistZnstand[[24]])
ols_coll_diag(modlistZnstand[[24]])
ols_correlations(modlistZnstand[[24]])

modlistZnstand[[25]] <- lm(Znstand ~ heatsubAADTlog100sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[25]])
#ols_plot_diagnostics(modlistZnstand[[25]])
ols_coll_diag(modlistZnstand[[25]])
ols_correlations(modlistZnstand[[25]])

modlistZnstand[[26]] <- lm(Znstand ~ heatsubAADTlog100sqrt + heatsubspdlpow100_3sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[26]])
#ols_plot_diagnostics(modlistZnstand[[26]])
ols_coll_diag(modlistZnstand[[26]])
ols_correlations(modlistZnstand[[26]])

modlistZnstand[[27]] <- lm(Znstand ~ heatsubAADTlog100sqrt*heatsubspdlpow100_3, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[27]])
#ols_plot_diagnostics(modlistZnstand[[27]])
ols_coll_diag(modlistZnstand[[27]])
ols_correlations(modlistZnstand[[27]])

modlistZnstand[[28]] <- lm(Znstand ~ heatsubAADTlog100sqrt*heatsubspdlpow100_3 + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[28]])
#ols_plot_diagnostics(modlistZnstand[[28]])
ols_coll_diag(modlistZnstand[[28]])
ols_correlations(modlistZnstand[[28]])

modlistZnstand[[29]] <- lm(Znstand ~ (heatsubAADTlog100sqrt)*heatsubspdlpow100_3 + heatbing1902log300proj, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[29]])
#ols_plot_diagnostics(modlistZnstand[[29]])
ols_coll_diag(modlistZnstand[[29]])
ols_correlations(modlistZnstand[[29]])

modlistZnstand[[30]] <- lm(Znstand ~ (heatsubAADTlog100sqrt)*heatsubspdlpow100_3 + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubAADTlog100sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[30]])
#ols_plot_diagnostics(modlistZnstand[[30]])
ols_coll_diag(modlistZnstand[[30]])
ols_correlations(modlistZnstand[[30]])


modlistZnstand[[31]] <- lm(Znstand ~ (heatsubAADTlog100sqrt) + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubAADTlog100sqrt, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[31]])
#ols_plot_diagnostics(modlistZnstand[[31]])
ols_coll_diag(modlistZnstand[[31]])
ols_correlations(modlistZnstand[[31]])

modlistZnstand[[32]] <- lm(Znstand ~ (heatsubspdlpow100_1) + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubspdlpow100_3, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[32]])
#ols_plot_diagnostics(modlistZnstand[[32]])
ols_coll_diag(modlistZnstand[[32]])
ols_correlations(modlistZnstand[[32]])

modlistZnstand[[33]] <- lm(Znstand ~ (heatsubAADTlog100sqrt) + heatbing1902log300proj +
                        heatbing1902log300proj:heatsubAADTlog100sqrt + nlcd_imp_ps, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[33]])
#ols_plot_diagnostics(modlistZnstand[[33]])
ols_coll_diag(modlistZnstand[[33]])
ols_correlations(modlistZnstand[[33]])

modlistZnstand[[34]] <- lm(Znstand ~ (heatsubAADTlog100sqrt) + heatbing1902log300proj + nlcd_imp_ps +
                        heatbustransitlog200, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[34]])
#ols_plot_diagnostics(modlistZnstand[[34]])
ols_coll_diag(modlistZnstand[[34]])
ols_correlations(modlistZnstand[[34]])

modlistZnstand[[35]] <- lm(Znstand ~ (heatsubAADTlog100sqrt) + heatbing1902log300proj + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[35]])
#ols_plot_diagnostics(modlistZnstand[[35]])
ols_coll_diag(modlistZnstand[[35]])
ols_correlations(modlistZnstand[[35]])

modlistZnstand[[36]] <- lm(Znstand ~ heatsubAADTlog100sqrt + nlcd_imp_ps_mean, 
                      data = pollutfieldclean_cast)
ols_regress(modlistZnstand[[36]])
#ols_plot_diagnostics(modlistZnstand[[36]])
ols_coll_diag(modlistZnstand[[36]])
ols_correlations(modlistZnstand[[36]])

#------ 4. Znstand - Make latex model summary table ----
vnum <- max(sapply(modlistZnstand, function(mod) {length(mod$coefficients)}))
model_summaryZnstand<- as.data.table(
  ldply(modlistZnstand, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summaryZnstand, -R2pred, AICc)  
cat(latex_format(model_summaryZnstand), file = file.path(moddir, 'modeltable_Znstand_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_Znstand_2019.tex')

#------ 5. Znstand - Make latex model summary table when excluding outliers----
vnum <- max(sapply(modlistZnstand, function(mod) {length(mod$coefficients)}))
model_summaryZnstand_nooutliers<- as.data.table(
  ldply(modlistZnstand, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summaryZnstand_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summaryZnstand_nooutliers), file = file.path(moddir, 'modeltable_Znstand_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_Znstand_nooutliers_2019.tex')

#------ 6. logZnstand - Single parameter models --------
modlistlogZnstand <- list() #List to hold models
modlistlogZnstand[[1]] <- lm(logZnstand ~ 1, data = pollutfieldclean_cast) #Null/Intercept model

modlistlogZnstand[[2]] <- lm(logZnstand ~ heatsubAADTlog100, data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[2]])
#ols_plot_diagnostics(modlistlogZnstand[[2]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog100^(1/4), y=Znstand, label=paste0(SiteID, Pair),
                                  color = heatbing1902log300proj)) + 
  geom_text() + 
  scale_color_distiller(palette='Spectral')

modlistlogZnstand[[3]] <- lm(logZnstand ~ heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[3]])
#ols_plot_diagnostics(modlistlogZnstand[[3]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902log300proj, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()
ggplot(pollutfieldclean_cast[as.numeric(SiteID) < 63,], aes(x=heatbing1902log300proj, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistlogZnstand[[4]] <- lm(logZnstand ~ nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[4]])
#ols_plot_diagnostics(modlistlogZnstand[[4]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps_mean, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistlogZnstand[[5]] <- lm(logZnstand ~ nlcd_imp_ps, data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[5]])
#ols_plot_diagnostics(modlistlogZnstand[[5]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistlogZnstand[[6]] <- lm(logZnstand ~ heatbustransitlog200, data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[6]])
#ols_plot_diagnostics(modlistlogZnstand[[6]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitlog200, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

modlistlogZnstand[[7]] <- lm(logZnstand ~ heatsubspdlpow100_3, data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[7]])
#ols_plot_diagnostics(modlistlogZnstand[[7]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow100_3, y=Znstand, 
                                  label=paste0(SiteID, Pair), color=nlcd_imp_ps_mean)) + 
  geom_text() +
  geom_smooth(method='lm') +
  scale_color_distiller(palette='Spectral')

modlistlogZnstand[[8]] <- lm(logZnstand ~ heatsubslopepow500_1, data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[8]])
#ols_plot_diagnostics(modlistlogZnstand[[8]])
ggplot(pollutfieldclean_cast, aes(x=heatsubslopepow500_1, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text()

#------ 7. logZnstand - Multiparameter models --------
modlistlogZnstand[[9]] <- lm(logZnstand ~ heatsubspdlpow100_3 + heatsubAADTlog100, 
                        data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[9]])
#ols_plot_diagnostics(modlistlogZnstand[[9]])
ols_coll_diag(modlistlogZnstand[[9]])
ols_correlations(modlistlogZnstand[[9]])

modlistlogZnstand[[10]] <- lm(logZnstand ~ heatsubspdlpow100_3 + heatbing1902log300proj, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[10]])
#ols_plot_diagnostics(modlistlogZnstand[[10]])
ols_coll_diag(modlistlogZnstand[[10]])
ols_correlations(modlistlogZnstand[[10]])

modlistlogZnstand[[11]] <- lm(logZnstand ~  heatsubspdlpow100_3 + nlcd_imp_ps_mean, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[11]])
#ols_plot_diagnostics(modlistlogZnstand[[11]])
ols_coll_diag(modlistlogZnstand[[11]])
ols_correlations(modlistlogZnstand[[11]])

modlistlogZnstand[[12]] <- lm(logZnstand ~ heatsubspdlpow100_3 + heatbustransitlog200, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[12]])
#ols_plot_diagnostics(modlistlogZnstand[[12]])
ols_coll_diag(modlistlogZnstand[[12]])
ols_correlations(modlistlogZnstand[[12]])

modlistlogZnstand[[13]] <- lm(logZnstand ~ heatsubspdlpow100_3 + heatsubslopepow500_1, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[13]])
#ols_plot_diagnostics(modlistlogZnstand[[13]])
ols_coll_diag(modlistlogZnstand[[13]])
ols_correlations(modlistlogZnstand[[13]])

modlistlogZnstand[[14]] <- lm(logZnstand ~ heatsubspdlpow100_3*heatsubAADTlog100, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[14]])
#ols_plot_diagnostics(modlistlogZnstand[[14]])
ols_coll_diag(modlistlogZnstand[[14]])
ols_correlations(modlistlogZnstand[[14]])

modlistlogZnstand[[15]] <- lm(logZnstand ~ heatsubspdlpow100_3*heatbing1902log300proj, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[15]])
#ols_plot_diagnostics(modlistlogZnstand[[15]])
ols_coll_diag(modlistlogZnstand[[15]])
ols_correlations(modlistlogZnstand[[15]])

modlistlogZnstand[[16]] <- lm(logZnstand ~ heatsubspdlpow100_3*heatbing1902log300proj + nlcd_imp_ps_mean , 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[16]])
#ols_plot_diagnostics(modlistlogZnstand[[16]])
ols_coll_diag(modlistlogZnstand[[16]])
ols_correlations(modlistlogZnstand[[16]])

modlistlogZnstand[[17]] <- lm(logZnstand ~ heatsubspdlpow100_3 + heatbing1902log300proj + nlcd_imp_ps_mean, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[17]])
#ols_plot_diagnostics(modlistlogZnstand[[17]])
ols_coll_diag(modlistlogZnstand[[17]])
ols_correlations(modlistlogZnstand[[17]])

modlistlogZnstand[[18]] <- lm(logZnstand ~ heatsubspdlpow100_3*nlcd_imp_ps_mean + heatbing1902log300proj,
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[18]])
#ols_plot_diagnostics(modlistlogZnstand[[18]])
ols_coll_diag(modlistlogZnstand[18])
ols_correlations(modlistlogZnstand[[18]])

modlistlogZnstand[[19]] <- lm(logZnstand ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[19]])
#ols_plot_diagnostics(modlistlogZnstand[[19]])
ols_correlations(modlistlogZnstand[[19]])

modlistlogZnstand[[20]] <- lm(logZnstand ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean + 
                           heatbustransitlog200, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[20]])
#ols_plot_diagnostics(modlistlogZnstand[[20]])
ols_coll_diag(modlistlogZnstand[20])
ols_correlations(modlistlogZnstand[[20]])

#------ 8. logZnstand - Transformed predictors --------
modlistlogZnstand[[21]] <- lm(logZnstand ~ heatsubAADTlog200frt, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[21]])
#ols_plot_diagnostics(modlistlogZnstand[[21]])
ols_coll_diag(modlistlogZnstand[[21]])
ols_correlations(modlistlogZnstand[[21]])

modlistlogZnstand[[22]] <- lm(logZnstand ~ heatsubspdlpow100_3frt, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[22]])
#ols_plot_diagnostics(modlistlogZnstand[[22]])
ols_coll_diag(modlistlogZnstand[[22]])
ols_correlations(modlistlogZnstand[[22]])

modlistlogZnstand[[23]] <- lm(logZnstand ~ heatbing1902log300projfrt, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[23]])
#ols_plot_diagnostics(modlistlogZnstand[[23]])
ols_coll_diag(modlistlogZnstand[[23]])
ols_correlations(modlistlogZnstand[[23]])

modlistlogZnstand[[24]] <- lm(logZnstand ~ heatsubAADTlog300frt, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[24]])
#ols_plot_diagnostics(modlistlogZnstand[[24]])
ols_coll_diag(modlistlogZnstand[[24]])
ols_correlations(modlistlogZnstand[[24]])

modlistlogZnstand[[25]] <- lm(logZnstand ~ heatsubAADTlog100frt, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[25]])
#ols_plot_diagnostics(modlistlogZnstand[[25]])
ols_coll_diag(modlistlogZnstand[[25]])
ols_correlations(modlistlogZnstand[[25]])

modlistlogZnstand[[26]] <- lm(logZnstand ~ heatsubAADTlog100frt + heatsubspdlpow100_3sqrt, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[26]])
#ols_plot_diagnostics(modlistlogZnstand[[26]])
ols_coll_diag(modlistlogZnstand[[26]])
ols_correlations(modlistlogZnstand[[26]])

modlistlogZnstand[[27]] <- lm(logZnstand ~ heatsubAADTlog100frt*heatsubspdlpow100_3, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[27]])
#ols_plot_diagnostics(modlistlogZnstand[[27]])
ols_coll_diag(modlistlogZnstand[[27]])
ols_correlations(modlistlogZnstand[[27]])

modlistlogZnstand[[28]] <- lm(logZnstand ~ heatsubAADTlog100frt*heatsubspdlpow100_3 + nlcd_imp_ps_mean, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[28]])
#ols_plot_diagnostics(modlistlogZnstand[[28]])
ols_coll_diag(modlistlogZnstand[[28]])
ols_correlations(modlistlogZnstand[[28]])

modlistlogZnstand[[29]] <- lm(logZnstand ~ (heatsubAADTlog100frt)*heatsubspdlpow100_3 + heatbing1902log300proj, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[29]])
#ols_plot_diagnostics(modlistlogZnstand[[29]])
ols_coll_diag(modlistlogZnstand[[29]])
ols_correlations(modlistlogZnstand[[29]])

modlistlogZnstand[[30]] <- lm(logZnstand ~ (heatsubAADTlog100frt)*heatsubspdlpow100_3 + heatbing1902log300proj +
                           heatbing1902log300proj:heatsubAADTlog100frt, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[30]])
#ols_plot_diagnostics(modlistlogZnstand[[30]])
ols_coll_diag(modlistlogZnstand[[30]])
ols_correlations(modlistlogZnstand[[30]])


modlistlogZnstand[[31]] <- lm(logZnstand ~ (heatsubAADTlog100frt) + heatbing1902log300proj +
                           heatbing1902log300proj:heatsubAADTlog100frt, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[31]])
#ols_plot_diagnostics(modlistlogZnstand[[31]])
ols_coll_diag(modlistlogZnstand[[31]])
ols_correlations(modlistlogZnstand[[31]])

modlistlogZnstand[[32]] <- lm(logZnstand ~ (heatsubspdlpow100_1) + heatbing1902log300proj +
                           heatbing1902log300proj:heatsubspdlpow100_3, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[32]])
#ols_plot_diagnostics(modlistlogZnstand[[32]])
ols_coll_diag(modlistlogZnstand[[32]])
ols_correlations(modlistlogZnstand[[32]])

modlistlogZnstand[[33]] <- lm(logZnstand ~ (heatsubAADTlog100frt) + heatbing1902log300proj +
                           heatbing1902log300proj:heatsubAADTlog100frt + nlcd_imp_ps, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[33]])
#ols_plot_diagnostics(modlistlogZnstand[[33]])
ols_coll_diag(modlistlogZnstand[[33]])
ols_correlations(modlistlogZnstand[[33]])

modlistlogZnstand[[34]] <- lm(logZnstand ~ (heatsubAADTlog100frt) + heatbing1902log300proj + nlcd_imp_ps +
                           heatbustransitlog200, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[34]])
#ols_plot_diagnostics(modlistlogZnstand[[34]])
ols_coll_diag(modlistlogZnstand[[34]])
ols_correlations(modlistlogZnstand[[34]])

modlistlogZnstand[[35]] <- lm(logZnstand ~ (heatsubAADTlog100frt) + heatbing1902log300proj + nlcd_imp_ps_mean, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[35]])
#ols_plot_diagnostics(modlistlogZnstand[[35]])
ols_coll_diag(modlistlogZnstand[[35]])
ols_correlations(modlistlogZnstand[[35]])

modlistlogZnstand[[36]] <- lm(logZnstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean, 
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[36]])
#ols_plot_diagnostics(modlistlogZnstand[[36]])
ols_coll_diag(modlistlogZnstand[[36]])
ols_correlations(modlistlogZnstand[[36]])

modlistlogZnstand[[37]] <- lm(logZnstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + heatbing1902log500proj,
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[37]])
#ols_plot_diagnostics(modlistlogZnstand[[37]])
ols_coll_diag(modlistlogZnstand[[37]])
ols_correlations(modlistlogZnstand[[37]])

modlistlogZnstand[[38]] <- lm(logZnstand ~ heatsubAADTlog200frt + nlcd_imp_ps_mean + heatbustransitlog200,
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[38]])
#ols_plot_diagnostics(modlistlogZnstand[[38]])
ols_coll_diag(modlistlogZnstand[[38]])
ols_correlations(modlistlogZnstand[[38]])

modlistlogZnstand[[39]] <- lm(logZnstand ~ heatsubAADTlog200frt + nlcd_imp_ps_mean + heatsubspdlpow300_3,
                         data = pollutfieldclean_cast)
ols_regress(modlistlogZnstand[[39]])
#ols_plot_diagnostics(modlistlogZnstand[[39]])
ols_coll_diag(modlistlogZnstand[[39]])
ols_correlations(modlistlogZnstand[[39]])

#------ 9. logZnstand - Make latex model summary table ----
vnum <- max(sapply(modlistlogZnstand, function(mod) {length(mod$coefficients)}))
model_summarylogZnstand<- as.data.table(
  ldply(modlistlogZnstand, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summarylogZnstand, -R2pred, AICc)  
cat(latex_format(model_summarylogZnstand), file = file.path(moddir, 'modeltable_logZnstand_2019.tex'))
setwd(moddir)

#Format table
modeltable_logZnstand_edit <- latex_format(model_summarylogZnstand) 
modeltable_logZnstand_edit <- gsub("pow100\\_1", "linear", modeltable_logZnstand_edit, fixed=T) 
modeltable_logZnstand_edit <- gsub("spdl", "speed", modeltable_logZnstand_edit) 
modeltable_logZnstand_edit <- gsub("nlcd\\_imp\\_ps", "imperviousness", modeltable_logZnstand_edit, fixed=T) 
modeltable_logZnstand_edit <- gsub("mean", "smooth", modeltable_logZnstand_edit, fixed=T) 
modeltable_logZnstand_edit <- gsub("bing1902", "congestion", modeltable_logZnstand_edit, fixed=T) 
modeltable_logZnstand_edit <- gsub("bustransit", "transit", modeltable_logZnstand_edit, fixed=T) 
modeltable_logZnstand_edit <- gsub("sub", "", modeltable_logZnstand_edit, fixed=T) 
modeltable_logZnstand_edit <- gsub("proj", "", modeltable_logZnstand_edit, fixed=T) 
modeltable_logZnstand_edit <- gsub("paperheight[=]10in[,]paperwidth[=]35in", "paperheight=8in,paperwidth=25in",
                                              modeltable_logZnstand_edit) 

cat(modeltable_logZnstand_edit, file = file.path(moddir, 'modeltable_logZnstand_2019.tex'))
texi2pdf('modeltable_logZnstand_2019.tex')

#------ 10. logZnstand - Make latex model summary table when excluding outliers----
vnum <- max(sapply(modlistlogZnstand, function(mod) {length(mod$coefficients)}))
model_summarylogZnstand_nooutliers<- as.data.table(
  ldply(modlistlogZnstand, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50,
                                                             remove_outliers = 'outliers & leverage',
                                                             labelvec = pollutfieldclean_cast[, SiteIDPair])}))
setorder(model_summarylogZnstand_nooutliers, -R2pred, AICc)  

#Format table
modeltable_logZnstand_nooutliers_edit <- latex_format(model_summarylogZnstand_nooutliers) 
modeltable_logZnstand_nooutliers_edit <- gsub("pow100\\_1", "linear", modeltable_logZnstand_nooutliers_edit, fixed=T) 
modeltable_logZnstand_nooutliers_edit <- gsub("spdl", "speed", modeltable_logZnstand_nooutliers_edit) 
modeltable_logZnstand_nooutliers_edit <- gsub("nlcd\\_imp\\_ps", "imperviousness", modeltable_logZnstand_nooutliers_edit, fixed=T) 
modeltable_logZnstand_nooutliers_edit <- gsub("mean", "smooth", modeltable_logZnstand_nooutliers_edit, fixed=T) 
modeltable_logZnstand_nooutliers_edit <- gsub("bing1902", "congestion", modeltable_logZnstand_nooutliers_edit, fixed=T) 
modeltable_logZnstand_nooutliers_edit <- gsub("bustransit", "transit", modeltable_logZnstand_nooutliers_edit, fixed=T) 
modeltable_logZnstand_nooutliers_edit <- gsub("sub", "", modeltable_logZnstand_nooutliers_edit, fixed=T) 
modeltable_logZnstand_nooutliers_edit <- gsub("proj", "", modeltable_logZnstand_nooutliers_edit, fixed=T) 
modeltable_logZnstand_nooutliers_edit <- gsub("paperheight[=]10in[,]paperwidth[=]35in", "paperheight=8in,paperwidth=27in",
                                              modeltable_logZnstand_nooutliers_edit) 

cat(modeltable_logZnstand_nooutliers_edit, file = file.path(moddir, 'modeltable_logZnstand_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_logZnstand_nooutliers_2019.tex')

#------ 11. logZnstand - test selected model ----
summary(modlistlogZnstand[[36]])
ols_plot_diagnostics(modlistlogZnstand[[36]])
regdiagnostic_customtab(mod=modlistlogZnstand[[36]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistlogZnstand[[36]]))
MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistlogZnstand[[36]]))
RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistlogZnstand[[36]]))

qplot(pollutfieldclean_cast[, heatbustransitlog200],modlistlogZnstand[[36]]$residuals) + 
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[, heatbing1902log500proj], modlistlogZnstand[[36]]$residuals) + 
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[, NLCD_reclass_final_PS], modlistlogZnstand[[36]]$residuals) + 
  geom_smooth(method='lm', color='red')


summary(modlistlogZnstand[[38]])
GAMrescheck(modlistlogZnstand[[38]])
regdiagnostic_customtab(mod=modlistlogZnstand[[38]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistlogZnstand[[38]]))
MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistlogZnstand[[38]]))
RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistlogZnstand[[38]]))

qplot(abs(fitted(modlistlogZnstand[[36]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand]), 
      abs(fitted(modlistlogZnstand[[38]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model 36 - aadt10frt + nlcd', 
       y='Absolute error for model 38 - aadt100frt + nlcd + transit200')

#------ 12. GLM Znstand - run all models for table ----
modlistglmZnstand <- list()

modlistglmZnstand[[1]] <- glm(Znstand ~ 1,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[2]] <- glm(Znstand ~ heatsubAADTlog200,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))
modlistglmZnstand[[3]] <- glm(Znstand ~ heatbing1902log300proj,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))
modlistglmZnstand[[4]] <- glm(Znstand ~ nlcd_imp_ps_mean,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[5]] <- glm(Znstand ~ nlcd_imp_ps,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))
modlistglmZnstand[[6]] <- glm(Znstand ~ heatbustransitlog200,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[7]] <- glm(Znstand ~ heatsubspdlpow100_3,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[8]] <- glm(Znstand ~ heatsubslopepow500_3,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[9]] <- glm(Znstand ~ heatsubspdlpow100_3,
                         data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[10]] <- glm(Znstand ~ heatsubspdlpow100_3 + heatsubAADTlog200,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[11]] <- glm(Znstand ~ heatsubspdlpow100_3+ heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[12]] <- glm(Znstand ~ heatsubspdlpow100_3 + nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[13]] <- glm(Znstand ~ heatsubspdlpow100_3 + heatbustransitlog200,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[14]] <- glm(Znstand ~ heatsubspdlpow100_3 + heatsubslopepow500_1,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[15]] <- glm(Znstand ~ heatsubspdlpow100_3*heatsubAADTlog200,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))


modlistglmZnstand[[16]] <- glm(Znstand ~ heatsubspdlpow100_3*heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[17]] <- glm(Znstand ~ heatsubspdlpow100_3*heatbing1902log300proj + nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[18]] <- glm(Znstand ~ heatsubspdlpow100_3*nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[19]] <- glm(Znstand ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[20]] <- glm(Znstand ~ heatsubspdlpow100_3 + heatbing1902log300proj*nlcd_imp_ps_mean + heatbustransitlog200,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[21]] <- glm(Znstand ~ heatsubspdlpow100_3*heatbustransitlog200 + heatbing1902log300proj*nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))


modlistglmZnstand[[22]] <- glm(Znstand ~ heatsubspdlpow100_3*heatbustransitlog200 + heatbing1902log300proj*nlcd_imp_ps,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[23]] <- glm(Znstand ~ heatsubAADTlog200frt,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[24]] <- glm(Znstand ~ heatsubAADTlog100frt,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[25]] <- glm(Znstand ~ heatsubAADTlog300frt,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))


modlistglmZnstand[[26]] <- glm(Znstand ~ heatsubAADTpow100_1frt,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[27]] <- glm(Znstand ~ heatsubAADTlog100frt + heatsubspdlpow100_3,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[28]] <- glm(Znstand ~ heatsubAADTlog100frt*heatsubspdlpow100_3,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[29]] <- glm(Znstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[29]] <- glm(Znstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[30]] <- glm(Znstand ~ heatsubAADTlog100frt*nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[31]] <- glm(Znstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[32]] <- glm(Znstand ~ heatsubAADTlog100frt*nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[33]] <- glm(Znstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean*heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[34]] <- glm(Znstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + heatbing1902log300proj + heatbustransitlog200,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[35]] <- glm(Znstand ~ heatsubAADTlog100frt*heatsubslopepow500_1 + nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[36]] <- glm(Znstand ~ heatsubAADTlog100frt*heatbustransitlog200 + nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[37]] <- glm(Znstand ~ heatsubAADTlog100frt*heatbustransitlog200 + heatsubspdllog100 + 
                            nlcd_imp_ps_mean + heatbing1902log300proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[38]] <- glm(Znstand ~ heatsubAADTlog100frt + heatbustransitlog200 + nlcd_imp_ps_mean,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[39]] <- glm(Znstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + heatbing1902log500proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

modlistglmZnstand[[40]] <- glm(Znstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean*heatbing1902log500proj,
                          data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))

# modlistglmZnstand[[41]] <- glm(Znstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + NLCD_reclass_final_PS_edit,
#                           data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))
# 
# modlistglmZnstand[[42]] <- glm(Znstand ~ heatsubAADTlog100frt + heatbing1902log500proj + NLCD_reclass_final_PS_edit,
#                           data = pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),], family=Gamma(link='log'))


#------ 13. GLM Znstand - Make latex model summary table ----
vnum <- max(sapply(modlistglmZnstand, function(mod) {length(mod$coefficients)}))

model_summaryglmZnstand<- as.data.table(
  ldply(modlistglmZnstand, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
model_summaryglmZnstand[, AICc := as.numeric(AICc)]
model_summaryglmZnstand[nvars==0, `:=`(
  MAEcv = as.character(round(
    MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistglmZnstand[[1]])),2)),
  RMSEcv = as.character(round(
    RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistglmZnstand[[1]])),2))
)]
setorder(model_summaryglmZnstand, AICc)  
cat(latex_format(model_summaryglmZnstand[, -c("BreuschPagan\\_fitp", "Score\\_fitp", "R2", "R2adj", "R2pred", "RMSE", "MAE", "VIF8")]),
    file = file.path(moddir, 'modeltable_glmZnstand_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_glmZnstand_2019.tex')

#------ 14. GLM Znstand - test selected model ----
summary(modlistglmZnstand[[29]])
GAMrescheck(modlistglmZnstand[[29]])
regdiagnostic_customtab(mod=modlistglmZnstand[[29]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistglmZnstand[[29]]))
MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistglmZnstand[[29]]))
RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistglmZnstand[[29]]))

qplot(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), heatbustransitlog200],
      modlistglmZnstand[[29]]$residuals) + 
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), heatbing1902log500proj],
      modlistglmZnstand[[29]]$residuals) + 
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), NLCD_reclass_final_PS],
      modlistglmZnstand[[29]]$residuals) + 
  geom_smooth(method='lm', color='red')


summary(modlistglmZnstand[[38]])
GAMrescheck(modlistglmZnstand[[38]])
regdiagnostic_customtab(mod=modlistglmZnstand[[38]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistglmZnstand[[38]]))
MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistglmZnstand[[38]]))
RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand], fitted(modlistglmZnstand[[38]]))

qplot(abs(fitted(modlistglmZnstand[[29]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand]), 
      abs(fitted(modlistglmZnstand[[38]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model 29 - aadt10frt + nlcd', 
       y='Absolute error for model 38 - aadt100frt + nlcd + transit200')


GAMrescheck(modlistglmZnstand[[39]])
qplot(abs(fitted(modlistglmZnstand[[29]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand]), 
      abs(fitted(modlistglmZnstand[[39]])-pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Znstand])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model 29 - aadt100frt + nlcd', 
       y='Absolute error for model 39 - aadt100frt + nlcd + bing500')

# plot_model(modlistglmZnstand[[33]], type='int', mdrt.values = 'meansd') +
#   scale_x_continuous(name = 'Imperviousness (%)', expand = c(0,0)) +
#   scale_y_continuous(name = 'Znstand index', expand=c(0,0), 
#                      breaks = c(0,0.25, 0.5, 0.75, 1.0), labels = c(0, 25, 50, 75, 100)) +
#   coord_cartesian(xlim=c(0,100), ylim=c(0,1)) +
#   theme_classic() +
#   theme(plot.title=element_blank(),
#         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
#         legend.position = c(0.8, 0.2), 
#         text =element_text(size=14))


#plot_model(modlistglmZnstand[[33]], type='pred', mdrt.values = 'meansd')

#------ 15. GAM Znstand - Single parameter models -------
modlistGAMZnstand <- list() #List to hold models
modlistGAMZnstand[[1]] <- lm(logZnstand ~ 1, data = pollutfieldclean_cast) #Null/Intercept model

modlistGAMZnstand[[2]] <- mgcv::gam(Znstand~s(heatsubAADTlog100, k=4), family=gaussian(link='log'), method='ML',
                               data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[2]])
plot(modlistGAMZnstand[[2]],residuals=TRUE,shade=T, cex=6)

pollutfieldclean_cast$predgamZnstand <- predict(modlistGAMZnstand[[2]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog200, y=log(Znstand), label=paste0(SiteID, Pair))) +
  geom_text() + 
  geom_line(aes(y=predgamZnstand))

modlistGAMZnstand[[3]] <- mgcv::gam(Znstand~s(heatsubAADTlog100, k=4), family=Gamma(link='log'), method='ML',
                               data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[3]])
plot(modlistGAMZnstand[[3]],residuals=TRUE,shade=T, cex=6)

pollutfieldclean_cast$predgamZnstand <- predict(modlistGAMZnstand[[3]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog200, y=log(Znstand), label=paste0(SiteID, Pair))) +
  geom_text() + 
  geom_line(aes(y=predgamZnstand), color='red', size=1.2) +
  theme_classic()

modlistGAMZnstand[[4]] <- mgcv::gam(Znstand~s(heatsubAADTlog200, k=4), family=gaussian(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[4]])
plot(modlistGAMZnstand[[4]],residuals=TRUE,shade=T, cex=6)

modlistGAMZnstand[[5]] <- mgcv::gam(Znstand~s(heatsubAADTlog200, k=4), family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[5]])
plot(modlistGAMZnstand[[5]],residuals=TRUE,shade=T, cex=6)

modlistGAMZnstand[[6]] <- mgcv::gam(Znstand~s(heatsubspdlpow100_3, k=4), family=gaussian(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[6]])
plot(modlistGAMZnstand[[6]],residuals=TRUE,shade=T, cex=6)

modlistGAMZnstand[[7]] <- mgcv::gam(Znstand~s(heatbing1902log300proj, k=4), 
                               family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[7]])
plot(modlistGAMZnstand[[7]],residuals=TRUE,shade=T, cex=6)

modlistGAMZnstand[[8]] <- mgcv::gam(Znstand~s(heatbustransitlog200, k=4), 
                               family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[8]])
plot(modlistGAMZnstand[[8]],residuals=TRUE,shade=T, cex=6)

modlistGAMZnstand[[9]] <- mgcv::gam(Znstand~s(nlcd_imp_ps, k=3), 
                               family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[9]])
plot(modlistGAMZnstand[[9]],residuals=TRUE,shade=T, cex=6)

modlistGAMZnstand[[10]] <- mgcv::gam(Znstand~s(nlcd_imp_ps_mean, k=3), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[10]])
plot(modlistGAMZnstand[[10]],residuals=TRUE,shade=T, cex=6)

modlistGAMZnstand[[11]] <- mgcv::gam(Znstand~s(heatsubslopepow500_1, k=4), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[11]])
plot(modlistGAMZnstand[[11]],residuals=TRUE,shade=T, cex=6)

modlistGAMZnstand[[12]] <- mgcv::gam(Znstand~s(heatbing1902log500proj, k=4), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[12]])
plot(modlistGAMZnstand[[12]],residuals=TRUE,shade=T, cex=6)

#------ 16. GAM Znstand - Multiparameter models -------
modlistGAMZnstand[[11]] <- mgcv::gam(Znstand~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3),
                                family=Gamma(link='log'), 
                                data=pollutfieldclean_cast,
                                select=TRUE)
GAMrescheck(modlistGAMZnstand[[11]])
GAMmultiplot(modlistGAMZnstand[[11]])
concurvity(modlistGAMZnstand[[11]])

pollutfieldclean_cast[, predgamZnstand := fitted(modlistGAMZnstand[[11]], type='response')]
ggplot(pollutfieldclean_cast, aes(x=predgamZnstand, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text() +
  geom_abline(intercept=0, slope=1) +
  #scale_x_continuous(limits=c(0,2)) +
  coord_fixed()

modlistGAMZnstand[[12]] <- mgcv::gam(Znstand~s(heatsubAADTlog100thd) + s(heatsubAADTlog100thd, k=4, by=nlcd_imp_ps_mean) + 
                                  s(nlcd_imp_ps, k=3), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[12]])
GAMmultiplot(modlistGAMZnstand[[12]])


modlistGAMZnstand[[13]] <- mgcv::gam(Znstand~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbing1902log300proj, k=4), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[13]])
GAMmultiplot(modlistGAMZnstand[[13]])
concurvity(modlistGAMZnstand[[13]])

modlistGAMZnstand[[14]] <- mgcv::gam(Znstand~s(heatsubAADTlog100,k=4) + s(heatsubAADTlog100, k=4, by=heatbing1902log500proj) + 
                                  s(nlcd_imp_ps, k=3), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[14]])
GAMmultiplot(modlistGAMZnstand[[14]])

modlistGAMZnstand[[15]] <- mgcv::gam(Znstand~s(heatsubAADTlog100thd,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[15]])
GAMmultiplot(modlistGAMZnstand[[15]])
#Plot
#Plot

modlistGAMZnstand[[16]] <- mgcv::gam(Znstand~s(heatsubAADTlog100, k=4) + s(heatsubspdlpow100_3, k=4)+ 
                                  s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[16]])
GAMmultiplot(modlistGAMZnstand[[16]])

modlistGAMZnstand[[17]] <- mgcv::gam(Znstand~s(heatsubAADTlog100, k=4) + s(heatsubAADTlog100, k=4, by=heatsubspdlpow100_3)+ 
                                  s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[17]])
GAMmultiplot(modlistGAMZnstand[[17]])
concurvity(modlistGAMZnstand[[17]])

pollutfieldclean_cast[, predgamZnstand := fitted(modlistGAMZnstand[[17]], type='response')]
ggplot(pollutfieldclean_cast, aes(x=predgamZnstand, y=Znstand, label=paste0(SiteID, Pair))) + 
  geom_text() +
  geom_abline(intercept=0, slope=1) +
  scale_x_continuous(limits=c(0,2)) +
  coord_fixed()

modlistGAMZnstand[[18]] <- mgcv::gam(Znstand~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200) +
                                  s(heatsubAADTlog100, k=4, by=heatsubslopepow500_1), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[18]])
GAMmultiplot(modlistGAMZnstand[[18]])
concurvity(modlistGAMZnstand[[18]])

modlistGAMZnstand[[19]] <- mgcv::gam(Znstand~s(heatsubAADTlog100thd,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200, k=4) +
                                  s(heatbing1902log500proj, k=3), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[19]])
GAMmultiplot(modlistGAMZnstand[[19]])
concurvity(modlistGAMZnstand[[19]])

modlistGAMZnstand[[20]] <- mgcv::gam(Znstand~s(heatsubAADTlog100thd,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog300, k=4) +
                                  s(heatbing1902log500proj, k=3), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast[])
GAMrescheck(modlistGAMZnstand[[20]])
GAMmultiplot(modlistGAMZnstand[[20]])
concurvity(modlistGAMZnstand[[20]])

modlistGAMZnstand[[21]] <- mgcv::gam(Znstand~s(heatsubspdllog100,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                                family=Gamma(link='log'), data=pollutfieldclean_cast)
GAMrescheck(modlistGAMZnstand[[21]])
GAMmultiplot(modlistGAMZnstand[[21]])
concurvity(modlistGAMZnstand[[21]])

#------ 17. GAM Znstand - Multiparameter models without Luwam data or 15th Ave outliers -------
GAMZnstand_sub11 <- mgcv::gam(Znstand~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3), 
                         family=Gamma(link='log'), data=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),])
GAMrescheck(GAMZnstand_sub11)
GAMmultiplot(GAMZnstand_sub11)
concurvity(GAMZnstand_sub11)

GAMZnstand_sub13 <- mgcv::gam(Znstand~s(heatsubAADTlog100, k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbing1902log300proj, k=4), 
                         family=Gamma(link='log'), data=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),])
GAMrescheck(GAMZnstand_sub13)
GAMmultiplot(GAMZnstand_sub13)
concurvity(GAMZnstand_sub13)

GAMZnstand_sub20 <- mgcv::gam(Znstand~s(heatsubAADTlog100,k=4) + s(nlcd_imp_ps_mean, k=3) + s(heatbustransitlog200), 
                         family=Gamma(link='log'),data=pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers),])
GAMrescheck(GAMZnstand_sub20)
GAMmultiplot(GAMZnstand_sub20)

#------ 18. GAM Znstand - Multiparameter models without outliers based on influence-------
#Model 11
influencemod11 <- GAMinfluence(modlistGAMZnstand[[11]], pollutfieldclean_cast)

ggplot(pollutfieldclean_cast, aes(x=fitted(modlistGAMZnstand[[11]]), y=Znstand, 
                                  color=influencemod11$influence, label=SiteIDPair)) + 
  geom_text() +
  scale_color_distiller(palette='Spectral')
ggplot(influencemod11, aes(x=influence, label=SiteIDPair)) + 
  geom_histogram()
pollutfieldclean_cast[!(SiteIDPair %in% 
                          influencemod11[influencemod11$influence>0.04, 'SiteIDPair'])]

GAMZnstand_sub11 <- mgcv::gam(Znstand~s(heatsubAADTlog100,k=4) + s(nlcd_imp_ps_mean, k=3), 
                         family=Gamma(link='log'), 
                         data=pollutfieldclean_cast[!(SiteIDPair %in% 
                                                        influencemod11[influencemod11$influence>0.04, 'SiteIDPair']),])
GAMrescheck(GAMZnstand_sub11)
GAMmultiplot(GAMZnstand_sub11)
GAMrescheck(modlistGAMZnstand[[11]])

#-------19. GLM Znstand - Compare final selected models -------
mod29_nooutliers <- regdiagnostic_customtab(modlistglmZnstand[[29]], maxpar=vnum, 
                                            remove_outliers = 'outliers',
                                            labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                            kCV = TRUE, k=10, cvreps=50)
#Determined by deleted studentized residuals
subdat29 <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                                      strsplit(gsub('\\\\', '', mod29_nooutliers['outliers']), ',')$outliers),]
#Simply excluding 15th Ave bus stop and Luwam's data
subdat29_2 <- pollutfieldclean_cast[!(paste0(SiteIDPair) %in% extraoutliers),]

mod29_nooutliersub <- glm(Znstand ~  heatsubAADTlog100frt + nlcd_imp_ps_mean, 
                          data = subdat29_2, family=Gamma('log'))
GAMrescheck(mod29_nooutliersub)
plot(mod29_nooutliersub)
#Plot
#Plot
#Plot
#Plot
MAPE(subdat29_2[, Znstand], fitted(mod29_nooutliersub))
MAE(subdat29_2[, Znstand], fitted(mod29_nooutliersub))
RMSE(subdat29_2[, Znstand], fitted(mod29_nooutliersub))

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
Znstandmod29plot <- ggplot(data=pollutfieldclean_cast[!(SiteIDPair %in% c('64A', '65A', '66A', '67A', '68A', '69A')),],
                      aes(x=mod29pred$fit,  y=Znstand)) +
  geom_point(alpha=0.75, size=2, color='red') + 
  geom_ribbon(aes(ymin=mod29pred$fit + 1.96*mod29pred$se.fit, ymax=mod29pred$fit - 1.96*mod29pred$se.fit), fill='orange', alpha=1/4) + 
  geom_point(data=subdat29_2, aes(x=fitted(mod29_nooutliersub), y=Znstand), color='black', alpha=1/2, size=2) +
  geom_abline(intercept=0, slope=1, size=1.3, color='red') +
  #geom_smooth(method='lm', se=FALSE) +
  #geom_text(aes(label=paste0(SiteID,Pair))) + 
  scale_x_continuous(limits=c(0,2.1), expand=c(0,0), breaks=seq(0,2,0.5)) +
  scale_y_continuous(limits=c(0,2.1), expand=c(0,0), breaks=seq(0,2,0.5)) +
  coord_fixed() +
  labs(x='Predicted Znstand index', y='Observed Znstand index') +
  theme_classic() + 
  theme(text= element_text(size=20))
png(file.path(moddir, 'scatterplot_Znstand_mod29.png'), width=9, height=9, units='in', res=300)
Znstandmod29plot
dev.off()

#Compare glm 29 with logZnstand 36
#Compare with predictions without transformation

qplot(abs(fitted(modlistglmZnstand[[29]])-pollutfieldclean_cast[!(paste0(SiteIDPair) %in% extraoutliers), Znstand]), 
      abs(exp(predict(modlistlogZnstand[[36]], newdata=pollutfieldclean_cast[!(paste0(SiteIDPair) %in% extraoutliers)]))-
            pollutfieldclean_cast[!(paste0(SiteIDPair) %in% extraoutliers), Znstand])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model glm 29 - glm  nlcd_imp_ps_mean + AADTlog100frt ', 
       y='Absolute error for model log 36 - nlcd_imp_ps_mean + AADTlog100frt')

#------ 20. Znstand - Check spatial and temporal autocorrelation of residuals for full dataset -------
#Other good resource: https://eburchfield.github.io/files/Spatial_regression_LAB.html
resnorm <- rstandard(modlistlogZnstand[[36]]) #Get standardized residuals from model
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
lm.morantest(modlistlogZnstand[[36]], listw = listw2U(weightmat_k[[1]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZnstand[[36]], listw = listw2U(weightmat_k[[2]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZnstand[[36]], listw = listw2U(weightmat_k[[3]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZnstand[[36]], listw = listw2U(weightmat_k[[4]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZnstand[[36]], listw = listw2U(weightmat_k[[5]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(modlistlogZnstand[[36]], listw = listw2U(weightmat_all)) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)

#Test for need for spatial regression model using Lagrange Multiplier (LM) tests
lm.LMtests(modlistlogZnstand[[36]], listw = listw2U(weightmat_k[[1]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(modlistlogZnstand[[36]], listw = listw2U(weightmat_k[[2]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(modlistlogZnstand[[36]], listw = listw2U(weightmat_k[[3]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))

#Spatial simultaneous autoregressive error model estimation with 1 nearest neighbors
sarlm_modlogZnstand <- errorsarlm(modlistlogZnstand[[36]]$call$formula, data = modlistlogZnstand[[36]]$model, 
                             listw = listw2U(weightmat_k[[1]]))
summary(sarlm_modlogZnstand)
bptest.sarlm(sarlm_modlogZnstand)
cor(modlistlogZnstand[[36]]$model$logZnstand, fitted(sarlm_modlogZnstand))^2

#Compare pseudo-R2
cor(pollutfieldclean_cast$Znstand, exp(predict(sarlm_modlogZnstand, pred.type='trend')))^2
cor(pollutfieldclean_cast$Znstand,
    with(pollutfieldclean_cast, exp(sarlm_modlogZnstand$coefficients['(Intercept)'] + 
                                      sarlm_modlogZnstand$coefficients['heatsubAADTlog100frt']*heatsubAADTlog100frt + 
                                      sarlm_modlogZnstand$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean)))^2

#Compare MAE
DescTools::MAPE(exp(modlistlogZnstand[[36]]$model$logZnstand), exp(predict(sarlm_modlogZnstand, pred.type='trend')))
DescTools::MAPE(exp(pollutfieldclean_cast$logZnstand), 
                with(pollutfieldclean_cast, exp(sarlm_modlogZnstand$coefficients['(Intercept)'] + 
                                                  sarlm_modlogZnstand$coefficients['heatsubAADTlog100frt']*heatsubAADTlog100frt + 
                                                  sarlm_modlogZnstand$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean)  
                       ))

#Compare observed~predicted for full-no outlier model and for aspatial and spatial model
spatial_comparisonplot <- ggplot(pollutfieldclean_cast, aes(x=exp(fitted(sarlm_modlogZnstand)), y=exp(logZnstand))) + 
  geom_point(aes(x=exp(fitted(modlistlogZnstand[[36]]))), size=2, alpha=1/2, color='orange') +
  geom_point(size=2, alpha=1/2, color='red') + 
  geom_abline(size=1, slope=1, intercept=0, color='red') + 
  #geom_text(aes(label=paste0(SiteID, Pair))) +
  coord_fixed() +
  theme_classic()
spatial_comparisonplot

resnorm_postsarlm <- residuals(sarlm_modlogZnstand) #Get standardized residuals from model
#Make bubble map of residuals
bubbledat_postsarlm <- data.frame(resnorm_postsarlm, pollutfieldclean_cast$coords.x1, pollutfieldclean_cast$coords.x2)
coordinates(bubbledat_postsarlm) <- c("pollutfieldclean_cast.coords.x1","pollutfieldclean_cast.coords.x2")
bubble(bubbledat_postsarlm, "resnorm_postsarlm", col = c("blue","red"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")

### FINAL VERDICT: GO WITH SARLM COEFFICIENTS FOR MODLOG36-----

#--------------- B. Custand ~ separate predictors ----
modlistCustand <- list() #List to hold models
modlistCustand[[1]] <- lm(Custand ~ 1, data = pollutfieldclean_cast) #Null/Intercept model
#------ 1. Custand - Single parameter models --------
modlistCustand[[2]] <- lm(Custand ~ heatsubAADTlog200, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[2]])
#ols_plot_diagnostics(modlistCustand[[2]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog200, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[3]] <- lm(Custand ~ heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[3]])
#ols_plot_diagnostics(modlistCustand[[3]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902log300proj, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[4]] <- lm(Custand ~ heatbing1902log500proj, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[4]])
#ols_plot_diagnostics(modlistCustand[[4]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902log500proj, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[5]] <- lm(Custand ~ nlcd_imp_ps, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[5]])
#ols_plot_diagnostics(modlistCustand[[5]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[6]] <- lm(Custand ~ nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[6]])
#ols_plot_diagnostics(modlistCustand[[6]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps_mean, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[7]] <- lm(Custand ~ heatbustransitlog200, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[7]])
#ols_plot_diagnostics(modlistCustand[[7]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitlog200, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[8]] <- lm(Custand ~ heatbustransitlog200sqrt, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[8]])
#ols_plot_diagnostics(modlistCustand[[8]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitlog200sqrt, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[9]] <- lm(Custand ~ heatsubspdlpow50_1, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[9]])
#ols_plot_diagnostics(modlistCustand[[9]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow50_1, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[10]] <- lm(Custand ~ heatsubslopepow500_1, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[10]])
#ols_plot_diagnostics(modlistCustand[[10]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow500_1, y=Custand, label=SiteIDPair)) + 
  geom_text()


modlistCustand[[11]] <- lm(Custand ~ heatsubslopepow500_1, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[11]])
#ols_plot_diagnostics(modlistCustand[[11]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow500_1, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[12]] <- lm(Custand ~ heatsubAADTlog100, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[12]])
#ols_plot_diagnostics(modlistCustand[[12]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog100, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[13]] <- lm(Custand ~ heatbustransitlog300sqrt, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[13]])
#ols_plot_diagnostics(modlistCustand[[13]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitlog300sqrt, y=Custand, label=SiteIDPair)) + 
  geom_text()

modlistCustand[[14]] <- lm(Custand ~ heatsubslopepow300_1, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[14]])
#ols_plot_diagnostics(modlistCustand[[14]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow300_1, y=Custand, label=SiteIDPair)) + 
  geom_text()


modlistCustand[[15]] <- lm(Custand ~ heatsubslopepow300_1, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[15]])
#ols_plot_diagnostics(modlistCustand[[15]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow300_1, y=Custand, label=SiteIDPair)) + 
  geom_text()

#------ 2. Custand - Multiparameter models --------
modlistCustand[[16]] <- lm(Custand ~ heatbustransitlog300sqrt + heatbing1902log500proj, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[16]])
#ols_plot_diagnostics(modlistCustand[[16]])
ols_coll_diag(modlistCustand[[16]])
ols_correlations(modlistCustand[[16]])

modlistCustand[[17]] <- lm(Custand ~ heatbustransitlog200sqrt + heatbing1902log500proj, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[17]])
#ols_plot_diagnostics(modlistCustand[[17]])
ols_coll_diag(modlistCustand[[17]])
ols_correlations(modlistCustand[[17]])

modlistCustand[[18]] <- lm(Custand ~ heatbustransitlog200sqrt + heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[18]])
#ols_plot_diagnostics(modlistCustand[[18]])
ols_coll_diag(modlistCustand[[18]])
ols_correlations(modlistCustand[[18]])

modlistCustand[[19]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[19]])
#ols_plot_diagnostics(modlistCustand[[19]])
ols_coll_diag(modlistCustand[[19]])
ols_correlations(modlistCustand[[19]])

modlistCustand[[20]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[20]])
#ols_plot_diagnostics(modlistCustand[[20]])
ols_coll_diag(modlistCustand[[20]])
ols_correlations(modlistCustand[[20]])

modlistCustand[[21]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubslopepow300_1, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[21]])
#ols_plot_diagnostics(modlistCustand[[21]])
ols_coll_diag(modlistCustand[[21]])
ols_correlations(modlistCustand[[21]])

modlistCustand[[22]] <- lm(Custand ~ heatsubspdlpow50_1 + nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[22]])
#ols_plot_diagnostics(modlistCustand[[22]])
ols_coll_diag(modlistCustand[[22]])
ols_correlations(modlistCustand[[22]])

modlistCustand[[23]] <- lm(Custand ~ heatsubspdlpow50_1 + heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[23]])
#ols_plot_diagnostics(modlistCustand[[23]])
ols_coll_diag(modlistCustand[[23]])
ols_correlations(modlistCustand[[23]])

modlistCustand[[24]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubspdlpow50_1, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[24]])
#ols_plot_diagnostics(modlistCustand[[24]])
ols_coll_diag(modlistCustand[[24]])
ols_correlations(modlistCustand[[24]])

modlistCustand[[25]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubspdlpow50_1, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[25]])
#ols_plot_diagnostics(modlistCustand[[25]])
ols_coll_diag(modlistCustand[[25]])
ols_correlations(modlistCustand[[25]])

modlistCustand[[26]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatbing1902log500proj, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[26]])
#ols_plot_diagnostics(modlistCustand[[26]])
ols_coll_diag(modlistCustand[[26]])
ols_correlations(modlistCustand[[26]])

modlistCustand[[27]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean*heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[27]])
#ols_plot_diagnostics(modlistCustand[[27]])
ols_coll_diag(modlistCustand[[27]])
ols_correlations(modlistCustand[[27]])

modlistCustand[[28]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog100, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[28]])
#ols_plot_diagnostics(modlistCustand[[28]])
ols_coll_diag(modlistCustand[[28]])
ols_correlations(modlistCustand[[28]])

modlistCustand[[29]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog200, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[29]])
#ols_plot_diagnostics(modlistCustand[[29]])
ols_coll_diag(modlistCustand[[29]])
ols_correlations(modlistCustand[[29]])

modlistCustand[[30]] <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog500, data = pollutfieldclean_cast)
ols_regress(modlistCustand[[30]])
#ols_plot_diagnostics(modlistCustand[[30]])
ols_coll_diag(modlistCustand[[30]])
ols_correlations(modlistCustand[[30]])

#------ 3. Custand - Make latex model summary table ----
vnum <- max(sapply(modlistCustand, function(mod) {length(mod$coefficients)}))
model_summaryCustand<- as.data.table(
  ldply(modlistCustand, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50,
                                                          labelvec = pollutfieldclean_cast[, SiteIDPair])}))
model_summaryCustand[, AICc := as.numeric(AICc)]
setorder(model_summaryCustand, AICc, -R2pred)  
cat(latex_format(model_summaryCustand), file = file.path(moddir, 'modeltable_Custand.tex'))
setwd(moddir)
texi2pdf('modeltable_Custand.tex')

#------ 4. Custand - Make latex model summary table when excluding outliers ----
vnum <- max(sapply(modlistCustand, function(mod) {length(mod$coefficients)}))
model_summaryCustand_nooutliers<- as.data.table(
  ldply(modlistCustand, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50,
                                                          remove_outliers = 'outliers',
                                                          labelvec = pollutfieldclean_cast[, SiteIDPair])}))
setorder(model_summaryCustand_nooutliers, -R2pred, AICc)  

model_summaryCustand_nooutliers_edit <- latex_format(model_summaryCustand_nooutliers) 
model_summaryCustand_nooutliers_edit <- gsub("pow100\\_1", "linear", model_summaryCustand_nooutliers_edit, fixed=T) 
model_summaryCustand_nooutliers_edit <- gsub("spdl", "speed", model_summaryCustand_nooutliers_edit) 
model_summaryCustand_nooutliers_edit <- gsub("nlcd\\_imp\\_ps", "imperviousness", model_summaryCustand_nooutliers_edit, fixed=T) 
model_summaryCustand_nooutliers_edit <- gsub("mean", "smooth", model_summaryCustand_nooutliers_edit, fixed=T) 
model_summaryCustand_nooutliers_edit <- gsub("bing1902", "congestion", model_summaryCustand_nooutliers_edit, fixed=T) 
model_summaryCustand_nooutliers_edit <- gsub("bustransit", "transit", model_summaryCustand_nooutliers_edit, fixed=T) 
model_summaryCustand_nooutliers_edit <- gsub("sub", "", model_summaryCustand_nooutliers_edit, fixed=T) 
model_summaryCustand_nooutliers_edit <- gsub("proj", "", model_summaryCustand_nooutliers_edit, fixed=T) 
model_summaryCustand_nooutliers_edit <- gsub("paperheight[=]10in[,]paperwidth[=]35in", "paperheight=8in,paperwidth=27in",
                                              model_summaryCustand_nooutliers_edit) 

cat(model_summaryCustand_nooutliers_edit, file = file.path(moddir, 'modeltable_Custand_nooutliers.tex'))
setwd(moddir)
texi2pdf('modeltable_Custand_nooutliers.tex')

#-------5. Custand - Check final selected model -------
mod19_nooutliers <- regdiagnostic_customtab(modlistCustand[[19]], maxpar=vnum, 
                                            remove_outliers = 'outliers',
                                            labelvec = pollutfieldclean_cast[, SiteIDPair],
                                            kCV = TRUE, k=10, cvreps=50)
subdat <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                       strsplit(gsub('\\\\', '', mod19_nooutliers['outliers']), ',')$outliers),]
mod19_nooutliersub <- lm(Custand ~ heatbustransitlog200sqrt +  nlcd_imp_ps_mean, 
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

ggplot(pollutfieldclean_cast, aes(y=predict(modlistCustand[[19]], pollutfieldclean_cast), x=Custand)) +
  geom_point(alpha=1/4, size=2) + 
  #geom_point(data=subdat, aes(y = predict(modlistCustand[[19]], pollutfieldclean_cast)), alpha=1/2, size=2) + 
  geom_point(data=subdat, aes(x=Custand, y = predict(mod19_nooutliersub, subdat)), color='red', alpha=1/2, size=2) + 
  geom_abline(intercept=0, slope=1) +
  theme_classic()


#------ 6. GLM Custand - run all models for table ----
modlistglmCustand <- list()

modlistglmCustand[[1]] <- glm(Custand ~ heatsubAADTlog200, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[2]] <- glm(Custand ~ heatsubAADTlog200, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[2]] <- glm(Custand ~ heatsubAADTlog200, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[3]] <- glm(Custand ~ heatbing1902log300proj, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[4]] <- glm(Custand ~ heatbing1902log500proj, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[5]] <- glm(Custand ~ nlcd_imp_ps, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[6]] <- glm(Custand ~ nlcd_imp_ps_mean, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[7]] <- glm(Custand ~ heatbustransitlog200, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[8]] <- glm(Custand ~ heatbustransitlog200sqrt, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[9]] <- glm(Custand ~ heatsubspdlpow50_1, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[10]] <- glm(Custand ~ heatsubslopepow500_1, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[11]] <- glm(Custand ~ heatsubslopepow500_1, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[12]] <- glm(Custand ~ heatsubAADTlog100, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[13]] <- glm(Custand ~ heatbustransitlog300sqrt, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[14]] <- glm(Custand ~ heatsubslopepow300_1, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[15]] <- glm(Custand ~ heatsubslopepow300_1, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[16]] <- glm(Custand ~ heatbustransitlog300sqrt + heatbing1902log500proj, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[17]] <- glm(Custand ~ heatbustransitlog200sqrt + heatbing1902log500proj, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[18]] <- glm(Custand ~ heatbustransitlog200sqrt + heatbing1902log300proj, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[19]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[20]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatbing1902log300proj, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[21]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubslopepow300_1, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[22]] <- glm(Custand ~ heatsubspdlpow50_1 + nlcd_imp_ps_mean, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[23]] <- glm(Custand ~ heatsubspdlpow50_1 + heatbing1902log300proj, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[24]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubspdlpow50_1, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[25]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubspdlpow50_1, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[26]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatbing1902log500proj, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[27]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean*heatbing1902log300proj, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[28]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog100, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[29]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog200, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[30]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog300, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[31]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + heatsubAADTlog500, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[32]] <- glm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean*heatsubAADTlog500, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[33]] <- glm(Custand ~ heatbustransitlog200sqrt + heatsubAADTlog300, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[34]] <- glm(Custand ~ heatbustransitlog200sqrt + heatsubAADTlog500, data=pollutfieldclean_cast, family=Gamma(link='log'))

modlistglmCustand[[35]] <- glm(Custand ~ nlcd_imp_ps_mean + heatsubAADTlog300, data=pollutfieldclean_cast, family=Gamma(link='log'))

#------ 7. GLM Custand - Make latex model summary table ----
vnum <- max(sapply(modlistglmCustand, function(mod) {length(mod$coefficients)}))

model_summaryglmCustand<- as.data.table(
  ldply(modlistglmCustand, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50)}))
model_summaryglmCustand[, AICc := as.numeric(AICc)]
model_summaryglmCustand[nvars==0, `:=`(
  MAEcv = as.character(round(
    MAE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Custand], fitted(modlistglmCustand[[1]])),2)),
  RMSEcv = as.character(round(
    RMSE(pollutfieldclean_cast[!(SiteIDPair %in% extraoutliers), Custand], fitted(modlistglmCustand[[1]])),2))
)]
setorder(model_summaryglmCustand, AICc)  
cat(latex_format(model_summaryglmCustand[, -c("BreuschPagan\\_fitp", "Score\\_fitp", "R2", "R2adj", "R2pred", "RMSE", "MAE", "VIF8")]),
    file = file.path(moddir, 'modeltable_glmCustand_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_glmCustand_2019.tex')

#Format table by hand
#pow100\_1 to linear (and other kernel size)
#spdl to speed
#nlcd\_imp\_ps to imperviousness
#heatbing1902 to congestion
#bustransit to transit
#remove sub
#reduce age size to 25 in width and 8 in height

texi2pdf('modeltable_glmCustand_2019edit.tex')

#------ 8. GLM Custand - Make latex model summary table when excluding outliers ----
vnum <- max(sapply(modlistglmCustand, function(mod) {length(mod$coefficients)}))
model_summaryglmCustand_nooutliers<- as.data.table(
  ldply(modlistglmCustand, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, kCV = TRUE, k=10, cvreps=50,
                                                             remove_outliers = 'outliers', # & leverage',
                                                             labelvec = pollutfieldclean_cast[, SiteIDPair])}))
setorder(model_summaryglmCustand_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summaryglmCustand_nooutliers), file = file.path(moddir, 'modeltable_glmCustand_nooutliers.tex'))
setwd(moddir)
texi2pdf('modeltable_glmCustand_nooutliers.tex')

#------ 9. GLM Custand - test selected model ----
summary(modlistglmCustand[[30]])
GAMrescheck(modlistglmCustand[[30]])
regdiagnostic_customtab(mod=modlistglmCustand[[30]], maxpar=vnum, kCV = TRUE, k=10, cvreps=50)
MAPE(pollutfieldclean_cast[, Custand], fitted(modlistglmCustand[[30]]))
MAE(pollutfieldclean_cast[, Custand], fitted(modlistglmCustand[[30]]))
RMSE(pollutfieldclean_cast[, Custand], fitted(modlistglmCustand[[30]]))

qplot(pollutfieldclean_cast[, heatbing1902log500proj],
      modlistglmCustand[[30]]$residuals) +
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[, heatsubspdlpow50_1],
      modlistglmCustand[[30]]$residuals) +
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[, heatsubslopepow300_1],
      modlistglmCustand[[30]]$residuals) +
  geom_smooth(method='lm', color='red')

qplot(pollutfieldclean_cast[, NLCD_reclass_final_PS],
      modlistglmCustand[[30]]$residuals) +
  geom_smooth(method='lm', color='red')

#Compare with predictions without transformation
qplot(abs(fitted(modlistglmCustand[[30]])-pollutfieldclean_cast[, Znstand]), 
      abs(fitted(modlistCustand[[19]])-pollutfieldclean_cast[, Znstand])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model glm 29 - glm transit200 + nlcd + AADT300 ', 
       y='Absolute error for model 19 - transit200 + nlcd')

#------ 10. Custand - Check spatial and temporal autocorrelation of residuals for dataset without outliers -------
par(mfrow=c(1,1))
resnorm <- rstandard(mod19_nooutliersub) #Get standardized residuals from model without outliers (even when removing points that are outliers and leverage, can't get normal residuals)
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
sarlm_modCustand <- errorsarlm(mod19_nooutliersub$call$formula, data = mod19_nooutliersub$model, 
                          listw = listw2U(weightmat_k[[1]]))
summary(sarlm_modCustand)
bptest.sarlm(sarlm_modCustand)
qplot(fitted(sarlm_modCustand), sarlm_modCustand$residuals)
shapiro.test(sarlm_modCustand$residuals)
cor(mod19_nooutliersub$model$Custand, fitted(sarlm_modCustand))^2

#Compare pseudo-R2
cor(mod19_nooutliersub$model$Custand, predict(sarlm_modCustand, pred.type='trend'))^2
cor(mod19_nooutliersub$model$Custand, fitted(mod19_nooutliersub))^2
cor(pollutfieldclean_cast$Custand, with(pollutfieldclean_cast, sarlm_modCustand$coefficients['(Intercept)'] + 
                                             sarlm_modCustand$coefficients['heatbustransitlog200sqrt']*heatbustransitlog200sqrt + 
                                             sarlm_modCustand$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean))^2

#Compare MAPE
DescTools::MAPE(mod19_nooutliersub$model$Custand, fitted(sarlm_modCustand))
DescTools::MAPE(mod19_nooutliersub$model$Custand, predict(sarlm_modCustand, pred.type='trend'))
DescTools::MAPE(mod19_nooutliersub$model$Custand, fitted(mod19_nooutliersub))

#Compare MAPE with all sites
DescTools::MAPE(pollutfieldclean_cast$Custand,with(pollutfieldclean_cast, sarlm_modCustand$coefficients['(Intercept)'] + 
                                                     sarlm_modCustand$coefficients['heatbustransitlog200sqrt']*heatbustransitlog200sqrt + 
                                                     sarlm_modCustand$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean)) #reported in paper
DescTools::MAPE(pollutfieldclean_cast$Custand, predict(mod19_nooutliersub, pollutfieldclean_cast))


#Compare observed~predicted for full-no outlier model and for aspatial and spatial model
spatial_comparisonplot <- ggplot(subdat, aes(x=fitted(sarlm_modCustand), y=Custand)) + 
  geom_point(aes(x=fitted(mod19_nooutliersub)), size=2, alpha=1/2, color='orange') +
  geom_point(size=2, alpha=1/2, color='red') + 
  geom_abline(size=1, slope=1, intercept=0, color='red') + 
  #geom_text(aes(label=paste0(SiteID, Pair))) +
  coord_fixed() +
  theme_classic()
spatial_comparisonplot

resnorm_postsarlm <- residuals(sarlm_modCustand) #Get standardized residuals from model
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
modlistPI[[1]] <- lm(pollution_index ~ 1, data = pollutfieldclean_cast) #Null/Intercept model

modlistPI[[2]] <- lm(pollution_index ~ nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistPI[[2]])
#ols_plot_diagnostics(modlistPI[[2]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps_mean, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[3]] <- lm(pollution_index ~ nlcd_imp_ps, data = pollutfieldclean_cast)
ols_regress(modlistPI[[3]])
#ols_plot_diagnostics(modlistPI[[3]])
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[4]] <- lm(pollution_index ~ heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistPI[[4]])
#ols_plot_diagnostics(modlistPI[[4]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902log300proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[5]] <- lm(pollution_index ~ heatbing1902log200proj, data = pollutfieldclean_cast)
ols_regress(modlistPI[[5]])
#ols_plot_diagnostics(modlistPI[[5]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902log200proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[6]] <- lm(pollution_index ~ heatbing1902log500proj, data = pollutfieldclean_cast)
ols_regress(modlistPI[[6]])
#ols_plot_diagnostics(modlistPI[[6]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902log500proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[7]] <- lm(pollution_index ~ heatbing1902pow200_1proj, data = pollutfieldclean_cast)
ols_regress(modlistPI[[7]])
#ols_plot_diagnostics(modlistPI[[7]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902pow200_1proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[8]] <- lm(pollution_index ~ heatbing1902pow300_1proj, data = pollutfieldclean_cast)
ols_regress(modlistPI[[8]])
#ols_plot_diagnostics(modlistPI[[8]])
ggplot(pollutfieldclean_cast, aes(x=heatbing1902pow200_1proj, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[9]] <- lm(pollution_index ~ heatbustransitpow200_1, data = pollutfieldclean_cast)
ols_regress(modlistPI[[9]])
#ols_plot_diagnostics(modlistPI[[9]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitpow200_1, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[10]] <- lm(pollution_index ~ heatbustransitpow200_2, data = pollutfieldclean_cast)
ols_regress(modlistPI[[10]])
#ols_plot_diagnostics(modlistPI[[10]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitpow200_2, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[11]] <- lm(pollution_index ~ heatsubslopelog500, data = pollutfieldclean_cast)
ols_regress(modlistPI[[11]])
#ols_plot_diagnostics(modlistPI[[11]])
ggplot(pollutfieldclean_cast, aes(x=heatsubslopelog500, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[12]] <- lm(pollution_index ~ heatsubspdlpow500_2, data = pollutfieldclean_cast)
ols_regress(modlistPI[[12]])
#ols_plot_diagnostics(modlistPI[[12]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdlpow500_2, y=pollution_index, label=SiteIDPair)) + 
  geom_text()


modlistPI[[13]] <- lm(pollution_index ~ heatsubspdllog500, data = pollutfieldclean_cast)
ols_regress(modlistPI[[13]])
#ols_plot_diagnostics(modlistPI[[13]])
ggplot(pollutfieldclean_cast, aes(x=heatsubspdllog500, y=pollution_index, label=SiteIDPair)) + 
  geom_text()


modlistPI[[14]] <- lm(pollution_index ~ heatsubAADTpow500_2, data = pollutfieldclean_cast)
ols_regress(modlistPI[[14]])
#ols_plot_diagnostics(modlistPI[[14]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog200^(1/3), y=pollution_index, label=SiteIDPair)) + 
  geom_text()


modlistPI[[15]] <- lm(pollution_index ~ heatsubAADTlog500, data = pollutfieldclean_cast)
ols_regress(modlistPI[[15]])
#ols_plot_diagnostics(modlistPI[[15]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog500, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[16]] <- lm(pollution_index ~ heatsubAADTlog300, data = pollutfieldclean_cast)
ols_regress(modlistPI[[16]])
#ols_plot_diagnostics(modlistPI[[16]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog300, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[17]] <- lm(pollution_index ~ heatbustransitpow200_1thd, data = pollutfieldclean_cast)
ols_regress(modlistPI[[17]])
#ols_plot_diagnostics(modlistPI[[17]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitpow200_1^(1/3), y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[18]] <- lm(pollution_index ~ heatbustransitpow100_1thd, data = pollutfieldclean_cast)
ols_regress(modlistPI[[18]])
#ols_plot_diagnostics(modlistPI[[18]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitpow100_1^(1/3), y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[19]] <- lm(pollution_index ~ heatbustransitpow100_2thd, data = pollutfieldclean_cast)
ols_regress(modlistPI[[19]])
#ols_plot_diagnostics(modlistPI[[19]])
ggplot(pollutfieldclean_cast, aes(x=heatbustransitpow100_2^(1/3), y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[20]] <- lm(pollution_index ~ heatsubAADTlog200thd, data = pollutfieldclean_cast)
ols_regress(modlistPI[[20]])
#ols_plot_diagnostics(modlistPI[[20]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog200thd, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

modlistPI[[21]] <- lm(pollution_index ~ heatsubAADTlog100thd, data = pollutfieldclean_cast)
ols_regress(modlistPI[[21]])
#ols_plot_diagnostics(modlistPI[[21]])
ggplot(pollutfieldclean_cast, aes(x=heatsubAADTlog100thd, y=pollution_index, label=SiteIDPair)) + 
  geom_text()

#------ 2. PI - Multiparameter models --------
modlistPI[[22]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1thd, data = pollutfieldclean_cast)
ols_regress(modlistPI[[22]])
#ols_plot_diagnostics(modlistPI[[22]])
ols_coll_diag(modlistPI[[22]])
ols_correlations(modlistPI[[22]])

modlistPI[[23]] <- lm(pollution_index ~ heatsubAADTlog200thd + heatbustransitpow100_1thd, data = pollutfieldclean_cast)
ols_regress(modlistPI[[23]])
#ols_plot_diagnostics(modlistPI[[23]])
ols_coll_diag(modlistPI[[23]])
ols_correlations(modlistPI[[23]])

modlistPI[[24]] <- lm(pollution_index ~ heatsubAADTlog100thd + nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistPI[[24]])
#ols_plot_diagnostics(modlistPI[[24]])
ols_coll_diag(modlistPI[[24]])
ols_correlations(modlistPI[[24]])

modlistPI[[25]] <- lm(pollution_index ~ heatsubAADTlog200thd + nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistPI[[25]])
#ols_plot_diagnostics(modlistPI[[25]])
ols_coll_diag(modlistPI[[25]])
ols_correlations(modlistPI[[25]])

modlistPI[[26]] <- lm(pollution_index ~ heatsubAADTlog100thd*nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistPI[[26]])
#ols_plot_diagnostics(modlistPI[[26]])
ols_coll_diag(modlistPI[[26]])
ols_correlations(modlistPI[[26]])

modlistPI[[27]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistPI[[27]])
#ols_plot_diagnostics(modlistPI[[27]])
ols_coll_diag(modlistPI[[27]])
ols_correlations(modlistPI[[27]])


modlistPI[[28]] <- lm(pollution_index ~ heatsubAADTlog100thd*heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistPI[[28]])
#ols_plot_diagnostics(modlistPI[[28]])
ols_coll_diag(modlistPI[[28]])
ols_correlations(modlistPI[[28]])

modlistPI[[29]] <- lm(pollution_index ~ heatbustransitpow100_1thd + nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistPI[[29]])
#ols_plot_diagnostics(modlistPI[[29]])
ols_coll_diag(modlistPI[[29]])
ols_correlations(modlistPI[[29]])

modlistPI[[30]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistPI[[30]])
#ols_plot_diagnostics(modlistPI[[30]])
ols_coll_diag(modlistPI[[30]])
ols_correlations(modlistPI[[30]])

modlistPI[[31]] <- lm(pollution_index ~ heatsubAADTlog100 + heatbustransitpow100_1thd + nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistPI[[31]])
#ols_plot_diagnostics(modlistPI[[30]])
ols_coll_diag(modlistPI[[30]])
ols_correlations(modlistPI[[30]])

modlistPI[[32]] <- lm(pollution_index ~ heatsubAADTlog100*heatbustransitpow100_1 + nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistPI[[32]])
#ols_plot_diagnostics(modlistPI[[32]])
ols_coll_diag(modlistPI[[32]])
ols_correlations(modlistPI[[32]])

modlistPI[[33]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean + heatbing1902log300proj, data = pollutfieldclean_cast)
ols_regress(modlistPI[[33]])
#ols_plot_diagnostics(modlistPI[[33]])
ols_coll_diag(modlistPI[[33]])
ols_correlations(modlistPI[[33]])

modlistPI[[34]] <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean + heatsubspdllog500, data = pollutfieldclean_cast)
ols_regress(modlistPI[[34]])
#ols_plot_diagnostics(modlistPI[[34]])
ols_coll_diag(modlistPI[[34]])
ols_correlations(modlistPI[[34]])

modlistPI[[35]] <- lm(pollution_index ~ heatsubAADTlog100thd*heatsubspdllog500 + heatbustransitpow100_1 + nlcd_imp_ps_mean, data = pollutfieldclean_cast)
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
                                                          remove_outliers = 'outliers',
                                                          labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                                          kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary_nooutliers, -R2pred, AICc)  

model_summarypollutionindex_nooutliers_edit <- latex_format(model_summary_nooutliers) 
model_summarypollutionindex_nooutliers_edit <- gsub("pow100\\_1", "linear", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub("spdl", "speed", model_summarypollutionindex_nooutliers_edit) 
model_summarypollutionindex_nooutliers_edit <- gsub("nlcd\\_imp\\_ps", "imperviousness", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub("mean", "smooth", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub("bing1902", "congestion", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub("bustransit", "transit", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub("sub", "", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub("proj", "", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub("paperheight[=]10in[,]paperwidth[=]35in", "paperheight=8in,paperwidth=27in",
                                             model_summarypollutionindex_nooutliers_edit) 

cat(model_summarypollutionindex_nooutliers_edit, file = file.path(moddir, 'modeltable_pollutionindex_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_pollutionindex_nooutliers_2019.tex')

#------ 5. log PI - Run models for table -----
modlistlogPI <- list() #List to hold models
modlistlogPI[[1]] <- lm(logPI ~ 1, data = pollutfieldclean_cast) #Null/Intercept model

modlistlogPI[[2]] <- lm(logPI ~ nlcd_imp_ps_mean, data = pollutfieldclean_cast)

modlistlogPI[[3]] <- lm(logPI ~ nlcd_imp_ps, data = pollutfieldclean_cast)

modlistlogPI[[4]] <- lm(logPI ~ heatbing1902log300proj, data = pollutfieldclean_cast)

modlistlogPI[[5]] <- lm(logPI ~ heatbing1902log200proj, data = pollutfieldclean_cast)

modlistlogPI[[6]] <- lm(logPI ~ heatbing1902log500proj, data = pollutfieldclean_cast)

modlistlogPI[[7]] <- lm(logPI ~ heatbing1902pow200_1proj, data = pollutfieldclean_cast)

modlistlogPI[[8]] <- lm(logPI ~ heatbing1902pow300_1proj, data = pollutfieldclean_cast)

modlistlogPI[[9]] <- lm(logPI ~ heatbustransitpow200_1, data = pollutfieldclean_cast)

modlistlogPI[[10]] <- lm(logPI ~ heatbustransitpow200_2, data = pollutfieldclean_cast)

modlistlogPI[[11]] <- lm(logPI ~ heatsubslopelog500, data = pollutfieldclean_cast)

modlistlogPI[[12]] <- lm(logPI ~ heatsubspdlpow500_2, data = pollutfieldclean_cast)

modlistlogPI[[13]] <- lm(logPI ~ heatsubspdllog500, data = pollutfieldclean_cast)

modlistlogPI[[14]] <- lm(logPI ~ heatsubAADTpow500_2, data = pollutfieldclean_cast)

modlistlogPI[[15]] <- lm(logPI ~ heatsubAADTlog500, data = pollutfieldclean_cast)

modlistlogPI[[16]] <- lm(logPI ~ heatsubAADTlog300, data = pollutfieldclean_cast)

modlistlogPI[[17]] <- lm(logPI ~ heatbustransitpow200_1thd, data = pollutfieldclean_cast)

modlistlogPI[[18]] <- lm(logPI ~ heatbustransitpow100_1thd, data = pollutfieldclean_cast)

modlistlogPI[[19]] <- lm(logPI ~ heatbustransitpow100_2thd, data = pollutfieldclean_cast)

modlistlogPI[[20]] <- lm(logPI ~ heatsubAADTlog200thd, data = pollutfieldclean_cast)

modlistlogPI[[21]] <- lm(logPI ~ heatsubAADTlog100thd, data = pollutfieldclean_cast)

modlistlogPI[[22]] <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1thd, data = pollutfieldclean_cast)

modlistlogPI[[23]] <- lm(logPI ~ heatsubAADTlog200thd + heatbustransitpow100_1thd, data = pollutfieldclean_cast)

modlistlogPI[[24]] <- lm(logPI ~ heatsubAADTlog100thd + nlcd_imp_ps_mean, data = pollutfieldclean_cast)

modlistlogPI[[25]] <- lm(logPI ~ heatsubAADTlog200thd + nlcd_imp_ps_mean, data = pollutfieldclean_cast)

modlistlogPI[[26]] <- lm(logPI ~ heatsubAADTlog100thd*nlcd_imp_ps_mean, data = pollutfieldclean_cast)

modlistlogPI[[27]] <- lm(logPI ~ heatsubAADTlog100thd + heatbing1902log300proj, data = pollutfieldclean_cast)

modlistlogPI[[28]] <- lm(logPI ~ heatsubAADTlog100thd*heatbing1902log300proj, data = pollutfieldclean_cast)

modlistlogPI[[29]] <- lm(logPI ~ heatbustransitpow100_1thd + nlcd_imp_ps_mean, data = pollutfieldclean_cast)

modlistlogPI[[30]] <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean, data = pollutfieldclean_cast)

modlistlogPI[[31]] <- lm(logPI ~ heatsubAADTlog100 + heatbustransitpow100_1thd + nlcd_imp_ps_mean, data = pollutfieldclean_cast)

modlistlogPI[[32]] <- lm(logPI ~ heatsubAADTlog100*heatbustransitpow100_1 + nlcd_imp_ps_mean, data = pollutfieldclean_cast)

modlistlogPI[[33]] <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean + heatbing1902log300proj, data = pollutfieldclean_cast)

modlistlogPI[[34]] <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean + heatsubspdllog500, data = pollutfieldclean_cast)

modlistlogPI[[35]] <- lm(logPI ~ heatsubAADTlog100thd*heatsubspdllog500 + heatbustransitpow100_1 + nlcd_imp_ps_mean, data = pollutfieldclean_cast)

#------ 3. log PI - Make latex model summary table ----
model_summary_nooutliers <- as.data.table(
  ldply(modlistPI, function(mod) {regdiagnostic_customtab(mod, maxpar=vnum, 
                                                          remove_outliers = 'outliers & leverage',
                                                          labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
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
                                                             labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                                             kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary_nooutliers, -R2pred, AICc)  
cat(latex_format(model_summary_nooutliers), file = file.path(moddir, 'modeltable_logpollutionindex_nooutliers_2019.tex'))
setwd(moddir)
texi2pdf('modeltable_logpollutionindex_nooutliers_2019.tex')

#------ 9. Compare final selected models ----
logmod34_nooutliers <- regdiagnostic_customtab(modlistlogPI[[34]], maxpar=vnum, 
                                               remove_outliers = 'outliers',
                                               labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                               kCV = TRUE, k=10, cvreps=50)
subdatPI <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
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
                                               labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                                               kCV = TRUE, k=10, cvreps=50)
subdatPI <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                                      strsplit(gsub('\\\\', '', logmod30_nooutliers['outliers']), ',')$outliers),]
logmod30_nooutliersub <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean, 
                            data = subdatPI)
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
modlistglmPI[[1]] <- glm(pollution_index ~ 1, data = pollutfieldclean_cast, family = Gamma('log')) #Null/Intercept model

modlistglmPI[[2]] <- glm(pollution_index ~ nlcd_imp_ps_mean, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[2]])

modlistglmPI[[3]] <- glm(pollution_index ~ nlcd_imp_ps, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[3]])

modlistglmPI[[4]] <- glm(pollution_index ~ heatbing1902log300proj, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[4]])

modlistglmPI[[5]] <- glm(pollution_index ~ heatbing1902log500proj, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[5]])

modlistglmPI[[6]] <- glm(pollution_index ~ heatbing1902log500projsqrt, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[6]])

modlistglmPI[[7]] <- glm(pollution_index ~ heatbing1902pow200_1proj, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[7]])

modlistglmPI[[8]] <- glm(pollution_index ~ heatbing1902pow500_1proj, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[8]])

modlistglmPI[[9]] <- glm(pollution_index ~ heatbustransitpow200_1, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[9]])

modlistglmPI[[10]] <- glm(pollution_index ~ heatbustransitpow200_2, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[10]])

modlistglmPI[[11]] <- glm(pollution_index ~ heatsubslopelog500, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[11]])

modlistglmPI[[12]] <- glm(pollution_index ~ heatsubspdllog500, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[12]])

modlistglmPI[[13]] <- glm(pollution_index ~ heatsubspdllog500sqrt, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[13]])

modlistglmPI[[14]] <- glm(pollution_index ~ heatsubAADTpow500_2, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[14]])

modlistglmPI[[15]] <- glm(pollution_index ~ heatsubAADTlog500, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[15]])

modlistglmPI[[16]] <- glm(pollution_index ~ heatsubAADTlog300, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[16]])

modlistglmPI[[17]] <- glm(pollution_index ~ heatbustransitpow200_1thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[17]])

modlistglmPI[[18]] <- glm(pollution_index ~ heatbustransitpow100_1thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[18]])

modlistglmPI[[19]] <- glm(pollution_index ~ heatbustransitpow100_2thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[19]])

modlistglmPI[[20]] <- glm(pollution_index ~ heatsubAADTlog200thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[20]])

modlistglmPI[[21]] <- glm(pollution_index ~ heatsubAADTlog100thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[21]])

#------ 6. GLM PI - Multiparameter models --------
modlistglmPI[[22]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[22]])

modlistglmPI[[23]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog200thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[23]])

modlistglmPI[[24]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbustransitpow100_1thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[24]])

modlistglmPI[[25]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbustransitpow200_1thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[25]])

modlistglmPI[[26]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubspdllog500sqrt, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[26]])

modlistglmPI[[27]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbing1902log500proj, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[27]])

modlistglmPI[[28]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbing1902log500projsqrt, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[28]])

modlistglmPI[[29]] <- glm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[29]])

modlistglmPI[[30]] <- glm(pollution_index ~ heatsubAADTlog200thd + heatbustransitpow100_1, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[30]])

modlistglmPI[[31]] <- glm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[31]])

modlistglmPI[[32]] <- glm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow200_1thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[32]])

modlistglmPI[[33]] <- glm(pollution_index ~ heatsubAADTlog100 + heatbustransitpow100_1thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[33]])

modlistglmPI[[34]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd + heatbustransitpow200_1thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[34]])

modlistglmPI[[35]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd + heatbing1902log500projsqrt, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[35]])

modlistglmPI[[36]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatbustransitpow200_1thd + heatbing1902log500projsqrt, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[36]])

modlistglmPI[[37]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd + heatsubspdllog500sqrt, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[37]])

modlistglmPI[[38]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd + heatbustransitpow200_1thd + heatsubspdllog500sqrt, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[38]])

modlistglmPI[[39]] <- glm(pollution_index ~ nlcd_imp_ps_mean*heatsubAADTlog100thd, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[39]])

modlistglmPI[[40]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd* heatbing1902log500projsqrt, data = pollutfieldclean_cast, family = Gamma('log'))
GAMrescheck(modlistglmPI[[40]])

modlistglmPI[[41]] <- glm(pollution_index ~ nlcd_imp_ps_mean + heatsubAADTlog100thd* heatbustransitpow200_1thd, data = pollutfieldclean_cast, family = Gamma('log'))
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
                                            labelvec = pollutfieldclean_cast[, SiteIDPair],
                                            kCV = TRUE, k=10, cvreps=50)
subdatPI <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
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
                                            labelvec = pollutfieldclean_cast[, SiteIDPair],
                                            kCV = TRUE, k=10, cvreps=50)
subdatPI <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                                    strsplit(gsub('\\\\', '', mod30_nooutliers['outliers']), ',')$outliers),]
mod30_nooutliersub <- lm(pollution_index ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean, 
                         data = subdatPI)
ols_regress(mod30_nooutliersub)
#ols_plot_diagnostics(mod30_nooutliersub)
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
resnorm <- rstandard(mod30_nooutliersub) #Get standardized residuals from model
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
lm.morantest(mod30_nooutliersub, listw = listw2U(weightmat_k[[1]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod30_nooutliersub, listw = listw2U(weightmat_k[[2]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod30_nooutliersub, listw = listw2U(weightmat_k[[3]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod30_nooutliersub, listw = listw2U(weightmat_k[[4]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod30_nooutliersub, listw = listw2U(weightmat_k[[5]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod30_nooutliersub, listw = listw2U(weightmat_all)) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)

#Test for need for spatial regression model using Lagrange Multiplier (LM) tests
lm.LMtests(mod30_nooutliersub, listw = listw2U(weightmat_k[[1]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(mod30_nooutliersub, listw = listw2U(weightmat_k[[2]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(mod30_nooutliersub, listw = listw2U(weightmat_k[[3]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(mod30_nooutliersub, listw = listw2U(weightmat_k[[4]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))                     
lm.LMtests(mod30_nooutliersub, listw = listw2U(weightmat_k[[5]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))  
lm.LMtests(mod30_nooutliersub, listw = listw2U(weightmat_all), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))

#Spatial simultaneous autoregressive error model estimation with 1 nearest neighbors
sarlm_modPI <- errorsarlm(mod30_nooutliersub$call$formula, data = mod30_nooutliersub$model, 
                          listw = listw2U(weightmat_k[[1]]))
summary(sarlm_modPI)
bptest.sarlm(sarlm_modPI)
cor(mod30_nooutliersub$model$pollution_index, fitted(sarlm_modPI))^2

#Compare pseudo-R2 when predicting the trend
cor(mod30_nooutliersub$model$pollution_index, predict(sarlm_modPI, pred.type='trend'))^2
cor(mod30_nooutliersub$model$pollution_index, fitted(mod30_nooutliersub))^2
cor(pollutfieldclean_cast$pollution_index , with(pollutfieldclean_cast, sarlm_modPI$coefficients['(Intercept)'] + 
                                                   sarlm_modPI$coefficients['heatsubAADTlog100thd']*heatsubAADTlog100thd  + 
                                                   sarlm_modPI$coefficients['heatbustransitpow100_1']*heatbustransitpow100_1  + 
                                                   sarlm_modPI$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean))^2
#Compare MAPE
DescTools::MAPE(mod30_nooutliersub$model$pollution_index, fitted(sarlm_modPI))
DescTools::MAPE(mod30_nooutliersub$model$pollution_index, fitted(mod30_nooutliersub))

#Compare MAPE with all 
DescTools::MAPE(pollutfieldclean_cast$pollution_index ,with(pollutfieldclean_cast, sarlm_modPI$coefficients['(Intercept)'] + 
                                                     sarlm_modPI$coefficients['heatsubAADTlog100thd']*heatsubAADTlog100thd  + 
                                                     sarlm_modPI$coefficients['heatbustransitpow100_1']*heatbustransitpow100_1  + 
                                                     sarlm_modPI$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean))
DescTools::MAPE(pollutfieldclean_cast$pollution_index, predict(mod30_nooutliersub, pollutfieldclean_cast))

#Compare observed~predicted for full-no outlier model and for aspatial and spatial model
spatial_comparisonplot <- ggplot(subdatPI, aes(x=fitted(sarlm_modPI), y=pollution_index)) + 
  geom_point(aes(x=fitted(mod30_nooutliersub)), size=2, alpha=1/2, color='orange') +
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
############################################################################################################################################

# 6. Export models and data for subsequent analysis, make tables and figures for manuscript
############################################################################################################################################
#--------------- A. Export models ----
exportlist <- list(logZnstandmod = sarlm_modlogZnstand, Custandmod = sarlm_modCustand, PImod = sarlm_modPI) 

#logCustandmod = modlistlogCustand[[23]])
saveRDS(exportlist, file = file.path(moddir, 'fieldXRFmodels.rds'))

#Export scaling variable for all pollution drivers
pollutmaxscale <- dcast(pollutfieldclean[!is.na(mean)], 
                        formula = castformula,
                        value.var= 'mean')[, lapply(.SD, max), .SDcols = heatcols]
saveRDS(pollutmaxscale, file = file.path(moddir, 'fieldXRFmodels_scaling.rds'))

#--------------- B. Export data for Luwam ----
write.csv(pollutfieldclean_cast, 
          file.path(resdir, 'map_forluwam/XRFsites_pollutiondata_forLuwam_20190530.csv'),
          row.names = F)

#--------------- C. Make table of metal values and pollutant predictors for each site -----
#Generate table of statistics
tablecols <- c('#', 'Pair','Br','Cd', 'Co', 'Cr', 'Cu', 'Fe', 'Mn', 'Ni', 'Pb', 'Se', 'Sr', 'Ti',
               'Zn', 'Zr', "pollution_index", "nlcd_imp_ps", "LU", "POINT_X", "POINT_Y")
digitvec <- c(0,0,2,3,2,3,2,1,2,3,2,3,2,2,2,2,0,0,0,4,4)
LUcodes <- data.table(code = as.character(c(21, 41, 43, 81, 95, 96, 97, 98, 99)), 
                      LU = c("devlpd open", "decid. forest", "mixed forest", 
                             "pasture", "wetlands", "road", "resid.", "commer.", "indus."))

pollutable <- pollutfieldclean_cast[LUcodes, on='NLCD_reclass_final_PS==code'] %>%
  .[order(-pollution_index), 
    `#` := as.character(as.integer(factor(SiteID, 
                             levels=unique(pollutfieldclean_cast[order(-pollution_index),SiteID]))))] %>%
  .[, tablecols, with=F]

setnames(pollutable, old = c('pollution_index', 'nlcd_imp_ps', 'POINT_X', 'POINT_Y'),
         new = c('pollut.i', 'imperv. (%)', 'longitude', 'latitude'))

tablestats <- dcast(melt(
  pollutable[, .(colname = colnames(pollutable),
                 mean = sapply(.SD, function(x) {ifelse(is.numeric(x), mean(x, na.rm=T), NA)}),
                 SD = sapply(.SD, function(x) {ifelse(is.numeric(x),sd(x, na.rm=T), NA)}),
                 min = sapply(.SD, function(x) {ifelse(is.numeric(x),min(x, na.rm=T), NA)}),
                 max = sapply(.SD, function(x) {ifelse(is.numeric(x),max(x, na.rm=T), NA)})), 
             .SDcols = colnames(pollutable)],
  id.vars = "colname"), variable ~ colname) %>%
  .[, `#`:= variable] %>%
  .[, variable := NULL]

#Format and output table
options(knitr.table.format = "html") 
rbind(tablestats, pollutable[order(-pollut.i),])[, colnames(pollutable), with=F] %>%
  kable(format = "html", escape = F, 
        digits = digitvec) %>%
  kable_styling("striped", full_width = F) %>%
  save_kable(file.path(resdir, 'XRFrestab_20181219.doc'), self_contained=T)
  

#--------------- D. Make scatterplot of three final results -----
pollutfieldclean_cast[, `:=`(
  predZn_sarlm = exp(sarlm_modlogZnstand$coefficients['(Intercept)'] + 
                       sarlm_modlogZnstand$coefficients['heatsubAADTlog100frt']*heatsubAADTlog100frt + 
                       sarlm_modlogZnstand$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean),
  predCu_sarlm = sarlm_modCustand$coefficients['(Intercept)'] + 
    sarlm_modCustand$coefficients['heatbustransitlog200sqrt']*heatbustransitlog200sqrt + 
    sarlm_modCustand$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean,
  predPI_sarlm = sarlm_modPI$coefficients['(Intercept)'] + 
    sarlm_modPI$coefficients['heatsubAADTlog100thd']*heatsubAADTlog100thd  + 
    sarlm_modPI$coefficients['heatbustransitpow100_1']*heatbustransitpow100_1  + 
    sarlm_modPI$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean
)]

MAE_Znstand <- pollutfieldclean_cast[, DescTools::MAE(Znstand,predZn_sarlm)]
MAPE_Znstand <- pollutfieldclean_cast[, DescTools::MAPE(Znstand,predZn_sarlm)]
MAE_Custand <- pollutfieldclean_cast[, DescTools::MAE(Custand,predCu_sarlm)]
MAPE_Custand <- pollutfieldclean_cast[, DescTools::MAPE(Custand,predCu_sarlm)]
MAE_PI <- pollutfieldclean_cast[, DescTools::MAE(pollution_index, predPI_sarlm)]
MAPE_PI <- pollutfieldclean_cast[, DescTools::MAPE(pollution_index, predPI_sarlm)]

eqZn <-paste0("Zn = exp(", format(unname(sarlm_modlogZnstand$coefficients[1]), digits = 3), " + ",
              format(unname(sarlm_modlogZnstand$coefficients[2]), digits = 3), " . traffic volumelog1001/4 + ",
              format(unname(sarlm_modlogZnstand$coefficients[3]), digits = 3), " . imperviousness)")
eqCu <-paste0("Cu = ", format(unname(sarlm_modCustand$coefficients[1]), digits = 3), " + ",
              format(unname(sarlm_modCustand$coefficients[2]), digits = 3), " . bus transitlog2001/2 + ",
              format(unname(sarlm_modCustand$coefficients[3]), digits = 3), " . imperviousness)")
eqPI <-paste0("PI = ", format(unname(sarlm_modPI$coefficients[1]), digits = 3), " + ",
              format(unname(sarlm_modPI$coefficients[2]), digits = 3), " . traffic volumelog1001/3 + ",
              format(unname(sarlm_modPI$coefficients[3]), digits = 3), " . bus transitlinear100 + ",
              format(unname(sarlm_modPI$coefficients[4]), digits = 3), " . imperviousness")

marginp <- 0.3
gplotZn <- ggplot(pollutfieldclean_cast, aes(x=predZn_sarlm, 
                                             y=Znstand)) +
  geom_point(alpha=1/2, size=2, color='#35978f') + 
  geom_abline(intercept=0, slope=1) +
  scale_y_continuous(limits=c(0,100), expand=c(0,0), name = 'Observed Zn') + 
  scale_x_continuous(limits=c(0,100), expand=c(0,0), name = 'Predicted Zn') + 
  annotate(geom='text', label=paste0('MAE: ', round(MAE_Znstand, 1), 
                                     '\n MAPE: ', round(100*MAPE_Znstand, 0), '%'), x=90, y=15) + 
  ggtitle(eqZn) + 
  coord_fixed() +
  theme_classic() + 
  theme(plot.margin = margin(marginp, marginp, marginp, marginp, "cm"),
        text = element_text(size=12),
        plot.title = element_text(size=10))

gplotCu <- ggplot(pollutfieldclean_cast[SiteIDPair %in% subdat$SiteIDPair,], 
                  aes(x=predCu_sarlm, y=Custand)) +
  geom_point(alpha=1/2, size=2, color = '#8c510a') + 
  geom_point(data=pollutfieldclean_cast[!(SiteIDPair %in% subdat$SiteIDPair),], color='#525252', size=2, alpha=0.75) + 
  geom_abline(intercept=0, slope=1) +
  scale_y_continuous(limits=c(0,100), expand=c(0,0), name = 'Observed Cu') + 
  scale_x_continuous(limits=c(0,100), expand=c(0,0), name = 'Predicted Cu') + 
  annotate(geom='text', label=paste0('MAE: ', round(MAE_Custand, 1), 
                                     '\n MAPE: ', round(100*MAPE_Custand, 0), '%'), x=90, y=15) + 
  ggtitle(eqCu) + 
  coord_fixed() +
  theme_classic() + 
  theme(plot.margin = margin(marginp, marginp, marginp, marginp, "cm"),
        text = element_text(size=12),
        plot.title = element_text(size=10))

gplotPI<- ggplot(pollutfieldclean_cast[SiteIDPair %in% subdatPI$SiteIDPair,], 
                 aes(x=predPI_sarlm, y=pollution_index)) +
  geom_point(alpha=1/2, size=2, color = '#762a83') +  
  geom_point(data=pollutfieldclean_cast[!(SiteIDPair %in% subdatPI$SiteIDPair),], color='#525252', size=2, alpha=0.75) + 
  geom_abline(intercept=0, slope=1) +
  scale_y_continuous(limits=c(0,75), expand=c(0,0), name = 'Observed pollution index') + 
  scale_x_continuous(limits=c(0,75), expand=c(0,0), name = 'Predicted pollution index') + 
  annotate(geom='text', label=paste0('MAE: ', round(MAE_PI, 1), 
                                     '\n MAPE: ', round(100*MAPE_PI, 0), '%'), x=60, y=10) + 
  ggtitle(eqPI) + 
  coord_fixed() +
  theme_classic() + 
  theme(plot.margin = margin(marginp, marginp, marginp, marginp, "cm"),
        text = element_text(size=12),
        plot.title = element_text(size=10))

#Allow to keep axes 0-100 while drawing points that fall right on the edge
gg_tableZn <- ggplot_gtable(ggplot_build(gplotZn))
gg_tableZn$layout$clip[gg_tableZn$layout$name=="panel"] <- "off"
gg_tableCu <- ggplot_gtable(ggplot_build(gplotCu))
gg_tableCu$layout$clip[gg_tableCu$layout$name=="panel"] <- "off"
gg_tablePI <- ggplot_gtable(ggplot_build(gplotPI))
gg_tablePI$layout$clip[gg_tablePI$layout$name=="panel"] <- "off"

pdf(file.path(moddir, 'scatterplots_ZnCuPI_20191220.pdf'), width=4, height=4)
grid.draw(cbind(arrangeGrob(gg_tableZn, gg_tablePI), arrangeGrob(gg_tableCu, rectGrob(gp=gpar(col=NA)))))
dev.off()

#--------------- C. Check predictions  -----
predtab<- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='sitescheck')) 
predtab <- predtab[,SiteIDPair := factor(paste0(F_, Pair))][pollutfieldclean_cast, on='SiteIDPair']

qplot(predtab$predPI_sarlm, predtab$predpi30)
qplot(predtab$predCu_sarlm, predtab$predcu19)
qplot(predtab$predZn_sarlm, predtab$predzn36)

#Export residuals
pollutfieldclean_cast[, `:=`(
  predZn_resid = (predZn_sarlm - Znstand)/Znstand,
  predCu_resid = (predCu_sarlm - Custand)/Custand,
  predPI_resid = (predPI_sarlm - pollution_index)/pollution_index)]
write.dbf(pollutfieldclean_cast, file.path(moddir, 'pollutfieldclean_cast.dbf'))


#--------------- E. Make graph of % area vs % total pollution -----
#This is a preliminary figure to show % area vs % total pollution for study area extent
predzntab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predzn36_tab'))
predcutab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predcu_tab'))
predpitab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predpi_tab'))

predtabs <- rbind(cbind(predzntab, 'predzn'), cbind(predcutab, 'predcu'), cbind(predpitab, 'predpi'))
predtabs[, Value_stand := (Value - min(Value))/(max(Value)-min(Value)),by=V2]

#Percentile # of cells
predtabs[order(-Value_stand), `:=`(cumcount = cumsum(Count),
                                   cumvalue = cumsum((Value_stand)*Count)), by=V2]

predzntab <- predtabs[V2 == 'predzn',]

#Function to compute cumulative area and pollution from raster attribute
cumpollutarea <- function(tab) {
  outcum <- unlist(
    sapply(round(seq(0, sum(tab$Count), sum(tab$Count)/1000)), function(x) { 
      tab[cumcount >= x, ][which.min(cumcount), cumvalue - (Value_stand)*(cumcount-x)]/tab[, max(cumvalue)]
    })
  )
  return(outcum)
}

cumsumdt <- data.table (cumpollution = 100*c(cumpollutarea(predtabs[V2 == 'predzn',]),
                                             cumpollutarea(predtabs[V2 == 'predcu',]),
                                             cumpollutarea(predtabs[V2 == 'predpi',])),
                        percarea = rep(seq(0,100, 0.1), 3),
                        var = factor(c(rep('cumzn', 1001), rep('cumcu', 1001), rep('cumpi', 1001)), levels=c('cumzn', 'cumcu', 'cumpi')))


cumplot <- ggplot(cumsumdt[var %in% c('cumzn', 'cumcu')], aes(x=percarea, y=cumpollution)) +
  #geom_bar(aes(fill=var), stat = 'identity', alpha=0.5, position = "identity") + 
  #geom_area(aes(fill=var), stat = 'identity', alpha=0.3, position = "identity") +
  geom_vline(xintercept = cumsumdt[cumpollution>=50 & var == 'cumzn', min(percarea)], size=1, color='#01665e') +
  geom_vline(xintercept = cumsumdt[cumpollution>=50 & var == 'cumcu', min(percarea)], size=1, color='#543005') +
  geom_hline(yintercept = 50, color='black') +
  geom_line(aes(color=var), size=1, alpha=1) +
  scale_y_sqrt(expand=c(0,0), breaks = c(0, 1, 5, 10, 25, 50, 100), name = 'Cumulative pollution (%)') +
  scale_x_sqrt(expand=c(0,0), breaks = c(0, 1, 5, 10, 25, 50), name = 'Cumulative area, from\n most to least polluted (%)', limits = c(0,50)) +
  coord_fixed() +
  theme_classic() +
  #scale_color_manual(values = c('#35978f', '#8c510a','#762a83')) + 
  scale_fill_manual(values = c('#35978f', '#8c510a','#762a83')) + 
  #facet_wrap(~var, ncol=1) +
  theme(text = element_text(size=14),
        legend.position = 'None', 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin= margin(0.25, 0.25, 0.25, 0.25, "cm"))
pdf(file.path(moddir, 'cumpollution_cumarea.pdf'), width=3, height=5)
cumplot
dev.off()


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




