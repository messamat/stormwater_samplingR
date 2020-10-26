#Modify rcompanion transformTukey to return lambda
#This Tukey ladder computation is designed to take in a raster attribute table and repeatedly sample it to determine whether a transformation 
#is required, and if so, performs the transformation
getlambda <- function(x, start = -2, end = 2, int = 0.025, verbose = FALSE, statistic = 1) {
  n = (end - start)/int
  lambda = as.numeric(rep(0, n))
  W = as.numeric(rep(0, n))
  Shapiro.p.value = as.numeric(rep(0, n))
  if (statistic == 2) {
    A = as.numeric(rep(1000, n))
    Anderson.p.value = as.numeric(rep(0, n))
  }
  for (i in (1:n)) {
    lambda[i] = signif(start + (i - 1) * int, digits = 4)
    if (lambda[i] > 0) {
      TRANS = x^lambda[i]
    }
    if (lambda[i] == 0) {
      TRANS = log(x)
    }
    if (lambda[i] < 0) {
      TRANS = -1 * x^lambda[i]
    }
    W[i] = NA
    if (statistic == 2) {
      A[i] = NA
    }
    if (any(is.infinite(TRANS)) == FALSE & any(is.nan(TRANS)) == 
        FALSE) {
      W[i] = signif(shapiro.test(TRANS)$statistic, digits = 4)
      Shapiro.p.value[i] = signif(shapiro.test(TRANS)$p.value, 
                                  digits = 4)
      if (statistic == 2) {
        A[i] = signif(ad.test(TRANS)$statistic, digits = 4)
        Anderson.p.value[i] = signif(ad.test(TRANS)$p.value, 
                                     digits = 4)
      }
    }
  }
  if (statistic == 2) {
    df = data.frame(lambda, W, Shapiro.p.value, A, Anderson.p.value)
  }
  if (statistic == 1) {
    df = data.frame(lambda, W, Shapiro.p.value)
  }
  if (verbose == TRUE) {
    print(df)
  }
  if (statistic == 1) {
    df2 = df[with(df, order(-W)), ]
  }
  if (statistic == 2) {
    df2 = df[with(df, order(A)), ]
  }
  cat("\n")
  print(df2[1, ])
  cat("\n")
  cat("if (lambda >  0){TRANS = x ^ lambda}", "\n")
  cat("if (lambda == 0){TRANS = log(x)}", "\n")
  cat("if (lambda <  0){TRANS = -1 * x ^ lambda}", "\n")
  cat("\n")
  return(df2[1, "lambda"])
}

transformTukey_lambda <- function (x,start = -2, end = 2, int = 0.025, rastertab=T, rep=10, verbose = FALSE, statistic = 1) {
  if (rastertab) { #If analyze a raster attribute table
    lambda_iter <- NULL
    for (iter in 1:rep) {
      x$prob <- x$Count/sum(x$Count)
      xsamp <- sample(x$Value, size=5000, replace=T, prob=x$prob)
      lambda_iter[iter] <- getlambda(x=xsamp)
      lambda <- mean(lambda_iter)
    }
    if (lambda > 0) {
      TRANS = x$Value^lambda
    }
    if (lambda == 0) {
      TRANS = log(x$Value)
    }
    if (lambda < 0) {
      TRANS = -1 * x$Value^lambda
    }
  } else{ #If processing a vector
    if (length(x)<=5000) { 
      lambda <- getlambda(x=x)
    } else {
      for (iter in 1:rep) {
        xsamp <- sample(x, size=5000)
        lambda_iter[iter] <- getlambda(x=xsamp)
        lambda <- mean(lambda_iter)
      }
    }
    if (lambda > 0) {
      TRANS = x^lambda
    }
    if (lambda == 0) {
      TRANS = log(x)
    }
    if (lambda < 0) {
      TRANS = -1 * x^lambda
    }
  }
  return(list(TRANS,lambda))
}

#Back transformTukey
backtransformTukey <- function (x, lambda) 
{
  if (lambda > 0) {
    BACKTRANS = x^(1/lambda)
  }
  if (lambda == 0) {
    BACKTRANS = exp(x)
  }
  if (lambda < 0) {
    BACKTRANS = (-x)^(1/lambda)
  }
  return(BACKTRANS)
}

#Function to bin raster attribute table
bin_rastertab<- function(tab, nbins, tukey=T, rep=10, rastertab) {
  df <- tab
  if (tukey) {
    dftrans <- transformTukey_lambda(df, rep=rep, rastertab=rastertab)
    df$valtuk <- dftrans[[1]]
    lambda <- dftrans[[2]]
  } else {
    df$valtuk <- df$Value
    lambda <- 1
  }
  binsize <- (max(df$valtuk)-min(df$valtuk))/nbins
  #minbin = min(df$valtuk)-min(df$valtuk)%%binsize+binsize
  #maxbin = max(df$valtuk)-max(df$valtuk)%%binsize+binsize
  bin <- with(df,seq(min(df$valtuk),max(df$valtuk), binsize))
  binback <- data.frame(bin,binmax=backtransformTukey(bin, lambda))
  setDT(binback)[,binmin:=shift(binmax, 1, type='lag'),]
  binback$id <- as.numeric(row.names(binback))-1
  df$bin <- .bincode(x=df$valtuk, breaks=bin, TRUE, include.lowest = T)
  if (rastertab) {
    df_hist <- setDT(df)[,sum(Count), .(bin)]
  } else {
    df_hist <- setDT(df)[,length(Value), .(bin)]
  }
  colnames(df_hist) <- c('id','count')
  df_hist <- merge(df_hist, binback, by='id', all.x=T, all.y=T)
  return(df_hist[-1,])
}

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