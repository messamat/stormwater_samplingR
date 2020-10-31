#Author: Mathis Messager
#Contact information: messamat@uw.edu
#Creation date: July 2018
#Purpose: import, format, merge, inspect, and model development of field XRF data

source('00_packages.R')
source('00_functions.R')
source('00_dirstructure.R')

#Import formatted data
datatab <- list.files(resdir, 'pollutfieldclean_cast.*[.]fst')
pollutfieldclean_cast <- read.fst(file.path(resdir, datatab[length(datatab)])) %>%
  setDT


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
#Check mixed effect model for temporal effect
pollutfieldclean_cast[, season := as.factor(season)]
AICc(modlistlogZnstand[[36]])

Znmod_seasoncheck <- lm(logZnstand ~ heatsubAADTlog100frt + nlcd_imp_ps_mean + season, 
                    data = pollutfieldclean_cast)
summary(Znmod_seasoncheck)
confint(Znmod_seasoncheck)
AICc(Znmod_seasoncheck)

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
Cusubdat <- pollutfieldclean_cast[!(paste0(SiteID, Pair) %in% 
                                    strsplit(gsub('\\\\', '', mod19_nooutliers['outliers']), ',')$outliers),]
mod19_nooutliersub <- lm(Custand ~ heatbustransitlog200sqrt +  nlcd_imp_ps_mean, 
                         data = Cusubdat)
summary(mod19_nooutliersub)
ols_plot_diagnostics(mod19_nooutliersub)
ols_regress(mod19_nooutliersub)
par(mfrow=c(2,2))
plot(mod19_nooutliersub)
ols_coll_diag(mod19_nooutliersub)
ols_correlations(mod19_nooutliersub)
AICc(mod19_nooutliersub)

qplot(Cusubdat$heatbing1902log300proj, mod19_nooutliersub$residuals) +
  geom_smooth(method='lm', color='red')
qplot(Cusubdat$heatsubspdllog300, mod19_nooutliersub$residuals) +
  geom_smooth(method='lm', color='red')
qplot(Cusubdat$heatsubAADTlog300, mod19_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(Cusubdat$heatsubslopelog500, mod19_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')

ggplot(Cusubdat, aes(x=NLCD_reclass_final_PS, y=mod19_nooutliersub$residuals)) +
  geom_boxplot() +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')

ggplot(pollutfieldclean_cast, aes(y=predict(modlistCustand[[19]], pollutfieldclean_cast), x=Custand)) +
  geom_point(alpha=1/4, size=2) + 
  #geom_point(data=Cusubdat, aes(y = predict(modlistCustand[[19]], pollutfieldclean_cast)), alpha=1/2, size=2) + 
  geom_point(data=Cusubdat, aes(x=Custand, y = predict(mod19_nooutliersub, Cusubdat)), color='red', alpha=1/2, size=2) + 
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
#Check for temporal effect
AICc(mod19_nooutliersub)
mod19_nooutliersub
Cumod_seasoncheck <- lm(Custand ~ heatbustransitlog200sqrt + nlcd_imp_ps_mean + season, 
                        data = Cusubdat)
summary(Cumod_seasoncheck)
confint(Cumod_seasoncheck)
AICc(Cumod_seasoncheck)

#Check for spatial autocorrelation
par(mfrow=c(1,1))
resnorm <- rstandard(mod19_nooutliersub) #Get standardized residuals from model without outliers (even when removing points that are outliers and leverage, can't get normal residuals)
#Make bubble map of residuals
bubbledat <- data.frame(resnorm, Cusubdat$coords.x1, Cusubdat$coords.x2)
coordinates(bubbledat) <- c("Cusubdat.coords.x1","Cusubdat.coords.x2")
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
  weightmat_IDW(Cusubdat[, .(coords.x1, coords.x2)], knb = i, mindist = 10)}) #Based on 1-10 nearest neighbors
weightmat_all <- weightmat_IDW(Cusubdat[, .(coords.x1, coords.x2)], knb = NULL, mindist = 10) #Based on all points

#Moran plots
#lag_resnorm <- lag.listw(weightmat_all, resnorm) #Can be used to create customized Moran plot by plotting residuals against matrix
moran.plot(resnorm, weightmat_all, labels=Cusubdat[,paste0(SiteID, Pair)], pch=19)
moran.plot(resnorm, weightmat_k[[2]], labels=Cusubdat[,paste0(SiteID, Pair)], pch=19)

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
spatial_comparisonplot <- ggplot(Cusubdat, aes(x=fitted(sarlm_modCustand), y=Custand)) + 
  geom_point(aes(x=fitted(mod19_nooutliersub)), size=2, alpha=1/2, color='orange') +
  geom_point(size=2, alpha=1/2, color='red') + 
  geom_abline(size=1, slope=1, intercept=0, color='red') + 
  #geom_text(aes(label=paste0(SiteID, Pair))) +
  coord_fixed() +
  theme_classic()
spatial_comparisonplot

resnorm_postsarlm <- residuals(sarlm_modCustand) #Get standardized residuals from model
#Make bubble map of residuals
bubbledat_postsarlm <- data.frame(resnorm_postsarlm, Cusubdat$coords.x1, Cusubdat$coords.x2)
coordinates(bubbledat_postsarlm) <- c("Cusubdat.coords.x1","Cusubdat.coords.x2")
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
ggplot(pollutfieldclean_cast, aes(x=nlcd_imp_ps_mean, y=pollution_index_old, label=SiteIDPair)) + 
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

modlistPI[[31]] <- lm(pollution_index ~ heatsubAADTlog100 + heatbustransitpow100_1thd +
                        nlcd_imp_ps_mean, data = pollutfieldclean_cast)
ols_regress(modlistPI[[31]])
#ols_plot_diagnostics(modlistPI[[31]])
ols_coll_diag(modlistPI[[31]])
ols_correlations(modlistPI[[31]])

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
cat(latex_format(model_summary), file = file.path(moddir, 'modeltable_pollutionindex_2020.tex'))
setwd(moddir)
texi2pdf('modeltable_pollutionindex_2020.tex')

#------ 4. PI - Make latex model summary table when excluding outliers ----
model_summary_nooutliers <- as.data.table(
  ldply(modlistPI, function(mod) {
    regdiagnostic_customtab(mod, maxpar=vnum, 
                            remove_outliers = 'outliers',
                            labelvec = pollutfieldclean_cast[, paste0(SiteID, Pair)],
                            kCV = TRUE, k=10, cvreps=50)}))
setorder(model_summary_nooutliers, -R2pred, AICc)  

model_summarypollutionindex_nooutliers_edit <- latex_format(model_summary_nooutliers) 
model_summarypollutionindex_nooutliers_edit <- gsub(
  "pow100\\_1", "linear", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub(
  "spdl", "speed", model_summarypollutionindex_nooutliers_edit) 
model_summarypollutionindex_nooutliers_edit <- gsub(
  "nlcd\\_imp\\_ps", "imperviousness", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub(
  "mean", "smooth", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub(
  "bing1902", "congestion", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub(
  "bustransit", "transit", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub(
  "sub", "", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub(
  "proj", "", model_summarypollutionindex_nooutliers_edit, fixed=T) 
model_summarypollutionindex_nooutliers_edit <- gsub(
  "paperheight[=]10in[,]paperwidth[=]35in", "paperheight=8in,paperwidth=27in",
  model_summarypollutionindex_nooutliers_edit) 

cat(model_summarypollutionindex_nooutliers_edit, 
    file = file.path(moddir, 'modeltable_pollutionindex_nooutliers_2020.tex'))
setwd(moddir)
texi2pdf('modeltable_pollutionindex_nooutliers_2020.tex')

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
logmodlistPI[[31]] <- lm(logPI ~ heatsubAADTlog100thd + heatbustransitpow100_1 + nlcd_imp_ps_mean, 
                            data = subdatPI)
summary(logmodlistPI[[31]])
logmod30_nooutliers['outliers']
ols_plot_diagnostics(logmodlistPI[[31]])
logmod30_nooutliers
AICc(logmodlistPI[[31]])
qplot(subdat$heatbing1902log500proj, logmodlistPI[[31]]$residuals) +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdat$heatsubslopelog500, logmodlistPI[[31]]$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')


ggplot(pollutfieldclean_cast, aes(x=exp(predict(logmod34_nooutliersub, pollutfieldclean_cast, type='response')), 
                                  y=pollution_index)) +
  geom_point(alpha=1/2, size=2) + 
  geom_point(data = subdatPI, aes(x = exp(predict(logmodlistPI[[31]], subdatPI))), color='red', alpha=1/2, size=2) + 
  geom_abline(intercept=0, slope=1) +
  theme_classic()

#Compare with predictions without transformation
qplot(abs(exp(fitted(modlistlogPI[[34]]))-subdatPI[, pollution_index]), 
      abs(exp(fitted(modlistlogPI[[30]]))-subdatPI[, pollution_index])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model glm 34 - glm  nlcd_imp_ps_mean + AADTlog100thd + transitpow200_1thd ', 
       y='Absolute error for model 30 - nlcd_imp_ps_mean + AADTlog100thd + transitpow100_1')

qplot(abs(fitted(modlistPI[[31]])-subdatPI[, pollution_index]), 
      abs(exp(predict(logmodlistPI[[31]], subdatPI))-subdatPI[, pollution_index])) + 
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


mod31_nooutliers <- regdiagnostic_customtab(modlistPI[[31]], maxpar=vnum, 
                                            remove_outliers = 'outliers',
                                            labelvec = pollutfieldclean_cast[, SiteIDPair],
                                            kCV = TRUE, k=10, cvreps=50)
subdatPI <- pollutfieldclean_cast[
  !(paste0(SiteID, Pair) %in% 
      strsplit(gsub('\\\\', '', mod31_nooutliers['outliers']), ',')$outliers),]
mod31_nooutliersub <- lm(pollution_index ~ heatsubAADTlog100 + 
                           heatbustransitpow100_1thd + nlcd_imp_ps_mean, 
                         data = subdatPI)
ols_regress(mod31_nooutliersub)
#ols_plot_diagnostics(mod31_nooutliersub)
ols_coll_diag(mod31_nooutliersub)
ols_correlations(mod31_nooutliersub)
AICc(mod31_nooutliersub)
qplot(subdatPI$heatbing1902log500proj, mod31_nooutliersub$residuals) +
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')
qplot(subdatPI$heatsubspdllog500sqrt, mod31_nooutliersub$residuals) + 
  geom_smooth(span=1) + geom_smooth(method='lm', color='red')


ggplot(subdatPI, aes(x=predict(mod34_nooutliersub, subdatPI, type='response'), 
                     y=pollution_index)) +
  geom_point(alpha=1/2, size=2) + 
  geom_point(aes(x = predict(mod31_nooutliersub, subdatPI)), color='red', alpha=1/2, size=2) + 
  geom_abline(intercept=0, slope=1) +
  theme_classic()

#Compare with predictions without transformation
qplot(abs(fitted(modlistglmPI[[34]])-subdatPI[, pollution_index]), 
      abs(fitted(modlistPI[[31]])-subdatPI[, pollution_index])) + 
  geom_abline(slope=1) + 
  labs(x='Absolute error for model glm 34 - glm  nlcd_imp_ps_mean + AADTlog100thd + transitpow200_1thd ', 
       y='Absolute error for model 30 - nlcd_imp_ps_mean + AADTlog100thd + transitpow100_1')

#------ 10. Check spatial and temporal autocorrelation of residuals for full and robust datasets -------
#Check for temporal effect
AICc(mod31_nooutliersub)
mod31_nooutliersub
PImod_seasoncheck <- lm(pollution_index ~ heatsubAADTlog100 + 
                          heatbustransitpow100_1thd + 
                          nlcd_imp_ps_mean + season, 
                        data = subdatPI)
summary(PImod_seasoncheck)
confint(PImod_seasoncheck)
AICc(PImod_seasoncheck)

"Fron Anselin 2006: ignoring spatially correlated errors is mostly a problem of efficiency, in the
sense that the OLS coefficient standard error estimates are biased, but the
coefficient estimates themselves remain unbiased. However, to the extent that
the spatially correlated errors mask an omitted variable, the consequences of
ignoring this may be more serious."
#Other good resource: https://eburchfield.github.io/files/Spatial_regression_LAB.html

par(mfrow=c(1,1))
resnorm <- rstandard(mod31_nooutliersub) #Get standardized residuals from model
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
lm.morantest(mod31_nooutliersub, listw = listw2U(weightmat_k[[1]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod31_nooutliersub, listw = listw2U(weightmat_k[[2]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod31_nooutliersub, listw = listw2U(weightmat_k[[3]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod31_nooutliersub, listw = listw2U(weightmat_k[[4]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod31_nooutliersub, listw = listw2U(weightmat_k[[5]])) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)
lm.morantest(mod31_nooutliersub, listw = listw2U(weightmat_all)) #lisw2U Make sure that distance matrix is symmetric (assumption in Moran's I)

#Test for need for spatial regression model using Lagrange Multiplier (LM) tests
lm.LMtests(mod31_nooutliersub, listw = listw2U(weightmat_k[[1]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(mod31_nooutliersub, listw = listw2U(weightmat_k[[2]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(mod31_nooutliersub, listw = listw2U(weightmat_k[[3]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))
lm.LMtests(mod31_nooutliersub, listw = listw2U(weightmat_k[[4]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))                     
lm.LMtests(mod31_nooutliersub, listw = listw2U(weightmat_k[[5]]), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))  
lm.LMtests(mod31_nooutliersub, listw = listw2U(weightmat_all), test=c("LMerr","RLMerr", "RLMlag", "SARMA"))

#Spatial simultaneous autoregressive error model estimation with 1 nearest neighbors
sarlm_modPI <- errorsarlm(mod31_nooutliersub$call$formula, data = mod31_nooutliersub$model, 
                          listw = listw2U(weightmat_k[[1]]))
summary(sarlm_modPI)
bptest.sarlm(sarlm_modPI)
cor(mod31_nooutliersub$model$pollution_index, fitted(sarlm_modPI))^2

#Compare pseudo-R2 when predicting the trend
cor(mod31_nooutliersub$model$pollution_index, predict(sarlm_modPI, pred.type='trend'))^2
cor(mod31_nooutliersub$model$pollution_index, fitted(mod31_nooutliersub))^2
cor(pollutfieldclean_cast$pollution_index , 
    with(pollutfieldclean_cast, sarlm_modPI$coefficients['(Intercept)'] + 
           sarlm_modPI$coefficients['heatsubAADTlog100']*heatsubAADTlog100  + 
           sarlm_modPI$coefficients['heatbustransitpow100_1thd']*heatbustransitpow100_1thd  + 
           sarlm_modPI$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean))^2
#Compare MAPE
DescTools::MAPE(mod31_nooutliersub$model$pollution_index, fitted(sarlm_modPI))
DescTools::MAPE(mod31_nooutliersub$model$pollution_index, fitted(mod31_nooutliersub))

#Compare MAPE with all 
DescTools::MAPE(pollutfieldclean_cast$pollution_index ,
                with(pollutfieldclean_cast, sarlm_modPI$coefficients['(Intercept)'] + 
                       sarlm_modPI$coefficients['heatsubAADTlog100']*heatsubAADTlog100  + 
                       sarlm_modPI$coefficients['heatbustransitpow100_1thd']*heatbustransitpow100_1thd  + 
                       sarlm_modPI$coefficients['nlcd_imp_ps_mean']*nlcd_imp_ps_mean))
DescTools::MAPE(pollutfieldclean_cast$pollution_index, predict(mod31_nooutliersub, pollutfieldclean_cast))

#Compare observed~predicted for full-no outlier model and for aspatial and spatial model
spatial_comparisonplot <- ggplot(subdatPI, aes(x=fitted(sarlm_modPI), y=pollution_index)) + 
  geom_point(aes(x=fitted(mod31_nooutliersub)), size=2, alpha=1/2, color='orange') +
  geom_point(aes(x=fitted(sarlm_modPI)), size=2, alpha=1/2, color='green') +
  geom_point(size=2, alpha=1/2, color='red') + 
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
exportlist <- list(logZnstandmod = sarlm_modlogZnstand, 
                   Custandmod = sarlm_modCustand,
                   PImod = sarlm_modPI) 

#logCustandmod = modlistlogCustand[[23]])
saveRDS(exportlist, file = file.path(moddir, 'fieldXRFmodels.rds'))

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
    sarlm_modPI$coefficients['heatsubAADTlog100']*heatsubAADTlog100  + 
    sarlm_modPI$coefficients['heatbustransitpow100_1thd']*heatbustransitpow100_1thd  + 
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


#--------------- E. Make graph of % area vs % total pollution + population + EJ -----
#Import pollution prediction tables
predzntab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predzn_tab'))
predcutab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predcu_tab'))
predpitab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predpi_tab'))

predznhightab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predznhighway_tab'))
predcuhightab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predcuhighway_tab'))
predpihightab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predpihighway_tab'))

#Get data on blockgroup dimensions and populations in PS watershed
EJWA <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer="EJWA_PS"))
EJWA[, quantile(Shape_Area, c(0.25,0.5,0.75))]/1000000
EJWA[, list(median(ACSTOTPOP), mean(ACSTOTPOP))]

#Make scatterplot of pollution vs EJ index
sampledpoints <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='randompoints_ps'))
sampledpoints_sub <- sampledpoints[!is.na(EJWApixeldipop10000sub), !"SHAPE"] %>%
  melt(id.vars=c('EJWApixeldipop10000sub', 'hpmstiger_PS_urbanhighways_buf500ras'),
       value.vars=c("predzn36indexland", "predcu19indexland")) %>%
  .[is.na(value), value:=0] %>%
  .[, highway := fifelse(is.na(hpmstiger_PS_urbanhighways_buf500ras),
                         'Beyond 500m from urban highway',
                         'Within 500m from urban highway')]

EJmetals_scatter <- ggplot(sampledpoints_sub, 
                           aes(x=100*EJWApixeldipop10000sub/max(sampledpoints_sub$EJWApixeldipop10000sub), 
                               y=value/100,)) + 
  geom_point(alpha=1/2, aes(color=as.factor(variable))) + 
  scale_color_manual(values=c("#81C2C7", "#C9A7CC"), labels=c('Zn', 'Cu')) +
  ggnewscale::new_scale_color() +
  geom_smooth(size=1.2, aes(color=as.factor(variable), fill=as.factor(variable))) +
  scale_color_manual(values=c("#22786F", "#7B4582"), labels=c('Zn', 'Cu')) +
  scale_fill_manual(values=c("#81C2C7", "#C9A7CC")) +
  scale_x_sqrt(expand=c(0,0), breaks=c(0,1,10,25,50,75,100), name= 'Environmental Justice index (scaled)') + 
  scale_y_sqrt(limits=c(0,100), expand=c(0,0), breaks=c(0,1,10,25,50,75,100), name='Metal index') + 
  coord_fixed() +
  facet_wrap(~highway) +
  theme_classic() + 
  theme(text = element_text(size=12),
        legend.title = element_blank())

png(file.path(moddir, 'EJindex_metalpollution_scatterplot_20200531.png'), width=10, height=8, units = 'in', res=600)
#pdf(file.path(moddir, 'EJindex_metalpollution_scatterplot_20200301.pdf'), width=6, height=5)
EJmetals_scatter
dev.off()

#This is a preliminary figure to show % area vs % total pollution for study area extent
EJWApixelpop10000tab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='EJWApixelpop10000sub_tab'))
EJWApixeldipop10000tab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='EJWApixeldipop10000sub_tab'))
predznpop10000tab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predznpop10000_tab'))
predcupop10000tab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predcupop10000_tab'))
predzndiratio10000tab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predzndiratio10000_tab'))
predcudiratio10000tab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predcudiratio10000_tab'))

EJWApixelpop10000hightab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='EJWApixelpop10000subhighway_tab'))
EJWApixeldipop10000hightab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='EJWApixeldipop10000subhighway_tab'))
predznpop10000hightab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predznpop10000highway_tab'))
predcupop10000hightab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predcupop10000highway_tab'))
predzndiratio10000hightab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predzndiratio10000highway_tab'))
predcudiratio10000hightab <- as.data.table(sf::st_read(dsn = file.path(resdir, 'PSpredictions.gdb'), layer='predcudiratio10000highway_tab'))

#Aggregate tables
weightingtab_zn <- rbind(cbind(predzntab, Target='predmetal'),
                         cbind(EJWApixelpop10000tab, Target='pop'), 
                         cbind(EJWApixeldipop10000tab, Target='diweight'),
                         cbind(predznpop10000tab, Target='predmetalpop'),
                         cbind(predzndiratio10000tab, Target='predmetaldi'))

weightingtab_cu <- rbind(cbind(predzntab, Target='predmetal'),
                         cbind(EJWApixelpop10000tab, Target='pop'), 
                         cbind(EJWApixeldipop10000tab, Target='diweight'),
                         cbind(predcupop10000tab, Target='predmetalpop'),
                         cbind(predcudiratio10000tab, Target='predmetaldi'))

weightingtab_join <- rbind(cbind(weightingtab_zn, Metal='Zn'),
                           cbind(weightingtab_cu, Metal='Cu'))


weightingtab_znhigh <- rbind(cbind(predznhightab, Target='predmetal'),
                         cbind(EJWApixelpop10000hightab, Target='pop'), 
                         cbind(EJWApixeldipop10000hightab, Target='diweight'),
                         cbind(predznpop10000hightab, Target='predmetalpop'),
                         cbind(predzndiratio10000hightab, Target='predmetaldi'))

weightingtab_cuhigh <- rbind(cbind(predznhightab, Target='predmetal'),
                         cbind(EJWApixelpop10000hightab, Target='pop'), 
                         cbind(EJWApixeldipop10000hightab, Target='diweight'),
                         cbind(predcupop10000hightab, Target='predmetalpop'),
                         cbind(predcudiratio10000hightab, Target='predmetaldi'))

weightingtab_joinhigh <- rbind(cbind(weightingtab_znhigh, Metal='Zn'),
                           cbind(weightingtab_cuhigh, Metal='Cu'))

weightingtab_joinall <- merge(weightingtab_join,
                              weightingtab_joinhigh, 
                              by=c('Value', 'Metal', 'Target'),
                              all.x=T, all.y=T, suffixes = c("", "_highway")) %>%
  .[is.na(Count_highway), Count_highway := 0] %>%
  .[, Count_nohighway := Count - Count_highway] %>%
  melt(id.vars=c('Metal', 'Target', 'Value'), 
       value.name = 'Count', variable.name = 'Variable')
  
#One graph for Cu, another for Zn
#Show accumulation by pollution only, population only, and DI ratio*population, pollution*population and population*DI ratio*population


#Percentile # of cells
weightingtab_joinall[order(-Value), `:=`(cumcount = cumsum(Count),
                                         cumvalue = cumsum((Value)*Count)), 
                     by=.(Target, Metal, Variable)]

#Function to compute cumulative area and pollution from raster attribute (messy set up but it works)
cumpollutarea <- function(tab) {
  outcum <- unlist(
    sapply(round(seq(0, sum(tab$Count), sum(tab$Count)/1000)), function(x) { 
      tab[cumcount >= x, ][which.min(cumcount), cumvalue - (Value)*(cumcount-x)]/tab[, max(cumvalue)]
    })
  )
  return(outcum)
}

cumcalc <- function(x, met, countvar) {
  data.table(
    cumval = cumpollutarea(
      weightingtab_joinall[!is.na(Count) & 
                             Target == x & 
                             Metal==met & 
                             Variable==countvar,]),
    percarea = seq(0,100, 0.1),
    Target = rep(x, 1001),
    Metal = rep(met, 1001),
    Variable = rep(countvar, 1001)
  )
}

valgrid <- weightingtab_joinall[, expand.grid(x=unique(Target), 
                                              met=unique(Metal))]

cumdt <- lapply(unique(weightingtab_joinall$Variable), function(countvar) {
  rbindlist(mapply(cumcalc, valgrid$x, valgrid$met, countvar,
                   SIMPLIFY = FALSE, USE.NAMES = FALSE))
}) %>%
  rbindlist()


targetsel <- c('predmetal', 'predmetaldi')
subtab <- cumdt[Target %in% targetsel,] %>%
  .[, Variable := factor(Variable, 
                         levels=c('Count', 
                                  'Count_nohighway', 
                                  'Count_highway'),
                         labels=c('All', 
                                  '> 500 m from highway', 
                                  '< 500 m for highway'))]

colorpal <- c('#D35F27','#0773B2', '#E69F25','#CC79A7', '#05A072') %>%
  .[seq_along(targetsel)]
vline_labels <- cbind(
  subtab[cumval>=0.5, min(percarea), by=.(Target, Metal)],
  color = rep(colorpal,2)) %>% 
  setorder(V1)

vlinetab <- subtab[cumval>=0.5, min(percarea), 
                   by=.(Target, Metal, Variable)]

#Only look at metal-only prioritization and metal+EJ*pop proritization
cumplot <- ggplot(subtab, aes(x=percarea, y=100*cumval, color=Target, linetype=Metal)) +
  #geom_bar(aes(fill=Target), stat = 'identity', alpha=0.5, position = "identity") + 
  #geom_area(aes(fill=Target), stat = 'identity', alpha=0.3, position = "identity") +
  geom_vline(data=vlinetab, alpha=1/2,
             aes(xintercept = V1, color=Target, linetype=Metal)) + #, color='#01665e') +
  geom_hline(yintercept = 50, color='black', alpha=1/2) +
  geom_line(aes(), size=1.2, alpha=1) +
  scale_y_sqrt(expand=c(0,0), breaks = c(0, 1, 5, 10, 25, 50, 100),
               name = 'Cumulative weight (%)') +
  scale_x_sqrt(expand=c(0,0), breaks = c(0, 1, 5, 10, 25, 50, 100), 
               name = 'Cumulative area, from\n most to least weight (%)',
               limits = c(0,100)) +
  coord_fixed() +
  theme_classic() +
  # annotate(geom='text',angle=45,
  #          color=vline_labels$color, 
  #          label=vline_labels$V1 , 
  #          x=vline_labels$V1, 
  #          y= 3*c(0.25, 0.75, 1.5, 1.5, 0.25, 0.75, 1.5, 1.5, 0.25, 0.25)) +
  scale_color_manual(values = colorpal[1:length(unique(subtab$Target))],
                     labels = c('pollution-only', 'pollution-EJ (v)')) +  #c('pollution-only', 'population-only (ii)', 'EJ-only (iii)', 'pollution-population (iv)', 'pollution-EJ (v)'))
  facet_wrap(~Variable, ncol=1) +
  theme(text = element_text(size=12),
        #legend.position = 'None', 
        strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.margin= margin(0.25, 0.25, 0.25, 0.25, "cm"))
cumplot

#pdf(file.path(moddir, 'cumpollution_cumarea.pdf'), width=6, height=5)
png(file.path(moddir, 'cumpollution_cumarea.png'), width=4, height=8, 
    units = 'in', res=600)
cumplot
dev.off()


kable(vlinetab) %>%
  kable_styling("striped", full_width = F) %>%
  save_kable(file.path(moddir, 'cumpollution_cumarea_50lines.doc'), self_contained=T)


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




