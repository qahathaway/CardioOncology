
##Load Dataset##
MESA <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/ALL_PRE.csv'))

##Impute Median Values##
detach("package:randomForestSRC", unload = TRUE)
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
dat <- as.data.frame(imputed$data)

#Load Packages
library(survival)
library(riskRegression)
library(randomForestSRC)
library(pec)
library(prodlim)

#COXPH Regression

coxOmicsM <- coxph(Surv(Mortality_Duration, Mortality) ~PW_PRE_CONVENTIONAL_min_ED+
                                  IVS_PRE_GLZLM_LZHGE_ED+
                                  PW_PRE_GLZLM_LZHGE_ED+
                                  IVS_PRE_DISCRETIZED_TLG+
                                  IVS_PRE_CONVENTIONAL_min_ES+
                                  IVS_PRE_CONVENTIONAL_Q1_ED+
                                  IVS_PRE_DISCRETIZED_TLG+
                                  IVS_PRE_GLZLM_HGZE_ED+
                                  PW_PRE_GLRLM_LRHGE_ED+
                                  IVS_PRE_GLZLM_SZHGE_ED+
                                  IVS_PRE_DISCRETIZED_Q1_ED+
                                  IVS_PRE_CONVENTIONAL_mean_ED+
                                  IVS_PRE_CONVENTIONAL_Q2_ED+
                                  IVS_PRE_DISCRETIZED_mean_ED+
                                  IVS_PRE_GLRLM_LRHGE_ED+
                                  IVS_PRE_GLRLM_HGRE_ED+
                                  IVS_PRE_GLRLM_SRHGE_ED+
                                  PW_PRE_GLRLM_HGRE_ED+
                                  IVS_PRE_CONVENTIONAL_Q3_ED+
                                  PW_PRE_GLRLM_SRHGE_ED+
                                  IVS_PRE_DISCRETIZED_Q2_ED+
                                  IVS_PRE_DISCRETIZED_Q3_ED+
                                  IVS_PRE_CONVENTIONAL_Q1_ES+
                                  PW_PRE_CONVENTIONAL_Q1_ED+
                                  PW_PRE_GLZLM_ZP_ED+
                                  PW_PRE_GLRLM_RP_ED+
                                  IVS_PRE_CONVENTIONAL_min_ED+
                                  IVS_PRE_DISCRETIZED_AUC_CSH_ED+
                                  PW_PRE_CONVENTIONAL_Q2_ED+
                                  PW_PRE_GLZLM_SZE_ED+
                                  PW_PRE_GLZLM_HGZE_ED+
                                  PW_PRE_DISCRETIZED_mean_ED+
                                  PW_PRE_CONVENTIONAL_mean_ED+
                                  PW_PRE_DISCRETIZED_Q2_ED+
                                  PW_PRE_DISCRETIZED_Q3_ED+
                                  IVS_PRE_CONVENTIONAL_mean_ES+
                                  IVS_PRE_CONVENTIONAL_Q2_ES+
                                  PW_PRE_DISCRETIZED_std_ES+
                                  PW_PRE_GLZLM_LZE_ED+
                                  PW_PRE_CONVENTIONAL_Skewness_ED+
                                  PW_PRE_DISCRETIZED_Skewness_ED+
                                  IVS_PRE_DISCRETIZED_Skewness_ED+
                                  IVS_PRE_DISCRETIZED_HISTO_Entropy_log10_ED+
                                  IVS_PRE_DISCRETIZED_HISTO_Entropy_log2_ED+
                                  IVS_PRE_CONVENTIONAL_Skewness_ED+
                                  PW_PRE_GLCM_Energy.AngularSecondMoment_ES+
                                  PW_PRE_CONVENTIONAL_Skewness_ES+
                                  PW_PRE_DISCRETIZED_Skewness_ES+
                                  PW_PRE_GLZLM_LZE_ES+
                                  PW_PRE_DISCRETIZED_Q1_ED+
                                  PW_PRE_GLZLM_ZP_ES+
                                  IVS_PRE_CONVENTIONAL_Q3_ES+
                                  PW_PRE_GLZLM_SZE_ES+
                                  PW_PRE_GLZLM_SZHGE_ED+
                                  IVS_PRE_GLRLM_LRHGE_ES+
                                  PW_PRE_CONVENTIONAL_Q3_ED+
                                  IVS_PRE_GLZLM_ZLNU_ES+
                                  IVS_PRE_DISCRETIZED_std_ED+
                                  IVS_PRE_GLRLM_HGRE_ES+
                                  IVS_PRE_GLZLM_GLNU_ED+
                                  IVS_PRE_DISCRETIZED_HISTO_Energy.Uniformity_ED+
                                  IVS_PRE_GLRLM_SRHGE_ES+
                                  PW_PRE_GLRLM_RP_ES,data=dat,x=TRUE,y=TRUE)

coxOmicsHF <- coxph(Surv(HF_Duration, HF) ~PW_PRE_CONVENTIONAL_min_ED+
                                IVS_PRE_CONVENTIONAL_std_ED+
                                PW_PRE_GLRLM_LRHGE_ED+
                                PW_PRE_GLRLM_HGRE_ED+
                                PW_PRE_GLRLM_SRHGE_ED+
                                PW_PRE_GLZLM_LZHGE_ED+
                                IVS_PRE_NGLDM_Coarseness_ED+
                                PW_PRE_CONVENTIONAL_Q1_ED+
                                PW_PRE_CONVENTIONAL_mean_ED+
                                IVS_PRE_DISCRETIZED_std_ED+
                                IVS_PRE_CONVENTIONAL_std_ES+
                                PW_PRE_CONVENTIONAL_Q2_ED+
                                IVS_PRE_CONVENTIONAL_Q3_ED+
                                PW_PRE_GLZLM_SZE_ED+
                                PW_PRE_CONVENTIONAL_Q3_ED+
                                PW_PRE_GLRLM_RP_ED+
                                IVS_PRE_CONVENTIONAL_max_ED+
                                IVS_PRE_DISCRETIZED_AUC_CSH_ED+
                                PW_PRE_GLZLM_ZP_ED+
                                PW_PRE_GLZLM_LZE_ED+
                                IVS_PRE_CONVENTIONAL_mean_ED+
                                PW_PRE_GLCM_Entropy_log10_ES+
                                IVS_PRE_DISCRETIZED_std_ES+
                                PW_PRE_GLRLM_LRHGE_ES+
                                PW_PRE_GLRLM_SRHGE_ES+
                                PW_PRE_GLRLM_HGRE_ES,data=dat,x=TRUE,y=TRUE)

coxDemM <- coxph(Surv(Mortality_Duration, Mortality) ~Age+BMI+HTN+HLD+DM+COPD+ASA+BB+Statin+
                    ACEi_ARB_Entresto+CCB+Diuretic+insulin+Metfromin+NSAID+smokinghx+currentsmoking+
                    Alcohol+Race,data=dat,x=TRUE,y=TRUE)

coxDemHF <- coxph(Surv(HF_Duration, HF) ~Age+BMI+HTN+HLD+DM+COPD+ASA+BB+Statin+
                    ACEi_ARB_Entresto+CCB+Diuretic+insulin+Metfromin+NSAID+smokinghx+currentsmoking+
                    Alcohol+Race,data=dat,x=TRUE,y=TRUE)

coxFuncM <- coxph(Surv(Mortality_Duration, Mortality) ~EF.Bi.Plane..Or.Visually.Estimated+IVSd+LVPWd+
                                 MV.E+MV.A+Lat.E.+Med.E.,data=dat,x=TRUE,y=TRUE)

coxFuncHF <- coxph(Surv(HF_Duration, HF) ~EF.Bi.Plane..Or.Visually.Estimated+IVSd+LVPWd+
                               MV.E+MV.A+Lat.E.+Med.E.,data=dat,x=TRUE,y=TRUE)

coxDrugM <- coxph(Surv(Mortality_Duration, Mortality) ~ATEZOLIZUMAB+CAPECITABINE+CARBOPLATIN+CELECOXIB+
                                 CYCLOPHOSPHAMIDE+DOCETAXEL+DOXORUBICIN+ERIBULIN+
                                 EVEROLIMUS+GEMCITABINE+METHOTREXATE+NERATINIB+
                                 PACLITAXEL+PALBOCICLIB+PEMBROLIZUMAB+PERTUZUMAB+
                                 TRASTUZUMAB,data=dat,x=TRUE,y=TRUE)

coxDrugHF <- coxph(Surv(HF_Duration, HF) ~ATEZOLIZUMAB+CAPECITABINE+CARBOPLATIN+CELECOXIB+
                               CYCLOPHOSPHAMIDE+DOCETAXEL+DOXORUBICIN+ERIBULIN+
                               EVEROLIMUS+GEMCITABINE+METHOTREXATE+NERATINIB+
                               PACLITAXEL+PALBOCICLIB+PEMBROLIZUMAB+PERTUZUMAB+
                               TRASTUZUMAB,data=dat,x=TRUE,y=TRUE)

Dem_M <- coxph(Surv(Mortality_Duration, Mortality) ~Age+BMI+ASA,data=dat,x=TRUE,y=TRUE)
Func_M <- coxph(Surv(Mortality_Duration, Mortality) ~Age+BMI+ASA+MV.E,data=dat,x=TRUE,y=TRUE)
Drug_M <- coxph(Surv(Mortality_Duration, Mortality) ~Age+BMI+ASA+DOCETAXEL,data=dat,x=TRUE,y=TRUE)
Omics_M <- coxph(Surv(Mortality_Duration, Mortality) ~Age+BMI+ASA+RiskScore_M,data=dat,x=TRUE,y=TRUE)

summary(Dem_M)
summary(Func_M)
summary(Drug_M)
summary(Omics_M)

Dem_HF <- coxph(Surv(HF_Duration, HF) ~Age+BMI+Metfromin,data=dat,x=TRUE,y=TRUE)
Func_HF <- coxph(Surv(HF_Duration, HF) ~Age+BMI+Metfromin+IVSd,data=dat,x=TRUE,y=TRUE)
Drug_HF <- coxph(Surv(HF_Duration, HF) ~Age+BMI+Metfromin+DOCETAXEL,data=dat,x=TRUE,y=TRUE)
Omics_HF <- coxph(Surv(HF_Duration, HF) ~Age+BMI+Metfromin+RiskScore_HF,data=dat,x=TRUE,y=TRUE)

summary(Dem_HF)
summary(Func_HF)
summary(Drug_HF)
summary(Omics_HF)


summary(coxDemM)
summary(coxFuncM)
summary(coxDrugM)
summary(coxOmicsM)
summary(coxDemHF)
summary(coxFuncHF)
summary(coxDrugHF)
summary(coxOmicsHF)

concordance(coxDemM)
concordance(coxFuncM)
concordance(coxDrugM)
concordance(coxOmicsM)
concordance(coxDemHF)
concordance(coxFuncHF)
concordance(coxDrugHF)
concordance(coxOmicsHF)

###Mortality###

AUC=Score(list(coxDemM, coxFuncM, coxDrugM, coxOmicsM), formula=Surv(Mortality_Duration, Mortality)~1,
                 data=dat, metrics="auc", method="bootcv", B=1000, null.model=TRUE, plots=c("calibration","ROC"), times=seq(0:3939))

Brier=Score(list(coxDemM, coxFuncM, coxDrugM, coxOmicsM), formula=Surv(Mortality_Duration, Mortality)~1,
         data=dat, metrics="brier", method="bootcv", B=1000, null.model=TRUE, plots=c("calibration","ROC"), times=seq(0:3939))

aucgraphDem <- plotAUC(AUC, models = "coxph", which = "score", xlim = c(0,4000), ylim = c(0,1), col = "#4682B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
aucgraphFunc <- plotAUC(AUC, models = "coxph.1", which = "score", xlim = c(0,4000), ylim = c(0,1), col = "#AF46B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
aucgraphDrug <- plotAUC(AUC, models = "coxph.2", which = "score", xlim = c(0,4000), ylim = c(0,1), col = "#B47846", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
aucgraphOmics <- plotAUC(AUC, models = "coxph.3", which = "score", xlim = c(0,4000), ylim = c(0,1), col = "#4BB446", lwd = 5, conf.int =TRUE, legend = FALSE)

write.csv(aucgraphDem, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/AUC_M_Dem.csv")
write.csv(aucgraphFunc, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/AUC_M_Func.csv")
write.csv(aucgraphDrug, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/AUC_M_Drug.csv")
write.csv(aucgraphOmics, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/AUC_M_Omics.csv")

briergraphDem <- plotBrier(Brier, models = "coxph", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#4682B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphFunc <- plotBrier(Brier, models = "coxph.1", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#AF46B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphDrug <- plotBrier(Brier, models = "coxph.2", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#B47846", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphOmics <- plotBrier(Brier, models = "coxph.3", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#4BB446", lwd = 5, conf.int =TRUE, legend = FALSE)

write.csv(briergraphDem, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Brier_M_Dem.csv")
write.csv(briergraphFunc, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Brier_M_Func.csv")
write.csv(briergraphDrug, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Brier_M_Drug.csv")
write.csv(briergraphOmics, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Brier_M_Omics.csv")

###HF###
AUCHF=Score(list(coxDemHF, coxFuncHF, coxDrugHF, coxOmicsHF), formula=Surv(HF_Duration, HF)~1,
          data=dat, metrics="auc", method="bootcv", B=1000, null.model=TRUE, plots=c("calibration","ROC"), times=seq(0:3754))

BrierHF=Score(list(coxDemHF, coxFuncHF, coxDrugHF, coxOmicsHF), formula=Surv(HF_Duration, HF)~1,
            data=dat, metrics="brier", method="bootcv", B=1000, null.model=TRUE, plots=c("calibration","ROC"), times=seq(0:3754))

aucgraphDemHF <- plotAUC(AUCHF, models = "coxph", which = "score", xlim = c(0,4000), ylim = c(0,1), col = "#4682B4", lwd = 5, conf.int =FALSE, legend = FALSE)
par(new=TRUE)
aucgraphFuncHF <- plotAUC(AUCHF, models = "coxph.1", which = "score", xlim = c(0,4000), ylim = c(0,1), col = "#AF46B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
aucgraphDrugHF <- plotAUC(AUCHF, models = "coxph.2", which = "score", xlim = c(0,4000), ylim = c(0,1), col = "#B47846", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
aucgraphOmicsHF <- plotAUC(AUCHF, models = "coxph.3", which = "score", xlim = c(0,4000), ylim = c(0,1), col = "#4BB446", lwd = 5, conf.int =TRUE, legend = FALSE)

write.csv(aucgraphDemHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/AUC_HF_Dem.csv")
write.csv(aucgraphFuncHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/AUC_HF_Func.csv")
write.csv(aucgraphDrugHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/AUC_HF_Drug.csv")
write.csv(aucgraphOmicsHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/AUC_HF_Omics.csv")

briergraphDemHF <- plotBrier(BrierHF, models = "coxph", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#4682B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphFuncHF <- plotBrier(BrierHF, models = "coxph.1", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#AF46B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphDrugHF <- plotBrier(BrierHF, models = "coxph.2", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#B47846", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphOmicsHF <- plotBrier(BrierHF, models = "coxph.3", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#4BB446", lwd = 5, conf.int =TRUE, legend = FALSE)

write.csv(briergraphDemHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Brier_HF_Dem.csv")
write.csv(briergraphFuncHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Brier_HF_Func.csv")
write.csv(briergraphDrugHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Brier_HF_Drug.csv")
write.csv(briergraphOmicsHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Brier_HF_Omics.csv")


###Mortalty###
ApparrentCindexM <- pec::cindex(list("COXPH Dem"=coxDemM,
                                    "COXPH Function"=coxFuncM,
                                    "COXPH Drugs"=coxDrugM,
                                    "COXPH Omics"=coxOmicsM),
              formula=Surv(Mortality_Duration, Mortality)~1,data=dat,
              eval.times=seq(0,3939,1), pred.times=seq(0,3939,1))

col = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")

plot(ApparrentCindexM, legend = FALSE, xlim=c(0,4000), ylim=c(0.5,1.0), lwd = 5, col = col)
write.csv(ApparrentCindexM$AppCindex, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Concordance_Mortality.csv")

###HF###
ApparrentCindexHF <- pec::cindex(list("COXPH Dem"=coxDemHF,
                                    "COXPH Function"=coxFuncHF,
                                    "COXPH Drugs"=coxDrugHF,
                                    "COXPH Omics"=coxOmicsHF),
                               formula=Surv(HF_Duration, HF)~1,data=dat,
                               eval.times=seq(0,3754,1), pred.times=seq(0,3754,1))

col = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")

plot(ApparrentCindexHF, legend = FALSE, xlim=c(0,4000), ylim=c(0.5,1.0), lwd = 5, col = col)
write.csv(ApparrentCindexHF$AppCindex, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Concordance_HF.csv")

##cumulative prediction error##

Models <-  list('CoxProb'=coxph(Surv(duration,event)~Prob,data=dat,x=TRUE,y=TRUE),
               'CoxALL'=coxph(Surv(duration,event)~.,data=dat,x=TRUE,y=TRUE))

# compute the apparent prediction error
PredError <- pec(object=Models,
                 formula=Surv(duration,event)~.,
                 data=dat,
                 exact=TRUE,
                 cens.model="marginal",
                 splitMethod="none",
                 B=0,
                 verbose=TRUE)
print(PredError,times=seq(189,189,1))
summary(PredError)
plot(PredError,xlim=c(0,200), legend = FALSE)
