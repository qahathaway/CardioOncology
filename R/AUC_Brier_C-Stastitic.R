
##Load Dataset##
CardioOnc <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/ALL_PRE.csv'))

##Impute Median Values##
detach("package:randomForestSRC", unload = TRUE)
library(mlr)
imputed = impute(CardioOnc, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
dat <- as.data.frame(imputed$data)

#Load Packages
library(survival)
library(riskRegression)
library(randomForestSRC)
library(pec)
library(prodlim)

#COXPH Regression

coxOmicsM <- coxph(Surv(Mortality_Duration, Mortality) ~Ultrasomics-Features_Here,data=dat,x=TRUE,y=TRUE)

coxOmicsHF <- coxph(Surv(HF_Duration, HF) ~Ultrasomics-Features_Here,data=dat,x=TRUE,y=TRUE)

coxDemM <- coxph(Surv(Mortality_Duration, Mortality) ~Demographic-Features_Here,data=dat,x=TRUE,y=TRUE)

coxDemHF <- coxph(Surv(HF_Duration, HF) ~Demographic-Features_Here,data=dat,x=TRUE,y=TRUE)

coxFuncM <- coxph(Surv(Mortality_Duration, Mortality) ~Functional-Features_Here,data=dat,x=TRUE,y=TRUE)

coxFuncHF <- coxph(Surv(HF_Duration, HF) ~Functional-Features_Here,data=dat,x=TRUE,y=TRUE)

coxDrugM <- coxph(Surv(Mortality_Duration, Mortality) ~Therapy-Features_Here,data=dat,x=TRUE,y=TRUE)

coxDrugHF <- coxph(Surv(HF_Duration, HF) ~Therapy-Features_Here,data=dat,x=TRUE,y=TRUE)

Dem_M <- coxph(Surv(Mortality_Duration, Mortality) ~Age+BMI+HTN,data=dat,x=TRUE,y=TRUE)
Func_M <- coxph(Surv(Mortality_Duration, Mortality) ~Age+BMI+HTN+MV.E,data=dat,x=TRUE,y=TRUE)
Drug_M <- coxph(Surv(Mortality_Duration, Mortality) ~Age+BMI+HTN+DOCETAXEL,data=dat,x=TRUE,y=TRUE)
Omics_M <- coxph(Surv(Mortality_Duration, Mortality) ~Age+BMI+HTN+RiskScore_M,data=dat,x=TRUE,y=TRUE)

summary(Dem_M)
summary(Func_M)
summary(Drug_M)
summary(Omics_M)

Dem_HF <- coxph(Surv(HF_Duration, HF) ~Age+BMI+HTN,data=dat,x=TRUE,y=TRUE)
Func_HF <- coxph(Surv(HF_Duration, HF) ~Age+BMI+HTN+IVSd,data=dat,x=TRUE,y=TRUE)
Drug_HF <- coxph(Surv(HF_Duration, HF) ~Age+BMI+HTN+DOCETAXEL,data=dat,x=TRUE,y=TRUE)
Omics_HF <- coxph(Surv(HF_Duration, HF) ~Age+BMI+HTN+RiskScore_HF,data=dat,x=TRUE,y=TRUE)

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

write.csv(aucgraphDem, file = "/path/to/csv")
write.csv(aucgraphFunc, file = "/path/to/csv")
write.csv(aucgraphDrug, file = "/path/to/csv")
write.csv(aucgraphOmics, file = "/path/to/csv")

briergraphDem <- plotBrier(Brier, models = "coxph", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#4682B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphFunc <- plotBrier(Brier, models = "coxph.1", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#AF46B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphDrug <- plotBrier(Brier, models = "coxph.2", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#B47846", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphOmics <- plotBrier(Brier, models = "coxph.3", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#4BB446", lwd = 5, conf.int =TRUE, legend = FALSE)

write.csv(briergraphDem, file = "/path/to/csv")
write.csv(briergraphFunc, file = "/path/to/csv")
write.csv(briergraphDrug, file = "/path/to/csv")
write.csv(briergraphOmics, file = "/path/to/csv")

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

write.csv(aucgraphDemHF, file = "/path/to/csv")
write.csv(aucgraphFuncHF, file = "/path/to/csv")
write.csv(aucgraphDrugHF, file = "/path/to/csv")
write.csv(aucgraphOmicsHF, file = "/path/to/csv")

briergraphDemHF <- plotBrier(BrierHF, models = "coxph", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#4682B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphFuncHF <- plotBrier(BrierHF, models = "coxph.1", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#AF46B4", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphDrugHF <- plotBrier(BrierHF, models = "coxph.2", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#B47846", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphOmicsHF <- plotBrier(BrierHF, models = "coxph.3", which = "score", xlim = c(0,4000), ylim = c(0,.5), col = "#4BB446", lwd = 5, conf.int =TRUE, legend = FALSE)

write.csv(briergraphDemHF, file = "/path/to/csv")
write.csv(briergraphFuncHF, file = "/path/to/csv")
write.csv(briergraphDrugHF, file = "/path/to/csv")
write.csv(briergraphOmicsHF, file = "/path/to/csv")


###Mortalty###
ApparrentCindexM <- pec::cindex(list("COXPH Dem"=coxDemM,
                                    "COXPH Function"=coxFuncM,
                                    "COXPH Drugs"=coxDrugM,
                                    "COXPH Omics"=coxOmicsM),
              formula=Surv(Mortality_Duration, Mortality)~1,data=dat,
              eval.times=seq(0,3939,1), pred.times=seq(0,3939,1))

col = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")

plot(ApparrentCindexM, legend = FALSE, xlim=c(0,4000), ylim=c(0.5,1.0), lwd = 5, col = col)
write.csv(ApparrentCindexM$AppCindex, file = "/path/to/csv")

###HF###
ApparrentCindexHF <- pec::cindex(list("COXPH Dem"=coxDemHF,
                                    "COXPH Function"=coxFuncHF,
                                    "COXPH Drugs"=coxDrugHF,
                                    "COXPH Omics"=coxOmicsHF),
                               formula=Surv(HF_Duration, HF)~1,data=dat,
                               eval.times=seq(0,3754,1), pred.times=seq(0,3754,1))

col = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")

plot(ApparrentCindexHF, legend = FALSE, xlim=c(0,4000), ylim=c(0.5,1.0), lwd = 5, col = col)
write.csv(ApparrentCindexHF$AppCindex, file = "/path/to/csv")
