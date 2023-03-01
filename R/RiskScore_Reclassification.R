
####Regression Predictions####
library(PredictABEL)
library(prodlim)
library(riskRegression)

##Load Dataset##
CardioOnc <- data.frame(read.csv(file = '/path/to/csv'))

##Impute Median Values##
library(mlr)
imputed = impute(CardioOnc, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

# detach mlr package because of 'plotCalibration' conflict
detach("package:mlr", unload = TRUE)

fgrDemM <- FGR(Hist(Mortality_Duration, Mortality) ~Demographics-Features-Here, data=final, cause=1)

fgrDemHF <- FGR(Hist(HF_Duration, HF) ~Demographics-Features-Here, data=final, cause=1)

fgrFuncM <- FGR(Hist(Mortality_Duration, Mortality) ~Functional-Features-Here, data=final, cause=1)

fgrFuncHF <- FGR(Hist(HF_Duration, HF) ~Functional-Features-Here, data=final, cause=1)

fgrDrugM <- FGR(Hist(Mortality_Duration, Mortality) ~Therapy-Features-Here, data=final, cause=1)

fgrDrugHF <- FGR(Hist(HF_Duration, HF) ~Therapy-Features-Here, data=final, cause=1)

fgrOmicsM <- FGR(Hist(Mortality_Duration, Mortality) ~Ultrasomics-Features-Here, data=final, cause=1)

fgrOmicsHF <- FGR(Hist(HF_Duration, HF) ~Ultrasomics-Features-Here, data=final, cause=1)

risk_score_DemM <- predictRisk(fgrDemM,times=c(3939),newdata=final)
risk_score_FuncM <- predictRisk(fgrFuncM,times=c(3939),newdata=final)
risk_score_DrugM <- predictRisk(fgrDrugM,times=c(3939),newdata=final)
risk_score_OmicsM <- predictRisk(fgrOmicsM,times=c(3939),newdata=final)

risk_score_DemHF <- predictRisk(fgrDemHF,times=c(3754),newdata=final)
risk_score_FuncHF <- predictRisk(fgrFuncHF,times=c(3754),newdata=final)
risk_score_DrugHF <- predictRisk(fgrDrugHF,times=c(3754),newdata=final)
risk_score_OmicsHF <- predictRisk(fgrOmicsHF,times=c(3754),newdata=final)

df_list <- list(risk_score_DemM, risk_score_FuncM, risk_score_DrugM, risk_score_OmicsM)
df <- data.frame(df_list)
write.csv(df, file = "/path/to/csv")

df_listHF <- list(risk_score_DemHF, risk_score_FuncHF, risk_score_DrugHF, risk_score_OmicsHF)
dfHF <- data.frame(df_listHF)
write.csv(dfHF, file = "/path/to/csv")

###Sample 1###
# specify column number of outcome variable
cOutcome <- 1
# specify column numbers of non-genetic predictors
cNonGenPred <- c(3)
# specify column numbers of non-genetic predictors that are categorical
cNonGenPredCat <- c(0)
# specify column numbers of genetic predictors
cGenPred <- c(0)
# specify column numbers of genetic predictors that are categorical
cGenPredCat <- c(0)
# fit logistic regression model
riskmodel1 <- fitLogRegModel(data=final, cOutcome=cOutcome,
                            cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
                            cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
summary(riskmodel1)

###Sample 2###
# specify column numbers of non-genetic predictors
cNonGenPred <- c(6)

# fit logistic regression model
riskmodel2 <- fitLogRegModel(data=final, cOutcome=cOutcome,
                            cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
                            cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
summary(riskmodel2)


# obtain multivariate OR(95% CI) for all predictors of the fitted model
ORmultivariate(riskModel=riskmodel1)
ORmultivariate(riskModel=riskmodel2)

# obtain predicted risks
predRisk1 <- predRisk(riskmodel1)
predRisk2 <- predRisk(riskmodel2)

# specify cutoff values for risk categories
cutoff <- c(0.0,0.15,1.0)
#cutoff <- c(0, 0.27, 0.36, 0.78, 1)
#compute reclassification measures
reclassification(data=final, cOutcome=cOutcome,
                 predrisk1=predRisk1, predrisk2=predRisk2, cutoff=cutoff)




# compute calibration measures and produce calibration plot
plotCalibration(data=final, cOutcome=cOutcome, predRisk=predRisk1, groups = 100)
plotCalibration(data=final, cOutcome=cOutcome, predRisk=predRisk2, groups = 100)


# specify labels for the groups without and with the outcome of interest
labels <- c("No Mortality", "Mortality")
# produce discrimination box plot
plotDiscriminationBox(data=final, cOutcome=cOutcome, predrisk=predRisk1,
                      labels=labels)
plotDiscriminationBox(data=final, cOutcome=cOutcome, predrisk=predRisk2,
                      labels=labels)
