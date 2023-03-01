
####Regression Predictions####
library(PredictABEL)
library(prodlim)
library(riskRegression)

##Load Dataset##
MESA <- data.frame(read.csv(file = '/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/ALL_PRE.csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

# detach mlr package because of 'plotCalibration' conflict
detach("package:mlr", unload = TRUE)

fgrDemM <- FGR(Hist(Mortality_Duration, Mortality) ~Age+BMI+HTN+HLD+DM+COPD+ASA+BB+Statin+
                   ACEi_ARB_Entresto+CCB+Diuretic+insulin+Metfromin+NSAID+smokinghx+currentsmoking+
                   Alcohol+Race, data=final, cause=1)

fgrDemHF <- FGR(Hist(HF_Duration, HF) ~Age+BMI+HTN+HLD+DM+COPD+ASA+BB+Statin+
                   ACEi_ARB_Entresto+CCB+Diuretic+insulin+Metfromin+NSAID+smokinghx+currentsmoking+
                   Alcohol+Race, data=final, cause=1)

fgrFuncM <- FGR(Hist(Mortality_Duration, Mortality) ~EF.Bi.Plane..Or.Visually.Estimated+IVSd+LVPWd+
                   MV.E+MV.A+Lat.E.+Med.E., data=final, cause=1)

fgrFuncHF <- FGR(Hist(HF_Duration, HF) ~EF.Bi.Plane..Or.Visually.Estimated+IVSd+LVPWd+
                   MV.E+MV.A+Lat.E.+Med.E., data=final, cause=1)

fgrDrugM <- FGR(Hist(Mortality_Duration, Mortality) ~ATEZOLIZUMAB+CAPECITABINE+CARBOPLATIN+CELECOXIB+
                   CYCLOPHOSPHAMIDE+DOCETAXEL+DOXORUBICIN+ERIBULIN+
                   EVEROLIMUS+GEMCITABINE+METHOTREXATE+NERATINIB+
                   PACLITAXEL+PALBOCICLIB+PEMBROLIZUMAB+PERTUZUMAB+
                   TRASTUZUMAB, data=final, cause=1)

fgrDrugHF <- FGR(Hist(HF_Duration, HF) ~ATEZOLIZUMAB+CAPECITABINE+CARBOPLATIN+CELECOXIB+
                   CYCLOPHOSPHAMIDE+DOCETAXEL+DOXORUBICIN+ERIBULIN+
                   EVEROLIMUS+GEMCITABINE+METHOTREXATE+NERATINIB+
                   PACLITAXEL+PALBOCICLIB+PEMBROLIZUMAB+PERTUZUMAB+
                   TRASTUZUMAB, data=final, cause=1)

fgrOmicsM <- FGR(Hist(Mortality_Duration, Mortality) ~PW_PRE_CONVENTIONAL_min_ED+
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
                     PW_PRE_GLRLM_RP_ED, data=final, cause=1)

fgrOmicsHF <- FGR(Hist(HF_Duration, HF) ~PW_PRE_CONVENTIONAL_min_ED+
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
                      PW_PRE_GLRLM_HGRE_ES, data=final, cause=1)

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
write.csv(df, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/RiskScore_M.csv")

df_listHF <- list(risk_score_DemHF, risk_score_FuncHF, risk_score_DrugHF, risk_score_OmicsHF)
dfHF <- data.frame(df_listHF)
write.csv(dfHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/RiskScore_HF.csv")



##Load RiskScore##
MESA <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/RiskScore_HF.csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

# detach mlr package because of 'plotCalibration' conflict
detach("package:mlr", unload = TRUE)


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

# specify range of y-axis
rangeyaxis <- c(0,1)
# specify labels of the predictiveness curves
labels <- c("B", "Framingham")
# produce predictiveness curves
plotPredictivenessCurve(predrisk=cbind(predRisk1,predRisk2),
                        rangeyaxis=rangeyaxis, labels=labels)


# specify the size of each interval
interval <- .05
# specify label of x-axis
xlabel <- "Predicted risk"
# specify label of y-axis
ylabel <- "Percentage"
# specify range of x-axis
xrange <- c(0,1)
# specify range of y-axis
yrange <- c(0,40)
# specify title for the plot
maintitle <- "Distribution of predicted risks"
# specify labels
labels <- c("No Mortality", "Mortality")
# produce risk distribution plot
plotRiskDistribution(data=final, cOutcome=cOutcome,
                     risks=predRisk1, interval=interval, plottitle=maintitle, rangexaxis=xrange,
                     rangeyaxis=yrange, xlabel=xlabel, ylabel=ylabel, labels=labels)

# specify range of x-axis
rangexaxis <- c(0,12)
# specify range of y-axis
rangeyaxis <- c(0,1)
# specify label of x-axis
xlabel <- "Risk score"
# specify label of y-axis
ylabel <- "Predicted risk"
# specify title for the plot
plottitle <- "Risk score versus predicted risk"
# function to compute unweighted genetic risk scores
riskScore1 <- riskScore(weights=riskmodel1, data=final,
                       cGenPreds=cGenPred, Type="unweighted")

# produce risk score-predicted risk plot
plotRiskscorePredrisk(data=final, riskScore=riskScore1, predRisk=predRisk1,
                      plottitle=plottitle, xlabel=xlabel, ylabel=ylabel, rangexaxis=rangexaxis,
                      rangeyaxis=rangeyaxis)


# specify cutoff values for risk categories
cutoff <- c(0,.10,.20,.30,.40,.50,.60,.70,.80,.90,1)
# compute reclassification measures
reclassification(data=final, cOutcome=cOutcome,
                 predrisk1=predRisk1, predrisk2=predRisk2, cutoff=cutoff)


#--- sample data (pbc in survival package) ---
D=subset(pbc, select=c("time","status","age","albumin","edema","protime","bili"))
D$status=as.numeric(D$status==2)
D=D[!is.na(apply(D,1,mean)),] ; dim(D)
mydata=D[1:100,]





####IDI Survival####
library(survIDINRI)

final <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/C-Statistic/NEW_PredScore.csv'))

library(mlr)
imputed = impute(final, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final1 <- as.data.frame(imputed$data)
final1$event=as.integer(final1$event)
final1$duration=as.integer(final1$duration)
detach("package:mlr", unload = TRUE)

covs0<-as.matrix(final1[5])
covs1<-as.matrix(final1[6])
#--- inference ---
t0=1140

x<-IDI.INF(final1[,2:1], covs0, covs1, t0, npert=10)
#--- results ---
IDI.INF.OUT(x)
#--- Graphical presentaion of the estimates###
IDI.INF.GRAPH(x)
