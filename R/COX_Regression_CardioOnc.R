
##Load Packages##
library(caret)
library(survival)
library(survminer)

##Load Dataset##
MESA <- data.frame(read.csv(
  file = '/path/to/csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

##Univariate##

covariates <- c(colnames(final)[6:273])

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(HF_Duration, HF)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = final)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res1 <- as.data.frame(res)
print(res1[order(res1$p.value, decreasing = FALSE), ]   )

write.csv(res1[order(res1$p.value, decreasing = FALSE), ], file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Univariate_HF.csv")

##Multivariate##
library(survcomp)

cOmicsMort <- concordance.index(x=final$Ultraosmics-Features-Here, method="noether",
                              surv.time=final$Mortality_Duration, surv.event=final$Mortality)

cOmicsHF <- concordance.index(x=final$Ultraosmics-Features-Here, method="noether",
                                  surv.time=final$HF_Duration, surv.event=final$HF)

cDemMort <- concordance.index(x=final$Demographics-Features-Here,
                            method="noether", surv.time=final$Mortality_Duration, surv.event=final$Mortality)

cDemHF <- concordance.index(x=final$Demographics-Features-Here,
                                method="noether", surv.time=final$HF_Duration, surv.event=final$HF)

cFuncMort <- concordance.index(x=final$Functional-Features-Here, method="noether",
                             surv.time=final$Mortality_Duration, surv.event=final$Mortality)

cFuncHF <- concordance.index(x=final$Functional-Features-Here, surv.event=final$HF)

cDrugMort <- concordance.index(x=final$Therapy-Features-Here, method="noether",
                             surv.time=final$Mortality_Duration, surv.event=final$Mortality)

cDrugHF <- concordance.index(x=final$Therapy-Features-Here, method="noether",
                             surv.time=final$HF_Duration, surv.event=final$HF)

cindex.comp(cDrugMort, cDemMort)

summary(cDrugMort)
concordance(cDrugMort)
anova(cDrugMort, cDemMort)


##Plot the baseline survival function##
fit <- surv_fit(Surv(HF_Duration, HF) ~Prob_Omics,
                data = final)

#ggsurvplot(fit, data = final, pval = TRUE, break.time.by = 500)

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  fun="event",
#  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",
  ylab = "Incidence",# customize X axis label.
  fontsize = 2,
  break.time.by = 500,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = TRUE,
  risk.table.font = 3,
  risk.table.title = "Number at Risk (%)",
  risk.table.pos = "out",
  cumevents.title = "Cumulative Events",
  cumevents = FALSE,
  cumevents.font = 3,
  font.x = 12,
  font.y = 12,
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
#  surv.median.line = "hv",  # add the median survival pointer.
  ylim = c(0, 1),
  legend.labs = 
    c("No Mortality", "Mortality"),    # change legend labels.
  palette = 
    c("black", "red") # custom color palettes.
)
