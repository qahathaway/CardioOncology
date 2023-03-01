# CardioOncology
## Background:
Breast cancer chemotherapy/immunotherapy can be associated with treatment-limiting cardiotoxicity. Radiomics techniques applied to ultrasound, known as ultrasomics, can be used in cardio-oncology to leverage echocardiography for added prognostic value.

## Objectives:
To utilize ultrasomics features collected prior to antineoplastic therapy to enhance prediction of mortality and heart failure (HF) in patients with breast cancer.

## Methods:
Patients were prospectively recruited to the study through the West Virginia University Cancer Institute. 134 of 455 identified patients met the inclusion criteria. Patients were imaged using echocardiography in the parasternal long-axis prior to receiving chemotherapy. All-cause mortality and HF, developed during treatment, were the primary outcomes. 268 features were assessed, grouped into four major classes: demographics (n=20), heart function (n=7), antineoplastic medication (n=17), and ultrasomics (n=224).

## Results:
Examining all features, mortality (n=40) and HF (n=18) revealed a significant association. Of these, 39/40 mortality and 14/18 HF features belonged to ultrasomics. There was an equal distribution of ultrasomics features coming from the interventricular septum (IVS) and posterior wall (PW) with most of these features derived from end-diastole. When developing a risk prediction score for each feature category, ultrasomics features were significantly associated with both mortality (P=0.019, log-rank test) and HF (P<0.0001, log-rank test). Further, only ultrasomics features provided significant improvement over demographic variables when predicting HF (C-Index: 0.93 [0.88-0.98]).

## Conclusions:
With further investigation, a clinical decision support tool could be developed utilizing routinely obtained patient data alongside ultrasomics variables to augment treatment regimens.

![alt text](https://github.com/qahathaway/CardioOncology/blob/main/Figure_1.jpg)

## Figure 1:
Study overview. (A) 134 patients met the inclusion criteria. Echocardiography was performed prior to antineoplastic therapy initiation. Following therapy, the two primary outcomes were measured: all-cause mortality and heart failure. (B) Echocardiographic images were selected in the end-diastolic and end-systolic phase for each patient. Using LifeX software, the interventricular septum and posterior wall were manually traced. This resulted in four unique segments for each patient: IVS-ED, IVS-ES, PW-ED, and PW-ES. Radiomics features were extracted using LifeX software, resulting in the collection of first and second order features and textural-based features. NGLDM = Neighboring Gray-Level Dependence Matrix, GLCM = Gray-Level Co-ccurrence Matrix, GLZLM = Gray-Level Zone Length Matrix, GLRLM = Gray-Level Run Length Matrix.
