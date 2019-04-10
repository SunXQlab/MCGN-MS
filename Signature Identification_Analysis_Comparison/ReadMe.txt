--- Multicellular gene network analysis identifies a macrophage-related gene signature predictive of therapeutic response and prognosis of gliomas

The codes were implemented in R (version 3.5.1).  The codes included two sections: (1) Multicellular Gene Network; and (2) Signature Identification_Analysis_Comparison.

(1) Multicellular Gene Network:

This file folder contains R codes of selecting DEGs and constructing PCC-based multicellular gene networks. (for Figures 2 and Figures S1-S7)

The gene expression data of mouse is saved in 'Gene_expression_DATA.csv' or 'Gene_expression_DATA.RData'. 

(1.1)  selected.R           --------       Analyze DEGs and select top 50 DEGs in TCs and TAMs respectively.

(1.2) pccEPTC.R           ---------       Construct PCC-based correlation network for DEGs in TCs of Ep group.

(1.3) pccEPTAM.R        ---------      Construct PCC-based correlation network for DEGs in TAMs of Ep group.

(1.4) pccEP_TAM_TC.R ---------      Construct PCC-based correlation network for DEGs between TAMs and TCs of Ep group.

(1.5) pccREBTC.R          ---------     Construct PCC-based correlation network for DEGs in TCs of Reb group.

(1.6) pccREBTAM.R       ---------     Construct PCC-based correlation network for DEGs in TAMs of Reb group.

(1.7) pccREB_TAM_TC.R ---------    Construct PCC-based correlation network for DEGs between TAMs and TCs of Reb group.

(1.8) compccdifPCC0.05.R --------   Construct differential network using network perturbation method.

(1.9) updown.R                ---------    Analysis of correlation-gain or -loss of the edges in the networks. 


(2) Signature Identification_Analysis_Comparison:

The clinical information of patient from TCGA and CGGA are saved in .csv files in this folder. The RNA-seq gene expression data of patients are in too large size (>25M) to be uploaded here, which are needed to be downloaded from TCGA (GBM and LGG) and CGGA and saved under the working path.

(2.1) Macrophage_Signature.R     ---------  Train a macrophage-related gene signature from CGGA set and validate on TCGA set; K-M analysis; plot ROC and compute AUC.  (for Figures 3-4)
 
(2.2) ROC_AUC_RandomGeneSet.R -------   Evaluating statistical significance of prognostic accuracy of the macrophage-related gene signature using a bootstrapping approach. (for Figure 5)

(2.3) Comparison with LASSO.R     -------    Build LASSO Cox model and compare the robustness with the macrophage-related gene signature.  (for Figure 6A-C, Figure S8A-B)

(2.4) Comparison with coorelation network.R     -------  Identify correlation network-based gene signature and compare the robustness with the macrophage-related gene signature.  (for Figure 6D-E, Figure S8C-D)

(2.5) Comparison with IGF1-PI3K pathway.R   --------  Build a Cox model for IGF1-PI3K pathways-based gene signature.

(2.6) KM_ROC_TargetedTherapy.R   ---------  Plot K-M survival curves for the patients who received targeted therapies. (for Figure 7A-B)

(2.7)  ROC_comparison_TCGA_Targeted_Therapy.R  ----  Compare ROCs of 4 signatures for the patients who received targeted therapies. (for Figure 7C-D)

(2.8) Multi-variate COX Model.R     --------   Multivariate COX regression analysis of clinicopathologic factors and four gene signatures for predicting overall survival and 5-year survival in the validation set.  (for Table 1)

(2.9) Stratified K-M Analysis.R        ---------  Prognostic significance of the macrophage-related gene signature in the stratified cohorts. (for Figures 8, 10)

(2.10) Combined_Signature.R         ---------  Compare the prognostic accuracy of the macrophage-related gene signature with clinicopathological risk factors and other existing gene signatures or their combination. (for Figure 9)

