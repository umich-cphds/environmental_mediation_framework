# environmental_mediation_framework

#There are four main scripts in this environmental mediaiton framework based on the approaches in the pipeline: 

#Approach 1a: 1a_pairwise_mediation.R
#expected run time is approximately 3-5 minutes per mediation model

#Approach 1b: 1b_exposure_reduction.R
#expected run time is approximately 3-5 minutes per mediation model

#Approach 2: 2_mediator_shrinkage_reduction.R and 2_apply_pathlasso.R
#expected run time is approximately 10-20 minutes per mediation model

#2_apply_pathlasso.R is located in the subfolder "Pathway_lasso" and will require source functions also located in that folder

#Approach 2 is chiefly focused on a setting where number of mediators is less than total sample size

#If the user has a dataset where the number of mediators are greater than sample size, refer to script: gPVD.R for example code on implementing a preliminary matrix decomposition step prior to high-dimensional mediation analysis

#Relevant packages needed for analyses are identified in each R script 
#Install time for packages is approximately 1-3 minutes per package
#Expected output from analytical framework will be data files with single mediation model results per row

#Sample_Data folder contains sample data file ("data_example.csv") available for demo and the data_dictionary.xlsx file to define variables
#data_example.csv contains true observed measurements of exposure and mediator variables and synthetic outcome and covariate variables based on the marginal association between the phthalate risk score and gestational age at delivery
#In total, the data_eample.csv includes: 42 exposures (including 4 environmental risk scores), 61 mediators, 6 covariates, and 1 outcome variable




