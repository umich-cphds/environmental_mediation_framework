# environmental_mediation_framework

#There are three main scripts in this environmental mediaiton framework based on the approaches in the pipeline: 

#Approach 1a: 1a_pairwise_mediation.R
##expected run time is approximately 3-5 minutes per mediation model

#Approach 1b: 1b_exposure_reduction.R
##expected run time is approximately 3-5 minutes per mediation model

#Approach 2: 2_mediator_shrinkage_reduction.R
##expected run time is approximately 10-20 minutes per mediation model

#Approach 2 is chiefly focused on a setting where number of mediators is less than total sample size

#If the user has a dataset where the number of mediators are greater than sample size, refer to script: gPVD.R for example code on implementing a preliminary matrix decomposition step prior to high-dimensional mediation analysis

#Relevant packages needed for analyses are identified in each R script 
#Install time for packages is approximately 1-3 minutes per package
#Expected output from analytical framework will be data files with single mediation model results per row

#sample data file ("sample.data.csv") is available for demo
##sample data frame columns in order include the following: 10 exposures, 2000 mediators, and outcome variable
##If desired, 2000 mediators can be sequentially divided into 50 synthetic biological groups of 20 mediators each



