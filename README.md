# environmental_mediation_pipeline

#There are three main scripts in this environmental mediaiton pipeline based on the steps in the pipeline: 

#Step 1: 1_pairwise_mediation.R

#Step 2: 2_mediator_shrinkage_reduction.R

#Step 3: 3_exposure_reduction.R

#Step 3 requires the user to return to script 2_mediator_srinkage_reduction after implementing exposure reduction

#Step 2 is chiefly focused on a setting where number of mediators is less than total sample size

#If the user has a dataset where the number of mediators are greater than sample size, refer to script: gPVD.R for example code on implementing a preliminary matrix decomposition step prior to high-dimensaional mediation analysis