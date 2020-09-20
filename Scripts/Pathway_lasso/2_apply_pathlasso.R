################################################################################
################ Pathway Lasso application (Zhao and Luo 2016) #################
################################################################################


source("functions_1.R")
source("functions_2_ADMM_adp.R")
source("functions_3_VSS_adp.R")

d1<-read.csv("dataset.csv")
d1<-na.omit(d1) #ensure your is complete with no missing observations
#order dataset as exposures, mediators, then covariates
#define:
nexp<- #define the number of exposure analytes in your dataset
nmed<- #define the number of mediators in your dataset
ncovars<- #define the number of covariates and outcome variable in your dataset
  
#pathlass.apply function asks you to input the main exposure of interest
pathlasso.apply<-function(exposure){
A<-exposure
M<-d1[,((nmed):(nexp+nmed+1))]
Y<-d1$outcome #repalce with your outcome variable
q=nmed
# approximate covariance matrix
Sigma10<-diag(rep(1,q))
Sigma20<-matrix(1,1,1)

m.A<-mean(A)
m.M<-apply(M,2,mean)
m.Y<-mean(Y)
sd.A<-sd(A)
sd.M<-apply(M,2,sd)
sd.Y<-sd(Y)

A<-scale(A) 
M<-scale(M) 
Y<-scale(Y) 

# Pathway Lasso method parameters
phi<-2
gamma<-0           # adpative lasso parameter
rho<-1             # ADMM parameter
max.itr<-10000     # change if needed for model convergence 
tol<-1e-6
thred<-1e-6
thred2<-1e-3
omega.p<-0.1      # omega = omega.p*lambda

# lambda for initial estimates
lambda.ini<-0.005
lambda<-c(10^c(seq(-5,-3,length.out=5),seq(-3,0,length.out=21)[-1],seq(0,2,length.out=11)[-1],seq(2,4,length.out=6)[-1]))

AB.est=B.est=A.est<-matrix(NA,q,length(lambda))

for (i in 1:length(lambda)) {

out<-mediation_net_ADMM_NC(A,M,Y,lambda=lambda[i],omega=omega.p*lambda[i],phi=phi,Phi1=diag(c(0,rep(1,q))),
                               Phi2=diag(rep(1,q+1)),rho=rho,rho.increase=FALSE,tol=tol,max.itr=max.itr,thred=thred,
                               Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE)

A.est[,i]<-out$A*(sd.M/sd.A) #exposure to mediator selection
B.est[,i]<-out$B*(sd.Y/sd.M) #mediator to outcome selection
AB.est[,i]<-A.est[,i]*B.est[,i] #mediation effect selection
print(AB.est[,i])
}

rownames(A.est)<-rownames(B.est)<-rownames(AB.est)<-colnames(d1[,((nmed):(nexp+nmed+1))])


colnames(A.est)<-lambda
colnames(B.est)<-lambda
colnames(AB.est)<-lambda

results<-list(A.est,B.est,AB.est)
return(results)

}

pathlasso_results<- pathlasso.apply(d1$exposure_variable_of_interest)

write.csv(pathlasso_results[[1]],"pathwaylasso_Alpha_selection_est.csv") #exposure to mediator selection
write.csv(pathlasso_results[[2]],"pathwaylasso_Beta_selection_est.csv") #mediator to outcome selection
write.csv(pathlasso_results[[3]],"pathwaylasso_Alpha-Beta_selection_est.csv") #mediation effect selection

