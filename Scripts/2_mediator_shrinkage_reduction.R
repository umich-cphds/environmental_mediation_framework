library(bama)
library(devtools)
install_github("oliverychen/PDM")
library(PDM)

d1<-read.csv("dataset.csv")

d1<-na.omit(d1) #ensure your is complete with no missing observations
#order your dataset's columns by: exposures, mediators, follwed by covariates and outcome
#transform your exposures and mediators as you see fit/necessary according to normality assumptions


#creating data frame to deposit results
nexp<- #define the number of exposure analytes in your dataset
nmed<- #define the number of mediators in your dataset
ncovars<- #define the number of covariates and outcome variable in your dataset
  
#read in a crosswalk data frame for differentiating mediators into groups
#cross walk should have at least a column for mediator group identity, and a column for the name of each mediator
crosswalk <- read.csv("crosswalk.csv", header=T, strip.white=T, stringsAsFactors=F)


#subset crosswalk for the groups you want to create
#this could be biological pathways (e.g. inflammation or oxidative stress), or gene grops, etc.
crosswalk.group1<-subset(crosswalk,crosswalk$group_name=="group 1")
crosswalk.group2<-subset(crosswalk,crosswalk$group_name=="group 2")
crosswalk.group3<-subset(crosswalk,crosswalk$group_name=="group 3")
crosswalk.group4<-subset(crosswalk,crosswalk$group_name=="group 4")
crosswalk.group5<-subset(crosswalk,crosswalk$group_name=="group 5")
crosswalk.group6<-subset(crosswalk,crosswalk$group_name=="group 6")
crosswalk.group7<-subset(crosswalk,crosswalk$group_name=="group 7")

#subset dataset so that each subset only contains mediators from one group
#the bama application function will require your columns to be ordered by exposures varaibles, then mediators, then covariates and outcome
d.group1<-cbind(d1[c("exposure variables")], d1[,colnames(d1) %in% unique(crosswalk.group1$mediator_name)],d1[c("covariates,outcome")])
d.group2<-cbind(d1[c("exposure variables")], d1[,colnames(d1) %in% unique(crosswalk.group2$mediator_name)],d1[c("covariates,outcome")])
d.group3<-cbind(d1[c("exposure variables")], d1[,colnames(d1) %in% unique(crosswalk.group3$mediator_name)],d1[c("covariates,outcome")])
d.group4<-cbind(d1[c("exposure variables")], d1[,colnames(d1) %in% unique(crosswalk.group4$mediator_name)],d1[c("covariates,outcome")])
d.group5<-cbind(d1[c("exposure variables")], d1[,colnames(d1) %in% unique(crosswalk.group5$mediator_name)],d1[c("covariates,outcome")])
d.group6<-cbind(d1[c("exposure variables")], d1[,colnames(d1) %in% unique(crosswalk.group6$mediator_name)],d1[c("covariates,outcome")])
d.group7<-cbind(d1[c("exposure variables")], d1[,colnames(d1) %in% unique(crosswalk.group7$mediator_name)],d1[c("covariates,outcome")])


################################################################################
### Mediator shrinkage using Bayesian mediation analysis (Song et al. 2020) ####
################################################################################

bama.apply<-function(dataset){
  d2<-dataset
  j=1
  k=ncol(d2[,(nexp+1):(ncol(d2)-ncovars)]) 
  l=ncol(d2[,(nexp+1):(ncol(d2)-ncovars)]) 
  results<-as.data.frame(matrix(nrow=(nexp*l),ncol=3)) 
  results[,2]<-rep(colnames(d2)[(nexp+1):(l+nexp)],nexp) 
  for (i in 1:nexp){
    Y<-as.numeric(d2$outcome) #replace with your outcome variable
    A<-d2[,i] 
    M<-as.matrix(d2[,(nexp+1):(l+nexp)])
    C<-as.matrix(d2[,(nexp+l+1):(nexp+l+ncovars)])
    beta.m<-rep(0,ncol(M))
    alpha.a<-rep(0,ncol(M))
    set.seed(111)
    bama.out<-bama(Y,A,M,C,beta.m,alpha.a,burnin=10000,ndraws=1000)
    pips<-summary(bama.out)[,4]
    results[j:k,3]<-pips
    results[j:k,1]<-colnames(d2)[i]
    print(results[k,])
    j=j+l
    k=k+l
  }
  return(results)
}
pips.allmediators<-bama.apply(d1)
pips.group1<-bama.apply(d.group1)
pips.group2<-bama.apply(d.group2)
pips.group3<-bama.apply(d.group3)
pips.group4<-bama.apply(d.group4)
pips.group5<-bama.apply(d.group5)
pips.group6<-bama.apply(d.group6)
pips.group7<-bama.apply(d.group7)

bama.results<-rbind(pips.group1,pips.group2,pips.group3,pips.group4,pips.group5,pips.group6,pips.group7)

write.csv(bama.grouped.results,'bama.grouped.results.csv')


#################################################################################
### High-dimensional multivariate mediation with mediator dimension reduction ###
### using population value decomposition (Chen et al. 2018)                   ###  
#################################################################################

# NOTE: this current application is for a setting where the number of mediators is less than
# the number of individuals in sample size (p<n scenario)
# if your data have more mediators than the sample size, refer to script "gPVD.R" for example code
# to implement an intial matrix decomposition step to reduce the dimensionality of your mediator matrix
# prior to running this population value decomposition code below

hdmm.apply<-function(dataset){
  d2<-dataset
  j=1
  k=ncol(d2[,(nexp+1):(ncol(d2)-ncovars)])
  results<-as.data.frame(matrix(nrow=nexp,ncol=42))
  results[,1]<-colnames(d2)[1:nexp]
  for (i in 1:nexp){
    Y<-as.numeric(d2$outcome) #replace with your outcome variable
    A<-d2[,i] 
    m1<-as.matrix(d2[,(nexp+1):(k+nexp)])
    output <- PDM_1(x=A,y=Y,m=m1,imax=100, tol=10^-{5}, theta=rep(1,5),w1=rep(1,ncol(m1)), interval=10^6, step=10^4)
    est_w <- output$w1
    dm1 <- m1%*%as.matrix(est_w, ncol = 1)
    m<-lm(dm1~A+d2$covariate1+d2$covariate2+d2$covariate3+d2$covariate4,data=d2) #replace with your covariates
    y<-lm(Y~A+dm1+d2$covariate1+d2$covariate2+d2$covariate3+d2$covariate4,data=d2) #replace with your covariates
    set.seed(111)
    med<-mediate(m,y,sims=2000,treat="A",mediator="dm1")
    results[j,2:42]<-cbind(nobs(y),med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p, med$d1, med$d1.ci[1], med$d1.ci[2], med$d1.p,
                           med$z0, med$z0.ci[1],med$z0.ci[2], med$z0.p, med$z1, med$z1.ci[1], med$z1.ci[2], med$z1.p, 
                           med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p, med$n1, med$n1.ci[1], med$n1.ci[2], med$n1.p, 
                           med$tau.coef, med$tau.ci[1], med$tau.ci[2], med$tau.p, med$d.avg, med$d.avg.ci[1], med$d.avg.ci[2], 
                           med$d.avg.p,med$z.avg, med$z.avg.ci[1], med$z.avg.ci[2], med$z.avg.p, med$n.avg, med$n.avg.ci[1],
                           med$n.avg.ci[2], med$n.avg.p)
    print(results[j,1])
    j=j+1 
  }
  return(results)
}
hdmm.allmediators<-hdmm.apply(d1)
hdmm.group1<-hdmm.apply(d.group1)
hdmm.group2<-hdmm.apply(d.group2)
hdmm.group3<-hdmm.apply(d.group3)
hdmm.group4<-hdmm.apply(d.group4)
hdmm.group5<-hdmm.apply(d.group5)
hdmm.group6<-hdmm.apply(d.group6)
hdmm.group7<-hdmm.apply(d.group7)

hdmm.all<-rbind(hdmm.allmediators,hdmm.group1,hdmm.group2,hdmm.group3,hdmm.group4,hdmm.group5,hdmm.group6,hdmm.group7,hdmm.group8)
colnames(hdmm.all)<-c('exp', 'nobs', 'ACME.C','ACME.C.lo','ACME.C.hi','ACME.C.Pval','ACME.T','ACME.T.lo','ACME.T.hi',
                         'ACME.T.pval','ADE.C','ADE.C.lo','ADE.C.hi','ADE.C.Pval','ADE.T', 'ADE.T.lo','ADE.T.hi','ADE.T.pval',
                         'PMed.C','PMed.C.lo','PMed.C.hi','PMed.C.pval','PMed.T','PMed.T.lo','PMed.T.hi','PMed.T.pval','TE',
                         'TE.lo','TE.hi','TE.pval','ACME.avg','ACME.avg.lo','ACME.avg.hi','ACME.avg.pval','ADE.avg',
                         'ADE.avg.lo','ADE.avg.hi','ADE.avg.pval','PMed.avg','PMed.avg.lo','PMed.avg.hi','PMed.avg.pval')
hdmm.all$group<-c(rep("allmediators",nexp), rep("group1",nexp),rep("group2",nexp),rep("group3",nexp),rep("group4",nexp),
                     rep("group5",nexp),rep("group6",nexp),rep("group7",nexp))

write.csv(hdmm.all,'hdmm.results.csv')


#################################################################################
### HIMA: high-dimensional mediation analysis (Zhang et al. 2016)             ###
#################################################################################

## 
install.packages("HIMA")
library(HIMA)

hima.apply<-function(dataset){
  results<-as.data.frame(matrix(nrow=nmed*nexp,ncol=9))
  j=1
  for(i in 1:nexp){
    d2<-dataset
    k=ncol(d2[,(nexp+1):(ncol(d2)-ncovars)])
    Xvar<-d2[,i] 
    Yvar<-d2$outcome #change with your outcome variable
    C<-as.matrix(d2[,(nexp+k+1):(nexp+k+ncovars)])
    Mvars<-as.matrix(d2[,(nexp+1):(k+nexp)])
    hima.model<-hima(Xvar,Yvar,Mvars,COV.XM=C,COV.MY = C,family="gaussian")
    print(hima.model)
    if(length(hima.model)>0){
      results[j:(j+dim(hima.model)[1]-1),1]<-colnames(d2)[i]
      results[j:(j+dim(hima.model)[1]-1),2]<-rownames(hima.model)
      results[j:(j+dim(hima.model)[1]-1),3:9]<-as.matrix(hima.model)
      print(results[j:(j+dim(hima.model)[1]-1),1])
      j=j+dim(hima.model)[1]
    }else{
      print("NULL")
    }
  }
  colnames(results)<-c("alpha","beta","gamma","alpha*beta","% total effect","adjusted.p","BH.FDR","mediator","exp")
  return(results)
}

hima.allmediators<-hima.apply(d1)
hima.group1<-hima.apply(d.group1)
hima.group2<-hima.apply(d.group2)
hima.group3<-hima.apply(d.group3)
hima.group4<-hima.apply(d.group4)
hima.group5<-hima.apply(d.group5)
hima.group6<-hima.apply(d.group6)
hima.group7<-hima.apply(d.group7)

write.csv(hima.all,"hima.allmediators.csv")
write.csv(hima.group1,"hima.group1.csv")
write.csv(hima.group2,"hima.group2.csv")
write.csv(hima.group3,"hima.group3.csv")
write.csv(hima.group4,"hima.group4.csv")
write.csv(hima.group5,"hima.group5.csv")
write.csv(hima.group6,"hima.group6.csv")
write.csv(hima.group7,"hima.group7.csv")


