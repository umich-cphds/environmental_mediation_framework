library(stats)  
library(epicalc)
library(Hmisc)
library(Epi)
library(mediation)


d1<-read.csv("dataset.csv")
d1<-na.omit(d1) #ensure your data is complete with no missing observations
#order your dataset's columns by: exposures, mediators, follwed by covariates and outcome
#transform your exposures and mediators as you see fit/necessary according normality assumptions


nexp<- #define number of exposure analytes in your dataset
nmed<- #define number of mediators in your dataset
ncovars<- #define number of covariates and outcome variable in your dataset
med.results<-as.data.frame(matrix(nrow=(nexp*nmed),ncol=41))
colnames(med.results)<-c('nobs', 'ACME.C','ACME.C.lo','ACME.C.hi','ACME.C.Pval','ACME.T','ACME.T.lo',
                         'ACME.T.hi','ACME.T.pval','ADE.C','ADE.C.lo','ADE.C.hi','ADE.C.Pval','ADE.T',
                         'ADE.T.lo','ADE.T.hi','ADE.T.pval','PMed.C','PMed.C.lo','PMed.C.hi','PMed.C.pval',
                         'PMed.T','PMed.T.lo','PMed.T.hi','PMed.T.pval','TE','TE.lo','TE.hi','TE.pval',
                         'ACME.avg','ACME.avg.lo','ACME.avg.hi','ACME.avg.pval','ADE.avg','ADE.avg.lo',
                         'ADE.avg.hi','ADE.avg.pval','PMed.avg','PMed.avg.lo','PMed.avg.hi','PMed.avg.pval')

#Loop to conduct pairwise mediation with multiple exposures and mediators
#Loop repeatedly subsets dataset d1 into d2 for individual pairs of exposures and biomarkers

k=1
for(i in 1:nexp){
  for(j in (nexp+1):(nexp+nmed)){
    d2<-d1[,c(i,j,((nexp+nmed+1):(nexp+nmed+ncovars)))]
    set.seed(111)
    #define your covariates as needed
    m<-lm(d2[,2]~d2[,1]+d2$covariate1+d2$covariate2+d2$covariate3+factor(d2$covariate4),data=d2)
    #define your outcome variable below
    y<-lm(d2$outcome~d2[,1]+d2[, 2]+d2$covariate1+d2$covariate2+d2$covariate3+factor(d2$covariate4),data=d2)
    med<-mediate(m,y,sims=2000,treat="d2[, 1]",mediator="d2[, 2]")
    med.results[k,]<-cbind(nobs(y),med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p, med$d1, med$d1.ci[1], 
                           med$d1.ci[2], med$d1.p, med$z0, med$z0.ci[1],med$z0.ci[2], med$z0.p, med$z1, 
                           med$z1.ci[1], med$z1.ci[2], med$z1.p, med$n0, med$n0.ci[1], med$n0.ci[2], 
                           med$n0.p, med$n1, med$n1.ci[1], med$n1.ci[2], med$n1.p, med$tau.coef, med$tau.ci[1], 
                           med$tau.ci[2], med$tau.p, med$d.avg, med$d.avg.ci[1], med$d.avg.ci[2], med$d.avg.p,
                           med$z.avg, med$z.avg.ci[1], med$z.avg.ci[2], med$z.avg.p, med$n.avg, med$n.avg.ci[1], 
                           med$n.avg.ci[2], med$n.avg.p)
    rownames(med.results)[k]<-paste(colnames(d2)[1],colnames(d2)[2],sep='.')  
    print(k) 
    print(rownames(med.results)[k])  
    k=k+1
  }
}
write.csv(med.results,'med.results.csv')






