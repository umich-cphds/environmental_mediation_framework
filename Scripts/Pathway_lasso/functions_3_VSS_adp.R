#######################################################################################
# Consistent selection of tuning parameters via variable selection stability
# Sun, Wang and Fang, Journal of Machine Learning Research (2013)
#######################################################################################


#######################################################################################
# Cohen's kappa coefficient
cohen.kappa<-function(S1,S2)
{
  if(length(S1)!=length(S2))
  {
    stop("Cardinality of two sets are different!")
  }else
  {
    p<-length(S1)
    
    n11<-length(which(S1==1&S2==1))
    n12<-length(which(S1==1&S2==0))
    n21<-length(which(S1==0&S2==1))
    n22<-length(which(S1==0&S2==0))
    
    re.table<-matrix(c(n11,n21,n12,n22),2,2)
    colnames(re.table)<-paste0("S2_",c("1","0"))
    rownames(re.table)<-paste0("S1_",c("1","0"))
    
    if(n11==p|n22==p)
    {
      kp<--1
    }else
    {
      pa<-(n11+n22)/p
      pe<-(n11+n12)*(n11+n21)/p^2+(n12+n22)*(n21+n22)/p^2
      
      kp<-(pa-pe)/(1-pe) 
    }
    
    return(list(table=re.table,kappa=kp))
  }
}

#######################################################################################
# Cauculate variable selection stability for each iteration
# data1<-list(Z=Z1,M=M1,R=R1)
# data2<-list(Z=Z2,M=M2,R=R2)
mediation_net_ADMM_NC_VSS<-function(data1,data2,zero.cutoff=1e-3,
                                    lambda=1,omega=0,phi=1,data1_Phi1=NULL,data1_Phi2=NULL,data2_Phi1=NULL,data2_Phi2=NULL,
                                    rho=1,rho.increase=FALSE,tol=1e-10,max.itr=10000,thred=1e-10,Sigma1=NULL,Sigma2=NULL,
                                    trace=FALSE,Theta0=NULL,D0=NULL,alpha0=NULL,beta0=NULL)
{
  #######################################################
  # data1
  Z1<-matrix(data1$Z,ncol=1)
  M1<-matrix(data1$M,ncol=ncol(data1$M))
  R1<-matrix(data1$R,ncol=1)
  
  re1<-mediation_net_ADMM_NC(Z1,M1,R1,lambda=lambda,omega=omega,phi=phi,Phi1=data1_Phi1,Phi2=data1_Phi2,
                             rho=rho,rho.increase=rho.increase,tol=tol,max.itr=max.itr,thred=thred,Sigma1=Sigma1,Sigma2=Sigma2,
                             trace=trace,Theta0=Theta0,D0=D0,alpha0=alpha0,beta0=beta0)
  
  # data2
  Z2<-matrix(data2$Z,ncol=1)
  M2<-matrix(data2$M,ncol=ncol(data2$M))
  R2<-matrix(data2$R,ncol=1)
  
  re2<-mediation_net_ADMM_NC(Z2,M2,R2,lambda=lambda,omega=omega,phi=phi,Phi1=data2_Phi1,Phi2=data2_Phi2,
                             rho=rho,rho.increase=rho.increase,tol=tol,max.itr=max.itr,thred=thred,Sigma1=Sigma1,Sigma2=Sigma2,
                             trace=trace,Theta0=Theta0,D0=D0,alpha0=alpha0,beta0=beta0)
  #######################################################
  
  #######################################################
  # active sets
  AB.est<-cbind(re1$AB,re2$AB)
  colnames(AB.est)<-c("data1","data2")
  
  S1<-as.numeric(abs(re1$AB)>zero.cutoff)
  S2<-as.numeric(abs(re2$AB)>zero.cutoff)
  ac.set<-cbind(S1,S2)
  colnames(ac.set)<-c("data1","data2")
  rownames(ac.set)<-paste0("M",1:length(S1))
  
  # Cohen's kappa coefficient (estimate of variable selection stability)
  re.kappa<-cohen.kappa(S1,S2)
  #######################################################
  
  re<-list(AB=AB.est,Active=ac.set,kappa=re.kappa$kappa)
  
  return(re)
}

# lambda and omega are vectors
mediation_net_ADMM_NC_VSS_rep<-function(data1,data2,zero.cutoff=1e-3,
                                        lambda=10^c(seq(-3,1,length.out=31)[-1],seq(1,5,length.out=11)[-1]),omega=0,phi=1,
                                        data1_Phi1=NULL,data1_Phi2=NULL,data2_Phi1=NULL,data2_Phi2=NULL,
                                        rho=1,rho.increase=FALSE,tol=1e-10,max.itr=10000,thred=1e-10,Sigma1=NULL,Sigma2=NULL,
                                        trace=FALSE,Theta0=NULL,D0=NULL,alpha0=NULL,beta0=NULL)
{
  if(length(lambda)==1)
  {
    lambda<-rep(lambda,length(omega))
  }
  if(length(omega)==1)
  {
    omega<-rep(omega,length(lambda))
  }
  
  n.lambda<-length(lambda)
  
  vss<-rep(NA,n.lambda)
  for(i in 1:length(lambda))
  {
    re.tmp<-mediation_net_ADMM_NC_VSS(data1,data2,zero.cutoff=zero.cutoff,
                                      lambda=lambda[i],omega=omega[i],phi=phi,
                                      data1_Phi1=data1_Phi1,data1_Phi2=data1_Phi2,data2_Phi1=data2_Phi1,data2_Phi2=data2_Phi2,
                                      rho=rho,rho.increase=rho.increase,tol=tol,max.itr=max.itr,thred=thred,Sigma1=Sigma1,Sigma2=Sigma2,
                                      trace=trace,Theta0=Theta0,D0=D0,alpha0=alpha0,beta0=beta0)
    
    vss[i]<-re.tmp$kappa
    
    # print(i)
  }
  
  re<-cbind(lambda,vss)
  colnames(re)<-c("lambda","vss")
  
  return(re)
}


mediation_net_ADMM_NC_KSC<-function(Z,M,R,zero.cutoff=1e-3,n.rep=5,vss.cut=0.1,
                                    lambda=10^c(seq(-3,1,length.out=31)[-1],seq(1,5,length.out=11)[-1]),
                                    omega=0,phi=1,Phi1=NULL,Phi2=NULL,rho=1,rho.increase=FALSE,
                                    tol=1e-10,max.itr=10000,thred=1e-10,Sigma1=NULL,Sigma2=NULL,
                                    trace=FALSE,Theta0=NULL,D0=NULL,alpha0=NULL,beta0=NULL)
{
  n<-nrow(R)
  k<-ncol(M)
  
  if(length(lambda)==1)
  {
    lambda<-rep(lambda,length(omega))
  }
  if(length(omega)==1)
  {
    omega<-rep(omega,length(lambda))
  }
  
  vss<-matrix(NA,length(lambda),n.rep)
  colnames(vss)<-paste0("sample",1:n.rep)
  
  for(b in 1:n.rep)
  {
    idx1<-sample(1:n,floor(n/2),replace=FALSE)
    idx2<-(1:n)[-idx1]
    
    Z1<-matrix(Z[idx1],ncol=1)
    M1<-matrix(M[idx1,],ncol=k)
    R1<-matrix(R[idx1],ncol=1)
    
    Z2<-matrix(Z[idx2],ncol=1)
    M2<-matrix(M[idx2,],ncol=k)
    R2<-matrix(R[idx2],ncol=1)
    
    data1<-list(Z=Z1,M=M1,R=R1)
    data2<-list(Z=Z2,M=M2,R=R2)
    
    re.tmp<-mediation_net_ADMM_NC_VSS_rep(data1,data2,zero.cutoff=zero.cutoff,
                                          lambda=lambda,omega=omega,phi=phi,
                                          data1_Phi1=Phi1,data1_Phi2=Phi2,data2_Phi1=Phi1,data2_Phi2=Phi2,
                                          rho=rho,rho.increase=FALSE,tol=tol,max.itr=max.itr,thred=thred,Sigma1=Sigma1,Sigma2=Sigma2,
                                          trace=trace,Theta0=Theta0,D0=D0,alpha0=alpha0,beta0=beta0)
    
    vss[,b]<-re.tmp[,2]
    
    # print(b)
  }
  
  vss.est<-apply(vss,1,mean,na.rm=TRUE)
  
  lambda.idx<-which((vss.est/max(vss.est))>=(1-vss.cut))
  
  if(length(unique(lambda))==1)
  {
    lambda.est<-min(omega[lambda.idx])
  }else
  {
    lambda.est<-min(lambda[lambda.idx]) 
  }
  
  re<-list(vss.mat=vss,vss=vss.est,lambda.idx=lambda.idx,lambda.est=lambda.est)
  return(re)
}
#######################################################################################
