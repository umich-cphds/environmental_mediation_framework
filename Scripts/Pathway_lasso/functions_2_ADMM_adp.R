#################################################################
# Adaptive pathway lasso functions
#################################################################

#################################################################
# log-likelihood function
mediation_net_ll<-function(Z,M,R,A,B,C,Sigma)
{
  n<-nrow(M)
  k<-ncol(M)
  
  Sigma1<-matrix(Sigma[1:k,1:k],k,k)
  Sigma2<-matrix(Sigma[(k+1),(k+1)],1,1)
  
  l1<--n*sum(log(diag(Sigma1)))/2
  l2<--n*log(Sigma2[1,1])/2
  if(k==1)
  {
    l3<--(t(M-Z%*%A)%*%(M-Z%*%A))[1,1]/(2*Sigma1[1,1])
  }else
  {
    l3<--sum(diag(diag(1/diag(Sigma1))%*%t(M-Z%*%A)%*%(M-Z%*%A)))/2
  }
  l4<--(t(R-Z%*%C-M%*%B)%*%(R-Z%*%C-M%*%B))[1,1]/(2*Sigma2[1,1])
  const.MR<-log(2*pi)*(-(n*(k+1))/2)  
  l.MR<-l1+l2+l3+l4+const.MR
  
  Sigma.C2<-(t(B)%*%Sigma1%*%B+Sigma2)[1,1]
  C2<-(C+A%*%B)[1,1]
  l.R<--n*log(2*pi*Sigma.C2)/2-(t(R-Z%*%C2)%*%(R-Z%*%C2))[1,1]/(2*Sigma.C2)
  
  re<-list(ll.MR=l.MR,ll.R=l.R,ll.AS=l1+l3,ll.BS=l2+l4,ll.A=l3,ll.B=l4)
  
  return(re)
}
#################################################################

#################################################################
# solution function

# soft-thresholding function
soft_thred<-function(mu,lambda)
{
  return(max(abs(mu)-lambda,0)*sign(mu))
}

# different omega's
solution<-function(lambda,omega1=0,omega2=0,phi1,phi2,mu1,mu2)
{
  if(lambda==0)
  {
    a<-soft_thred(mu1,omega1)/phi1
    b<-soft_thred(mu2,omega2)/phi2
    
    return(data.frame(a=a,b=b))
  }else
  {
    de<-phi1*phi2-lambda^2
    
    x1<-phi2*(mu1-omega1)-lambda*(mu2-omega2)
    x2<-phi2*(mu1-omega1)+lambda*(mu2+omega2)
    x3<-phi2*(mu1+omega1)+lambda*(mu2-omega2)
    x4<-phi2*(mu1+omega1)-lambda*(mu2+omega2)
    
    y1<-phi1*(mu2-omega2)-lambda*(mu1-omega1)
    y2<-phi1*(mu2+omega2)+lambda*(mu1-omega1)
    y3<-phi1*(mu2-omega2)+lambda*(mu1+omega1)
    y4<-phi1*(mu2+omega2)-lambda*(mu1+omega1)
    
    s1<-data.frame(a=x1/de,b=y1/de)
    s2<-data.frame(a=x2/de,b=y2/de)
    s3<-data.frame(a=x3/de,b=y3/de)
    s4<-data.frame(a=x4/de,b=y4/de)
    s5<-data.frame(a=soft_thred(mu1,omega1)/phi1,b=0)
    s6<-data.frame(a=0,b=soft_thred(mu2,omega2)/phi2)
    s7<-data.frame(a=0,b=0)
    
    if(x1>0&y1>0)
    {
      return(s1)
    }else
      if(x2>0&y2<0)
      {
        return(s2)
      }else
        if(x3<0&y3>0)
        {
          return(s3)
        }else
          if(x4<0&y4<0)
          {
            return(s4)
          }else
            if(abs(mu1)>omega1&((phi1*abs(mu2)-lambda*abs(mu1))<=(-lambda*omega1+phi1*omega2)))
            {
              return(s5)
            }else
              if(abs(mu2)>omega2&((phi2*abs(mu1)-lambda*abs(mu2))<=(-lambda*omega2+phi2*omega1)))
              {
                return(s6)
              }else
              {
                return(s7)
              }
  }
}
#################################################################

#################################################################
# objective function
opt_func<-function(Z,M,R,Theta,D,Sigma1,Sigma2,lambda=1,phi=1,omega=0,Phi1=NULL,Phi2=NULL)
{
  n<-nrow(M)
  k<-ncol(M)
  
  if(is.null(Phi1)==TRUE)
  {
    Phi1<-diag(c(0,rep(1,k)))
  }
  if(is.null(Phi2)==TRUE)
  {
    Phi2<-diag(rep(1,k+1))
  }
  
  X<-cbind(Z,M)
  Omega1<-matrix(0,k+1,k+1)
  Omega1[2:(k+1),2:(k+1)]<-solve(Sigma1)/n  ##Omega1 and w2 are essentially the inversed covariance matrix
  w2<-1/(Sigma2*n)
  
  J<-c(0,rep(1,k))
  
  u<-sum(diag(Omega1%*%t(X-Z%*%Theta)%*%(X-Z%*%Theta)))/2+w2*t(R-X%*%D)%*%(R-X%*%D)/2
  v1<-lambda*(abs(Theta)%*%Phi1%*%Phi2%*%abs(D)+phi*Theta%*%Phi1^2%*%t(Theta)+phi*t(D)%*%Phi2^2%*%D)
  v2<-omega*(abs(Theta)%*%Phi1%*%J+t(abs(D))%*%Phi2%*%J)
  
  f<-as.numeric(u+v1+v2)  ## only u and v, the function needed to be minimized, no constraints included
  
  return(f)
}

#################################################################
# ADMM for (adptive) Pathway Lasso without constraint
mediation_net_ADMM_NC<-function(Z,M,R,lambda=1,omega=0,phi=1,Phi1=NULL,Phi2=NULL,rho=1,rho.increase=FALSE,tol=1e-10,max.itr=10000,thred=1e-10,est_thred=FALSE,Sigma1=NULL,Sigma2=NULL,
                                trace=FALSE,Theta0=NULL,D0=NULL,alpha0=NULL,beta0=NULL)
{
  n<-nrow(M)
  k<-ncol(M)
  
  if(is.null(Phi1)==TRUE)
  {
    Phi1<-diag(c(0,rep(1,k)))
  }
  if(is.null(Phi2)==TRUE)
  {
    Phi2<-diag(rep(1,k+1))
  }
  
  X<-cbind(Z,M)
  e1<-c(1,rep(0,k))
  
  # C'
  C2.hat<-solve(t(Z)%*%Z)%*%t(Z)%*%R
  tau2.hat<-(t(R-Z%*%C2.hat)%*%(R-Z%*%C2.hat)/n)[1,1]
  
  # OLS for A
  A.hat<-solve(t(Z)%*%Z)%*%t(Z)%*%M
  Theta.tilde<-cbind(1,A.hat)
  Sigma1.tilde<-diag(diag(t(M-Z%*%A.hat)%*%(M-Z%*%A.hat)/n))
  
  # Initial values
  if(is.null(Theta0)==TRUE)
  {
    Theta0<-matrix(rep(0,k+1),nrow=1)
  }
  if(is.null(alpha0)==TRUE)
  {
    alpha0<-matrix(rep(0,k+1),nrow=1)
  }
  if(is.null(D0)==TRUE)
  {
    D0<-matrix(rep(0,k+1),ncol=1)
  }
  if(is.null(beta0)==TRUE)
  {
    beta0<-matrix(rep(0,k+1),ncol=1)
  }
  if(is.null(Sigma1)==TRUE)
  {
    Sigma10<-diag(diag(t(M-Z%*%Theta0[-1])%*%(M-Z%*%Theta0[-1])/n))
  }else
  {
    Sigma10<-Sigma1
  }
  if(is.null(Sigma2)==TRUE)
  {
    Sigma20<-t(R-X%*%D0)%*%(R-X%*%D0)
  }else
  {
    Sigma20<-Sigma2
  }
  
  nu1=nu2<-rep(0,k+1)
  nu3<-0
  
  if(rho.increase==FALSE)
  {
    rho0<-0
  }else
  {
    rho0<-rho
  }
  
  # Trace
  if(trace==TRUE)
  {
    Theta.trace=D.trace=alpha.trace=beta.trace<-NULL
    Theta.sg=D.sg=alpha.sg=beta.sg<-NULL
    f<-NULL
    
    Theta.trace<-rbind(Theta.trace,Theta0)
    D.trace<-cbind(D.trace,D0)
    alpha.trace<-rbind(alpha.trace,alpha0)
    beta.trace<-cbind(beta.trace,beta0)
  }
  
  diff<-100
  s<-0
  
  time<-system.time(
    while(diff>=tol&s<=max.itr)
    {
      s<-s+1
      Omega1<-matrix(0,k+1,k+1)
      Omega1[2:(k+1),2:(k+1)]<-solve(Sigma10)/n
      w2<-1/(Sigma20*n)
      
      # Update Theta
      de.Theta<-solve((t(Z)%*%Z)[1,1]*Omega1+2*rho*diag(rep(1,k+1))+2*rho*e1%*%t(e1))
      nu.Theta<-t(Z)%*%X%*%Omega1-t(nu1)+2*rho*alpha0+(2*rho-nu3)*t(e1)
      Theta.new<-nu.Theta%*%de.Theta
      if(est_thred==TRUE)
      {
        Theta.idx<-which(abs(Theta.new)<thred)
        Theta.new[Theta.idx]<-0
      }
      
      # Update D
      de.D<-solve(w2[1,1]*t(X)%*%X+2*rho*diag(rep(1,k+1)))
      nu.D<-w2[1,1]*t(X)%*%R-nu2+2*rho*beta0
      D.new<-de.D%*%nu.D
      if(est_thred==TRUE)
      {
        D.idx<-which(abs(D.new)<thred)
        D.new[D.idx]<-0
      }
      
      # Update alpha and beta
      alpha.new<-matrix(rep(NA,k+1),nrow=1)
      beta.new<-matrix(rep(NA,k+1),ncol=1)
      J<-c(0,rep(1,k))
      if(lambda==0)
      {
        J[1]<-1
      }
      for(j in 1:(k+1))
      {
        pp1<-2*lambda*phi*Phi1[j,j]^2+2*rho
        pp2<-2*lambda*phi*Phi2[j,j]^2+2*rho
        mu1<-nu1[j]+2*rho*Theta.new[1,j]
        mu2<-nu2[j]+2*rho*D.new[j,1]
        re<-solution(lambda=lambda*Phi1[j,j]*Phi2[j,j],omega1=omega*Phi1[j,j]*J[j],omega2=omega*Phi2[j,j]*J[j],
                     phi1=pp1,phi2=pp2,mu1=mu1,mu2=mu2)
        alpha.new[1,j]<-re$a
        beta.new[j,1]<-re$b
      }
      
      nu1<-nu1+2*rho*t(Theta.new-alpha.new)
      nu2<-nu2+2*rho*(D.new-beta.new)
      nu3<-nu3+2*rho*(Theta.new%*%e1-1)[1,1]
      
      d.Theta<-max(abs(Theta.new-Theta0))
      d.D<-max(abs(D.new-D0))
      
      diff<-max(d.Theta,d.D)
      
      Theta0<-Theta.new
      D0<-D.new
      alpha0<-alpha.new
      beta0<-beta.new
      
      if(is.null(Sigma1)==TRUE)
      {
        Sigma10<-diag(diag(t(M-Z%*%Theta0[-1])%*%(M-Z%*%Theta0[-1])/n))
      }else
      {
        Sigma10<-Sigma1
      }
      if(is.null(Sigma2)==TRUE)
      {
        Sigma20<-t(R-X%*%D0)%*%(R-X%*%D0)
      }else
      {
        Sigma20<-Sigma2
      }
      
      if(trace==TRUE)
      {
        # Trace
        Theta.trace<-rbind(Theta.trace,Theta0)
        D.trace<-cbind(D.trace,D0)
        alpha.trace<-rbind(alpha.trace,alpha0)
        beta.trace<-cbind(beta.trace,beta0)
        
        # subgradient
        Theta.sg<-rbind(Theta.sg,Theta0%*%solve(de.Theta)-nu.Theta)
        D.sg<-cbind(D.sg,solve(de.D)%*%D0-nu.D)
        alpha.sg<-rbind(alpha.sg,lambda*(t(abs(beta.new))*sign(alpha.new))%*%(Phi1*Phi2)+2*lambda*phi*alpha.new%*%(Phi1^2)+
                          omega*(sign(alpha.new)*t(J))%*%(Phi1))
        beta.sg<-cbind(beta.sg,lambda*(Phi1*Phi2)%*%(t(abs(alpha.new))*sign(beta.new))+2*lambda*phi*(Phi2^2)%*%beta.new+
                         omega*Phi2%*%(J*sign(beta.new)))
        
        f<-c(f,opt_func(Z,M,R,Theta0,D0,Sigma10,Sigma20,lambda=lambda,phi=phi,omega=omega,Phi1=Phi1,Phi2=Phi2))
      }
    }
  )
  if(s>max.itr)
  {
    warning("Method does not converge!")
  }
  
  constraint1=constraint2<-matrix(NA,1,3)
  colnames(constraint1)=colnames(constraint2)<-c("Theta=alpha","D=beta","Theta[1]")
  constraint1[1,1]<-(max(abs(Theta0-alpha0))<thred)
  constraint1[1,2]<-(max(abs(D0-beta0))<thred)
  constraint1[1,3]<-(abs(Theta0[1,1]-1)<thred)
  constraint2[1,1]<-max(abs(Theta0-alpha0))
  constraint2[1,2]<-max(abs(D0-beta0))
  constraint2[1,3]<-Theta0[1,1]
  constraint<-cbind(data.frame(t(constraint1)),data.frame(t(constraint2)))
  colnames(constraint)<-c("Satisfied","value")
  
  if(constraint[1,1]==TRUE)
  {
    A.hat<-matrix(alpha0[-1],nrow=1)
    colnames(A.hat)<-paste0("M",1:k)
    rownames(A.hat)<-"Z"
  }else
  {
    A.hat<-matrix(Theta0[-1],nrow=1)
    colnames(A.hat)<-paste0("M",1:k)
    rownames(A.hat)<-"Z"
  }
  if(constraint[2,1]==TRUE)
  {
    D.hat<-beta0
    C.hat<-matrix(beta0[1],1,1)
    colnames(C.hat)<-"R"
    rownames(C.hat)<-"Z"
    B.hat<-matrix(beta0[-1],ncol=1)
    colnames(B.hat)<-"R"
    rownames(B.hat)<-paste0("M",1:k) 
  }else
  {
    D.hat<-D0
    C.hat<-matrix(D0[1],1,1)
    colnames(C.hat)<-"R"
    rownames(C.hat)<-"Z"
    B.hat<-matrix(D0[-1],ncol=1)
    colnames(B.hat)<-"R"
    rownames(B.hat)<-paste0("M",1:k)
  }
  Sigma1.hat<-Sigma10
  colnames(Sigma1.hat)=rownames(Sigma1.hat)<-paste0("M",1:k)
  Sigma2.hat<-Sigma20
  colnames(Sigma2.hat)=rownames(Sigma2.hat)<-"R"
  Sigma.hat<-cbind(rbind(Sigma1.hat,matrix(0,1,k)),rbind(matrix(0,k,1),Sigma2.hat))
  colnames(Sigma.hat)=rownames(Sigma.hat)<-c(colnames(Sigma1.hat),colnames(Sigma2.hat))
  
  ll<-mediation_net_ll(Z,M,R,A.hat,B.hat,C.hat,Sigma.hat)
  
  d<-sum(abs(c(t(A.hat)*B.hat,C.hat))>thred)
  net.BIC<--2*ll$ll.MR+d*log(n)
  
  if(trace==TRUE)
  {
    re<-list(lambda=lambda,omega=omega,phi=phi,Phi1=Phi1,Phi2=Phi2,rho=rho,A=A.hat,C=C.hat,B=B.hat,AB=t(A.hat)*B.hat,Sigma=Sigma.hat,
             C2=C2.hat[1,1],tau2=tau2.hat,logLik=data.frame(ZMR=ll$ll.MR,ZR=ll$ll.R,lA=ll$ll.A,lB=ll$ll.B),
             BIC=net.BIC,converge=(s<=max.itr),alpha=alpha0,beta=beta0,Theta.1=Theta0[1,1],constraint=constraint,time=time,
             Theta.trace=Theta.trace,D.trace=D.trace,alpha.trace=alpha.trace,beta.trace=beta.trace,
             Theta.sg=Theta.sg,D.sg=D.sg,alpha.sg=alpha.sg,beta.sg=beta.sg,opt_func.sg=f)
  }else
  {
    re<-list(lambda=lambda,omega=omega,phi=phi,Phi1=Phi1,Phi2=Phi2,rho=rho,A=A.hat,C=C.hat,B=B.hat,AB=t(A.hat)*B.hat,Sigma=Sigma.hat,
             C2=C2.hat[1,1],tau2=tau2.hat,logLik=data.frame(ZMR=ll$ll.MR,ZR=ll$ll.R,lA=ll$ll.A,lB=ll$ll.B),
             BIC=net.BIC,converge=(s<=max.itr),alpha=alpha0,beta=beta0,Theta.1=Theta0[1,1],constraint=constraint,time=time)
  }
  
  return(re)
}






