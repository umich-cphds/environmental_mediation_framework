#=================================
# Generate simulation data: causally dependent
sim.data_dep<-function(n,Z,a,b,c,Delta,Xi1,Sigma2)
{
  k<-nrow(b)
  
  M.coef<-rbind(a,Delta)
  D<-rbind(c,b)
  
  eps1<-matrix(NA,n,k)
  M<-matrix(NA,n,k)
  for(j in 1:k)
  {
    eps1[,j]<-rnorm(n,sd=sqrt(Xi1[j,j]))
    if(j==1)
    {
      M[,j]<-Z*A[1,1]+eps1[,j]
    }else
    {
      M[,j]<-cbind(Z,M[,1:(j-1)])%*%M.coef[1:j,j]+eps1[,j]
    }
  }
  eps2<-rnorm(n,sd=sqrt(Sigma2))
  R<-cbind(Z,M)%*%D+eps2
  
  return(list(M=M,R=R))
}
#=================================

#Generate data based on full model:
#a = (a1, a2, ..., ak)
#Delta: k by k
#M[,j] <- cbind(Z, M[,1:(j-1)]) %*% (a, Delta) + esp1
#R = cbind(Z, M) %*% (c, b) + eps2


