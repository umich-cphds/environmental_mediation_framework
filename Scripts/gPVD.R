##########################################################################################################################################
##########################################################################################################################################
### Generalized Populaton Value Decomposition
### Last Update: Oct 18, 2017
### Author: Oliver Y. Chén
### Department of Engineering Science, University of Oxford
### Description: The codes input a large data matrices from n subjects: M_1, M_2, ... M_n, each of which is of n_i by p, where p is high-dimensional
### 
### Depends: Unix, Linux
### License: The MIT License (MIT)
###
### References: [1] Crainiceanu et al. (2011) Population value decomposition, a framework for the analysis of image populations. JASA 106, no. 495 (2011): 775-790.
###             [2] Chén et al. (2015) High-dimensional multivariate mediation with application to neuroimaging data. Biostatistics.
##########################################################################################################################################
##########################################################################################################################################

#rm(list = ls())
#setwd("")

#####################################################################
#####################################################################
#####################################################################
### gPVD (generalized Population Value Decomposition)
#####################################################################
#####################################################################
#####################################################################

library(R.matlab)
library(MASS)
library(data.table)

#####################################################################
### Input
#####################################################################
grp <- 1  ### number of groups
r <- 25  ### remaining dimension after data reduction; namely, reducing from n by p to n by r, where r << p.

#M <- list(log(d2[,40:100]))
for (i in 1: grp) {
  M[[i]] =  m1                                 ### Read data M[[i]], each is a n_i by p matrix containing data 
                                               ### from subject i, where n_i refers to the number repeated measurement
                                               ### and p is the number of features (or data measured at p time points).  
                                               ### If data consists groups, then M[[i]] is data from the group i, and
                                               ### n_i is the number of subjects in that group.
}

##################
# Output
#################
### There will be two elements in the output, namely [D_r] and [M_tilde_r_final].
### Consider the generalized PVD [see Chén et al. (2015) High-dimensional multivariate mediation with application to neuroimaging data],
### where M = PVD \approximate M_tilde_r_final %*% D_r.
### M_tilde_r_final is the dimension-reduced matrix of n by r, where r << p, and D_r containing the rest information (variance) of M.



##################
# Begin
#################

##################
# SVD of each M_i
#################
U=list()
Sigma=list()
V=list()

for (i in 1: grp ) {
  
  U[[i]] = svd(M[[i]])$u
  Sigma[[i]] = svd(M[[i]])$d
  V[[i]] = svd(M[[i]])$v
}

for (i in 1: grp ) {
  
  #write.csv(U[[i]], file = paste0("U", i, ".csv"),row.names = FALSE)
  #write.csv(Sigma[[i]], file = paste0("Sigma", i, ".csv"),row.names = FALSE)
  #write.csv(V[[i]], file = paste0("V", i, ".csv"),row.names = FALSE)
}

#####################################################################
### Forming V  
#####################################################################
require(data.table)
#V=list()

for (i in 1:grp) {
   #V[[i]] =  as.matrix ( fread(paste0("V", i, ".csv")) )  #[, 1: n_pca]
   V[[i]] = V[[i]]

}

#V_final=cbind(V[[1]], V[[2]], V[[3]])
V_final=cbind(V[[1]])

############################################################################
############################################################################
m <- floor(nrow(V_final)/1000)   ### break the data row-wise into m chunks each is of 1000 by p

for (i in 1: m) {
  
  #write.csv(V_final[ (1000*(i-1)+1): (1000*i), ],  file=  paste0("V_tilde_", i, ".csv"), row.names = FALSE ) 
  
}

#write.csv(V_final[m*1000: nrow(V_final), ],  file=  paste0("V_tilde_", (m+1), ".csv"), row.names = FALSE ) 

################################################################################
################################################################################

VV=list()     ### a list of V_k ^T V_k

for (i in 1: m) {
  
  VV[[i]] = tcrossprod(t(V_final[ (1000*(i-1)+1): (1000*i), ]))

}


VV_final=Reduce('+', VV)     

#####################################################################
### SVD on VV_final to obtain B and C
#####################################################################

svd_VV=svd(VV_final)
C_tilde=svd_VV$u

B=sqrt(svd_VV$d)

B_tilde = svd_VV$d

B_tilde_inv=diag(1/B)

CB_inv=C_tilde %*%B_tilde_inv   # C_tilde * B_tilde_inv

#######################################################################
### Obtain A_k tilde, i.e. A_k tilde = V_k tilde C_tilde B_tilde ^ {-1}
#######################################################################
A_tilde = list()

for (i in 1: m) {
  
  A_tilde[[i]] = as.matrix( V_final[ (1000*(i-1)+1): (1000*i), ] ) %*% CB_inv
  
}


###  
library(plyr)
A_tilde_final=rbind.fill.matrix(A_tilde)  ### binding all A_tilde_k into A_tilde_final

#####################################################################
### Form D: take the first (20) columns of A_tilde_final
#####################################################################
D_r = t ( A_tilde_final[,1:r] )
#####################################################################
### Get U_k
#####################################################################


#####################################################################
### Get V_rk   Here we take R_k as all columns of V_k
#####################################################################
VD_r=list()

for (i in 1:grp) {
  VD_r[[i]] =  tcrossprod( (as.matrix(t(V[[i]]))[1:r,]), D_r )   ## here the V_i is from M_k = U_k Sigma_k V_k
}

#####
V_rk_r =list()

for (i in 1:grp) {
  V_rk_r[[i]] =  diag( Sigma[[i]] ) [, 1:r] %*% VD_r[[i]]     ### choose the first n_pca number of sigma
  ## here the Sigma_i is from M_k = U_k Sigma_k V_k
}

#####################################################################
### Form M_tilde
#####################################################################
M_tilde_r=list()

for (i in 1:grp) {
  M_tilde_r[[i]] = as.matrix( U[[i]] ) %*% V_rk_r[[i]]
}

M_tilde_r_final <- do.call(rbind,M_tilde_r)    ### binding all groups.  M_tilde_r_final is the reduced data of n by r, where r << p.


#M_tilde_r_final <- rbind(M_tilde_r[[1]], M_tilde_r[[2]], M_tilde_r[[3]])


#####################################################################
### output
#####################################################################

#write.csv(D_r, file = "D_r.csv",row.names = FALSE)
#write.csv(M_tilde_r_final, file = "M_tilde_r.csv",row.names = FALSE)

########################################################################################################################
########################################################################################################################
########################################################################################################################
#### END 
########################################################################################################################
########################################################################################################################
########################################################################################################################
