########################################################################
## Group SCCA+ algorithm
## handle group lasso with overlap
## 
## Alissa Castleberry
## Written 05/08/18
## Modified from Lei Du, leidu@iu.edu
## Date created: DEC-19-2014
## Date updated: Jan-17-2015
## @Indiana University School of Medicine.
#######################################################################
#Note: all code within ### blocks is mine, all code within --- blocks is original


##--------------------------------------------------------------------
##Input:
#       - X, geno matrix
#       - Y, pheno matrix
#       - paras, parameters: gamma1, gamma2, beta1, beta2 (unknown, tuned from other functions)
# Output:
#       - u, weight of X
#       - v, weight of Y
#       - corrs, all corr of every iteration
#       - U, all u of every iteration
#       - V, all v of every iteration.
#---------------------------------------

######################################################################
##Input
##Geno matrix
genoMat <- read.csv("CRC_Metabolites.pone.0152126.s002.csv", header = TRUE)
##Pheno matrix
phenoMat <- read.csv("CRC_Metabolites.pone.0152126.s002.csv", header = TRUE)

A <- matrix(c(12, 4, 2, 14, 5, 6), nrow = 2, ncol = 3)
B <- matrix(c(3, 2, 1, 7, 1, 2), nrow = 2, ncol = 3)
######################################################################

##-------------------------
##set parameters
##not necessary in R?
#alpha1 = paras.alpha1;
#alpha2 = paras.alpha2;
#lambda1 = paras.lambda1;
#lambda2 = paras.lambda2;
#beta1 = paras.beta1;
#beta2 = paras.beta2;
#n_snp = size(X,2);
#n_pheno = size(Y,2);
##-------------------------

######################################################################
##The function
######################################################################

agn_scca <- function(X, Y, paras){
######################################################################
##Calculate covariance within geno and pheno 
##matrix multiplciation of transposed X and X
XX <- H(X) %*% X 
YY <- H(Y) %*% Y
##direct calculation of covariance of X 
#XX <- cov(X)
#YY <- cov(Y)
######################################################################

######################################################################
##Calculate covariance between geno and pheno
XY <- H(X) %*% Y
YX <- H(Y) %*% X

#need to load pracma to use bsxfun 
#library('pracma')
#XY <- H(bsxfun("-", X, colMeans(X))) * (bsxfun("-", Y, colMeans(Y))/(ncol(X)-1))

#-----------------------------------------------------------------------
#original code
#XY = bsxfun(@minus,X,mean(X))'*bsxfun(@minus,Y,mean(Y))/(size(X,1)-1);
#YX = bsxfun(@minus,Y,mean(Y))'*bsxfun(@minus,X,mean(X))/(size(Y,1)-1);
#-----------------------------------------------------------------------
######################################################################


######################################################################
##Identify matrix 
d10 <- matrix(1, nrow = nrow(X), ncol = 1)
d20 <- matrix(1, nrow = nrow(Y), ncol = 1)
d11 <- matrix(1, nrow = nrow(X), ncol = 1)
d12 <- matrix(1, nrow = nrow(X), ncol = 1)
d21 <- matrix(1, nrow = nrow(Y), ncol = 1)
d22 <- matrix(1, nrow = nrow(Y), ncol = 1)

##pre-calculate for loop 
##set group information 
##cannot find documentation on this function??? 

#%-------------------------
#% pre calculate for loop
#% % set group information
#H1 = getGroupInfo(X,'lp');
#H2 = getGroupInfo(Y,'lp');
#%-------------------------
######################################################################


######################################################################
##Initialization
#u0 <- matrix(1, nrow = nrow(X), ncol = 1)
#v0 <- matrix(1, nrow = nrow(Y), ncol = 1)
#scale1 <- sqrt(H(u0)*XX*u0)
#u <- u0 / scale1
#scale2 <- sqrt(H(v0)*YY*v0)
#v <- v0 / scale2

u <- matrix(1, nrow = nrow(X), ncol = 1) / nrow(X)
v <- matrix(1, nrow = nrow(Y), ncol = 1) / nrow(Y)

#u <- u0 / nrow(X)
#v <- v0 / nrow(Y)
#u0 <- runif(nrow(X))
#v0 <- runif(nrow(Y))
#u <- u0 / Norm(u0)
#v <- v0 / Norm(v0)
######################################################################

max_iter = 50 #pre-set max # of iterations, default is 50
err = 1e-5 # 0.01 ~ 0.05
diff_u = err*10
diff_v = err*10
diff_obj = err*10
i = 0

#while loop to fix v, solve u
while (i < max_iter && diff_u > err && diff_v > err){
  i = i+1 #figure out how to initialize this? 
  iter1 = 0;
  while ((iter1 < max_iter) && (diff_u > err)){
    iter1 = iter1+1
    D11 = diag(d11)
    M1 = alpha1*XX+lambda1*diag((H1*abs(u)))*D11+beta1*D11
    u_new = mldivide(M1,(XY*v))
    #scale1 = sqrt(H(u_new)*XX*u_new)
    #u_new = u_new/scale1
    if (sum(is.nan(u_new))){
      ##eps not defined???
      u = u+eps
      v = v+eps
    }
  }
    diff_u = max(abs(u_new - u))
    u = u_new
    [d11, u_group1] <- updateD(u) ##need to code updateD function 
}

#while loop to fix u, solve v
while (i < max_iter && diff_u > err && diff_v > err){
  i = i+1 #figure out how to initialize this? 
  iter2 = 0
  while ((iter2 < max_iter) && (diff_v > err)){
    iter2 = iter2+1
    D21 = diag(d21)
    M2 = alpha2*YY+lambda2*diag((H2*abs(v)))*D21+beta2*D21
    v_new = mldivide(M2,(YX*u))
    # scale2 = sqrt(v_new'*YY*v_new)
    # v_new = v_new./scale2
    if (sum(is.nan(v_new))){
      u = u+eps
      v = v+eps
    }
  }
  diff_v = max(abs(v_new - v))
  v = v_new
  [d21, v_group1] = updateD(v)
}
}
######################################################################



######################################################################
##updateD function needed for while loop 
##updates the diagnoal matrix
##also written by Lei Du

updateD <- function(beta, struct_in,lnorm){
# group = 0;
# if nargin == 1
#     d = 1 ./ sqrt(beta.^2+eps)
#     group = sum(abs(beta))
# else
  #     [nrow,ncol] = size(group_in)
#     for g_i = 1:nrow
#         idx = group_in(g_i,:)~=0
#         wc1 = beta(idx, :)
#         group = sqrt(sum(wc1.*wc1))+group # for calculate objective function
#         d_gi = sqrt(sum(wc1.*wc1)+eps)
#         beta_i(idx) = d_gi
#     end
#     d = 1 ./ beta_i
# end
# group_out = group

group = 0;
if (nargin == 1){
  d <- 1 / sqrt(beta^2+eps)
  group <- sum(abs(beta))
} elseif (strcmpi(lnorm,'group')){
  for (g_i in 1:nrow(struct_in)){
  idx <- struct_in(g_i,:)~=0
  wc1 <- beta(idx)
  group <- sqrt(sum(wc1.*wc1))+group # for calculate objective function
  d_gi <- sqrt(sum(wc1.*wc1)+eps)
  beta_i(idx) <- d_gi
  }
  d <- 1 / beta_i 
}elseif (strcmpi(lnorm,'graph')){
  coef <- zeros(nrow(struct_in),ncol(struct_in))
  for (g_i in 1:nrow){
    idx0 <- struct_in(g_i,)==0
    wc1 <- beta
    wc1(idx0) <- 0
    group <- sqrt(sum(wc1.*wc1))+group # for calculate objective function
    d_gi <- sqrt(sum(wc1.*wc1)+eps)
    coef(g_i,idx0) <- d_gi
    #beta_i(idx) <- d_gi
  }
beta_i <- sum(coef,1)
d <- 1 / beta_i
}
struct_out <- group
return(d)
}







######################################################################





