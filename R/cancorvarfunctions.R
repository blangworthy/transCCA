kernelfun <-function(i,j,var1,var2,data) {
  sign(data[i,var1] - data[j,var1])* sign(data[i,var2] - data[j,var2])
}

transcorKvar <- function(x,y){
  if(sum(is.na(x)) >0 | sum(is.na(x)) > 0){stop("Data includes missing values")}
  if(nrow(x)!=nrow(y)){stop("x and y are not same length")}
  data <- cbind(x,y)
  kendallmat <- pcaPP::cor.fk(data)
  dimx <- ncol(x)
  dimy <- ncol(y)
  n <- nrow(x)
  kendallmatvec <- c(as.vector(kendallmat[1:dimx,1:dimx]),
                     as.vector(kendallmat[1:dimx,(dimx+1):(dimx+dimy)]),
                     as.vector(kendallmat[(dimx+1):(dimx+dimy),(dimx+1):(dimx+dimy)]))
  
  vecindices <- cbind(rbind(rep(1:dimx,dimx),rep(1:dimx,each=dimx)),
                      rbind(rep(1:dimx,dimy),rep((dimx+1):(dimx+dimy),each = dimx)),
                      rbind(rep((dimx+1):(dimx+dimy),dimy),rep((dimx+1):(dimx+dimy),each =dimy)))
  
  kernelests <- lapply(as.list(1:n), function(u) apply(vecindices,2,function(v) mean(kernelfun(setdiff(1:n,u),u,v[1],v[2],data))))
  
  ####Gamma estimate from Rublik paper
  gamma <- Reduce('+',lapply(kernelests,function(u) u%*%t(u)))/n
  
  ####Estimate of variance
  vhat <- 4*(gamma - kendallmatvec%*%t(kendallmatvec))
  
  ####Jacobian for transformation
  jac <- diag(pi/2*cos(pi*kendallmatvec/2))
  
  ###Don't need to transpose jac since it is diagonal matrix
  vhattrans <- jac%*%vhat%*%jac
  return(vhattrans)
}


####Finds asymptotic variance of A%*%PXX%*%t(A), A%*%PXY%*%t(B) and B%*%PYY%*%t(B)
####Our A, Ahat, B, Bhat are transpose of what they are in paper
vhatuuvv <- function(Ahat,Bhat,vhattrans){
  jacuuvv<-rbind(cbind(Ahat%x%Ahat,matrix(0,nrow = nrow(Ahat)^2,ncol = ncol(Ahat)*ncol(Bhat)),matrix(0,nrow = nrow(Ahat)^2,ncol = ncol(Bhat)^2)),
                 cbind(matrix(0,nrow = nrow(Ahat)*nrow(Bhat),ncol = ncol(Ahat)^2),Bhat%x%Ahat,matrix(0,nrow = nrow(Bhat)*nrow(Ahat),ncol = ncol(Bhat)^2)),
                 cbind(matrix(0,nrow = nrow(Bhat)^2,ncol = ncol(Ahat)^2),matrix(0,nrow = nrow(Bhat)^2,ncol = ncol(Bhat)*ncol(Ahat)),Bhat%x%Bhat))
  
  return(jacuuvv%*%vhattrans%*%t(jacuuvv))
}

#####Row indices for matrix for the i/jth entry of PUU, PUV or PVV in the vectorized version of these matrices
#####mat should be character string PUU, PUV or PVV indicating which matrix
#####it is
matvecindx <- function(i,j,mat,dimx,dimy){
  if (mat=="PUU"){
    return((j-1)*dimx + i)
  }
  else if (mat == "PUV"){
    return((j-1)*dimx + i+dimx^2)
  }
  else if (mat == "PVV"){
    return((j-1)*dimy + i + dimx^2 + dimy*dimx)
  }
  else { return("Indicate which matrix to use")}
}

#####returns variance of the i/jth entry of the matrix specified by mat
varvecentry <- function(varmat,i,j,mat,dimx,dimy){
  varmat[matvecindx(i,j,mat,dimx,dimy),matvecindx(i,j,mat,dimx,dimy)]
}

covvecentry <- function(varmat,i,j,k,l,mat1,mat2,dimx,dimy){
  varmat[matvecindx(i,j,mat1,dimx,dimy),matvecindx(k,l,mat2,dimx,dimy)]
}

cov.var.ii <- function(varmat,i,dimx,dimy){
  cov.mat <- diag(c(varvecentry(varmat,i,i,"PUU",dimx,dimy),varvecentry(varmat,i,i,"PVV",dimx,dimy),varvecentry(varmat,i,i,"PUV",dimx,dimy)))
  cov.mat[1,2] <- covvecentry(varmat,i,i,i,i,"PUU","PVV",dimx,dimy)
  cov.mat[1,3] <- covvecentry(varmat,i,i,i,i,"PUU","PUV",dimx,dimy)
  cov.mat[2,3] <- covvecentry(varmat,i,i,i,i,"PVV","PUV",dimx,dimy)
  cov.mat[lower.tri(cov.mat)] <- t(cov.mat)[lower.tri(t(cov.mat))]
  return(cov.mat)
}

cov.var.ij <- function(varmat,i,j,dimx,dimy){
  cov.mat <- diag(c(varvecentry(varmat,i,j,"PUU",dimx,dimy),varvecentry(varmat,i,j,"PVV",dimx,dimy),varvecentry(varmat,i,j,"PUV",dimx,dimy),varvecentry(varmat,j,i,"PUV",dimx,dimy)))
  cov.mat[1,2] <- covvecentry(varmat,i,j,i,j,"PUU","PVV",dimx,dimy)
  cov.mat[1,3] <- covvecentry(varmat,i,j,i,j,"PUU","PUV",dimx,dimy)
  cov.mat[1,4] <- covvecentry(varmat,i,j,j,i,"PUU","PUV",dimx,dimy)
  cov.mat[2,3] <- covvecentry(varmat,i,j,i,j,"PVV","PUV",dimx,dimy)
  cov.mat[2,4] <- covvecentry(varmat,i,j,j,i,"PVV","PUV",dimx,dimy)
  cov.mat[3,4] <- covvecentry(varmat,i,j,j,i,"PUV","PUV",dimx,dimy)
  cov.mat[lower.tri(cov.mat)] <- t(cov.mat)[lower.tri(t(cov.mat))]
  return(cov.mat)
}




####variance for gij and hij, can get same thing from covariance functions with i=k and j=l
g.ij.var <- function(varmat,i,j,Lambda,dimx,dimy){
  cov.mat <- cov.var.ij(varmat,i,j,dimx,dimy)
  mat <- as.matrix(c(-Lambda[i]*Lambda[j],-Lambda[j]^2,Lambda[j],Lambda[i]))
  var.out <- (t(mat)%*%cov.mat%*%mat)/((Lambda[j]^2-Lambda[i]^2)^2)
  return(var.out)
}

h.ij.var <- function(varmat,i,j,Lambda,dimx,dimy){
  cov.mat <- cov.var.ij(varmat,i,j,dimx,dimy)
  mat <- as.matrix(c(-Lambda[j]^2,-Lambda[i]*Lambda[j],Lambda[i],Lambda[j]))
  var.out <- (t(mat)%*%cov.mat%*%mat)/((Lambda[j]^2-Lambda[i]^2)^2)
  return(var.out)
}

####Covariance between g_ij and g_kl
g.ij.kl.cov <- function(varmat,i,j,k,l,Lambda,dimx,dimy,ndir){
  if(i<=ndir ){
    indx.ij <- c(sapply(c("PUU","PVV","PUV"),function(x) matvecindx(i,j,x,dimx,dimy)),matvecindx(j,i,"PUV",dimx,dimy))
    mat1 <- as.matrix(c(-Lambda[j]^2,-Lambda[j]*Lambda[i],Lambda[j],Lambda[i]))
  }
  if (i >ndir){
    indx.ij <- c(sapply(c("PUU","PUV"),function(x) matvecindx(i,j,x,dimx,dimy)))
    mat1 <- as.matrix(c(-Lambda[j]^2,Lambda[j]))
  }
  if(k<=ndir){
    indx.kl <- c(sapply(c("PUU","PVV","PUV"),function(x) matvecindx(k,l,x,dimx,dimy)),matvecindx(l,k,"PUV",dimx,dimy))
    mat2 <- as.matrix(c(-Lambda[l]^2,-Lambda[k]*Lambda[l],Lambda[l],Lambda[k]))
  }
  if (k >ndir){
    indx.kl <- c(sapply(c("PUU","PUV"),function(x) matvecindx(k,l,x,dimx,dimy)))
    mat2 <- as.matrix(c(-Lambda[l]^2,Lambda[l]))
  }
  cov.mat <- varmat[indx.ij,indx.kl]
  var.out <- (t(mat1)%*%cov.mat%*%mat2)/((Lambda[j]^2-Lambda[i]^2)*(Lambda[l]^2-Lambda[k]^2))
  return(var.out)
}

####Covariance between h_ij and h_kl
h.ij.kl.cov <- function(varmat,i,j,k,l,Lambda,dimx,dimy,ndir){
  if (i <= ndir){
    indx.ij <- c(sapply(c("PUU","PVV","PUV"),function(x) matvecindx(i,j,x,dimx,dimy)),matvecindx(j,i,"PUV",dimx,dimy))
    mat1 <- as.matrix(c(-Lambda[i]*Lambda[j],-Lambda[j]^2,Lambda[i],Lambda[j]))
  }
  if( i > ndir){
    indx.ij <- c(sapply(c("PVV"),function(x) matvecindx(i,j,x,dimx,dimy)),matvecindx(j,i,"PUV",dimx,dimy))
    mat1 <- as.matrix(c(-Lambda[j]^2,Lambda[j]))
  }
  if (k <= ndir){
    indx.kl <- c(sapply(c("PUU","PVV","PUV"),function(x) matvecindx(k,l,x,dimx,dimy)),matvecindx(l,k,"PUV",dimx,dimy))
    mat2 <- as.matrix(c(-Lambda[k]*Lambda[l],-Lambda[l]^2,Lambda[k],Lambda[l]))
  }
  if (k > ndir){
    indx.kl <- c(sapply(c("PVV"),function(x) matvecindx(k,l,x,dimx,dimy)),matvecindx(l,k,"PUV",dimx,dimy))
    mat2 <- as.matrix(c(-Lambda[l]^2,Lambda[l])) 
  }
  cov.mat <- varmat[indx.ij,indx.kl]
  
  
  var.out <- (t(mat1)%*%cov.mat%*%mat2)/((Lambda[j]^2-Lambda[i]^2)*(Lambda[l]^2-Lambda[k]^2))
  return(var.out)
}

####Variance for g_ii and h_ii
g.ii.var <- function(varmat,i,dimx,dimy){
  varvecentry(varmat,i,i,"PUU",dimx,dimy)/4
}

h.ii.var <- function(varmat,i,dimx,dimy){
  varvecentry(varmat,i,i,"PVV",dimx,dimy)/4
}

####covariance between g_ii and g_kl
g.ii.kl.cov <- function(varmat,i,k,l,Lambda,dimx,dimy,ndir){
  indx.i <- matvecindx(i,i,"PUU",dimx,dimy)
  if(k<=ndir){
    indx.kl <- c(sapply(c("PUU","PVV","PUV"),function(x) matvecindx(k,l,x,dimx,dimy)),matvecindx(l,k,"PUV",dimx,dimy))
    mat <- as.matrix(c(-Lambda[l]^2,-Lambda[k]*Lambda[l],Lambda[l],Lambda[k]))
  }
  if (k >ndir){
    indx.kl <- c(sapply(c("PUU","PUV"),function(x) matvecindx(k,l,x,dimx,dimy)))
    mat <- as.matrix(c(-Lambda[l]^2,Lambda[l]))
  }
  cov.mat <- varmat[indx.i,indx.kl]
  var.out <- t(mat)%*%cov.mat/(2*(Lambda[l]^2-Lambda[k]^2))
  return(var.out)
}

####covariance between h_ii and h_kl
h.ii.kl.cov <- function(varmat,i,k,l,Lambda,dimx,dimy,ndir){
  indx.i <- matvecindx(i,i,"PVV",dimx,dimy)
  if (k <= ndir){
    indx.kl <- c(sapply(c("PUU","PVV","PUV"),function(x) matvecindx(k,l,x,dimx,dimy)),matvecindx(l,k,"PUV",dimx,dimy))
    mat <- as.matrix(c(-Lambda[k]*Lambda[l],-Lambda[l]^2,Lambda[k],Lambda[l]))
  }
  if (k > ndir){
    indx.kl <- c(sapply(c("PVV"),function(x) matvecindx(k,l,x,dimx,dimy)),matvecindx(l,k,"PUV",dimx,dimy))
    mat <- as.matrix(c(-Lambda[l]^2,Lambda[l])) 
  }
  cov.mat <- varmat[indx.i,indx.kl]
  var.out <- t(mat)%*%cov.mat/(2*(Lambda[l]^2-Lambda[k]^2))
  return(var.out)
}

