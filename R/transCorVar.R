#'Transelliptical Canonical correlation asymptotic variance
#'
#'Estimates the asymptotic variance for the ith transelliptical canonical correlation and direction estimate
#'@param x The first data matrix to be included in the calculation of the transformed Kendall scatter matrix
#'@param y The second data matrix to be included in the calculation of the transformed Kendall scatter matrix
#'@param i The index for the canonical correlation variance to estimate
#'@param eigenmin The minimum eigenvalue when transforming matrix to be positive definite
#'@return The estimate of the asymptotic variance for the ith transelliptical canonical correlation and direction
#'@importFrom pcaPP cor.fk
#'@importFrom expm sqrtm
#'
#'@export
transCorVar <- function(x,y,i,eigenmin = 0.001){
  if(sum(is.na(x)) >0 | sum(is.na(x)) > 0){stop("Data includes missing values")}
  if(nrow(x)!=nrow(y)){stop("x and y are not same length")}
  data <- cbind(x,y)
  dimx <- ncol(x)
  dimy <- ncol(y)
  ndir <- min(dimx,dimy)
  if (i >min(dimx,dimy)){stop("subscript for canonical correlation out of bounds")}
  n <- nrow(x)
  Lambda <- c(transCCA(x,y,eigenmin)$cancor,rep(0,max(dimx,dimy) - min(dimx,dimy)))
  Ahat <- transCCAalldir(x,y,eigenmin)$xcoef
  Bhat <- transCCAalldir(x,y,eigenmin)$ycoef
  
  
  vhatkend <- transcorKvar(x,y)
  
  vhatuv <- vhatuuvv(Ahat,Bhat,vhatkend)
  
  covmat <- cov.var.ii(vhatuv,i,dimx,dimy)
  mat <- as.matrix(c(-Lambda[i]^2,-Lambda[i]^2,2*Lambda[i]))
  varcor <- (t(mat)%*%covmat%*%mat)/(4*Lambda[i]^2)
  
  covmatsa <- vector(mode ="list",length = dimx*dimx)
  counter <- 1
  for(j in 1:dimx){
    for (k in 1:dimx){
      if (j==i & k == i){
        covmatsa[[counter]]<- Ahat[i,]%*%t(Ahat[i,])*as.numeric(g.ii.var(vhatuv,i,dimx,dimy))
      }
      else if(j==i){
        covmatsa[[counter]] <- Ahat[i,]%*%t(Ahat[k,])*as.numeric(g.ii.kl.cov(vhatuv,i,k,i,Lambda,dimx,dimy,ndir))
      }
      else if(k==i){
        covmatsa[[counter]] <- Ahat[j,]%*%t(Ahat[i,])*as.numeric(g.ii.kl.cov(vhatuv,i,j,i,Lambda,dimx,dimy,ndir))
      }
      else{
        covmatsa[[counter]] <- Ahat[j,]%*%t(Ahat[k,])*as.numeric(g.ij.kl.cov(vhatuv,j,i,k,i,Lambda,dimx,dimy,ndir))
      }
      counter <- counter + 1
      
    }
  }
  
  xdirvar <- Reduce("+",covmatsa)
  
  covmatsb <- vector(mode ="list",length = dimy*dimy)
  counter <- 1
  for(j in 1:dimy){
    for (k in 1:dimy){
      if (j==i & k == i){
        covmatsb[[counter]]<- Bhat[i,]%*%t(Bhat[i,])*as.numeric(h.ii.var(vhatuv,i,dimx,dimy))
      }
      else if(j==i){
        covmatsb[[counter]] <- Bhat[i,]%*%t(Bhat[k,])*as.numeric(h.ii.kl.cov(vhatuv,i,k,i,Lambda,dimx,dimy,ndir))
      }
      else if(k==i){
        covmatsb[[counter]] <- Bhat[j,]%*%t(Bhat[i,])*as.numeric(h.ii.kl.cov(vhatuv,i,j,i,Lambda,dimx,dimy,ndir))
      }
      else{
        covmatsb[[counter]] <- Bhat[j,]%*%t(Bhat[k,])*as.numeric(h.ij.kl.cov(vhatuv,j,i,k,i,Lambda,dimx,dimy,ndir))
      }
      counter <- counter + 1
    }
  }
  
  
  ydirvar <- Reduce("+",covmatsb)
  outputlist <- list("XCoefVar" = xdirvar,
                     "YCoefVar" = ydirvar,
                     "corVar" = varcor)
  return(outputlist)
}
