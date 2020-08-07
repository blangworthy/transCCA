#'Output transelliptical CCA directions and correlations
#'
#'This takes two data sets and returns the transelliptical canonical directions and correlations
#'@param x The first data matrix to be included in the calculation of the transformed Kendall scatter matrix
#'@param y The second data matrix to be included in the calculation of the transformed Kendall scatter matrix
#'@param eigenmin The minimum eigenvalue when transforming matrix to be positive definite
#'@return The estimate of transelliptical canonical correlations and directions
#'@importFrom pcaPP cor.fk
#'@importFrom expm sqrtm
#'
#'@export
transCCA <- function(x,y,eigenmin=0.001){
  cormat<- transcorK(x,y,eigenmin)
  dimx <- ncol(x)
  dimy <- ncol(y)
  ndir <- min(dimx,dimy)
  
  sigma.xx <- cormat[(1:dimx),(1:dimx)]
  sigma.yy <- cormat[(dimx+1):(dimx+dimy),(dimx+1):(dimx+dimy)]
  sigma.xy <- cormat[1:dimx,(dimx+1):(dimx+dimy)]
  sigma.yx <- cormat[(dimx+1):(dimx+dimy),1:dimx]
  can.eig.vals.y <- eigen(expm::sqrtm(solve(sigma.yy))%*%sigma.yx%*%solve(sigma.xx)%*%sigma.xy%*%expm::sqrtm(solve(sigma.yy)))
  can.eig.vals.x <- eigen(expm::sqrtm(solve(sigma.xx))%*%sigma.xy%*%solve(sigma.yy)%*%sigma.yx%*%expm::sqrtm(solve(sigma.xx)))
  xcoef <- t(can.eig.vals.x$vectors)%*% expm::sqrtm(solve(sigma.xx))
  ycoef <- t(can.eig.vals.y$vectors)%*% expm::sqrtm(solve(sigma.yy))
  for (i in 1:ndir){
    if (xcoef[i,] %*% sigma.xy %*% ycoef[i,] < 0){
      xcoef[i,] <- xcoef[i,]*-1
    }
  }
  cancor <- sqrt(can.eig.vals.x$values[1:ndir])
  outputlist <- list("xcoef"= xcoef[1:ndir,],
                     "ycoef" = ycoef[1:ndir,],
                     "cancor" = cancor)
  class(outputlist) <- 'transcca'
  return(outputlist)
}


transCCAalldir <- function(x,y,eigenmin=0.001){
  ##This outputs all directions to be used for variance estimation
  cormat<- transcorK(x,y,eigenmin)
  dimx <- ncol(x)
  dimy <- ncol(y)
  ndir <- min(dimx,dimy)
  
  
  sigma.xx <- cormat[(1:dimx),(1:dimx)]
  sigma.yy <- cormat[(dimx+1):(dimx+dimy),(dimx+1):(dimx+dimy)]
  sigma.xy <- cormat[1:dimx,(dimx+1):(dimx+dimy)]
  sigma.yx <- cormat[(dimx+1):(dimx+dimy),1:dimx]
  can.eig.vals.y <- eigen(expm::sqrtm(solve(sigma.yy))%*%sigma.yx%*%solve(sigma.xx)%*%sigma.xy%*%expm::sqrtm(solve(sigma.yy)))
  can.eig.vals.x <- eigen(expm::sqrtm(solve(sigma.xx))%*%sigma.xy%*%solve(sigma.yy)%*%sigma.yx%*%expm::sqrtm(solve(sigma.xx)))
  xcoef <- t(can.eig.vals.x$vectors)%*% expm::sqrtm(solve(sigma.xx))
  ycoef <- t(can.eig.vals.y$vectors)%*% expm::sqrtm(solve(sigma.yy))
  
  for (i in 1:ndir){
    if (xcoef[i,] %*% sigma.xy %*% ycoef[i,] < 0){
      xcoef[i,] <- xcoef[i,]*-1
    }
  }
  outputlist <- list("xcoef"= xcoef,
                     "ycoef" = ycoef)
  return(outputlist)
}