#'Create transformed Kendall Scatter Matrix
#'
#'Creates transformed Kendall's scatter matrix between
#'two sets of variables, x and y, to be used in scale invariant CCA
#'@param x The first data matrix to be included in the calculation of the transformed Kendall scatter matrix
#'@param y The second data matrix to be included in the calculation of the transformed Kendall scatter matrix
#'@param eigenmin The minimum eigenvalue when transforming matrix to be positive definite
#'@return The estimate of the scatter matrix for data using the transformation of Kendall's tau
#'@importFrom pcaPP cor.fk
#'
#'@export
transcorK <- function(x,y,eigenmin=0.001){
  if(sum(is.na(x)) >0 | sum(is.na(y)) > 0){stop("Data includes missing values")}
  if(nrow(x)!=nrow(y)){stop("x and y are not same length")}
  data <- cbind(x,y)
  kendallmat <- pcaPP::cor.fk(data)
  cormat <-sin(pi*kendallmat/2) 
  cormateigen <- eigen(cormat)
  
  if(min(cormateigen$values)>0){
    return(cormat)
  }else{
    return(cormateigen$vectors%*%diag(pmax(eigenmin,cormateigen$values))%*%t(cormateigen$vectors))
  }
  
}

