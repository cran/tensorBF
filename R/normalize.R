sd.na <- function(v) { sd(v,na.rm=TRUE) }
mean.na <- function(v) { mean(v,na.rm=TRUE) }
mean.fun <- function(x) { mean.na(x) }
get.slab.var <- function(Y,d)
{
  v <- apply(Y,d,function(x) {var(as.vector(x),na.rm=TRUE)})
  return(v)
}

#' Preprocessing: Slab Scaling
#'
#' \code{normSlabScaling} scales the slabs of the \eqn{o^{th}} mode of the tensor to unit variance.
#'
#' @param Y the tensor data. See function \code{\link{tensorBF}} for
#'   details.
#' @param o the \eqn{o^{th}} (default: 2) mode of the tensor in which the slabs are to be scaled to unit variance.
#' @return a list containing the following elements:
#' \item{data}{The data after performing the required scaling operation.}
#' \item{pre}{The scale's used for preprocessing.}
#' @references
#' Kolda, Tamara G., and Brett W. Bader. "Tensor decompositions and applications." SIAM review 51.3 (2009): 455-500.
#' @export
#' @examples
#' #Data generation
#' K <- 3
#' X <- matrix(rnorm(20*K),20,K)
#' W <- matrix(rnorm(30*K),30,K)
#' U <- matrix(rnorm(3*K),3,K)
#' Y = 0
#' for(k in 1:K) Y <- Y + outer(outer(X[,k],W[,k]),U[,k])
#'  Y <- Y + array(rnorm(20*30*3),dim=c(20,30,3))
#'
#' #scale the slabs in second mode of tensor Y
#' res <- normSlabScaling(Y=Y,o=2)
#' dim(res$data) #the scaled data
normSlabScaling <- function(Y,o=2)
{
  if(o<1 || o>length(dim(Y)) || !is.numeric(o))
    stop(paste0("In-appropriate value of slab scaling mode specified. The value should be between 1 and ",length(dim(Y))))

  pten <- array(0,dim=dim(Y),dimnames=dimnames(Y))
  sten <- sqrt(get.slab.var(Y,o))
  if(sum(sten == 0))
  {
    warning("normSlabScaling: Some Slabs have ZERO variance. Ignoring these in the normalization. You should consider removing these slabs from the data.")
    sten[sten == 0] = 1
  }

  if(o == 1)
    for(j in 1:dim(Y)[1])
      pten[j,,] = Y[j,,]/sten[j]

  if(o == 2)
    for(j in 1:dim(Y)[2])
      pten[,j,] = Y[,j,]/sten[j]

  if(o == 3)
    for(j in 1:dim(Y)[3])
      pten[,,j] = Y[,,j]/sten[j]

  return(list(data = pten, pre=list(scales=sten,scalingMode=o)))
}

#' Preprocessing: fiber Centering
#'
#' \code{normFiberCentering} center the fibers of the \eqn{o^{th}} mode of the tensor to zero mean.
#'
#' @param Y the tensor data. See function \code{\link{tensorBF}} for
#'   details.
#' @param o the \eqn{o^{th}} (default: 1) mode of the tensor in which the fibers are to be centered to zero mean.
#' @return a list containing the following elements:
#' \item{data}{The data after performing the required centering operation.}
#' \item{pre}{The centering values used for preprocessing.}
#' @references
#' Kolda, Tamara G., and Brett W. Bader. "Tensor decompositions and applications." SIAM review 51.3 (2009): 455-500.
#' @export
#' @examples
#' #Data generation
#' K <- 3
#' X <- matrix(rnorm(20*K),20,K)
#' W <- matrix(rnorm(30*K),30,K)
#' U <- matrix(rnorm(3*K),3,K)
#' Y = 0
#' for(k in 1:K) Y <- Y + outer(outer(X[,k],W[,k]),U[,k])
#'  Y <- Y + array(rnorm(20*30*3),dim=c(20,30,3))
#'
#' #center the fibers in first mode of tensor Y
#' res <- normFiberCentering(Y=Y,o=1)
#' dim(res$data) #the centered data
normFiberCentering <- function(Y,o)
{
  if(o<1 || o>length(dim(Y)) || !is.numeric(o))
    stop(paste0("In-appropriate value of fiber centering mode specified. The value should be between 1 and ",length(dim(Y))))

  pten <- array(0,dim=dim(Y),dimnames=dimnames(Y))
  mten = 0
  if(o == 1)
  {
    mten = apply(Y,c(2,3),mean.fun)
    for(i in 1:dim(Y)[1])
      pten[i,,] = Y[i,,] - mten
  }

  if(o == 2)
  {
    mten = apply(Y,c(1,3),mean.fun)
    for(i in 1:dim(Y)[2])
      pten[,i,] = Y[,i,] - mten
  }

  if(o == 3)
  {
    mten = apply(Y,c(1,2),mean.fun)
    for(i in 1:dim(Y)[3])
      pten[,,i] = Y[,,i] - mten
  }

  return(list(data = pten, pre = list(centers=mten, centeringMode=o)))
}

#' Postprocessing: Undo Slab Scaling
#'
#' \code{undoSlabScaling} reverts the slabs of the \eqn{o^{th}} mode to undo the scaling effect.
#'
#' @param Yn the normalized tensor data. This can be, for example, the output of \code{\link{reconstructTensorBF}}.
#' @param pre The scaling values and mode used for preprocessing in the format as produced by \code{\link{normSlabScaling}}.
#' @return The data tensor after reversing the scaling operation.
#' @references
#' Kolda, Tamara G., and Brett W. Bader. "Tensor decompositions and applications." SIAM review 51.3 (2009): 455-500.
#' @export
#' @examples
#' #Given tensor Y
#' \dontrun{Yscaled <- normSlabScaling(Y=Y,o=2)}
#' \dontrun{Yunscaled <- undoSlabScaling(Yscaled$data,Yscaled$pre)}
undoSlabScaling <- function(Yn,pre)
{
  if(dim(Yn)!=3)
    stop("Yn should be a tensor.")
  sten = pre$scales
  o = pre$scalingMode

  if(o == 1)
    for(j in 1:dim(Yn)[1])
      Yn[j,,] = Yn[j,,]*sten[j]

  if(o == 2)
    for(j in 1:dim(Yn)[2])
      Yn[,j,] = Yn[,j,]*sten[j]

  if(o == 3)
    for(j in 1:dim(Yn)[3])
      Yn[,,j] = Yn[,,j]*sten[j]

  return(Yn)
}

#' Postprocessing: Undo fiber Centering
#'
#' \code{undoFiberCentering} reverts the fiber's of the \eqn{o^{th}} mode to undo the centering effect.
#'
#' @param Yn the normalized tensor data. This can be, for example, the output of \code{\link{reconstructTensorBF}}.
#' @param pre The centering parameters used for preprocessing in the format as produced by \code{\link{normFiberCentering}}.
#' @return The data tensor after reversing the centering operation.
#' @references
#' Kolda, Tamara G., and Brett W. Bader. "Tensor decompositions and applications." SIAM review 51.3 (2009): 455-500.
#' @export
#' @examples
#' #Given tensor Y
#' \dontrun{Ycentered <- normFiberCentering(Y=Y,o=1)}
#' \dontrun{Yuncentered <- undoFiberCentering(Ycentered$data,Ycentered$pre)}
undoFiberCentering <- function(Yn,pre)
{
  if(dim(Yn)!=3)
    stop("Yn should be a tensor.")

  mten = pre$centers
  o = pre$centeringMode
  if(o == 1)
  {
    for(i in 1:dim(Yn)[1])
      Yn[i,,] = Yn[i,,] + mten
  }

  if(o == 2)
  {
    for(i in 1:dim(Yn)[2])
      Yn[,i,] = Yn[,i,] + mten
  }

  if(o == 3)
  {
    for(i in 1:dim(Yn)[3])
      Yn[,,i] = Yn[,,i] + mten
  }

  return(Yn)
}

KhatriRao.reshaped <- function(Z,U)
{
  if(ncol(U) != ncol(Z))
  {
    stop("Number of columns not same. Khatri Rao product is not possible.")
  }

  N = nrow(Z)
  L = nrow(U)
  K = ncol(U)
  KR <- array(0, dim=c(N,L,K))

  for(k in 1:K)
    KR[,,k] <- Z[,k] %*% t(U[,k])

  return(KR)
}
