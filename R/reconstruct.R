#' @title Reconstruct the data based on posterior samples
#'
#' @description
#' \code{reconstructTensorBF} returns the reconstruction of the data based on
#' posterior samples of a given run. The function reconstructs the tensor for
#' each posterior sample and then computes the expected value.
#' The reconstruction is returned in the un-normalized space if \code{res$pre}
#' contains appropriate preprocessing information.
#'
#' @param res The model object from function \code{\link{tensorBF}}.
#' @return The reconstructed data, a tensor of the size equivalent to the
#' data on which the model was run.
#' @export
#' @examples
#' #Data generation
#' K <- 3
#' X <- matrix(rnorm(20*K),20,K)
#' W <- matrix(rnorm(30*K),30,K)
#' U <- matrix(rnorm(3*K),3,K)
#' Y = 0
#' for(k in 1:K) Y <- Y + outer(outer(X[,k],W[,k]),U[,k])
#'  Y <- Y + array(rnorm(20*30*3,0,0.25),dim=c(20,30,3))
#'
#' #Run the method with default options and reconstruct the model's representation of the tensor
#' \dontrun{res <- tensorBF(Y)}
#' \dontrun{recon = reconstructTensorBF(res)}
#' \dontrun{inds = sample(prod(dim(Y)),100)}
#' \dontrun{plot(Y[inds],recon[inds],xlab="obs",ylab="recon",main=round(cor(Y[inds],recon[inds]),2))}

reconstructTensorBF <- function(res)
{
  recon = remake.tensor.samples.EV(res)[[1]]
  if(is.null(res$pre)) {
    #print("reconstructTensorBF: Preprocessing information not specified in res$pre, returning reconstruction without un-normalization.")
  }else{
    if(!is.null(res$pre$scales) && !is.null(res$pre$scalingMode)) {
      recon <- undoSlabScaling(recon,res$pre);
    }

    if(!is.null(res$pre$centers) && !is.null(res$pre$centeringMode)){
      recon <- undoFiberCentering(recon,res$pre);
    }
  }
  return(recon)
}

#
# Remake the tensor for each posterior sample and then take the expected value. Uses all posterior samples.
#
remake.tensor.samples.EV <- function(model,IsMat=-999)
{
  samples = length(model$posterior$X)
  Rec <- list()
  for(s in 1:samples)
  {
    md = getPosteriorSample(model$posterior,s)
    if(s == 1)
      Rec = remake.tensor.mat(md,IsMat)
    else
    {
      RecSamp = remake.tensor.mat(md,IsMat)
      for(m in 1:length(Rec))
        Rec[[m]] = Rec[[m]] + RecSamp[[m]]
    }
    if((s %% 100) == 0)
      cat(".",append=TRUE)
  }
  for(m in 1:length(Rec))
    Rec[[m]] = Rec[[m]]/samples
  return(Rec)
}

remake.tensor.mat <- function(model,IsMat,kk="all")
{
  return(make.tensor.mat(model,IsMat,kk))
}

make.tensor.mat <- function(model,IsMat,kk="all")
{
  K = ncol(model$X)
  if(kk=="all")
    kk = 1:K

  M <- length(model$W)
  if(IsMat[1]==-999) IsMat = rep(0,M)
  Yestim <- list()
  ## Can be speeded up - see proof tensor
  for(m in 1:M)
  {
    Yestim[[m]] <- 0

    if(IsMat[m]==0)
    {
      for(k in kk)
        Yestim[[m]] <- Yestim[[m]] + outer(outer(model$X[,k],model$W[[m]][,k]),model$U[,k])
    }
    else
    {
      for(k in kk)
        Yestim[[m]] <- Yestim[[m]] + outer(model$X[,k],model$W[[m]][,k])
    }
  }

  return(Yestim)
}

getPosteriorSample <- function(post,s)
{
  samp <- list()
  samp$X <- post$X[[s]]
  samp$U <- post$U[[s]]
  samp$W <- post$W[[s]]
  samp$tau <- post$tau[[s]]
  samp$alpha <- post$alpha[[s]]
  samp$alphaU <- post$alphaU[[s]]
  samp$Z <- post$Z[[s]]
  return(samp)
}
