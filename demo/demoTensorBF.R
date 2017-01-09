
require("tensorBF")
#' Demo of Bayesian Tensor Factorization
#'
#' \code{demo} demonstrates an example application of tensorBF package on simulated data
#'
#' @return res A list containing the following :
#'   \item{Y}{the simulated data}
#'   \item{res}{the model fit using \code{\link{tensorBF}}}
#'
demoTensorBF <- function()
{
  #Data generation
  K <- 3
  N = 20; D = 30; L = 3
  print(""); print("")
  print("*************************************************")
  print(paste0("Generate Data with N:",N," D:",D," L:",L," K:",K))

  X <- matrix(rnorm(N*K),N,K)
  W <- matrix(rnorm(D*K),D,K)
  U <- matrix(rnorm(L*K),L,K)
  Y = 0
  for(k in 1:K) Y <- Y + outer(outer(X[,k],W[,k]),U[,k])
   Y <- Y/sd(Y)
   Y <- Y + array(rnorm(N*D*L,0,0.5),dim=c(N,D,L))
  dimnames(Y) = list(paste0("row.",1:N),paste0("col.",1:D),paste0("lat.",1:L))

  print(""); print("")
  print("*************************************************")
  print("Inserting Missing Values in the tensor")
  missing.inds = sample(prod(dim(Y)),100)
  Yobs = Y[missing.inds]
  Y[missing.inds] = NA

  print(""); print("")
  print("*************************************************")
  print("Running the model with default options")
  res <- tensorBF(Y=Y)
  print(""); print("")
  print("*************************************************")
  print("Model run completed")

  print(""); print("")
  print("*************************************************")
  print("Predict Missing Values and plot the predictions")
  pred = predictTensorBF(Y=Y,res=res)
  plot(Yobs,pred[missing.inds],xlab="obs",ylab="pred",main=round(cor(Yobs,pred[missing.inds]),2))
  print(paste0("Prediction Correlation: ",round(cor(Yobs,pred[missing.inds]),2)))

  print(""); print("")
  print("*************************************************")
  if(res$K>=1){
  print("Plot Components")
  for(k in 1:res$K)
    plotTensorBF(res = res,Y=Y,k=k)
  } else print("No active component detected.")
}

demoTensorBF()
