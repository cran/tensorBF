noiseProportion <- function(Y,opts,prop,conf,fiberCentering=1,slabScaling=2)
{
  if(!is.numeric(prop) || prop < 0.001 || prop > 1)
    stop("In-appropriate value of noise prop specified. Value should be between 0.001 and 1.")
  if(!is.numeric(conf) || conf < 0.001 || conf > 100)
    stop("In-appropriate value of conf specified in function noiseProportion. Value should be between 0.001 and 100.")

  if(!is.null(fiberCentering))
    Y <- normFiberCentering(Y,fiberCentering)$data;
  if(!is.null(slabScaling))
    Y <- normSlabScaling(Y,slabScaling)$data

  prop.to.be.explained.by.noise <- prop;
  tau_setting = "S";
  Y = list(Y);

  D <- unlist(lapply(Y,function(x){dim(x)[2]}))
  N <- unlist(lapply(Y,function(x){dim(x)[1]}))
  L <- unlist(lapply(Y,function(x){dim(x)[3]})); L[is.na(L)] = 1
  M <- length(Y)

  if(tau_setting == "S")
    prior.alpha_0t <- conf*N*L*D/2
  else
    prior.alpha_0t <- rep(conf*N*L/2, M)

  view.var.modes <- rep(NA, M)
  for (i in 1:M) {
    total.variance.in.view <- (N*L*D[i]*var(Y[[i]][!is.na(Y[[i]])]))/2 #sum(Y[[i]][!is.na(Y[[i]])]^2)/2
    if(tau_setting != "S")
      total.variance.in.view <- total.variance.in.view/D[i]

    sigma.hat <- prop.to.be.explained.by.noise * total.variance.in.view
    view.var.modes[i] <- sigma.hat
  }

  prior.beta_0t <- conf * view.var.modes
  opts$prior.alpha_0t=prior.alpha_0t
  opts$prior.beta_0t=prior.beta_0t
  return(opts)
}
