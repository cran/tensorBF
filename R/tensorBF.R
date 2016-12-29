library(tensor)

#' Bayesian Factorization of a Tensor
#'
#' \code{tensorBF} implements the Bayesian factorization of a tensor.
#'
#'
#' Bayesian Tensor Factorization performs tri-linear (CP) factorization of a tensor.
#' The method automatically identifies the number of components,
#' given K is initialized to a large enough value, see arguments.
#' Missing values are supported and should be set as NA's in the data.
#' They will not affect the model parameters, and can be predicted
#' with function \code{\link{predictTensorBF}}, based on the observed values.
#'
#' @param Y is a three-mode tensor to be factorized.
#' @param method the factorization method. Currently only "CP" (default) is supported.
#' @param opts List of model options; see function \code{\link{getDefaultOpts}} for details and default.
#' @param K The number of components (i.e. latent variables or factors). Recommended to be
#'   set somewhat higher than the expected component number, so that the method
#'   can determine the model complexity by prunning excessive components
#'   (default: 20\% of the sum of lower two dimensions).
#'   High values result in high CPU time.
#'
#'   NOTE: Adjust parameter noiseProp if sufficiently
#'   large values of K do not lead to a model with pruned components.
#' @param fiberCentering the mode for which fibers are to be centered at zero (default = NULL).
#' Fiber is analogous to a vector in a particular mode.
#' Fiber centering and Slab scaling are the recommended normalizations for a tensor.
#' For details see the provided normalization functions and the references therein.
#' @param slabScaling the mode for which slabs are to be scaled to unit variance (default = NULL).
#' Slab is analogous to a matrix in a particular mode.
#' Alternativly, you can preprocess the data using the provided normalization functions.
#' @param noiseProp c(prop,conf); sets an informative noise prior for tensorBF.
#' The model sets the noise prior such that the expected proportion of
#' variance explained by noise is defined by this parameter. It is recommended when
#' the standard prior from \code{\link{getDefaultOpts}} seems to overfit the model
#' by not prunning any component with high initial K. Use NULL to switch off
#' informative noise prior.
#'
#' - prop defines the proportion of total variance to be explained by noise (between 0.1 and 0.9),
#'
#' - conf defines the confidence in the prior (between 0.1 and 10).
#'
#'
#' We suggest a default value of c(0.5,0.5) for real data sets.
#' @return A list containing model parameters.
#'   For key parameters, the final posterior sample ordered w.r.t. component variance is provided to aid in
#'   initial checks; all the posterior samples should be used for model
#'   analysis. The list elements are:
#' \item{K}{The number of learned components. If this value is not less then the input argument K, the model should be rerun with a larger K or use the noiseProp parameter.}
#' \item{X}{a matrix of \eqn{N \times K} dimensions, containing the last Gibbs sample of the first-mode latent variables.}
#' \item{W}{a matrix of \eqn{D \times K} dimensions, containing the last Gibbs sample of the second-mode latent variables.}
#' \item{U}{a matrix of \eqn{L \times K} dimensions, containing the last Gibbs sample of the third-mode latent variables.}
#' \item{tau}{The last sample of noise precision.}
#' and the following elements:
#' \item{posterior}{the posterior samples of model parameters (X,U,W,Z,tau).}
#' \item{cost}{The likelihood of all the posterior samples.}
#' \item{opts}{The options used to run the model.}
#' \item{conv}{An estimate of the convergence of the model, based on reconstruction
#'  of data using the Geweke diagnostic. Values significantly above 0.05 occur when
#'  model has not converged and should therefore be rerun with a higher value of iter.burnin in \code{\link{getDefaultOpts}}.}
#' \item{pre}{A list of centering and scaling values used to transform the data, if any. Else an empty list.}
#' @export
#' @examples
#' #Data generation
#' K <- 2
#' X <- matrix(rnorm(20*K),20,K)
#' W <- matrix(rnorm(30*K),30,K)
#' U <- matrix(rnorm(3*K),3,K)
#' Y = 0
#' for(k in 1:K) Y <- Y + outer(outer(X[,k],W[,k]),U[,k])
#'  Y <- Y + array(rnorm(20*30*3,0,0.25),dim=c(20,30,3))
#'
#' #Run the method with default options
#' \dontrun{res2 <- tensorBF(Y=Y)}
#'
#' #Run the method with K=3 and iterations=1000
#' opts <- getDefaultOpts(); opts$iter.burnin = 1000
#' res1 <- tensorBF(Y=Y,K=3,opts=opts)
#'
#' #Vary the user defined expected proportion of noise variance
#' #explained. c(0.2,1) represents 0.2 as the noise proportion
#' #and confidence of 1
#' \dontrun{res3 <- tensorBF(Y=Y,noiseProp=c(0.2,1))}
#'
tensorBF <- function(Y,method="CP",K=NULL,opts=NULL,fiberCentering=NULL,slabScaling=NULL,noiseProp=c(0.5,0.5)){

  if(is.null(Y)) stop("Please specify the input tensor Y to factorize")
  if(length(dim(Y))!=3) stop(paste0("Y should be a 3-mode tensor. Current Y is a ",length(dim(Y)),"-mode tensor."))
  if(mean(apply(Y,1:length(dim(Y)),is.numeric))!=1) print("Y contains non-numeric values")

  if("CP"!=method) {
    print("Only CP factorization supported currently. Running with CP factorization."); method = "CP"
  }

  if(is.null(K))
    K <- ceiling(0.2*sum(sort(dim(Y))[1:2]))

  if(is.null(opts) || is(opts)[1] != "list"){
    print("Initializing with default options.")
    opts <- getDefaultOpts()
  }

  pre <- list()
  if(!is.null(fiberCentering)){
    tmp <- normFiberCentering(Y,fiberCentering);
    Y <- tmp$data; pre$centers <- tmp$pre$centers;
    pre$centeringMode <- tmp$pre$centeringMode;
  }
  if(!is.null(slabScaling)) {
    tmp <- normSlabScaling(Y,slabScaling);
    Y <- tmp$data; pre$scales <- tmp$pre$scales;
    pre$scalingMode <- tmp$pre$scalingMode;
  }
  if(!is.null(noiseProp)){
    if(length(noiseProp)==2)
      opts <- noiseProportion(Y,opts,prop=noiseProp[1],conf=noiseProp[2],fiberCentering=NULL,slabScaling=NULL)
    else
      stop(paste0("noiseProp should contain prop and conf in the format c(prop,conf). User specified ",length(noiseProp)," values."))
  }

  res <- tensorBF.compute(Y=list(Y),0,K,opts);
  res$pre <- pre
  return(res)
}

#' Predict Missing Values using the Bayesian tensor factorization model
#'
#' \code{predictTensorBF} predicts the missing values in the data \code{Y} using the learned model \code{res}.
#'
#' If the original data \code{Y} contained missing values (NA's),
#' this function predicts them using the model. The predictions are
#' returned in the un-normalized space if \code{res$pre} contains appropriate
#' preprocessing information.
#'
#' @param Y is a 3-mode tensor containing missing values as NA's. See function \code{\link{tensorBF}} for details.
#' @param res the model object returned by the function \code{\link{tensorBF}}.
#' @return A tensor of the same size as \code{Y} containing predicted values in place of NA's.
#' @export
#' @examples
#' #Data generation
#' \dontrun{K <- 2}
#' \dontrun{X <- matrix(rnorm(20*K),20,K)}
#' \dontrun{W <- matrix(rnorm(30*K),30,K)}
#' \dontrun{U <- matrix(rnorm(3*K),3,K)}
#' \dontrun{Y = 0}
#' \dontrun{for(k in 1:K) Y <- Y + outer(outer(X[,k],W[,k]),U[,k])}
#' \dontrun{ Y <- Y + array(rnorm(20*30*3,0,0.25),dim=c(20,30,3))}
#'
#' #insert missing values
#' \dontrun{m.inds = sample(prod(dim(Y)),100)}
#' \dontrun{Yobs = Y[m.inds]}
#' \dontrun{Y[m.inds] = NA}
#'
#' #Run the method with default options and predict missing values
#' \dontrun{res <- tensorBF(Y)}
#' \dontrun{pred = predictTensorBF(Y=Y,res=res)}
#' \dontrun{plot(Yobs,pred[m.inds],xlab="obs",ylab="pred",main=round(cor(Yobs,pred[m.inds]),2))}

predictTensorBF <- function(Y,res){
  recon <- reconstructTensorBF(res)
  recon[!is.na(Y)] = Y[!is.na(Y)]
  return(recon)
}

tensorBF.compute <- function(Y,IsMat,K,opts){
  if(!is.numeric(K)) stop("K must be a numeric value.")
  if(K<1) stop(paste0("K =",K,". Too small value specified. Increase K and rerun."))

  M <- length(Y)
  N <- dim(Y[[1]])[1]
  if(sum(IsMat==0) == 0)
	  L <- 1
  else
	  L <- dim(Y[[which(!IsMat)[1]]])[3]
  D <- unlist(lapply(Y,function(x){dim(x)[2]}))
  Ls <-unlist(lapply(Y,function(x){dim(x)[3]})); Ls[is.na(Ls)] = 1

  const <- - N*sum(Ls*D)*log(2*pi)/2
  id <- rep(1,K)
  X <- matrix(rnorm(N*K),N,K)
  covX <- diag(K)
  XX <- crossprod(X)
  XXtmp <- XX
  XXdiag <- rep(1,K)

  if(sum(IsMat==0) == 0)
  	U <- matrix(1,L,K)
  else
  	U <- matrix(rnorm(L*K),L,K)
  covU <- diag(K)
  UU <- crossprod(U)
  UUdiag <- rep(1,K)

  alpha_0 <- opts$prior.alpha_0 # ARD hyperparameters
  beta_0 <- opts$prior.beta_0
  Z <- matrix(1,M,K) # binary mask
  alpha_0t <- opts$prior.alpha_0t
  beta_0t <- opts$prior.beta_0t
  if(length(alpha_0t)<M) alpha_0t = rep(alpha_0t,M)
  if(length(beta_0t)<M) beta_0t = rep(beta_0t,M)
  prior.betaW1 <- opts$prior.betaW1
  prior.betaW2 <- opts$prior.betaW2
  ARDW <- opts$ARDW
  ARDU <- opts$ARDU
  ARDX <- opts$ARDX
  TAU_SETTING <- "S"
  DEBUGMODE <- FALSE
  VERBOSE <- opts$verbose

  alphaU <- matrix(K*ARDU+1,L,K)
  alphaX <- matrix(K*ARDX+1,N,K)
  alpha <- vector("list")
  for(m in 1:M)
    alpha[[m]] <- matrix(K*ARDW+1,D[m],K)

  tau <- vector("list")
  b_tau <- vector("list")
  for(m in 1:M)
  {
    	tau[[m]] <- matrix(rep(opts$init.tau,Ls[m]*D[m]),Ls[m],D[m],byrow=TRUE)
    	b_tau[[m]] <- matrix(rep(1,Ls[m]*D[m]),Ls[m],D[m],byrow=TRUE) # initialization
  }
  a_tau <- N/2#(N*L)/2 #N/2#

  W <- vector("list",length=M)
  for(m in 1:M)
    W[[m]] <- matrix(0,D[m],K)

  r <- rep(0.5,K)

  iter.burnin <- opts$iter.burnin
  iter.sampling <- opts$iter.sampling
  iter.thinning <- opts$iter.thinning
  maxiter <- iter.burnin + iter.sampling*iter.thinning
  if(VERBOSE >= 0) print(paste0("Initializing the model with ",K," components."))
  posterior <- list();
  if(iter.sampling>0)# && length(postVar)>0)
  {
	  posterior$tau <- list(); length(posterior$tau) <- ceiling(iter.sampling/iter.thinning)
	  posterior$U <- list(); length(posterior$U) <- ceiling(iter.sampling/iter.thinning)
	  posterior$W <- list(); length(posterior$W) <- ceiling(iter.sampling/iter.thinning)
	  posterior$X <- list(); length(posterior$X) <- ceiling(iter.sampling/iter.thinning)
	  posterior$Z <- list(); length(posterior$Z) <- ceiling(iter.sampling/iter.thinning)
	  #posterior$cost <- rep(NA,ceiling(iter.sampling/iter.thinning))
  }

  # the component numbers and loglikelihood
  cost <- rep(0,maxiter)

  if(K>300 & VERBOSE>0) {
    print("Computational Costs may be high for K>300. Consider noiseProp parameter to set an informative noise prior and lower values of K.")
  }
  ##Missing Values
  missingValues = FALSE
  na.inds = list()
  for(m in 1:M)
  {
  	inds = which(is.na(Y[[m]]))
  	if(length(inds)>0)
  	{
  		na.inds[[m]] = inds
  		missingValues = TRUE
  	}
  }
  na.inds[M+1] = 0
  missingValues.InitIter = round(iter.burnin/2,0) #5000
  if(missingValues && VERBOSE>0){
    str="Missing Values Detected"
    if(missingValues.InitIter<=iter.burnin)
      str = paste0(str,", Learning using EM type approximation after ",missingValues.InitIter," iterations.")
      print(str)
  }
  ## Missing Values end

  for(iter in 1:maxiter){
  if(VERBOSE==2)
    print(paste("iter: ",iter,sep=""))
    #
    # sample Z,U and W
    #

  if(iter > 1){
	XXUU <- XX*UU#DONE
	XXUUdiag <- diag(XXUU)
	diag(XXUU) <- 0 #for j!=k matrix summations
	XXtmp <- XX
	XXdiag <- diag(XXtmp)
	diag(XXtmp) <- 0 #for j!=k matrix summations

	for(m in 1:M){ #DONE
	for(k in 1:K){

	if(IsMat[m]) lambda <- tau[[m]][1,]*XXdiag[k] + alpha[[m]][,k] #DONE
	else lambda <- tau[[m]][1,]*XXUUdiag[k] + alpha[[m]][,k]

	if(missingValues && (missingValues.InitIter >= iter) )
	{
		mu_sub <- 0
		if(IsMat[m])
		{
			UX = X; UXK = UX[,k];  UX[,k] = 0
			ss <- tcrossprod(UX,W[[m]])
			tmp = Y[[m]]-ss; tmp[is.na(tmp)] = 0;
			mu_sub <- mu_sub + crossprod(tmp,UXK)
		}
		else{
		for(l in 1:L)
		{
			UX = X*matrix(rep(U[l,],N),N,K,byrow=T); UXK = UX[,k];  UX[,k] = 0
			ss <- tcrossprod(UX,W[[m]])
			tmp = Y[[m]][,,l]-ss; tmp[is.na(tmp)] = 0;
			mu_sub <- mu_sub + crossprod(tmp,UXK)
		}}
		mu = mu_sub*tau[[m]][1,1]/lambda
	}
	else
	{
		if(IsMat[m]) mu <- (tau[[m]][1,]/lambda)*(crossprod(Y[[m]],X[,k]) - W[[m]]%*%(XXtmp[k,]))
		else mu <- (tau[[m]][1,]/lambda)*(crossprod(tensor::tensor(Y[[m]],U[,k],3,1),X[,k]) - W[[m]]%*%(XXUU[k,]))
	}
	#logpr of activating the component
      	logpr <- 0.5*( sum(log(alpha[[m]][,k]) - log(lambda)) + crossprod(mu*sqrt(lambda)) ) + log(r[k]) - log(1-r[k]) # m - k
	zone <- 1/(1+exp(-logpr))
	if(DEBUGMODE)
		debug$zone[m,k] = logpr

	if(iter > 500){
		a = Z[m,k]
	  	Z[m,k] <- as.double((runif(1) < zone))
	}

	if(Z[m,k]==1){
		W[[m]][,k] <- mu + 1/sqrt(lambda)*rnorm(D[m])
	}else{
		W[[m]][,k] <- 0
	}
	}
	}

	r <- rep( (prior.betaW1+sum(Z))/(M*K+prior.betaW1+prior.betaW2), K)
    }else{
	# alternative is to use method similar to VB
	# note now this is more connected to initialization
	for(m in 1:M){
	if(missingValues && (missingValues.InitIter < iter))
	{
	      ##Could give bad init when high missing values. therefore skip from here.
	      if(length(na.inds[[m]])>0)
	      Y[[m]][na.inds[[m]]] = mean(Y[[m]][-na.inds[[m]]])
	}
	if(IsMat[m]) covW <- diag(1,K) + opts$init.tau*XX
	else covW <- diag(1,K) + opts$init.tau*XX*UU
        eSW <- eigen(covW)
        covW <- tcrossprod(eSW$vectors*outer(id,1/eSW$values),eSW$vectors)
	estimW = matrix(0,D[m],K) #equivalent of crossprod(Y[[m]],X)
	for(k in 1:K){
		if(missingValues && (missingValues.InitIter >= iter))
		{
			tmp = Y[[m]]; tmp[is.na(tmp)] = 0;
			if(IsMat[m]) estimW[,k] = crossprod(tmp,X[,k])
			else estimW[,k] = crossprod(tensor::tensor(tmp,U[,k],3,1),X[,k])
		}
		else{
			if(IsMat[m]) estimW[,k] = crossprod(Y[[m]],X[,k])
			else estimW[,k] = crossprod(tensor::tensor(Y[[m]],U[,k],3,1),X[,k])
		}
	}
	W[[m]] <- estimW%*%covW*opts$init.tau + matrix(rnorm(D[m]*K),D[m],K)%*%t( eSW$vectors*outer(id,1/sqrt(eSW$values)) )
	}

  	#Initialize alphaU here
  	if(TRUE == ARDU)
  	{
    	##Elementwise Alpha
    	alphaU <- apply((U^2)/2,c(1,2),function(x){rgamma(1,shape=alpha_0+1/2,rate=beta_0+x)})
  	}
  	#else
  	#{
    	##Viewwise Alpha
    	#tmp <- apply(U^2,2,function(x){rgamma(1,shape=alpha_0+1/2,rate=beta_0+sum(x)/2)})
    	#alphaU <- matrix(rep(tmp,L),L,K,byrow=TRUE)
  	#}

    if(TRUE == ARDX)
    {
      alphaX <- apply((X^2)/2,c(1,2),function(x){rgamma(1,shape=alpha_0+1/2,rate=beta_0+x)})
    }
    #else
    #{
      #tmp <- apply(X^2,2,function(x){rgamma(1,shape=alpha_0+1/2,rate=beta_0+sum(x)/2)})
      #alphaX <- matrix(rep(tmp,N),N,K,byrow=TRUE)
    #}

  }# end initialization else

################################################################################
	#
	# sample X (latent variables)
	#
################################################################################

	#Sample X_n Matrix Version
	X <- 0
	if(missingValues && (missingValues.InitIter >= iter) )
	{
		for(m in 1:M){
		tmp = Y[[m]]; tmp[is.na(tmp)] = 0;
		if(IsMat[m]) X <- X + tmp %*% (W[[m]]*tau[[m]][1,1])
		else{
			for(l in 1:L) X <- X + ((tmp[,,l]*matrix(rep(tau[[m]][1,],N),N,D[m],byrow=TRUE)) %*% (W[[m]]*matrix(rep(U[l,],D[m]),D[m],ncol(U),byrow=TRUE)))
		}
		}
	}
	else
	{
		for(m in 1:M){
		if(IsMat[m]) X <- X + Y[[m]] %*% (W[[m]]*tau[[m]][1,1])
		else{
			for(l in 1:L) X <- X + ((Y[[m]][,,l]*matrix(rep(tau[[m]][1,],N),N,D[m],byrow=TRUE)) %*% (W[[m]]*matrix(rep(U[l,],D[m]),D[m],ncol(U),byrow=TRUE)))
		}
		}
	}
	covX.bk <- 0
	for(m in 1:M){
	  if(IsMat[m]) covX.bk <- covX.bk + crossprod(W[[m]]*sqrt(tau[[m]][1,]))
	  else covX.bk <- covX.bk + crossprod(W[[m]]*sqrt(tau[[m]][1,]))*UU
	}
	if(TRUE == ARDX)
	{
	  for(n in 1:N)
	  {
	    covX <- covX.bk + diag(alphaX[n,]) #ARD Prior
	    eS <- eigen(covX)
	    covX <- tcrossprod( eS$vectors*outer(id,1/sqrt(eS$values)) )
	    X[n,] <- X[n,]%*%covX + matrix(rnorm(K),1,K)%*%t( eS$vectors*outer(id,1/sqrt(eS$values)) )
	  }
	} else {
	  #Gaussian Prior with variance I on X
	  covX <- covX.bk + diag(K)
	  eS <- eigen(covX)
	  covX <- tcrossprod( eS$vectors*outer(id,1/sqrt(eS$values)) )
	  X <-  X%*%covX + matrix(rnorm(N*K),N,K)%*%t( eS$vectors*outer(id,1/sqrt(eS$values)) )
	}

	# do scaling
	Rx <- apply(X,2,sd)
	X <- X*outer(rep(1,N),1/Rx)
	for(m in 1:M)
		W[[m]] <- W[[m]]*outer(rep(1,D[m]),Rx)
	XX <- crossprod(X)

################################################################################
	#
	# sample U
	#
################################################################################

  if(sum(IsMat==0) != 0) #Update if there is atleast one tensor, else ignore U
  {
	U <- 0
	if(missingValues && (missingValues.InitIter >= iter) )
	{
		for(m in 1:M){
		if(!IsMat[m]){
			tmp = Y[[m]]; tmp[is.na(tmp)] = 0;
			for(n in 1:N){
		      		U <- U + ((t(tmp[n,,])*matrix(rep(tau[[m]][1,],L),L,D[m],byrow=TRUE)) %*% (W[[m]]*matrix(rep(X[n,],D[m]),D[m],ncol(X),byrow=TRUE)))
			}
		}}
	}
	else
	{
		for(m in 1:M){
			if(!IsMat[m]) {
			for(n in 1:N) U <- U + ((t(Y[[m]][n,,])*matrix(rep(tau[[m]][1,],L),L,D[m],byrow=TRUE)) %*% (W[[m]]*matrix(rep(X[n,],D[m]),D[m],ncol(X),byrow=TRUE))) }
		}
	}
	covU.bk <- 0
	for(m in 1:M){
	  if(!IsMat[m]) covU.bk <- covU.bk + crossprod(W[[m]]*sqrt(tau[[m]][1,]))*XX
	}
	if(TRUE == ARDU)
	{
		for(l in 1:L)
		{
		covU <- covU.bk + diag(alphaU[l,]) #ARD Prior
		eS <- eigen(covU)
		covU <- tcrossprod( eS$vectors*outer(id,1/sqrt(eS$values)) )
		U[l,] <-  U[l,]%*%covU + matrix(rnorm(K),1,K)%*%t( eS$vectors*outer(id,1/sqrt(eS$values)) )
		}
	} else {
	    #Gaussian Prior with variance I on U
	    covU <- covU.bk + diag(K)
	    eS <- eigen(covU)
	    covU <- tcrossprod( eS$vectors*outer(id,1/sqrt(eS$values)) )
	    U <-  U%*%covU + matrix(rnorm(L*K),L,K)%*%t( eS$vectors*outer(id,1/sqrt(eS$values)) )
	}


	# do scaling
	if(L == 1)
	{
		Ru <- as.vector(U)
		U <- U*outer(rep(1,L),1/Ru)
		U[is.na(U)] = 0
		for(m in 1:M) #DONE, not in all M
			if(!IsMat[m]) W[[m]] <- W[[m]]*outer(rep(1,D[m]),Ru)
	}
	else
	{
		Ru <- apply(U,2,sd)
		U <- U*outer(rep(1,L),1/Ru)
		U[is.na(U)] = 0
		for(m in 1:M)
			if(!IsMat[m]) W[[m]] <- W[[m]]*outer(rep(1,D[m]),Ru)
	}
	UU <- crossprod(U)
   }
################################################################################
	#
	# Sample Hyperparameters
	#
################################################################################

	if(TRUE == ARDU)
	{
	  #Elementwise alpha
 	  alphaU <- U^2/2
  	keep <- which(alphaU!=0)
  	for(n in keep) alphaU[n] <- rgamma(1,shape=alpha_0+0.5,rate=beta_0+alphaU[n])+ 1e-7 #+ 1e-7 to prevent bad values of gamma
  	alphaU[-keep] <- rgamma( length(alphaU[-keep]),shape=alpha_0,rate=beta_0)+ 1e-7
	}
	#else
	#{
	  #View-Wise alpha
	  #tmp <- apply(U^2,2,function(x){rgamma(1,shape=alpha_0+L/2,rate=beta_0+sum(x)/2)})
    #alphaU <- matrix(rep(tmp,L),L,K,byrow=TRUE)+ 1e-7
	#}

	if(TRUE == ARDX)
	{
	  alphaX <- X^2/2
	  keep <- which(alphaX!=0)
	  for(n in keep) alphaX[n] <- rgamma(1,shape=alpha_0+0.5,rate=beta_0+alphaX[n])+ 1e-7 #+ 1e-7 to prevent bad values of gamma
	  alphaX[-keep] <- rgamma( length(alphaX[-keep]),shape=alpha_0,rate=beta_0)+ 1e-7
	}

	for(m in 1:M)
	{
	  if(TRUE == ARDW)
	  {
      #Elementwise alpha
  	  alpha[[m]] <- W[[m]]^2/2
  	  keep <- which(alpha[[m]]!=0)
  	  for(n in keep) alpha[[m]][n] <- rgamma(1,shape=alpha_0+0.5,rate=beta_0+alpha[[m]][n])+ 1e-7
  	  alpha[[m]][-keep] <- rgamma( length(alpha[[m]][-keep]),shape=alpha_0,rate=beta_0)+ 1e-7
	  }
	  #else
	  #{
  	  #View-Wise alpha
	    #tmp <- apply(W[[m]]^2,2,function(x){rgamma(1,shape=alpha_0+(as.double(sum(x)!=0)*D[m])/2,rate=beta_0+sum(x)/2)})
      #alpha[[m]] <- matrix(rep(tmp,D[m]),D[m],K,byrow=TRUE)+ 1e-7
	  #}
	}


	#
	# sample noise
	#
	XU <- KhatriRao.reshaped(X,U)
	for(m in 1:M){
	if(IsMat[m]){ #DONE
		XW <- tcrossprod(X,W[[m]])
		if(missingValues && (missingValues.InitIter > iter) )
		{
			XW <- Y[[m]] - XW
			XW[is.na(XW)] = 0;
			b_tau[[m]][1,] <- 0.5*crossprod((XW)^2,rep(1,dim(Y[[m]])[1]))
		}
		else{
			if(missingValues)
			{
				if(length(na.inds[[m]])>0)
					Y[[m]][na.inds[[m]]] = XW[na.inds[[m]]]
			}
			b_tau[[m]][1,] <- 0.5*crossprod((Y[[m]] - XW)^2,rep(1,dim(Y[[m]])[1]))
		}

		#Vector b_tau
		if("S" == TAU_SETTING)
		{
			a_tau <- (N*Ls[m]*D[m])/2 #DONE
			b_tau2 <- sum(b_tau[[m]][1,])
			b_tau[[m]] <- matrix(b_tau2,Ls[m],D[m])
			tau[[m]] <- matrix(rgamma(1,shape= alpha_0t[m]+ a_tau,rate= beta_0t[m] + b_tau2),Ls[m],D[m])
		}
	}
	else {
		#single valued tau
		#b_tau[m] <- ( tau_tmp[m] - 2*sum( W[[m]]*crossprod(Y[[m]],X) ) +
		#                sum((W[[m]]%*%XX)*W[[m]]) )/2

		#Matrix b_tau generation of D values (not LD) is correct, verified computationally
		XWU <- aperm(tensor::tensor(XU,W[[m]],3,2),c(1,3,2))
		if(missingValues && (missingValues.InitIter > iter) ) #in the last iter of ignoring missing values, update Y[[m]] with the model updates
		{
			XWU <- Y[[m]] - XWU
			XWU[is.na(XWU)] = 0;
			b_tau[[m]][1,] <- 0.5*tensor::tensor( tensor::tensor((XWU)^2,rep(1,dim(Y[[m]])[3]),3,1) ,rep(1,dim(Y[[m]])[1]),1,1)
		}
		else{
			if(missingValues)
			{
				if(length(na.inds[[m]])>0)
					Y[[m]][na.inds[[m]]] = XWU[na.inds[[m]]]
			}
			b_tau[[m]][1,] <- 0.5*tensor::tensor( tensor::tensor((Y[[m]] - XWU)^2,rep(1,dim(Y[[m]])[3]),3,1) ,rep(1,dim(Y[[m]])[1]),1,1)
		}

		#Vector b_tau
		if("S" == TAU_SETTING)
		{
			a_tau <- (N*Ls[m]*D[m])/2
			b_tau2 <- sum(b_tau[[m]][1,])
			b_tau[[m]] <- matrix(b_tau2,Ls[m],D[m])
			tau[[m]] <- matrix(rgamma(1,shape= alpha_0t[m]+ a_tau,rate= beta_0t[m] + b_tau2),Ls[m],D[m])
		}
   }
  }

    # calculate log likelihood
    cost[iter] <- const
    for(m in 1:M)
    {
      if("S" == TAU_SETTING)
      	cost[iter] <- cost[iter] + N*L*D[m]*0.5*log(tau[[m]][1,1]) - b_tau[[m]][1,1]*tau[[m]][1,1] # Will change when tau reduced to D vector
      if("D" == TAU_SETTING)
      	cost[iter] <- cost[iter] + N*L*0.5*sum(log(tau[[m]][1,])) - sum(b_tau[[m]][1,]*tau[[m]][1,])
      if("LD" == TAU_SETTING)
      	cost[iter] <- cost[iter] + N*0.5*sum(log(tau[[m]])) - sum(b_tau[[m]]*tau[[m]])
    }

    # print component number and collect numbers
    if(VERBOSE==2 || ((VERBOSE==1) && ((iter %% 500) == 0)))
    {
	    str = paste0("Iter: ",iter," Components: ",rowSums(Z)," Empty: ",sum(colSums(Z) == 0)," Cost: ",round(cost[iter],0))
	    print(str)
    }

  if(iter > iter.burnin && (((iter - iter.burnin) %% iter.thinning) == 0))
  {
  	ind <- (iter - iter.burnin)/iter.thinning
  	posterior$tau[[ind]] <- list(tau[[1]][1,1]) #for Bayesian CP
  	posterior$U[[ind]] <- U
  	posterior$W[[ind]] <- W
  	posterior$X[[ind]] <- X
  	posterior$Z[[ind]] <- Z

  	if(ind==1 && (sum(colSums(Z) == 0) == 0) )
  	{
  	  warning("Burn-in period finished, but no components have been shut down, preventing automatic model complexity selection. Consider higher inital K, if computational resources permit. Alternativly, refer to noiseProp parameter of the tensorBF function for defining an informative prior for residual noise.")
  	}

  	if(ind==1 && (sum(Z) == 0))
  	{
  	  warning("Burn-in period finished, and all components have been shut down. Returning a NULL model.")
  	}

  	if(ind==1 && (sum(colSums(Z) == 0) > 0.7*K))
  	{
  	  warning("Burn-in period finished, and more than 70% components have been shut down. Very high values of K could result in a poorly initialized model, therefore, it is recommended to re-run the model around K=", ceiling((sum(colSums(Z) != 0)+2)*1.25 ) )
  	}
  }

  }

  if(opts$checkConvergence & length(posterior$X)>=8) {

    if((sum(Z) == 0)){
      conv <- NA
      if(VERBOSE>0) print("No active components. Can not estimate convergence.")
      break;
    }

    #Estimate the convergence of the data reconstruction, based on the Geweke diagnostic
    ps <- floor(length(posterior$X)/4)
    start <- 1:ps
    end <- (-ps+1):0+length(posterior$X)
    Df <- N*D*L

    y <- matrix(NA,0,Df)
    for(ps in c(start,end)) {
      tmp <- rep(0,N*D*L)
      for(k in 1:K)
        tmp <- tmp + c(outer(outer(posterior$X[[ps]][,k],posterior$W[[ps]][[1]][,k]),posterior$U[[ps]][,k]))
      y <- rbind(y, tmp)
    }

    foo <- rep(NA,Df)
    for(j in 1:Df)
      foo[j] <- t.test(y[start,j],y[-start,j])$p.value
    conv <- mean(foo<=0.05)
    if(VERBOSE>0) {
      print(paste0("Convergence diagnostic: ",round(conv,4),". Values significantly greater than 0.05 imply a non-converged model."))
    }
  } else {
    conv <- NA
  }

  emptyK = which(colSums(Z) == 0)
  if(length(emptyK)>0 && length(emptyK)<K){
    if(nrow(Z)==1)
       Z = matrix(Z[,-emptyK],1,K-length(emptyK))
    else
       Z = Z[,-emptyK]
    U = U[,-emptyK]
    X = X[,-emptyK]
    for(m in 1:M)
      W[[m]] = W[[m]][,-emptyK]
    K = ncol(X)
  }
  if(K>1){
    sdv = apply(X,2,sd)*apply(U,2,sd)
    for(m in 1:M)
      sdv = sdv*apply(W[[m]],2,sd)
    sdv = order(sdv,decreasing=T)
    if(nrow(Z)==1)
      Z = matrix(Z[,sdv],1,K)
    else
      Z = Z[,sdv]
    U = U[,sdv]
    X = X[,sdv]
    for(m in 1:M)
      W[[m]] = W[[m]][,sdv]
  }
  if(VERBOSE>=0) {
    print(paste0("Run completed with ",K," active components."))
  }
  tau = tau[[1]][1,1] #for Bayesian CP
  return(list(X=X,W=W,U=U,K=K,tau=tau,posterior=posterior,cost=cost,opts=opts,conv=conv))
}

#' A function for generating a default set of parameters for Bayesian Tensor Factorization methods
#'
#' \code{getDefaultOpts} returns the default choices for model parameters.
#'
#' This function returns options for defining the model's high-level structure
#' (sparsity priors), the hyperparameters, and the uninformative
#' priors. We recommend keeping these as provided.
#'
#' @param method the factorization method for which options are required.
#' Currently only "CP" (default) is supported.
#' @return A list with the following model options:
#' \item{ARDX}{TRUE: use elementwise ARD prior for X, resulting in sparse X's.
#' FALSE: use guassian prior for a dense X (default).}
#' \item{ARDW}{TRUE: use elementwise ARD prior for W, resulting in sparse W's (default).
#' FALSE: use guassian prior for a dense W.}
#' \item{ARDU}{TRUE: use elementwise ARD prior for U, resulting in sparse U's.
#' FALSE: use guassian prior for a dense U (default).}
#' \item{iter.burnin}{The number of burn-in samples (default 5000).}
#' \item{iter.sampling}{The number of saved posterior samples (default 50).}
#' \item{iter.thinning}{The thinning factor to use in saving posterior samples (default 10).}
#' \item{prior.alpha_0t}{The shape parameter for residual noise (tau's) prior (default 1).}
#' \item{prior.beta_0t}{The rate parameter for residual noise (tau's) prior (default 1).}
#' \item{prior.alpha_0}{The shape parameter for the ARD precisions (default 1e-3).}
#' \item{prior.beta_0}{The rate parameter for the ARD precisions (default 1e-3).}
#' \item{prior.betaW1}{Bernoulli prior for component activiations, prior.betaW1 < prior.betaW2: sparsity inducing (default: 1).}
#' \item{prior.betaW2}{Bernoulli prior for component activation, (default: 1).}
#' \item{init.tau}{The initial value for noise precision (default 1e3).}
#' \item{verbose}{The verbosity level. 0=no printing, 1=moderate printing,
#'  2=maximal printing (default 1).}
#' \item{checkConvergence}{Check for the convergence of the data reconstruction,
#'  based on the Geweke diagnostic (default TRUE).}
#' @export
#' @examples
#' #To run the algorithm with other values:
#' opts <- getDefaultOpts()
#' opts$ARDW <- FALSE #Switch off Feature-level Sparsity on W's
#'  \dontrun{res <- tensorBF(Y=Y,opts=opts)}
getDefaultOpts <- function(method = "CP")
{
  if("CP" != method)
    print("Only CP factorization supported in this version.")

  return(list(prior.alpha_0=1e-3,prior.beta_0=1e-3,
              prior.alpha_0t=1,prior.beta_0t=1,
              verbose=1, init.tau=10^3,
              prior.betaW1=1,prior.betaW2=1,
              iter.burnin=5000,iter.sampling=50,iter.thinning=10,
	            checkConvergence=TRUE,
              ARDW = TRUE, ARDU = FALSE, ARDX = FALSE))
}
