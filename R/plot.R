#' Plot Tensor Components
#'
#' \code{plotTensorBF} shows the heatmap of components inferred by \code{\link{tensorBF}}.
#'
#' @param res The learned tensorBF model.
#' @param Y The original input data to be plotted. If specified NULL,
#' the function plots the data reconstruction using \code{\link{reconstructTensorBF}} (default: NULL).
#' @param k the component number to visualize (default: 1).
#' @param modesOnAxis which mode to plot on each axis c(Yaxis,Xaxis,lateral). Defaults to c(1,2,3).
#' @param nTopFeatures The number of most relevant features to show for the data space
#'        visualizations in each of the modes. Defaults to c(5,15,3) for displaying top 10 features
#'        of \eqn{1^{st}} mode, 20 of \eqn{2^{nd}} mode and 5 of \eqn{3^{rd}} mode.
#' @param margins numeric vector of length 4 containing the margins (see par(mar= *))
#' @param cex.axis positive numbers, used as cex.axis (default: 1)
#' @param cols colors used for the image. Defaults to a blue-white-red color scale.
#' @param key logical indicating whether a color-key should be drawn.
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
#' #Run the method with default options
#' \dontrun{res1 <- tensorBF(Y)}
#' \dontrun{plotTensorBF(res = res1,Y=Y,k=1)}
plotTensorBF <- function(res, Y=NULL, k=1, modesOnAxis=c(1,2,3), nTopFeatures=c(5,15,3), margins=c(4,4,4,12), cex.axis=1, cols=colorRampPalette(c("blue","white","red"))(101), key=TRUE)
{
  if(is.null(res)) stop("Please specify a correct model object.")
  if(is.null(k)) stop("Please specify a correct value of component k to plot.")
  if(!is.numeric(k)) stop("k should be a positive numeric number.")
  if(k<1 || k>ncol(res$X))
    stop(paste0("K should be between: ",1," and ",ncol(res$X)))
  if(!is.logical(key)) stop("key should be a logical value of TRUE/FALSE.")
  if(!is.numeric(cex.axis) || cex.axis < 0.01) {
    warning("cex.axis should be a positive numeric value."); cex.axis = 1;
  }

  if(length(dim(Y)) != length(modesOnAxis))
    stop(paste0("Number of modes specified is ",length(modesOnAxis),". Should be ",length(dim(Y)),"."))
  if(length(dim(Y)) != length(nTopFeatures))
    stop(paste0("Number of nTop feature count is ",length(nTopFeatures),". Should be ",length(dim(Y)),"."))
  if(4 != length(margins))
    stop(paste0("Number of margins count is ",length(margins),". Should be 4."))

  if(any( !(modesOnAxis %in% 1:length(dim(Y))) ))
    stop(paste0("Wrong value of mode specified in modesOnAxis. Should be unique values in the range of 1 to ",length(dim(Y))))

  if(!is.numeric(nTopFeatures) || any(is.na(nTopFeatures)) || any(nTopFeatures<1))
    stop("Wrong value of top features specified in nTopFeatures Should be positive numbers.")

  inds = which(nTopFeatures > dim(Y));
  if(length(inds)>0) {
    nTopFeatures[inds] = dim(Y)[inds]
    warning("Too high values in nTopFeatures. Plotting only the maximum number of available dimensions.")
  }
  if(any(nTopFeatures>50))
    print("Plot for more than 50 features could be cluttered. Try reducing the values in nTopFeatures if the image is not clear.")

  side="A"
  if(is.null(Y)) {
    Y = reconstructTensorBF(res)
    print("Y not supplied. Plotting the reconstruction.")
  }
  # adjusting for modes on different axis
  md = list()
  md$X = res$X
  md$U = res$U
  md$W = res$W#[[1]]
  if(2==modesOnAxis[1]) md$X = res$W#[[1]]
  if(3==modesOnAxis[1]) md$X = res$U
  if(1==modesOnAxis[2]) md$W = res$X
  if(3==modesOnAxis[2]) md$W = res$U
  if(1==modesOnAxis[3]) md$U = res$X
  if(2==modesOnAxis[3]) md$U = res$W#[[1]]
  Y = aperm(Y,modesOnAxis)
  nTopFeatures = nTopFeatures[modesOnAxis]
  mnames = c("X","W","U")[modesOnAxis]

  md <- normX.max1(md)
  md <- normU.max1(md)

	if(side=="A")
		x.inds <- order(md$X[,k],decreasing=TRUE)[1:nTopFeatures[1]]
	else
		x.inds <- order(md$X[,k],decreasing=FALSE)[1:nTopFeatures[1]]

	x.inds <- rev(x.inds) #reveresed for ploting purpose.
	w.inds.1 <- order(abs(md$W[,k]),decreasing=TRUE)[1:nTopFeatures[2]]

	dd = Y[x.inds,w.inds.1,1]
	if(nTopFeatures[3]>1)
	  for(n3 in 2:nTopFeatures[3])
	    dd = rbind(dd,Y[x.inds,w.inds.1,n3])
	dd <- t(dd)

	w.inds.1 <- w.inds.1[hclust(dist(dd))$order] #reorder for clustering of genes
	cex.a = cex.axis*1.3
	Y <- plotlimitdata.tensor(Y,limit=3)
	o.op <- par(mfcol=c(nTopFeatures[3]+key,1), oma=c(1,0,0,0.3),mar=c(0,margins[2],margins[3],margins[4])+0.2)
	breaks <- 1:(length(cols) + 1)
	extreme <- max(abs(Y), na.rm = TRUE)
	breaks <- seq(-extreme, extreme, length = length(breaks))
	cell.label.pos = c(-8.3,-21.2,-29.4)
	m3Names = dimnames(Y)[[3]]
	for(o in 1:nTopFeatures[3])
	{
  	if(o < nTopFeatures[3] && o > 1) { par(mar=c(margins[1]/2,margins[2],margins[3]/2,margins[4])+0.2); }
  	if(o == nTopFeatures[3]) par(mar=c(margins[1],margins[2],0,margins[4])+0.2)
  	dat <- Y[x.inds,w.inds.1,o]
  	image(x=t(dat),col=cols,axes=FALSE,breaks=breaks)
  	axis(2,at=0.5,labels=paste(m3Names[o],"\n",mnames[3],": ",format(round(md$U[o,k],2),width=4,nsmall=2),sep=""),las=3,cex.axis=cex.a,tick=0)

  	if(o == 1){
  		axis(3,at=seq(0,1,l=ncol(dat)),labels=colnames(dat),las=2,cex.axis=cex.a)
  	  axis(1,at=seq(0,1,l=ncol(dat)),labels=rep("",ncol(dat)),las=1,cex.axis=cex.a)
  	}
  	if(o < nTopFeatures[3] && o > 1){
  		axis(3,at=seq(0,1,l=ncol(dat)),labels=rep("",ncol(dat)),las=2,cex.axis=cex.a)
  	  axis(1,at=seq(0,1,l=ncol(dat)),labels=rep("",ncol(dat)),las=1,cex.axis=cex.a)
  	}
  	axis(4,at=seq(0,1,l=nrow(dat)),labels=paste(mnames[1],": ",format(round(md$X[x.inds,k],2),width=4)," ",rownames(dat),sep=""),las=1,cex.axis=cex.a)
	}
	axis(3,at=seq(0,1,l=ncol(dat)),labels=rep("",ncol(dat)),las=1,cex.axis=cex.a)
	axis(1,at=seq(0,1,l=ncol(dat)),labels=paste(mnames[2],":",format(round(md$W[w.inds.1,k],2),width=4),sep=""),las=2,cex.axis=cex.a)
	if(key) plotKey(cols,breaks)
	par(o.op);
}#EndFunction

plotKey <- function(col,v)
{
	op <- par(mar = par("mar")+0.5+c(-3,3,3,3))
	z=matrix(1:100, ncol=1)
	image((z), col=col, xaxt="n", yaxt="n")
	title(main="Key",outer=FALSE,line=0.5)
	axis(1,at=c(0,0.5,1),labels=round(c(min(v),mean(v),max(v)),1),las=1)
	par(op)
}

plotlimitdata.tensor <- function(Y,limit=2)
{
	inds = which(Y > limit); if(length(inds)>0) Y[inds] <- limit
	inds = which(Y < (-limit)); if(length(inds)>0) Y[inds] <- (-limit)
	return(Y)
}

normU.max1 <- function(model)
{
	ll <- normK.max1(model$U,model$W)
	model$U <- ll$U
	model$W <- ll$W
	return(model)
}

normX.max1 <- function(model)
{
	ll <- normK.max1(model$X,model$W)
	model$X <- ll$U
	model$W <- ll$W
	return(model)
}

normK.max1 <- function(U,W)
{
	L <- nrow(U)
	D <- nrow(W)
	if(L == 1)
	{
		Ru <- as.vector(U)
		U <- U*outer(rep(1,L),1/Ru)
		W <- W*outer(rep(1,D),Ru)
	}
	else
	{
		Ru.x <- apply(U,2,max)
		Ru.n <- apply(U,2,min)
		Ru <- Ru.x
		inds <- which(abs(Ru.n)  > Ru.x)
		if(length(inds)>0) Ru[inds] <- Ru.n[inds]
		U <- U*outer(rep(1,L),1/Ru)
		W <- W*outer(rep(1,D),Ru)
	}
	return(list(U=U,W=W))
}
