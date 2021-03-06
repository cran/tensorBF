% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tensorBF.R
\name{predictTensorBF}
\alias{predictTensorBF}
\title{Predict Missing Values using the Bayesian tensor factorization model}
\usage{
predictTensorBF(Y, res)
}
\arguments{
\item{Y}{is a 3-mode tensor containing missing values as NA's. See function \code{\link{tensorBF}} for details.}

\item{res}{the model object returned by the function \code{\link{tensorBF}}.}
}
\value{
A tensor of the same size as \code{Y} containing predicted values in place of NA's.
}
\description{
\code{predictTensorBF} predicts the missing values in the data \code{Y} using the learned model \code{res}.
}
\details{
If the original data \code{Y} contained missing values (NA's),
this function predicts them using the model. The predictions are
returned in the un-normalized space if \code{res$pre} contains appropriate
preprocessing information.
}
\examples{
#Data generation
\dontrun{K <- 2}
\dontrun{X <- matrix(rnorm(20*K),20,K)}
\dontrun{W <- matrix(rnorm(30*K),30,K)}
\dontrun{U <- matrix(rnorm(3*K),3,K)}
\dontrun{Y = 0}
\dontrun{for(k in 1:K) Y <- Y + outer(outer(X[,k],W[,k]),U[,k])}
\dontrun{ Y <- Y + array(rnorm(20*30*3,0,0.25),dim=c(20,30,3))}

#insert missing values
\dontrun{m.inds = sample(prod(dim(Y)),100)}
\dontrun{Yobs = Y[m.inds]}
\dontrun{Y[m.inds] = NA}

#Run the method with default options and predict missing values
\dontrun{res <- tensorBF(Y)}
\dontrun{pred = predictTensorBF(Y=Y,res=res)}
\dontrun{plot(Yobs,pred[m.inds],xlab="obs",ylab="pred",main=round(cor(Yobs,pred[m.inds]),2))}
}
