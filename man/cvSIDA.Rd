\name{cvSIDA}
\alias{cvSIDA}
\title{
Cross validation for Sparse Integrative Discriminant Analysis for Multi-view Data
}
\description{
Peforms nfolds cross validation  to select optimal tuning parameters for sida based on training data, which are then used with the training or testing  data
 to predict class membership. Allows for inclusion of covariates which are not penalized. If you want to apply optimal tuning parameters to testing data, you may also use sida.
}
\usage{
cvSIDA(Xdata,Y,withCov=FALSE,plotIt=FALSE, Xtestdata=NULL,Ytest=NULL,
        isParallel=TRUE, ncores=NULL,gridMethod='RandomSearch',
        AssignClassMethod='Joint',nfolds=5,ngrid=8,standardize=TRUE,
        maxiteration=20,weight=0.5,thresh=1e-03)
}
\arguments{
  \item{Xdata}{
A list with each entry containing training views of size \eqn{n \times p_d}, where \eqn{d =1,...,D}. Rows are samples and columns are variables. If covariates are available, they  should be included as a separate view, and set as the last dataset. For binary  or categorical covariates (assumes no ordering), we suggest the use of indicator variables.
}
  \item{Y}{
\eqn{n \times 1} vector of class membership.
}
  \item{withCov}{
TRUE or FALSE if covariates are available. If TRUE, please set all covariates as one dataset and should be the last dataset. For binary and categorical variables, use indicator matrices/vectors. Default is FALSE.
}
\item{plotIt}{
 TRUE or FALSE.  If TRUE, produces discriminants and correlation plots. Default is FALSE.
}
  \item{Xtestdata}{
A list with each entry containing testing views of size \eqn{ntest \times p_d}, where \eqn{d =1,...,D}. Rows are samples and columns are variables. The order of the list should be the same as the order for the training data, Xdata. Use if you want to predict on a testing dataset. If no Xtestdata, set to NULL.
}
   \item{Ytest}{
\eqn{ntest \times 1} vector of test class membership. If no testing data provided, set to NULL.
}
 \item{isParallel}{
TRUE or FALSE for parallel computing. Default is TRUE.
}
\item{ncores}{
Number of cores to be used for parallel computing. Only used if isParallel=TRUE. If isParallel=TRUE and ncores=NULL, defaults to half the size of the number of system cores.
}
 \item{gridMethod}{
GridSearch or RandomSearch. Optimize tuning parameters over full grid or random grid. Default is RandomSearch.
}

  \item{AssignClassMethod}{
  Classification method. Either Joint or Separate. Joint uses all discriminant vectors from D datasets to predict class membership. Separate predicts       class membership separately for each dataset. Default is Joint
}
\item{nfolds}{
  Number of cross validation folds. Default is 5.
}
\item{ngrid}{
   Number of grid points for tuning parameters. Default is 10 for each view if \eqn{D=2}. If \eqn{D>2}, default is 8.
}
 \item{standardize}{
 TRUE or FALSE. If TRUE, data will be normalized to have mean zero  and variance one for each variable. Default is TRUE.
}
  \item{maxiteration}{
  Maximum iteration for the algorithm if not converged. Default is 20.
}
  \item{weight}{
Balances separation and association. Default is 0.5.
}
  \item{thresh }{
  threshold for convergence. Default is 0.001.
}
}
\details{
The function will return several R objects, which can be assigned to a variable. To see the results, use the ``$" operator.
}
\value{
 \item{sidaerror}{ estimated classication error. If testing data provided, this will be test classification error, otherwise, training error }
		\item{sidacorrelation}{ sum of pairwise RV coefficients. Normalized to be within 0 and 1, inclusive.
	}
	\item{hatalpha}{ A list of estimated sparse discriminant vectors for each view. }
	\item{PredictedClass}{ Predicted class. If AssignClassMethod='Separate', this will be a \eqn{ntest\times D} matrix, with each column the predicted class for each data. }
	\item{optTau}{ Optimal tuning parameters for each view, not including covariates, if available. }
\item{gridValues}{ Grid values used for searching optimal tuning paramters. }
\item{AssignClassMethod}{ Classification method used. Joint or Separate. }
\item{gridMethod}{ Grid method used. Either GridSearch or RandomSearch}

}
\references{
Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019) , \emph{Sparse Linear Discriminant Analysis for Multi-view Structured Data}, \emph{submitted}%.\url{}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sida},\link{CorrelationPlots},\link{DiscriminantPlots} }
}
\examples{
library(SIDA)
##---- read in sample data
data(DataExample)

Xdata=DataExample[[1]]
Y=DataExample[[2]]
Xtestdata=DataExample[[3]]
Ytest=DataExample[[4]]


##---- call cross validation

mycv=cvSIDA(Xdata,Y,withCov=FALSE,plotIt=FALSE, Xtestdata=Xtestdata,Ytest=Ytest,
            isParallel=TRUE,ncores=NULL,gridMethod='RandomSearch',
            AssignClassMethod='Joint',nfolds=5,ngrid=8,standardize=TRUE,
            maxiteration=20, weight=0.5,thresh=1e-03)

#check output
test.error=mycv$sidaerror

test.correlation=mycv$sidacorrelation

optTau=mycv$optTau

hatalpha=mycv$hatalpha

#---------Discriminant plot
mydisplot=DiscriminantPlots(Xtestdata,Ytest,mycv$hatalpha)

mycorrplot=CorrelationPlots(Xtestdata,Ytest,mycv$hatalpha)
}
