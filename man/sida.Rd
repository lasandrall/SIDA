\name{sida}
\alias{sida}
\title{
Sparse Integrative Discriminant Analysis for Multi-view Data
}
\description{
Performs sparse integrative disdcriminant analysis of multi-view data  to 1) obtain discriminant vectors that are associated and optimally separate subjects into different classes 2) estimate misclassification rate, and total correlation coefficient. Allows for the inclusion of other covariates which are not penalized in the algorithm. It is recommended to use cvSIDA to choose best tuning parameter.
}
\usage{
sida(Xdata=Xdata,Y=Y,Tau=Tau,withCov=FALSE,Xtestdata=Xtestdata,Ytest=Ytest,
    AssignClassMethod='Joint',plotIt=FALSE, standardize=TRUE,
    maxiteration=20,weight=0.5,thresh= 1e-03)
}
\arguments{
  \item{Xdata}{
A list with each entry containing training views of size \eqn{n \times p_d}, where \eqn{d =1,...,D} views. Rows are samples and columns are variables. If covariates are available, they  should be included as a separate view, and set as the last dataset. For binary  or categorical covariates (assumes no ordering), we suggest the use of indicator variables.
}
  \item{Y}{
\eqn{n \times 1} vector of class membership.
}
  \item{Tau}{
\eqn{d \times 1} vector of tuning parameter. It is recommended to use sidatunerange to obtain lower and upper bounds for the tuning parameters since too large a tuning parameter will result in a trivial solution vector (all zeros) and too small may result in non-sparse vectors.
}
  \item{withCov}{
TRUE or FALSE if covariates are available. If TRUE, please set all covariates as one dataset and should be the last dataset. For binary and categorical variables, use indicator matrices/vectors. Default is FALSE.
}
  \item{Xtestdata}{
A list with each entry containing testing views of size \eqn{ntest \times p_d}, where \eqn{d =1,...,D}. Rows are samples and columns are variables. The order of the list should be the same as the order for the training data, Xdata. Use if you want to predict on a testing dataset. If no Xtestdata, set to NULL.
}
   \item{Ytest}{
\eqn{ntest \times 1} vector of test class membership. If no testing data provided, set to NULL.
}
  \item{AssignClassMethod}{
  Classification method. Either Joint or Separate. Joint uses all discriminant vectors from D datasets to predict class membership. Separate predicts       class membership separately for each dataset. Default is Joint
}

\item{plotIt}{
 TRUE or FALSE.  If TRUE, produces discriminants and correlation plots. Default is FALSE
}
 \item{standardize}{
 TRUE or FALSE. If TRUE, data will be normalized to have mean zero  and variance one for each variable. Default is TRUE.
}

  \item{maxiteration}{
  Maximum iteration for the algorithm if not converged.Default is 20.
}
  \item{weight}{
Balances separation and association. Default is 0.5.
}
  \item{thresh }{
  Threshold for convergence. Default is 0.001.
}
}
\details{
The function will return several R objects, which can be assigned to a variable. To see the results, use the ``$" operator.
}
\value{
 \item{sidaerror}{ Estimated classication error. If testing data provided, this will be test classification error, otherwise, training error }
		\item{sidacorrelation}{ Sum of pairwise RV coefficients. Normalized to be within 0 and 1, inclusive.
	}
	\item{hatalpha}{ A list of estimated sparse discriminant vectors for each view. }
	\item{PredictedClass}{ Predicted class. If AssignClassMethod='Separate', this will be a \eqn{ntest\times D} matrix, with each column the predicted class for each data. }
%% ...
}
\references{
Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019), \emph{Sparse Linear Discriminant Analysis for Multi-view Structured Data}, \emph{submitted}%.\url{}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cvSIDA},\link{sidatunerange}, \link{CorrelationPlots},\link{DiscriminantPlots} }
}
\examples{
library(SIDA)
##---- read in data
data(DataExample)

Xdata=DataExample[[1]]
Y=DataExample[[2]]
Xtestdata=DataExample[[3]]
Ytest=DataExample[[4]]


##---- call sida algorithm to estimate discriminant vectors, and predict on testing data

#call sidatunerange to get range of tuning paramater
ngrid=10
mytunerange=sidatunerange(Xdata,Y,ngrid,standardize=TRUE,weight=0.5,withCov=FALSE)

# an example with Tau set as the lower bound
Tau=c(mytunerange$Tauvec[[1]][1], mytunerange$Tauvec[[2]][1])

mysida=sida(Xdata,Y,Tau,withCov=FALSE,Xtestdata=Xtestdata,Ytest=Ytest,
            AssignClassMethod='Joint',plotIt=TRUE, standardize=TRUE,
            maxiteration=20,weight=0.5,thresh= 1e-03)

test.error=mysida$sidaerror

test.correlation=mysida$sidacorrelation

hatalpha=mysida$hatalpha

predictedClass=mysida$PredictedClass


##----plot discriminant and correlation plots

#---------Discriminant plot
mydisplot=DiscriminantPlots(Xtestdata,Ytest,mysida$hatalpha)

mycorrplot=CorrelationPlots(Xtestdata,Ytest,mysida$hatalpha)

}
