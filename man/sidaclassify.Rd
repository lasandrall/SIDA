\name{sidaclassify}
\alias{sidaclassify}
\title{
Classification approach for Sparse Integrative Discriminant Analysis
}
\description{
Performs classification using nearest centroid on separate or combined estimated discriminant vectors, and predicts class membership.
}
\usage{
sidaclassify(hatalpha=hatalpha,Xtestdata=Xtestdata,Xdata=Xdata,Y=Y,
            AssignClassMethod='Joint', standardize=TRUE)
}
\arguments{
	\item{hatalpha}{ A list of estimated sparse discriminant vectors for each view. This may be obtained from sida or cvSIDA }
	 \item{Xtestdata}{
A list with each entry containing testing views of size \eqn{ntest \times p_d}, where \eqn{d =1,...,D} views. Rows are samples and columns are variables. The order of the list should be the same as the order for the training data, Xdata. If covariates are available, they  should be included as a separate view, and set as the last dataset. For binary  or categorical covariates (assumes no ordering), we suggest the use of indicator variables.
If you want to obtain training error, set as Xdata.}
  \item{Xdata}{
A list with each entry containing training views of size \eqn{n \times p_d}, where \eqn{d =1,...,D} views. Rows are samples and columns are variables. If covariates are available, they  should be included as a separate view, and set as the last dataset. For binary  or categorical covariates (assumes no ordering), we suggest the use of indicator variables.
}
  \item{Y}{
\eqn{n \times 1} vector of class membership. Same size as the number of training samples.
}
  \item{AssignClassMethod}{
  Classification method. Either Joint or Separate. Joint uses all discriminant vectors from D datasets to predict class membership. Separate predicts       class membership separately for each dataset. Default is Joint
}
 \item{standardize}{
 TRUE or FALSE. If TRUE, data will be normalized to have mean zero  and variance one for each variable. Default is TRUE.
}
}
\details{
The function will return an R object, showing the predicted class and the classification method. To see the results, use the ``$" operator.
}
\value{
 	\item{PredictedClass}{ Predicted class. If AssignClassMethod='Separate', this will be a \eqn{ntest\times D} matrix, with each column the predicted class for each data. }
\item{AssignClassMethod}{ Classification method used. Joint or Separate. }
}
\references{
Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019) , \emph{Sparse Linear Discriminant Analysis for Multi-view Structured Data}, \emph{submitted}%.\url{}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cvSIDA},\link{sida},\link{cvSIDANet},\link{sidanet} }
}
\examples{
##---- read in data
data(DataExample)

Xdata=DataExample[[1]]
Y=DataExample[[2]]
Xtestdata=DataExample[[3]]
Ytest=DataExample[[4]]


#call sidatunerange to get range of tuning paramater
ngrid=10
mytunerange=sidatunerange(Xdata,Y,ngrid,standardize=TRUE,weight=0.5,withCov=FALSE)

# an example with Tau set as the lower bound
Tau=c(mytunerange$Tauvec[[1]][1], mytunerange$Tauvec[[2]][1])

mysida=sida(Xdata,Y,Tau,withCov=FALSE,Xtestdata=Xtestdata,Ytest=Ytest)


#classification with combined estimated vectors
mysida.classify.Joint=sidaclassify(mysida$hatalpha,Xtestdata,Xdata,Y,
                        AssignClassMethod='Joint')

mysida.PredClass.Joint=mysida.classify.Joint$PredictedClass

#classification with separate estimated vectors
mysida.classify.Separate=sidaclassify(mysida$hatalpha,Xtestdata,Xdata,Y,
                        AssignClassMethod='Separate')

mysida.PredClass.Separate=mysida.classify.Separate$PredictedClass

}
