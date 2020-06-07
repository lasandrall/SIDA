\name{DiscriminantPlots}
\alias{DiscriminantPlots}
\title{
Discriminant Plots
}
\description{
Plots discriminant vectors for visualizing class separation
}
\usage{
DiscriminantPlots(Xtestdata=Xtestdata,Ytest=Ytest,hatalpha=hatalpha)
}

\arguments{
  \item{Xtestdata}{
A list with each entry containing views of size \eqn{ntest \times p_d}, where \eqn{d =1,...,D}. Rows are samples and columns are variables. Can use testing or training data.
}
   \item{Ytest}{
\eqn{ntest \times 1} vector of class membership.
}
  \item{hatalpha}{
  A list of estimated sparse discriminant vectors for each view.
}
}
\details{
The function will return discriminant plots.
}
\value{
 \item{NULL}{ }
	 }
%% ...

\references{
Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019) , \emph{Sparse Linear Discriminant Analysis for Multi-view Structured Data}, \emph{submitted}%.\url{}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cvSIDA},\link{sidatunerange}, \link{CorrelationPlots} }
}
\examples{
library(SIDA)
##---- read in data
data(DataExample)

##---- call sida algorithm to estimate discriminant vectors, and predict on testing data

Xdata=DataExample[[1]]
Y=DataExample[[2]]
Xtestdata=DataExample[[3]]
Ytest=DataExample[[4]]


#call sidatunerange to get range of tuning paramater
ngrid=10
mytunerange=sidatunerange(Xdata,Y,ngrid,standardize=TRUE,weight=0.5,withCov=FALSE)

# an example with Tau set as the lower bound
Tau=c(mytunerange$Tauvec[[1]][1], mytunerange$Tauvec[[2]][1])

mysida=sida(Xdata,Y,Tau,withCov=FALSE,Xtestdata=Xtestdata,Ytest=Ytest,
            AssignClassMethod='Joint', plotIt=TRUE, standardize=TRUE,
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

