\docType{data}
\name{SIDANetDataExample}
\alias{SIDANetDataExample}
\title{Data example for SIDANet}
\format{A list with 6 elements:
\describe{
  \item{XdataNet}{ A list with each entry containing two views of training data with dimension \eqn{240\times 1000} each. Rows are samples and columns are variables.}
  \item{YNet}{ \eqn{240 \times 1} vector of training class membership. There are three classes each with size 80.}
 \item{XtestdataNet}{ A list with each entry containing two views of testing data with dimension \eqn{480\times 1000} each. Rows are samples and columns are variables.}
  \item{YtestNet}{ \eqn{480 \times 1} vector of testing class membership. There are three classes each with size 160.}
 \item{myedges}{ A list with each entry containing a \eqn{36\times 2} matrix of edge information for each view. Assumes variable 1 is connected to variables 2 to 10, variable 11 is connected to variables 12 to 20, variable 21 is connected to variables 22 to 30 and variable 31 is connected to variables 32 to 40. All remaining variables are singletons.}
\item{myedgeweight}{ A list with each entry containing edgeweight. In this example, views 1 and 2 have edge weights so the Laplacian of a weighted graph will be used.}
}}

\usage{
data(SIDANetDataExample)
}
\description{
Simulated data to demonstrate the use of SIDANet.
}
\examples{

# see cvSIDANet or sidanet examples

}
\references{
Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019) , \emph{Sparse Linear Discriminant Analysis for Multi-view Structured Data}, \emph{submitted}
}




