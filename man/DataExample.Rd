\docType{data}
\name{DataExample}
\alias{DataExample}
\title{
Data example for SIDA
}

\format{A list with 4 elements
\describe{
  \item{Xdata}{ A list with each entry containing two views of training data with dimension \eqn{160\times 2000} each. Rows are samples and columns are variables.}
  \item{Y}{ \eqn{160 \times 1} vector of training class membership. There are three classes each with size 80.}
 \item{Xtestdata}{ A list with each entry containing two views of testing data with dimension \eqn{320\times 2000} each. Rows are samples and columns are variables.}
  \item{Ytest}{ \eqn{320 \times 1} vector of testing class membership. There are three classes each with size 160.}
}}

\usage{
data(DataExample)
}
\description{
Simulated data to demonstrate the use of SIDA.
}
\examples{

# see cvSIDA or sida examples

}
\references{
Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019) , \emph{Sparse Linear Discriminant Analysis for Multi-view Structured Data}, \emph{submitted}
}




