\name{sidatunerange}
\alias{sidatunerange}
\title{
Tuning paramter grid values for sida
}
\description{
Sida function to provide tuning parameter grid values for each view, not including covariates, if available. It is recommended to use this to get lower and upper bounds of tuning parameters for each view that can be used in sida. This function is called by cvSIDA to select optimal tuning parameters.
}
\usage{
sidatunerange(Xdata=Xdata,Y=Y,ngrid=8,standardize=TRUE,weight=0.5,withCov=FALSE)
}
\arguments{
  \item{Xdata}{
A list with each entry containing each \eqn{n\times p_d} training view, where \eqn{d =1,...,D} views. Rows are samples and columns are variables. If covariates are available, they  should be included as a separate view, and set as the last dataset. For binary  or categorical covariates (assumes no ordering), we suggest the use of indicator variables.
}
  \item{Y}{
\eqn{n \times 1} vector of class membership. Same size as the number of training samples.
}
\item{ngrid}{
   Number of grid points for tuning parameters.
}
 \item{standardize}{
 TRUE or FALSE. If TRUE, data will be normalized to have mean zero  and variance one for each variable. Default is TRUE.
}
\item{weight}{
Balances separation and association. Default is 0.5.
}
 \item{withCov}{
TRUE or FALSE if covariates are available. If TRUE, set all covariates as one dataset and should be the last dataset. For binary and categorical variables, use indicator matrices/vectors. Default is FALSE.
}
}
\details{
The function will return an R object with grid values for each data, not including covariates, if available. To see the results, use the ``$" operator.
}
\value{
 	\item{Tauvec}{grid values for each data, not including covariates, if available. }
 	}

\references{
Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019) , \emph{Sparse Linear Discriminant Analysis for Multi-view Structured Data}, \emph{submitted}%.\url{}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sida} }
}
\examples{
## see examples in sida
}
