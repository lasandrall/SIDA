\name{sidanettunerange}
\alias{sidanettunerange}
\title{
Tuning paramter grid values for sidanet
}
\description{
Sidanet function to provide tuning parameter grid values for each view, not including covariates, if available. It is recommended to use this to get lower and upper bounds of tuning parameters for each view that can be used in sidanet. This function is called by cvSIDANet to select optimal tuning parameters.
}
\usage{
sidanettunerange(Xdata,Y,ngrid,standardize=TRUE,weight=0.5,eta=0.5,
                myedges,myedgeweight,withCov=FALSE)
}
\arguments{
  \item{Xdata}{
A list with each entry containing training views of size \eqn{n \times p_d}, where \eqn{d =1,...,D} views. Rows are samples and columns are variables. If covariates are available, they  should be included as a separate view, and set as the last dataset. For binary  or categorical covariates (assumes no ordering), we suggest the use of indicator variables.
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
\item{eta}{
Balances the selection of network, and variables within network. Default is 0.5.
}
\item{myedges}{
A list with each entry containing a \eqn{M_d\times 2} matrix of edge information for each view. If a view has no edge information, set to 0; this will default to SIDA. If covariates are available as a view (\eqn{D}th view), the edge information should be set to 0.
}
\item{myedgeweight}{
A list with each entry containing a \eqn{M_d\times 1} vector of weight information for each view. If a view has no weight information, set to 0; this will use the Laplacian of an unweighted graph. If covariates are available as a view (\eqn{D}th view), the weight information should be set to 0.
}
 \item{withCov}{
TRUE or FALSE if covariates are available. If TRUE,  set all covariates as one dataset and should be the last dataset. For binary and categorical variables, use indicator matrices/vectors. Default is FALSE.
}
}
\details{
The function will return an R object with grid values for each data, not including covariates, if available. To see the results, use the ``$" operator.
}
\value{
 	\item{Tauvec}{Grid values for each data, not including covariates, if available. }
 	}

\references{
Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019) , \emph{Sparse Linear Discriminant Analysis for Multi-view Structured Data}, \emph{submitted}%.\url{}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sidanet} }
}
\examples{
# see examples in sidanet

}
