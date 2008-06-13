\name{srLSE}
\alias{srLSE}
\title{Compute the LSE of a convex hazard with fixed antimode}
\description{
	This function computes the LSE of a convex hazard function on 
	\code{[0,TT]} with antimode \code{a}.  
}
\usage{
srLSE(x, a, TT, M = 100, GRIDLESS = 0, tol = 1e-08, max.loop = 250, print = 0, 
solve.tol = .Machine$double.eps) 
}
\arguments{
  \item{x}{vector containing the data.}
  \item{a}{value of the antimode (the location of the minimum of the hazard rate).}
  \item{TT}{value specifiying \code{TT} in the interval \code{[0,TT]} over which the LSE if found.}
  \item{M}{size of the grid used in the implementation of the algorithm.  The default value is \code{M=100}.}
  \item{GRIDLESS}{boolean variable.  If \code{GRIDLESS=1}, then a 'gridless' implementation of the algorithm is used. The default value is \code{GRIDLESS=0}.}
  \item{tol}{tolerance value of the algorithm.  The default is \code{tol=1e-08}.}
  \item{max.loop}{maximum number of iterations in the algorithm.  The default is \cr 
\code{max.loop=250}.  To be used for troubleshooting.}
  \item{print}{boolean variable.  If \code{print=1}, then the results of each iteration are printed.  The default value is \code{print=0}.  To be used for troubleshooting.}
  \item{solve.tol}{tolerance value used in \code{solve()} invoked by the algorithm.  To be used for troubleshooting.  The default value is the same as for the \code{solve()} function.}
}
\details{
Let HHn be the empirical cumulative hazard function for the data x.  Then the LSE is the convex hazard function, with minimum at \code{a}, which is closest to the hazard function whose integral is given by HHn, in the least squares sense.  See Jankowski and Wellner (2007) for further details.

The function implements the support reduction algorithm Groeneboom, Jongbloed and Wellner (2007).

Increasing \code{M}, decreasing \code{tol}, and setting \code{GRIDLESS=1} will each increase the accuracy of the algorithm.  However, these will also increase the computing time, especially increasing \code{M}.  
}
\value{
A list which includes
  \item{lse}{
	list containing the support \code{supp} and mixing measure \code{mix} (further lists) for the 
	computed LSE.  The function \code{\link{h}} can translate these into a 
	hazard function.
	}
  \item{ls}{
	value of the least squares criterion function at the LSE.  
	}
  \item{iter}{
	number of iterations of the algorithm.
	}
  \item{conv}{
	boolean indicating if the algorithm has converged.
	}
}
\section{Troubleshooting}{
If the algorithm fails to converge, then there are several options available.  Setting \code{print=1} can reveal the reason behind the failure.  If the maximal number of iterations was reached, then increasing \code{max.loop} may solve the problem.  If a \code{linear solve problem} occured, then one can try to decrease the tolerance of the \code{solve()} function, \code{solve.tol}.
}
\references{
   Groeneboom, Jongbloed and Wellner (2008). The support reduction 
   algorithm for computing nonparametric function estimates in mixture
   models.  \emph{Scan. J. Statist.} \bold{35}, 385--399.

   Jankowski and Wellner (2007). Nonparametric estimation of a convex
   bathtub-shaped hazard function.  \emph{Technical Report} \bold{521}, Department of Statistics, University of Washington.

   Jankowski and Wellner (2008). Computation of nonparametric convex hazard estimators via profile methods.  \emph{Technical Report} \bold{542}, Department of Statistics, University of Washington.
}
\author{Hanna Jankowski: \email{hkj@mathstat.yorku.ca}, Ivy Wang, Hugh McCague, Jon A. Wellner}
\seealso{\code{\link{convexMLE}} \code{\link{convexLSE}} \code{\link{hazard}} \code{\link{srMLE}}}
\examples{
# Generate sample data:
set.seed(1111, kind="default")
x 	<- rweibull(50, 3) 

# Find the LSE with antimode a=0 over the range [0,1]:
TT	<- 1	
lse 	<- srLSE(x, a=0, TT=TT)

# plot the true hazard function and the lse
h.true 	<- function(x) 3*x^2
tt	<- c(0,sort(lse$lse$supp$tau), sort(lse$lse$supp$eta), TT)
yy	<- sapply(tt, h, supp=lse$lse$supp, hpar=lse$lse$mix)


plot(h.true, xlim=c(0,TT), lwd=2, col="blue", ylab="hazard", xlab="time", main="LSE (a=0)")
lines(tt,yy, lwd=2, col="red")
legend("topleft", c("true","estimated"),col=c("blue", "red"), lwd=2, inset = .05)
}
\keyword{nonparametric}
\keyword{survival}
\keyword{iteration}
\keyword{optimize}
