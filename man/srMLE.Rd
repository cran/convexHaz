\name{srMLE}
\alias{srMLE}
\title{Compute the MLE of a convex hazard with fixed antimode}
\description{
	This function computes the MLE of a convex hazard function with antimode \code{a}.  
}
\usage{
srMLE(x, a, M = 100, GRIDLESS = 0, ini.LSE = 0, tol = 1e-08, max.inner.loop = 250, 
max.outer.loop = 50, print = 0, solve.tol = .Machine$double.eps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector containing the data.}
  \item{a}{value of the antimode (location of the minimum of the hazard rate). }
  \item{M}{size of the grid used in the implementation of the algorithm.  The default value is \code{M=100}.}
  \item{GRIDLESS}{boolean variable.  \code{If GRIDLESS=1}, then a 'gridless' implementation of the algorithm is used.  The default value is \code{GRIDLESS=0}.}
  \item{ini.LSE}{boolean variable.  If \code{ini.LSE=1} then the initial value for the support reduction algorithm for the MLE is taken from the least squares estimator.  The default value is \code{ini.LSE=0}.}
  \item{tol}{tolerance value of the algorithm.  The default value is \code{tol=1e-08}.}
  \item{max.inner.loop}{maximum number of iterations in the inner loop of the support reduction algorithm.  The default value is \code{max.inner.loop=250}.  To be used for troubleshooting.}
  \item{max.outer.loop}{maximum number of iterations in the outer loop of the support reduction algorithm.  The default value is \code{max.outer.loop=50}.  To be used for troubleshooting.}
  \item{print}{boolean variable.  If print=1, then the results of each iteration are printed.  The default value is \code{print=0}.  To be used for troubleshooting.}
  \item{solve.tol}{tolerance value used in \code{solve()} invoked by the algorithm.  The default value is the same as that of the \code{solve()} function. To be used for troubleshooting.}
}
\details{
This function implements the support reduction algorithm to find the nonparametric MLE of a convex hazard function with antimode at \code{a} under IID sampling.  The likelihood in this case is actually a modified likelihood, which excludes the term h(max(x)). The full likelihood can be made arbitrarily large by increasing h(max(x)), and therefore this term is omitted. The estimated MLE returned by the algorithm is thus the function on [0, max(x)), and we assume that the MLE is equal to infinity beyond max(x). 

Increasing \code{M} and/or decreasing \code{tol} will increase the accuracy of the algorithm, but increase the computation time.  Setting \code{GRIDLESS=1} will also increase the accuracy, but this may cause the algorithm not to converge (see 'Troubleshooting').  
}
\value{
A list which includes
  \item{mle}{
	list containing the support \code{supp} and mixture \code{mix} (further lists) for the 
	computed MLE.  The function \code{\link{h}} can translate these into a 
	hazard function.
	}
  \item{llh}{
	value of the (modified) log likelihood at the MLE.
	}
  \item{conv}{
	boolean indicating if the algorithm has converged.
	}
  \item{iter}{
	number of iterations (of the outer loop) of the algorithm.
	}
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
\section{Troubleshooting}{
If the algorithm fails to converge, then there are several options available.  Setting \code{print=1} in the options may reveal to cause of the problem.  

The support reduction algorithm consists of inner and outer loops, if the algorithm reaches the maximum in either one of these, then increasing \code{max.outer.loop} or \code{max.inner.loop} may fix the problem, though this is unlikely.

At each iteration, the support reduction algorithm preforms a finite dimensional quadratic optimization using the function \code{solve()}.  Occasionally, the matrix passed to \code{solve()} is computationally singular.  Several built-in catches exist in the code to handle this, but if the problem occurs too frequently, the algorithm will abort.  In fact, this is by far the most common cause of nonconvergence.  The larger the value of \code{M}, the more likely the problem is, but this also depends on the data set.  The problem is also likely to occur if \code{GRIDLESS=1}.  Setting \code{solve.tol} to a smaller value may fix the problem.  

In addition, changing the starting point of the algorithm may help, which is done by setting \code{sr.LSE=1}.  This may also increase the speed of the algorithm. 
}
\seealso{\code{\link{convexLSE}} \code{\link{convexMLE}} \code{\link{hazard}} \code{\link{srLSE}}} 
\examples{
# Generate sample data:
set.seed(1111, kind="default")
x 	<- rweibull(50, 3) 

# Find the LSE with antimode a=0 over the range [0,1]:
TT	<- 1	
mle 	<- srMLE(x, a=0)

# plot the true hazard function and the lse
h.true 	<- function(x) 3*x^2
tt	<- c(0,sort(mle$mle$supp$tau), sort(mle$mle$supp$eta), max(x))
yy	<- sapply(tt, h, supp=mle$mle$supp, hpar=mle$mle$mix)


plot(h.true, xlim=c(0,max(x)), lwd=2, col="blue", ylab="hazard", xlab="time", main="MLE (a=0)")
lines(tt,yy, lwd=2, col="red")
legend("topleft", c("true","estimated"),col=c("blue", "red"), lwd=2, inset = .05)
}
\keyword{nonparametric}
\keyword{survival}
\keyword{iteration}
\keyword{optimize}
