\name{convexMLE}
\alias{convexMLE}
\title{Compute the nonparametric MLE of a convex hazard}
\description{
This function computes the nonparametric MLE of a convex hazard function.
}
\usage{
convexMLE(x, M = 100, GRIDLESS = 0, ini.LSE = 0, tol = 1e-05, tol.SR = 1e-08, type = 1, 
max.loop = 100, max.inner.loop.SR = 250, max.outer.loop.SR = 50, range = c(0, max(x)), 
solve.tol = .Machine$double.eps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector containing the data.}
  \item{M}{size of the grid used in the implementation of the (support reduction) algorithm.  The default value is \code{M=100}.}
  \item{GRIDLESS}{boolean variable.  If \code{GRIDLESS=1}, then a 'gridless' implementation of the algorithm is used.  The default value is \code{GRIDLESS=0}.}
  \item{ini.LSE}{boolean variable.  If \code{ini.LSE=1} then the initial value for each call of the support reduction algorithm for the MLE is taken from the least squares estimator.  The default value is \code{ini.LSE=0}.}
  \item{tol}{tolerance value of the bisection algorithm.  The default value is \code{tol=1e-05}.}
  \item{tol.SR}{tolerance value of the support reduction algorithm.  The default value is \code{tol.SR=1e-08}.}
  \item{type}{integer between 1 and 3 indicating the exit criterion function used by the bisection algorithm.  The default is \code{type=1}.}
  \item{max.loop}{maximum number of iterations in the bisection algorithm.  The default value is \code{max.loop=100}.  To be used for troubleshooting.}
  \item{max.inner.loop.SR}{maximum number of iterations in the inner loop of the support reduction algorithm.  The default value is \code{max.inner.loop.SR=250}.  To be used for troubleshooting.}
  \item{max.outer.loop.SR}{maximum number of iterations in the outer loop of the support reduction algorithm.  The default value is \code{max.outer.loop.SR=50}. To be used for troubleshooting.}
  \item{range}{vector of length two describing the range of antimodes over which the (bisection) algorithm searches.  The default is \code{range=c(0,max(x))}.  The range should always be contained in \code{[0,max(x)]}.  To be used in troubleshooting.}
  \item{solve.tol}{tolerance value used in \code{solve()} invoked by the algorithm.  The default is \code{solve.tol=.Machine$double.eps}, which is also the default of the \code{solve()} function.  To be used in troubleshooting. }
}
\details{
This function gives an iterative procedure to find the nonparametric MLE of a convex hazard function under IID sampling.  The likelihood in this case is actually a modified likelihood, which excludes the term h(\code{max(x)}).  The full likelihood can be made arbitrarily large by increasing h(\code{max(x)}), and therefore this term is omitted.  The estimated MLE returned by the algorithm is thus the function on [0, \code{max(x)}), with the assumption that the MLE is equal to infinity beyond \code{max(x)}.

The function consists of an outer, bisection algorithm loop and an inner, support reduction algorithm loop which calls upon \code{\link{srMLE}}.  The goal is to maximize the likelihood over the space of positive convex hazards.  This is done in two steps:  the support reduction algorithm maximizes the log likelihood (or minimizes its negative, \code{negllh}) over the space of positive convex hazard with minimum (antimode) at \code{a};  the bisection algorithm searches over the possible values of \code{a} in \code{[0,max(x)]}.   

The bisection algorithm ends once the exit criterion function is smaller than \code{tol}.  The exit criterion function may be specified via the option \code{type}:

\code{type=1}
\code{sum((cur.phi-min(old.phi))^2)}

\code{type=2} 		
\code{sqrt(sum((cur.phi - old.phi)^2))}

\code{type=3} 		
\code{sort(cur.phi)[2]+sort(cur.phi)[3]-2*sort(cur.phi)[1]}

where \code{cur.phi} is the vector of the least squares criterion function \code{phi} along different values of antimode \code{a}, and \code{old.phi} is the same but in the previous iteration.

The support reduction algorithm was first described in Groeneboom, Jongbloed and Wellner (2007).  It is an iterative optimization procedure which exits when a sufficient tolerance level \code{tol.SR} is reached.   The algorithm uses a gridded implementation (see Jankowski and Wellner (2008) for the details) described by the parameter \code{M}.  The idea is that any convex function may be written as a mixture of elbow functions, and the grid describes the possible support of the mixing distribution.  (The documentation for the function \code{h} provides a little more information; see also Jankowski and Wellner (2007, 2008).)  The larger \code{M} is the more accurate, and slower, the algorithm.  Setting \code{GRIDLESS=1} will also increase the accuracy, but this often results in non-convergence of the support reduction algorithm (see 'Troubleshooting').    
}
\value{
  \item{mle}{
	list containing the support \code{supp} and mixing measure \code{mix} (further lists) for the 
	computed MLE.  The function \code{\link{h}} can translate these into a 
	hazard function.
	}
  \item{llh}{
	value of the (modified) log likelihood at the MLE.
	}
  \item{antimode}{
	value of the antimode \code{a} at which the MLE of the hazard has its minimum.
        }
  \item{conv}{
	boolean indicating if the (bisection) algorithm has converged.
	}
  \item{iter}{
	number of iterations of the bisection algorithm.
	}
  \item{"convexMLE.out"}{
	The algorithm also creates a file called \code{convexMLE.out} in R's current working directory.  The file contains three variables:  (1) \code{negllh} the value of the smallest least squares criterion function at the antimode found by the support reduction algorithm, (2) \code{antimode} the value of the antimode, and (3) \code{conv} a boolean indicating if the support reduction algorithm converged.  If \code{convexMLE.out} already exists, then the algorithm will overwrite it.  
	}
}
\references{
   Groeneboom, Jongbloed and Wellner (2008). The support reduction 
   algorithm for computing nonparametric function estimates in mixture
   models.  \emph{Scan. J. Statist.} \bold{35}, 385--399.

   Jankowski and Wellner (2007). Nonparametric estimation of a convex
   bathtub-shaped hazard function.  \emph{Technical Report} \bold{521}, Department of Statistics,
   University of Washington.

   Jankowski and Wellner (2008). Computation of nonparametric convex hazard estimators via profile methods.  \emph{Technical Report} \bold{542}, Department of Statistics, University of Washington.
}
\author{Hanna Jankowski: \email{hkj@mathstat.yorku.ca}, Ivy Wang, Hugh McCague, Jon A. Wellner}
\section{Diagnostics}{
The file \code{"convexMLE.out"} contains information which may be used to determine how well the algorithm has worked.  Plotting \code{negllh} against \code{antimode} should result in a U-shaped (first decreasing, then increasing) function, with minimum lying at the antimode of the MLE.  This is only true though for the values of \code{negllh} (and \code{antimode}) where the \code{\link{srMLE}} algorithm converged (ie. \code{conv=TRUE}).  Plotting the output in \code{"convexMLE.out"} can be used to check this (see 'Examples'), and to check where convergence of \code{\link{srMLE}} failed.  
}
\section{Troubleshooting}{
If the bisection algorithm fails to converge, then \code{max.loop} may be increased.  Also, changing the \code{range} may used temporarily to diagnose the issue.

Alternatively, the support reduction algorithm may fail to converge.  The file \cr \code{"convexMLE.out"} shows which of the calls of the support reduction algorithm by the bisection loop converged.  If too many of these failed, then changing the values of \cr \code{max.inner.loop.SR}, \code{max.outer.loop.SR}, \code{solve.tol} or \code{ini.LSE} may fix this.  Increasing \code{max.inner.loop.SR} or \code{max.outer.loop.SR} gives the algorithms more time to converge, though it is rare that this is the cause of the problem.  Setting \code{ini.LSE} changes the starting point of the support reduction algorithm, and may even increase its speed.

At each iteration, the support reduction algorithm does a finite dimensional quadratic optimization using the function \code{solve()}.  Occasionally, the matrix passed to \code{solve()} is computationally singular.  Several built-in catches exist in the code to handle this, but if the problem occurs too frequently, the algorithm will abort.  In fact, this is by far the most common cause of nonconvergence.  The larger the value of \code{M}, the more likely the problem is, but this also depends on the data set.  Also, choosing \code{GRIDLESS=1} highly increases the probability that this type of problem will occur.   Setting \code{solve.tol} to a smaller value may fix the problem.  
}
\seealso{\code{\link{convexLSE}} \code{\link{hazard}} \code{\link{srLSE}} \code{\link{srMLE}}}
\examples{
# Generate sample data:
set.seed(1111, kind="default")
x 	<- rbeta(100, 0.5, 0.5) 

# Find the mle (this takes a few minutes)	
mle 	<- convexMLE(x, M=1000)

# And the mle has minimum at
mle$antimode

# Check the bisection algorithm
diag	<- read.table("convexMLE.out", header=TRUE)
diag

plot(diag$antimode, diag$negllh, pch=16, col="blue", main="diagnostics", ylab="phi", xlab="antimode")


# plot the true hazard function and the lse
h.true 	<- function(x) dbeta(x, 0.5, 0.5)/(1-pbeta(x, 0.5, 0.5))
tt	<- c(0,sort(mle$mle$supp$tau), sort(mle$mle$supp$eta), max(x))
yy	<- sapply(tt, h, supp=mle$mle$supp, hpar=mle$mle$mix)


plot(h.true, xlim=c(0, 0.95), lwd=2, col="blue", ylab="hazard", xlab="time", main="MLE of convex hazard")
lines(tt,yy, lwd=2, col="red")
legend("topleft", c("true","estimated"),col=c("blue", "red"), lwd=2, inset = .05)

%# plot the cumulative hazard functions
%H.true	<- function(x) x^3

%tt.HH	<- sort(c(0,x,x))
%ytemp	<- sapply(sort(x), HHn, xx=sort(x))
%yy.HH	<- sort(c(0,0, ytemp[1:(length(ytemp)-1)], ytemp))
%yy.H	<- sapply(tt, H, supp=mle$mle$supp, hpar=mle$mle$mix)

%plot(H.true, xlim=c(0,max(x)), lwd=2, col="blue", ylab="cumulative hazard", xlab="time", main="comparing cumulative hazards")
%lines(tt.HH, yy.HH, col="black", lwd=2)
%lines(tt, yy.H, col="red", lwd=2)
%legend("topleft", c("true","estimated", "empirical"),col=c("blue", "red", "black"), lwd=2, inset = .05)
}
\keyword{nonparametric}
\keyword{survival}
\keyword{iteration}
\keyword{optimize}
