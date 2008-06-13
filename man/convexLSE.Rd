\name{convexLSE}
\alias{convexLSE}
\title{Compute the nonparamateric LSE of a convex hazard}
\description{
  This function computes the nonparametric LSE of a convex hazard function on \code{[0,TT]}.
}
\usage{
convexLSE(x, TT, M = 100, GRIDLESS = 0, tol = 1e-05, tol.SR = 1e-08, type = 1, 
max.loop = 100, max.loop.SR = 250, range = c(0, TT), solve.tol = .Machine$double.eps)
}
\arguments{
  \item{x}{vector containing the data.}
  \item{TT}{value specifying \code{TT} in the interval \code{[0,TT]} over which the LSE is found.}
  \item{M}{size of the grid used in the implementation of the (support reduction) algorithm.  The default value is \code{M=100}.}
  \item{GRIDLESS}{boolean variable.  If \code{GRIDLESS=1}, then a 'gridless' implementation of the algorithm is used.  The default is \code{GRIDLESS=0}.}
  \item{tol}{tolerance value of the bisection algorithm.  The default value is \code{tol=1e-05}.}
  \item{tol.SR}{tolerance value of the support reduction algorithm.  The default value is \code{tol.SR=1e-08}.}
  \item{type}{integer between 1 and 3 indicating the exit criterion function used by the bisection algorithm.  The default is \code{type=1}.}
  \item{max.loop}{maximum number of iterations in the bisection algorithm. The default is \code{max.loop=100}.  To be used for troubleshooting.}
  \item{max.loop.SR}{maximum number of iteration in the support reduction algorithm. The default is \code{max.loop.SR=250}.  To be used for troubleshooting.}
  \item{range}{vector of length two describing the range of antimodes over which the (bisection) algorithm searches.  The default is \code{range=[0,TT]}.  The range should always be contained within \code{[0,TT]}.  To be used in troubleshooting.}
  \item{solve.tol}{tolerance value used in \code{solve()} invoked by the algorithm.  The default is \code{solve.tol=.Machine$double.eps}, which is also the default of the \code{solve()} function.  To be used in troubleshooting. }
}
\details{
Let HHn be the empirical cumulative hazard function for the data x.  Then the LSE is the convex hazard function which is closest to the hazard function whose integral is given by HHn, in the least squares sense.  See Jankowski and Wellner (2007) for further details.

The function consists of an outer, bisection algorithm loop and an inner, support reduction algorithm loop.  The goal is to minimize the least squares criterion function over the space of positive convex hazards.  This is done in two steps:  the support reduction algorithm minimizes the least squares criterion function (\code{phi}) over the space of positive convex hazard with minimum (antimode) at \code{a} using the function \code{\link{srLSE}};  the bisection algorithm searches over the possible values of \code{a} in \code{[0,TT]}.   

The bisection algorithm ends once the exit criterion function is smaller than \code{tol}.  The exit criterion function may be specified via the option \code{type}:

\code{type=1}
\code{sum((cur.phi-min(old.phi))^2)}

\code{type=2} 		
\code{sqrt(sum((cur.phi - old.phi)^2))}

\code{type=3} 		
\code{sort(cur.phi)[2]+sort(cur.phi)[3]-2*sort(cur.phi)[1]}

%\tabular{rl}{
% \code{type=1} \tab \preformatted{sum((cur.phi-min(old.phi))^2)} \cr
% \code{type=2} \tab sqrt(sum((cur.phi - old.phi)^2)) \cr
% \code{type=3} \tab sort(cur.phi)[2]+sort(cur.phi)[3]-2*sort(cur.phi)[1]
%}
%\itemize{
% \item{\code{type=1}} {sum((cur.phi-min(old.phi))^2)}
% \item{\code{type=2}} {sqrt(sum((cur.phi - old.phi)^2))}
% \item{\code{type=3}} {sort(cur.phi)[2]+sort(cur.phi)[3]-2*sort(cur.phi)[1]}
%}
where \code{cur.phi} is the vector of the least squares criterion function \code{phi} along different values of antimode \code{a}, and \code{old.phi} is the same but in the previous iteration.

The support reduction algorithm was first described in Groeneboom, Jongbloed and Wellner (2007).  It is an iterative optimization procedure which exits when a sufficient tolerance level \code{tol.SR} is reached.   The algorithm uses a gridded implementation (see Jankowski and Wellner (2008) for the details) described by the parameter \code{M}.  The idea is that any convex function may be written as a mixture of elbow functions, and the grid describes the possible support of the mixing distribution.  (The documentation for the function \code{h} provides a little more information; see also Jankowski and Wellner (2007, 2008).)  The larger \code{M} is, the more accurate, but slower, the algorithm.  Setting \code{GRIDLESS=1} invokes an implementation which tries to overcome the grid, and increases the accuracy again.  Decreasing \code{tol} and/or \code{tol.SR} may also increase the accuracy of the algorithm, but not as much as manipulating the grid size.
}
\value{
A list which includes
  \item{lse}{
	list containing the support \code{supp} and mixing measure \code{mix} for the 
	computed LSE.  The function \code{\link{h}} can translate these into a 
	hazard function.
	}
  \item{ls}{
	value of the least squares criterion function at the LSE.  See Jankowski
	and Wellner (2007) for its exact form.
	}
  \item{antimode}{
	value of the antimode \code{a} at which the LSE of the hazard has its minimum.
        }
  \item{conv}{
	boolean indicating if the (bisection) algorithm has converged.
	}
  \item{iter}{
	number of iterations of the bisection algorithm.
	}
  \item{"convexLSE.out"}{
	The algorithm also creates a file called \code{convexLSE.out} in R's current working directory.  The file contains three variables:  (1) \code{phi} the value of the smallest least squares criterion function at the antimode found by the support reduction algorithm, (2) \code{antimode} the value of the antimode, and (3) \code{conv} a boolean indicating if the support reduction algorithm converged.  If \code{convexLSE.out} already exists, then the algorithm simply appends to it.
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
\section{Diagnostics}{
The file \code{"convexLSE.out"} contains information which may be used to determine how well the algorithm has worked.  Plotting \code{phi} against \code{antimode} should result in a u-shaped (first decreasing, then increasing) function, with minimum lying at the antimode of the LSE.  This is only true though for the values of \code{phi} (and \code{antimode}) where the \code{\link{srLSE}} algorithm converged (ie. \code{conv=TRUE}).  Plotting the output in \code{"convexLSE.out"} can be used to check this (see 'Examples'), and to check where convergence of \code{\link{srLSE}} failed.  
}
\section{Troubleshooting}{
If the bisection algorithm fails to converge, then \code{max.loop} may be increased, or changing the \code{range} may be used temporarily to diagnose the issue.  

Alternatively, the support reduction algorithm may fail to converge.  The file \code{"convexLSE.out"} shows which of the calls of the support reduction algorithm by the bisection loop converged.  If too many of these failed, then changing the values of \code{max.loop}, \code{max.loop.SR}, or \code{solve.tol} may fix this.

At each iteration, the support reduction algorithm does a finite dimensional quadratic optimization using the function \code{solve()}.  Occasionally (quite rarely for the LSE), the matrix passed to \code{solve()} is computationally singular.  Several built-in catches exist in the code to handle this, but if the problem occurs too frequently, the algorithm will abort.  In fact, this is by far the most common cause of nonconvergence.  The larger the value of \code{M}, the more likely the problem is, but this also depends on the data set.  Setting \code{solve.tol} to a smaller value may fix the problem.  

If the algorithm is interrupted, one should type \code{sink()} on the R command line, so that output will no longer be sent to the file \code{convexLSE.out}.
}
\seealso{\code{\link{convexMLE}} \code{\link{srLSE}} \code{\link{srMLE}} \code{\link{h}}}
\examples{
# Generate sample data:
set.seed(1111, kind="default")
x 	<- rweibull(50, 3) 

# Set the value of TT:
TT	<- 1

# Find the lse	
lse 	<- convexLSE(x, TT=TT)

# And the lse has minimum at
lse$antimode

# Check the bisection algorithm
diag	<- read.table("convexLSE.out", header=TRUE)
diag

plot(diag$antimode, diag$phi, pch=16, col="blue", cex=1, main="diagnostics", ylab="phi", xlab="antimode")


# plot the true hazard function and the lse
h.true 	<- function(x) 3*x^2
tt	<- c(0,sort(lse$lse$supp$tau), sort(lse$lse$supp$eta), TT)
yy	<- sapply(tt, h, supp=lse$lse$supp, hpar=lse$lse$mix)


plot(h.true, xlim=c(0,TT), lwd=2, col="blue", ylab="hazard", xlab="time", main="LSE of convex hazard")
lines(tt,yy, lwd=2, col="red")
legend("topleft", c("true","estimated"),col=c("blue", "red"), lwd=2, inset = .05)
}
\keyword{nonparametric}
\keyword{survival}
\keyword{iteration}
\keyword{optimize}
