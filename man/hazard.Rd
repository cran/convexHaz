\name{hazard}
\alias{hazard}
\alias{H}
\alias{h}
\title{Create a hazard function (and its integral) given its support and mixing measure.}
\description{
	Create a hazard function (and its integral) given its support and mixing measure.
}
\usage{
h(tt, supp, hpar)

H(tt, supp, hpar)
}
\arguments{
  \item{tt}{real positive number.}
  \item{supp}{list containing the values in the support: e.g. \code{list(constant = 0, tau = c(0.5), eta = c(1,2))}.}
  \item{hpar}{list containing the values in the mixing measure: e.g. \code{list(alpha = numeric(), nu = c(3), mu = c(2,0.5))}.}
}
\details{
The functions \code{\link{convexMLE}} and \code{\link{convexLSE}}, as well as \code{\link{srMLE}} and \code{\link{srLSE}} return the MLE/LSE in terms of the support and mixture of the hazard function.  The function \code{h} translates these into a hazard function, and the function \code{H} gives the cumulative hazard function (integral of \code{h}).

	Any convex hazard function can be written as a mixture of three types of basis functions: a constant function \code{e(t)=1}, decreasing elbow fuctions \code{e(t)=max(tau-t,0)}, and increasing elbow functions \code{e(t)=max(t-eta,0)}.  If the minimum of the hazard is located at a point \code{a}, then the possible values of tau lie in the set \code{[0,a]} and the possible values of eta must be greater than \code{a}.  We call the point \code{a} the antimode.  \code{supp} describes the list of basis functions to use: if constant is 0 then there is no constant function, and \code{tau} and \code{eta} list the values of tau and eta used, respectively.  \code{hpar} describes the mixing weights associated to these functions, with \code{alpha} giving the weight of the constant, and \code{nu} and \code{mu} giving the weights of tau and eta, respectively.  Note that the weights need not total to one, and that for a discrete mixture the hazard rate will be a piecewise constant function.

For example, if \code{supp=list(constant=1, tau=c(0.5), eta=numeric())} and \cr
\code{hpar=list(constant=3, nu=c(2), mu=numeric())}, then the hazard function is\cr  \code{h(t)=3*1+2*max(0.5-t,0)}.
}
\value{
value of the hazard or cumulative hazard function at (time) t.
}
\references{
  Jankowski and Wellner (2007). Nonparametric estimation of a convex
   bathtub-shaped hazard function.  \emph{Technical Report} \bold{521}, Department of Statistics, University of Washington.

   Jankowski and Wellner (2008). Computation of nonparametric convex hazard estimators via profile methods.  \emph{Technical Report} \bold{542}, Department of Statistics, University of Washington.
}
\author{Hanna Jankowski: \email{hkj@mathstat.yorku.ca}}
\seealso{\code{\link{convexLSE}} \code{\link{convexMLE}} \code{\link{srLSE}} \code{\link{srMLE}}}
\examples{
# Generate sample data:
set.seed(1111, kind="default")
x 	<- rweibull(50, 3) 

# Find the LSE with antimode a=0 over the range [0,1]:
TT	<- 1
lse 	<- srLSE(x, a=0, TT=TT)

# create simpler function to evaluate hazard and cumulative hazard
h.lse 	<- function(t){ return(h(t, lse$lse$supp, lse$lse$mix))}
H.lse 	<- function(t){ return(H(t, lse$lse$supp, lse$lse$mix))}


# hazard function at t=0
h.lse(0)

# cumulative hazard function at t=1
H.lse(1)

# plot the hazard function h.lse
tt	<- c(0,sort(lse$lse$supp$tau), sort(lse$lse$supp$eta), TT)  	# where h.lse changes slope
yy	<- sapply(tt, h.lse)						# values of h.lse at tt
plot(tt, yy, xlim=c(0,TT), type="l", lwd=2, col="red", ylab="hazard", xlab="time", main="LSE (a=0)")
}
\keyword{survival}
\keyword{nonparametric}
