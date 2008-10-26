\name{convexHaz-package}
\alias{convexHaz-package}
\alias{convexHaz}
\docType{package}
\title{
Find the nonparametric MLE/LSE of convex hazard
}
\description{
This package contains functions to compute the nonparametric
maximum likelihood estimator (MLE) and the nonparametric least squares
estimator (LSE) of a convex hazard function, assuming that the data is IID.
}
\details{
\tabular{ll}{
Package: \tab convexHaz\cr
Version: \tab 0.1\cr
Date: \tab 2008-10-27\cr
Depends: \tab R (>= 2.4)\cr
License: \tab GPL (version 2 or later)\cr
Built: \tab R 2.7.1; ; 2008-09-02 13:05:46; windows\cr
}

Index:
\preformatted{
convexLSE               Compute the nonparamateric LSE of a convex
                        hazard
convexMLE               Compute the nonparametric MLE of a convex
                        hazard
hazard                  Create a hazard function (and its integral)
                        given its support and mixing measure.
srLSE                   Compute the LSE of a convex hazard with fixed
                        antimode
srMLE                   Compute the MLE of a convex hazard with fixed
                        antimode
}
}
\author{
Hanna Jankowski, Ivy Wang, Hugh McCague, Jon A. Wellner

Maintainer: Hanna Jankowski <hkj@mathstat.yorku.ca>
}
\references{
   Groeneboom, Jongbloed and Wellner (2008). The support reduction 
   algorithm for computing nonparametric function estimates in mixture
   models.  \emph{Scan. J. Statist.} \bold{35}, 385--399.

   Jankowski and Wellner (2007). Nonparametric estimation of a convex
   bathtub-shaped hazard function.  \emph{Technical Report} \bold{521}, Department of Statistics, University of Washington.

   Jankowski and Wellner (2008). Computation of nonparametric convex hazard estimators via profile methods.  \emph{Technical Report} \bold{542}, Department of Statistics, University of Washington.
}
\keyword{package}
\keyword{nonparametric}
\keyword{survival}
\keyword{iteration}
\keyword{optimize}
