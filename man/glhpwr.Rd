\name{glhpwr}
\alias{glhpwr}
\title{
Power Function for Tests of the General Linear Hypothesis
}
\description{
Use a noncentral F approximation to compute the power of a test of
the general linear hypothesis.
}
\usage{
glhpwr(alpha, x, bt1, sigma, cc, u, th0, test=1, tol=1e-08)
}
\arguments{
\item{alpha}{
type I error probability of the test.
}
\item{x}{
N-by-q design matrix.
}
\item{bt1}{
q-by-p alternative hypothesis value of the regression coefficient.
}
\item{sigma}{
p-by-p variance matrix of a single observation.
}
\item{cc}{
a-by-q between-rows contrast matrix.
}
\item{u}{
p-by-b between-columns contrast matrix.
}
\item{th0}{
a-by-b null value of cc\%*\%beta\%*\%u.
}
\item{test}{
test for which power is desired.  test=1 (the default) is the Wilks
lambda test; test=2 is the Hotelling-Lawley trace test; test=3 is
the Pillai-Bartlett trace test.  See Muller and Peterson (1984) for
details.
}
\item{tol}{
tolerance for computing the rank of x (using qr()).  The default is
1.e-8.
}}
\value{
the approximate power of a test of the general linear hypothesis
cc\%*\%beta\%*\%u=th0 under the alternative beta=bt1.
The model is y(Nxp)=x(Nxq)\%*\%beta(qxp)+e(Nxp), where x is the design,
beta is the matrix of multivariate regression parameters, and e is
the error matrix, whose rows are assumed to be independent draws
from a multivariate normal with mean 0 and p-by-p variance matrix
sigma.
}

\author{Daniel F. Heitjan <dheitjan@peter.cpmc.columbia.edu>,  
R-port by Stefan Funke <funke@attglobal.net>}

\references{
Muller, K.E. and Peterson, B.L. (1984).  Practical methods for
computing power in testing the multivariate general linear
hypothesis.  Computational Statistics and Data Analysis 2, 143--158.
}
\section{NOTE}{
If the computed approximate degrees of freedom are negative, a
power of 0 is returned; this may mean that the proposed design
has no degrees of freedom for error.
}
\seealso{
\code{\link{pfnc}}, \code{\link{sscomp}}
}
% Converted by Sd2Rd version 0.3-3.
