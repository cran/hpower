\name{sscomp}
\alias{sscomp}
\title{
Sample Size for Tests of the General Linear Hypothesis
}
\description{
Compute the minimum sample size required to have a specified power
for rejecting the null hypothesis in a test of the general linear
hypothesis.
}
\usage{
sscomp(alpha, power, x, bt1, sigma, cc, u, th0, test=1, ninit=1)
}
\arguments{
\item{alpha}{
type I error probability of the test.
}
\item{power}{
desired power of the test.
}
\item{x}{
N0-by-q design matrix.
This matrix should represent the design for a sample of size 1.
For example, if there are to be two groups, x should contain two
rows, with the first row representing the design for group 1 and
the second row the design for group 2.
The sample size returned will then be the number per group.
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
test to use.  test=1 (the default) is the Wilks lambda test; test=2
is the Hotelling-Lawley trace test; test=3 is the Pillai-Bartlett
trace test.  See Muller and Peterson (1984) for details.
}
\item{ninit}{
number of replications of x to use as the starting sample size.
The sample size algorithm is simple -- it starts at ninit and keeps
adding or subtracting 1 (that is, 1 copy of x) until the desired
power is achieved.
The default is ninit=1.
}}
\value{
list containing two elements -- n, the smallest number of copies of
x that guarantee the desired power for the test described, and
power, the power of the design with that sample size.
}

\author{Daniel F. Heitjan <dheitjan@peter.cpmc.columbia.edu>,  
R-port by Stefan Funke <funke@attglobal.net>}

\references{
Muller, K.E. and Peterson, B.L. (1984).  Practical methods for
computing power in testing the multivariate general linear
hypothesis.  Computational Statistics and Data Analysis 2, 143--158.
}
\section{NOTE}{
The power function glhpwr is an approximation, so there may be
slight differences when compared with more exact methods.
This function increases (or decreases) n one unit at a time, so
well chosen starting values can save time.
Some experimentation will often be helpful.
}
\seealso{
\code{\link{pfnc}}, \code{\link{glhpwr}}.
}
% Converted by Sd2Rd version 0.3-3.
