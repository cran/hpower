***
\name{pchisqnc}
\alias{pchisqnc}
\title{
Noncentral Chi-square Distribution
}
\description{
Using the first terms in the series expansion 26.4.25 in
Abramowitz and Stegun, compute the probability integral of the
noncentral chi-square distribution.
}
\usage{
pchisqnc(q, df, lm, iprec=c(6))
}
\arguments{
\item{q}{
vector of quantiles.
}
\item{df}{
vector of degrees of freedom.
}
\item{lm}{
vector of noncentrality parameters.
}
\item{iprec}{
a parameter governing the precision of the answer; higher values
lead to greater precision.
The default is iprec=6.
For large lm, the true answer will exceed the returned value by
no more than 1-pnorm(iprec).
}}
\value{
probability that a noncentral chi-square variable with degrees of
freedom df and noncentrality parameter lm is less than q.
}

\author{Daniel F. Heitjan <dheitjan@peter.cpmc.columbia.edu>,  
R-port by Stefan Funke <funke@attglobal.net>}

\references{
Milton Abramowitz and Irene A. Stegun, Handbook of Mathematical
Functions, p. 942, Dover, 1970.
}
\section{NOTE}{
Let n be the length of the longest of q, df and lm.
If any of q, df and lm is of length less than n, then all
values after the first in that vector are ignored.
Only the first element of iprec is used.
}
\seealso{
\code{\link{pfnc}} , \code{\link{pchisq}}
}
% Converted by Sd2Rd version 0.3-3.
