\name{pfnc}
\alias{pfnc}
\title{
Noncentral F Distribution
}
\description{
Using the first terms in the series expansion 26.6.20 in
Abramowitz and Stegun, compute the probability integral of
the noncentral F distribution.
}
\usage{
pfnc(q, df1, df2, lm, iprec=c(6))
}
\arguments{
\item{q}{
vector of quantiles.
}
\item{df1}{
vector of numerator degrees of freedom.
}
\item{df2}{
vector of denominator degrees of freedom.
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
probability that a noncentral F variable with numerator
degrees of freedom df1, denominator degrees of freedom df2 and
noncentrality parameter lm is less than q.
}
\references{
Milton Abramowitz and Irene A. Stegun, Handbook of Mathematical
Functions, p. 947, Dover, 1970.
}
\section{NOTE}{
Let n be the length of the longest of q, df1, df2 and lm.
If any of q, df1, df2 and lm is of length less than n, then all
values after the first in that vector are ignored.
Only the first element of iprec is used.
}
\seealso{
\code{\link{pf}}, \code{\link{pchisqnc}}
}
% Converted by Sd2Rd version 0.3-3.
