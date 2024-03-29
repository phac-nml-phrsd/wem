% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{broken_line}
\alias{broken_line}
\title{Piecewise linear function}
\usage{
broken_line(x, b, v)
}
\description{
For a given \code{x}, returns the value on a piecewise 
linear function as defined by 
vector \code{b} of its break points and associated 
values \code{v}.
It is typically used for relative factor.for time-dependent 
parameters. For example: if n = 2, broken line returns below v values:
t < t1: returns v1
t > t2: returns v2
t1 < t < t2: returns the line from (t1,v1) to (t2,v2)
}
