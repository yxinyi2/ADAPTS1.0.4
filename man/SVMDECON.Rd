% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onlyDeconAlgorithms.R
\name{SVMDECON}
\alias{SVMDECON}
\title{Support vector machine deconvolution}
\usage{
SVMDECON(m, B)
}
\arguments{
\item{m}{a matrix represenging the mixture (genes X 1 sample)}

\item{B}{a matrix representing the references (genes X cells), m should be subset to match B}
}
\value{
A matrix with cell type estimates for each samples
}
\description{
Use SVMDECONV to estimate the cell count percentage
David L Gibbs, dgibbs@systemsbiology.org
June 9, 2017

v-SVR is applied with a linear kernel to solve for f,
and the best result from three values of v = {0.25, 0.5, 0.75}
is saved, where ‘best’ is defined as the lowest root mean squared error
between m and the deconvolution result, f x B.

Our current implementation executes v-SVR using the
‘svm’ function in the R package, ‘e1071’.

w2 <- SVMDECON(m, B)
}
