% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onlyDeconAlgorithms.R
\name{estCellPercent.spillOver}
\alias{estCellPercent.spillOver}
\title{Estimate cell percentage from spillover}
\usage{
estCellPercent.spillOver(spillExpr, refExpr, geneExpr, method = "DCQ", ...)
}
\arguments{
\item{spillExpr}{A spill over matrix, as calculated by buildSpilloverMat(). (e.g. LM22.spillover.csv.gz)}

\item{refExpr}{a data frame representing immune cell expression profiles. Each row represents an expression of a gene, and each column represents a different immune cell type. colnames contains the name of each immune cell type and the rownames includes the genes' symbol. The names of each immune cell type and the symbol of each gene should be unique. Any gene with missing expression values must be excluded.}

\item{geneExpr}{a data frame representing RNA-seq or microarray gene-expression profiles of a given complex tissue. Each row represents an expression of a gene, and each column represents a different experimental sample. colnames contain the name of each sample and rownames includes the genes' symbol. The name of each individual sample and the symbol of each gene should be unique. Any gene with missing expression values should be excluded.}

\item{method}{One of 'DCQ', 'SVMDECON', 'DeconRNASeq', 'proportionsInAdmixture', 'nnls' (DEFAULT: DCQ)}

\item{...}{Parameters for estCellPercent.X (e.g. number_of_repeats for .DCQ)}
}
\value{
a matrix of estimate cell type percentages in samples
}
\description{
Use a spillover matrix to deconvolve a samples
}
\examples{
#This toy example 
library(ADAPTS)
fullLM22 <- ADAPTS::LM22[1:30, 1:4]
smallLM22 <- fullLM22[1:25,] 

spillover <- buildSpilloverMat(refExpr=smallLM22, geneExpr=fullLM22) 
cellEst <- estCellPercent.spillOver(spillExpr=spillover, refExpr=smallLM22, geneExpr=fullLM22)
}
