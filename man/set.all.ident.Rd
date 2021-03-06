% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/seurat.R
\name{set.all.ident}
\alias{set.all.ident}
\title{Switch identity class definition to another variable}
\usage{
set.all.ident(object, id = NULL)
}
\arguments{
\item{object}{Seurat object}

\item{id}{Variable to switch identity class to (for example, 'DBclust.ident', the output
of density clustering) Default is orig.ident - the original annotation pulled from the cell name.}
}
\value{
A Seurat object where object@ident has been appropriately modified
}
\description{
Switch identity class definition to another variable
}

