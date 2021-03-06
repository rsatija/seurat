% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/seurat.R
\name{doKMeans}
\alias{doKMeans}
\title{K-Means Clustering}
\usage{
doKMeans(object, genes.use = NULL, k.genes = NULL, k.cells = NULL,
  k.seed = 1, do.plot = TRUE, data.cut = 2.5, k.cols = pyCols,
  pc.row.order = NULL, pc.col.order = NULL, rev.pc.order = FALSE,
  use.imputed = FALSE, set.ident = TRUE, ...)
}
\arguments{
\item{object}{Seurat object}

\item{genes.use}{Genes to use for clustering}

\item{k.genes}{K value to use for clustering genes}

\item{k.cells}{K value to use for clustering cells (default is NULL, cells
are not clustered)}

\item{k.seed}{Random seed}

\item{do.plot}{Draw heatmap of clustered genes/cells (default is TRUE)}

\item{data.cut}{Clip all z-scores to have an absolute value below this.
Reduces the effect of huge outliers in the data.}

\item{k.cols}{Color palette for heatmap}

\item{pc.row.order}{Order gene clusters based on the average PC score within
a cluster. Can be useful if you want to visualize clusters, for example,
based on their average score for PC1.}

\item{pc.col.order}{Order cell clusters based on the average PC score within
a cluster}

\item{rev.pc.order}{Use the reverse PC ordering for gene and cell clusters
(since the sign of a PC is arbitrary)}

\item{use.imputed}{Cluster imputed values (default is FALSE)}

\item{set.ident}{If clustering cells (so k.cells>0), set the cell identity
class to its K-means cluster (default is TRUE)}

\item{\dots}{Additional parameters passed to doHeatMap for plotting}
}
\value{
Seurat object where the k-means results for genes is stored in
object@kmeans.obj[[1]], and the k-means results for cells is stored in
object@kmeans.col[[1]]. The cluster for each cell is stored in object@data.info[,"kmeans.ident"]
and also object@ident (if set.ident=TRUE)
}
\description{
Perform k=means clustering on both genes and single cells
}
\details{
K-means and heatmap are calculated on object@scale.data
}

