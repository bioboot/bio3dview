% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/view_ngl.R
\name{pdbs2string}
\alias{pdbs2string}
\title{Convert PDBS to Multi-Element Character Strings}
\usage{
pdbs2string(pdbs, collapse = TRUE)
}
\arguments{
\item{pdbs}{a multi-structure \code{pdbs} object as obtained from \code{pdbaln()},
\code{read.fasta.pdb()}, etc.}

\item{collapse}{logical, if TRUE a single element vector is returned. If
FALSE a multi-element vector with a element per \code{pdb} structure is
returned. The later is required for setting distinct viewing options per
structure - such as user defined colors (e.g. one per structure) etc.}
}
\value{
a character vector of structure data.
}
\description{
Convert a \strong{bio3d} multi-structure \code{pdbs} object to a multi-element
character vector that can be used as input for \strong{NGLVieweR} and r3dmol.
}
\examples{
  pth <- "~/Desktop/courses/BIMM143/class10/pdbs/split_chain/"
  files <- list.files(path=pth, full.names = TRUE)
  pdbs <- bio3d::pdbaln(files, fit=TRUE, exefile="msa")

 NGLVieweR::NGLVieweR( bio3d::pdbs2pdb(pdbs), format="pdb") |>
    NGLVieweR::addRepresentation("cartoon")

 # Or more simpley...
 view.pdbs(pdbs)
 # Trace, tube, line, cartoon, ball+stick
 view.pdbs(pdbs, representation = "trace")
 view.pdbs(pdbs, cols = c("red","blue") )
 view.pdbs(pdbs, colorScheme = "residueindex")

}
\seealso{
\code{view.pdbs()}, \code{pdb2string()}, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}, \code{bio3d::read.fasta.pdb()}, \code{bio3d::pdbaln()}.
}
\author{
Barry Grant, \email{bjgrant@ucsd.edu}
}
