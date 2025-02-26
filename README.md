
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bio3dview

<!-- badges: start -->
<!-- badges: end -->

Interactive biomolecular structure visualization of
[bio3d](http://thegrantlab.org/bio3d/) objects in R.

## Installation

You can install the development version of bio3dview from
[GitHub](https://github.com/bioboot/bio3dview) with:

``` r
# install.packages("pak")
pak::pak("bioboot/bio3dview")
```

Dependencies include the R CRAN packages
[bio3d](https://cran.r-project.org/web/packages/bio3d/index.html) and
[NGLVieweR](https://cran.r-project.org/web/packages/NGLVieweR/), which
can be installed from CRAN with:

``` r
install.packages("bio3d")
install.packages("NGLVieweR")
```

## Example

Generate a quick NGL (webGL based) structure overview of bio3d pdb class
objects with a number of simple defaults. The returned NGLVieweR object
can be further added to for custom interactive visualizations.

``` r
library(bio3dview)
library(bio3d)
library(NGLVieweR)

pdb <- read.pdb("1hsg")
#>   Note: Accessing on-line PDB file
#view.pdb(pdb)
```

``` r
#htmlwidgets::saveWidget(v, "man/figures/temp.html")
#webshot2::webshot("man/figures/temp.html", "man/figures/README-screenshot.png", delay = 0.5)
```

``` r
#view.pdb(pdb, ligand=FALSE, cols=c("pink","aquamarine"))
```

``` r
#view.pdb(pdb, colorScheme = "sstruc", backgroundColor = "black")
```

``` r
sele <- atom.select(pdb, resno=c(25, 50))
#view.pdb(pdb, highlight = sele,
#        cols = c("navy","orange"),
#        backgroundColor = "pink",
#        highlight.style = "spacefill")
```

Quick interactive multi-structure ensemble bio3d pdbs object viewing
using the NGLVieweR package.

``` r
#pth <- "~/Desktop/courses/BIMM143/class10/pdbs/split_chain/"
#files <- list.files(path=pth, full.names = TRUE)
#pdbs <- pdbaln(files, fit=TRUE, exefile="msa")

#view.pdbs(pdbs, representation = "cartoon")
```

``` r
#view.pdbs(pdbs, colorScheme = "residueindex")
```

And a PCA analysis result

``` r
#pc <- pca(pdbs)

# Plot a single conformer plot of PC1 v PC2
#plot(pc, pc.axes = 1:2)
```

``` r
## Plot atom wise loadings
#plot.bio3d(pc$au[,1], ylab="PC1 (A)")
```

and a trajectory along PC1

``` r
#view.pca(pc)
```

Major functions of the package include:
