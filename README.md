
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Bio3DView <img src="man/figures/bio3d-logo.png" style="width: 50%;" align="right" />

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

## Example 1

First let’s load up the packages and generate a quick NGL (webGL based)
structure overview of a bio3d pdb class object with a number of simple
defaults. The returned NGLVieweR object can be further added to for
custom interactive visualizations:

``` r
library(bio3dview)
library(bio3d)
library(NGLVieweR)

pdb <- read.pdb("5p21")
view.pdb(pdb) |>
  setSpin()
```

<figure>
<img src="man/figures/fig1a.gif"
alt="Figure 1. Structure of HRas PDB code: 5p21. Note that the image here is not interactive due to restrictions with GitHub GFM format." />
<figcaption aria-hidden="true"><strong>Figure 1</strong>. Structure of
HRas PDB code: 5p21. <em>Note that the image here is not interactive due
to restrictions with GitHub GFM format.</em></figcaption>
</figure>

## Example 2.

Here we generate a quick interactive multi-structure ensemble view of a
bio3d `pdbs` object:

``` r
data(transducin)

view.pdbs(transducin$pdbs, colorScheme = "res") 
```

<figure>
<img src="man/figures/fig2.png"
alt="Figure 2. All 53 PDB structures of Transducin colored by residue index" />
<figcaption aria-hidden="true"><strong>Figure 2.</strong> All 53 PDB
structures of Transducin colored by residue index</figcaption>
</figure>

## Example 3.

As a final example let’s perform a quick Normal Mode Analysis (NMA) and
view the predicted large scale motions:

``` r
adk <- read.pdb("6s36")
m <- nma(adk)
view.nma(m, pdb=adk) |>
  setPlay()
```

<figure>
<img src="man/figures/fig3.gif"
alt="Figure 3. Predicted large scale domain motions of Adenalate kinase from a bio3d based Normal Mode Analysis." />
<figcaption aria-hidden="true"><strong>Figure 3.</strong> Predicted
large scale domain motions of Adenalate kinase from a bio3d based Normal
Mode Analysis.</figcaption>
</figure>

## Going further

There are many, many more options for customized viewing options
available. Further examples can be found on our [getting started
article]() with full details documented within in the individual
functions [help pages]().
