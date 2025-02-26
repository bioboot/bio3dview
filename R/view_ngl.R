# Functions for NGL based viewing of bio3d objects
# 2025-02-14   (11:44:38 PST on Fri, Feb 14)
#
# ToDo: Improve view nma/pca  function that makes use of pca2string

#' Convert a PDB Object to a Character String
#'
#' Function to take a **bio3d** `pdb` structure object and convert it to a single
#' element character vector that can then be used as input for the **NGLVieweR**
#' and **r3dmol** packages.
#'
#' @param pdb a bio3d pdb object as obtained from `read.pdb()`.
#' @param ... Extra arguments that are passed to `write.pdb()`.
#'
#' @return a single element character vector with return characters as required by `NGLVieweR::NGLVieweR()` function with `format=pdb` option. Note that no intermediate files are pro
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#'
#' @seealso \code{view.pdb()} which uses this function internally, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}.
#'
#' @export
#'
#' @examples
#' pdb <- bio3d::read.pdb("5p21")
#'
#' NGLVieweR::NGLVieweR(pdb2string(pdb), format="pdb") |>
#'   NGLVieweR::addRepresentation("cartoon")
#'
#' # Or more simply
#' view.pdb(pdb)
#'
pdb2string <- function(pdb, ...) {
  return( paste(utils::capture.output(
    bio3d::write.pdb(pdb, file="", ...) ),
      collapse = "\n") )
}



#' Convert PDBS to Multi-Element Character Strings
#'
#' Convert a **bio3d** multi-structure `pdbs` object to a multi-element
#' character vector that can be used as input for **NGLVieweR** and r3dmol.
#'
#' @param pdbs a multi-structure `pdbs` object as obtained from `pdbaln()`,
#' `read.fasta.pdb()`, etc.
#' @param collapse logical, if TRUE a single element vector is returned. If
#'   FALSE a multi-element vector with a element per `pdb` structure is
#'   returned. The later is required for setting distinct viewing options per
#'   structure - such as user defined colors (e.g. one per structure) etc.
#'
#' @returns a character vector of structure data.
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#'
#' @seealso \code{view.pdbs()}, \code{pdb2string()}, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}, \code{bio3d::read.fasta.pdb()}, \code{bio3d::pdbaln()}.
#'
#' @export
#'
#' @examples
#'   pth <- "~/Desktop/courses/BIMM143/class10/pdbs/split_chain/"
#'   files <- list.files(path=pth, full.names = TRUE)
#'   pdbs <- bio3d::pdbaln(files, fit=TRUE, exefile="msa")
#'
#'  NGLVieweR::NGLVieweR( bio3d::pdbs2pdb(pdbs), format="pdb") |>
#'     NGLVieweR::addRepresentation("cartoon")
#'
#'  # Or more simpley...
#'  view.pdbs(pdbs)
#'  # Trace, tube, line, cartoon, ball+stick
#'  view.pdbs(pdbs, representation = "trace")
#'  view.pdbs(pdbs, cols = c("red","blue") )
#'  view.pdbs(pdbs, colorScheme = "residueindex")
#'
pdbs2string <- function(pdbs, collapse=TRUE) {
  z <- bio3d::pdbs2pdb(pdbs)#, rm.gaps = TRUE)
  x <- NULL
  for(i in 1:length(pdbs$id)) {
    x <- c(x, pdb2string( pdb=z[[i]]) )
  }
  # Rtn all in one element or element per pdb
  if(collapse) {
    return( paste(x, collapse = "\n") )
  } else{ return(x) }
}


#pdbs2string <- function(pdbs) {
#  z <- pdbs2pdb(pdbs, inds = 1)
#  x <- pdb2string(pdb=z[[1]], xyz=pdbs$xyz)
#  return(x)
#}


#' Convert a bio3d PCA or NMA object to a character vector for NGLVieweR input
#'
#' @param pc the results of principal component analysis or normal mode analysis as obtained with `pca()` or `nma()` and friends.
#' @param ... Extra arguments passed to `mktrj()`.
#'
#' @returns a character vector of structure data.
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#'
#' @seealso \code{view.pdb()}, \code{pdb2string()}, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::nma()}, \code{bio3d::pca()}, \code{bio3d::mktrj()}.
#'
#' @export
#'
#' @examples
#'  #pdb <- bio3d::read.pdb("6s36")
#'  #n <- bio3d::nma(pdb)
#'  #NGLVieweR::NGLVieweR( pca2string(n), format="pdb") |>
#'  #   NGLVieweR::addRepresentation("cartoon")
#'
#'  # Or more simpley...
#'  # view.pdb(n) # Does not work yet!!!
#'  # view.pdb(n, colorScheme = "residueindex")
#'
pca2string <- function(pc, ...) {
  return( paste(utils::capture.output(
    bio3d::mktrj(pc, file="", ...)),
      collapse = "\n") )
}


#' Interactive PDBS Viewing
#'
#' Quick interactive multi-structure ensemble **bio3d** `pdbs` object viewing
#'  using the **NGLVieweR** package.
#'
#' Just like the related `view.pdb()` function the objective here is to speed up
#' the inspection of multi structure ensembles without having to write out many
#' aligned and superposed individual PDB files for viewing in a external
#' molecular viewer (like PyMol, Mol-star, or VMD - all of which are excellent
#' but often overkill for a quick inspection) or write many lines of **NGLVieweR**
#' code that loop through all structures in a given `pdbs` object.
#'
#' @param pdbs a multi-structure `pdbs` object as obtained from `pdbaln()`, `read.fasta.pdb()`, etc.
#' @param cols a vector of colors, typically one entry per structure. If `NULL` then the output of `vmd.colors()` is used.
#' @param colorScheme if not NULL then this over-rides the previous `cols` input argument. Possible values include "residueindex", "modelindex", "sstruc", "bfactor", "chainid", "chainindex", "atomindex", "occupancy"
#' @param representation the representation style, useful values are line, tube, cartoon, trace, and backbone, ball+stick.
#' @param backgroundColor set the display area background color.
#'
#' @returns an **NGLVieweR** display object that can be displayed or further added to using `NGLVieweR::addRepresentation()` and friends.
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#'
#' @seealso \code{bio3d::pdbaln()}, \code{bio3d::read.fasta.pdb()},
#'   \code{pdbs2string()}, \code{view.pdb()},
#'   \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}
#'
#' @export
#'
#' @examples
#'   pth <- "~/Desktop/courses/BIMM143/class10/pdbs/split_chain/"
#'   files <- list.files(path=pth, full.names = TRUE)
#'   pdbs <- bio3d::pdbaln(files, fit=TRUE, exefile="msa")
#'
#'  view.pdbs(pdbs, representation = "cartoon")
#'  # Trace, tube, line, cartoon, ball+stick
#'  view.pdbs(pdbs, representation = "trace")
#'  view.pdbs(pdbs, cols = c("red","blue") )
#'  # Perhaps this should be the default?
#'  view.pdbs(pdbs, colorScheme = "residueindex")
#'
view.pdbs <- function(pdbs,
                      cols=NULL,
                      colorScheme=NULL,
                      representation="cartoon",
                      backgroundColor = "white"){

  # Convert to multi-element character vector
  x <- pdbs2string(pdbs, collapse=FALSE)

  # Make sure styles are valid for NGL
  representation <- style_key_match(representation)

  n.pdbs <- length(pdbs$id)

  if(is.null(colorScheme)) {
    if(is.null(cols)) {
      cols=bio3d::vmd_colors( length(x) )
    }
    if(length(cols) < n.pdbs) {
      warning("Not enough distinct cols for each structure, recycling")
      cols <- rep(cols, length.out=n.pdbs)
    }
    # Names cause JSON/NGL problems
    names(cols) <- NULL
    params <- list(color = cols[1])
  } else{
    cols <- NULL
    colorScheme <- color_key_match(colorScheme)
    #params <- list(colorScheme = "residueindex")
    params <- list(colorScheme = colorScheme)
  }

  model <- NGLVieweR::NGLVieweR( x[1], format = "pdb") |>
    NGLVieweR::stageParameters(backgroundColor = backgroundColor) |>
    NGLVieweR::addRepresentation(representation, param=params)
    #param = list(color = colors[1]))

  for(j in 2:length(x)) {
    if(!is.null(cols)) {
      params <- list(color = cols[j])
    }
    model <- model |>
      NGLVieweR::addStructure(x[j], format="pdb") |>
      NGLVieweR::addRepresentation(representation, param=params)
      #param = list(color = colors[j]))
  }

  return(model)
}


#' Convert a bio3d PDB object for NGLVieweR
#'
#' Function to take a bio3d structure and use in the NGLVieweR package.
#'
#' @param pdb a bio3d pdb object as obtained from `read.pdb()`.
#'
#' @return a single element character vector with return characters as required by `NGLVieweR::NGLVieweR()` function with `format=pdb` option.
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#'
#' @seealso \code{view.pdb()} which uses this function internally, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}.
#'
#' @export
#'
#' @examples
#' pdb <- bio3d::read.pdb("5p21")
#' NGLVieweR::NGLVieweR(pdb2ngl(pdb), format="pdb") |>
#'   NGLVieweR::addRepresentation("cartoon")
#'
pdb2ngl <- function(pdb) {
  temp <- tempfile()
  on.exit(unlink(temp))
  bio3d::write.pdb(pdb = pdb, file = temp)
  return( paste(readLines(temp), collapse = "\n") )
}

#' Internal utility function for matching available NGL colorSchemes
#'
#' @description
#' This function is for internal use only. No guarantee is provided for its
#' stability across package versions.
#'
#' @keywords internal
#'
#'
color_key_match <- function(colorScheme) {
  # Make sure only valid NGL colorScheme options are used
  # See http://nglviewer.org/ngl/api/classes/colormaker.html
  keyword <- switch(tolower(colorScheme),
                    res = ,
                    resno = ,
                    residue = ,
                    resindex = ,
                    residueindex = "residueindex",
                    resid = ,
                    resname = "resname",
                    model = ,
                    modelindex = "modelindex",
                    ss = ,
                    sse = ,
                    structure = ,
                    secondary = ,
                    sstruc = "sstruc",
                    b = ,
                    bfac = ,
                    tempfactor = ,
                    bfactor = "bfactor",
                    chain = ,
                    chainid = "chainid",
                    chainnum = ,
                    chainindex = "chainindex",
                    chainname = "chainname",
                    eleno = ,
                    atomindex = "atomindex",
                    o = ,
                    occ = "occupancy",
                    hydrophobicity = "hydrophobicity",
                    random = "random",
                    atom = ,
                    elety = ,
                    element = "element",
                    valadation = "valadation",
                    #stop("Invalid colorScheme: Must be one of the allowed types")
                    "residueindex" # Default value
  )
  return(keyword)
}

#' Internal utility function for matching available NGL representations
#'
#' @description
#' This function is for internal use only. No guarantee is provided for its
#' stability across package versions.
#'
#' @keywords internal
#'
#'
style_key_match <- function(representation) {
  # Make sure only valid NGL representation options are used
  # See http://nglviewer.org/ngl/api/
  keyword <- switch(tolower(representation),
                    richardson = ,
                    cartoon = "cartoon",
                    ribbon = "ribbon",
                    tube = "tube",
                    smooth = ,
                    trace = "trace",
                    cpk = ,
                    ballandstick = ,
                    'ball+stick' = "ball+stick",
                    lines = ,
                    line = "line",
                    ca = ,
                    simple = ,
                    backbone = "backbone",
                    licorice = "licorice",
                    surface = "surface",
                    point= "point",
                    rope = "rope",
                    vdw = ,
                    spacefill = "spacefill",
                    #"rocket" # helix cylinder
                    #contact
                    #axis
                    #dot
                    #"helixorient"
                    #"hyperball"
                    #label
                    #stop("Invalid representation: Must be one of the allowed types")
                    "ball+stick" # Default value line
  )
  return(keyword)
}


#' Interactive PDB Viewing
#'
#' Generate a quick NGL (webGL based) structure overview of **bio3d** `pdb` class
#'  objects with a number of simple defaults. The returned **NGLVieweR** object
#'  can be further added to for custom interactive visualizations.
#'
#'  The purpose of this function is to quickly view a given PDB
#'  structure object without having to write output files for input to molecular
#'  graphics programs or write multiple lines of NGLVieweR code.
#'
#'  The extra argument `highlight` takes a **bio3d** `atom.select()` object to highlight
#'  atom subsets as a given `highlight.style`. Useful `highlight.style` options include:
#'  *"ball+stick"*, *"spacefill"*, *"licorice"*, *"line"*, *"surface"* and
#'  *"ribbon"*.
#'
#' The full set of supported `representation` styles include:
#'  - "cartoon",
#'  - "ribbon",
#'  - "tube",
#'  - "trace",
#'  - "ball+stick",
#'  - "line",
#'  - "backbone",
#'  - "licorice",
#'  - "surface",
#'  - "point",
#'  - "rope",
#'  - "spacefill".
#'
#'  Available `colorScheme` options include:
#'  - "residueindex",
#'  - "resname",
#'  - "modelindex",
#'  - "sstruc",
#'  - "bfactor",
#'  - "chainid",
#'  - "chainindex",
#'  - "chainname",
#'  - "atomindex",
#'  - "occupancy",
#'  - "hydrophobicity",
#'  - "random",
#'  - "element",
#'  - "validation",
#'
#'  Limitations:
#'  Currently the function attempts to map common abbreviations and acronyms of
#'  `colorScheme` and `representation` values (e.g. "cpk" for "ballandstick",
#'  "sse" for "sstruct", etc.). If a mapping is not available then the default
#'  values of "residueindex" and "ball+stick" are used respectively.
#'  Currently the function does not check for bio3d or NGLVieweR
#'  availability. It also does not work well with trajectory
#'  or `pdbs` objects. For the later use `view.pdbs()`.
#'
#'
#' @param pdb a PDB structure object as obtained from `read.pdb()`.
#' @param cols a vector of colors, typically one entry per chain. If `NULL` then `colorScheme` is used. If numeric then output of `vmd.colors(cols)` is used.
#' @param colorScheme keyword based coloring used only if `cols` input is NULL. Possible values include "residueindex", "modelindex", "sstruc", "bfactor", "chainid", "chainindex", "atomindex", and "occupancy".
#' @param representation the representation style, useful values are "line", "tube", "cartoon", "trace", "backbone", and "ball+stick".
#' @param backgroundColor a single element color vector that set the display area background color.
#' @param ligand logical, if TRUE ligands will be rendered as atom colored `ligand.style`.
#' @param ligand.style the representation style to use for ligands.
#' @param highlight an optional `atom.select()` object for highlighting.
#' @param highlight.style the representation style to use for selected `highlight` atoms.
#' @param water.rm logical, if TRUE water molecules are removed pior to viewing.
#'
#' @return an **NGLVieweR** display object that can be displayed or further added to using `NGLVieweR::addRepresentation()` and friends.
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#'
#' @seealso \code{pdb2string()}, \code{view.pdbs()},
#'   \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}
#'
#' @export
#'
#' @examples
#'  pdb <- bio3d::read.pdb("1hsg")
#'  view.pdb(pdb)
#'  view.pdb(pdb, ligand=FALSE, cols=c("pink","aquamarine"))
#'  view.pdb(pdb, colorScheme = "sstruc")
#'
#'  #ras <- read.pdb("5p21")
#'  #view.pdb(ras)
#'
#'  sele <- bio3d::atom.select(pdb, resno=c(25, 50))
#'  view.pdb(pdb, highlight = sele,
#'         cols = c("navy","orange"),
#'         backgroundColor = "pink",
#'         highlight.style = "spacefill")
#'
view.pdb <- function(pdb,
                     cols = NULL,
                     colorScheme = "residueindex",
                     representation = "cartoon", # mv "style"?
                     ligand=TRUE,
                     ligand.style = "licorice",
                     backgroundColor = "white",
                     highlight=NULL,
                     highlight.style="ball+stick",
                     water.rm=FALSE) {

  # Should we strip water by default?
  if(water.rm) {
    pdb <- bio3d::atom.select(pdb, string="water", inverse=TRUE, value=TRUE)
  }

  # Convert to character vector
  x <- pdb2string(pdb)


  # Make sure styles are valid for NGL
  representation  <- style_key_match(representation)
  highlight.style <- style_key_match(highlight.style)
  ligand.style    <- style_key_match(ligand.style)

  # Find ligand resid/resn
  lig <- bio3d::atom.select(pdb, "ligand", value=TRUE)
  lig.resid <- unique(lig$atom$resid)

  # Do we have multiple protein chains
  #chains <- unique(pdb$atom$chain)
  chains <- paste(":", unique(bio3d::atom.select(pdb, "protein", value=TRUE)$atom$chain),sep="")
  multi.chain <- length( chains ) > 1

  # Setup the initial viewer and stage
  model <- NGLVieweR::NGLVieweR(x, format="pdb") |>
    NGLVieweR::stageParameters(backgroundColor = backgroundColor)

  # Two options for Color setup:
  #   1) by `colorScheme` keyword
  #   2) by chain from `cols` vector

  col.chain   <- FALSE
  colorScheme <- color_key_match(colorScheme)
  params <- list(colorScheme = colorScheme)

  if(!is.null(cols)) {
    # Override colorScheme for chain coloring
    col.chain <- TRUE

    if(is.numeric(cols)) {
      # Color by number using VMD palette
      cols <- bio3d::vmd_colors()[cols]
    }

    # Names cause JSON/NGL problems
    names(cols) <- NULL

    params <- list(color = cols[1])
  }

  if(!col.chain) {
    # Option 1. keyword coloring (colorScheme)
    model <- model |>
      NGLVieweR::addRepresentation(representation, param=params)
  } else {
    # Option 2. chain coloring
    if(length(cols) < length(chains)) {
      warning("Not enough distinct cols for each chain, recycling")
      cols <- rep(cols, length.out=length(chains))
    }
    # Loop through each chain
    for(j in 1:length(chains)) {
      model <- model |>
        NGLVieweR::addRepresentation(representation,
          param = list(color = cols[j], sele=chains[j]))
    }

  }

  # Add ligand as licorice/"spacefill"
  if(length(lig.resid) > 0) {
    if(ligand) {
      model <- model |> NGLVieweR::addRepresentation(ligand.style,
        param=list(sele=paste(lig.resid, collapse=" "), radius=0.3) )
    } else {
      if(!is.null(lig.resid)) {
        message("Potential ligands found but not displayed, use ligand=T to view")
      }
    }
  }

  # Highlight atoms based on atom.select() object
  if(!is.null(highlight)) {
    #highlight <- atom.select(pdb, resno=90)
    highlight.eleno <- pdb$atom[highlight$atom,]$eleno
    eleno <- paste(paste("@", highlight.eleno, sep=""),collapse = " ")
    model <- model |>
      NGLVieweR::addRepresentation(highlight.style, param=list(sele=eleno))
  }

  return(model)
}



#' Quick molecular viewing of bio3d PCA results
#'
#' @param pc an object of class "pca" as obtained with function `pca.xyz` or `pca`.
#' @param colorScheme keyword based coloring used only if `cols` input is NULL. Possible values include "residueindex", "modelindex", "sstruc", "bfactor", "chainid", "chainindex", "atomindex", and "occupancy".
#' @param representation the representation style, useful values are "line", "tube", "cartoon", "trace", "backbone", and "ball+stick".
#' @param backgroundColor a single element color vector that set the display area background color.
#' @param ... additional arguments passed to and from functions (e.g. to function `write.pdb`).
#'
#' @returns an **NGLVieweR** display object that can be displayed or further added to using `NGLVieweR::addRepresentation()` and friends.
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#'
#' @export
#'
#' @examples
#' #pdb <- bio3d::read.pdb("6s36")
#' #n <- bio3d::nma(pdb)
#' #view.nma(n)
#' #view.nma(n, colorScheme = "model")
view.pca <- function(pc,
                     colorScheme = "residueindex",
                     representation = "cartoon",
                     backgroundColor = "white", ...) {

  x <- pca2string(pc, ...)

  # Make sure styles are valid for NGL
  representation <- style_key_match(representation)
  colorScheme    <- color_key_match(colorScheme)
  params <- list(colorScheme = colorScheme)

  # Setup the initial viewer and stage
  model <- NGLVieweR::NGLVieweR(x, format="pdb") |>
    NGLVieweR::stageParameters(backgroundColor = backgroundColor)|>
    NGLVieweR::addRepresentation(representation, param=params)

  return(model)
}

#' @describeIn view.pca Quick viewing of NMA modes
#' @export
view.nma <- view.pca

#' @describeIn view.pdb An alternative pdb viewer function
#' @param chain.colors Color vector parameter specific to coloring by chain
#' @export
view.old <- function(pdb,
                     ligand=TRUE,
                     chain.colors=bio3d::vmd_colors(),
                     backgroundColor = "white", # model = NULL
                     highlight=NULL,
                     highlight.style="ball+stick") {

  # Find ligand resid/resn
  lig <- bio3d::atom.select(pdb, "ligand", value=TRUE)
  lig.resid <- unique(lig$atom$resid)

  # Do we have multiple protein chains
  #chains <- unique(pdb$atom$chain)
  chains <- paste(":", unique(bio3d::atom.select(pdb, "protein", value=TRUE)$atom$chain),sep="")
  multi.chain <- length( chains ) > 1

  # Setup the initial viewer and stage
  model <- NGLVieweR::NGLVieweR( pdb2string(pdb), format="pdb") |>
    NGLVieweR::stageParameters(backgroundColor = backgroundColor)

  if(!multi.chain) {
    # Color N-C term spectrum
    model <- model |> NGLVieweR::addRepresentation("cartoon",
                                                   param = list(colorScheme = "residueindex") )
  } else {
    ## Color by chain
    ## Not clear how to set custom colors when using NGLs colorScheme = "chainid"
    ## so we instead loop through each chain with a for loop here:
    #model <- model |> NGLVieweR::addRepresentation("cartoon",
    #    param = list(colorScheme = "chainid") )

    # Ensure chain.colors has no names set
    names(chain.colors) <- NULL

    for(j in 1:length(chains)) {
      model <- model |>
        NGLVieweR::addRepresentation("cartoon",
                                     param = list(color = chain.colors[j], sele=chains[j]))
    }

  }

  # Add ligand as licorice/"spacefill"
  if(length(lig.resid) > 0) {
    if(ligand) {
      model <- model |> NGLVieweR::addRepresentation("licorice",
                                                     param=list(sele=paste(lig.resid, collapse=" "),
                                                                radius=0.3) )
    } else {
      if(!is.null(lig.resid)) {
        message("Potential ligands found but not displayed, use ligand=T to view")
      }
    }
  }

  # Highlight atoms based on atom.select() object
  if(!is.null(highlight)) {
    #highlight <- atom.select(pdb, resno=90)
    highlight.eleno <- pdb$atom[highlight$atom,]$eleno
    eleno <- paste(paste("@", highlight.eleno, sep=""),collapse = " ")
    model <- model |>
      NGLVieweR::addRepresentation(highlight.style, param=list(sele=eleno))
  }

  return(model)
}






#' Quick interactive PDBS object viewing using NGLVieweR
#'
#' @param pdbs a multi-structure `pdbs` object as obtained from `pdbaln()`, `read.fasta.pdb()`, etc.
#' @param colors a vector of colors for each structure. If NULL then the output of `vmd.colors()` is used.
#' @param representation the representation style, usefull values are lines, tube and cartoon.
#' @param backgroundColor set the display area background color.
#'
#' @returns an **NGLVieweR** display object that can be displayed or further added to using `NGLVieweR::addRepresentation()` and friends.
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#'
#' @seealso \code{view.pdb()}, \code{pdb2ngl()}, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}
#'
#' @export
#'
#' @examples
#'   pth <- "~/Desktop/courses/BIMM143/class10/pdbs/split_chain/"
#'   files <- list.files(path=pth, full.names = TRUE)
#'   pdbs <- bio3d::pdbaln(files, fit=TRUE, exefile="msa")
#'
#'  view.pdbs2(pdbs, representation = "cartoon")
#'  view.pdbs2(pdbs, colors = c("red","blue") )
#'
view.pdbs2 <- function(pdbs,
                      colors=NULL,
                      representation="line",
                      backgroundColor = "white"){

  # CAN DELETE THIS ONE - it is superceeded by view.pdbs()
  ## Convert to individual PDB objects
  all.pdbs <- bio3d::pdbs2pdb(pdbs)
  n.pdbs <- length(pdbs$id)

  # Setup default function params
  #colors=NULL #vmd_colors()
  #backgroundColor = "white"
  #highlight=NULL
  #highlight.style="ball+stick"
  #representation="line" # "cartoon"

  if(is.null(colors)) {
    colors <- bio3d::vmd_colors()
  }
  if(length(colors) < n.pdbs) {
    warning("Not ennough distinct colors for each structure, recycling")
    colors <- rep(colors, length.out=n.pdbs)
  }
  # Names cause JSON/NGL problems
  names(colors) <- NULL

  # Setup stage and viewer with first pdb
  model <- NGLVieweR::NGLVieweR( pdb2ngl( all.pdbs[[1]] ), format="pdb") |>
    NGLVieweR::stageParameters(backgroundColor = backgroundColor) |>
    NGLVieweR::addRepresentation(representation, param = list(color = colors[1]))

  # Work through each remaining structure/pdb
  for(k in 2:n.pdbs) {
    model <- model |>
      NGLVieweR::addStructure(pdb2ngl(all.pdbs[[k]]), format="pdb") |>
      NGLVieweR::addRepresentation(representation, param = list(color = colors[k]))
  }
  return(model)
}



