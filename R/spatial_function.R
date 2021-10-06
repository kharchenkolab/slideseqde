# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' create KDE boundary for a pagoda object
#' @param p.obj is a pagoda2 object
#' @param cell.type a cell type default value is Tumor
#' @return a modified object with boundary stored in
#'         an object contour misc$cont
#'         and create primary contexts
createKDEBoundary = function(p.obj, cell.type = "Tumor"){
  if (cell.type %in% levels(p.obj$clusters$PCA$annot)){
    cnames <- names(p.obj$clusters$PCA$annot[
      p.obj$clusters$PCA$annot == cell.type]
    )
    if(is.null(cnames)){
      stop("cell names with celltype is 0")
    }
  }else{
    stop(paste(cell.type,"not in annotation"))
  }

  linetype <- 1
  color="red"
  coord = p.obj$embeddings$PCA$physical[cnames,]
  kd <- ks::kde(coord, compute.cont=TRUE)
  lnc <- with(
    kd,
    grDevices::contourLines(
      x= kd$eval.points[[1]],
      y= kd$eval.points[[2]],
      z= kd$estimate,
      levels = kd$cont["10%"]
    )[[1]]
  )

  cont <- ggplot2::geom_path(
    aes(x, y),
    data=data.frame(lnc),
    linetype = linetype ,
    color=color
  )

  p.obj$misc$cont <- cont

  # create tumor vs non tumor context
  emb = p.obj$embeddings$PCA$physical
  tim = sp::point.in.polygon(emb[,1], emb[,2], lnc$x, lnc$y)
  names(tim)=rownames(emb)

  p.obj$clusters$PCA$contexts = factor(tim)
  p.obj
}

#' computer joint context
#' @param p.obj a pagoda object
#' @param annotation.type.name name of the annotation object within pagoda object
#' @param joint.type.name name of the object that would context specific annotation
#' @return modified pagoda2 object
computeJointContext <- function(
  p.obj,
  annotation.type.name = "annot",
  joint.type.name = "jcontext"
){
  if( ((annotation.type.name %in% names(p.obj$clusters$PCA) ) &
       ("contexts" %in% names(p.obj$clusters$PCA)) )
  ){
    context.anot <- data.frame(context = p.obj$clusters$PCA$contexts)
    context.anot$anot <- p.obj$clusters$PCA[[annotation.type.name]][
      rownames(context.anot)
    ]

    context.anot <- context.anot %>%
      mutate(context = str_replace(context, "0","Non_Tumor")) %>%
      mutate(context = str_replace(context, "1","Tumor"))

    context.anot$compound <- paste(context.anot$anot, context.anot$context, sep="_")
    compound.anot <- as.factor(context.anot$compound)
    names(compound.anot) <- rownames(context.anot)
    # Make it a part of tumor.2 as joined context
    p.obj$clusters$PCA[[joint.type.name]] <- compound.anot

    p.obj
  }else{
    stop("error")
  }
}
