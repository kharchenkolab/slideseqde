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

#' inner embedding from the pagoda2 object
#' @param pagoda2 object
#' @param tumor.center coordinate for tentative center for the
#'                     tumor context
#' @param nontumor.center coordinate for tentative center for the
#'                        non tumor context
#' @param ctype the cell type
#' @param mode if we _don't_ fix the centers of the contexts manually
#'             they will be calculated
#' @return a dataframe with the distance from the boundary
get.inner.embedding <- function(
  r,
  tumor.center,
  nontumor.center,
  ctype = NULL,
  mode = "manual"
) {

  if (!is.null(ctype)){
    cnames <- names( r$clusters$PCA$annot[r$clusters$PCA$annot == ctype] )
    ctype.emb.df <- data.frame(r$embeddings$PCA$physical[cnames,])
  }else{
    ctype.emb.df <- data.frame(r$embeddings$PCA$physical)
  }
  if (mode == "auto"){
    tumor.center <- as.vector( cluster::pam(
      ctype.emb.df[ctype.emb.df$context == 1, c('xcoord','ycoord') ],1)$medoids
    )
    nontumor.center <- as.vector( cluster::pam(
      ctype.emb.df[ctype.emb.df$context == 0, c('xcoord','ycoord') ],1)$medoids
    )
  }

  ctype.emb.df$context <- r$clusters$PCA$contexts[rownames(ctype.emb.df)]

  tumor.p1 <- ctype.emb.df %>% filter(context == 1) %>%
    dplyr::select(xcoord, ycoord) %>%
    set_names(c('x','y')) %>%
    as.matrix

  nontumor.p1 <- ctype.emb.df %>% filter(context == 0) %>%
    dplyr::select(xcoord, ycoord) %>%
    set_names(c('x','y')) %>%
    as.matrix

  p1 <- ctype.emb.df %>%
    dplyr::select(xcoord, ycoord) %>%
    set_names(c('x','y')) %>%
    as.matrix

  p2 <- r$misc$cont$data %>% dplyr::select(x,y)

  ctype.emb.df$mindist <- rowMin( raster::pointDistance(p1, p2, lonlat=F) )
  tumor.dist <- raster::pointDistance(tumor.p1, tumor.center,lonlat=F)
  non.tumor.dist <- raster::pointDistance(nontumor.p1, nontumor.center,lonlat=F)

  ctype.emb.df$dist <- 0
  ctype.emb.df[ names(tumor.dist),'dist' ] = tumor.dist
  ctype.emb.df[ names(non.tumor.dist),'dist' ] = non.tumor.dist
  ctype.emb.df
}

#' plot the beads with the distance boundary
#' @param r a pagoda2 object with all the information
#' @param ctype.emb.df the dataframe to be generated from
#'                     `get.inner.embedding`
#' @param distance.from.boundary the distance from the boundary drawn from KDE
#' @param pattern TODO: we have to get rid of multiple patterns
#' @return a ggplot object
plot.with.distance.boundary <- function(
  r,
  ctype.emb.df,
  distance.from.boundary,
  distance.from.center=NULL,
  pattern="pattern1"
){
  if(pattern == "pattern1"){
    p <- ggplot2::ggplot(
      ctype.emb.df %>% dplyr::filter(mindist > distance.from.boundary ) ,
      aes(x=xcoord, y=ycoord, color=mindist)) +
      ggplot2::geom_point(size=0.3)+
      r$misc$cont
  }else{
    x.min <- colMins(r$embeddings$PCA$physical)[[1]]
    y.min <- colMins(r$embeddings$PCA$physical)[[2]]
    x.max <- colMaxs(r$embeddings$PCA$physical)[[1]]
    y.max <- colMaxs(r$embeddings$PCA$physical)[[2]]
    p <- ctype.emb.df %>%
      dplyr::filter((mindist > distance.from.boundary)&(mindist < distance.from.center) ) %>%
      ggplot2::ggplot(aes(x=xcoord, y = ycoord, color=mindist)) +
      ggplot2::geom_point(size=0.1) +
      xlim(x.min, x.max) +
      ylim(y.min, y.max) +
      r$misc$cont
  }
  p
}

#' plot the beads with the distance boundary
#' @param r a pagoda2 object with all the information
#' @param ctype.emb.df the dataframe to be generated from
#'                     `get.inner.embedding`
#' @param distance.from.boundary the distance from the boundary drawn from KDE
#' @param pattern TODO: we have to get rid of multiple patterns
#' @return a ggplot object
get.subset <- function(
  r,
  ctype.emb.df,
  distance.from.boundary,
  distance.from.center=NULL,
  pattern = "pattern1"
){
  if(pattern == "pattern1"){
    cnames <- ctype.emb.df %>%
      dplyr::filter(mindist > distance.from.boundary ) %>%
      rownames
    t <- table(r$clusters$PCA$jcontext.L2[cnames]) %>%
      data.frame %>%
      rename(cell.type=Var1) %>%
      dplyr::mutate(context = ifelse( grepl("_Non_Tumor",cell.type), "Non_Tumor", "Tumor" ) ) %>%
      dplyr::mutate(cell.type = ifelse( grepl("_Non_Tumor",cell.type), gsub('_Non_Tumor','',cell.type),
                                 gsub('_Tumor','',cell.type)) ) %>%
      reshape2::dcast(cell.type ~ context, value.var = "Freq") %>%
      dplyr::filter( (Non_Tumor > 10) & (Tumor > 10) ) %>%
      tidyr::pivot_longer(!cell.type, names_to="context",values_to="Freq") %>%
      ggplot2::ggplot(aes(x=cell.type,y=Freq,fill=as.factor(context))) +
      ggplot2::geom_bar(stat='identity') +
      ggplot2::theme_classic() +
      ggplot2::theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }else{
    cnames <- ctype.emb.df %>%
      dplyr::filter((mindist > distance.from.boundary)&(mindist < distance.from.center) ) %>%
      rownames
    t <- table(r$clusters$PCA$jcontext.L2[cnames]) %>%
      data.frame %>%
      rename(cell.type=Var1) %>%
      dplyr::mutate(context = ifelse( grepl("_Non_Tumor",cell.type), "Non_Tumor", "Tumor" ) ) %>%
      dplyr::mutate(cell.type = ifelse( grepl("_Non_Tumor",cell.type), gsub('_Non_Tumor','',cell.type),
                                 gsub('_Tumor','',cell.type)) ) %>%
      reshape2::dcast(cell.type ~ context, value.var = "Freq") %>%
      dplyr::filter( (Non_Tumor > 10) & (Tumor > 10) ) %>%
      tidyr::pivot_longer(!cell.type, names_to="context",values_to="Freq") %>%
      ggplot2::ggplot(aes(x=cell.type,y=Freq,fill=as.factor(context))) +
      ggplot2::geom_bar(stat='identity') +
      ggplot2::theme_classic() +
      ggplot2::theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  }

  list(cnames=cnames,barplot=t)
}


#' make `dgCMatrix` for a pseudo bulk profiles
#' @param counts a sparse matrix instance where the rows are genes
#'               and the columns are cell names
#' @param cell.groups a list of cell names with cell type with
#'                    annotations
#' @return return the collapsed `dgCMatrix` with gene x cell-types
collapse.clusters <- function(counts,cell.groups) {
  cn <- intersect(names(cell.groups),colnames(counts))
  cell.groups <- as.factor(cell.groups[cn])
  if(length(cn)<2) stop("the matrix and the factor have too few common cells")
  x <- t(sccore::colSumByFactor(t(counts[,cn]),cell.groups)[-1,,drop=FALSE]) # remove the first NA column
  colnames(x) <- levels(cell.groups); rownames(x) <- rownames(counts);
  return(as(x,'dgCMatrix'))
}

#' construct the profiles within a construct
#' @param r the pagoda object
#' @param r.raw the SparseSpatial object with raw counts
#' @param cnames the cell names
#' @param cutoff the cutoff number of cells per gene type
#' @param special.genes the special genes that are needed to be included
#' @return
construct.subset.profile <- function(
  r,
  r.raw,
  cnames,
  cutoff,
  special.genes = NULL
  ){
  ad.mix.ctypes <- table(r$clusters$PCA$jcontext.L2[cnames]) %>%
    data.frame %>%
    rename(cell.type=Var1) %>%
    mutate(context = ifelse( grepl("_Non_Tumor",cell.type), "Non_Tumor", "Tumor" ) ) %>%
    mutate(cell.type = ifelse( grepl("_Non_Tumor",cell.type), gsub('_Non_Tumor','',cell.type),
                               gsub('_Tumor','',cell.type)) ) %>%
    reshape2::dcast(cell.type ~ context, value.var = "Freq") %>%
    filter( (Non_Tumor > cutoff) & (Tumor > cutoff) ) %>%
    dplyr::select(cell.type) %>% unlist %>% as.vector

  ad.mix.ctypes <- union(
    paste(ad.mix.ctypes,"Tumor",sep="_"),
    paste(ad.mix.ctypes,"Non_Tumor",sep="_")
  )

  cname.anot <- r$clusters$PCA$jcontext.L2[cnames]
  counts = r.raw$cm[colnames(r$counts),]
  counts = counts[rowSums(counts) > 50,]
  final.gene.names <- rownames(counts)
  cname.anot <- cname.anot[cname.anot %in% ad.mix.ctypes]

  rprofile <- collapse.clusters(
    r.raw$cm[ final.gene.names, names(cname.anot)],
    cname.anot
  )

  if(!is.null(special.genes)){
    special.genes = intersect(rownames(rprofile), special.genes)
  }
  valid.genes = names(sort(rowSums( rprofile ),decreasing = TRUE))[1:7e3]
  if(!is.null(special.genes)){
    valid.genes = union(valid.genes, special.genes)
    # print(length(valid.genes))
  }
  rprofile <- rprofile[valid.genes,colSums(rprofile) > 0]
  rprofile
}


#' Prepare aggregated profiles
#' @param r the pagoda object
#' @param r.raw the SparseSpatial object with raw counts
#' @param distance.from.boundary the cell names
#' @param tumor.center the coordinates for the tumor center
#' @param nontumor.center the coordinates for the nontumor center
#' @param special.genes the special genes that are needed to be included
#' @return the collapsed profiles and the dataframe, distance, plots and
#'         the cell names
preparePseudoBulkProfile = function(
  r,
  r.raw,
  distance.from.boundary,
  tumor.center,
  nontumor.center,
  special.genes=NULL
){
  ctype.emb.df <- get.inner.embedding(
    r,
    tumor.center,
    nontumor.center
  )

  p <- plot.with.distance.boundary(r, ctype.emb.df, distance.from.boundary)
  subset <- get.subset(r, ctype.emb.df, distance.from.boundary)
  cnames <- subset$cnames
  p2 <- subset$barplot

  rprofile <- construct.subset.profile(r, r.raw, cnames, 10, special.genes)
  ctype.emb.df$cell.type <- r$clusters$PCA$annot.L2[rownames(ctype.emb.df)]
  list(rprofile=rprofile,
       ctype.emb.df = ctype.emb.df,
       dist=distance.from.boundary ,
       plots = p | p2,
       cnames=cnames
       )
}
