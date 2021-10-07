#' @import Matrix
NULL

#' Convert csv file to the sparse matrix
#'
#' @param d_table csv file to be read
#' @return a sparse matrix
to_sparse <- function(d_table){
  # create list (one element for each column) of row indices where values are not zero
  i_list <- lapply(d_table, function(x) which(x != 0))

  #count nonzero values in each column
  counts <- unlist(lapply(i_list, length), use.names = F)
  # create sparse matrix
  sparseMatrix(
    i = unlist(i_list, use.names = F), # row indices column by column
    j = rep(1:ncol(d_table), counts), # column indices (each repeated as many times as needed
    x = unlist(lapply(d_table, function(x) x[x != 0]), use.names = F), # take nonzero values from dt
    dims = dim(d_table), # dimensions do not change
    dimnames = list(NULL, names(d_table)))
}

#' convert the csv file and spatial info
#' @param pw   path received from the RCTD file should
#'             contain the `alignment/MappedDGEForR.csv` which
#'             is a gene to cell count matrix
#'             and another file
#'             `barcode_matching/BeadLocationsForR.csv`
#'             which contains the spatial information
#' @return a list containing a sparse matrix and
#'         a dataframe
#'
#' @export
createSparseSpatialObjFromRCTD = function(pw){
  cm <- data.table::fread(paste(pw,'alignment/MappedDGEForR.csv',sep='/'),header=T)
  gn <- cm$GENE; cm <- cm[,-1]
  cm <- to_sparse(cm);
  rownames(cm) <- gn;
  cm <- as(cm,'dgCMatrix')

  bp <- as.data.frame(data.table::fread(paste(pw,'barcode_matching/BeadLocationsForR.csv',sep='/'),header=T))
  cm <- cm[,colnames(cm) %in% bp$barcodes,drop=F]
  return(list(cm=cm,pos=bp[match(colnames(cm),bp$barcodes),]))
}

#' convert the csv file and spatial info
#' @param gene.matrix.file Direct csv file for gene to cell matrix
#' @param location.file file containing the the dataframe
#' @return a list containing a sparse matrix and
#'         a dataframe
#'
#' @export
createSparseSpatialObj = function(gene.matrix.file, location.file){
  cm <- data.table::fread(gene.matrix.file,header=T)
  gn <- cm$GENE; cm <- cm[,-1]
  cm <- to_sparse(cm);
  rownames(cm) <- gn;
  cm <- as(cm,'dgCMatrix')

  bp <- as.data.frame(data.table::fread(location.file,header=T))
  cm <- cm[,colnames(cm) %in% bp$barcodes,drop=F]
  return(list(cm=cm,pos=bp[match(colnames(cm),bp$barcodes),]))
}

#' Convert the Spatial objects to pagoda
#' call `createSparseObj` or `createSparseSpatialObjFromRCTD`
#' before calling this function
#' @param cdl a list of SparseSpatialObj objects with names
#'
#' @return a list of two objects
#'         object 1 is a list of SparseSpatial objects with modified
#'         names, it contains the raw counts
#'         object 2 is a list of pagoda2 objects
convertSpatialObj2Pagoda = function(cdl){
  cdl <- lapply(cdl,function(z) {
    x <- pagoda2::gene.vs.molecule.cell.filter(z$cm,min.cell.size = 100)
    #z$cm <- x[rowSums(x)>50,]
    z$cm <- x
    z
  })

  cdl <- mapply(function(name, x) {
    colnames(x$cm) <- paste(name,colnames(x$cm),sep='_');
    x$pos$barcodes <- paste(name,x$pos$barcodes,sep='_');
    return(x)
  },
  names(cdl),
  cdl,
  SIMPLIFY=F
  )

  cdl.p2 <- lapply(cdl,
   function(x)
      pagoda2::basicP2proc(
       x$cm,nPcs=30,
       n.cores=30,
       get.tsne = F,
       make.geneknn = F,
       #min.cells.per.gene =20,
       min.transcripts.per.cell = 100
  ))

  list(cdl = cdl, cdl.p2 = cdl.p2)
}

#' adding the cell type annotation in the pagoda objects
#' call `convertSpatialObj2Pagoda` before calling this function
#' @param cdl a list of SparseSpatial
#' @param cdl.p2 a list of pagoda2 objects made from `convertSpatialObj2Pagoda`
#' @param anno an R object containing the annotation file
#'        it assumes a specific format for now and
#'        TODO: we need to change this and make more generic
#' @return a modified pagoda object with the positions and
#'         spatial positions
updateSpatialInfo = function(cdl,cdl.p2,anno){
  mapply(
    function(p2,x) {
      p2$embeddings$PCA$physical <- as.matrix(x$pos[,c(2,3)]);
      rownames(p2$embeddings$PCA$physical) <- x$pos[,1];
      p2$embeddings$PCA$physical <- p2$embeddings$PCA$physical[
        rownames(p2$embeddings$PCA$physical) %in% rownames(p2$counts),
      ];
      tmp.anot <- as.factor(anno$cell.type.L1[intersect(rownames(p2$counts),
                                                        names(anno$cell.type.L1))]);
      p2$clusters$PCA$annot <- tmp.anot;

      tmp.anot <- as.factor(anno$cell.type.L2[intersect(rownames(p2$counts),
                                                        names(anno$cell.type.L2))]);
      p2$clusters$PCA$annot.L2 <- tmp.anot;
      p2
    },
    cdl.p2,
    cdl
  )
}
