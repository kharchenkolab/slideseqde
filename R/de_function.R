#' @importFrom methods as
#' @importFrom stats lm.fit
#' @import ggplot2
#' @import edgeR
#' @import DESeq2
#' @import ggrepel
#' @importFrom stats model.matrix
#' @importFrom graphics par
#' @importFrom stats qnorm
NULL

#' Context-dependent DE approximation using a linear model
#'
#' @param profiles either a matrix of profiles (columns - cell types, rows - genes)
#'                or a list of such matrices, each corresponding to a replicate
#' @param target1 name of the cell type in context 1 (e.g. Fibroblast_Tumor)
#' @param target2 name of the cell type in context 2 (e.g. Firboblast_Non_Tumor)
#' @param de.method 'edgeR' (default) or 'DESeq2'
#' @param default.bcv default value of
#' @param gene.scale.function function used in library scaling
#' @param library.size.scale library size to scale default is 1e6
#' @param remove.zeroes a flag to remove zeros
#' @return
#' @export
getDifferentialExpression <- function(profiles, target1, target2, 
                   de.method='edgeR.qlf', 
                   default.bcv=0.1 , 
                   gene.scale.function=sqrt , 
                   library.size.scale=1e6,
                   remove.zeroes=FALSE) {
  
  if(!is.list(profiles)) { profiles <- list(profiles)} # convert to a general list representation
  if(is.null(names(profiles))) names(profiles) <- paste0('rep',c(1:length(profiles)))
  for(p in profiles){
      if(!(target1 %in% colnames(p)) | !(target2 %in% colnames(p))){
          print(paste(target1,'or',target2,'present in matrix'))
          return(NULL)
      }
  }

  # TODO: remove genes that were 0s in both target1 and target2
  profiles <- lapply(profiles, function(x){
      x[which(rowSums(x[,c(target1,target2)]) > 0),]
  })
  
  genes <- Reduce(intersect,lapply(profiles,rownames))
  profiles <- lapply(profiles, function(x){x[genes,]})
  # determine gene scales
  gene.scaling <- rowSums(do.call(cbind,lapply(profiles,function(p) rowSums(p[genes,]))))
  gene.scaling <- gene.scale.function(gene.scaling)/gene.scaling
  
  # run linear models
  ml <- lapply(profiles,function(p) {
    # normalize
    # take common genes 
    p <- t(t(p)/colSums(p)*library.size.scale)*gene.scaling
    m <- lm.fit( as.matrix(p[,colnames(p)!=target1]), p[,target1])
  })
  
  
  # construct "adjusted profiles" for DE
  ap <- lapply(names(ml), function(rep) {
     m <- ml[[rep]];
     pos <- names(m$coefficients[m$coefficients>0])
     neg <- names(m$coefficients[!m$coefficients>0])
     p <- profiles[[rep]]
     # take the intersecting genes

     libsizes <- colSums(p)
     t1ap <- p[,target1] + colSums(t(p[,neg])/libsizes[neg]*libsizes[target1]*(-1*m$coefficients[neg]))
     t2ap <- colSums(t(p[,pos])/libsizes[pos]*libsizes[target2]*(m$coefficients[pos]))
     am <- cbind(t2ap,t1ap);
     colnames(am) <- paste(rep,c('t1','t2'),sep='.')
     am
  })
  

  if(de.method %in% c('edgeR.qlf','edgeR.lrt')) {
     mat <- as.matrix(do.call(cbind,ap)); storage.mode(mat) <- 'integer'
     mat.na <- as.matrix(do.call(cbind,lapply(profiles, function(p) {p[,c(target2,target1)]}))); colnames(mat.na) <- colnames(mat);
     storage.mode(mat.na) <- 'integer'
     meta <- data.frame(sample.id=colnames(mat),group=as.factor(rep(c('t1','t2'),length(ap))))
     group <- meta$group
     y <- DGEList(counts=mat,group=group)
     y <- calcNormFactors(y)
     y.na <- DGEList(counts=mat.na,group=group)
     y.na <- calcNormFactors(y.na)
     design <- model.matrix(~group)
     if(length(ap)==1) {
       result <- exactTest(y,dispersion=default.bcv^2)
       tt <- result$table;
       result.na <- exactTest(y.na,dispersion=default.bcv^2)
       tt.na <- result.na$table;
     } else {
       y <- edgeR::estimateDisp(y,design)
       y.na <- edgeR::estimateDisp(y.na,design)
       if(de.method=='edgeR.qlf') {
         fit <- glmQLFit(y,design)
         result <- glmQLFTest(fit,coef=2)
         tt <- result$table
         # non-adjusted
         fit.na <- glmQLFit(y.na,design)
         result.na <- glmQLFTest(fit.na,coef=2)
         tt.na <- result.na$table
       } else {
         fit <- glmFit(y,design)
         result <- glmLRT(fit,coef=2)
         tt <- result$table
         fit.na <- glmFit(y.na,design)
         result.na <- glmLRT(fit.na,coef=2)
         tt.na <- result.na$table
       }
     }
  } else if(de.method %in% c("DESeq2.Wald", "DESeq2.LRT")) {
    mat <- as.matrix(do.call(cbind,ap)); storage.mode(mat) <- 'integer'
    mat.na <- as.matrix(do.call(cbind,lapply(profiles, function(p) {p[,c(target2,target1)]})));
    
    storage.mode(mat.na) <- 'integer'
    meta <- data.frame(sample.id=colnames(mat),group=as.factor(rep(c('t1','t2'),length(ap))))
    group <- meta$group  
    result <- DESeq2::DESeqDataSetFromMatrix(mat,meta,design=~group)
    result.na <- DESeq2::DESeqDataSetFromMatrix(mat.na,meta,design=~group)
    if(de.method=='DESeq2.LRT') {
      result <- DESeq2::DESeq(result,test="LRT", reduced = ~ 1,quiet=TRUE)
      result.na <- DESeq2::DESeq(result.na,test="LRT", reduced = ~ 1,quiet=TRUE)
    } else {
      result <- DESeq2::DESeq(result,test="Wald",quiet=TRUE)
      result.na <- DESeq2::DESeq(result.na,test="Wald",quiet=TRUE)
    }
    tt <- as.data.frame(DESeq2::results(result));
    tt.na <- as.data.frame(DESeq2::results(result.na)); 
    colnames(tt) <- colnames(tt.na) <- c('logCPM','logFC','lfcSE','stat','PValue','pjad');
  } else {
    stop(paste('DE method',de.method,'is not supported'))
  }
  return(list(tt=tt,result=result,tt.na=tt.na,result.na=result.na,p=ap,ml=ml))
}

#' get z values from edgeR DE
#' @param d the data.frame
#' @return updated data.frame
get.z.edge.R = function(d) {
  res.table <- d;
  res.table$Z <- -qnorm(pmin(res.table$PValue,0.5))* sign(res.table$logFC)
  #res.table$Z <- -qnorm(res.table$PValue,0.5)
  # res.table$Z[is.na(res.table$Z)] <- 0
  # res.table$Z[is.infinite(res.table$Z)] <- 0
  res.table
}

#' Scatter plot with some top DE genes along
#' with the tumor markers to show the effect of regression
#' @param res the de object with both the corrected and uncorrected DE results
#' @param red.marker marker genes from the cell-type that we control for (tumor)
#' @param target the target cell-type to create the title of the plot
#' @param top.genes how many DE genes to annotate
#' @return two ggplots appended side by side
#' @export
scatterPlotWithCorrection = function(res, red.marker, target, top.genes){
  etuc <- res$tt.na; et <- res$tt
  #genes <- rownames(et)[order(et$PValue)[1:top.genes]]
  par(mfrow=c(1,2), mar=c(3.5,3.5,0.5,0.5),mgp=c(2,0.65,0),
      oma=c(0,0,2,0))
  sdf <- data.frame(uc=etuc$logFC,c=et$logFC); rownames(sdf) <- rownames(et)
  genes <- sdf[order(-abs(sdf$c))[1:top.genes],] %>% rownames

  p1 <- ggplot(sdf, aes(x=uc,y=c)) +
    geom_point(size=0.05) +
    geom_hline(yintercept=0,linetype="dashed")+
    labs(x='Uncorrected M value', y='Corrected M value') +
    geom_abline(slope = 1, intercept = 0, linetype="dashed") +
    geom_point(data=sdf[ red.marker ,],
               pch=21, fill=NA, size=2, colour="red", stroke=0.6,na.rm=TRUE)+
    theme_classic() +
    ggtitle(target)+
    geom_text_repel(
      data=sdf[ genes ,], # Filter data first
      aes(label = genes ),
      size = 4,
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.5, "lines"),
      max.overlaps = Inf,
      na.rm=TRUE
    ) + theme(plot.title = element_text(size = 8))

  etuc <- get.z.edge.R(etuc); et <- get.z.edge.R(et)
  sdf <- data.frame(uc=etuc$Z,c=et$Z); rownames(sdf) <- rownames(et)
  genes2 <- sdf[order(-abs(sdf$c))[1:top.genes],] %>% rownames

  #genes2 <- rownames(et)[order(et$PValue)[1:top.genes]]
  p2 <- ggplot(sdf, aes(x=uc,y=c)) +
    geom_point(size=0.05) +
    geom_hline(yintercept=0,linetype="dashed")+
    geom_abline(slope = 1, intercept = 0, linetype="dashed") +
    geom_point(data=sdf[ red.marker ,],
               pch=21, fill=NA, size=2, colour="red", stroke=0.6,na.rm=TRUE)+
    theme_classic() +
    ggtitle(target)+
    labs(x='Uncorrected Z value', y='Corrected Z value') +
    geom_text_repel(
      data=sdf[ genes2 ,], # Filter data first
      aes(label = genes2 ),
      size = 4,
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.5, "lines"),
      max.overlaps = Inf,
      na.rm=TRUE
    ) + theme(plot.title = element_text(size = 8))

  p1+p2
}
