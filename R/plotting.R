#' plotKdModel
#'
#' Plots the summary of an affinity model. Requires the packages `seqLogo` and 
#' `cow_plot`.
#'
#' @param mod A `KdModel`
#' @param what Either 'seeds', 'logo', or 'both' (default)
#'
#' @return If `what="logo"`, returns nothing a plots a position weight matrix.
#' Otherwise returns a ggplot.
#' @import ggplot2
#' @export
plotKdModel <- function(mod, what=c("both","seeds","logo"), n=10){
  what <- match.arg(what)
  if(what=="seeds"){
    mer8 <- getSeed8mers(mod$canonical.seed)
    wA <- which(substr(mer8,8,8)=="A")
    mer7 <- substr(mer8,1,7)
    As <- mod$mer8[wA]
    names(As) <- mer7[wA]
    mer.mean <- rowsum(mod$mer8[-wA],mer7[-wA])[,1]/3
    As <- As-mer.mean[names(As)]
    d <- data.frame(seed=names(mer.mean), base=mer.mean/-1000, "A"=As[names(mer.mean)]/-1000,
                    type=getMatchTypes(names(mer.mean),mod$canonical.seed), row.names=NULL)
    d <- d[head(order(d$base+d$A, decreasing=TRUE),n=n),]
    d$seed <- factor(as.character(d$seed), rev(as.character(d$seed)))
    levels(d$type) <- .matchLevels(FALSE)
    d2 <- data.frame(seed=rep(d$seed,2), log_kd=c(d$base,d$A), type=c(as.character(d$type), rep("+A",n)))
    p <- ggplot(d2, aes(seed, log_kd, fill=type)) + geom_col() + coord_flip() + 
      ylab("-log_kd") + ggtitle(mod$name)
    if(mod$name != mod$mirseq) p <- p + labs(subtitle=gsub("T","U",mod$mirseq))
    return( p )
  }
  
  if(what=="logo") return(seqLogo::seqLogo(mod$pwm, xfontsize=12, yfontsize=12, xaxis=FALSE))
  cowplot::plot_grid( plotKdModel(mod, "seeds"),
                      grid::grid.grabExpr(plotKdModel(mod, "logo")),
                      nrow=2, rel_heights=c(6,4))
}


#' viewTargetAlignment
#'
#' @param m A GRanges of length 1 giving the information for a given match, as
#' produced by \link{\code{findSeedMatches}}.
#' @param miRNA A miRNA sequence, or a \link{\code{KdModel}} object of the miRNA
#' corresponding to the match in `m`; alternatively, a \link{\code{KdModelList}}
#' including the model.
#' @param seqs The sequences corresponding to the seqnames of `m`. Not needed if
#' `m` contains the target sequences.
#' @param flagBulgeMatches Logical; whether to flag matches inside the bulge
#' (default FALSE)
#' @param maxBulgeSize The maximum bulge size to consider (default none)
#' @param maxBulgeDiff The maximum difference between miRNA and target bulges
#'
#' @return A data.frame containing the aligned sequences.
#' @importFrom stringi stri_reverse
#' @export
viewTargetAlignment <- function(m, miRNA, seqs=NULL, flagBulgeMatches=FALSE,
                                maxBulgeSize=Inf, maxBulgeDiff=Inf){
  stopifnot(is(m,"GRanges"))
  stopifnot(length(m)==1)
  if(is.list(miRNA) && is(miRNA[[1]],"KdModelList") && !is.null(m$miRNA)){
    stopifnot(as.character(m$miRNA) %in% names(miRNA))
    miRNA <- miRNA[[as.character(m$miRNA)]]
  }
  if(is(miRNA,"KdModel")){
    miRNA <- miRNA$mirseq
  }else{
    stopifnot(is.character(miRNA) && length(miRNA)==1)
  }
  if(!("p3.mir.bulge" %in% colnames(mcols(m))))
    stop("`m` does not contain detailed 3' alignment information. Re-run the ",
         "scan with p3.extra=TRUE to be able to use this features.")
  if(!("sequence" %in% colnames(mcols(m)))){
    if(is.null(seqs))
      stop("`m` does not contain target RNA sequences. Please provide them ",
         "via the `seqs` argument.")
    stopifnot(as.character(seqnames(m)) %in% names(seqs))
    seqs <- DNAString(seqs[[as.character(seqnames(m))]])
    bulgeDiff <- max(m$p3.target.bulge-m$p3.mir.bulge,0L)
    r <- IRanges(max(start(m)-bulgeDiff-nchar(mirseq)+8L,1L), end(m)+2L)
    m$sequence <- as.character(unlist(extractAt(seqs, r)))
  }
  if(m$p3.mir.bulge>maxBulgeSize | m$p3.target.bulge>maxBulgeSize |
     abs(m$p3.mir.bulge-m$p3.target.bulge) > maxBulgeDiff ){
    m$p3.mir.bulge <- m$p3.target.bulge <- 3L
  }
  mirseq <- gsub("T","U",miRNA)
  target <- stringi::stri_reverse(gsub("T","U",as.character(m$sequence)))
  mirseq2 <- as.character(complement(RNAString(mirseq)))
  mirseq2 <- paste0("A",substr(mirseq2, 2, nchar(mirseq2)))
  if(grepl("g-bulged",as.character(m$type),fixed=TRUE)){
    # TO DO!
  }
  if(flagBulgeMatches){
    minBulge <- min(m$p3.mir.bulge,m$p3.target.bulge)
  }else{
    minBulge <- 0
  }
  mm <- ifelse(strsplit(substr(mirseq2,1,8+minBulge),"")[[1]]==
                 strsplit(substr(target,3,10+minBulge),"")[[1]], "|", " ")
  if(!flagBulgeMatches && m$p3.mir.bulge==m$p3.target.bulge){
    mm <- c(mm,rep(" ",m$p3.mir.bulge))
  }else if(m$p3.mir.bulge<m$p3.target.bulge){
    mm <- c(mm,rep(" ",m$p3.target.bulge-minBulge))
    mirseq <- paste0(
      paste(substr(mirseq,1,8+m$p3.mir.bulge),collapse=""),
      paste(rep("-",m$p3.target.bulge-m$p3.mir.bulge),collapse=""),
      paste(substr(mirseq,9+m$p3.mir.bulge,nchar(mirseq)),collapse=""))
  }else{
    mm <- c(mm,rep(" ",m$p3.mir.bulge-minBulge))
    target <- paste0(
      paste(substr(target,1,10+m$p3.target.bulge),collapse=""),
      paste(rep("-",m$p3.mir.bulge-m$p3.target.bulge),collapse=""),
      paste(substr(target,11+m$p3.target.bulge,nchar(target)),collapse=""))
  }
  minl <- min(nchar(target),nchar(mirseq2)-8-m$p3.mir.bulge+10+m$p3.target.bulge)
  mirseq2 <- substr(mirseq2,9+m$p3.mir.bulge,nchar(mirseq2))
  target2 <- substr(target,11+m$p3.target.bulge,nchar(target))
  minl <- min(nchar(mirseq2),nchar(target2))
  mirseq2 <- substr(mirseq2,1,minl)
  target2 <- substr(target2,1,minl)
  mm2 <- ifelse(strsplit(mirseq2,"")[[1]]==strsplit(target2,"")[[1]], "|", " ")
  mm <- paste(c(mm,mm2),collapse="", sep="")
  sp5 <- paste0(rep(" ",5),collapse="")
  d <- data.frame(
    row.names=c("miRNA 3'  ",paste0(sp5,sp5),"target 5' "),
    "alignment"=c(stringi::stri_reverse(paste0(sp5,mirseq,"   ")),
                  stringi::stri_reverse(paste0(sp5,mm,"   ")),
                  stringi::stri_reverse(paste0("...",target,"..."))))
  d$alignment <- paste0(sapply(max(nchar(d$alignment))-nchar(d$alignment),
                               FUN=function(x) paste0(rep(" ",x),collapse="")),
                        d$alignment)
  d
}