#' plotKdModel
#'
#' Plots the summary of an affinity model.
#'
#' @param mod A `KdModel`
#' @param what Either 'seeds', 'logo', or 'both' (default). 'logo' and 'both
#' @param n The number of top 7-mers to plot (when `what='seeds'`)
#'
#' @return If `what="logo"`, returns nothing a plots a position weight matrix.
#' Otherwise returns a ggplot.
#'
#' @details
#' `what='seeds'` plots the -log(Kd) values of the top `n` 7-mers (including
#' both canonical and non-canonical sites), with or without the final "A"
#' vis-a-vis the first miRNA nucleotide.
#' `what='logo'` plots a `seqLogo` (requires the
#' [seqLogo](https://bioconductor.org/packages/release/bioc/html/seqLogo.html)
#' package) showing the nucleotide-wise information content and preferences for
#' all 12-mers (centered around the seed). `what="both"` plots both.
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom grid grid.grabExpr
#' @importFrom seqLogo seqLogo
#' @export
#' @examples
#' data(SampleKdModel)
#' plotKdModel(SampleKdModel, what="seeds")
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
    d <- data.frame(seed=names(mer.mean), base=mer.mean/-1000,
                    "A"=As[names(mer.mean)]/-1000,
                    type=getMatchTypes(names(mer.mean),mod$canonical.seed),
                    row.names=NULL)
    d <- d[head(order(d$base+d$A, decreasing=TRUE),n=n),]
    d$seed <- factor(as.character(d$seed), rev(as.character(d$seed)))
    levels(d$type) <- .matchLevels(FALSE)
    d2 <- data.frame(seed=rep(d$seed,2), log_kd=c(d$base,d$A),
                     type=c(as.character(d$type), rep("+A",n)))
    p <- ggplot(d2, aes(seed, log_kd, fill=type)) + geom_col() + coord_flip() +
      ylab("-log(KD)") + xlab("7-mer") + ggtitle(mod$name)
    if(mod$name != mod$mirseq) p <- p + labs(subtitle=gsub("T","U",mod$mirseq))
    return( p )
  }

  if(what=="logo")
    return(seqLogo::seqLogo(mod$pwm, xfontsize=12, yfontsize=12, xaxis=FALSE))
  gridExtra::grid.arrange(plotKdModel(mod, "seeds"),
                          grid::grid.grabExpr(plotKdModel(mod, "logo")),
                          nrow=2, heights=c(6,4))
}


#' viewTargetAlignment
#'
#' @param m A GRanges of length 1 giving the information for a given match, as
#' produced by \code{\link{findSeedMatches}}.
#' @param miRNA A miRNA sequence, or a \code{\link{KdModel}} object of the miRNA
#' corresponding to the match in `m`; alternatively, a \code{\link{KdModelList}}
#' including the model.
#' @param seqs The sequences corresponding to the seqnames of `m`. Not needed if
#' `m` contains the target sequences.
#' @param flagBulgeMatches Logical; whether to flag matches inside the bulge
#' (default FALSE)
#' @param maxBulgeSize The maximum bulge size to consider (default none)
#' @param maxBulgeDiff The maximum difference between miRNA and target bulges
#' @param min3pMatch The minimum 3' alignment for any to be plotted
#' @param UGsub Logical; whether to show U-G matches
#' @param outputType Either 'print' (default, prints to console), 'data.frame',
#' or 'plot'.
#' @param ... Passed to `text` if `outputType="plot"`.
#'
#' @return Returns nothing `outputType="print"`. If `outputType="data.frame"`,
#' returns a data.frame containing the alignment strings; if
#' `outputType="plot"` returns a `ggplot` object.
#' @importFrom stringi stri_reverse
#' @importFrom graphics par text
#' @export
#' @examples
#' data(SampleKdModel)
#' seq <- c(seq1="CGACCCCTATCACGTCCGCAGCATTAAAT")
#' m <- findSeedMatches(seq, SampleKdModel)
#' viewTargetAlignment(m, miRNA=SampleKdModel, seqs=seq)
viewTargetAlignment <- function(m, miRNA, seqs=NULL, flagBulgeMatches=FALSE,
                                maxBulgeSize=9L, maxBulgeDiff=4L,
                                min3pMatch=3L, UGsub=TRUE, ...,
                                outputType=c("print","data.frame",
                                             "plot","ggplot")){
  stopifnot(is(m,"GRanges"))
  stopifnot(length(m)==1)
  outputType <- match.arg(outputType)
  if(is.list(miRNA) && is(miRNA[[1]],"KdModelList") && !is.null(m$miRNA)){
    stopifnot(as.character(m$miRNA) %in% names(miRNA))
    miRNA <- miRNA[[as.character(m$miRNA)]]
  }
  mod <- miRNA
  if(is(miRNA,"KdModel")){
    miRNA <- miRNA$mirseq
  }else{
    stopifnot(is.character(miRNA) && length(miRNA)==1)
    if(nchar(miRNA)<15) stop("`miRNA` should be a full mature miRNA sequence")
  }
  if(is.null(seqs) && ( !("p3.mir.bulge" %in% colnames(mcols(m))) ||
                        !("sequence" %in% colnames(mcols(m))) ) ){
    stop("`m` does not contain target RNA sequences. Please provide them ",
           "via the `seqs` argument.")
  }
  if(!("p3.mir.bulge" %in% colnames(mcols(m)))){
    # re-scan to get additional data
    seq2 <- subseq(DNAStringSet(seqs[[as.character(seqnames(m))]]),
                   max(1L,start(m)-(nchar(miRNA)+maxBulgeSize-8L)), end(m)+2L)
    m <- findSeedMatches(seq2, mod, keepMatchSeq=TRUE, p3.extra=TRUE,
                         p3.params = list(maxMirLoop=maxBulgeSize,
                                          maxTargetLoop=maxBulgeSize,
                                          maxLoopDiff=maxBulgeDiff),
                         verbose=FALSE)
  }
  if(!("sequence" %in% colnames(mcols(m)))){
    stopifnot(as.character(seqnames(m)) %in% names(seqs))
    seqs <- DNAString(seqs[[as.character(seqnames(m))]])
    bulgeDiff <- max(m$p3.target.bulge-m$p3.mir.bulge,0L)
    r <- IRanges(max(start(m)-bulgeDiff-nchar(miRNA)+8L,1L), end(m)+2L)
    m$sequence <- as.character(unlist(extractAt(seqs, r)))
  }
  if(m$p3.mir.bulge>maxBulgeSize | m$p3.target.bulge>maxBulgeSize |
     abs(m$p3.mir.bulge-m$p3.target.bulge) > maxBulgeDiff ){
    m$p3.mir.bulge <- m$p3.target.bulge <- 3L
  }
  mirseq <- gsub("T","U",miRNA)
  if(bulged <- grepl("g-bulged",m$type))
    mirseq <- paste0(substr(mirseq,1,5),"-",substr(mirseq,6,nchar(mirseq)))
  target <- stringi::stri_reverse(gsub("T","U",as.character(m$sequence)))
  target2 <- target
  mirseq2 <- as.character(complement(RNAString(mirseq)))
  mirseq2 <- paste0("A",substr(mirseq2, 2, nchar(mirseq2)))
  if(flagBulgeMatches){
    minBulge <- min(m$p3.mir.bulge,m$p3.target.bulge)
  }else{
    minBulge <- 0
  }
  mm <- .matchStrings(substr(mirseq2,1,8+bulged+minBulge),
                      substr(target,3,10+bulged+minBulge), FALSE )
  if(!flagBulgeMatches && m$p3.mir.bulge==m$p3.target.bulge){
    mm <- c(mm,rep(" ",m$p3.mir.bulge))
  }else if(m$p3.mir.bulge<m$p3.target.bulge){
    mm <- c(mm,rep(" ",m$p3.target.bulge-minBulge))
    mirseq <- paste0(
      paste(substr(mirseq,1,8+bulged+m$p3.mir.bulge),collapse=""),
      paste(rep("-",m$p3.target.bulge-m$p3.mir.bulge),collapse=""),
      paste(substr(mirseq,9+bulged+m$p3.mir.bulge,nchar(mirseq)),collapse=""))
  }else{
    mm <- c(mm,rep(" ",m$p3.mir.bulge-minBulge))
    target <- paste0(
      paste(substr(target,1,10+m$p3.target.bulge),collapse=""),
      paste(rep("-",m$p3.mir.bulge-m$p3.target.bulge),collapse=""),
      paste(substr(target,11+m$p3.target.bulge,nchar(target)),collapse=""))
  }
  minl <- min(nchar(target),
              nchar(mirseq2)-8-m$p3.mir.bulge+10+m$p3.target.bulge)
  mirseq2 <- substr(mirseq2,9+m$p3.mir.bulge,nchar(mirseq2))
  target2 <- substr(target2,11+m$p3.target.bulge,nchar(target2))
  minl <- min(nchar(mirseq2),nchar(target2))
  mirseq2 <- substr(mirseq2,1,minl)
  target2 <- substr(target2,1,minl)
  mm2 <- .matchStrings(mirseq2, target2, UGsub)
  if(min3pMatch>1L &&
     !grepl(paste(rep("|",min3pMatch),collapse=""),
            gsub("-","|",paste(mm2,collapse=""),fixed=TRUE), fixed=TRUE)){
    mm2 <- rep(" ",sum(nchar(mm2)))
  }
  mm <- paste(c(mm,mm2),collapse="", sep="")
  sp <- function(x) paste0(rep(" ",x),collapse="")
  d <- data.frame(
    row.names=c("miRNA ",sp(6),"target"),
    "alignment"=c(paste0(sp(3),"3'-",stringi::stri_reverse(mirseq),"-5'",sp(5)),
                  paste0(sp(8),stringi::stri_reverse(mm),sp(8)),
                  paste0("5'-...",stringi::stri_reverse(target),"...-3'")))
  d$alignment <- paste0(vapply(max(nchar(d$alignment))-nchar(d$alignment),
                               FUN.VALUE=character(1),
                               FUN=function(x) paste0(rep(" ",x),collapse="")),
                        d$alignment)
  if(outputType=="data.frame") return(d)
  d2 <- paste0(paste(c("miRNA ","      ","target"), d$alignment), collapse="\n")
  if(outputType=="ggplot"){
    p <- ggplot(data.frame(x=1,y=1,label=d2), aes(x,y,label=label)) +
      theme_void() + geom_text(family="mono", fontface="bold")
    return(p)
  }else if(outputType=="plot"){
    par(xpd = NA, mar=c(0,0,0,0))
    plot(1, 1, ann=FALSE, bty='n', type='n', xaxt='n', yaxt='n')
    text(x=1, y=1, d2, family="mono", font=2, ...)
  }else{
    cat(paste0("\n",d2,"\n"))
  }
}

.matchStrings <- function(s1, s2, UGsub=TRUE){
  s1 <- strsplit(s1,"")[[1]]
  s2 <- strsplit(s2,"")[[1]]
  mm <- ifelse(s1==s2, "|", " ")
  if(UGsub){
    mm[s1!=s2 & s1 %in% c("U","C") & s2 %in% c("U","C")] <- "-"
  }
  mm <- gsub("\\-$"," ",paste(mm,collapse=""))
  mm <- gsub("\\-[^|]|[^|]\\-","  ", mm)
  mm
}

