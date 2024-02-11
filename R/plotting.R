#' plotKdModel
#'
#' Plots the summary of an affinity model.
#'
#' @param mod A `KdModel`
#' @param what Either 'seeds', 'logo', or 'both' (default).
#' @param n The number of top 7-mers to plot (when `what='seeds'`)
#'
#' @return If `what="logo"`, returns nothing and plots a position weight matrix.
#' Otherwise returns a ggplot.
#'
#' @details
#' `what='seeds'` plots the -$log(K_d)$ values of the top `n` 7-mers (including
#' both canonical and non-canonical sites), with or without the final "A"
#' vis-a-vis the first miRNA nucleotide.
#' `what='logo'` plots a `seqLogo` (requires the
#' [seqLogo]{https://bioconductor.org/packages/release/bioc/html/seqLogo.html}
#' package) showing the nucleotide-wise information content and preferences for
#' all 12-mers (centered around the seed, oriented in the direction of the target
#' mRNA). `what="both"` plots both.
#' Note that if the package `ggseqlogo` is installed, this will be used instead
#' to plot the logo, resulting in more detailed plot annotation.
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom seqLogo seqLogo
#' @export
#' @examples
#' data(SampleKdModel)
#' plotKdModel(SampleKdModel, what="seeds")
plotKdModel <- function(mod, what=c("both","seeds","logo"), n=10){
  stopifnot(is(mod,"KdModel"))
  what <- match.arg(what)
  if(what=="seeds"){
    type_cols <- .typeColors()
    mirseq2 <- strsplit(gsub("T","U",mod$mirseq),"")[[1]]
    mirseq2 <- paste0("3'-",stringi::stri_reverse(gsub("T","U",mod$mirseq)),"-5'")
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
    d$seed <- gsub("T","U",d$seed)
    d$seed <- factor(as.character(d$seed), rev(as.character(d$seed)))
    levels(d$type) <- .matchLevels(FALSE)
    d2 <- data.frame(seed=rep(d$seed,2), log_kd=c(d$base,d$A),
                     type=c(as.character(d$type), rep("+A",n)))
    p <- ggplot(d2, aes(seed, log_kd, fill=type)) + geom_col() + coord_flip() +
      ylab(bquote("-"*log(K[D]))) + xlab("7-mer") + ggtitle(mod$name) + 
      theme_minimal() +
      theme(axis.text.x = element_text(size=11),
            axis.text.y = element_text(size=11, family="mono", face = "bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
    p <- p + scale_fill_manual(values = type_cols[p$data$type])
    if(mod$name != mod$mirseq) p <- p + labs(subtitle=mirseq2)
    return( p )
  }

  if(what=="logo"){
#    if(requireNamespace("ggseqlogo", quietly=TRUE)){
#      p <- suppressWarnings(ggseqlogo( mod$pwm )) + 
#        labs(y="Information (bits)", x="miRNA 5' position")
#      p$scales$scales[[1]] <- scale_x_continuous(breaks=1:12,labels=c(10:1,"",""))
#      p <- p + 
#      theme(axis.text.x=element_text(size=11), axis.text.y=element_text(size=11),
#            axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))
#    }else{
      p <- plot_grid(grid::grid.grabExpr(seqLogo( mod$pwm )))
#    }
    return(p)
  }
  plot_grid(plotKdModel(mod, "seeds"),plotKdModel(mod, "logo"),
            ncol=1, rel_heights = c(6,4))
}

.typeColors <- function(){
  c("8mer"="darkred","7mer-m8"="#44AA99","7mer-a1"="#117733", "7mer"="#117733",
    "6mer"="#4CAF50","6mer-m8"="#449d48", "6mer-a1"="#3c8c40",
    "g-bulged 8mer"="#5bb9e8", "g-bulged 7mer"="#88CCEE",
    "g-bulged 6mer"="#b5dff4", "wobbled 7mer"="grey85",
    "wobbled 8mer"="grey75", "non-canonical"="grey55", "+A"="darkred")
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
#' @param p3.params See \code{\link{findSeedMatches}}.
#' @param min3pMatch The minimum 3' alignment for any to be plotted
#' @param UGsub Logical; whether to show U-G matches
#' @param hideSingletons Logical; whether to hide isolated single base-pair
#' matches
#' @param outputType Either 'print' (default, prints to console), 'data.frame',
#' or 'plot'.
#' @param ... Passed to `text` if `outputType="plot"`.
#'
#' @return Returns nothing `outputType="print"`. If `outputType="data.frame"`,
#' returns a data.frame containing the alignment strings; if
#' `outputType="ggplot"` returns a `ggplot` object.
#' @importFrom stringi stri_reverse
#' @importFrom graphics par text
#' @export
#' @examples
#' data(SampleKdModel)
#' seq <- c(seq1="CGACCCCTATCACGTCCGCAGCATTAAAT")
#' m <- findSeedMatches(seq, SampleKdModel, verbose=FALSE)
#' viewTargetAlignment(m, miRNA=SampleKdModel, seqs=seq)
viewTargetAlignment <- function(m, miRNA, seqs=NULL, flagBulgeMatches=FALSE,
                                p3.params=list(), min3pMatch=3L,
                                hideSingletons=FALSE, UGsub=TRUE, ...,
                                outputType=c("print","data.frame",
                                             "plot","ggplot")){
  stopifnot(is(m,"GRanges"))
  stopifnot(length(m)==1)
  outputType <- match.arg(outputType)
  p3.params <- .check3pParams(p3.params)
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
    maxBulgeSize <- max(p3.params$maxMirLoop, p3.params$maxTargetLoop)
    internStart <- max(1L,start(m)-(nchar(miRNA)+maxBulgeSize-8L))
    seq2 <- subseq(DNAStringSet(seqs[as.character(seqnames(m))]),
                   internStart, end(m)+2L)
    tx_start <- start(m)
    m <- findSeedMatches(seq2, mod, keepMatchSeq=TRUE, p3.extra=TRUE,
                         p3.params=p3.params,
                         maxLogKd=0, minDist=-Inf, verbose=FALSE)
    m <- m[tx_start==start(m)+internStart-1L]
  }
  if(!("sequence" %in% colnames(mcols(m)))){
    stopifnot(as.character(seqnames(m)) %in% names(seqs))
    seqs <- DNAString(seqs[[as.character(seqnames(m))]])
    bulgeDiff <- max(m$p3.target.bulge-m$p3.mir.bulge,0L)
    r <- IRanges(max(start(m)-bulgeDiff-nchar(miRNA)+8L,1L), end(m)+2L)
    m$sequence <- as.character(unlist(extractAt(seqs, r)))
  }
  do3p <- TRUE
  if(m$p3.mir.bulge>p3.params$maxMirLoop |
     m$p3.target.bulge>p3.params$maxTargetLoop |
     abs(m$p3.mir.bulge-m$p3.target.bulge) > p3.params$maxLoopDiff ){
    m$p3.mir.bulge <- m$p3.target.bulge <- 3L
    do3p <- FALSE
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
                      substr(target,3,10+bulged+minBulge), UGsub=UGsub )
  # prevent wobble at first two positions:
  mm <- gsub("^(.)-","\\1-",gsub("^-"," ",mm))
  if(!flagBulgeMatches && m$p3.mir.bulge==m$p3.target.bulge){
    mm <- c(mm,rep(" ",m$p3.mir.bulge))
  }else if(m$p3.mir.bulge<m$p3.target.bulge){
    mm <- c(mm,rep(" ",m$p3.target.bulge-minBulge))
    di <- m$p3.target.bulge-m$p3.mir.bulge
    di1 <- paste(rep("-",ifelse(flagBulgeMatches,0,floor(di/2))),collapse="")
    di2 <- paste(rep("-",ifelse(flagBulgeMatches,di,ceiling(di/2))),collapse="")
    mirseq <- paste0(
      paste(substr(mirseq,1,8+bulged),collapse=""), di1,
      paste(substr(mirseq,9+bulged,8+bulged+m$p3.mir.bulge),collapse=""), di2,
      paste(substr(mirseq,9+bulged+m$p3.mir.bulge,nchar(mirseq)),collapse=""))
  }else{
    mm <- c(mm,rep(" ",m$p3.mir.bulge-minBulge))
    di <- m$p3.mir.bulge-m$p3.target.bulge
    target <- paste0(
      paste(substr(target,1,10),collapse=""),
      paste(rep("-",floor(di/2)),collapse=""),
      paste(substr(target,11,10+m$p3.target.bulge),collapse=""),
      paste(rep("-",ceiling(di/2)), collapse=""),
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
  if(hideSingletons) d$alignment <- gsub(" | ", paste(rep(" ",3),collapse=""),
                                         d$alignment, fixed=TRUE)
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
    mm[s1!=s2 & ((s1=="A" & s2=="G") | (s1=="C" & s2=="U"))] <- "-"
  }
  mm <- gsub("\\-$", " ", paste(mm, collapse=""))
  mm <- gsub("\\-[^|]|[^|]\\-", "  ", mm)
  mm
}

