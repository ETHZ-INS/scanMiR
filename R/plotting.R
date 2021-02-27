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


viewTargetAlignment <- function(m, mod){
  stopifnot(is(m,"GRanges"))
  stopifnot(length(m)==1)
  if(all(c("sequence","p3.mir.bulge") %in% colnames(mcols(m))))
    return(.targetAlignment_internal(m, mod))
  stop("Not yet implemented")
}

.targetAlignment_internal <- function(m, mod){
  mirseq <- gsub("T","U",mod$mirseq)
  target <- stringi::stri_reverse(gsub("T","U",as.character(m$sequence)))
  mirseq2 <- as.character(complement(RNAString(mirseq)))
  mirseq2 <- paste0("A",substr(mirseq2, 2, nchar(mirseq2)))
  mm <- ifelse(strsplit(substr(mirseq2,1,8),"")[[1]]==
                 strsplit(substr(target,3,10),"")[[1]], "|", " ")
  cat(paste(paste0("--",mirseq),paste(c("  ",mm),collapse=""),target,sep="\n"))
}