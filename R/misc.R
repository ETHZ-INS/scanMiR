#' getKmers
#'
#' Returns all combinations of `n` elements of `from`
#'
#' @param n Number of elements
#' @param from Letters sampled
#'
#' @return A character vector
#' @export
#'
#' @examples
#' getKmers(3)
getKmers <- function(n=4, from=c("A", "C", "G", "T")){
  apply(expand.grid(lapply(seq_len(n), FUN=function(x) from)),
        1,collapse="",FUN=paste)
}


#' getSeed8mers
#' 
#' Generates all possible 8mers with 4 consecutive and positioned matches to a 
#' given seed.
#'
#' @param seed The miRNA seed (target DNA sequence), a character vector of 
#' length 8 (if of length 7, a "A" will be added on the right)
#' @param addNs Logical; whether to include 8mers with one flanking N
#'
#' @return A vector of 1024 8mers.
#' @export
#'
#' @examples
#' head(getSeed8mers("ACACTCCA"))
getSeed8mers <- function(seed, addNs=FALSE){
  a <- strsplit(as.character(seed),"")[[1]]
  if(length(a)==7) a <- c(a, "A")
  nts <- c("A","C","G","T")
  kmers <- expand.grid(nts,nts,nts,nts)
  y <- unique(unlist(lapply(0:4, FUN=function(x){
    paste0(
      do.call(paste0,kmers[,seq_len(x),drop=FALSE]),
      paste(a[x+1:4],collapse=""),
      do.call(paste0,kmers[,x+seq_len(4-x),drop=FALSE])
    )
  })))
  if(addNs){
    kmers <- expand.grid(nts,nts,nts)
    y <- c(y,unique(unlist(lapply(1:4, FUN=function(x){
      paste0("N",
             do.call(paste0,kmers[,seq_len(x-1),drop=FALSE]),
             paste(a[x+1:4],collapse=""),
             do.call(paste0,kmers[,x+seq_len(4-x)-1,drop=FALSE])
      )
    }))))
    y <- c(y,unique(unlist(lapply(0:3, FUN=function(x){
      paste0(do.call(paste0,kmers[,seq_len(x),drop=FALSE]),
             paste(a[x+1:4],collapse=""),
             do.call(paste0,kmers[,x+seq_len(3-x),drop=FALSE]),
             "N"
      )
    }))))
  }
  y
}

.build4mersRegEx <- function(seed){
  if(is(seed,"KdModel")){
    a <- strsplit(as.character(seed$canonical.seed),"")[[1]]
  }else{
    a <- strsplit(as.character(seed),"")[[1]]
  }
  if(length(a)==7) a <- c(a, "A")
  pats <- sapply(0:4, FUN=function(x){
    paste0(paste(rep("[^N]",x), collapse=""),
           paste(a[x+1:4],collapse=""),
           paste(rep(".",4-x), collapse=""))
  })
  paste(pats,collapse="|")
}


.guessSeqType <- function(x, use.subset=TRUE){
  seqs <- x[sample.int(length(x),min(length(x),10))]
  u <- any(grepl("U",seqs))
  t <- any(grepl("T",seqs))
  if(t && u) stop("Sequences contain both T and U!")
  if(t || u) return(ifelse(u,"RNA","DNA"))
  if(length(seqs)>10) return(.guessSeqType(x, FALSE))
  warning("Sequences contain neither T or U - assuming they are DNA...")
  return("DNA")
}


#' @import Matrix
.matches2sparse <- function(x){
  as(as.matrix(table(as.factor(seqnames(x)), as.factor(x$miRNA))), 
     "sparseMatrix")
}
.sparse2df <- function(x, content="value"){
  dimn <- dimnames(x)
  x <- summary(x)
  w <- which(x$x!=ifelse(is.logical(x$x),FALSE,0))
  xout <- data.frame( feature=factor(x$i[w], levels=seq_len(length(dimn[[1]])), labels=dimn[[1]]),
                      set=factor(x$j[w], levels=seq_len(length(dimn[[2]])), labels=dimn[[2]]) )
  if(!is.logical(x$x)) xout[[content]] <- x$x[w]
  xout
}

#' enrichedMirTxPairs
#'
#' Identifies pairs of miRNA and target transcripts that have an unexpectedly
#' high number of sites.
#'
#' @param m A GRanges of matches, as produced by \link{\code{findSeedMatches}}.
#' It is recommended to filter this to have only canonical sites.
#' @param minSites The minimum number of sites for a given miRNA-target pair to
#' be considered.
#' @param max.binom.p The maximum binomial p-value of miRNA-target pairs.
#'
#' @return A data.frame of top combinations, including number of sites and 
#' the log-transformed binomial p-value.
#' @export
#'
#' @import Matrix
#' @importFrom stats pbinom
enrichedMirTxPairs <- function(m, minSites=5, max.binom.p=0.001){
  m <- m[as.integer(m$type) %in% grep("8mer|7mer",.matchLevels())]
  b <- .matches2sparse(m)
  b <- b[rowSums(b>=minSites)>0,]
  rs <- rowSums(b)
  cs <- colSums(b)
  p <- as.matrix(rs/sum(rs)) %*% t(cs/sum(cs))
  S <- pbinom(as.matrix(b)-1, prob=p, size=rs, lower.tail=FALSE, log.p=TRUE)
  mode(S) <- "integer"
  S <- as(round(S), "sparseMatrix")
  S <- .sparse2df(S, "logp.binom")
  b <- .sparse2df(b, "sites")
  S$sites <- b$sites
  rm(b)
  S <- S[S$logp.binom < as.integer(log(max.binom.p)) & 
           Sd$sites>=as.integer(minSites),]
  S[order(S$logp.binom),]
}


#' getTranscriptSequence
#' 
#' Utility wrapper to extracts the sequence of a given transcript (UTR or 
#' CDS+UTR).
#'
#' @param tx The ensembl ID of the transcript
#' @param species The species, either 'hsa', 'rno', or 'mmu'; ignored if 
#' `ensdbs` is given and is an `EnsDb` object.
#' @param ensdbs An \code{\link[ensembldb]{EnsDb}} object, or a named list of 
#' such objects.
#' @param genome The genome sequence (e.g. 
#' \code{\link[BSgenome]{BSgenome-class}}). If one of the three pre-defined 
#' `species`, will attempt to fetch the corresponding packages (if isntalled).
#' @param UTRonly Logical; whether to fetch only the UTR sequences.
#' @param ... Passed to \code{\link{AnnotationHub}}
#'
#' @return A \link{\code{DNAStringSet}}.
#' @export
#'
#' @importFrom GenomicFeatures threeUTRsByTranscript cdsBy extractTranscriptSeqs
#' @importFrom ensembldb metadata organism seqlevelsStyle
#' @importFrom AnnotationHub AnnotationHub query
#' @import Biostrings
#' @examples
#' # not run:
#' # getTranscriptSequence("ENST00000641515", species="hsa")
getTranscriptSequence <- function(tx, species=NULL, ensdbs=NULL, genome=NULL,
                                  UTRonly=TRUE, ...){
  if(is.null(ensdbs)){
    species <- match.arg(species, c("hsa","rno","mmu"))
    ah <- AnnotationHub(...)
    ahid <- switch(species,
      hsa=rev(query(ah, c("EnsDb", "Homo sapiens"))$ah_id)[1],
      mmu=rev(query(ah, c("EnsDb", "Mus musculus"))$ah_id)[1],
      rno=rev(query(ah, c("EnsDb", "Rattus norvegicus"))$ah_id)[1],
      stop("Species not among the pre-defined one, please provide `ensdbs` ",
           "and `genome` manually.")
    )
    ensdb <- ah[[ahid]]
    em <- metadata(ensdb)
    em <- setNames(em$value, em$name)
    message("Using ", em[["genome_build"]], ", Ensembl version ",
            em[["ensembl_version"]])
  }else{
    if(is.null(species) && length(ensdbs)==1){
      if(is.list(ensdbs)){
        stopifnot(!is.null(names(ensdbs)))
        species <- names(ensdbs)
        ensdb <- ensdbs[[1]]
      }else if(is(ensdbs,"EnsDb")){
        species <- organism(ensdbs)
        ensdb <- ensdbs
      }else{
        stop("`ensdbs` should either by a object of class EnsDb or a named ",
             "list of such objects.")
      }
    }else{
      species <- match.arg(species, names(ensdbs))
      ensdb <- ensdbs[[species]]
    }
  }
  if(is.null(genome)){
    em <- metadata(ensdb)
    em <- setNames(em$value, em$name)
    gbuild <- em[["genome_build"]]
    genome <- switch(gbuild,
      GRCh38=BSgenome.Hsapiens.UCSC.hg38:::BSgenome.Hsapiens.UCSC.hg38,
      GRCm38=BSgenome.Mmusculus.UCSC.mm10:::BSgenome.Mmusculus.UCSC.mm10,
      "Rnor_6.0"=BSgenome.Rnorvegicus.UCSC.rn6:::BSgenome.Rnorvegicus.UCSC.rn6,
      stop("Genome not among the pre-defined one, please provide `genome` ",
           "manually.")
    )
  }
  seqlevelsStyle(genome) <- "Ensembl"

  gr <- suppressWarnings(threeUTRsByTranscript(ensdb, filter=~tx_id==tx))
  gr <- gr[seqnames(gr) %in% seqlevels(genome)]
  if(length(gr)==0) stop("Transcript not found!")
  seqs <- extractTranscriptSeqs(genome, gr)
  if(!UTRonly){
    grl_ORF <- cdsBy(ensdb, by="tx", filter=~tx_id==tx)
    seqs_ORF <- extractTranscriptSeqs(genome, grl_ORF)
    orf.len <- setNames(lengths(seqs_ORF), names(seqs_ORF))
    seqs_ORF[names(seqs)] <- xscat(seqs_ORF[names(seqs)],seqs)
    seqs <- seqs_ORF
    rm(seqs_ORF)
    mcols(seqs)$ORF.length <- orf.len[names(seqs)]
  }
  seqs
}


.mirTargetAlignment <- function(mirseq, target){
  mirseq <- paste0("A",substr(gsub("T","U",mirseq),2,nchar(mirseq)))
  target <- RNAString(gsub("T","U",target))
  mirseq <- reverseComplement(RNAString(mirseq))
  al <- pairwiseAlignment(mirseq, target, type="local-global")
}
