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

#' getRandomSeq
#' 
#' Produces a random sequence of the given letters
#'
#' @param length Length of the sequence
#' @param alphabet Letters from which to sample
#'
#' @return A character vector of length 1
#' @export
#'
#' @examples
#' getRandomSeq()
getRandomSeq <- function(length=3000, alphabet=c("A","C","G","T"), xs=FALSE){
  x <- paste(sample(alphabet, size = length, replace=TRUE), collapse="")
  names(x) <- "seq1"
  x
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
           S$sites>=as.integer(minSites),]
  S[order(S$logp.binom),]
}


#' Create dummy log_kd per 12-mer data
#'
#' @param mod Optional model from which to create the dummy data
#'
#' @return A data.frame with 12-mers and log_kds
#' @export
#'
#' @examples
#' kd <- dummyKdData()
dummyKdData <- function(mod=NULL){
  if(is.null(mod)){
    data("SampleKdModel")
    mod <- SampleKdModel
  }
  mer12 <- paste0(getKmers(2),getSeed8mers(mod$canonical.seed),getKmers(2))
  data.frame(X12mer=mer12, log_kd=mod$mer8/1000)
}
