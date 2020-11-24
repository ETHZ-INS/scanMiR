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
#' Generates all possible 8mers with 4 consecutive and positioned matches to a given
#' seed.
#'
#' @param seed The miRNA seed (target DNA sequence), a character vector of length 8 (if 
#' of length 7, a "A" will be added on the right)
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