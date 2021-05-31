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
        1, collapse="",FUN=paste)
}

#' getRandomSeq
#'
#' Produces a random sequence of the given letters
#'
#' @param length Length of the sequence
#' @param alphabet Letters from which to sample
#' @param n The number of sequences to generate
#'
#' @return A character vector of length 1
#' @export
#'
#' @examples
#' getRandomSeq(100)
getRandomSeq <- function(length=3000, alphabet=c("A","C","G","T"), n=1){
  seqs <- vapply(seq_len(n), FUN.VALUE=character(1), FUN=function(i){
    paste(sample(alphabet, size = length, replace=TRUE), collapse="")
  })
  names(seqs) <- paste0("seq",seq_len(n))
  seqs
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
  pats <- vapply(0:4, FUN.VALUE=character(1), FUN=function(x){
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

#' Create dummy log_kd per 12-mer data
#'
#' @param mod Optional model from which to create the dummy data
#'
#' @return A data.frame with 12-mers and log_kds
#' @export
#' @importFrom utils data
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

                    
#' get8merRange
#'
#' Returns the minimum and maximum 8-mer log-kd values
#'
#' @param mod A `KdModel`
#'
#' @return A numeric vector of length two
#' @export
#' @examples
#' data("SampleKdModel")
#' get8merRange(SampleKdModel)
get8merRange <- function(mod){
  stopifnot(is(mod,"KdModel"))
  mer8 <- getSeed8mers(mod$canonical.seed)
  mer8 <- which(mer8==mod$canonical.seed)
  fl <- rowSums(apply(.flankingValues(), 2, FUN=range))
  fl*mod$fl[mer8]+mod$mer8[mer8]
}
