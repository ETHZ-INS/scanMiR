#' findSeedMatches
#'
#' @param seqs A character vector or `XStringSet`` of sequences in which to look.
#' @param seeds A character vector of 7-nt seeds to look for. If RNA, will be 
#' reversed and complemented before matching. If DNA, they are assumed to be
#' the target sequence to look for. Alternatively, a list of objects of class
#' `KdModel` or an object of class `CompressedKdModelList` can be given.
#' @param seedtype Either RNA, DNA or 'auto' (default)
#' @param shadow Integer giving the shadow, i.e. the number of nucleotides
#'  hidden at the beginning of the sequence (default 0)
#' @param minLogKd Minimum log_kd value to keep (default 1). Set to 0 to disable.
#' @param types The type of sites to return (default all)
#' @param max.noncanonical.motifs The maximum number of non-canonical motifs to search 
#' for (default all)
#' @param keepMatchSeq Logical; whether to keep the sequence (including flanking
#' dinucleotides) for each seed match (default FALSE).
#' @param minDist Integer specifying the minimum distance between matches of the same 
#' miRNA (default 1). Closer matches will be reduced to the highest-affinity. To 
#' disable the removal of overlapping features, use `minDist=-Inf`.
#' @param fastRemoveOverlaps Whether to use the fast (rather than exact) method to remove 
#' overlaps. Will improve speed, but some matching sites near stronger ones (but beyond 
#' minDist) might be lost when more than two sites overlap. Not recommended.
#' @param BP Pass `BiocParallel::MulticoreParam(ncores, progressbar=TRUE)` to enable 
#' multithreading.
#' @param verbose Logical; whether to print additional progress messages (default on if 
#' not multithreading)
#'
#' @return A GRanges of all matches
#' 
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @import BiocParallel Biostrings GenomicRanges
#' @export
#'
#' @examples
#' # we create mock RNA sequences and seeds:
#' seqs <- sapply(1:10, FUN=function(x) paste(sample(strsplit("ACGU", "")[[1]], 
#'                                      1000, replace=TRUE),collapse=""))
#' names(seqs) <- paste0("seq",1:length(seqs))
#' seeds <- c("AAACCAC", "AAACCUU")
#' findSeedMatches(seqs, seeds)
findSeedMatches <- function( seqs, seeds, seedtype=c("auto", "RNA","DNA"), shadow=0L, 
                             minLogKd=1, keepMatchSeq=FALSE, maxLoop=10L, mir3p.nts=8L, 
                             minDist=7L, fastRemoveOverlaps=FALSE, BP=NULL, verbose=NULL){
  
  if(is.null(verbose)) verbose <- is(seeds,"KdModel") || length(seeds)==1 || is.null(BP)
  if(verbose) message("Preparing sequences...")
  args <- .prepSeqs(seqs, seeds, seedtype, shadow=shadow, pad=c(maxLoop+mir3p.nts+6L,6L))
  seqs <- args$seqs
  if("seeds" %in% names(args)) seeds <- args$seeds
  offset <- args$offset
  rm(args)

  if(is(seeds,"KdModel") || length(seeds)==1){
    if(is.list(seeds[[1]])) seeds <- seeds[[1]]
    if(is.null(verbose)) verbose <- TRUE
    m <- .find1SeedMatches(seqs, seeds, keepMatchSeq=keepMatchSeq, minDist=minDist, 
                           minLogKd=minLogKd, maxLoop=maxLoop, mir3p.nts=mir3p.nts,
                           fastRemoveOverlaps=fastRemoveOverlaps, verbose=verbose)
    if(length(m)==0) return(m)
  }else{
    if(is.null(BP)) BP <- SerialParam()
    if(is.null(verbose)) verbose <- !(bpnworkers(BP)>1 | length(seeds)>5)
    m <- bplapply( seeds, seqs=seqs, keepMatchSeq=keepMatchSeq, verbose=verbose, 
                   minDist=minDist, minLogKd=minLogKd, maxLoop=maxLoop, 
                   mir3p.nts=mir3p.nts, fastRemoveOverlaps=fastRemoveOverlaps,
                   BPPARAM=BP, FUN=.find1SeedMatches)
    m <- GRangesList(m)
    if(is.null(names(m))){
      if(!is.character(seeds)) seeds <- sapply(seeds, FUN=function(x){
        if(is.null(x$name)) return(x$canonical.seed)
        x$name
      })
      names(m) <- seeds
    }
    mirs <- Rle(as.factor(names(m)),lengths(m))
    m <- unlist(m)
    m$miRNA <- mirs
    m
  }

  gc(verbose = FALSE, full = TRUE)
  
  names(m) <- row.names(m) <- NULL
  m <- shift(m, offset)

  metadata(m)$call.params <- list(
    shadow=shadow,
    minDist=minDist,
    maxLoop=maxLoop,
    mir3p.nts=mir3p.nts
  )
  m
}


#' removeOverlappingMatches
#' 
#' Removes elements from a GRanges that overlap (or are within a given distance of) other 
#' elements higher up in the list (i.e. assumes that the ranges are sorted in order of
#' priority). The function handles overlaps between more than two ranges by successively
#' removing those that overlap higher-priority ones
#'
#' @param x A GRanges, sorted by (decreasing) importance
#' @param minDist Minimum distance between ranges
#'
#' @return A filtered GRanges.
#' @export
removeOverlappingMatches <- function(x, minDist=7L, method=c("exact","fast")){
  method <- match.arg(method)
  red <- GenomicRanges::reduce(x, with.revmap=TRUE, min.gapwidth=minDist)$revmap
  red <- red[lengths(red)>1]
  if(length(red)==0) return(x)
  if(method=="fast"){
    toRemove <- setdiff(unlist(red),min(red))
  }else{
    red2 <- red[lengths(red)==2]
    toRemove <- setdiff(unlist(red2), min(red2))
    toRemove <- c(toRemove, unlist(lapply( red[lengths(red)>2], FUN=function(i){
      w <- j <- i <- sort(i)
      y <- ranges(x)[i]
      while(length(y)>1 && length(w)>0){
        # remove anything overlappnig the top one
        w <- 1+which(overlapsAny(y[-1],y,maxgap=minDist))
        if(length(w)>0){
          j <- j[-w]
          y <- y[-w]
        }
      }
      return(setdiff(i,j))
    })))
  }
  if(length(toRemove)>0) x <- x[-toRemove]
  x
}


.prepSeqs <- function(seqs, seeds, seedtype=c("auto", "RNA","DNA"), shadow=0, pad=c(0,0)){
  if(is.null(names(seqs))) names(seqs) <- paste0("seq",seq_along(seqs))
  seedtype <- match.arg(seedtype)
  seqtype <- .guessSeqType(seqs)
  ret <- list()
  if( is(seeds, "KdModel") || 
      (is.list(seeds) && all(sapply(seeds, is.list))) ){
    if(is.null(names(seeds)))
      stop("If `seeds` is a list of kd models, it should be named.")
    if(seedtype=="RNA" || seqtype=="RNA") 
      stop("If `seeds` is a list of kd models, both the seeds and the target
sequences should be in DNA format.")
  }else{
    if(is.null(names(seeds))) n <- names(seeds) <- seeds
    if(seedtype=="auto") seedtype <- .guessSeqType(seeds)
    if(seedtype=="RNA"){
      message("Matching reverse complements of the seeds...")
      seeds <- as.character(reverseComplement(RNAStringSet(seeds)))
    }else{
      message("Matching the given seeds directly...")
    }
    if(seqtype=="RNA"){
      seeds <- gsub("T", "U", seeds)
    }else{
      seeds <- gsub("U", "T", seeds)
    }
    names(seeds) <- n
    ret$seeds <- seeds
  }
  seqnms <- names(seqs)
  if(is.character(seqs)) seqs <- DNAStringSet(seqs)
  names(seqs) <- seqnms
  seqs <- seqs[lengths(seqs)>=(shadow+7)]
  if(shadow>0){
    seqs <- subseq(seqs, shadow+1-min(c(shadow,pad[1],2)))
    pad[1] <- pad[1]-min(shadow,pad[1])
  }
  ret$offset <- shadow-pad[1]
  seqs <- padAndClip(seqs, views=IRanges(start=1-pad[1], width=lengths(seqs)+pad[1]+pad[2]), 
                     Lpadding.letter = "N", Rpadding.letter = "N")
  c(ret, list(seqs=seqs))
}


.find1SeedMatches <- function(seqs, seed, keepMatchSeq=FALSE, minLogKd=0, maxLoop=10, 
                              mir3p.nts=8, minDist=1, fastRemoveOverlaps=FALSE, 
                              verbose=FALSE){
  if(verbose) message("Scanning for matches...")

  if(isPureSeed <- is.character(seed)){
    pos <- gregexpr(paste0("(?=.",substr(seed,2,7),".)"), seqs, perl=TRUE)
  }else{
    mod <- seed
    seed <- mod$canonical.seed
    pos <- gregexpr(paste0("(?=",.build4mersRegEx(seed),")"), seqs, perl=TRUE)
  }
  pos <- lapply(lapply(pos, as.numeric), y=-1, setdiff)
  if(sum(lengths(pos))==0){
    if(verbose) message("Nothing found!")
    return(GRanges())
  }
  m <- GRanges( rep(names(seqs), lengths(pos)), IRanges( start=unlist(pos), width=8 ) )
  m <- keepSeqlevels(m, seqlevelsInUse(m))
  m <- m[order(seqnames(m))]
  
  if(verbose) message("Exctracting sequences and characterizing matches...")
  seqs <- seqs[seqlevels(m)]
  r <- ranges(m)

  if(isPureSeed){
    r <- split(r, seqnames(m))
    names(r) <- NULL
    ms <- as.factor(unlist(extractAt(seqs, r)))
    if(keepMatchSeq) mcols(m)$sequence <- ms
    mcols(m)$type <- characterizeSeedMatches(ms, substr(seed,1,7))[,"type"]
  }else{
    start(r) <- start(r)-1-maxLoop-mir3p.nts
    end(r) <- end(r)+2
    r <- split(r, seqnames(m))
    names(r) <- NULL
    ms <- unlist(extractAt(seqs, r))
    names(ms) <- NULL
    mir.3p <- as.character(reverseComplement(DNAString(
      substr(mod$mirseq, 12, min(c(11+mir3p.nts, nchar(mod$mirseq)))) )))
    al <- pairwiseAlignment(subseq(ms,1,maxLoop+mir3p.nts), mir.3p, 
                            type="local", scoreOnly=TRUE)
    mcols(m)$align.3p <- as.integer(round(1000*al))
    rm(al)
    ms <- subseq(ms, maxLoop+mir3p.nts, 11+maxLoop+mir3p.nts)
    if(keepMatchSeq) mcols(m)$sequence <- as.factor(ms)
    fl.score <- as.numeric(.getFlankingScore(ms)$score)
    mer8 <- factor(as.character(subseq(ms, 3,10)), levels=getSeed8mers(seed))
    type <- .getMatchTypes(levels(mer8), seed)
    mer8 <- as.integer(mer8)
    mcols(m)$type <- type[mer8]
    mcols(m)$log_kd <- as.integer(round(mod$mer8[mer8] + fl.score*mod$fl[mer8]))
  }
  rm(ms)
  
  if(!is.null(minLogKd) && "log_kd" %in% colnames(mcols(m))){
    if(minLogKd>0) minLogKd <- -minLogKd
    if(minLogKd < -10) minLogKd <- minLogKd*1000
    m <- m[which(m$log_kd <= as.integer(round(minLogKd)))]
  }
  
  if(minDist>-Inf){
    if(verbose) message("Removing overlaps...")
    m <- removeOverlappingMatches(m, minDist=minDist, 
                                  method=ifelse(fastRemoveOverlaps,"fast","exact"))
  }
  
  names(m) <- NULL
  m
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

#' characterizeSeedMatches
#'
#' @param x A factor or character vector of matches, or a data.frame or GRanges
#' containing a `sequence` column.
#' @param seed The seed which sequences of `x` should match, or a 'KdModel' object.
#'
#' @return A data.frame, or a `GRanges` if `x` is a `GRanges`
#' @export
#'
#' @examples
#' characterizeSeedMatches(c("UAAACCACCC","CGAACCACUG"), "AAACCAC")
characterizeSeedMatches <- function(x, seed=NULL){
  if(is(seed, "KdModel")){
    kd.model <- seed
    sseed <- seed$canonical.seed
  }else{
    kd.model <- NULL
    sseed <- seed
  }
  sseed <- as.character(sseed)
  if(length(sseed)>1 || nchar(sseed)!=7) 
    stop("`seed` must contain a single string of 7 characters.")
  if(!is.character(x) && !is.factor(x) && !is(x,"Rle")){
    y <- characterizeSeedMatches(x$sequence, seed)
    if(is(x,"GRanges")){
      mcols(x) <- cbind(mcols(x), y)
    }else{
      x <- cbind(x,y)
    }
    return(x)
  }
  if(is(x,"Rle") || any(duplicated(x))){
    # we resolve type for unique sequences
    if(!is.null(levels(x))){
      levels(x) <- gsub("N","x",levels(x))
      y <- characterizeSeedMatches(levels(x), seed)
    }else{
      x <- gsub("N","x",x)
      y <- characterizeSeedMatches(unique(x), seed)
    }
    y <- y[as.character(x),,drop=FALSE]
    row.names(y) <- NULL
    return(y)
  }
  d <- data.frame( row.names=x, type=.getMatchTypes(as.character(x), seed=sseed) )
  if(!is.null(kd.model)) d$log_kd <- predictKD(row.names(d), kd.model)
  d
}

.getMatchTypes <- function(x, seed){
  x <- as.character(x)
  y <- rep(1L,length(x))
  if(nchar(seed)==7) seed <- paste0(seed,"A")
  seed6 <- substr(seed,2,7)
  y[grep(seed6,x,fixed=TRUE)] <- 2L # offset 6mer
  y[grep(paste0("[ACGT]","[ACGT]",substr(seed,3,8)),x)] <- 3L # 6mer-a1
  y[grep(paste0(substr(seed,1,6),"[ACGT][ACGT]"),x)] <- 4L # 6mer-m8
  y[grep(paste0("[ACGT]",substr(seed,2,7)),x)] <- 5L # 6mer
  y[grep(paste0("[ACGT]",substr(seed,2,8)),x)] <- 6L # 7mer-a1
  y[grep(substr(seed,1,7),x,fixed=TRUE)] <- 7L # 7mer-m8
  y[grep(seed,x,fixed=TRUE)] <- 8L # 8mer
  factor(y, levels=8:1, labels=c("8mer","7mer-m8","7mer-a1","6mer","6mer-m8",
                                 "6mer-a1","offset 6mer","non-canonical"))
}

# deprecated, to be removed
.getMatchType <- function(x, seed){
  if(grepl(paste0(seed,"A"),x,fixed=TRUE)) return("8mer")
  if(grepl(seed,x,fixed=TRUE)) return("7mer-m8")
  seed6 <- substr(seed,2,7)
  if(grepl(paste0("[ACGT]",seed6,"A"),x)) return("7mer-a1")
  if(grepl(paste0("[ACGT]",seed6),x)) return("6mer")
  seed5m <- substr(seed,1,6)
  if(grepl(paste0(seed5m,"[ACGT]","[ACGT]"),x)) return("6mer-m8")
  seed5 <- substr(seed,3,7)
  if(grepl(paste0("[ACGT]","[ACGT]",seed5,"A"),x)) return("6mer-a1")
  if(grepl(seed6,x,fixed=TRUE)) return("offset 6mer")
  "non-canonical"
}

.datatable.aware = TRUE

#' aggregateMatches
#'
#' @param e A GRanges object as produced by `findSeedMatches`.
#'
#' @return An aggregated data.frame
#' @importFrom data.table data.table as.data.table dcast
#' @importFrom GenomicRanges mcols
#' @export
aggregateMatches <- function(e, fn=agg.repr){
  d <- as.data.frame(mcols(e))
  d$transcript <- as.factor(seqnames(e))
  d <- as.data.table(d)
  d2 <- subset(d, type!="non-canonical")
  ag1a <- d2[,.( log_kd.canonical=log10(1/sum(1/10^log_kd)), repr.canonical=fn(log_kd)),
                 by=c("transcript","seed")]
  ag2 <- d[,.( log_kd=log10(1/sum(1/10^log_kd)), repr=fn(log_kd)),
             by=c("transcript","seed")]  
  ag1b <- dcast( d2[,.(N=.N), by=c("transcript","seed","type")],
                 formula=transcript+seed~type, value.var="N", fill=0)
  rm(d,d2)
  m <- merge(ag1b, ag1a, by=c("transcript","seed"), all=TRUE)
  m <- merge(m, ag2, by=c("transcript","seed"), all=TRUE)
  m <- as.data.frame(m)
  for(f in c("log_kd.canonical", "log_kd", "repr.canonical", "repr")){
    m[[f]][is.na(m[[f]])] <- 0
  }
  m$`offset 6mer` <- NULL
  for( f in c("8mer","7mer-m8","7mer-A1","6mer") ) m[[f]][is.na(m[[f]])] <- 0
  colnames(m)[3:6] <- gsub("-","",paste0("n.",colnames(m)[3:6]))
  m
}

agg.repr <- function(x, b=1.8, ag=10^-2){
  log(1+b*sum(ag/(ag+rep(1,length(x)))))-log(1+b*sum(ag/(ag+10^x)))
}
