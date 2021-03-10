#' findSeedMatches
#'
#' @param seqs A character vector or `XStringSet` of sequences in which to look.
#' @param seeds A character vector of 7-nt seeds to look for. If RNA, will be 
#' reversed and complemented before matching. If DNA, they are assumed to be
#' the target sequence to look for. Alternatively, a list of objects of class
#' `KdModel` or an object of class `KdModelList` can be given.
#' @param seedtype Either RNA, DNA or 'auto' (default)
#' @param shadow Integer giving the shadow, i.e. the number of nucleotides
#'  hidden at the beginning of the sequence (default 0)
#' @param maxLogKd Maximum log_kd value to keep (default 0). Set to Inf to disable.
#' @param keepMatchSeq Logical; whether to keep the sequence (including flanking
#' dinucleotides) for each seed match (default FALSE).
#' @param onlyCanonical Logical; whether to restrict the search only to canonical
#' binding sites.
#' @param minDist Integer specifying the minimum distance between matches of the same 
#' miRNA (default 1). Closer matches will be reduced to the highest-affinity. To 
#' disable the removal of overlapping features, use `minDist=-Inf`.
#' @param p3.maxLoop The maximum loop size for the 3' alignment
#' @param p3.mismatch Logical; whether to allow mismatches in 3' alignment
#' @param p3.params Named list of parameters for the 3' alignment.
#' @param agg.params a named list with slots `ag`, `b` and `c` indicating the 
#' parameters for the aggregation. Ignored if `ret!="aggregated"`.
#' @param ret The type of data to return, either "GRanges" (default), 
#' "data.frame" (lighter weight than GRanges), or "aggregated" (aggregated per 
#' transcript).
#' @param BP Pass `BiocParallel::MulticoreParam(ncores, progressbar=TRUE)` to enable 
#' multithreading.
#' @param verbose Logical; whether to print additional progress messages (default on if 
#' not multithreading)
#'
#' @return A GRanges of all matches. If `seeds` is a `KdModel` or `KdModelList`, the 
#' `log_kd` column will report the ln(Kd) multiplied by 1000, rounded and saved as an 
#' integer.
#' 
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @import Biostrings GenomicRanges
#' @export
#'
#' @examples
#' # we create mock RNA sequences and seeds:
#' seqs <- sapply(1:10, FUN=function(x) paste(sample(strsplit("ACGU", "")[[1]], 
#'                                      1000, replace=TRUE),collapse=""))
#' names(seqs) <- paste0("seq",1:length(seqs))
#' seeds <- c("AAACCAC", "AAACCUU")
#' findSeedMatches(seqs, seeds)
findSeedMatches <- function( seqs, seeds, shadow=0L, onlyCanonical=FALSE, 
                             maxLogKd=c(-0.3,-0.3), keepMatchSeq=FALSE, 
                             minDist=7L, p3.extra=FALSE, 
                             p3.params=list(maxMirLoop=5L, maxTargetLoop=9L, 
                                            maxLoopDiff=4L, mismatch=TRUE),
                             agg.params=.defaultAggParams(),
                             ret=c("GRanges","data.frame","aggregated"), 
                             BP=NULL, verbose=NULL, ...){
  p3.params <- .check3pParams(p3.params)
  # This might not be most efficient:
  length.seqs <- width(seqs)
  if(is.character(seqs) || is.null(mcols(seqs)$ORF.length)){
    utr_len <- length.seqs
    orf_len <- rep(0, length.out = length(utr_len))
  }else{
    orf_len <- mcols(seqs)[,"ORF.length"] - 15
    utr_len <- length.seqs - orf_len
  }
  length.info <- cbind(orf_len,utr_len)
  if(!is.null(names(seqs))) row.names(length.info) <- names(seqs)
  ###
  ret <- match.arg(ret)
  if(ret=="aggregated"){
    if(!is.list(agg.params)) agg.params <- as.list(agg.params)
    if(!all(c("ag","b","c","p3","coef_utr","coef_orf") %in% names(agg.params)))
      stop("`agg.params` should be a named list with slots ",
           "`ag`, `b`, `c`, `p3`, `coef_utr` and `coef_orf`.")
  }
  seedInputType <- .checkSeedsInput(seeds)
  if(is.null(verbose)) verbose <- is(seeds,"KdModel") || length(seeds)==1 || is.null(BP)
  if(verbose) message("Preparing sequences...")
  maxLoop <- max(unlist(p3.params[c("maxMirLoop","maxTargetLoop")]))
  args <- .prepSeqs(seqs, seeds, shadow=shadow, pad=c(maxLoop+30L,8L))
  seqs <- args$seqs
  if("seeds" %in% names(args)) seeds <- args$seeds
  offset <- args$offset
  rm(args)

  params <- list(
    shadow=shadow,
    minDist=minDist,
    maxLogKd=maxLogKd,
    p3.params=p3.params
  )
  if(ret=="aggregated") params$agg.params <- agg.params
  
  if(is(seeds,"KdModel") || length(seeds)==1){
    if(is.list(seeds[[1]])) seeds <- seeds[[1]]
    params$miRNA <- ifelse(is(seeds, "KdModel"), seeds$name, seeds)
    if(is.null(verbose)) verbose <- TRUE
    m <- .find1SeedMatches(seqs, seeds, keepMatchSeq=keepMatchSeq, 
                           minDist=minDist, maxLogKd=maxLogKd, 
                           onlyCanonical=onlyCanonical, p3.extra=p3.extra,
                           p3.params=p3.params, offset=offset, 
                           verbose=verbose, ret=ret, ...)
    if(length(m)==0) return(m)
    if(ret=="aggregated"){
      if(verbose) message("Aggregating...")
      ll <- as.data.frame(length.info)
      ll$transcript <- row.names(ll)
      m <- .aggregate_miRNA(m,ll, ag=agg.params$ag, b=agg.params$b, 
                            c=agg.params$c, p3 = agg.params$p3, coef_utr = agg.params$coef_utr,
                            coef_orf = agg.params$coef_orf, toInt=TRUE)
    }
  }else{
    if(is.null(BP)) BP <- SerialParam()
    if(is.null(verbose)) verbose <- !(bpnworkers(BP)>1 | length(seeds)>5)
    m <- bplapply( seeds, BPPARAM=BP, FUN=function(oneseed){
      m <- .find1SeedMatches(seqs=seqs, seed=oneseed, keepMatchSeq=keepMatchSeq,
                 minDist=minDist, maxLogKd=maxLogKd, p3.extra=p3.extra,
                 onlyCanonical=onlyCanonical, p3.params=p3.params, ret=ret, 
                 offset=offset, verbose=verbose, ...)
      if(length(m)==0) return(m)
      if(ret=="aggregated"){
        if(verbose) message("Aggregating...")
        ll <- as.data.frame(length.info)
        ll$transcript <- row.names(ll)
        m <- .aggregate_miRNA(m,ll, ag=agg.params$ag, b=agg.params$b, 
                              c=agg.params$c,p3 = agg.params$p3, coef_utr = agg.params$coef_utr,
                              coef_orf = agg.params$coef_orf, toInt=TRUE)
      }
      m
    } )
      
    if(ret=="GRanges") m <- GRangesList(m)
    
    if(is.null(names(m))){
      if(!is.character(seeds)) seeds <- sapply(seeds, FUN=function(x){
        if(is.null(x$name)) return(x$canonical.seed)
        x$name
      })
      names(m) <- seeds
    }
    
    if(ret=="GRanges"){
      mirs <- Rle(as.factor(names(m)),lengths(m))
      m <- unlist(m)
      metadata(m)$call.params <- params
      metadata(m)$length.info <- length.info
      names(m) <- row.names(m) <- NULL
      m$miRNA <- mirs
    }else{
      mirs <- rep(as.factor(names(m)),lengths(m))
      m <- dplyr::bind_rows(m, .id="miRNA")
      m$miRNA <- as.factor(m$miRNA)
      attr(m, "call.params") <- params
      row.names(m) <- NULL
      m$transcript <- as.factor(m$transcript)
    }
  }

  gc(verbose = FALSE, full = TRUE)

  m
}

# scan for a single seed
.find1SeedMatches <- function(seqs, seed, keepMatchSeq=FALSE, maxLogKd=0, 
                              minDist=1L, onlyCanonical=FALSE, p3.extra=FALSE,
                              p3.params=list(), offset=0L, 
                              ret=c("GRanges","data.frame","aggregated"), 
                              verbose=FALSE){
  ret <- match.arg(ret)
  p3.params <- .check3pParams(p3.params)

  if(is.null(maxLogKd)) maxLogKd <- c(Inf,Inf)
  if(length(maxLogKd)==1) maxLogKd <- rep(maxLogKd,2)
  
  if(verbose) message("Scanning for matches...")
  
  if(isPureSeed <- is.character(seed)){
    stopifnot(nchar(seed)>6)
    if(nchar(seed) %in% 7:8){
      mirseq <- NULL
    }else{
      mirseq <- gsub("U","T",seed)
      seed <- as.character(reverseComplement(DNAStringSet(substr(mirseq,2,8))))
      seed <- paste0(seed,"A")
      if(verbose && (nchar(mirseq)<18 | nchar(mirseq)>24))
        warning("The `seed` given seems to be neither a miRNA target seed nor",
                "a mature miRNA sequence! Scanning for ", seed)
    }
    pat <- substr(seed,2,7)
    if(!onlyCanonical)
      pat <- paste0(pat,"|",substr(seed,2,3),"G",substr(seed,4,7))
    pos <- gregexpr(paste0("(?=.",substr(seed,2,7),".)"), seqs, perl=TRUE)
  }else{
    mod <- seed
    seed <- mod$canonical.seed
    mirseq <- mod$mirseq
    if(onlyCanonical){
      patt <- paste0(".",substr(seed,2,7),".")
    }else{
      patt <- .build4mersRegEx(seed)
    }
    pos <- gregexpr(paste0("(?=",patt,")"), seqs, perl=TRUE)
  }
  pos <- lapply(lapply(pos, as.integer), y=-1L, setdiff)
  if(sum(lengths(pos))==0){
    if(verbose) message("Nothing found!")
    return(GRanges())
  }
  m <- GRanges( rep(names(seqs), lengths(pos)), IRanges( start=unlist(pos), width=8 ) )
  m <- keepSeqlevels(m, seqlevelsInUse(m))
  m <- m[order(seqnames(m))]
  
  if(verbose) message("Extracting sequences and characterizing matches...")
  seqs <- seqs[seqlevels(m)]
  r <- ranges(m)
  
  if(isPureSeed && is.null(mirseq)){
    r <- split(r, seqnames(m))
    names(r) <- NULL
    ms <- as.factor(unlist(extractAt(seqs, r)))
    if(keepMatchSeq) mcols(m)$sequence <- ms
    mcols(m)$type <- getMatchTypes(ms, substr(seed,1,7))
    m <- m[order(seqnames(m), m$type)]
  }else{
    maxLoop <- max(unlist(p3.params[c("maxMirLoop","maxTargetLoop")]))
    plen <- maxLoop+nchar(mirseq)-8L
    start(r) <- start(r)-1L-plen
    end(r) <- end(r)+2L
    r <- split(r, seqnames(m))
    names(r) <- NULL
    ms <- unlist(extractAt(seqs, r))
    names(ms) <- NULL
    p3 <- get3pAlignment( subseq(ms,1L,1+plen), mirseq, 
                          allow.mismatch=p3.params$mismatch,
                          maxMirLoop=p3.params$maxMirLoop, 
                          maxLoopDiff=p3.params$maxLoopDiff,
                          maxTargetLoop=p3.params$maxTargetLoop)
    if(p3.extra){
      mcols(m) <- cbind(mcols(m), p3)
      if(keepMatchSeq) mcols(m)$sequence <- ms
    }else{
      mcols(m)$p3.score <- p3$p3.score
    }
    ms <- subseq(ms, width(r)[[1]]-11L, width(r)[[1]])
    if(keepMatchSeq && !p3.extra) mcols(m)$sequence <- as.factor(ms)
    if(isPureSeed){
      mcols(m)$type <- getMatchTypes(ms, substr(seed,1,7))
      mcols(m)$note <- .TDMD(cbind(type=mcols(m)$type, p3), mirseq=mirseq)
      m <- m[order(seqnames(m), m$type)]
    }else{
      mcols(m) <- cbind(mcols(m), assignKdType(ms, mod))
      mcols(m)$note <- .TDMD(cbind(type=mcols(m)$type, p3), mirseq=mirseq)
      if(maxLogKd[[1]]!=Inf){
        if(all(maxLogKd>=0)) maxLogKd <- -maxLogKd
        if(all(maxLogKd > -10)) maxLogKd <- maxLogKd*1000L
        m <- m[which(m$log_kd <= as.integer(round(maxLogKd[1])))]
      }else{
        m <- m[!is.na(m$log_kd)]
      }
      m <- m[order(seqnames(m), m$log_kd, m$type)]
    }
  }
  rm(ms)
  if(!is.null(mcols(seqs)$ORF.length)){
    mcols(m)$ORF <- start(m) <= mcols(seqs)[as.integer(seqnames(m)),"ORF.length"]
    if(!isPureSeed && maxLogKd[2]!=Inf){
      m <- m[which(!m$ORF | m$log_kd <= as.integer(round(maxLogKd[2])))]
    }
  }
  if(minDist>-Inf){
    if(verbose) message("Removing overlaps...")
    m <- removeOverlappingRanges(m, minDist=minDist, ignore.strand=TRUE)
  }
  names(m) <- NULL
  if(offset!=0) m <- IRanges::shift(m, -offset)
  if(ret=="GRanges") return(m)
  .gr2matchTable(m)
}

.checkSeedsInput <- function(seeds){
  if(is.character(seeds)){
    if(all(nchar(seeds) %in% c(7,8))) return("seed")
    if(all(nchar(seeds) >= 16 & nchar(seeds) <= 26)) return("mirseq")
    stop("If `seeds` is a character vector, it should either be vector of ",
         "seeds (corresponding DNA sequence, 7 or 8 nucleotides), or of mature",
         " miRNA sequences (RNA sequence).\n",
         "The size of the characters in `seeds` is compatible with neither.")
  }
  if(!is(seeds,"KdModel") & !is(seeds,"KdModelList") &
     !(is.list(seeds) && all(sapply(seeds, class)=="KdModel")))
    stop("`seeds` should either be a character vector or an object of class ",
         "`KdModel` or `KdModelList`")
  NULL
}


.gr2matchTable <- function(m, include_name=FALSE, include_ORF=TRUE, p3=TRUE){
  d <- data.frame(transcript=as.factor(seqnames(m)), start=start(m))
  if(include_name) d$miRNA <- as.factor(m$miRNA)
  if(include_ORF && !is.null(m$ORF)) d$ORF <- m$ORF
  d$type <- m$type
  d$log_kd <- m$log_kd
  d$TDMD <- m$TDMD
  for(f in c("p3.mir.bulge","p3.target.bulge","p3.mismatch","p3.matches","p3.score")){
    if(p3 && f %in% colnames(mcols(m)))
      d[[f]] <- mcols(m)[[f]]
  }
  d
}

.check3pParams <- function(p3.params){
  stopifnot(is.list(p3.params))
  def <- list(maxMirLoop=5L, maxTargetLoop=9L, maxLoopDiff=4L, mismatch=TRUE)
  for(f in names(def)){
    if(!(f %in% names(p3.params))) p3.params[[f]] <- def[[f]]
  }
  stopifnot(is.logical(p3.params$mismatch))
  stopifnot(is.numeric(unlist(p3.params[names(def)[1:3]])))
  for(f in names(def)[1:3]) p3.params[[f]] <- as.integer(p3.params[[f]])
  p3.params
}

#' get3pAlignment
#'
#' @param seqs A set of sequences in which to look for matches
#' @param mirseq The sequence of the miRNA
#' @param mir3p.start The position in `mirseq` in which to start looking
#' @param allow.mismatch Logical; whether to allow mismatches
#' @param TGsub Logical; whether to allow T/G substitutions
#' @param maxMirLoop Maximum miRNA loop size
#' @param maxTargetLoop Maximum target loop size
#' @param maxLoopDiff Maximum size difference between miRNA and target loops
#'
#' @return A data.frame with one row for each element of `seqs`.
#' @export
#'
#' @examples
#' get3pAlignment(seqs="NNAGTGTGCCATNN", mirseq="TGGAGTGTGACAATGGTGTTTG")
get3pAlignment <- function(seqs, mirseq, mir3p.start=9L, allow.mismatch=TRUE, 
                           maxMirLoop=5L, maxTargetLoop=9L, maxLoopDiff=4L,
                           TGsub=TRUE){
  mir.3p <- as.character(DNAString(substr(x=mirseq, start=mir3p.start, stop=nchar(mirseq))))
  seqs <- reverseComplement(seqs)
  subm <- .default3pSubMatrix(ifelse(allow.mismatch,-3,-Inf), TG=TGsub)
  al <- pairwiseAlignment(seqs, mir.3p, type="local", substitutionMatrix=subm)
  df <- data.frame( p3.mir.bulge=start(subject(al))-1L,
                    p3.target.bulge=start(pattern(al))-1L )
  df$p3.mismatch <- nchar(mir.3p)-width(pattern(al))-df$p3.mir.bulge
  df$p3.score <- as.integer(score(al))
  diff <- abs(df$p3.mir.bulge-df$p3.target.bulge)
  df$p3.score <- ifelse(diff > 2,df$p3.score - (diff - 2),df$p3.score)
  df[which(df$p3.mir.bulge>maxMirLoop | df$p3.target.bulge>maxTargetLoop | 
             diff>maxLoopDiff), 
     c("p3.mir.bulge","p3.target.bulge","p3.score")] <- 0L
  df
}

.TDMD <- function(m,mirseq){
  is78 <- which(m$type %in% c("8mer","7mer-m8","7mer-a1"))
  m2 <- m[is78,]
  m2$TDMD <- rep(1L,nrow(m2))
  absbulgediff <- abs(m2$p3.mir.bulge-m2$p3.target.bulge)
  w <- which(m2$p3.mismatch==0L & m2$p3.mir.bulge <= 5L & 
             m2$p3.mir.bulge>0L & absbulgediff <= 4L & m2$p3.score >= 6L)
  m2$TDMD[w] <- 2L
  w <- which(m2$p3.mismatch==0L & m2$p3.mir.bulge < 5L & 
             m2$p3.mir.bulge>0L & absbulgediff <= 2L & m2$p3.score >= 6L)
  m2$TDMD[w] <- 3L
  is8 <- m2$type == "8mer"
  m2$not.bound <- nchar(mirseq) - 8L - m2$p3.score
  w <- which(is8 & m2$p3.mismatch<=1L & m2$p3.mir.bulge == 0L & 
               absbulgediff == 0L & m2$not.bound <= 6L)
  m2$TDMD[w] <- 4L
  w <- which(is8 & m2$p3.mismatch == 0L & m2$p3.mir.bulge == 0L & 
               absbulgediff == 0L & m2$not.bound == 0L)
  m2$TDMD[w] <- 5L
  TDMD <- rep(1L,nrow(m))
  TDMD[is78] <- m2$TDMD
  factor(TDMD, levels = 1L:5L, labels = c("-","TDMD?","TDMD","Slicing?","Slicing"))
}

.default3pSubMatrix <- function(mismatch=-3, TG=TRUE){
  subm <- Biostrings::nucleotideSubstitutionMatrix(match=1, mismatch=mismatch)
  l <- c("A","C","G","T","N")
  subm <- subm[l,l]
  if(TG) subm["A", "G"] <- subm["C", "T"] <- 0
  subm
}


#' removeOverlappingRanges
#' 
#' Removes elements from a GRanges that overlap (or are within a given distance of) other 
#' elements higher up in the list (i.e. assumes that the ranges are sorted in order of
#' priority). The function handles overlaps between more than two ranges by successively
#' removing those that overlap higher-priority ones.
#'
#' @param x A GRanges, sorted by (decreasing) importance
#' @param minDist Minimum distance between ranges
#' @param retIndices Logical; whether to return the indices of entries to remove, rather
#' than the filtered GRanges.
#'
#' @return A filtered GRanges, or an integer vector of indices to be removed if 
#' `retIndices==TRUE`.
#' @export
#' @examples
#' gr <- GRanges(seqnames=rep("A",4), IRanges(start=c(10,25,45,35), width=6))
#' removeOverlappingRanges(gr, minDist=7)
removeOverlappingRanges <- function(x, minDist=7L, retIndices=FALSE, ignore.strand=FALSE){
  red <- GenomicRanges::reduce(x, with.revmap=TRUE, min.gapwidth=minDist, ignore.strand=ignore.strand)$revmap
  red <- red[lengths(red)>1]
  if(length(red)==0){
    if(retIndices) return(c())
    return(x)
  }
  i <- seq_along(x)
  toRemove <- c()
  while(length(red)>0){
    ## for each overlap set, we flag the index (relative to i) of the maximum
    ## (i.e. lowest in the list)
    top <- min(red) ## indexes of the top entry per overlap set, relative to i
    ## overlap of non-top entries to the top entries:
    o <- overlapsAny(x[i[-top]],x[i[top]],maxgap=minDist)
    torem <- i[-top][which(o)] ## entries to remove, relative to x
    toRemove <- c(toRemove, torem) ## relative to x
    i <- setdiff(i,torem)
    ## and check again overlaps among this subset (revmap indexes are relative to i)
    red <- GenomicRanges::reduce(x[i], with.revmap=TRUE, min.gapwidth=minDist, ignore.strand=ignore.strand)$revmap
    red <- red[lengths(red)>1]
  }
  if(retIndices) return(toRemove)
  if(length(toRemove)>0) x <- x[-toRemove]
  x
}

# determines target and seed sequence type, converts if necessary, and adds padding/shadow
.prepSeqs <- function(seqs, seeds, shadow=0, pad=c(0,0)){
  if(is.null(names(seqs))) names(seqs) <- paste0("seq",seq_along(seqs))
  seqtype <- .guessSeqType(seqs)
  if(seqtype=="RNA")
    stop("Both the seeds and the target sequences should be in DNA format.")
  ret <- list()
  if( is(seeds, "KdModel") || 
      (is.list(seeds) && all(sapply(seeds, is.list))) ){
    if(is.null(names(seeds)))
      stop("If `seeds` is a list of KdModels, it should be named.")
  }else{
    if(is.null(names(seeds))) n <- names(seeds) <- seeds
    seeds <- gsub("U", "T", seeds)
    names(seeds) <- n
    ret$seeds <- seeds
  }
  seqnms <- names(seqs)
  if(is.character(seqs)) seqs <- DNAStringSet(seqs)
  names(seqs) <- seqnms
  shadow <- max(c(0,shadow-1))
  ret$offset <- max(c(0,pad[1]-max(0,shadow)))
  seqs <- seqs[lengths(seqs)>=(shadow+8)]
  seqs <- subseq(seqs,1+shadow,lengths(seqs))
  seqs <- padAndClip(seqs, views=IRanges( start=1-shadow-ret$offset,
                                          width=lengths(seqs)+shadow+ret$offset+pad[2] ),
                     Lpadding.letter = "N", Rpadding.letter = "N")
  if(!is.null(mcols(seqs)$ORF.length))
    mcols(seqs)$ORF.length <- mcols(seqs)$ORF.length + ret$offset + shadow
  c(ret, list(seqs=seqs))
}

#' getMatchTypes
#' 
#' Given a seed and a set of sequences mathcing it, returns the type of match.
#'
#' @param x A character vector of short sequences.
#' @param seed A 7 or 8 nucleotides string indicating the seed (5' to 3' 
#' sequence of the target RNA). If of length 7, an "A" will be appended.
#'
#' @return A factor of match types.
#' @export
#'
#' @examples
#' x <- c("AACACTCCAG","GACACTCCGC","GTACTCCAT","ACGTACGTAC")
#' getMatchTypes(x, seed="ACACTCCA")
getMatchTypes <- function(x, seed){
  if(is.factor(x)) return(getMatchTypes(levels(x), seed)[as.integer(x)])
  x <- as.character(x)
  y <- rep(1L,length(x))
  seed <- as.character(seed)
  stopifnot(length(seed)==1)
  stopifnot(nchar(seed) %in% 7:8)
  if(nchar(seed)==7) seed <- paste0(seed,"A")
  seed6 <- substr(seed,2,7)
  seedGb6 <- paste0(substr(seed,2,3),"G",substr(seed,4,7))
  seedGb7 <- paste0(substr(seed,1,3),"G",substr(seed,4,7))
  y[grep(paste0("[ACGT][ACGT]",substr(seed,3,8)),x)] <- 2L # 6mer-a1
  y[grep(paste0(substr(seed,1,6),"[ACGT][ACGT]"),x)] <- 3L # 6mer-m8
  if(substr(seedGb6,2,7)!=seed6)
    y[grep(seedGb6,x,fixed=TRUE)] <- 4L # g-bulged 6mer
  y[grep(paste0("[ACGT]",substr(seed,2,7)),x)] <- 5L # 6mer
  if(seedGb6!=substr(seed,1,7)){
    y[grep(paste0(seedGb6,"A|",seedGb7),x)] <- 6L # g-bulged 7mer
    y[grep(paste0(seedGb7,"A"),x,fixed=TRUE)] <- 7L # g-bulged 8mer
  }
  y[grep(paste0("[ACGT]",substr(seed,2,8)),x)] <- 8L # 7mer-a1
  y[grep(substr(seed,1,7),x,fixed=TRUE)] <- 9L # 7mer-m8
  y[grep(seed,x,fixed=TRUE)] <- 10L # 8mer
  factor(y, levels=10L:1L, labels=.matchLevels())
}

.matchLevels <- function(withA=TRUE){
  if(withA) return(c("8mer","7mer-m8","7mer-a1","g-bulged 8mer","g-bulged 7mer",
                     "6mer","g-bulged 6mer","6mer-m8","6mer-a1",
                     "non-canonical"))
  c("7mer","7mer","6mer","g-bulged 7mer","g-bulged 6mer","6mer","g-bulged 6mer",
    "6mer-m8","non-canonical","non-canonical")
}
  

#' runFullScan
#' 
#' @export
runFullScan <- function(species, mods=NULL, UTRonly=TRUE, shadow=15, cores=8, maxLogKd=c(-0.3,-0.3), save.path=FALSE, ...){
  message("Loading annotation")
  suppressPackageStartupMessages({
    library(ensembldb)
    library(AnnotationHub)
    library(BSgenome)
    library(BiocParallel)
  })
  ah <- AnnotationHub()
  species <- match.arg(species, c("mmu","hsa","rno"))
  if(species=="hsa"){
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    if(is.null(mods)) mods <- readRDS(file = "/mnt/schratt/miRNA_KD/Data_Output/mods_hsa_comp.rds")
    ahid <- rev(query(ah, c("EnsDb", "Homo sapiens"))$ah_id)[1]
  }else if(species=="mmu"){
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if(is.null(mods)) mods <- readRDS(file = "/mnt/schratt/miRNA_KD/Data_Output/mods_mmu_comp.rds")
    ahid <- rev(query(ah, c("EnsDb", "Mus musculus"))$ah_id)[2]
  }else if(species=="rno"){
    genome <- BSgenome.Rnorvegicus.UCSC.rn6::BSgenome.Rnorvegicus.UCSC.rn6
    if(is.null(mods)) mods <- readRDS(file = "/mnt/schratt/miRNA_KD/Data_Output/mods_rno_comp.rds")
    ahid <- rev(query(ah, c("EnsDb", "Rattus norvegicus"))$ah_id)[1]
  }
  ensdb <- ah[[ahid]]
  seqlevelsStyle(genome) <- "Ensembl"
  
  # restrict to canonical chromosomes
  canonical_chroms <- seqlevels(genome)[!grepl('_', seqlevels(genome))]
  filt <- SeqNameFilter(canonical_chroms)
  
  message("Extracting transcripts")
  grl_UTR <- suppressWarnings(threeUTRsByTranscript(ensdb, filter=filt))
  seqs <- extractTranscriptSeqs(genome, grl_UTR)
  utr.len <- lengths(seqs)
  if(!UTRonly){
    grl_ORF <- cdsBy(ensdb, by="tx", filter=filt)
    seqs_ORF <- extractTranscriptSeqs(genome, grl_ORF)
    tx_info <- data.frame(strand=unlist(unique(strand(grl_ORF))))
    orf.len <- lengths(seqs_ORF)
    names(orf.len) <- names(grl_ORF)
    tx_info$ORF.length <- orf.len[row.names(tx_info)]
    seqs_ORF[names(seqs)] <- xscat(seqs_ORF[names(seqs)],seqs)
    seqs <- seqs_ORF
    rm(seqs_ORF)
    mcols(seqs)$ORF.length <- orf.len[names(seqs)]
  }else{
    tx_info <- data.frame(strand=unlist(unique(strand(grl_UTR))))
  }
  tx_info$UTR.length <- utr.len[row.names(tx_info)]
  
  message("Scanning with ", cores, " cores")
  if(cores>1){
    BP <- MulticoreParam(cores, progress=TRUE)
  }else{
    BP <- SerialParam()
  }
  m <- findSeedMatches(seqs, mods, shadow=shadow, maxLogKd=maxLogKd, BP=BP, ...)

  metadata(m)$tx_info <- tx_info
  metadata(m)$ah_id <- ahid
  if(!isFALSE(save.path)) save.path <- paste(species, ifelse(UTRonly,"utrs","full"), "matches.rds", sep=".")
  if(isFALSE(save.path)) return(m)
  saveRDS(m, file=save.path)
  rm(m)
  gc()
  message("Saved in: ", save.path)
}


.defaultAggParams <- function(){
  c(ag=-4.863126 , b=0.5735, c=-1.7091, p3=0.04403, 
    coef_utr = -0.28019, coef_orf = -0.08622)
}