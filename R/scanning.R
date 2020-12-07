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
#' @param extra.3p Logical; whether to provide more detail information about the
#' 3' alignment (default FALSE).
#' @param p3.params a named list of parameters for the 3' alignment (see 
#' \code{\link{get3pAlignment}})
#' parameters for the aggregation. Ignored if `ret!="aggregated"`.
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
findSeedMatches <- function( seqs, seeds, seedtype=c("auto", "RNA","DNA"), 
                             shadow=0L, maxLogKd=c(-0.3,-0.3), keepMatchSeq=FALSE, minDist=7L, 
                             onlyCanonical=FALSE, extra.3p=FALSE, maxLoop=15L, mir3p.nts=8L,
                             p3.params=c(maxLoop=15L, mir.nts=8L, minS=1L, maxS=7L, minDist=1L, maxDist=12L),
                             agg.params=c(ag=-5.5, b=0.8656, c=-1.8488, p3=0.2733),
                             ret=c("GRanges","data.frame","aggregated"), 
                             BP=NULL, verbose=NULL, ...){
  ret <- match.arg(ret)
  if(ret=="aggregated"){
    if(!is.list(agg.params)) agg.params <- as.list(agg.params)
    if(!all(c("ag","b","c") %in% names(agg.params)))
      stop("`agg.params` should be a named list with slots `ag`, `b` and `c`.")
  }
  if(is.null(verbose)) verbose <- is(seeds,"KdModel") || length(seeds)==1 || is.null(BP)
  if(verbose) message("Preparing sequences...")
  args <- .prepSeqs(seqs, seeds, seedtype, shadow=shadow, pad=c(maxLoop+mir3p.nts+6L,6L))
  seqs <- args$seqs
  if("seeds" %in% names(args)) seeds <- args$seeds
  offset <- args$offset
  rm(args)

  params <- list(
    shadow=shadow,
    minDist=minDist,
    maxLoop=maxLoop,
    mir3p.nts=mir3p.nts,
    maxLogKd=maxLogKd
  )
  if(ret=="aggregated") params$agg.params <- agg.params
  
  if(is(seeds,"KdModel") || length(seeds)==1){
    if(is.list(seeds[[1]])) seeds <- seeds[[1]]
    params$miRNA <- ifelse(is(seeds, "KdModel"), seeds$name, seeds)
    if(is.null(verbose)) verbose <- TRUE
    m <- .find1SeedMatches(seqs, seeds, keepMatchSeq=keepMatchSeq, minDist=minDist, 
                           maxLogKd=maxLogKd, onlyCanonical=onlyCanonical, 
                           extra.3p=extra.3p, p3.params=p3.params, 
                           offset=offset, verbose=verbose, ret=ret, ...)
    if(length(m)==0) return(m)
    if(ret=="aggregated"){
      if(verbose) message("Aggregating...")
      m <- .aggregate_miRNA(m, ag=agg.params$ag, b=agg.params$b, 
                            c=agg.params$c, toInt=TRUE)
    }
  }else{
    if(is.null(BP)) BP <- SerialParam()
    if(is.null(verbose)) verbose <- !(bpnworkers(BP)>1 | length(seeds)>5)
    m <- bplapply( seeds, BPPARAM=BP, FUN=function(oneseed){
      m <- .find1SeedMatches(seqs=seqs, seed=oneseed, keepMatchSeq=keepMatchSeq,
                   minDist=minDist, maxLogKd=maxLogKd, ret=ret, 
                   onlyCanonical=onlyCanonical, p3.params=p3.params, 
                   offset=offset, extra.3p=extra.3p, verbose=verbose, ...)
      if(length(m)==0) return(m)
      if(ret=="aggregated"){
        if(verbose) message("Aggregating...")
        m <- .aggregate_miRNA(m, ag=agg.params$ag, b=agg.params$b, 
                              c=agg.params$c, toInt=TRUE)
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
                              minDist=1L, onlyCanonical=FALSE, extra.3p=FALSE,
                              p3.params=c(), offset=0L, 
                              ret=c("GRanges","data.frame","aggregated"), 
                              verbose=FALSE){
  ret <- match.arg(ret)
  p3.params <- .check3pParams(p3.params)

  if(is.null(maxLogKd)) maxLogKd <- c(Inf,Inf)
  if(length(maxLogKd)==1) maxLogKd <- rep(maxLogKd,2)
  
  if(verbose) message("Scanning for matches...")
  
  if(isPureSeed <- is.character(seed)){
    pos <- gregexpr(paste0("(?=.",substr(seed,2,7),".)"), seqs, perl=TRUE)
  }else{
    mod <- seed
    seed <- mod$canonical.seed
    if(onlyCanonical){
      patt <- paste0(".",substr(seed,2,7),".")
    }else{
      patt <- .build4mersRegEx(seed)
    }
    pos <- gregexpr(paste0("(?=",patt,")"), seqs, perl=TRUE)
  }
  pos <- lapply(lapply(pos, as.numeric), y=-1, setdiff)
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
  
  if(isPureSeed){
    r <- split(r, seqnames(m))
    names(r) <- NULL
    ms <- as.factor(unlist(extractAt(seqs, r)))
    if(keepMatchSeq) mcols(m)$sequence <- ms
    mcols(m)$type <- getMatchTypes(levels(ms), substr(seed,1,7))[as.integer(ms)]
    m <- m[order(seqnames(m), m$type)]
  }else{
    plen <- p3.params$maxLoop+p3.params$mir.nts
    start(r) <- start(r)-1-plen
    end(r) <- end(r)+2
    r <- split(r, seqnames(m))
    names(r) <- NULL
    ms <- unlist(extractAt(seqs, r))
    names(ms) <- NULL
    mcols(m) <- cbind(mcols(m),
                      get3pAlignment( subseq(ms,1,plen),
                                      mod$mirseq, p3.params=p3.params,
                                      extra.3p=extra.3p ) )
    ms <- subseq(ms, plen, 11+plen)
    if(keepMatchSeq) mcols(m)$sequence <- as.factor(ms)
    mcols(m) <- cbind(mcols(m), assignKdType(ms, mod))
    if(maxLogKd[[1]]!=Inf){
      if(all(maxLogKd>=0)) maxLogKd <- -maxLogKd
      if(all(maxLogKd > -10)) maxLogKd <- maxLogKd*1000
      m <- m[which(m$log_kd <= as.integer(round(maxLogKd[1])))]
    }else{
      m <- m[!is.na(m$log_kd)]
    }
    m <- m[order(seqnames(m), m$log_kd, m$type)]
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

.check3pParams <- function(p3.params){
  p3.pn <- c("maxLoop", "mir.nts", "minS", "maxS", "minDist", "maxDist")
  p3.params <- as.list(p3.params)
  if(any(lengths(p3.params)>1) || !all(p3.pn %in% names(p3.params)))
    stop("`p3.params` should be a list with a single numeric value for each of",
         paste(p3.pn, collapse=", "))
  p3.params
}

.gr2matchTable <- function(m, include_name=FALSE, include_ORF=TRUE){
  d <- data.frame(transcript=as.factor(seqnames(m)), start=start(m))
  if(include_name) d$miRNA <- as.factor(m$miRNA)
  if(include_ORF && !is.null(m$ORF)) d$ORF <- m$ORF
  d$type <- m$type
  d$log_kd <- m$log_kd
  d
}

#' get3pAlignment
#'
#' @param seqs A set of sequences in which to look for matches
#' @param mirseq The sequence of the miRNA
#' @param mir3p.nts The number of miRNA nucelotide in which to look 
#' for matches
#' @param mir3p.start The position in `mirseq` in which to start looking
#' @param extra.3p Logical; whether to provide more detail information about the
#' alignment.
#' @param subm An optional substitution matrix; by default a binary diagonal 
#' matrix with additional 0.65 on G/T is used.
#'
#' @return A data.frame with one row for each element of `seqs`.
#' @export
#'
#' @examples
#' get3pAlignment(target="NNAGTGTGCCATNN", mirseq="TGGAGTGTGACAATGGTGTTTG")
get3pAlignment <- function(seqs, mirseq, mir3p.start=12L, extra.3p=TRUE, 
                           p3.params=c(maxLoop=15L, mir.nts=8L, minS=1L, 
                                       maxS=7L, minDist=1L, maxDist=12L),
                           subm=NULL){
  p3 <- .check3pParams(p3.params)
  mir3p.nts <- as.integer(p3$mir.nts)
  target.len <- width(seqs[1])
  mir.3p <- as.character(reverseComplement(DNAString(
    substr(x=mirseq, start=mir3p.start, 
           stop=min(c(mir3p.start-1+mir3p.nts, nchar(mirseq))))
    )))
  if(is.null(subm)){
    subm <- diag(1,nrow=5,ncol=5)
    colnames(subm) <- row.names(subm) <- c("A","C","G","T","N")
    subm["G", "T"] <- subm["T", "G"] <- 0.65
  }
  al <- pairwiseAlignment(seqs, mir.3p, type="local", substitutionMatrix=subm)
  df <- data.frame( mir.pos.3p=end(subject(al)),
                    target.pos.3p=end(pattern(al)) )
  df$dist.3p <- start(pattern(al))+nchar(mir.3p)-start(subject(al))-target.len
  al <- score(al)
  al[al<p3$minS] <- 0L
  al[al>p3$maxS] <- 0L #or maxS?
  al[-df$dist.3p > p3$maxDist | -df$dist.3p < p3$minDist] <- 0L
  df$align.3p <- al
  if(!extra.3p){
    df$mir.pos.3p <- df$target.pos.3p <- df$dist.3p <- NULL
  }
  df
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
#' @param seed A 7 or 8 nucleotides string indicating the seed (5' to 3' sequence of the
#' target RNA). If of length 7, an "A" will be appended.
#'
#' @return A factor of match types.
#' @export
#'
#' @examples
#' x <- c("AACACTCCAG","GACACTCCGC","GTACTCCAT","ACGTACGTAC")
#' getMatchTypes(x, seed="ACACTCCA")
getMatchTypes <- function(x, seed){
  x <- as.character(x)
  y <- rep(1L,length(x))
  seed <- as.character(seed)
  if(length(seed)!=1 || !(nchar(seed) %in% c(7,8)))
    stop("'seed' should be a string of 7 or 8 characters")
  if(nchar(seed)==7) seed <- paste0(seed,"A")
  seed6 <- substr(seed,2,7)
  y[grep(paste0("[ACGT]","[ACGT]",substr(seed,3,8)),x)] <- 2L # 6mer-a1
  y[grep(paste0(substr(seed,1,6),"[ACGT][ACGT]"),x)] <- 3L # 6mer-m8
  y[grep(paste0("[ACGT]",substr(seed,2,7)),x)] <- 4L # 6mer
  y[grep(paste0("[ACGT]",substr(seed,2,8)),x)] <- 5L # 7mer-a1
  y[grep(substr(seed,1,7),x,fixed=TRUE)] <- 6L # 7mer-m8
  y[grep(seed,x,fixed=TRUE)] <- 7L # 8mer
  factor(y, levels=7:1, labels=c("8mer","7mer-m8","7mer-a1","6mer","6mer-m8",
                                 "6mer-a1","non-canonical"))
}

#' runFullScan
#' 
#' @export
runFullScan <- function(species, mods=NULL, UTRonly=TRUE, shadow=15, cores=8, minLogKd=c(-0.3,-1), save.path=NULL, ...){
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
    ahid <- rev(query(ah, c("EnsDb", "Mus musculus"))$ah_id)[1]
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
  m <- findSeedMatches(seqs, mods, shadow=shadow, minLogKd=minLogKd, BP=BP, ...)

  metadata(m)$tx_info <- tx_info
  metadata(m)$ah_id <- ahid
  if(!is.null(save.path)) save.path <- paste(species, ifelse(UTRonly,"utrs","full"), "matches.rds", sep=".")
  if(isFALSE(save.path)) return(m)
  saveRDS(m, file=save.path)
  rm(m)
  gc()
  message("Saved in: ", save.path)
}
