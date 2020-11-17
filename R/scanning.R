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
#' @param minLogKd Minimum log_kd value to keep (default 0). Set to Inf to disable.
#' @param keepMatchSeq Logical; whether to keep the sequence (including flanking
#' dinucleotides) for each seed match (default FALSE).
#' @param onlyCanonical Logical; whether to restrict the search only to canonical
#' binding sites.
#' @param minDist Integer specifying the minimum distance between matches of the same 
#' miRNA (default 1). Closer matches will be reduced to the highest-affinity. To 
#' disable the removal of overlapping features, use `minDist=-Inf`.
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
findSeedMatches <- function( seqs, seeds, seedtype=c("auto", "RNA","DNA"), shadow=0L, 
                             minLogKd=0, keepMatchSeq=FALSE, maxLoop=10L, mir3p.nts=6L, 
                             minDist=7L, onlyCanonical=FALSE, BP=NULL, verbose=NULL){
  
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
                           onlyCanonical=onlyCanonical, verbose=verbose)
    if(length(m)==0) return(m)
  }else{
    if(is.null(BP)) BP <- SerialParam()
    if(is.null(verbose)) verbose <- !(bpnworkers(BP)>1 | length(seeds)>5)
    m <- bplapply( seeds, seqs=seqs, keepMatchSeq=keepMatchSeq, verbose=verbose, 
                   minDist=minDist, minLogKd=minLogKd, maxLoop=maxLoop, 
                   mir3p.nts=mir3p.nts, onlyCanonical=onlyCanonical, 
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
  m <- IRanges::shift(m, -offset)

  metadata(m)$call.params <- list(
    shadow=shadow,
    minDist=minDist,
    maxLoop=maxLoop,
    mir3p.nts=mir3p.nts
  )
  m
}

# scan for a single seed
.find1SeedMatches <- function(seqs, seed, keepMatchSeq=FALSE, minLogKd=0, maxLoop=10, 
                              mir3p.nts=8, minDist=1, onlyCanonical=FALSE, verbose=FALSE){
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
    mcols(m) <- cbind(mcols(m), assignKdType(ms, mod))
    if(!is.null(minLogKd) && minLogKd!=Inf){
      if(minLogKd>0) minLogKd <- -minLogKd
      if(minLogKd < -10) minLogKd <- minLogKd*1000
      m <- m[which(m$log_kd <= as.integer(round(minLogKd)))]
    }else{
      m <- m[!is.na(m$log_kd)]
    }
    m <- m[order(seqnames(m), m$log_kd, m$type)]
  }
  rm(ms)
  if(minDist>-Inf){
    if(verbose) message("Removing overlaps...")
    m <- removeOverlappingRanges(m, minDist=minDist)
  }
  names(m) <- NULL
  m
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
#'
#' @return A filtered GRanges.
#' @export
#' @examples
#' gr <- GRanges(seqnames=rep("A",4), IRanges(start=c(10,25,45,35), width=6))
#' removeOverlappingRanges(gr, minDist=7)
removeOverlappingRanges <- function(x, minDist=7L){
  red <- GenomicRanges::reduce(x, with.revmap=TRUE, min.gapwidth=minDist)$revmap
  red <- red[lengths(red)>1]
  if(length(red)==0) return(x)
  i <- seq_along(x)
  toRemove <- c()
  while(length(red)>0){
    ## for each overlap set, we flag the index (relative to i) of the maximum
    ## (i.e. lowest in the list)
    top <- min(red) ## indexes of the top entry per overlap set, relative to i
    ## overlap of non-top entries to the top entries:
    o <- GenomicRanges::overlapsAny(x[i[-top]],x[i[top]],maxgap=minDist)
    torem <- i[-top][which(o)] ## entries to remove, relative to x
    toRemove <- c(toRemove, torem) ## relative to x
    i <- setdiff(i,torem)
    ## and check again overlaps among this subset (revmap indexes are relative to i)
    red <- GenomicRanges::reduce(x[i], with.revmap=TRUE, min.gapwidth=minDist)$revmap
    red <- red[lengths(red)>1]
  }
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
runFullScan <- function(species, UTRonly=TRUE, shadow=15, cores=8, minLogKd=0, mods=NULL, ...){
  start <- Sys.time()
  message("Loading annotation")
  suppressPackageStartupMessages({
    library(ensembldb)
    library(AnnotationHub)
    library(BSgenome)
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
  
  message("Preparing sequences")
  grl_UTR <- suppressWarnings(threeUTRsByTranscript(ensdb, filter=filt))
  seqs <- extractTranscriptSeqs(genome, grl_UTR)
  utr.len <- lengths(seqs)
  if(!UTRonly){
    grl_ORF <- cdsBy(ensdb, by="tx", filter=filt)
    seqs_ORF <- extractTranscriptSeqs(genome, grl_ORF)
    tx_info <- data.frame(strand=unlist(unique(strand(grl_ORF))))
    orf.len <- lengths(seqs_ORF)
    tx_info$ORF.length <- orf.len[row.names(tx_info)]
    seqs_ORF[names(seqs)] <- xscat(seqs_ORF[names(seqs)],seqs)
    seqs <- seqs_ORF
    rm(seqs_ORF)
  }else{
    tx_info <- data.frame(strand=unlist(unique(strand(grl_UTR))))
  }
  tx_info$UTR.length <- utr.len[row.names(tx_info)]
  
  message("Scanning")
  m <- findSeedMatches(seqs, mods, shadow=ifelse(UTRonly,shadow,0), minLogKd=minLogKd,
                       BP=MulticoreParam(cores, progress=TRUE), ...)
  if(!UTRonly)
    m$ORF <- start(m) < orf.len[seqlevels(m)][as.integer(seqnames(m))] + shadow

  metadata(m)$tx_info <- tx_info
  metadata(m)$ah_id <- ahid
  saveRDS(m, file=paste(species, ifelse(UTRonly,"utrs","full"), "matches.rds", sep="."))
  end <- Sys.time()
  return(end-start)
}
