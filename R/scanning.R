#' Predicting and characterizing miRNA binding sites
#'
#' `findSeedMatches` takes a set of sequences and a set of miRNAs (given either
#' as target seeds, mature miRNA sequences, or a \code{\link{KdModelList}}).
#'
#' @param seqs A character vector or `DNAStringSet` of DNA sequences in which to
#' look.
#' @param seeds A character vector of 7-nt seeds to look for. If RNA, will be
#' reversed and complemented before matching. If DNA, they are assumed to be
#' the target sequence to look for. Alternatively, a list of objects of class
#' `KdModel` or an object of class `KdModelList` can be given.
#' @param shadow Integer giving the shadow, i.e. the number of nucleotides
#'  hidden at the beginning of the sequence (default 0).
#' @param onlyCanonical Logical; whether to restrict the search only to
#' canonical binding sites.
#' @param maxLogKd Maximum log_kd value to keep. This has a major impact on the
#' number of sites returned, and hence on the memory requirements. Set to Inf
#' to disable (_not_ recommended when running large scans!).
#' @param keepMatchSeq Logical; whether to keep the sequence (including flanking
#' dinucleotides) for each seed match (default FALSE).
#' @param minDist Integer specifying the minimum distance between matches of the
#' same miRNA (default 7). Closer matches will be reduced to the
#' highest-affinity. To disable the removal of overlapping features, use
#' `minDist=-Inf`.
#' @param p3.extra Logical; whether to keep extra information about 3'
#' alignment. Disable (default) this when running large scans, otherwise you
#' might hit your system's memory limits.
#' @param p3.params Named list of parameters for 3' alignment with slots
#' `maxMirLoop` (integer, default = 5), `maxTargetLoop` (integer, default = 9),
#' `maxLoopDiff` (integer, default = 4), and `mismatch`
#' (logical, default = TRUE).
#' @param agg.params A named list with slots `a`, `b`, `c`, `p3`, `coef_utr`,
#' `coef_orf` and `keepSiteInfo` indicating the parameters for the aggregation.
#' Ignored if `ret!="aggregated"`. For further details see documentation of
#' `aggregateMatches`.
#' @param ret The type of data to return, either "GRanges" (default),
#' "data.frame", or "aggregated" (aggregates affinities/sites for each
#' seed-transcript pair).
#' @param BP Pass `BiocParallel::MulticoreParam(ncores, progressbar=TRUE)` to
#' enable multithreading.
#' @param verbose Logical; whether to print additional progress messages
#' (default on if not multithreading)
#' @param n_seeds Integer; the number of seeds that are processed in parallel to
#' avoid memory issues.
#' @param useTmpFiles Logical; whether to write results for single miRNAs in
#' temporary files (ignored when scanning for a single seed). Alternatively,
#' `useTmpFiles` can be a character vector of length 1 indicating the path to
#' the directory in which to write temporary files.
#' @param keepTmpFiles Logical; whether to keep the temporary files at the end
#' of the process; ignored if `useTmpFiles=FALSE`. Temporary files are removed
#' only upon successful completion of the function, meaning that they will not
#' be deleted in case of errors.
#'
#' @return A GRanges of all matches. If `seeds` is a `KdModel` or `KdModelList`,
#' the `log_kd` column will report the ln(Kd) multiplied by 1000, rounded and
#' saved as an integer. If `ret!="GRanges`, returns a data.frame.
#'
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @importFrom GenomeInfoDb seqlevels
#' @import Biostrings GenomicRanges
#' @importFrom S4Vectors mcols mcols<- metadata metadata<- Rle DataFrame
#' @importFrom IRanges IRanges
#' @export
#'
#' @examples
#' # we create mock RNA sequences and seeds:
#' seqs <- getRandomSeq(n=10)
#' seeds <- c("AAACCAC", "AAACCUU")
#' findSeedMatches(seqs, seeds)
findSeedMatches <- function( seqs, seeds, shadow=0L, onlyCanonical=FALSE,
                             maxLogKd=c(-1,-1.5), keepMatchSeq=FALSE,
                             minDist=7L, p3.extra=FALSE,
                             p3.params=list(maxMirLoop=5L, maxTargetLoop=9L,
                                            maxLoopDiff=4L, mismatch=TRUE),
                             agg.params=.defaultAggParams(),
                             ret=c("GRanges","data.frame","aggregated"),
                             BP=NULL, verbose=NULL, n_seeds=NULL,
                             useTmpFiles=FALSE, keepTmpFiles=FALSE){
  p3.params <- .check3pParams(p3.params)
  if(length(maxLogKd)==1) maxLogKd <- rep(maxLogKd,2)
  seqtype <- .guessSeqType(seqs)
  if(seqtype=="RNA")
    stop("Target sequences should be in DNA format.")
  if(is.character(seqs)) seqs <- DNAStringSet(seqs)
  mcols(seqs)$length <- length.seqs <- width(seqs)
  tx_info <- mcols(seqs)
  if(is.null(tx_info$ORF.length)){
    hasORF <- FALSE
    utr.length <- length.seqs
    orf.length <- 0L
    mcols(seqs)$ORF.length <- orf.length
    mcols(seqs)$C.length <- orf.length
  }else{
    hasORF <- TRUE
    orf.length <- tx_info$ORF.length
    utr.length <- ifelse(length.seqs > orf.length, length.seqs - orf.length, 0L)
    mcols(seqs)$ORF.length <- orf.length
    mcols(seqs)$C.length <- orf.length + shadow
    shadow <- 0L
  }
  tx_info$UTR.length <- utr.length
  length.info <- cbind(orf.length, utr.length)
  if(!is.null(names(seqs))) row.names(length.info) <- names(seqs)

  ret <- match.arg(ret)
  if(ret=="aggregated"){
    if(!is.list(agg.params)) agg.params <- as.list(agg.params)
    if(!all(names(agg.params) %in% names(.defaultAggParams())))
      stop("`agg.params` should be a named list with slots among ",
           "`a`, `b`, `c`, `p3`, `coef_utr`, `coef_orf` and `keepSiteInfo.")
    for(name in names(.defaultAggParams())){
      if(!name %in% names(agg.params)){
        agg.params[[name]] <- .defaultAggParams()[[name]]
      }
    }
  }
  seedInputType <- .checkSeedsInput(seeds)
  if(is.null(verbose)) {
    verbose <- is(seeds,"KdModel") || length(seeds)==1 || is.null(BP)
  }
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
    useTmpFiles <- FALSE
    if(is.list(seeds[[1]])) seeds <- seeds[[1]]
    params$miRNA <- ifelse(is(seeds, "KdModel"), seeds$name, seeds)
    if(is.null(verbose)) verbose <- TRUE
    m <- .find1SeedMatches(seqs, seeds, keepMatchSeq=keepMatchSeq,
                           minDist=minDist, maxLogKd=maxLogKd,
                           onlyCanonical=onlyCanonical, p3.extra=p3.extra,
                           p3.params=p3.params, offset=offset,
                           verbose=verbose, ret=ret)
    if(length(m)==0) return(m)
    if(ret=="aggregated"){
      if(verbose) message("Aggregating...")
      ll <- as.data.frame(length.info)
      ll$transcript <- row.names(ll)
      m <- .aggregate_miRNA(m, ll, a=agg.params$a, b=agg.params$b,
                            c=agg.params$c, p3 = agg.params$p3,
                            coef_utr = agg.params$coef_utr,
                            coef_orf = agg.params$coef_orf,
                            keepSiteInfo = agg.params$keepSiteInfo, toInt=TRUE)
    }
  }else{
    if(is.null(BP)) BP <- SerialParam()
    if(is.null(verbose)) verbose <- !(bpnworkers(BP)>1 | length(seeds)>5)
    if(is.null(n_seeds)) n_seeds <- length(seeds)
    if(isTRUE(useTmpFiles) ||
       (is.character(useTmpFiles) && length(useTmpFiles)==1)){
      if(is.logical(useTmpFiles)){
        tmpDir <- tempdir()
      }else{
        tmpDir <- useTmpFiles
        useTmpFiles <- TRUE
      }
      stopifnot(!file.access(tmpDir, mode=2))
    }
    if(verbose) {
      if(is.numeric(BP$workers)) {
        message("Scanning with ", n_seeds, "seeds at a time on ", BP$workers,
                " cores...")
        if(useTmpFiles) message("Temporary writing results in:\n", tmpDir)
      } else {
        message("Scanning with ", n_seeds, " seeds at a time...")
      }
    }
    split_seeds <- split(seeds, ceiling(seq_along(seeds)/n_seeds))
    m <- lapply(split_seeds, function(seeds) {
      m <- bplapply(seeds, BPPARAM=BP, FUN=function(oneseed){
        if(ret=="aggregated" && is(oneseed, "KdModel") && agg.params$p3 == 0) {
          oneseed$mirseq <- NULL
        }
        m <- .find1SeedMatches(seqs=seqs, seed=oneseed,
                               keepMatchSeq=keepMatchSeq, minDist=minDist,
                               maxLogKd=maxLogKd, p3.extra=p3.extra,
                               onlyCanonical=onlyCanonical, p3.params=p3.params,
                               ret=ret, offset=offset, verbose=verbose)
        if(ret=="aggregated"){
          if(verbose) message("Aggregating...")
          if(length(m)==0) return(data.frame())
          ll <- as.data.frame(length.info)
          ll$transcript <- row.names(ll)
          m <- .aggregate_miRNA(m, ll, a=agg.params$a, b=agg.params$b,
                                c=agg.params$c, p3 = agg.params$p3,
                                coef_utr = agg.params$coef_utr,
                                coef_orf = agg.params$coef_orf,
                                keepSiteInfo = agg.params$keepSiteInfo,
                                toInt=TRUE)
        }
        if(!useTmpFiles) return(m)
        f <- paste0(ifelse(is.character(oneseed), oneseed, oneseed$name), "_")
        f <- tempfile(f, tmpdir=tmpDir)
        saveRDS(m, f, compress=FALSE)
        rm(m)
        gc(verbose = FALSE)
        f
      })
    })
    m <- do.call(c, m)
    names(m) <- NULL

    if(!is.character(seeds)) {
      seeds <- vapply(seeds, FUN=function(x){
        if(is.null(x$name)) return(x$canonical.seed)
        x$name
      }, character(1))
    }
    names(m) <- seeds

    if(useTmpFiles){
      if(verbose) message("Scan completed, reading and aggregating results ",
                          "from temporary files...")
      ff <- unlist(m,use.names=FALSE)
      m <- lapply(m, FUN=function(x) readRDS(x))
    }

    if(ret=="GRanges"){
      m <- .unlistGRL(m, .id="miRNA")
      names(m) <- row.names(m) <- NULL
      mcols(m)$miRNA <- Rle(as.factor(mcols(m)$miRNA))
    }else{
      m <- lapply(m, as.data.frame)
      m <- data.table::rbindlist(m, fill=TRUE, use.names=TRUE,
                                               idcol="miRNA")
      for(f in colnames(m)){
        if(is.numeric(m[[f]])) m[[f]][is.na(m[[f]])] <- 0L
      }
      m$miRNA <- as.factor(m$miRNA)
      row.names(m) <- NULL
      m$transcript <- as.factor(m$transcript)
    }
  }
  if(ret=="GRanges"){
    metadata(m)$call.params <- params
    metadata(m)$tx_info <- as.data.frame(tx_info)
  }else{
    m <- as.data.frame(m)
    attr(m, "call.params") <- params
    attr(m, "tx_info") <- as.data.frame(tx_info)
  }
  if(ret=="aggregated" && onlyCanonical) m[["non-canonical"]] <- NULL
  if(useTmpFiles && !keepTmpFiles) unlink(ff)
  gc(verbose = FALSE, full = TRUE)
  return(m)
}

# scan for a single seed
.find1SeedMatches <- function(seqs, seed, keepMatchSeq=FALSE, maxLogKd=-1,
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
      if(verbose && p3.extra)
        warning("`p3.extra` ignored when input is only a seed")
      p3.extra <- FALSE
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
  m <- GRanges(rep(names(seqs), lengths(pos)),
               IRanges(start=unlist(pos), width=8))
  rm(pos)
  m <- m[order(seqnames(m))]

  if(verbose) message("Extracting sequences and characterizing matches...")
  r <- ranges(m)
  if(isPureSeed){
    r <- split(r, seqnames(m))
    names(r) <- NULL
    ms <- as.factor(unlist(extractAt(seqs[seqlevels(m)], r))) # 8mers
    mcols(m)$type <- getMatchTypes(as.factor(ms), substr(seed,1,7))
    if(keepMatchSeq && !p3.extra) mcols(m)$sequence <- ms
    m <- m[order(seqnames(m), m$type)]
  }else{
    start(r) <- start(r)-2L
    end(r) <- end(r)+2L
    r <- split(r, seqnames(m))
    names(r) <- NULL
    ms <- unlist(extractAt(seqs[seqlevels(m)], r))  # 12mers
    mcols(m) <- cbind(mcols(m), assignKdType(ms, mod))
    if(keepMatchSeq) mcols(m)$sequence <- ms
    if(maxLogKd[[1]]!=Inf){
      if(all(maxLogKd>=0)) maxLogKd <- -maxLogKd
      if(all(maxLogKd > -10)) maxLogKd <- maxLogKd*1000L
      m <- m[which(m$log_kd <= as.integer(round(maxLogKd[1])))]
    }else{
      m <- m[!is.na(m$log_kd)]
    }
    m <- m[order(seqnames(m), m$log_kd, m$type)]
  }
  rm(ms)
  if(minDist>-Inf){
    if(verbose) message("Removing overlaps...")
    m <- removeOverlappingRanges(m, minDist=minDist, ignore.strand=TRUE)
  }
  if(!is.null(mirseq)){
    if(verbose) message("Performing 3' alignment...")
    maxLoop <- max(unlist(p3.params[c("maxMirLoop","maxTargetLoop")]))
    plen <- maxLoop+nchar(mirseq)-8L
    r <- ranges(m)
    start(r) <- start(r)-plen-1L
    end(r) <- start(r)+plen
    if(isPureSeed && keepMatchSeq && p3.extra) end(r) <- end(r) + 10L
    r <- split(r, seqnames(m))
    names(r) <- NULL
    ms <- unlist(extractAt(seqs[seqlevels(m)], r))
    rm(r)
    names(ms) <- NULL
    p3 <- get3pAlignment( subseq(ms,1,plen+1L), mirseq, siteType=mcols(m)$type,
                          allow.mismatch=p3.params$mismatch,
                          maxMirLoop=p3.params$maxMirLoop,
                          maxLoopDiff=p3.params$maxLoopDiff,
                          maxTargetLoop=p3.params$maxTargetLoop )
    if(p3.extra){
      mcols(m) <- cbind(mcols(m), p3)
      if(keepMatchSeq){
        if(isPureSeed){
          mcols(m)$sequence <- ms
        }else{
          mcols(m)$sequence <- xscat(ms, subseq(mcols(m)$sequence,3L))
        }
      }
    }else{
      mcols(m)$p3.score <- p3$p3.score
      mcols(m)$note <- p3$note
    }
    rm(ms)
    mcols(m)$note <- Rle(mcols(m)$note)
  }
  if(!is.null(mcols(seqs)$C.length) && !all(mcols(seqs)$ORF.length == 0)) {
    mcols(m)$ORF <-
      start(m) <= mcols(seqs[seqlevels(m)])[as.integer(seqnames(m)),"C.length"]
    if(!isPureSeed && maxLogKd[2]!=Inf && maxLogKd[2]!=maxLogKd[1]){
      m <- m[which(!m$ORF | m$log_kd <= as.integer(round(maxLogKd[2])))]
    }
    mcols(m)$ORF <- Rle(mcols(m)$ORF)
  }
  names(m) <- NULL
  if(offset!=0) m <- IRanges::shift(m, -offset)
  if(ret=="data.frame") m <- .gr2matchTable(m)
  return(m)
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
  if(!is(seeds,"KdModel") && !is(seeds,"KdModelList") &&
     !(is.list(seeds) && all(vapply(seeds, class, character(1))=="KdModel")))
    stop("`seeds` should either be a character vector or an object of class ",
         "`KdModel` or `KdModelList`")
  NULL
}


.gr2matchTable <- function(m, include_name=FALSE, include_ORF=TRUE, p3=TRUE){
  d <- DataFrame(transcript=as.factor(seqnames(m)), start=start(m))
  if(include_name) d$miRNA <- as.factor(m$miRNA)
  if(include_ORF && !is.null(m$ORF)) d$ORF <- m$ORF
  d$type <- m$type
  d$log_kd <- m$log_kd
  d$TDMD <- m$TDMD
  d$note <- m$note
  for(f in c("p3.mir.bulge",
             "p3.target.bulge",
             "p3.mismatch",
             "p3.matches",
             "p3.score")){
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
  stopifnot(is.numeric(unlist(p3.params[names(def)[seq_len(3)]])))
  for(f in names(def)[seq_len(3)]) p3.params[[f]] <- as.integer(p3.params[[f]])
  p3.params
}

#' Finds 3' complementary binding of a miRNA
#'
#' Performs a local alignment of the miRNA 3' sequence (determined by
#' `mir3p.start`) on given the given sequences.
#'
#' @param seqs A set of sequences in which to look for 3' matches (i.e. upstream
#' of the seed match)
#' @param mirseq The sequence of the mature miRNA
#' @param siteType The optional type of seed-complementarity, as returned by 
#' \code{\link{getMatchTypes}}. This is needed to identify slicing/TDMD sites.
#' If given, should be a vector of the same length as `seqs`.
#' @param mir3p.start The position in `mirseq` in which to start looking
#' @param allow.mismatch Logical; whether to allow mismatches
#' @param TGsub Logical; whether to allow T/G substitutions.
#' @param maxMirLoop Maximum miRNA loop size
#' @param maxTargetLoop Maximum target loop size
#' @param maxLoopDiff Maximum size difference between miRNA and target loops
#'
#' @return A data.frame with one row for each element of `seqs`, indicating the
#' size of the miRNA bulge, the size of the target mRNA bulge, the number of
#' mismatches at the 3' end, and the partial 3' alignment score (i.e. roughly
#' the number of consecutive matching nucleotides)
#'
#' @export
#'
#' @examples
#' get3pAlignment(seqs="NNAGTGTGCCATNN", mirseq="TGGAGTGTGACAATGGTGTTTG")
get3pAlignment <- function(seqs, mirseq, mir3p.start=9L, allow.mismatch=TRUE,
                           maxMirLoop=5L, maxTargetLoop=9L, maxLoopDiff=4L,
                           TGsub=TRUE, siteType=NULL){
  if(!is.null(siteType)) stopifnot(length(seqs)==length(siteType))
  mir.3p <- as.character(DNAString(substr(x=mirseq,
                                          start=mir3p.start,
                                          stop=nchar(mirseq))))
  
  if(!is(seqs, "XStringSet") && !is(seqs, "XString")) seqs <- DNAString(seqs)
  seqs <- reverseComplement(seqs)
  subm <- .default3pSubMatrix(ifelse(allow.mismatch,-3L,-Inf), TG=TGsub)
  al <- pairwiseAlignment(seqs, mir.3p, type="local", substitutionMatrix=subm)
  df <- data.frame( p3.mir.bulge=start(subject(al))-1L,
                    p3.target.bulge=start(pattern(al))-1L )
  df$p3.mismatch <- nchar(mir.3p)-width(pattern(al))-df$p3.mir.bulge
  df$p3.score <- as.integer(score(al))
  diff <- abs(df$p3.mir.bulge-df$p3.target.bulge)
  df$p3.score <- ifelse(diff > 2L, df$p3.score - (diff - 2L), df$p3.score)
  df[which(df$p3.mir.bulge>maxMirLoop | df$p3.target.bulge>maxTargetLoop |
             diff>maxLoopDiff),
     c("p3.mir.bulge","p3.target.bulge","p3.score")] <- 0L
  if(!is.integer(df$p3.score)) df$p3.score <- as.integer(round(df$p3.score))
  if(!is.null(siteType)){
    df$type <- siteType
    df$note <- .TDMD(df, mirseq, acceptWobble=TGsub)
    n_9_11 <- substr(as.character(mirseq),nchar(mirseq)-10,nchar(mirseq)-8)
    if(length(w <- which(df$note %in% c("Slicing","Slicing?")))>0 &&
       sum(w2 <- !grepl(paste0(n_9_11,"$"), as.character(seqs)[w]))>0){
      # for slicing sites, ensure that positions 9-11 are complementary
      df$note[w[which(!w2)]] <- "-"
    }
    df$type <- NULL
  }
  df
}

.TDMD <- function(m, mirseq, acceptWobble=TRUE){
  tt <- c("8mer","7mer-m8","7mer-a1","g-bulged 8mer")
  if(acceptWobble) tt <- c(tt, "wobbled 8mer","wobbled 7mer")
  is78 <- which(m$type %in% tt)
  m2 <- m[is78,]
  isWobble <- m2$type %in% c("wobbled 8mer","wobbled 7mer")
  isBulged <- m2$type=="g-bulged 8mer"
  m2$TDMD <- rep(1L,nrow(m2))
  absbulgediff <- abs(m2$p3.mir.bulge-m2$p3.target.bulge)
  w <- which(m2$p3.mismatch<=1L & m2$p3.mir.bulge<=5L & !isBulged & !isWobble &
             m2$p3.mir.bulge > 0L & absbulgediff <= 4L & m2$p3.score >= 6L)
  m2$TDMD[w] <- 2L
  w <- which(m2$TDMD == 2L & m2$p3.mir.bulge < 5L & absbulgediff <= 2L)
  m2$TDMD[w] <- 3L
  w <- which( m2$type != "7mer-a1"  &
              m2$p3.score >= 7L & m2$p3.mir.bulge == 0L & absbulgediff == 0L )
  m2$TDMD[w] <- 4L
  m2$not.bound <- nchar(mirseq) - 8L - m2$p3.score
  w <- which(m2$TDMD == 4L & m2$p3.mismatch <= 1L & m2$not.bound <= 1L &
               !isWobble & !isBulged)
  m2$TDMD[w] <- 5L
  TDMD <- rep(1L,nrow(m))
  TDMD[is78] <- m2$TDMD
  factor(TDMD, levels = 1L:5L, labels = c("-",
                                          "TDMD?",
                                          "TDMD",
                                          "Slicing?",
                                          "Slicing"))
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
#' Removes elements from a GRanges that overlap (or are within a given distance
#' of) other elements higher up in the list (i.e. assumes that the ranges are
#' sorted in order of priority). The function handles overlaps between more than
#' two ranges by successively removing those that overlap higher-priority ones.
#'
#' @param x A GRanges, sorted by (decreasing) importance.
#' @param minDist Minimum distance between ranges.
#' @param retIndices Logical; whether to return the indices of entries to
#' remove, rather than the filtered GRanges.
#' @param ignore.strand Logical. Whether the strand of the input ranges should
#' be ignored or not.
#'
#' @return A filtered GRanges, or an integer vector of indices to be removed if
#' `retIndices==TRUE`.
#' @export
#' @import GenomicRanges
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames=rep("A",4), IRanges(start=c(10,25,45,35), width=6))
#' removeOverlappingRanges(gr, minDist=7)
removeOverlappingRanges <- function(x, minDist=7L, retIndices=FALSE,
                                    ignore.strand=FALSE){
  red <- GenomicRanges::reduce(x, with.revmap=TRUE, min.gapwidth=minDist,
                               ignore.strand=ignore.strand)$revmap
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
    o <- IRanges::overlapsAny(x[i[-top]],x[i[top]],maxgap=minDist)
    torem <- i[-top][which(o)] ## entries to remove, relative to x
    toRemove <- c(toRemove, torem) ## relative to x
    i <- setdiff(i,torem)
    ## and check again overlaps among this subset (revmap ind are relative to i)
    red <- GenomicRanges::reduce(x[i], with.revmap=TRUE, min.gapwidth=minDist,
                                 ignore.strand=ignore.strand)$revmap
    red <- red[lengths(red)>1]
  }
  if(retIndices) return(toRemove)
  if(length(toRemove)>0) x <- x[-toRemove]
  x
}

## determines target and seed sequence type, converts if necessary,
## and adds padding/shadow
.prepSeqs <- function(seqs, seeds, shadow=0, pad=c(0,0)){
  if(is.null(names(seqs))) names(seqs) <- paste0("seq",seq_along(seqs))
  seqtype <- .guessSeqType(seqs)
  if(seqtype=="RNA")
    stop("Both the seeds and the target sequences should be in DNA format.")
  ret <- list()
  if( is(seeds, "KdModel") ||
      (is.list(seeds) && all(vapply(seeds, is.list, logical(1)))) ){
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
  seqs <- padAndClip(seqs, views=IRanges(start=1-shadow-ret$offset,
                                         width=lengths(seqs)+
                                           shadow+ret$offset+pad[2]),
                     Lpadding.letter = "N", Rpadding.letter = "N")
  if(!is.null(mcols(seqs)$C.length))
    mcols(seqs)$C.length <- mcols(seqs)$C.length + ret$offset + shadow
  c(ret, list(seqs=seqs))
}

#' getMatchTypes
#'
#' Given a seed and a set of sequences matching it, returns the type of match.
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
getMatchTypes <- function(x, seed, checkWobble=TRUE){
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
  y[grep(paste0("[ACGTN][ACGTN]",substr(seed,3,8)),x)] <- 2L # 6mer-a1
  y[grep(paste0(substr(seed,1,6),"[ACGTN][ACGTN]"),x)] <- 3L # 6mer-m8
  if(substr(seedGb6,2,7)!=seed6)
    y[grep(seedGb6,x,fixed=TRUE)] <- 4L # g-bulged 6mer
  y[grep(paste0("[ACGTN]",seed6),x)] <- 5L # 6mer
  if(seedGb6!=substr(seed,1,7)){
    y[grep(paste0(seedGb6,"A|",seedGb7),x)] <- 6L # g-bulged 7mer
    y[grep(paste0(seedGb7,"A"),x,fixed=TRUE)] <- 7L # g-bulged 8mer
  }
  if(checkWobble){
    y[.isWobble(x,seed)] <- 8L # wobbled 7-mer
    y[.isWobble(x,seed,FALSE)] <- 9L # wobbled 8-mer
  }
  y[grep(paste0("[ACGTN]",substr(seed,2,8)),x)] <- 10L # 7mer-a1
  y[grep(substr(seed,1,7),x,fixed=TRUE)] <- 11L # 7mer-m8
  y[grep(seed,x,fixed=TRUE)] <- 12L # 8mer
  factor(y, levels=12L:1L, labels=.matchLevels())
}

.isWobble <- function(x, seed, allow7mer=TRUE, positions=2:5){
  seed <- strsplit(seed,"")[[1]]
  if(allow7mer) seed <- head(seed,7)
  wo <- c(A="G", C="T")
  seeds <- vapply(intersect(which(seed %in% names(wo)), positions), 
                  FUN.VALUE=character(1), FUN=function(i){
    seed[i] <- wo[seed[i]]
    paste(seed, collapse="")
  })
  if(length(seeds)==0) return(rep(FALSE, length(x)))
  grepl(paste(seeds, collapse="|"), x)
}

.matchLevels <- function(withA=TRUE){
  if(withA) return(c("8mer","7mer-m8","7mer-a1","wobbled 8mer","wobbled 7mer",
                     "g-bulged 8mer","g-bulged 7mer","6mer","g-bulged 6mer",
                     "6mer-m8","6mer-a1","non-canonical"))
  c("7mer","7mer","6mer","wobbled 8mer","wobbled 7mer","g-bulged 7mer",
    "g-bulged 6mer","6mer","g-bulged 6mer","6mer-m8","non-canonical",
    "non-canonical")
}

.defaultAggParams <- function(){
  list(a=0.007726,
    b=0.5735,
    c=0.1810,
    p3=0.04403,
    coef_utr = -0.28019,
    coef_orf = -0.08622,
    keepSiteInfo = TRUE
    )
}

#' @importFrom IRanges FactorList
.unlistGRL <- function(m, .id=NULL, tryStandard=TRUE, verbose=FALSE){
  # to avoid c-stack errors on some systems
  if(tryStandard){
    gr <- try(unlist(GRangesList(m)), silent=TRUE)
    if(!is(gr,"try-error")){
      if(!is.null(.id))
        mcols(gr)[[.id]] <- rep(as.factor(names(m)), lengths(m))
      return(gr)
    }
  }
  if(verbose) message("Ranges...")
  gr <- GRanges(unlist(FactorList(lapply(m,seqnames)), use.names=FALSE),
                IRanges(unlist(lapply(m, start)), unlist(lapply(m,end))))
  for(f in colnames(mcols(m[[1]]))){
    if(verbose) message(f)
    x <- lapply(m, FUN=function(x) mcols(x)[[f]])
    if(is.factor(x[[1]])){
      x <- unlist(FactorList(x), use.names=FALSE)
    }else{
      x <- unlist(x)
    }
    mcols(gr)[[f]] <- x
  }
  if(verbose) message("Names...")
  if(!is.null(.id))
    mcols(gr)[[.id]] <- rep(as.factor(names(m)), lengths(m))
  gr
}
