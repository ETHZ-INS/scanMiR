#' aggregateMatches
#'
#' @param m A GRanges or data.frame of matches as returned by `findSeedMatches`.
#' @param a The relative concentration of unbound AGO-miRNA complexes.
#' @param b Factor specifying the additional repression by a single bound AGO.
#' @param c Penalty for sites that are found within the ORF region.
#' @param p3 Factor specifying additional repression due to 3p alignment.
#' @param coef_utr Factor specifying additional repression due to UTR length.
#' @param coef_orf Factor specifying additional repression due to ORF length.
#' @param p3.range Range used for 3p alignment.
#' @param keepSiteInfo Logical; whether to return information about site types
#' (default = TRUE). Ignored if `m` does not contain `log_kd` values
#' @param toInt Logical; whether to convert repression scores to integers
#' (default = FALSE).
#' @param BP Pass `BiocParallel::MulticoreParam(ncores, progressbar=TRUE)` to
#' enable multithreading. Note that in addition, `aggregateMatches` uses the
#' \link{data.table} package, which is often set to use multi-threading by
#' default (which would be multiplied by threads determined by `BP`). See
#' \code{\link[data.table]{setDTthreads}} for more information.
#'
#' @return a data.frame containing aggregated repression values and/or
#' information about the numbers and types of matches
#' @export
#' @importFrom stats quantile
#' @importFrom data.table as.data.table .N := dcast rbindlist
#'
#' @examples
#' # we create mock RNA sequences and seeds:
#' seqs <- getRandomSeq(n=10)
#'
#' # load sample KdModel
#' data(SampleKdModel)
#'
#' # find matches
#' matches <- findSeedMatches(seqs, SampleKdModel)
#'
#' # aggregate matches
#' aggregateMatches(matches)
aggregateMatches <- function(m, a=-4.863126 , b=0.5735, c=-1.7091, p3=0.04403,
                           coef_utr = -0.28019, coef_orf = -0.08622,
                           p3.range=c(3L,8L), keepSiteInfo = TRUE, toInt=FALSE,
                           BP=NULL){
  if(is.null(BP)) BP <- BiocParallel::SerialParam()
  if(is(m,"GRanges")){
    m$transcript <- as.factor(seqnames(m))
    m <- mcols(m)
    if(!is.null(m$miRNA)) m$miRNA <- as.factor(m$miRNA)
    m <- as.data.frame(m)
  }
  # if(is.null(m$ORF)) m$ORF <- 0L
  if(!is.null(m$miRNA)){
    # m <- m[,c("miRNA","transcript","ORF","log_kd","p3.score","type")]
    m <- split(m, m$miRNA)
    m <- bplapply(m, BPPARAM=BP, FUN=function(x){
      .aggregate_miRNA(x, a=a, b=b, c=c, p3=p3,coef_utr = coef_utr,
                       coef_orf = coef_orf, keepSiteInfo = keepSiteInfo,
                       toInt=toInt, p3.range=p3.range)
    })
    m <- as.data.frame(data.table::rbindlist(m, use.names=TRUE, fill=TRUE))
    for(f in colnames(m)){
      if(is.numeric(m[[f]])) m[[f]][is.na(m[[f]])] <- 0L
    }
  }else{
    # m <- m[,c("transcript","ORF","log_kd","p3.score","type")]
    m <- .aggregate_miRNA(m, a=a, b=b, c=c, p3=p3,coef_utr = coef_utr,
                          coef_orf = coef_orf, keepSiteInfo = keepSiteInfo,
                          toInt=toInt, p3.range=p3.range)
  }
  return(m)
}


.aggregate_miRNA <- function(m, ll = NULL, a=-4.863126 , b=0.5735, c=-1.7091,
                             p3=0.04403, coef_utr = -0.28019,
                             coef_orf = -0.08622, p3.range=c(3L,8L),
                             keepSiteInfo = FALSE, toInt=FALSE){
  if(is(m,"GRanges")){
    m$transcript <- as.factor(seqnames(m))
    m <- mcols(m)
    if(!is.null(m$miRNA)) m$miRNA <- as.factor(m$miRNA)
    m <- as.data.frame(m)
  }
  if(is.null(m$log_kd)){
    m <- .aggregateSiteInfo(as.data.table(m))
    return(as.data.frame(m, stringsAsFactor=TRUE))
  }
  for(col in c("ORF", "p3.score", "type")) if(is.null(m[[col]])) m[[col]] <- 0L
  if(!is.null(m$miRNA)){
    m <- m[,c("miRNA","transcript","ORF","log_kd","p3.score","type")]
  }else{
    m <- m[,c("transcript","ORF","log_kd","p3.score","type")]
  }
  m <- as.data.table(m)
  m[, ORF:=as.integer(ORF)]
  m[, log_kd:=-log_kd/1000]
  m <- m[log_kd>0]
  if(keepSiteInfo){
    m_type_table <- .aggregateSiteInfo(m)
  }
  m$p3.score <- ifelse(m$type == "non-canonical" , 0L, m$p3.score)
  m$p3.score[m$p3.score>max(p3.range)] <- as.integer(max(p3.range))
  m$p3.score[m$p3.score<min(p3.range)] <- 0L
  m$N <- 1 / (1 + exp(-1 * (a + m$log_kd + c*m$ORF + p3*m$p3.score) ))
  m$log_kd <- NULL
  m$N_bg <- 1 / (1 + exp(-1 * (a  + c*m$ORF) ))
  m <- as.data.frame(rowsum(as.matrix(m[,c("N","N_bg")]), group=m$transcript))
  m <- data.frame( transcript=row.names(m),
                   repression=log2(1+exp(b)*m$N_bg) - log2(1 + exp(b)*m$N))
  if(!is.null(ll) && nrow(m) > 1){
    m <- merge(m,ll,by = "transcript", all.x = TRUE)

    # get the utr score
    m$utr_len <- log10(m$utr_len)
    m$utr_len[is.infinite(m$utr_len) || is.na(m$utr_len)] <- 0
    qu_un <- m[!duplicated(m$transcript),"utr_len"]
    qu <- quantile(qu_un, probs = c(0.05,0.95), na.rm = TRUE)
    m$utr_score <- (m$utr_len - qu[1]) / (qu[2] - qu[1])
    m$utr_score[is.na(m$utr_score)] <- 0

    # get the orf score
    if(sum(m$orf_len, na.rm = TRUE) > 0){
      m$orf_len <- log10(m$orf_len)
      m$orf_len[is.infinite(m$orf_len) || is.na(m$orf_len)] <- 0
      qu_un <- m[!duplicated(m$transcript),"orf_len"]
      qu <- quantile(qu_un, probs = c(0.05,0.95), na.rm = TRUE)
      m$orf_score <- (m$orf_len - qu[1]) / (qu[2] - qu[1])
      m$orf_score[is.na(m$orf_score)] <- 0
    }else{
      m$orf_score <- 0
    }
  m$repression <- m$repression + coef_utr*m$utr_score*m$repression +
    coef_orf*m$orf_score*m$repression
  m <- subset(m, select = - c(orf_len,utr_len,utr_score,orf_score))
  }
  if(toInt) m$repression <- as.integer(round(1000*m$repression))
  m$repression <- ifelse(m$repression >= 0, 0, m$repression)
  if(keepSiteInfo){
    if(!is.null(m$miRNA)){
      m <- merge(m, m_type_table, by=c("transcript","miRNA"), all=TRUE)
    }else{
      m <- merge(m, m_type_table, by=c("transcript"), all=TRUE)
    }
  }
  m
}

.aggregateSiteInfo <- function(m) {
  if(is.null(m$miRNA)){
    sites <- dcast( m[,.(N=.N), by=c("transcript","type")],
                    formula=transcript~type, value.var="N", fill=0L )    
  }else{
    sites <- dcast( m[,.(N=.N), by=c("transcript","miRNA","type")],
                    formula=transcript+miRNA~type, value.var="N",
                    fill=0L )
  }
  ind_bulged <- grep("bulged", names(sites))
  if(length(ind_bulged)>0) {
    ind_nc <- c(ind_bulged, grep("non-canonical", names(sites)))
    sites[,"non-canonical":= rowSums(.SD, na.rm = TRUE), .SDcols = ind_nc]
    sites[, grep("bulged", names(sites)) := NULL]
  }
  ind_6mer <- grep("6mer", names(sites))
  sites[, "6mer" := rowSums(.SD, na.rm = TRUE), .SDcols = ind_6mer]
  ind_7mer <- grep("7mer", names(sites))
  sites[, "7mer" := rowSums(.SD, na.rm = TRUE), .SDcols = ind_7mer]
  for(col in setdiff(c("8mer", "7mer", "6mer", "non-canonical"), names(sites))){
    sites[, (col):=0L]
  }
  if(is.null(m$miRNA)) {
    cols <- c("transcript", "8mer", "7mer", "6mer", "non-canonical")
    return(sites[,..cols])
  } else{
    cols <- c("transcript", "miRNA", "8mer", "7mer", "6mer", "non-canonical")
    return(sites[,..cols])
  }
}


.datatable.aware = TRUE
