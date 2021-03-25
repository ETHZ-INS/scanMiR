#' aggregateSites
#'
#' @param m A GRanges or data.frame of matches.
#' @param ag 
#' @param b 
#' @param c 
#' @param p3 
#' @param toInt Logical; whether to convert repression scores to integers
#' @param BP BPPARAM argument for multithreading
#'
#' @return a data.frame
#' @export
aggregateSites <- function(m, ag=-4.863126 , b=0.5735, c=-1.7091, p3=0.04403, 
                           coef_utr = -0.28019, coef_orf = -0.08622, p3.range=c(3L,8L), 
                           keepSiteInfo = FALSE, toInt=FALSE, BP=NULL){
  if(is.null(BP)) BP <- BiocParallel::SerialParam()
  if(is(m,"GRanges")){
    m$transcript <- as.factor(seqnames(m))
    m <- mcols(m)
    if(!is.null(m$miRNA)) m$miRNA <- as.factor(m$miRNA)
    m <- as.data.frame(m)
  }
  if(is.null(m$ORF)) m$ORF <- 0L
  if(!is.null(m$miRNA)){
    m <- m[,c("miRNA","transcript","ORF","log_kd","p3.score","type")]
    m <- split(m, m$miRNA)
    m <- bplapply(m, BPPARAM=BP, FUN=function(x){
      .aggregate_miRNA(x, ag=ag, b=b, c=c, p3=p3,coef_utr = coef_utr, coef_orf = coef_orf, 
                       keepSiteInfo = keepSiteInfo, toInt=toInt, p3.range=p3.range)
    })
    dplyr::bind_rows(m, .id="miRNA")
    m <- dplyr::mutate_if(m, is.numeric, tidyr::replace_na, 0L)
  }else{
    m <- m[,c("transcript","ORF","log_kd","p3.score","type")]
    m <- .aggregate_miRNA(m, ag=ag, b=b, c=c, p3=p3,coef_utr = coef_utr, coef_orf = coef_orf, 
                          keepSiteInfo = keepSiteInfo, toInt=toInt, p3.range=p3.range)
    m
  }
}


.aggregate_miRNA <- function(m, ll = NULL, ag=-4.863126 , b=0.5735, c=-1.7091, p3=0.04403, 
                             coef_utr = -0.28019,coef_orf = -0.08622, p3.range=c(3L,8L), 
                             keepSiteInfo = FALSE, toInt=FALSE){
  if(is(m,"GRanges")){
    m$transcript <- as.factor(seqnames(m))
    m <- mcols(m)
    if(!is.null(m$miRNA)) m$miRNA <- as.factor(m$miRNA)
    m <- as.data.frame(m)
  }
  if(is.null(m$ORF)) m$ORF <- 0L
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
    if(!is.null(m$miRNA)){
      m_type_table <- dcast( m[,.(N=.N), by=c("transcript","miRNA","type")],
                           formula=transcript+miRNA~type, value.var="N", fill=0L)
    }else{
      m_type_table <- dcast( m[,.(N=.N), by=c("transcript","type")],
                             formula=transcript~type, value.var="N", fill=0L)
    }
  }
  if(is.null(m$p3.score)) m$p3.score <- 0L
  m$p3.score <- ifelse(m$type == "non-canonical" , 0L, m$p3.score)
  m$p3.score[m$p3.score>max(p3.range)] <- as.integer(max(p3.range))
  m$p3.score[m$p3.score<min(p3.range)] <- 0L
  m$N <- 1 / (1 + exp(-1 * (ag + m$log_kd + c*m$ORF + p3*m$p3.score) ))
  m$log_kd <- NULL
  m$N_bg <- 1 / (1 + exp(-1 * (ag  + c*m$ORF) ))
  m <- as.data.frame(rowsum(as.matrix(m[,c("N","N_bg")]), group=m$transcript))
  m <- data.frame( transcript=row.names(m),
                   repression=log(1+exp(b)*m$N_bg) - log(1 + exp(b)*m$N) )
  
  if(!is.null(ll) && nrow(m) > 1){
    m <- merge(m,ll,by = "transcript", all.x = TRUE)
    
    # get the utr score
    m$utr_len <- log10(m$utr_len)
    m$utr_len[is.infinite(m$utr_len)] <- 0
    qu_un <- m[!duplicated(m$transcript),"utr_len"]
    qu <- quantile(qu_un, probs = c(0.05,0.95), na.rm = TRUE)
    m$utr_score <- (m$utr_len - qu[1]) / (qu[2] - qu[1])
    
    # get the orf score
    if(sum(m$orf_len) > 0){
      m$orf_len <- log10(m$orf_len)
      m$orf_len[is.infinite(m$orf_len)] <- 0
      qu_un <- m[!duplicated(m$transcript),"orf_len"]
      qu <- quantile(qu_un, probs = c(0.05,0.95), na.rm = TRUE)
      m$orf_score <- (m$orf_len - qu[1]) / (qu[2] - qu[1])
    }else{
      m$orf_score <- 0
    }
  m$repression <- m$repression + coef_utr*m$utr_score*m$repression + coef_orf*m$orf_score*m$repression
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



.datatable.aware = TRUE




# deprecated, to be removed

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



#' aggregateMatches_Biochem
#'
#' Aggregates Matches of the findSeedMatches function according to the
#' "Biochemical Model" of McGeary et al., 2020, Science
#'
#' @param e A GRanges object as produced by `findSeedMatches`.
#' @param kd_cut_off A cutoff value for log_kd values
#' @param ag The 'ag' value for the aggregation, corresponding to the free
#' concentration of AGO
#' @param keepSiteInfo An option on whether to also keep info concerning the
#' number of Binding Sites
#' @param ORF Option indicating whether sites in the open reading frame are
#' included in the scan. If they are included, a column named 'ORF' with
#' 'TRUE' / 'FALSE' entries should characterize each site.
#'
#' @return An aggregated data.frame
#' @importFrom data.table data.table as.data.table dcast
#' @importFrom GenomicRanges mcols
#' @export
aggregateMatches_Biochem <- function(e, kd_cut_off = 0, ag = -6.5, keepSiteInfo=FALSE){
  b <- 0.8655766248703003
  c <- -1.848806619644165
  m <- as.data.frame(mcols(e))
  m$transcript <- as.factor(seqnames(e))
  m <- as.data.table(m[,intersect(colnames(m), 
                                  c("transcript","miRNA","type","log_kd","ORF"))])

  if(keepSiteInfo)
    m_agg2 <- dcast( m[,.(N=.N), by=c("transcript","miRNA","type")],
                     formula=transcript+miRNA~type, value.var="N", fill=0)
  m$type <- NULL
                       
  if("ORF" %in% colnames(m)){
    m$ORF <- as.integer(m$ORF)
  }else{
    m$ORF <- 0L
  }

  m$log_kd <- m$log_kd / 1000
  m <- m[m$log_kd < kd_cut_off,]
  m$log_kd <- -m$log_kd
  m$N <- 1 / (1 + exp(-1 * (ag + m$log_kd + c*m$ORF)))
  m$log_kd <- NULL
  m$N_bg <- 1 / (1 + exp(-1 * (ag  + c*m$ORF)))
  
  m <- m[,.(N=sum(N), N_bg=sum(N_bg)), by=c("transcript","miRNA")]
  
  if(keepSiteInfo)
    m <- merge(m, m_agg2, by=c("transcript","miRNA"), all=TRUE)
  
  m$repression <- log(1 + exp(b)*m$N_bg) - log(1 + exp(b)*m$N)
  m$N <- m$N_bg <- NULL
  m
} 
