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
aggregateSites <- function(m, ag=-5.5, b=0.8656, c=-1.8488, p3=NULL, toInt=FALSE, BP=NULL){
  if(is.null(BP)) BP <- BiocParallel::SerialParam()
  if(is(m,"GRanges")){
    m$transcript <- as.factor(seqnames(m))
    m <- mcols(m)
    m$miRNA <- as.factor(m$miRNA)
    m <- as.data.frame(m)
  }
  m <- m[,c("miRNA","transcript","ORF","log_kd",intersect("align.3p",colnames(m)))]
  m <- split(m, m$miRNA)
  m <- bplapply(m, BPPARAM=BP, FUN=function(x){
    .aggregate_miRNA(m, ag=ag, b=b, c=c, p3=p3, toInt=toInt)
  })
  dplyr::bind_rows(m, .id="miRNA")
}


.aggregate_miRNA <- function(m, ag=-5.5, b=0.8656, c=-1.8488, p3=NULL, toInt=FALSE){
  if(is.null(m$ORF)) m$ORF <- 0L
  m <- m[,c("transcript","ORF","log_kd",intersect("align.3p",colnames(m)))]
  m$ORF <- as.integer(m$ORF)
  m$log_kd <- m$log_kd / 1000
  m <- m[m$log_kd < 0,]
  m$log_kd <- -m$log_kd
  m$N <- 1 / (1 + exp(-1 * (ag + m$log_kd + c*m$ORF)))
  m$log_kd <- NULL
  m$N_bg <- 1 / (1 + exp(-1 * (ag  + c*m$ORF)))
  m <- rowsum(m[,c("N","N_bg")], group=m$transcript)
  m <- data.frame( transcript=row.names(m),
                   repression=log(1+exp(b)*m$N_bg) - log(1 + exp(b)*m$N) )
  if(toInt) m$repression <- as.integer(round(1000*m$repression))
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
