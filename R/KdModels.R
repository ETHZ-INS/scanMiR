#' @exportClass KdModel
setClass(
  "KdModel",
  contains="list",
  validity=function(object){
    for(f in c("name","canonical.seed","mirseq")){
      if(is.null(object[[f]]) || !is.character(object[[f]]) || length(object[[f]])!=1)
        stop("The model should have a `",f,"` slot (character of length 1).")
    }
    if(is.null(object$mer8) || length(object$mer8) != 1024 || !is.integer(object$mer8)){
      stop("The `mer8` slot should be an integer vector of length 1024.")
    }
    if(is.null(object$fl) || length(object$fl) != 1024 || !is.integer(object$fl)){
      stop("The `fl` slot should be an integer vector of length 1024.")
    }
  }
)

#' @export
setMethod("show", "KdModel", function(object){
  cat(paste0("A `KdModel` for ", object$name, "\n  Sequence: ", object$mirseq,
             "\n  Canonical seed: ", object$canonical.seed))
})

#' @export
setMethod("summary", "KdModel", function(object){
  c( name=object$name, sequence=object$mirseq, canonical.seed=object$canonical.seed )
})

getKdModel <- function(kd, mirseq=NULL, name=NULL, conservation=NA, ...){
  if(is.character(kd) && length(kd)==1){
    if(is.null(name)) name <- gsub("\\.txt$","",gsub("_kds","",basename(kd)))
    kd <- read.delim(kd, header=TRUE, stringsAsFactors=FALSE)[,c(1,2,4)]
  }
  if(is.null(mirseq)) mirseq <- as.character(kd$mirseq[1])
  if(is.null(name)) name <- as.character(kd$mir[1])
  kd <- kd[,c("X12mer","log_kd")]
  seed <- paste0(as.character(reverseComplement(DNAString(substr(mirseq, 2,8)))),"A")
  w <- grep("X|N",kd$X12mer)
  pwm <- Biostrings::consensusMatrix(
    as.character(rep(kd$X12mer[-w], floor( (exp(-kd$log_kd[-w]))/3 ))),
    as.prob=TRUE
  )
  fields <- c("mer8","fl.score")
  if(!all(fields %in% colnames(kd))) kd <- .prep12mers(kd, seed=seed)
  fields <- c(fields, "log_kd")
  if(!all(fields %in% colnames(kd))) stop("Malformed `kd` data.frame.")
  co <- t(sapply(split(kd[,c("log_kd","fl.score")], kd$mer8), FUN=function(x){
    .lm.fit(cbind(1,x$fl.score),x$log_kd)$coefficients
  }))
  fitted <- co[kd$mer8,1]+co[kd$mer8,2]*kd$fl.score
  new("KdModel", list(mer8=as.integer(round(co[,1]*1000)), 
                      fl=as.integer(round(co[,2]*1000)), 
                      name=name, mirseq=mirseq, canonical.seed=seed,
                      pwm=pwm, conservation=conservation,
                      cor=cor(fitted, kd$log_kd),
                      mae=median(abs(kd$log_kd-fitted)), ... ))
}

.prep12mers <- function(x, seed){
  if(is.data.frame(x)){
    if(!all(c("X12mer","log_kd") %in% colnames(x)))
      stop("`x` should be a character vector or a data.frame with the columns ",
           "'X12mer' and 'log_kd'")
    x <- x[grep("N|X",substr(x$X12mer, 3,10),invert=TRUE),]
    x <- cbind(x[,"log_kd",drop=FALSE], .prep12mers(x$X12mer, seed=seed))
    return(x[!is.na(x$mer8),])
  }
  x <- gsub("X","N",as.character(x))
  y <- .getFlankingScore(x)
  data.frame(mer8=as.integer(factor(substr(x, 3,10), levels=getSeed8mers(seed))), 
             fl.score=y$score, fl.ratio=y$ratio)
}

.getFlankingScore <- function(x){
  fl.s <- matrix(c(-0.24, -0.14, 0, 0.1, 0.28, -0.24, -0.3, 0, 0.13, 0.42, -0.075, -0.18, 
                   0, 0, 0.25, -0.1, -0.1, 0, 0, 0.26), 
                 nrow=5, dimnames=list(c("A","T","N","C","G")))
  fl.m <- cbind(substr(x,1,1), substr(x,2,2), substr(x,11,11),substr(x,12,12))
  fl.m <- matrix(as.integer(factor(fl.m, row.names(fl.s))), ncol=4)
  fl.score <- vapply(1:4, FUN.VALUE=numeric(length(x)), 
                     FUN=function(i) fl.s[fl.m[,i,drop=FALSE],i,drop=FALSE])
  if(is.null(dim(fl.score))) fl.score <- matrix(fl.score, ncol=4)
  fl.score <- rowSums(fl.score)
  fl.ratio <- rowSums( (fl.m-3)>0 ) - rowSums( (fl.m-3)<0 )
  return( list(score=fl.score, ratio=fl.ratio) )
}

.add8merN <- function(mod, mer8=NULL){
  if(is.null(mer8)) mer8 <- getSeed8mers(mod$canonical.seed, addNs=TRUE)
  i1 <- split(1:1024, substr(mer8[1:1024],2,8))
  i2 <- split(1:1024, substr(mer8[1:1024],1,7))
  fl <- co <- rep(0,416)
  m1 <- grep("^N",mer8)
  m2 <- grep("N$",mer8)
  co[m1-1024] <- sapply(i1[gsub("N","",mer8[m1])], FUN=function(i) median(mod$mer8[i]))
  fl[m1-1024] <- sapply(i1[gsub("N","",mer8[m1])], FUN=function(i) median(mod$fl[i]))
  co[m2-1024] <- sapply(i2[gsub("N","",mer8[m2])], FUN=function(i) median(mod$mer8[i]))
  fl[m2-1024] <- sapply(i2[gsub("N","",mer8[m2])], FUN=function(i) median(mod$fl[i]))
  mod$mer8 <- c(mod$mer8, co)
  mod$fl <- c(mod$fl, fl)
  mod
}

#' assignKdType
#'
#' Assigns a log_kd and match type to a set of matched sequences.
#'
#' @param x A vector of matched sequences, each of 12 nucleotides
#' @param mod An object of class `KdModel`
#' @param mer8 The optional set of 8mers included in the model (if omitted, will be 
#' reconstructed from the model)
#'
#' @return A data.frame with one row for each element of `x`, and the columns `type` and
#' `log_kd`. To save space, the reported log_kd is multiplied by 1000, rounded and saved
#' as an integer.
#' @export
assignKdType <- function(x, mod, mer8=NULL){
  if(is.null(mer8)) mer8 <- getSeed8mers(mod$canonical.seed, addNs=TRUE)
  mod <- .add8merN(mod, mer8)
  fl.score <- as.numeric(.getFlankingScore(x)$score)
  mer8 <- factor(as.character(subseq(x, 3,10)), levels=mer8)
  data.frame(type=getMatchTypes(levels(mer8), mod$canonical.seed)[as.integer(mer8)],
             log_kd=as.integer(round(mod$mer8[mer8] + fl.score*mod$fl[mer8])))
}
