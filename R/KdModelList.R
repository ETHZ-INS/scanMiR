#' @exportClass KdModelList
setClass(
  "KdModelList",
  contains="list",
  validity=function(object){
    if(!all(sapply(object,FUN=function(x) is(x,"KdModel"))))
      stop("A KdModelList should be a list of objects of class 'KdModel'")
    return(TRUE)
  }
)

KdModelList <- function(x, description=NULL){
  names(x) <- vapply(x, FUN.VALUE=character(1), FUN=function(x) x$name)
  x <- new("KdModelList", x)
  attr(x, "created") <- Sys.Date()
  if(!is.null(description)) attr(x, "description") <- description
  x
}

#' @export
setMethod("as.data.frame", "KdModelList", function(x){
  x <- t(sapply(x, summary))
  as.data.frame(x[,-1,drop=FALSE], row.names=x[,1])
})

#' @export
setMethod("summary", "KdModelList", function(object){
  d <- attr(object, "created")
  cat(paste0("A `KdModelList` object", 
            ifelse(is.null(d), "", paste0(" created on ",d,",\n")),
            " containing binding affinity models from ", length(object), " miRNAs.\n"))
  if(!is.null(desc <- attr(object, "description"))) cat(desc)
  cons <- conservation(object)
  if(!all(is.na(cons))){
    print(table(conservation(object)))
  }
})

  #' conservation
#'
#' @param x A KdModelList, or a KdModel
#'
#' @return A vector of the conservation status for each miRNA
#' @export
conservation <- function(x){
  lvls <- c("-1"="Low-confidence","0"="Poorly conserved","1"="Conserved across mammals",
            "2"="Conserved across vertebrates")
  if(is(x,"KdModelList")){
    y <- factor(sapply(x, FUN=function(x){
      if(is.null(x$conservation)) return(NA_integer_)
      x$conservation
    }), levels=names(lvls))
  }else if(is(x,"KdModel")){
    if(is.null(x$conservation)) return(NA)
    y <- factor(x$conservation, levels=names(lvls))
  }else{
    stop("Undefined for an object of class ", class(x))
  }
  levels(y) <- as.character(lvls)
  y
}
