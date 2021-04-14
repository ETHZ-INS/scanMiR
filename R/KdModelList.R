#' @exportClass KdModelList
setClass(
  "KdModelList",
  contains="list",
  validity=function(object){
    if(!all(unlist(lapply(object,FUN=function(x) is(x,"KdModel")))))
      stop("A KdModelList should be a list of objects of class 'KdModel'")
    return(TRUE)
  }
)

#' KdModelList
#'
#' @param x A list of KdModels
#' @param description A description for the collection.
#'
#' @return A KdModelList
#' @export
#' @examples
#' data(SampleKdModel)
#' KdModelList(list(SampleKdModel,SampleKdModel))
KdModelList <- function(x, description=NULL){
  names(x) <- vapply(x, FUN.VALUE=character(1), FUN=function(x) x$name)
  x <- new("KdModelList", x)
  if(!is.null(attr(x, "created"))) attr(x, "created") <- Sys.Date()
  if(!is.null(description)) attr(x, "description") <- description
  x
}

#' @export
setMethod("summary", "KdModelList", function(object){
  d <- attr(object, "created")
  cat(paste0("A `KdModelList` object", 
            ifelse(is.null(d), "", paste0(" created on ",d,",\n")),
            " containing binding affinity models from ", length(object), " miRNAs.\n"))
  if(!is.null(desc <- attr(object, "description"))) cat(paste0(desc,"\n"))
  cons <- conservation(object)
  if(!all(is.na(cons))){
    print(table(conservation(object)))
  }
})

#' @export
setMethod("[", "KdModelList", function(x, i, j=NULL, ..., drop = TRUE){
  xo <- new("KdModelList", unclass(x)[i])
  if(!is.null(attr(x, "created"))) attr(xo, "created") <- attr(x, "created")
  if(!is.null(attr(x, "description")))
    attr(xo, "description") <- attr(x, "description")
  xo
})


#' conservation
#'
#' @param x A KdModelList, or a KdModel
#'
#' @return A vector of the conservation status for each miRNA
#' @export
#' @examples 
#' data(SampleKdModel)
#' conservation(SampleKdModel)
conservation <- function(x){
  lvls <- .conservation_levels()
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

.conservation_levels <- function(){
  c("-1"="Low-confidence", "0"="Poorly conserved",
    "1"="Conserved across mammals", "2"="Conserved across vertebrates")
}

