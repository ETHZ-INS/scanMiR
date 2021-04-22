#' @exportClass KdModelList
#' @rdname KdModelList
setClass(
  "KdModelList",
  contains="list",
  validity=function(object){
    if(!all(vapply(object, is, logical(1), class2="KdModel")))
      stop("A KdModelList should be a list of objects of class 'KdModel'")
    if(any(duplicated(names(object))))
      stop("There are multiple KdModels with the same name!")
    mn <- as.character(unlist(lapply(object, FUN=function(x) x$name)))
    if(!identical(names(object),mn)) stop("Mismatched names!")
    return(TRUE)
  }
)

#' KdModelList
#'
#' @param ... Any number of \code{\link{KdModel}} objects or lists thereof.
#' @param description A description for the collection.
#' @param makeUnique Logical; whether to rename models if names are duplicated.
#'
#' @return A KdModelList
#' @export
#' @examples
#' data(SampleKdModel)
#' mods <- KdModelList(SampleKdModel, SampleKdModel, makeUnique = TRUE)
#' mods
KdModelList <- function(..., description=NULL, makeUnique=FALSE){
  x <- list(...)
  isKdm <- vapply(x, is, logical(1), class2="KdModel")
  if(any(!isKdm)){
    warn <- FALSE
    for(f in which(!isKdm)){
      y <- x[[f]]
      x <- x[-f]
      if(is(y,"KdModelList") ||
         (is.list(y) && vapply(y, is, logical(1), class2="KdModel"))){
        x <- c(x,as.list(y))
      }else{
        warn <- TRUE
      }
    }
    if(warn) warning("Some objects were not KdModels and were discarded.")
  }
  nn <- vapply(x, FUN.VALUE=character(1), FUN=function(x) x$name)
  if(length(wdup <- which(duplicated(nn)))>0){
    if(makeUnique){
      wdup <- which(nn %in% unique(nn[wdup]))
      nn <- make.unique(nn)
      for(i in wdup) x[[i]]$name <- nn[i]
    }else{
      stop("There are multiple KdModels with the same name. ",
           "Use `KdModelList(..., makeUnique=TRUE)` to automatically rename ",
           "them.")
    }
  }
  names(x) <- nn
  x <- new("KdModelList", x)
  if(!is.null(attr(x, "created"))) attr(x, "created") <- Sys.Date()
  if(!is.null(description)) attr(x, "description") <- description
  x
}

#' Methods for the \code{\link{KdModelList}} classes
#' @name KdModelList-methods
#' @aliases KdModelList-methods
#' @seealso \code{\link{KdModel}}, \code{\link{KdModelList}}
#' @param object,x An object of class \code{\link{KdModelList}}
#' @return Depends on the method.
#' @examples
#' # create a KdModelList :
#' data(SampleKdModel)
#' kml <- KdModelList( SampleKdModel, SampleKdModel, makeUnique=TRUE )
#'
#' summary(kml)
#' kml[1] # returns a KdModelList
#' kml[[2]] # returns a KdModel
#' conservation(kml)
#' @rdname KdModelList-methods
#' @export
setMethod("summary", "KdModelList", function(object){
  d <- attr(object, "created")
  cat(paste0("A `KdModelList` object",
            ifelse(is.null(d), "", paste0(" created on ",d,",\n")),
            " containing binding affinity models from ", length(object),
            " miRNAs.\n"))
  if(!is.null(desc <- attr(object, "description"))) cat(paste0(desc,"\n"))
  cons <- conservation(object)
  if(!all(is.na(cons))){
    print(table(conservation(object)))
  }
})

#' @rdname KdModelList-methods
#' @export
#' @param i the index of item(s) to select
#' @param j,drop,... ignored
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
    y <- factor(vapply(x, FUN.VALUE=integer(1), FUN=function(x){
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

