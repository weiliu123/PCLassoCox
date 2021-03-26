#' Convert selected features into ComplexID
#'
#' Convert selected features into ComplexID at protein complex level
#' @param x Names of selected vars by \code{PCLasso} or \code{cv.PCLasso}.
#'
#' @return Selected protein complexes.
#' @export
#'
ext2Group <-
function(x){
    x.Group <- rep(0, length = length(x))
    for(i in 1:length(x)){
        x.Group[i] <- unlist(strsplit(x[i], split = "[.]"))[1]
    }
    unique(x.Group)
}
