#' Transform selection features into EntrezID
#'
#' @description
#' Transform selection features into EntrezID at the gene level.
#' @param x Names of selected vars by \code{PCLasso} or \code{cv.PCLasso}.
#'
#' @return Selected genes.
#' @export
ext2EntrezID <-
function(x){
    x.EntrezID <- rep(0, length = length(x))
    for(i in 1:length(x)){
        x.EntrezID[i] <- unlist(strsplit(x[i], split = "[.]"))[2]
    }
    x.EntrezID
}
