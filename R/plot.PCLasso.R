#' Plot coefficients from a PCLasso object
#'
#' @description Produces a plot of the coefficient paths for a fitted
#' \code{PCLasso} object.
#' @param x Fitted \code{PCLasso} model.
#' @param norm If TRUE, plot the norm of each group, rather than the individual
#' coefficients.
#' @param ... Other graphical parameters to \code{plot}.
#'
#' @seealso \code{\link{PCLasso}}
#' @export
#'
#' @examples
#' # load data
#' data(GBM)
#' data(PCGroup)
#'
#' x = GBM$GBM.train$Exp
#' y = GBM$GBM.train$survData
#'
#' # fit the PCLasso model
#' fit1 <- PCLasso(x, y, group = PCGroup, penalty = "grLasso")
#'
#' # plot the norm of each group
#' plot(fit1, norm = TRUE)
#'
#' # plot the individual coefficients
#' plot(fit1, norm = FALSE)
plot.PCLasso <-
function(x, norm = TRUE, ...){
    plot(x$fit, norm = norm, ...)
}
