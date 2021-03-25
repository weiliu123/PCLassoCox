#' Plot the cross-validation curve from a \code{cv.PCLasso} object
#'
#' @description Plot the cross-validation curve from a \code{cv.PCLasso} object,
#' along with standard error bars.
#' @param x Fitted \code{cv.PCLasso} model.
#' @param type What to plot on the vertical axis. "cve" plots the
#' cross-validation error (deviance); "rsq" plots an estimate of the fraction of
#' the deviance explained by the model (R-squared); "snr" plots an estimate of
#' the signal-to-noise ratio; "all" produces all of the above.
#' @param norm If TRUE, plot the norm of each group, rather than the individual
#' coefficients.
#' @param ... Other graphical parameters to \code{plot}
#' @details Error bars representing approximate +/- 1 SE (68\% confidence
#' intervals) are plotted along with the estimates at value of lambda. See
#' \code{plot.cv.grpreg} in the R package \code{grpreg} for details.
#'
#' @export
#' @seealso \code{\link{cv.PCLasso}}
#' @examples
#' # load data
#' data(GBM)
#' data(PCGroup)
#' x = GBM$GBM.train$Exp
#' y = GBM$GBM.train$survData
#'
#' cv.fit1 <- cv.PCLasso(x, y, group = PCGroup, penalty = "grLasso",nfolds = 10)
#'
#' # plot the norm of each group
#' plot(cv.fit1, norm = TRUE)
#'
#' # plot the individual coefficients
#' plot(cv.fit1, norm = FALSE)
#'
#' # plot the cross-validation error (deviance)
#' plot(cv.fit1, type = "cve")
plot.cv.PCLasso <-
function(x, type = c("cve", "rsq", "snr", "all"),
    norm = NULL, ...){
    if(is.null(norm)){
        type <- match.arg(type)
        plot(x$cv.fit, type = type, ...)
    }else{
        plot(x$cv.fit$fit, norm = norm, ...)
    }
}
