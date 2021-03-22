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