predict.PCLasso <-
function(object, x = NULL,
    type = c("link", "response", "survival", "median", "norm", "coefficients",
    "vars", "nvars","vars.unique", "nvars.unique", "groups", "ngroups"),
    lambda, ...){

    type <- match.arg(type)
    if(type == "vars.unique"){
        vars.tmp <- predict(object = object$fit,
                            type = "vars", lambda = lambda, ...)

        if(is.list(vars.tmp)){
            vars.list <- vector(mode = "list", length = length(vars.tmp))
            names(vars.list) <- names(vars.tmp)
            for(vars.list.i in 1:length(vars.tmp)){
                if(length(vars.tmp[[vars.list.i]]) > 0){
                    vars.list[[vars.list.i]] <-
                        unique(ext2EntrezID(rownames(object$fit$beta)[vars.tmp[[vars.list.i]]]))
                }else{
                    vars.list[[vars.list.i]] <- vars.tmp[[vars.list.i]]
                }

            }
            vars.list
        }else{if(length(lambda) > 1){
                vars.vector <- rep(NA, length = length(vars.tmp))
                names(vars.vector) <- names(vars.tmp)
                for(ii in 1:length(vars.tmp)){
                    vars.vector[ii] <-
                        ext2EntrezID(rownames(object$fit$beta)[vars.tmp[ii]])
                }
                vars.vector
            }else{
                unique(ext2EntrezID(rownames(object$fit$beta)[vars.tmp]))
            }
        }

    }else if(type == "nvars.unique"){
        vars.tmp <- predict(object = object$fit,
                            type = "vars", lambda = lambda, ...)

        if(is.list(vars.tmp)){
            vars.list <- vector(mode = "list", length = length(vars.tmp))
            names(vars.list) <- names(vars.tmp)
            nvars.vector <- rep(0, length = length(vars.tmp))
            names(nvars.vector) <- names(vars.tmp)
            for(vars.list.i in 1:length(vars.tmp)){
                if(length(vars.tmp[[vars.list.i]]) > 0){
                    vars.list[[vars.list.i]] <-
                        unique(ext2EntrezID(rownames(object$fit$beta)[vars.tmp[[vars.list.i]]]))
                    nvars.vector[vars.list.i] <-
                        length(vars.list[[vars.list.i]])
                }
            }
            nvars.vector
        }else{
            if(length(lambda) > 1){
                nvars.vector <- rep(NA, length = length(vars.tmp))
                names(nvars.vector) <- names(vars.tmp)
                for(ii in 1:length(vars.tmp)){
                    nvars.vector[ii] <-
                        length(ext2EntrezID(rownames(object$fit$beta)[vars.tmp[ii]]))
                }
                nvars.vector
            }else{
                length(unique(ext2EntrezID(rownames(object$fit$beta)[vars.tmp])))
            }

        }


    }else{
        if(is.null(x)){
            predict(object = object$fit, type = type,
                    lambda = lambda, ...)
        }else{
            # extended genes
            commonFeat.ext <- unlist(object$group.dt)

            # New names of extended genes
            # The new name consists of "group+.+gene name"
            commonFeat.extName <- c()
            for(i in 1:length(object$group.dt)){
                names.i <- paste0(names(object$group.dt)[i], ".",
                                  object$group.dt[[i]])
                commonFeat.extName <- c(commonFeat.extName, names.i)
            }

            # extended dataset
            x.ext <- x[, commonFeat.ext]
            colnames(x.ext) <- commonFeat.extName

            predict(object = object$fit, X = x.ext,
                    type = type, lambda = lambda, ...)
        }
    }
}
