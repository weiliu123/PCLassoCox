#' Make predictions from a PCLasso model
#'
#' @description Similar to other predict methods, this function returns
#' predictions from a fitted \code{PCLasso} object.
#'
#' @param object Fitted \code{PCLasso} model object.
#' @param x Matrix of values at which predictions are to be made. The features
#' (genes) contained in \code{x} should be consistent with those contained in
#' \code{x} in the \code{PCLasso} function.  Not used for type="coefficients"
#' or for some of the type settings in \code{predict}.
#' @param type Type of prediction: "link" returns the linear predictors;
#'   "response" gives the risk (i.e., exp(link)); "vars" returns the indices for
#'   the nonzero coefficients; "vars.unique" returns unique features (genes)
#'   with nonzero coefficients (If a feature belongs to multiple groups and
#'   multiple groups are selected, the feature will be repeatedly selected.
#'   Compared with "var", "var.unique" will filter out repeated features.);
#'   "groups" returns the groups with at least one nonzero coefficient; "nvars"
#'   returns the number of nonzero coefficients; "nvars.unique" returens the
#'   number of unique features (genes) with nonzero coefficients; "ngroups"
#'   returns the number of groups with at least one nonzero coefficient; "norm"
#'   returns the L2 norm of the coefficients in each group."survival" returns
#'   the estimated survival function; "median" estimates median survival times.
#' @param lambda Values of the regularization parameter \code{lambda} at which
#'   predictions are requested. For values of \code{lambda} not in the sequence
#'   of fitted models, linear interpolation is used.
#' @param ... Arguments to be passed to \code{predict.grpsurv} in the R package
#' \code{grpreg}.
#' @details
#' See \code{predict.grpsurv} in the R package \code{grpreg} for details.
#' @return The object returned depends on \code{type}.
#' @seealso \code{\link{PCLasso}}
#' @export
#'
#' @examples
#' # load data
#' data(GBM)
#' data(PCGroup)
#'
#' fit1 <- PCLasso(x = GBM$GBM.train$Exp, y = GBM$GBM.train$survData, group =
#' PCGroup)
#'
#' # predict risk scores of samples in x.test
#' s <- predict(object = fit1, x = GBM$GBM.test$Exp, type="link",
#' lambda=fit1$fit$lambda)
#'
#' s <- predict(object = fit1, x = GBM$GBM.test$Exp, type="link",
#' lambda=fit1$fit$lambda[10])
#'
#' s <- predict(object = fit1, x = GBM$GBM.test$Exp, type="link", lambda=c(0.1,
#' 0.01))
#'
#' # Nonzero coefficients
#' sel.groups <- predict(object = fit1, type="groups",
#'                       lambda = fit1$fit$lambda)
#' sel.ngroups <- predict(object = fit1, type="ngroups",
#'                        lambda = fit1$fit$lambda)
#' sel.vars.unique <- predict(object = fit1, type="vars.unique",
#'                           lambda = fit1$fit$lambda)
#' sel.nvars.unique <- predict(object = fit1, type="nvars.unique",
#'                             lambda = fit1$fit$lambda)
#' sel.vars <- predict(object = fit1, type="vars",
#'                     lambda=fit1$fit$lambda)
#' sel.nvars <- predict(object = fit1, type="nvars",
#'                      lambda=fit1$fit$lambda)
#'
#' # For values of lambda not in the sequence of fitted models,
#' # linear interpolation is used.
#' sel.groups <- predict(object = fit1, type="groups",
#'                       lambda = c(0.1, 0.01))
#' sel.ngroups <- predict(object = fit1, type="ngroups",
#'                        lambda = c(0.1, 0.01))
#' sel.vars.unique <- predict(object = fit1, type="vars.unique",
#'                            lambda = c(0.1, 0.01))
#' sel.nvars.unique <- predict(object = fit1, type="nvars.unique",
#'                             lambda = c(0.1, 0.01))
#' sel.vars <- predict(object = fit1, type="vars",
#'                     lambda=c(0.1, 0.01))
#' sel.nvars <- predict(object = fit1, type="nvars",
#'                      lambda=c(0.1, 0.01))
#'
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
