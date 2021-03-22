cv.PCLasso <-
function(x, y, group, 
    penalty = c("grLasso", "grMCP", "grSCAD","gel", "cMCP"), nfolds = 5,
    standardize = TRUE,...){
    penalty <- match.arg(penalty)
    
    if(standardize){
        x <- scale(x, center = TRUE, scale = TRUE)			
    }
    
    # feature set in all groups
    featureSet <- unique(unlist(group))
    
    # common features in groups and expression matrix x
    commonFeat <- intersect(colnames(x), featureSet)
    
    # filter undetected genes in expression matrix x
    x <- x[,commonFeat]
    
    # filter undetected genes in groups
    # Construct groups whose expressions are detected
    group.dt <- vector(mode = "list", length = 0)
    idx <- 0
    for(i in 1:length(group)){
        group.i <- intersect(group[[i]], commonFeat)
        if(length(group.i) > 0){
            idx <- idx + 1
            group.dt[[idx]] <- group.i
            names(group.dt)[idx] <- names(group)[i]
        }
    }
    
    # Filter duplicate groups (generated due to undetected genes) 
    group.dt <- group.dt[!duplicated(group.dt)]
    
    # extended genes
    commonFeat.ext <- unlist(group.dt)
    
    # New names of extended genes
    # The new name consists of "group+.+gene name" 
    commonFeat.extName <- c()
    for(i in 1:length(group.dt)){
        names.i <- paste0(names(group.dt)[i], ".", group.dt[[i]])
        commonFeat.extName <- c(commonFeat.extName, names.i)
    }
    
    # group of extended genes
    groupOfFeats <- c()
    for(i in 1:length(group.dt)){
        group.i <- rep(names(group.dt)[i], length = length(group.dt[[i]]))
        groupOfFeats <- c(groupOfFeats, group.i)
    }
    
    # extended dataset
    x.ext <- x[, commonFeat.ext]
    colnames(x.ext) <- commonFeat.extName
    
    # grpsurv
    cv.fit <- cv.grpsurv(X=x.ext,
                      y=y,
                      group = groupOfFeats,
                      penalty = penalty, 
                      nfolds = nfolds,...)
    
    res <- list(cv.fit = cv.fit, group.dt = group.dt)
    class(res) <- "cv.PCLasso"
    
    return(res)
}
