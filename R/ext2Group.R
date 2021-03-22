ext2Group <-
function(x){
    x.Group <- rep(0, length = length(x))
    for(i in 1:length(x)){
        x.Group[i] <- unlist(strsplit(x[i], split = "[.]"))[1]
    }
    unique(x.Group)
}
