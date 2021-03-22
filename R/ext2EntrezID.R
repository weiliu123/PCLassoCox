ext2EntrezID <-
function(x){
    x.EntrezID <- rep(0, length = length(x))
    for(i in 1:length(x)){
        x.EntrezID[i] <- unlist(strsplit(x[i], split = "[.]"))[2]
    }
    x.EntrezID
}
