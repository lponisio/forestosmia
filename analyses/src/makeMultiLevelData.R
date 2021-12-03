
makeDataMultiLevel <- function(indiv.data, site.col){
    ##
    site.ids <- unlist(tapply(indiv.data[, site.col],
                              indiv.data[, site.col],
                              function(x) 1:length(x)))
    names(site.ids) <- NULL
    indiv.data$SiteIDs <- site.ids
    indiv.data$Weights <- indiv.data$SiteIDs
    indiv.data$Weights[indiv.data$Weights > 1] <- 0
    return(indiv.data)
}
