



## This functions takes site-species-abundance data and creates a
## matrix where the sites are columns and the rows are species.

samp2site.spp <- function(site, spp, abund, FUN=sum) {
  x <- tapply(abund, list(site = site, spp = spp), FUN)
  x[is.na(x)] <- 0
  return(x)
}


## write to a pdf
pdf.f <- function(f, file, ...) {
  cat(sprintf('Writing %s\n', file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}


## does the reverse of samp2site

comm.mat2sample <-  function (z) {
  temp <- data.frame(expand.grid(dimnames(z))[1:2],
                     as.vector(as.matrix(z)))
  temp <- temp[sort.list(temp[, 1]), ]
  data.frame(Site = temp[, 1], Samp = temp[, 3],
             Date = temp[, 2])
}


## function to clean up white-space in a column of data (replaces all
## instances of white-space with " " and empty cells with ""
fix.white.space <- function(d) {
    d <- as.character(d)
    remove.first <- function(s) substr(s, 2, nchar(s))
    d <- gsub("        ", " ", d, fixed=TRUE)
    d <- gsub("       ", " ", d, fixed=TRUE)
    d <- gsub("      ", " ", d, fixed=TRUE)
    d <- gsub("     ", " ", d, fixed=TRUE)
    d <- gsub("    ", " ", d, fixed=TRUE)
    d <- gsub("   ", " ", d, fixed=TRUE)
    d <- gsub("  ", " ", d, fixed=TRUE)

    tmp <- strsplit(as.character(d), " ")
    d <- sapply(tmp, function(x) paste(x, collapse=" "))

    first <- substr(d, 1, 1)
    d[first==" "] <- remove.first(d[first==" "])

    trim <- function (x) gsub("^\\s+|\\s+$", "", x)

    d <- trim(d)
    return(d)
}


## return sorted unique values
id <- function(x) unique(sort(x))
