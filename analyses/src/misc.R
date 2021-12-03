
pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3],
                  alpha=alpha))
}


standardize <- function(x)
(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

unstandardize <- function(x, orig){
    (x*sd(orig, na.rm=TRUE)) + mean(orig, na.rm=TRUE)
}




standardize.axis <- function(x, orig)
(x-mean(orig, na.rm=TRUE))/sd(orig, na.rm=TRUE)
