


calcPcoa <- function(par.dat, parasites, nperm,
                     sites){
    ## function for calculating pcoas by bee species or site, running
    ## permutations to determine dignificance. par.dat is the main
    ## dataset of parasite and pathogen exposure for each individual,
    ## path.or.parasites is a vector of the columns of either the
    ## parasites or pathogens, nperm is the number of permutations to
    ## determine significance, genSp is the vector of bee genuses,
    ## sites is a vector of sites
    com <- par.dat[, parasites]
    no.p <- rowSums(com) == 0 | apply(com, 1, function(x) any(is.na(x)))
    com <- com[!no.p,]
    sites <- as.character(sites)[!no.p]

    com.dist  <- vegdist(com, method="jaccard")

    ## are the site different in parasite composition?
    beta.disper.site <- betadisper(com.dist, sites,
                                   type="centroid")
    perm.test.site <- permutest(beta.disper.site,
                                control = permControl(nperm = nperm),
                                pairwise = TRUE)

    return(list(tests= perm.test.site,
                dists=list(dist=com.dist, sites=sites)))
}


plotCommDist  <- function(dist.mat, sites, par.or.path){
    ## plotting function for pcoas
    f.pcoa <- function(){
        cols <- rainbow(length(unique(sites)))
        names(cols) <- unique(sites)

        dist.mat <- as.matrix(dist.mat)
        pcoa.comm <- cmdscale(dist.mat)

        plot(NA, asp=1,  cex=1.5,
             ylim=range(pcoa.comm[,2]),
             xlim=range(pcoa.comm[,1]),
             xlab='',
             ylab='',
             xaxt='n',
             yaxt='n',
             cex.lab=1.5)
        for(site in unique(sites)){
            points(jitter(pcoa.comm[sites == site,]),
                   col=cols[site], pch=16, cex=1.5)
            ## points(pcoa.comm[sites == site,],
            ##        col="black", pch=1, cex=1.5)
        }
        ordihull(pcoa.comm, sites)

        legend("topright", legend=unique(sites),
               bty="n", cex=0.6, col=cols[unique(sites)],
               pch=16)
        legend("topright", legend=unique(sites),
               bty="n", cex=0.6, col="black", pch=1)

        mtext('PCoA1', 1, line=2, cex=1.5)
        mtext('PCoA2', 2, line=2, cex=1.5)
    }
    ## function for plotting PcoA axes
    path <- 'figures'
    pdf.f(f.pcoa, file= file.path(path,
                                  sprintf("%s_pcoa.pdf", par.or.path)),
          width=7, height=7)
}
