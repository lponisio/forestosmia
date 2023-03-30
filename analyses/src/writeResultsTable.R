write.ms.table <- function(mod.output, mod.name){
    ci.hdi <- hdi(mod.output)
    sum.mod <- as.data.frame(round(summary(mod.output)$fixed,2))
    sum.ci.hdi <- data.frame(Parameter=ci.hdi$Parameter,
                             HDI.lb=round(ci.hdi$CI_low, 2),
                             HDI.ub=round(ci.hdi$CI_high, 2))
    coeffs <- c(paste0("b_",
                       rownames(sum.mod)),
                paste0("bs_",
                       rownames(sum.mod)))
    samps.mod <- posterior_samples(mod.output)
    coeffs <- coeffs[coeffs %in% colnames(samps.mod)]
    samps.mod <- samps.mod[, coeffs]
    coeff.samps <- colnames(samps.mod)
    coeff.samps.sum <- sub("[a-z]*_", "", coeff.samps)
    sum.ci.hdi$Parameter <- sub("[a-z]*_", "",
                                sum.ci.hdi$Parameter)
    samps.mod <- samps.mod[order(match(coeff.samps.sum,
                                       rownames(sum.mod)))]
    sum.mod$Pgt0  <- round(apply(samps.mod, 2, function(x)
        sum(x > 0)/length(x)), 2)
    sum.mod$Plt0  <- round(apply(samps.mod, 2, function(x)
        sum(x < 0)/length(x)),2)
    sum.mod$Parameter <- rownames(sum.mod)
    rownames(sum.mod) <- NULL
    sum.mod <- merge(sum.mod, sum.ci.hdi)

    sum.mod <- sum.mod[, c("Parameter", "Estimate",  "Est.Error",
                           "HDI.lb", "HDI.ub", "Rhat", "Bulk_ESS",
                           "Tail_ESS", "Pgt0", "Plt0")]
    write.table(sum.mod,
                file=sprintf("saved/tables/%s.txt", mod.name),
                sep="&", row.names=FALSE)

    write.csv(sum.mod,
              file=sprintf("saved/tables/%s.csv", mod.name), row.names=FALSE)
}
