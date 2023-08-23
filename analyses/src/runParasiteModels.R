
runParasiteModels <- function(spec.data,
                              species.group, parasite,
                              xvars,
                              iter = 10^4,
                              chains = 1,
                              thin=1,
                              init=0,
                              SEM = TRUE){

  if(SEM){
    formula.parasite  <- as.formula(paste(
      paste(parasite, "| weights(WeightsPar)"),
      paste(xvars,
            collapse=" + "),
      sep=" ~ "))

    bf.par <- bf(formula.parasite, family="bernoulli")  

    bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.par +
      set_rescor(FALSE)

  } else{

    spec.data <- spec.data[spec.data$WeightsPar == 1,]
    formula.parasite  <- as.formula(paste(
      paste(parasite),
      paste(xvars,
            collapse=" + "),
      sep=" ~ "))

    bform <- bf(formula.parasite, family="bernoulli")  

  }
  
  fit.parasite <- brm(bform, spec.data,
                      cores=ncores,
                      iter = iter,
                      chains = chains,
                      thin=thin,
                      init=init,
                      control = list(adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE))
  
  write.ms.table(fit.parasite,
                 sprintf("parasitism_%s_%s",
                         species.group, parasite))

  r2 <- bayes_R2(fit.parasite)
  print(r2)

  save(fit.parasite, spec.data, r2,
       file=sprintf("saved/parasiteFit_%s_%s.Rdata",
                    species.group, parasite))
  
  plot.res(fit.parasite,  sprintf("%s_%s",
                                  species.group, parasite))
  
  return(list(fit=fit.parasite, r2=r2))
}




runCombinedParasiteModels <- function(spec.data,
                                      species.group,
                                      parasites,
                                      xvars,
                                      iter = 10^4,
                                      chains = 1,
                                      thin=1,
                                      init=0){

  browser()
  bf.parasite.formulas <- vector(mode="list",
                                 lenght=length(parasites))
  names(bf.parasite.formulas) <- parasites
  
  for(parasite in parasites){
    formula.parasite  <- as.formula(paste(
      paste(parasite, "| weights(WeightsPar)"),
      paste(xvars,
            collapse=" + "),
      sep=" ~ "))
    bf.parasite.formulas[[parasite]] <-  bf(formula.parasite, family="bernoulli")  
  }

  bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv +
    bf.parasite.formulas[[1]]+
    bf.parasite.formulas[[2]]+ bf.parasite.formulas[[3]] +
    set_rescor(FALSE)

  fit.parasite <- brm(bform, spec.data,
                      cores=ncores,
                      iter = iter,
                      chains = chains,
                      thin=thin,
                      init=init,
                      control = list(adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE))

  write.ms.table(fit.parasite,
                 sprintf("parasitism_%s_%s",
                         species.group, parasite))

  r2 <- bayes_R2(fit.parasite)
  print(r2)

  save(fit.parasite, spec.data, r2,
       file=sprintf("saved/parasiteFit_%s_%s.Rdata",
                    species.group, parasite))

  plot.res(fit.parasite,  sprintf("%s_%s",
                                  species.group, parasite))

  return(list(fit=fit.parasite, r2=r2))
}
