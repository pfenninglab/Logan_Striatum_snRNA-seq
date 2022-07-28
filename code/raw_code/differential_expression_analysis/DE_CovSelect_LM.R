#expr vector
#diagnosis vector
#coefficient dataframe
#need cpm transformed

Run_LM_CovSelect = function(x = list(expressionMatrix = dat,
                                     Main = covariables.main,
                                     Covariables = covariables.potential,
                                     MainContrasts = NULL),
                            GeneSymbol = GeneSymbol,
                            parallel.cores = 20,
                            filename = "OUD_comorbid_NAc",
                            out.dir = job.dir,
                            n.permutation = 100){

  expressionMatrix = x$expressionMatrix
  Main = x$Main
  Covariables = x$Covariables
  MainContrasts = x$MainContrasts
  GeneIndeces = 1:nrow(expressionMatrix)

  VariableListOne = combn(colnames(Covariables),1)
  VariableListTwo = combn(colnames(Covariables),2)

  library(snowfall)
  set.seed(15213)
  StartTime = Sys.time()
  #parallel computing
  sfInit(parallel=TRUE,type="SOCK", cpus=parallel.cores)
  sfExport("StartTime")
  sfExport("filename")
  sfExport("expressionMatrix")
  sfExport("Main")
  sfExport("Covariables")
  sfExport("VariableListOne")
  sfExport("VariableListTwo")
  sfExport("GeneIndeces")
  sfExport("MainContrasts")
  sfExport("bestModelSelection_LM")

  indLM<-c()
  indLM<-sfLapply(GeneIndeces,function(x) try(bestModelSelection_LM(x,expressionMatrix,covariables3 = Covariables,
                                                                    Main,VariableListOne,VariableListTwo, MainContrasts),silent=TRUE ) )
  # indLM.vars = lapply(indLM, `[[`, 2)
  # indLM.race = sapply(indLM.vars, function(a){"Race"%in%a})
  # which(indLM.race)[1]
  # indLM[[20]]

  sfStop()
  save(indLM,file=paste0(out.dir, "/LM_true", filename, ".rdata"))

  a.null.dir = paste0(out.dir, "/null_", filename)
  dir.create(a.null.dir, recursive = TRUE)

  B<-100
  sfInit(parallel=TRUE,type="SOCK", cpus=parallel.cores)
  #permutation to calculate the Corrected.pvalue
  sfExport("StartTime")
  sfExport("filename")
  sfExport("expressionMatrix")
  sfExport("Main")
  sfExport("Covariables")
  sfExport("VariableListOne")
  sfExport("VariableListTwo")
  sfExport("GeneIndeces")
  sfExport("MainContrasts")
  sfExport("bestModelSelection_LM")
  sfExport("a.null.dir")

  permLM = sfLapply(1:B,function(b){
    print(b)
    set.seed(b)
    expr.b <- expressionMatrix[,sample(ncol(expressionMatrix))]
    result.b <- replicate(nrow(expr.b),list())
    names(result.b) <- rownames(expr.b)
    result.b<-lapply(1:nrow(expr.b),function(x) bestModelSelection_LM(x,expr.b,covariables3 = Covariables,
                                                                      Main,VariableListOne,VariableListTwo,MainContrasts))
    afile <- paste(a.null.dir, '/result_null_', filename, "_", b,'.rdata',sep='')
    save(result.b,file=afile)
  })
  sfStop()

  #calculate p-value
  #since main effect might have multiple levels, we need to find out which are main effects
  # for debugging usage ---
  # indLM = get(load(paste0(out.dir, "/LM_true", filename, ".rdata")))
  # a.null.dir = paste0(out.dir, "/null_", filename)
  # permLM = parallel::mclapply(1:B, function(b){
  #   get(load(paste(a.null.dir, '/result_null_', filename, "_", b,'.rdata',sep='')))
  # }, mc.cores = 20)
  # for debugging usage ---
  
  All.coefs = rownames(indLM[[1]]$coefficients)
  # intercept.ind = grep("Intercept", All.coefs)
  main.ind = unlist(sapply(colnames(Main), function(a.main){
    grep(a.main, All.coefs)
  }))
  main.res = lapply(main.ind, function(a.ind){
    a.coef = All.coefs[a.ind]
    est.true = do.call(c, lapply(1:length(indLM), function(a){indLM[[a]]$coefficients[a.coef, "Estimate"]}))
    sd.true = do.call(c, lapply(1:length(indLM), function(a){indLM[[a]]$coefficients[a.coef, "Std. Error"]}))
    t.true = do.call(c, lapply(1:length(indLM), function(a){indLM[[a]]$coefficients[a.coef, "t value"]}))
    p.true = do.call(c, lapply(1:length(indLM), function(a){indLM[[a]]$coefficients[a.coef, "Pr(>|t|)" ]}))
    p.null = do.call(c, lapply(1:length(permLM), function(b){
      do.call(c, lapply(1:length(permLM[[b]]), function(a){permLM[[b]][[a]]$coefficients[a.coef, "Pr(>|t|)" ]}))
    }))
    p.pool = c(p.true,p.null)
    p.corrected= rank(p.pool)[1:length(p.true)]/length(p.pool)
    p.BH = p.adjust(p.corrected, "BH")
    a.tab = data.frame(est = est.true, sd = sd.true, t_statistics = t.true, p.observed = p.true, p.corrected = p.corrected, p.BH = p.BH)
    colnames(a.tab)[1] = a.coef
    rownames(a.tab) = rownames(expressionMatrix)
    a.tab$GeneSymbol = GeneSymbol
    write.csv(a.tab, paste0(out.dir, "/", filename, "_", a.coef, "_LM_CovSelect.csv"))
    return(a.tab)
  })

  if(!is.null(MainContrasts)){
    con.names = rownames(indLM[[1]]$res.contrasts)
    con.res = lapply(1:length(con.names), function(a.ind){
      a.coef = con.names[a.ind]
      est.true = do.call(c, lapply(1:length(indLM), function(a){indLM[[a]]$res.contrasts[a.coef, "est"]}))
      sd.true = do.call(c, lapply(1:length(indLM), function(a){indLM[[a]]$res.contrasts[a.coef, "sd"]}))
      t.true = do.call(c, lapply(1:length(indLM), function(a){indLM[[a]]$res.contrasts[a.coef, "t.statistics"]}))
      p.true = do.call(c, lapply(1:length(indLM), function(a){indLM[[a]]$res.contrasts[a.coef, "pvalue"]}))
      p.null = do.call(c, lapply(1:length(permLM), function(b){
        do.call(c, lapply(1:length(permLM[[b]]), function(a){permLM[[b]][[a]]$res.contrasts[a.coef, "pvalue"]}))
      }))
      p.pool = c(p.true,p.null)
      p.corrected= rank(p.pool)[1:length(p.true)]/length(p.pool)
      p.BH = p.adjust(p.corrected, "BH")
      a.tab = data.frame(est = est.true, sd = sd.true, t_statistics = t.true, p.observed = p.true, p.corrected = p.corrected, p.BH = p.BH)
      colnames(a.tab)[1] = a.coef
      rownames(a.tab) = rownames(expressionMatrix)
      a.tab$GeneSymbol = GeneSymbol
      a.tab = a.tab[order(a.tab$p.corrected), ]
      write.csv(a.tab, paste0(out.dir, "/", filename, "_", a.coef, "_LM_CovSelect.csv"))
      return(a.tab)
    })
  }

  #extracts coefficients
  main.tab = do.call(cbind.data.frame, lapply(main.res, function(a){a[, 1, drop = FALSE]}))
  if(!is.null(MainContrasts)){
    con.tab = do.call(cbind.data.frame, lapply(con.res, function(a){a[, 1, drop = FALSE]}))
    main.tab = cbind.data.frame(main.tab, con.tab)
  }

  Covariables.extend = colnames(model.matrix(~., data = Covariables))[-1]
  Cov.count = sapply(colnames(Covariables), function(a){
    sum(grepl(a, Covariables.extend))
  })
  Cov.Factors = names(Cov.count)[Cov.count>1]

  n.main = 1+length(main.ind) # 1 is for intercept
  cov.lists = lapply(colnames(Covariables), function(a.cov.name){
    a.cov.res = do.call(rbind, lapply(1:length(indLM), function(a.ind){
      if(a.cov.name%in%indLM[[a.ind]]$bestBICVariables){
        match.ind = which(indLM[[a.ind]]$bestBICVariables==a.cov.name)
        if(match.ind==1){
          if(a.cov.name%in%Cov.Factors){
            n.match = Cov.count[a.cov.name]
            a.tab = indLM[[a.ind]]$coefficients[(n.main+1):(n.main+n.match), "Estimate"]
            a.tab = t(a.tab)
            colnames(a.tab) = names(unlist(sapply(Covariables.extend, function(a.extend.name){grep(a.cov.name, a.extend.name)})))
            return(a.tab)
          }else{
            return(indLM[[a.ind]]$coefficients[(n.main+1), "Estimate"])
          }
        }else{
          if(indLM[[a.ind]]$bestBICVariables[1]%in%Cov.Factors){
            cov1 = indLM[[a.ind]]$bestBICVariables[1]
            n.match1 = Cov.count[cov1]
            return(indLM[[a.ind]]$coefficients[(n.main+n.match1), "Estimate"])
          }else{
            return(indLM[[a.ind]]$coefficients[(n.main+2), "Estimate"])
          }
        }
      }else{
        a.coef = 0
      }
      return(a.coef)
    }))
  })
  cov.tab = do.call(cbind.data.frame, cov.lists)
  colnames(cov.tab) = Covariables.extend

  coef.tab = cbind.data.frame(main.tab, cov.tab)
  rownames(coef.tab) = rownames(expressionMatrix)
  coef.tab$GeneSymbol = GeneSymbol
  write.csv(coef.tab, paste0(out.dir, "/", filename, "_", "Coefficients", "_LM_CovSelect.csv"))
}

Individual_Coeff_adj = function(index,data,diagnosis,coefficient){
  expr = as.numeric(data[index,])
  result = rep(0,ncol(coefficient))
  names(result) = colnames(coefficient)
  for (i in 1:ncol(coefficient)) {
    variable = unlist(coefficient[,i])
    fit = lm(expr~diagnosis+variable)
    result[i] = summary(fit)$coefficients[3,4]
  }
  result
}


bestModelSelection_LM <- function(aGeneIndex,expressionMatrix,covariables3,Main,VariableListOne,VariableListTwo, MainContrasts){
  curGeneExpression = as.numeric(expressionMatrix[aGeneIndex,])
  Main2 = cbind.data.frame(curGeneExpression = curGeneExpression, Main)

  ## linear model with adjusting for covariates
  basicFormula = update(curGeneExpression ~ .,Main2,.~.)
  basicFormula1 = update(basicFormula,.~.+unlist(covariables3[,aVariable]))
  basicFormula2 = update(basicFormula,.~.+unlist(covariables3[,aVariable1]) + unlist(covariables3[,aVariable2]))

  bestNullFormula = update(curGeneExpression~1, .~.)

  fm0 = lm(basicFormula, data = Main2)
  bestfm = fm0
  bestBIC = BIC(fm0)
  bestBICVariables = NULL

  ## adjusting for one covariate
  for(i in 1:length(VariableListOne)){
    aVariable = VariableListOne[i]
    fm1 = lm(basicFormula1, data = Main2)
    fm1BIC = tryCatch(BIC(fm1), error=function(e) e)

    if(!sum(class(fm1BIC)=="numeric")){
      cat(i)
      cat("th, one variable")
      cat("\n")
      cat("gene index: ")
      cat(aGeneIndex)
      cat("\n")
      next
    }

    if(fm1BIC<bestBIC){
      bestBIC = fm1BIC
      bestBICVariables = aVariable
      bestfm = fm1
    }
  }

  ## mixed model with adjusting for two covariates
  for(i in 1:dim(VariableListTwo)[2]){
    aVariable1 = VariableListTwo[1,i]
    aVariable2 = VariableListTwo[2,i]
    fm2 = lm(basicFormula2, data = Main2)

    fm2BIC = tryCatch(BIC(fm2), error=function(e) e)
    if(!sum(class(fm2BIC)=="numeric")){
      cat(i)
      cat("th, two variable")
      cat("\n")
      cat("gene index: ")
      cat(aGeneIndex)
      cat("\n")
      next
    }
    if(fm2BIC<bestBIC){
      bestBIC = fm2BIC
      bestBICVariables = VariableListTwo[,i]
      bestfm = fm2
    }
  }

  if(length(bestBICVariables)==0){
    bestNullFormula = update(bestNullFormula,.~.)
  } else if(length(bestBICVariables)==1){
    bestNullFormula = update(bestNullFormula,.~.+unlist(covariables3[,bestBICVariables]))
  } else if(length(bestBICVariables)==2){
    bestNullFormula = update(bestNullFormula,.~.+unlist(covariables3[,bestBICVariables[1]])+unlist(covariables3[,bestBICVariables[2]]))
  } else {
    stop('length of best BIC variables error')
  }



  bestNullfm = lm(bestNullFormula, data = Main2)
  lrt.pvalue = anova(bestfm,bestNullfm)[2,'Pr(>F)']

  coefficients = summary(bestfm)$coefficients

  if(!is.null(MainContrasts)){
    #add sum of interaction terms into the full model
    vcov.main = vcov(bestfm)
    res.contrasts = lapply(MainContrasts, function(a){
      eff = coefficients[a$var0, "Estimate"]
      for(i in 1:length(a$fun)){
        eff = a$fun[[i]](eff, coefficients[a$var1[i], "Estimate"])
      }
      eff.sd = sqrt(sum(vcov.main[c(a$var0, a$var1), c(a$var0, a$var1)]))
      eff.t = eff/eff.sd
      eff.p = pt(abs(eff.t), df = bestfm$df.residual, lower.tail = FALSE)*2
      one.row = data.frame(est = eff, sd = eff.sd, t.statistics = eff.t, pvalue = eff.p)
      return(one.row)
    })
    res.contrasts.tab = do.call(rbind.data.frame, res.contrasts)
    rownames(res.contrasts.tab) = names(MainContrasts)
  }else{
    res.contrasts.tab = NULL
  }

  # if(is.null(bestBICVariables)){
  #  residuals = curGeneExpression
  #} else if(length(bestBICVariables)==1){
  # residuals = curGeneExpression - covariables3[,bestBICVariables]*coefficients['covariables3[, aVariable]','Estimate']
  #} else if(length(bestBICVariables)==2){
  #residuals = curGeneExpression - as.matrix(covariables3[,bestBICVariables])%*%as.matrix(coefficients[c('covariables3[, aVariable1]','covariables3[, aVariable2]'),'Estimate'])
  #}
  result = list(coefficients=coefficients,bestBICVariables=bestBICVariables,lrt.pvalue=lrt.pvalue,
                res.contrasts = res.contrasts.tab) #residuals=residuals
  return(result)
}
