#DE analysis with covariate selection
#packages

#source("/Users/xiangningxue/OneDrive - University of Pittsburgh/Research/Analysis/DE_Analysis/DE_CovSelect_LM.R")
source("/home/xix66/code/DE_CovSelect_LM.R")
hdir = "/home/xix66/Seney/Collaboration_HumanPyramidalNeuron"
setwd(hdir)
out.dir = paste0(hdir, "/DE_CovSelect/LM_rm")
dir.create(out.dir, recursive = TRUE)
#load data
dat = read.csv(paste0("data/HumanPyramidalNeuron_data_combat.csv"), row.names = 1)
meta = read.csv(paste0("data/HumanPyramidalNeuron_meta.csv"), row.names = 1, stringsAsFactors = FALSE)
gene_annotation = read.csv(paste0("data/HumanPyramidalNeuron_Gene_annotations.csv"), row.names = 1, stringsAsFactors = FALSE)
all(meta$Sample == colnames(dat))

meta = meta[meta$Race!="As", ]
meta$Race = factor(meta$Race, levels = c("W", "B"))
dat = dat[, meta$Sample]
rownames(meta) = meta$Sample

Regions = as.character(unique(meta$Region))
 for(a in 1:length(Regions)){
#for(a in 3){
  a.Region = Regions[a]
  a.meta = meta[meta$Region==a.Region, ]
  a.dat = dat[, a.meta$Sample]

  Dx = a.meta$Sex
  coefficient = a.meta[, c("Age", "Race", "PMI", "pH", "RIN")]

  # #individual covariate check
  # Res = lapply(1:nrow(a.dat),function(x) Individual_Coeff_adj(index = x,data = a.dat,diagnosis = Dx, coefficient = coefficient))
  # Pvalue_matrix = as.data.frame(do.call(rbind,Res))

  VariableListOne = combn(colnames(coefficient),1)
  VariableListTwo = combn(colnames(coefficient),2)

  library(snowfall) #parallel computing to speed up process
  set.seed(15213)
  sfInit(parallel=TRUE,type="SOCK", cpus=24)
  GeneIndeces = 1:nrow(a.dat)
  sfExportAll()

  indLM<-c()
  sfExport("bestModelSelection_LM")
  indLM<-sfLapply(GeneIndeces,function(x) try(bestModelSelection_LM(x,a.dat,covariables3 = coefficient,
                                                                    Diagnosis = Dx,VariableListOne,VariableListTwo),silent=TRUE ) )
  # indLM.vars = lapply(indLM, `[[`, 2)
  # indLM.race = sapply(indLM.vars, function(a){"Race"%in%a})
  # which(indLM.race)[1]
  # indLM[[20]]

  sfStop()
  save(indLM,file=paste0(out.dir, "/LM_true", a.Region, ".rdata"))

  a.null.dir = paste0(out.dir, "/null_", a.Region)
  dir.create(a.null.dir, recursive = TRUE)

  #setwd("/net/wong02/home/wez97/collaboration/Logan/Human_OUD/output/DE_0905/quantileFilter/null_NAC")
  expr = a.dat
  B<-100
  sfInit(parallel=TRUE,type="SOCK", cpus=24)
  #Corrected
  sfExportAll()
  result.B = sfLapply(1:B,function(b){
    print(b)
    set.seed(b)
    expr.b <- expr[,sample(ncol(expr))]
    result.b <- replicate(nrow(expr.b),list())
    names(result.b) <- rownames(expr.b)
    result.b<-lapply(1:nrow(expr.b),function(x) bestModelSelection_LM(x,expr.b,covariables3 = coefficient,Diagnosis = Dx,
                                                                        VariableListOne,VariableListTwo))
    afile <- paste(a.null.dir, '/result_null_', a.Region, "_", b,'.rdata',sep='')
    save(result.b,file=afile)
  })
}

for(a in 1:length(Regions)){
#for(a in 3){
  a.Region = Regions[a]
  a.meta = meta[meta$Region==a.Region, ]
  a.dat = dat[, a.meta$Sample]

  Dx = a.meta$Sex
  coefficient = a.meta[, c("Age", "Race", "PMI", "pH", "RIN")]
  COVs.covariates = c("Age", "Race", "PMI", "pH", "RIN")

  VariableListOne = combn(colnames(coefficient),1)
  VariableListTwo = combn(colnames(coefficient),2)

  ##P-Value Corrections #Run permutation test to correct for biased pval from LRT (see paper)
  ##Full
  a.null.dir = paste0(out.dir, "/null_", a.Region)
  LM.NULL<-c()
  for(i in 1:100){
    print(i)
    aNull<-get(load(paste(a.null.dir, "/result_null_", a.Region, "_", i,".rdata", sep = "")))
    LM.NULL[[i]]<-aNull
  }

  unlistLM_null<-unlist(LM.NULL)
  names.pval<-grep("lrt.pvalue", names(unlistLM_null))
  pval.null<-unname(as.numeric(unlistLM_null[names.pval]))

  indLM<-get(load(paste0(out.dir, "/LM_true", a.Region, ".rdata")))
  unlistLM_true<-unlist(indLM)
  names.pval<-grep("lrt.pvalue", names(unlistLM_true))
  pval.true<-unname(as.numeric(unlistLM_true[names.pval]))

  es.Dx = c()
  st.Dx = c()
  for( i in 1:length(indLM)){
    es.Dx[i]<-indLM[[i]]$coefficients[2,1]
    st.Dx[i]<-indLM[[i]]$coefficients[2,2]
  }

  p.pool = c(pval.true,pval.null)
  p.corr = rank(p.pool)[1:length(pval.true)]/length(p.pool)
  a.Region_res<-cbind(es.Dx,st.Dx, pval.true, p.corr)
  genes<-row.names(a.dat)
  row.names(a.Region_res)<-make.names(genes, unique = TRUE)

  colnames(a.Region_res)<-c("coefficient","sd","obs.p","corrected.p")
  a.Region_res<-data.frame(a.Region_res)
  a.Region_res$p.BH<-p.adjust(a.Region_res$corrected.p, "BH")
  a.Region_res$Gene_Symbols = gene_annotation$Gene[match(rownames(a.Region_res), gene_annotation$ENSEMBL)]
  a.Region_res_Sorted<-a.Region_res[order(a.Region_res$corrected.p, decreasing = FALSE),]
  head(a.Region_res_Sorted)
  write.csv(a.Region_res_Sorted, paste0(out.dir, "/Male_vs_Female_", "LMcovariateSelection_",  a.Region, "_rm.csv"))


  #Pull out coefficient of adjusted variables for all genes
  #"sex","age","race","PMI","avgpH","RPF_RIN"
  n.main = 1+1
  n.cov = ncol(coefficient)
  n.total = n.main+n.cov
  est.vecs = replicate(n.total, vector())
  names(est.vecs) = c("Intercept", "GenderM", "Age", "Race", "PMI",  "pH", "RIN")
  new.colnames = c("Intercept", "GenderM", "Age", "RaceB", "PMI",  "pH", "RIN")

  for(i in 1:length(indLM)){
    indLM.names = rownames(indLM[[i]]$coefficients)
    indLM.names[2] = "GenderM"

    if(is.null(indLM[[i]]$bestBICVariables)){
      for(m in 1:n.total){
        one.name = names(est.vecs)[m]
        if(m<=n.main){
          est.vecs[[m]][i] = indLM[[i]]$coefficients[grepl(one.name, indLM.names), "Estimate"]
        }else{
          est.vecs[[m]][i] = 0
        }
      }
    }else{
      for(m in 1:n.total){
        one.name = names(est.vecs)[m]
        if(m<=n.main){
          est.vecs[[m]][i] = indLM[[i]]$coefficients[grepl(one.name, indLM.names), "Estimate"]
        }else{
          if(one.name %in% indLM[[i]]$bestBICVariables){
            est.vecs[[m]][i] = indLM[[i]]$coefficients[which(indLM[[i]]$bestBICVariables==one.name)+n.main, "Estimate"]
          }else{
            est.vecs[[m]][i] = 0
          }
        }
      }

    }
  }

  coeff.df = do.call(cbind.data.frame, est.vecs)
  colnames(coeff.df) = new.colnames
  coeff.df$Gene_Symbols = gene_annotation$Gene[match(rownames(a.dat), gene_annotation$ENSEMBL)]
  colnames(coeff.df)= names(est.vecs)
  rownames(coeff.df) = rownames(a.dat)
  write.csv(coeff.df, paste0(out.dir, "/Male_vs_Female__CovariateSelected_coeff_", a.Region, "_rm.csv"))

}



