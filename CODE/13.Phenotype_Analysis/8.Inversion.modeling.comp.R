###model effect of inversion status 2L on DGRP phenotypes
###inversion modeling
library(lme4)
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel) ### AOB
registerDoParallel(4)

#setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/")
#system("cp /project/berglandlab/Yang_Adam/inversion.csv ./")

##load in inversion data
inversions = fread("inversion.csv")

#bring in phenotype data
#system("cp /project/berglandlab/Yang_Adam/wideform.fixed.phenotable.RDS ./")

full.data = readRDS("wideform.fixed.phenotable.RDS")
#daria's issue
x = as.vector(colSums(is.na(full.data)) < 50)
narrow.data =  full.data[, ..x]
#need to melt this table into long form
melt.data = melt(full.data, id.vars = "ral_id", variable.name = "pheno_env", value.name = "avg")

#
#fix up inversions data table- fix column names, make line names universal at RAL_ID
columnnames = names(inversions)

columnnames <- gsub(pattern = "\\)", replacement = "_", x = columnnames)
columnnames <- gsub(pattern = "\\(", replacement = "_", x = columnnames)
colnames(inversions) = columnnames
inversions <- inversions[-c(1,2),]
inversions$`DGRP Line` = gsub("DGRP_","", inversions$`DGRP Line`)
colnames(inversions)[1] = "ral_id"

#combine data
total.data <- merge(melt.data, inversions, by = "ral_id")
head(total.data)
# ###graph effect of in2lt on basal activity
# act.data = total.data[pheno_env == "WalkingActivity_standard_Female"]
# 
# ggplot(data = act.data, aes(x = In_2L_t, y = avg)) + geom_boxplot()

#### Run a loop recording model values for each pheno X inversion possibility
#create reference table of different combinations
ref.table = expand.grid(unique(total.data$pheno_env) , columnnames[-1]) #issue, some combinations have no inversions. 
###for this analysis, we are only concerned with 2L. remove the following line if interested in other inversions
ref.table = ref.table %>% filter(Var2 == "In_2L_t")
###create a foreach loop that runs through each combination using the reftable for reference
out <- foreach(f = c(1:dim(ref.table)[1]) ) %do% {
 #f = 461
  #assign pheno and inversion using ref table
  pheno = as.character(ref.table[f,1])
  inversion = as.character(ref.table[f,2])
  guess = inversion[1]
  #filter data to match pheno of interest
  filtered.data = total.data[pheno_env == pheno]
  #use data to create model
  #use reformulate to write formula
  form <- reformulate(termlabels = paste("0 +", inversion, sep = ""),response= "avg")
  testmodel = lm(form, data = filtered.data)
  
  #remove intercept and beta values
  inversiontype = unlist(testmodel$xlevels)
  beta = unlist(testmodel$coefficients)
  
  #create a null model and use anova to compare
  nullmodel= lm(avg ~ 1, filtered.data)
  anova = anova(testmodel, nullmodel)
  
  #return a dataframe with phenotype, inversion, inversion type, beta, P value, and F statistic data
  o <- data.frame( 
    phenotype = pheno,
    inversion = inversion,
    inversion.type = inversiontype,
    beta = beta,
    P.value = unlist(anova[2,6]),
    F.stat = unlist(anova[2,5])
  )
}


short_out <- rbindlist(out)
saveRDS(short_out, "fullpheno.inversions.summary")
short_out = readRDS("fullpheno.inversions.summary")

short_out$p.adjust = p.adjust(short_out$P.value, "fdr")
short_out$permutation.number = 0
short_out$perm.status = "observed"


####repeat analysis with permutations###
# use a forloop and shuffling inversion status 
#shuffle lineids in inversion file, then merge onto phenotypes


perm.out = foreach(i = c(1:10)) %do% {
 # i = 1
  print(i)
  set.seed(i)

  #shuffle inversion status
  inversions$ral_id <- inversions$ral_id[sample(nrow(inversions))]
  
 
  #combine and melt data
  
  melt.data = melt(full.data, id.vars = "ral_id", variable.name = "pheno_env", value.name = "avg")
  total.data <- merge(melt.data, inversions, by = "ral_id")
  #redo modeling # with only ln2lt
  ref.table = expand.grid(unique(total.data$pheno_env) , columnnames[-1])
  ref.table = ref.table %>% filter(Var2 == "In_2L_t")
  ###create a foreach loop that runs through each combination using the reftable for reference
  out <- foreach(f = c(1:dim(ref.table)[1]), .packages = "data.table") %dopar% {
    # f = 1
    #assign pheno and inversion using ref table
    pheno = as.character(ref.table[f,1])
    inversion = as.character(ref.table[f,2])
    guess = inversion[1]
    #filter data to match pheno of interest
    total.data = as.data.table(total.data)
    filtered.data = total.data[pheno_env == pheno]
    #use data to create model
    #use reformulate to write formula
    form <- reformulate(termlabels = paste("0 +", inversion, sep = ""),response= "avg")
    testmodel = lm(form, data = filtered.data)
    
    #remove intercept and beta values
    inversiontype = unlist(testmodel$xlevels)
    beta = unlist(testmodel$coefficients)
    
    #create a null model and use anova to compare
    nullmodel= lm(avg ~ 1, filtered.data)
    anova = anova(testmodel, nullmodel)
    
    #return a dataframe with phenotype, inversion, inversion type, beta, P value, and F statistic data
    o <- data.table( 
      phenotype = pheno,
      inversion = inversion,
      inversion.type = inversiontype,
      beta = beta,
      P.value = unlist(anova[2,6]),
      F.stat = unlist(anova[2,5]),
      permutation.number = i
    )
    o
  }
  permutationdata = rbindlist(out)
  
  # saveRDS(permutationdata,paste0("permutated.model.stats.", i))
  return(permutationdata)
}

fullout = rbindlist(perm.out)
# saveRDS(fullout, "100perms.inversion.allchrom.modeling")
# fullout = readRDS("1000perms.inversion.modeling")

fullout$perm.status = ifelse(fullout$permutation.number == 0, "observed", "permutation")
permutationdata$perm.status = ifelse(permutationdata$permutation.number == 0, "observed", "permutation")

fullout.o = rbind(fullout,permutationdata)


#bind to observed
#fullout.o = fullout.o[,-7]
fullout.o$perm.status = "permutation"
fullout.o =fullout.o  %>% 
  group_by(permutation.number) %>% 
  mutate(p.adjust = p.adjust(P.value, method = "fdr")) %>% 
  as.data.table(.)

#rbind
comp.data = rbind(fullout.o, short_out)
saveRDS(comp.data, "complete.inversion.data")
comp.data = readRDS("complete.inversion.data")

comp.data$chr = tstrsplit(comp.data$inversion, "_")[2]
#lets make data simpler by removing inversion type data we're not using
shortout = comp.data[,-c(3,4)]
shortout = unique(shortout)
#add in groups
#now load in phenotypes with groups
#system("cp /project/berglandlab/Yang_Adam/phenogroups3.22.csv ./")

groups = read_csv("phenogroups3.22.csv")
groups = groups[,c(1,5)]
#remove special character
groups$phenotype = gsub("μ","",groups$phenotype)
shortout$phenotype = gsub("Î¼","",shortout$phenotype)

mergedata = merge(shortout, groups, by = "phenotype")
unique(mergedata$phenotype)

#find number of phenotype models for each classification with pvalue < 0.05
mergedata %>%
  filter(P.value < 0.05) %>%
  group_by(chr,perm.status, permutation.number, group.general) %>%
  summarise(N=n()) %>% 
  as.data.table(.) -> out_count

