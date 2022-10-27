# Outputs compiled genomic data from SLiM VCF
# By: Connor Murray on 6.16.2022
# ijob -A berglandlab_standard --mem=2G -p standard -c 4
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(forcats)
library(gmodels)
library(gtools)
library(poolSeq)
library(poolfstat)
library(FactoMineR)
library(tidyverse)
library(adegenet)

# Executable in command line
arg <- commandArgs(TRUE)
path.name <- arg[1] # pathway to output
file <- arg[2] # seed for vcf
out.name <- arg[3] # the output file name

# 25k to 25k 0% bottleneck
#path.name="/scratch/csm6hg/bottleneck/tmp/19/100/24223"
#file="24223_.vcf"
#out.name="/project/berglandlab/connor/bottleneck/data/freq.output.new"

# 25k to 2.5k 90%
#path.name="/scratch/csm6hg/bottleneck/tmp/3/67/4941"
#file="4941_.vcf"
#out.name="/project/berglandlab/connor/bottleneck/data/freq.output.new"

# 25k to 250 99%
#path.name="/scratch/csm6hg/bottleneck/tmp/4/25/14399"
#file="14399_.vcf"
#out.name="/project/berglandlab/connor/bottleneck/data/freq.output.new"

# VCF output folder
setwd(path.name)

# All frequency data
filenames <- list.files(pattern = paste(file, "$", sep=""))

## define Euclidean distances
euc.dist.3d <- function(coord_vec) {
  #x1, y1, z1, x2, y2, z2
  x1 = coord_vec[1] 
  y1 = coord_vec[2] 
  z1 = coord_vec[3] 
  x2 = coord_vec[4]  
  y2 = coord_vec[5]
  z2 = coord_vec[6]
  sqrt( ((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))
} 

## define Euclidean distances
euc.dist.2d <- function(coord_vec) {
  #x1, y1, z1, x2, y2, z2
  x1 = coord_vec[1] 
  y1 = coord_vec[2] 
  x2 = coord_vec[3]  
  y2 = coord_vec[4]
  sqrt( ((x1 - x2)^2) + ((y1 - y2)^2) )
} 

# Forloop - read files
out <- foreach(i = 1:length(filenames), .combine="rbind", .errorhandling="remove") %do% {
  
  # Read VCF
  vcfi <- fread(filenames[i])
  
  # Read in VCF and extract info
  dt <- data.table(data.table(data.table(CHROM=vcfi[,1],
                   POS=as.numeric(vcfi$POS),
                   ID=as.numeric(tstrsplit(tstrsplit(vcfi$INFO, ";")[[1]], 
                                           "MID=")[[2]]),
                   REF=vcfi$REF,
                   ALT=vcfi$ALT,
                   ALT_CTS=as.numeric(tstrsplit(tstrsplit(vcfi$INFO, ";")[[7]], 
                                                "AC=")[[2]]),
                   OBS_CTS=as.numeric(tstrsplit(filenames[i], "_")[[4]])*2) %>% 
            mutate(REF_CTS=(OBS_CTS-ALT_CTS),
                   AF=(ALT_CTS/(OBS_CTS))) %>% 
            mutate(MAF = case_when(
                   AF >= 0.5 ~ 1-AF,
                   AF < 0.5 ~ AF)), 
                   vcf = filenames[i],
                   slurm_ID = tstrsplit(filenames[i], "_")[[3]],
                   nSamp = tstrsplit(filenames[i], "_")[[4]],
                   sim.gen = tstrsplit(filenames[i], "_")[[5]],
                   fin.gen = tstrsplit(filenames[i], "_")[[6]],
                   nMax = tstrsplit(filenames[i], "_")[[7]],
                   nMin = tstrsplit(filenames[i], "_")[[8]],
                   replicate = tstrsplit(filenames[i], "_")[[9]],
                   seed = tstrsplit(filenames[i], "_")[[10]],
                   iteration = i) %>% 
                mutate(Year = case_when(
                              sim.gen %in% c(1:16) ~ 1,
                              sim.gen %in% c(17:33) ~ 2,
                              sim.gen %in% c(34:50) ~ 3)) %>% 
                mutate(Year.samp = paste(Year, nMax, nMin, replicate, sep="_")))
  
  # Fix column name
  colnames(dt)[1] <- "CHROM"
 
  # Progress message
  print(paste("SLURM_ID-", 
              tstrsplit(filenames[i], "_")[[3]], 
              "-Iteration-", i, ": ", 
              round((i/length(filenames))*100, digits=2), 
              "%", " complete", sep=""))

  # Finish
  return(dt)
}

# Standardize number of SNPs - choose 3000 random SNPs
pos.samp <- sample(unique(out$POS), 3000)
out <- out[POS %in% pos.samp]

# VCF list to make pools
vcf.list <- unique(out$vcf)

# Pairwise comparisons across and within years
compare <- data.frame(expand.grid(vcf.list, vcf.list))
colnames(compare) <- c("Sample1", "Sample2")

# Remove pairwise duplicates
compare <- data.table(compare[!duplicated(t(apply(compare, 1, sort))),])
rownames(compare) <- 1:dim(compare)[1]
compare <- data.table(compare[!Sample1==Sample2])

# Generate outfile object
outfile = data.frame(
  samp1 = rep(NA, dim(compare)[1]),
  samp2 = rep(NA, dim(compare)[1]),
  FST = rep(NA, dim(compare)[1]))

# Extract allele freqeuncy info
alefeq_fd <- data.table(out %>% 
                          group_by(vcf) %>% 
                          dplyr::select(POS, AF, vcf, ALT, REF))

# Get allele frequencies  
a <- data.table(data.table(alefeq_fd %>% 
                             group_by(vcf) %>%
                             summarize(table(POS))) %>% 
                  mutate(POS=as.numeric(table.POS..V1),
                         n=as.numeric(table.POS..N)) %>% 
                  dplyr::select(vcf, POS, n))

# Merge and Create pool
alefeq_fd1 <- data.table(a[n == 1] %>% 
               left_join(alefeq_fd, by=c("vcf", "POS")) %>%
                 group_by(vcf) %>% 
                 mutate(sample.alleles(AF, 
                                       size=60, 
                                       mode="coverage")))

# Fix column names
colnames(alefeq_fd1)[7:8] <- c("AF.pool", "coverage")

# Rbind pools
base_obj <- data.table(alefeq_fd1 %>% 
                         mutate(count = as.numeric(AF.pool)*as.numeric(coverage),
                                SNP_id = paste("sim", POS, sep = "_"))) 

# Long to wide coverages matrix
coverages <- data.table(base_obj %>% 
                          dplyr::select(vcf, POS, coverage) %>%
                          pivot_wider(names_from = POS, values_from = coverage))

# Replace NAs with coverage of 60
#coverages[is.na(coverages)]=60

# Get counts matrix
counts <- data.table(base_obj %>% 
                       dplyr::select(vcf, POS, count) %>%
                       pivot_wider(names_from = POS, values_from = count))

#counts[is.na(counts)]=0

# SNP info object
snp.info <- data.table(Chromosome=1,
                       Position=unique(alefeq_fd1$POS),
                       RefAllele=unique(alefeq_fd1$REF),
                       AltAllele=unique(alefeq_fd1$ALT))

#counts[is.na(counts)]=0

# Make pool
pool <- new("pooldata",
            npools = dim(counts)[1],
            nsnp = dim(counts)[2]-1,
            refallele.readcount = counts[,-"vcf", with=F]%>%as.matrix%>%t,
            readcoverage = coverages[,-"vcf", with=F]%>%as.matrix%>%t,
            snp.info = snp.info,
            poolsizes = rep(50*2, 50),
            poolnames = as.character(paste("pool",1:50, sep="")))

# Allele frequency spectrum
af.dt <- unlist(data.table(pool@refallele.readcount/100) %>% 
       mutate(across(V1:V50, funs(case_when(. > 0.5 ~ 1-., 
                                            . <= 0.5 ~ .)))) %>%
       select(V1_case_when:V50_case_when))

# Pairwise FST Function 
pair <- compute.pairwiseFST(pool, method = "Anova", verbose = TRUE)

# Compress FST
compi <- data.table(vcf=unique(out$vcf), 
                    gen=unique(out$sim.gen)) %>% 
                    mutate(sample=paste("pool", gen, sep=""))

# Join VCF with pairwise FST
pair.dt <- data.table(comparison=rownames(pair@values),
                      pair@values) %>% 
  mutate(sample1=tstrsplit(comparison, ";")[[1]],
         sample2=tstrsplit(comparison, ";")[[2]]) %>% 
  left_join(compi %>% dplyr::select(samp1=vcf, sample), by=c("sample1"="sample")) %>% 
  left_join(compi %>% dplyr::select(samp2=vcf, sample), by=c("sample2"="sample"))

# Change names of columns
colnames(pair.dt)[2:3] <- c("FST", "Q2")

# Join with comparison output
total <- data.table(data.table(pair.dt %>% 
                             dplyr::select(-c(comparison,Q2,sample1,sample2,Nsnp)) %>% 
                          left_join(out %>% 
                             dplyr::select(vcf, nMax, nMin, replicate, seed, YearS1=Year, sim.gen1=sim.gen) %>% 
                                    distinct(), 
                                           by=c("samp1"="vcf")) %>%
                           left_join(out %>% 
                             dplyr::select(vcf, YearS2=Year, sim.gen2=sim.gen) %>% 
                                    distinct(), 
                                           by=c("samp2"="vcf"))) %>% 
                      mutate(Year.diff=abs(as.numeric(YearS1)-as.numeric(YearS2))) %>% 
                      mutate(Year=case_when(Year.diff==0 ~ "Within",
                                            Year.diff==1 ~ "Overwinter",
                                            Year.diff==2 ~ "Multi")))

# Fix FST output
total <- data.table(total %>% 
                      mutate(iteration=row_number()) %>% 
                      dplyr::select(samp1,samp2,FST,iteration, 
                             nMax,nMin,replicate,seed,
                             YearS1,sim.gen1,YearS2,sim.gen2,
                             Year))

#total %>% group_by(Year) %>% summarize(mean.fst=mean(FST))

# Get AF pool matrix
af <- data.table(base_obj %>% 
                 dplyr::select(vcf, POS, AF.pool) %>% 
                 pivot_wider(names_from = POS, values_from = AF.pool))

### Run pca - impute average ###
pca <- PCA(af %>% dplyr::select(-c(vcf)), 
           graph = F)

# Include sample numbers
pca.dt <- data.table(pca$ind$coord, 
                     pc_perc=t(pca$eig[c(1:5),2]),
                     Sample=af$vcf, 
                     Sample.gen=tstrsplit(af$vcf, "_")[[5]]) %>% 
  mutate(Year = case_when(Sample.gen %in% c(1:16) ~ "Year1",
                          Sample.gen %in% c(17:33) ~ "Year2",
                          Sample.gen %in% c(34:50) ~ "Year3")) %>% 
  mutate(samp.pop=paste(Sample.gen, Year, sep="_"))
colnames(pca.dt)[6:10] <- str_remove_all(colnames(pca.dt)[6:10], " ")

#pca.dt %>% ggplot(.,aes(x=Dim.1,y=Dim.2, color=Year)) +geom_point(size=3)

### Euclidean distances ###
permutations(n = length(pca.dt$samp.pop), r = 2, repeats.allowed = F, v = pca.dt$samp.pop) %>%
  as.data.frame() -> perm_samps

perm_samps %>%
  dplyr::select(samp.pop = V1) %>%
  left_join(pca.dt[,c("samp.pop","Dim.1", "Dim.2", "Dim.3")]) %>%
  dplyr::select(sampleId.s1 = samp.pop,
                Dim.1.s1 = Dim.1,
                Dim.2.s1 = Dim.2,
                Dim.3.s1 = Dim.3) -> S1_dat

perm_samps %>%
  dplyr::select(samp.pop = V2) %>%
  left_join(pca.dt[,c("samp.pop","Dim.1", "Dim.2", "Dim.3")]) %>%
  dplyr::select(sampleId.s2 = samp.pop,
                Dim.1.s2 = Dim.1,
                Dim.2.s2 = Dim.2,
                Dim.3.s2 = Dim.3) -> S2_dat

cbind(S1_dat, S2_dat) -> euc.in.dat

eu_dist.2d.D.pca = foreach(i = 1:nrow(euc.in.dat), .combine = c )%do%{
  in_coords = unlist(c(euc.in.dat[i, c("Dim.1.s1",  "Dim.2.s1")],
                       euc.in.dat[i, c("Dim.1.s2",  "Dim.2.s2")]))
  euc.dist.2d(as.vector(in_coords))
}

# Summarize data - 2D
cbind(euc.in.dat, eu_dist.2d.D.pca) %>%
  as.data.table %>% 
  separate(sampleId.s1, into = c("gen1", "year1"), sep = "_Year", remove = F) %>%
  separate(sampleId.s2, into = c("gen2", "year2"), sep = "_Year", remove = F) %>% 
  mutate(year_diff = abs(as.numeric(year1) - as.numeric(year2))) -> Euclidean.analysis.results.2d.pca

euc.dist2 <- data.table(Euclidean.analysis.results.2d.pca %>%
  group_by(year.diff=as.factor(year_diff)) %>%
  summarize(Mean.euc = mean(eu_dist.2d.D.pca, na.rm = T),
            Median.euc = quantile(eu_dist.2d.D.pca, 0.5),
            IQR05.euc = quantile(eu_dist.2d.D.pca, 0.05),
            IQR95.euc = quantile(eu_dist.2d.D.pca, 0.95)),
  data="Euc2D")

# Rbind euclidean distance - final
fin.euc <- data.table(euc.dist2,
                      replicate=unique(total$replicate),
                      nMax=unique(total$nMax),
                      nMin=unique(total$nMin),
                      model=paste(unique(total$nMax), 
                                  unique(total$nMin), sep="_"))

#fin.euc %>% ggplot(., aes(x=year.diff, y=Mean.euc, ymin=IQR05.euc, ymax=IQR95.euc)) + facet_wrap(~data) + geom_pointrange() +theme_bw()

# Allele frequency partitioning
af.dt <- data.frame(af %>% 
  mutate(Sample.gen=tstrsplit(vcf, "_")[[5]]) %>% 
  mutate(Year = case_when(Sample.gen %in% c(1:16) ~ "Year1",
                          Sample.gen %in% c(17:33) ~ "Year2",
                          Sample.gen %in% c(34:50) ~ "Year3")) %>% 
  mutate(samp.pop=paste(Sample.gen, Year, sep="_")) %>%
  select(-c(vcf, Sample.gen, samp.pop, Year)))

af.meta <- data.frame(af %>% 
      mutate(Sample.gen=tstrsplit(vcf, "_")[[5]]) %>% 
      mutate(Year = case_when(Sample.gen %in% c(1:16) ~ "Year1",
                              Sample.gen %in% c(17:33) ~ "Year2",
                              Sample.gen %in% c(34:50) ~ "Year3")) %>% 
      mutate(samp.pop=paste(Sample.gen, Year, sep="_")) %>%
      select(c(Year, samp.pop, vcf)))

### AF variance ###
af.var <- data.table(af.dt %>% 
              cbind(af.meta) %>% 
              summarize(across(contains("X"), ~var(.x, na.rm = T))) %>% 
              summarize(Mean.var=mean(c_across(cols=everything()), na.rm = T),
                        Median.var=quantile(c_across(cols=everything()), probs=0.5, na.rm = T),
                        IQR95.var=quantile(c_across(cols=everything()), probs=0.95, na.rm = T),
                        IQR05.var=quantile(c_across(cols=everything()), probs=0.05, na.rm = T)),
              replicate=unique(total$replicate),
              nMax=unique(total$nMax),
              nMin=unique(total$nMin),
              model=paste(unique(total$nMax), 
                          unique(total$nMin), sep="_"))

# Impute average
for(i in 1:ncol(af.dt)) {
  af.dt[,i][is.na(af.dt[,i])] <- mean(af.dt[,i], na.rm = TRUE)
}

# Perform DAPC
af.dt %>%  
  dapc(., grp=as.factor(af.meta$Year), n.pca=9, n.da=5) ->
  dapc_optim

dapc_optim$ind.coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  cbind(af.meta) ->
  dapc_optim_dat

#scatter(dapc_optim)

### Output information ###

# Save fst output
write.csv(total, 
          file = paste(out.name, 
                       "fst",
                       unique(total$nMax), 
                       unique(total$nMin),
                       unique(total$replicate),
                       unique(total$seed),
                       "csv", 
                       sep="."), 
          quote = F, 
          row.names = F)

# Combine PCA and DAPC information
pca.dt2 <- data.table(pca.dt %>% 
             left_join(dapc_optim_dat %>% 
                         select(-c('sampleId')), 
                       by=c('samp.pop', "Year", "Sample"="vcf")))

# Save PCA output
write.csv(pca.dt2, 
          file = paste(out.name,
                       "pca",
                       unique(total$nMax), 
                       unique(total$nMin),
                       unique(total$replicate),
                       unique(total$seed),
                       "csv", 
                       sep="."), 
          quote = F, 
          row.names = F)

# Write Euclidean distance
write.csv(fin.euc, 
          file = paste(out.name,
                       "euc",
                       unique(total$nMax), 
                       unique(total$nMin),
                       unique(total$replicate),
                       unique(total$seed),
                       "csv", 
                       sep="."), 
          quote = F, 
          row.names = F)

# Write AF Variance
write.csv(af.var, 
          file = paste(out.name,
                       "afvar",
                       unique(total$nMax), 
                       unique(total$nMin),
                       unique(total$replicate),
                       unique(total$seed),
                       "csv", 
                       sep="."), 
          quote = F, 
          row.names = F)
