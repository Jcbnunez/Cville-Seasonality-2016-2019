# Outputs compiled genomic data from SLiM and PLINK
# Connor Murray 11.16.2021
# ijob -A berglandlab_standard --mem=2G -p standard -c 4
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(forcats)
library(gmodels)
library(poolSeq)
library(poolfstat)
library(FactoMineR)
library(tidyverse)

# Executable in command line
arg <- commandArgs(TRUE)
path.name <- arg[1] # pathway to output
file <- arg[2] # seed for vcf
out.name <- arg[3] # the output file name

path.name="/scratch/csm6hg/slim_bottleneck/tmp/1/503/20193"
file='20193_.vcf'
out.name="/project/berglandlab/connor/BACKUP_scratch/slim_bottleneck/data/freq.output.new"

# VCF output folder
setwd(path.name)

# All frequency data
filenames <- list.files(pattern = paste(file, "$", sep=""))

# Forloop - read files
out <- foreach(i = 1:length(filenames), .combine="rbind", .errorhandling="remove") %do% {
  
  # Read VCF
  vcfi <- fread(filenames[i])
  
  # Read in VCF and extract info
  dt <- data.table(data.table(data.table(CHROM=vcfi[,1],
                   POS=as.numeric(vcfi$POS),
                   ID=as.numeric(tstrsplit(tstrsplit(vcfi$INFO, ";")[[1]], 
                                           "MID=")[[2]]),
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
                              sim.gen %in% c(1:16) ~ "Year_1",
                              sim.gen %in% c(17:33) ~ "Year_2",
                              sim.gen %in% c(34:50) ~ "Year_3")) %>% 
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
                          select(POS, AF, vcf, ALT, REF))

# Get allele frequencies  
a <- data.table(data.table(alefeq_fd %>% 
                             group_by(vcf) %>%
                             summarize(table(POS))) %>% 
                  mutate(POS=as.numeric(table.POS..V1),
                         n=as.numeric(table.POS..N)) %>% 
                  select(vcf, POS, n))

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
                          select(vcf, POS, coverage) %>%
                          pivot_wider(names_from = POS, values_from = coverage))

# Get counts matrix
counts <- data.table(base_obj %>% 
                       select(vcf, POS, count) %>%
                       pivot_wider(names_from = POS, values_from = count))

# SNP info object
snp.info <- data.table(Chromosome=1,
                       Position=unique(alefeq_fd1$POS),
                       RefAllele=unique(alefeq_fd1$REF),
                       AltAllele=unique(alefeq_fd1$ALT))

# Make pool
pool <- new("pooldata",
            npools = dim(counts)[1],
            nsnp = dim(counts)[2]-1,
            refallele.readcount = counts[,-"vcf", with=F]%>%as.matrix%>%t,
            readcoverage = coverages[,-"vcf", with=F]%>%as.matrix%>%t,
            snp.info = snp.info,
            poolsizes = rep(50*2, 50),
            poolnames = as.character(paste("pool",1:50, sep="")))

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
  left_join(compi %>% select(samp1=vcf, sample), by=c("sample1"="sample")) %>% 
  left_join(compi %>% select(samp2=vcf, sample), by=c("sample2"="sample"))

# Change names of columns
colnames(pair.dt)[2:3] <- c("FST", "Q2")

# Join with comparison output
total <- data.table(data.table(pair.dt %>% select(-c(comparison,Q2,sample1,sample2,Nsnp)) %>% 
                          left_join(out %>% 
                             select(vcf, nMax, nMin, replicate, seed, YearS1=Year, sim.gen1=sim.gen) %>% 
                                    distinct(), 
                                           by=c("samp1"="vcf")) %>%
                           left_join(out %>% 
                             select(vcf, YearS2=Year, sim.gen2=sim.gen) %>% 
                                    distinct(), 
                                           by=c("samp2"="vcf"))) %>% 
                      mutate(Year = case_when(YearS1==YearS2 ~ "Within",
                                              !YearS1==YearS2 ~ "Overwinter")))
# Fix FST output
total <- data.table(total %>% mutate(iteration=row_number()) %>% 
                      select(samp1,samp2,FST,iteration, 
                             nMax,nMin,replicate,seed,
                             YearS1,sim.gen1,YearS2,sim.gen2,
                             Year))

# Get AF pool matrix
af <- data.table(base_obj %>% 
                 select(vcf, POS, AF.pool) %>% 
                 pivot_wider(names_from = POS, values_from = AF.pool))

# Run pca - impute average
pca <- PCA(af %>% select(-c(vcf)), 
           graph = FALSE)

# Include sample numbers
pca.dt <- data.table(pca$ind$coord, 
                     Sample=af$vcf, 
                     Sample.gen=tstrsplit(af$vcf, "_")[[5]]) %>% 
  mutate(Year = case_when(Sample.gen %in% c(1:16) ~ "Year1",
                          Sample.gen %in% c(17:33) ~ "Year2",
                          Sample.gen %in% c(34:50) ~ "Year3"))

# Save fst output
write.csv(total, 
          file = paste(out.name, 
                       unique(total$nMax), 
                       unique(total$nMin),
                       unique(total$replicate),
                       unique(total$seed),
                       "csv", 
                       sep="."), 
          quote = F, 
          row.names = F)

# Save pca output
write.csv(pca.dt, 
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
