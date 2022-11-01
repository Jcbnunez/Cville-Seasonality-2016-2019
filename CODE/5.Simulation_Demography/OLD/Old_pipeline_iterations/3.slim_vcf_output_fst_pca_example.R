# Outputs compiled genomic data from SLiM and PLINK
# Connor Murray 11.15.2021
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
library(doParallel)
library(tidyverse)

# Executable in command line
arg <- commandArgs(TRUE)
#path.name <- arg[1] # pathway to output
#file <- arg[2] # seed for vcf
#out.name <- arg[3] # the output file name

path.name="/project/berglandlab/connor/Shared_w_Alan"
file='31638_.acount'
out.name="/project/berglandlab/connor/Shared_W_Alan/freq.output"

# VCF output folder
setwd(path.name)

# All frequency data
filenames <- list.files(pattern = paste(file, "$", sep=""))

# Register cores
registerDoParallel(cores = 4)

# Forloop - read files
out <- foreach(i = 1:length(filenames), .combine="rbind", .errorhandling="remove") %do% {
  
  # Load allele frequency information
  dt <- data.table(data.table(fread(filenames[i]) %>% 
                   mutate(REF_CTS=(OBS_CT-ALT_CTS),
                          AF=(ALT_CTS/(OBS_CT))) %>% 
                   mutate(MAF = case_when(
                          AF >= 0.5 ~ 1-AF,
                          AF < 0.5 ~ AF
                   )), 
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
                                SNP_id = as.character(paste("sim", POS, sep = "_")))) 

# Long to wide coverages matrix
coverages <- as.matrix(base_obj %>% 
                          select(vcf, SNP_id, coverage) %>%
                          pivot_wider(names_from = SNP_id, 
                                      values_from = coverage, 
                                      values_fill = 0))

coverages <- apply(coverages[,-1], 1, as.numeric)
colnames(coverages) <- as.character(paste("pool",1:50, sep=""))
rownames(coverages) <- unique(paste("SNP", alefeq_fd1$POS, sep="_"))

# Get counts matrix
counts <- as.matrix(base_obj %>% 
                       select(vcf, SNP_id, count) %>%
                       pivot_wider(names_from = SNP_id, 
                                   values_from = count,
                                   values_fill = 0))

counts <- apply(counts[,-1], 1, as.numeric)
colnames(counts) <- as.character(paste("pool",1:50, sep=""))
rownames(counts) <- unique(paste("SNP", alefeq_fd1$POS, sep="_"))

# SNP info object
snp.info <- data.table(Chromosome=1,
                       Position=unique(paste("SNP", alefeq_fd1$POS, sep="_")),
                       RefAllele=unique(alefeq_fd1$REF),
                       AltAllele=unique(alefeq_fd1$ALT))

# Make pool
pool <- new("pooldata",
            npools = dim(counts)[1],
            nsnp = dim(counts)[2], 
            refallele.readcount = counts, 
            readcoverage = coverages,
            snp.info = snp.info,
            poolsizes = 50*2,
            poolnames = as.character(paste("pool",1:50, sep="")))


is.pooldata(pool) <- TRUE
#### Pairwise FST Function ####

t <- compute.pairwiseFST(pool,
                         method = "Anova", 
                         verbose = TRUE)
t@values



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

# Make linear model for PCA
lm.pc1year <- summary(lm(data = pca.dt, Dim.1 ~ Year))
lm.pc2year <- summary(lm(data = pca.dt, Dim.2 ~ Year))
lm.pc3year <- summary(lm(data = pca.dt, Dim.3 ~ Year))

# Forloop through comparisons across and within Years
outfile <- foreach(j=1:dim(outfile)[1], .errorhandling = "pass", .combine = "rbind") %dopar% {

  # Progress message
  print(paste("Comparison set:", j, sep=" "))
  
  # Comparison set  
  samps <- c(as.character(compare$Sample1[j]), 
             as.character(compare$Sample2[j]))
  
  # Rbind pools
  base_obj <- data.table(rbind(data.table(alefeq_fd1[vcf==samps[1]], pool = "p1"),
                               data.table(alefeq_fd1[vcf==samps[2]], pool = "p2")) %>%
                mutate(count = as.numeric(AF.pool)*as.numeric(coverage),
                       SNP_id = paste("sim", POS, sep = "_"))) 
  
  # Long to wide coverages matrix
  coverages <- data.table(base_obj %>% 
                   select(pool, SNP_id, coverage) %>%
                   pivot_wider(names_from = SNP_id, values_from = coverage))
  
  # Get counts matrix
  counts <- data.table(base_obj %>% 
                   select(pool, SNP_id, count) %>%
                   pivot_wider(names_from = SNP_id, values_from = count))
  
  # Make pool
  pool <- new("pooldata",
              npools = dim(counts)[1],
              nsnp = dim(counts)[2], 
              refallele.readcount = t(data.matrix(counts)), 
              readcoverage = t(data.matrix(coverages)),
              poolsizes = 50*2,
              poolnames = c(unique(alefeq_fd[vcf==samps[1]]$vcf),
                            unique(alefeq_fd[vcf==samps[2]]$vcf)))

  # Calculate FST
  fst.out <- computeFST(pool, method = "Anova")
  
  # Output to file
  dt <- data.table(samp1=as.character(compare$Sample1[[j]]),
                   samp2=as.character(compare$Sample2[[j]]),
                   FST=fst.out$FST,
                   iteration=j)
  
  # Finish j
  return(dt)
  
}

# Compile and summarize data
total <- data.table(data.table(outfile %>%
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
