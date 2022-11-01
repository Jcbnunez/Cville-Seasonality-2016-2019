# Outputs compiled genomic data from SLiM and PLINK
# Connor Murray 10.5.2021
# ijob -A berglandlab_standard --mem=10G -p standard -c 4
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(PopGenome))
suppressMessages(library(foreach))
suppressMessages(library(gmodels))
suppressMessages(library(poolSeq))
suppressMessages(library(poolfstat))

# Executable in command line
arg <- commandArgs(TRUE)
path.name <- "/dev/shm/csm6hg/1/100" # pathway to output
file <- '83197612_.acount' # seed for vcf
out.name <- "/project/berglandlab/connor/slim_bottleneck/freq.bottleneck.csv" # the output file name

# VCF output folder
setwd(path.name)

# All frequency data
filenames <- list.files(pattern = paste(file, "$", sep=""))

# Register cores
# doParallel::registerDoParallel(cores = 4)

# ReadVCF forloop - calculate diversity statistics.
out <- foreach(i = 1:length(filenames), .combine="rbind", .errorhandling="remove") %do% {
  
  i=1
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
                   bottle = tstrsplit(filenames[i], "_")[[9]],
                   seed = tstrsplit(filenames[i], "_")[[10]],
                   iteration = i) %>% 
                mutate(Year = case_when(
                              sim.gen %in% c(1:18) ~ "Year_1",
                              sim.gen %in% c(19:35) ~ "Year_2",
                              sim.gen %in% c(36:51) ~ "Year_3")) %>% 
                mutate(Year.samp = paste(Year, nMax, nMin, bottle, sep="_")))
 
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





# Forloop through comparisons across and within Years
foreach(j=1:length(vcf.list), .combine = "rbind") %do% {

  j=20
  
  print(j)
  
  # Simulate pool-seq coverage
  simTraj <- out[vcf %in% samps]$MAF, 
                     Ne=(0.001/(4*2.89e-09)), 
                     t=c(0))
  
  af <- sample.alleles(simTraj, 
                       size=60, 
                       mode="coverage")
  
  simTraj <- matrix(af$p.smpld, 
                    nrow=nrow(simTraj), 
                    dimnames=dimnames(simTraj))
  
  simCov <- matrix(af$size, 
                   nrow=nrow(simTraj), 
                   dimnames=dimnames(simTraj))
  
  # Allele counts matrix
  ad.matrix = as.matrix(out[vcf %in% samps] %>% 
                          group_by(vcf) %>% 
                          mutate(row=row_number()) %>% 
                          select(seed, C2, vcf, row) %>% 
                          pivot_wider(names_from = vcf, 
                                      values_from = C2) %>% 
                          select(-row, -seed))
  
  # Add rownames to matrix
  rownames(ad.matrix) <- 1:dim(ad.matrix)[1]
  
  # Coverage matrix
  rd.matrix = as.matrix(out[vcf %in% samps] %>% 
                          cbind(simCov) %>% 
                          group_by(vcf) %>% 
                          mutate(row=row_number()) %>% 
                          select(seed, F0, vcf, row) %>% 
                          pivot_wider(names_from = vcf, 
                                      values_from = F0) %>% 
                          select(-row, -seed))
  
  # Add rownames to matrix
  rownames(rd.matrix) <- 1:dim(rd.matrix)[1]
  
  # Make pool
  pool <- new("pooldata",
              npools=dim(ad.matrix)[2],
              nsnp=dim(ad.matrix)[1], 
              refallele.readcount=ad.matrix, 
              readcoverage=rd.matrix,
              poolsizes=pool_sizes * 2)
  
  # Calculate FST
  fst.out <- computeFST(pool, method = "Anova")
  
  # Output to file
  outfile$samp1[j] = as.character(compare$Sample1[[j]])
  outfile$samp2[j] = as.character(compare$Sample2[[j]])
  outfile$FST[j] = fst.out$FST
  
  # Finish j
  return(fst.out)
  
}

# Pairwise comparisons across and within years
samps <- c("Year_1", "Year_2", "Year_3")
compare <- data.frame(expand.grid(samps, samps))
colnames(compare) <- c("Sample1", "Sample2")

# Remove pairwise duplicates
compare <- data.table(compare[!duplicated(t(apply(compare, 1, sort))),])
rownames(compare) <- 1:dim(compare)[1]

# Forloop through comparisons across and within Years
foreach(j=1:length(vcf.list), .combine = "rbind") %do% {
  
  # Simulate pool-seq coverage
  simTraj <- wf.traj(out[vcf %in% vcf.list[i]]$MAF, 
                     Ne=(0.001/(4*2.89e-09)), 
                     t=c(0))
  
  af <- sample.alleles(simTraj, 
                       size=60, 
                       mode="coverage")
  
  simTraj <- matrix(af$p.smpld, 
                    nrow=nrow(simTraj), 
                    dimnames=dimnames(simTraj))
  
  simCov <- matrix(af$size, 
                   nrow=nrow(simTraj), 
                   dimnames=dimnames(simTraj))
  
  # Allele counts matrix
  ad.matrix = as.matrix(out[vcf %in% vcf.list[i]] %>% 
                          group_by(vcf) %>% 
                          mutate(row=row_number()) %>% 
                          select(seed, C2, vcf, row) %>% 
                          pivot_wider(names_from = vcf, 
                                      values_from = C2) %>% 
                          select(-row, -seed))
  
  # Add rownames to matrix
  rownames(ad.matrix) <- 1:dim(ad.matrix)[1]
  
  # Coverage matrix
  rd.matrix = as.matrix(out[vcf %in% vcf.list[i]] %>% 
                          cbind(simCov) %>% 
                          group_by(vcf) %>% 
                          mutate(row=row_number()) %>% 
                          select(seed, F0, vcf, row) %>% 
                          pivot_wider(names_from = vcf, 
                                      values_from = F0) %>% 
                          select(-row, -seed))
  
  # Add rownames to matrix
  rownames(rd.matrix) <- 1:dim(rd.matrix)[1]
  
  # Make pool
  pool <- new("pooldata",
              npools=dim(ad.matrix)[2],
              nsnp=dim(ad.matrix)[1], 
              refallele.readcount=ad.matrix, 
              readcoverage=rd.matrix,
              poolsizes=pool_sizes * 2)
  
  # Calculate FST
  fst.out <- computeFST(pool, method = "Anova")
  
  # Output to file
  outfile$samp1[j] = samps[1]
  outfile$samp2[j] = samps[2]
  outfile$FST[j] = fst.out$FST
  
  # Finish j
  return(fst.out)
  
}



# Forloop through comparisons across and within Years
foreach(j=1:dim(compare)[1], .combine = "rbind") %do% {

  # Samples to calculate FST
  samps <- c(as.character(compare$Sample1[[j]]), 
             as.character(compare$Sample2[[j]]))
  
  # Simulate pool-seq coverage
  simTraj <- wf.traj(out[Year %in% samps]$MAF, 
                       Ne=(0.001/(4*2.89e-09)), 
                       t=c(0))
    
  af <- sample.alleles(simTraj, 
                         size=60, 
                         mode="coverage")
    
  simTraj <- matrix(af$p.smpld, 
                      nrow=nrow(simTraj), 
                      dimnames=dimnames(simTraj))
    
  simCov <- matrix(af$size, 
                     nrow=nrow(simTraj), 
                     dimnames=dimnames(simTraj))
    
  # Allele counts matrix
  ad.matrix = as.matrix(out[Year %in% samps] %>% 
                             group_by(vcf) %>% 
                             mutate(row=row_number()) %>% 
                             select(seed, C2, vcf, row) %>% 
                             pivot_wider(names_from = vcf, 
                                         values_from = C2) %>% 
                             select(-row, -seed))
  
  # Add rownames to matrix
  rownames(ad.matrix) <- 1:dim(ad.matrix)[1]
    
  # Coverage matrix
  rd.matrix = as.matrix(out[Year %in% samps] %>% 
                            cbind(simCov) %>% 
                            group_by(vcf) %>% 
                            mutate(row=row_number()) %>% 
                            select(seed, F0, vcf, row) %>% 
                            pivot_wider(names_from = vcf, 
                                        values_from = F0) %>% 
                            select(-row, -seed))
    
  # Add rownames to matrix
  rownames(rd.matrix) <- 1:dim(rd.matrix)[1]
    
  # Make pool
  pool <- new("pooldata",
              npools=dim(ad.matrix)[2],
              nsnp=dim(ad.matrix)[1], 
              refallele.readcount=ad.matrix, 
              readcoverage=rd.matrix,
              poolsizes=pool_sizes * 2)
    
  # Calculate FST
  fst.out <- computeFST(pool, method = "Anova")

  # Output to file
  outfile$samp1[j] = samps[1]
  outfile$samp2[j] = samps[2]
  outfile$FST[j] = fst.out$FST
  
# Finish j
return(fst.out)

}


# Compile and summarize data
total <- data.table(dt %>% 
                      group_by() %>% 
                      summarise(varAF=var(AF),
                                meanAF=mean(AF),
                                uciAF=ci(AF, confidence = 0.975)[3],
                                lciAF=ci(AF, confidence = 0.975)[2]))

# Write final output
write_csv(x=out, path=out.name, append=T, quote_escape=F)
