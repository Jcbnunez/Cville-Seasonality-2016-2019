# Outputs compiled genomic data from SLiM and PLINK
# Connor Murray 10.5.2021
# ijob -A berglandlab_standard --mem=10G -p standard -c 4
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(PopGenome))
suppressMessages(library(foreach))
suppressMessages(library(forcats))
suppressMessages(library(gmodels))
suppressMessages(library(poolSeq))
suppressMessages(library(poolfstat))

# Executable in command line
arg <- commandArgs(TRUE)
path.name <- "/dev/shm/csm6hg/1/20" # pathway to output
file <- '135959205_.acount' # seed for vcf
out.name <- "/project/berglandlab/connor/slim_bottleneck/freq.bottleneck.csv" # the output file name

# VCF output folder
setwd(path.name)

# All frequency data
filenames <- list.files(pattern = paste(file, "$", sep=""))

# Register cores
doParallel::registerDoParallel(cores = 4)

# ReadVCF forloop - calculate diversity statistics.
out <- foreach(i = 1:length(filenames), .combine="rbind", .errorhandling="remove") %dopar% {
  
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

# Forloop through comparisons across and within Years
outfile <- foreach(j=1:dim(outfile)[1], .errorhandling = "pass", .combine = "rbind") %dopar% {
  
  # Progress message
  print(paste("Comparison set:", j, sep=" "))
  
  # Comparison set  
  samps <- c(as.character(compare$Sample1[j]), 
             as.character(compare$Sample2[j]))
  
  # Extract allele freqeuncy info
  alefeq_fd <- data.table(out[vcf %in% c(samps)] %>% 
              group_by(vcf) %>% 
              select(POS, AF, vcf))

  # Sample 1 - get allele frequencies  
  alefeq_fd1 <- alefeq_fd[vcf==samps[1]] 
  alefeq_fd1 %>% .$POS %>% table -> a1
  a1[which(a1 == 2)] -> names_of_dup
  dt1 <- alefeq_fd1[which(!alefeq_fd1$POS %in% names(names_of_dup) ),]
  alefeq_fd1 <- dt1[which(!dt1$POS %in% names(names_of_dup) ),] %>% select(POS, AF)
  
  # Create sample 1 pool
  pool1 <- sample.alleles(alefeq_fd1$AF, size=60, mode="coverage")
  
  # Sample 2 - get allele frequencies  
  alefeq_fd2 <- alefeq_fd[vcf==samps[2]] 
  alefeq_fd2 %>% .$POS %>% table -> a2
  a2[which(a2 == 2)] -> names_of_dup
  dt2 <- alefeq_fd2[which(!alefeq_fd2$POS %in% names(names_of_dup) ),]
  alefeq_fd2 <- dt2[which(!dt2$POS %in% names(names_of_dup) ),] %>% select(POS, AF)
  
  # Create sample 2 pool
  pool2 <- sample.alleles(alefeq_fd2$AF, size=60, mode="coverage")
  
  # Rbind pools
  rbind(
    data.frame(pool1, pool = "p1", alefeq_fd1),
    data.frame(pool2, pool = "p2", alefeq_fd2)) %>%
    mutate(count = p.smpld*size,
           SNP_id = paste("sim", POS, sep = "_")) ->
    base_obj
  
  # Long to wide coverages matrix
  coverages = base_obj %>%
    dcast(pool~SNP_id, value.var = "size")
  
  # Get counts matrix
  counts  = base_obj %>%
    dcast(pool~SNP_id, value.var = "count" )
  
  # Make pool
  pool <- new("pooldata",
              npools=dim(counts)[1],
              nsnp=dim(counts)[2], 
              refallele.readcount=as.matrix(apply(counts, 1, as.numeric)), 
              readcoverage=as.matrix(apply(coverages, 1, as.numeric)),
              poolsizes=50 * 2,
              poolnames=c(unique(alefeq_fd[vcf==samps[1]]$vcf),
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

# Save output
saveRDS(object = outfile, file = out.name)

# Compile and summarize data
total <- data.table(outfile %>%
                      left_join(out %>% 
                                select(vcf, nMax, nMin, bottle, seed, sim.gen1=sim.gen) %>% 
                                distinct(), 
                              by=c("samp1"="vcf")) %>%
                      left_join(out %>% 
                                select(vcf, sim.gen2=sim.gen) %>% 
                                distinct(), 
                              by=c("samp2"="vcf")))


png("/project/berglandlab/connor/slim_bottleneck/fst.pairwise.png")

ggplot(total,aes(x=as.numeric(sim.gen1),
                   #fct_reorder(sim.gen1, as.numeric(sim.gen1)), 
                 y=fct_reorder(sim.gen2, as.numeric(sim.gen2)), 
                 fill=FST)) +
  geom_raster() +
  labs(x="Generation 1", 
       y="Generation 2", 
       title="Constant population size - 100,000") +
  theme_classic() +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        #legend.position = "none", 
        legend.background = element_blank(),
        axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

dev.off()
  
png("/project/berglandlab/connor/slim_bottleneck/fst.pairwise.walk.png")

ggplot(total,aes(x=as.numeric(sim.gen1),
                 y=as.numeric(FST))) +
  geom_line() +
  geom_point() +
  facet_wrap(~sim.gen2) +
  labs(x="Generation",
       y="FST",
       title="Constant population size - 100,000") +
  theme_classic() +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        #legend.position = "none", 
        legend.background = element_blank(),
        axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

dev.off()

# Write final output
write_csv(x=total, 
          path="/project/berglandlab/connor/slim_bottleneck/pairwise.fst.csv", 
          append=T, 
          quote_escape=F)
