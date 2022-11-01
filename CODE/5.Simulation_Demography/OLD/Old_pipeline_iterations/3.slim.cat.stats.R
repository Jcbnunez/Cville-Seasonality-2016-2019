# Outputs compiled genomic data from SLiM and PLINK
# Connor Murray 10.29.2021
# ijob -A berglandlab_standard --mem=10G -p standard -c 4
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(foreach))

# Executable in command line
arg <- commandArgs(TRUE)
path.name <- arg[1] # pathway to output
file <- arg[2] # seed for vcf
out.name <- arg[3] # the output file name

#path.name="/dev/shm/csm6hg/1/1"
#file="7976955_."
#out.name="/dev/shm/csm6hg/1/1/output"

# VCF output folder
setwd(path.name)

# All frequency data
filenames <- list.files(pattern = file)

# Stats used
stats <- c("het", "LROH", "hwe", "Tajima.D", "windowed.pi")

# Register cores
doParallel::registerDoParallel(cores = 4)

# Read forloop 
out <- foreach(i = 1:length(stats), .errorhandling="remove") %dopar% {
  
  # Filenames
  files <- filenames[which(filenames %like% stats[i]==TRUE)]
  
  # Progress message
  print(paste("Processing: ", stats[i], sep=""))
  
  # Nested - individuals
  out.cat <- foreach(j = 1:length(files), .combine="rbind", .errorhandling="remove") %do% {
  
    # Load allele frequency information
    dt <- data.table(fread(files[j], header = TRUE), 
                     vcf = files[j],
                     slurm_ID = tstrsplit(files[j], "_")[[3]],
                     nSamp = tstrsplit(files[j], "_")[[4]],
                     sim.gen = tstrsplit(files[j], "_")[[5]],
                     fin.gen = tstrsplit(files[j], "_")[[6]],
                     nMax = tstrsplit(files[j], "_")[[7]],
                     nMin = tstrsplit(files[j], "_")[[8]],
                     replicate = tstrsplit(files[j], "_")[[9]],
                     seed = tstrsplit(files[j], "_")[[10]],
                     iteration = j)
    
    # Progress message
    print(paste("SLURM_ID-", 
                tstrsplit(files[j], "_")[[3]], 
                "-Iteration-", i, ": ", 
                round((j/length(files))*100, digits=2), 
                "%", " complete", sep=""))
  
    # Finish
    return(dt)
  }
  
  # Write output
  write_delim(out.cat,
              file=paste(out.name,
                         stats[i],
                         sep="."), 
              delim = "\t")
  
}
