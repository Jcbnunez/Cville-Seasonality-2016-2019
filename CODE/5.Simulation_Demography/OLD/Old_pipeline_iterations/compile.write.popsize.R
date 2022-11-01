# Compile population data from SLiM
# 8.23.2020

# Libraries
suppressMessages(library(foreach))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

# Executable commands
arg <- commandArgs(TRUE)
path.name <- arg[1]
file <- arg[2]
out.name <- arg[3]

# Output folder
setwd(path.name)

# Extract all population files
tbl <- list.files(pattern = paste(file, "$", sep="")) 

# Read in files as data table
out <- foreach(i=1:length(tbl), .combine="rbind", .errorhandling="remove") %do% {
  out <- data.table(read.table(tbl[i]), file=tbl[i])
  print(paste(round(i/length(tbl)*100, digits = 3), "% Complete", sep=""))
  return(out)
}

# Fix header
# colnames(out)[1:8] <- c("pop", "generation", "gen.fin", "pop.size", "seed", "EG", "K", "nSamp")

# Output
write_csv(x=out, path=out.name, col_names=F, append=T, quote_escape=F)
