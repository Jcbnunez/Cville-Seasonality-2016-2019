### Create the guide file for boot and resamp
### 
library(permute)

####
id_vec = 1:250

SNP_sampler <- c(seq(from = 100, to = 1000, by = 100),
                 seq(from = 1000, to = 20000, by = 1000))

expand.grid(id_vec, SNP_sampler, 28) -> guide_file

##shuffle to avoid batch effect
write.table(guide_file[shuffle(guide_file[1]),], 
            file = "guide_file_2022_boot_resamp.28Cov.varCov.txt",
            append = FALSE,
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)

