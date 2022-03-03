### JCBN
total_snps=739157
step=1200

seq(from=1, to=total_snps, by= step) %>%
  data.frame(start=.) %>%
  mutate(end = start+step-1 ) ->
  guide_file

guide_file[dim(guide_file)[1], "end"] = total_snps

write.table(guide_file,
            file = "./3.annotation_guide_file.txt", 
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")
