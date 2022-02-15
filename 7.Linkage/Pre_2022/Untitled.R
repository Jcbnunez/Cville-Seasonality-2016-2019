#####
## haplovalidate
library(haploReconstruct)  
library(psych)
library(stringr)
library(data.table)
library(haplovalidate)

load("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")

glm.out %>% 
  filter(
    rnp.clean<0.05,
    mod == "aveTemp",
    chr == "2L")  %>%
  select(chr, pos, b) ->
  glm05

names(glm05)[3] = "score"
##
cands <- merge(comp_cases_id_dedupl_sort,glm05,by=.(chr,pos))
saveRDS(cands,"cands.rds")




names(out_snps) = c("chr", "pos", "score")

parameters <- get.mncs.win(cands.all,
                           out_snps,
                           wins=seq(0.1,10,0.05),
                           mncs=0.01)

print(parameters)

parameters$win = "0.1"

happy <- haplovalidate(cands.all,
                       out_snps,
                       parameters,
                       repl,
                       gens)

plot.haplovalidate(blocks=happy$dominant_haplotypes,cmh,title="My beautiful haplotype blocks",label=TRUE)

