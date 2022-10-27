library(SeqArray)

args = commandArgs(trailingOnly=TRUE)

vcf.fn=args[[1]]

gds.fn=gsub(".vcf", ".gds", vcf.fn)

vcf.fn=paste(vcf.fn, ".gz", sep="")

seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA")