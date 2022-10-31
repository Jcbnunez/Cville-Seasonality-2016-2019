# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

fl <- list.files("/scratch/aob2x/summarized_dest_glm/inversionEnrichment", ".Rdata", full.names=T)

invEn <- foreach(fl.i=fl, .combine="rbind")%do%{
  #fl.i <- fl[1]
  message(fl.i)
  load(fl.i)
  fet.invs
}

save(invEnv, file="~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/2.InversionEnrichmentAnalysis/inversion_enrichment_output.Rdata")
