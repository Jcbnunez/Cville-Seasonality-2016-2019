# ijob -A berglandlab_standard -c2 -p dev --mem=10G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])

#jobId=455 # "temp.max;4;5.Cville"
# jobId=2

### libraries
  library(data.table)
  #library(gdata)
  library(foreach)
  library(doMC)
  #registerDoMC(10)

  registerDoMC(30)

### load guide file
   # load(fl[1])
   # mod_var <- glm.out[,.N,list(variable, mod, cluster)]
   # save(mod_var, file="/project/berglandlab/alan/environmental_ombibus_global_permV2/mod_var.redoPerm.Rdata")

   load("/project/berglandlab/alan/environmental_ombibus_global/mod_var.redoPerm.Rdata")
   mainDir <- "/scratch/aob2x/environmental_ombibus_global_permV2"
   mvi <- paste(mod_var[jobId]$variable, mod_var[jobId]$mod, mod_var[jobId]$cluster, sep=";")
   mvi

### check to see if output file already exists
 outDir <- paste("/project/berglandlab/alan/environmental_ombibus_global_permV2/glmOut/", mvi, sep="")
 dir.create(file.path(outDir), showWarnings = FALSE)

 if(file.exists(file=paste(outDir, "/", mvi, ".glmRNP.Rdata", sep=""))) {
   message( fileInfo <- file.info(file=paste(outDir, "/", mvi, ".glmRNP.Rdata", sep="")))
   message("done")
   q("no")
 }

### cycle through raw data and save mvi
### combine with snp identity
  load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
  setkey(snp.dt, id)

  use <- apply(snp.dt[,c("VA_ch"), with=F],
                1, any)
  snp.dt <- snp.dt[use]
  setkey(snp.dt, id)


### split output and get summaries
  fl <- list.files("/scratch/aob2x/environmental_ombibus_global_permV2",full.names=T)
  #fl <- fl[1:100]

  glm.out <- foreach(fl.i=fl, .errorhandling="remove")%do%{
    # fl.i <- fl[1]
    message(paste(which(fl.i==fl), length(fl), sep=" / "))

    load(fl.i)
    #glm.out[cluster=="2.North_America_E", cluster:="2.North_America_I95"]
    glm.out <- merge(glm.out, snp.dt, by.x="variant.id", by.y="id", all.x=T)
    setkey(glm.out, variable, mod, cluster)

    tmp <- glm.out[J(mod_var[jobId])]

    ### return mvi slice
      return(tmp)
  }
  #glm.out <- rbindlist(glm.out)
  #dim(glm.out)
  length(glm.out)


  glm.out.ag <- glm.out[,list(rnp_r=rank(p_lrt, ties="random")/length(p_lrt), rnp=rank(p_lrt)/length(p_lrt), chr, pos, inv=invName!="none"), list(perm)]

  setkey(glm.out, chr, pos, perm)
  setkey(glm.out.ag, chr, pos, perm)

  glm.out_new <- merge(glm.out, glm.out.ag)

### load in old p-values just for comparison
  load(paste("/project/berglandlab/alan/environmental_ombibus_global/", mvi, "/", mvi, ".glmRNP.Rdata", sep=""))
  setkey(glm.out, variant.id, perm)
  setkey(glm.out_new, variant.id, perm)

  glm.out <- merge(glm.out_new, glm.out[,c("variant.id", "perm", "p_lrt")])
  rm(glm.out_new)
  dim(glm.out)

### output raw GLM results
  outDir <- paste("/project/berglandlab/alan/environmental_ombibus_global_permV2/glmOut/", mvi, sep="")
  dir.create(file.path(outDir), showWarnings = FALSE)

  save(glm.out, file=paste(outDir, "/", mvi, ".glmRNP.Rdata", sep=""))
  #save(glm.out2, file=paste("~", "/", mvi, ".glmRNP.Rdata", sep=""))
  # load(paste("~", "/", mvi, ".glmRNP.Rdata", sep="")) ; glm.out <- glm.out2

#### summarize
  p_thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-1)))$Var1

  glm.enrich <- foreach(thr=p_thrs, .errorhandling="remove", .combine="rbind")%dopar%{
    # thr <- .05
    tmp <- glm.out[,list(new=mean(p_lrt.x<=thr), orig=mean(p_lrt.y<=thr), thr, .N), list(perm, mod, variable, chr, invName)]

  }

  glm.enrich[perm==0][order(thr)][which.max(new/thr)]
  glm.enrich.ag <- glm.enrich[,list(rr_new=mean(new[perm==0]/new[perm!=0]), rr_old=mean(orig[perm==0]/orig[perm!=0]),
                                    perm_new=mean(new[perm!=0]/thr), perm_old=mean(orig[perm!=0]/thr),
                                    obs_new=mean(new[perm==0]/thr), obs_old=mean(orig[perm==0]/thr)), list(thr, chr, inv=invName!="none")]
  glm.enrich.ag[order(rr_new)][thr==.005]

  outDir <-"/project/berglandlab/alan/environmental_ombibus_global_permV2/glmEnrich/"
  dir.create(file.path(outDir), showWarnings = FALSE)
  save(glm.enrich, glm.enrich.ag, file=paste(outDir, "/", mvi, ".glm_enrich.Rdata", sep=""))


### summarize v2
  grid <- 0.001
  my_seq = data.table(min_p=seq(from=0, to=1-grid, by=grid), max_p=seq(from=grid, to=1, by=grid))
  glm.out[,inv:=invName!="none"]
  setkey(glm.out, perm, chr, inv)
  glm.bins <- foreach(perm.i=unique(glm.out$perm), .combine="rbind")%dopar%{
    foreach(chr.i=unique(glm.out$chr), .combine="rbind")%do%{
      foreach(inv.i=unique(glm.out$inv), .combine="rbind")%do%{
        message(paste(chr.i, inv.i, perm.i, sep=" / "))
        # chr.i <- "2L"; inv.i <- TRUE; perm.i<-0
        ### new
          new_tmp <- glm.out[perm==perm.i][chr==chr.i][inv==inv.i][my_seq, .(N = .N), on = .(p_lrt.x > min_p, p_lrt.x <= max_p), by = .EACHI]
          new_tmp[,chr:=chr.i]
          new_tmp[,inv:=inv.i]
          new_tmp[,perm:=perm.i]
          new_tmp[,perm_method:="new"]
          setnames(new_tmp, 1:2, c("p_lrt.low", "p_lrt.high"))

        ### old
          old_tmp <- glm.out[perm==perm.i][chr==chr.i][inv==inv.i][my_seq, .(N = .N), on = .(p_lrt.y > min_p, p_lrt.y <= max_p), by = .EACHI]
          old_tmp[,chr:=chr.i]
          old_tmp[,inv:=inv.i]
          old_tmp[,perm:=perm.i]
          old_tmp[,perm_method:="old"]
          setnames(old_tmp, 1:2, c("p_lrt.low", "p_lrt.high"))

        rbind(new_tmp, old_tmp)
      }
    }
  }
  outDir <-"/project/berglandlab/alan/environmental_ombibus_global_permV2/glmBins/"
  dir.create(file.path(outDir), showWarnings = FALSE)
  save(glm.bins, file=paste(outDir, "/", mvi, ".glm_bins.Rdata", sep=""))





# #### how many jobs are missing?
#   fl <- list.files("/scratch/aob2x/environmental_ombibus_global_permV2",full.names=T)
#   fld <- data.table(job=tstrsplit(fl, "/")[[5]])
#   fld[,num:=gsub("job", "", job)]
#   fld[,real:=as.numeric(gsub(".Rdata", "", num))]
#
#   c(1:5000)[!c(1:5000)%in%fld$real]
#
# paste(c(1:5000)[!c(1:5000)%in%fld$real] , collapse=",")
#
# prop.table(table(glm.out[,list(pr=mean(p_lrt[perm==0]<p_lrt[perm!=0]), p=p_lrt[perm==0]), list(variant.id, variable, cluster)]$pr>=.95))
#
#
# "12,14,24,26,85,87,127,128,131,133,135,136,138,140,141,143,144,145,146,147,148,149,150,153,154,155,156,157,158,161,163,165,167,168,169,172,173,174,175,176,1
# 77,179,180,181,182,183,184,185,186,188,189,190,191,192,193,195,196,199,200,201,202,204,205,206,207,209,210,211,212,214,215,216,217,218,219,221,222,223,224,225,2
# 26,228,229,230,231,233,234,235,236,237,238,239,240,241,242,243,244,245,246,248,249,251,252,253,254,255,256,257,258,277,278,279,288,289,290,291,293,307,392,405,4
# 16,417,418,422,423,424,435,436,439,440,441,443,445,446,452,453,473,476,478,479,481,482,483,485,520,570,571,572,574,575,577,578,579,580,582,585,588,590,592,593,5
#
# 99,600,601,604,606,607,608,609,610,611,613,614,615,694,708,709,711,712,717,734,736,737,738,739,740,741,743,795,814,826,827,850,987,1039,1042,1043,1044,1045,1046
# ,1049,1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1061,1062,1063,1064,1070,1072,1094,1096,1107,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193
# ,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1220,1221,1222,1223,1224,1225,1226,1227,1228,1229,1231,1234,1237,1239,1241,1244,1245,1246,1247,1248
# ,1249,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1264,1265,1266,1268,1269,1270,1271,1273,1274,1276,1277,1278,1279,1282,1283,1284,1285,1287,1289
# ,1291,129
