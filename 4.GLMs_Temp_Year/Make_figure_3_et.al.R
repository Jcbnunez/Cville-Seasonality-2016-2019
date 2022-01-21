# Load modules in Rivanna
#module load gcc/7.1.0
#module load openmpi/3.1.4
#module load R/4.1.1
## load R
#R


### libraries
  library(tidyverse)
  library(data.table)
  library(tidyr)
  library(viridis)
  library(patchwork)
  library(magrittr)
  library(RColorBrewer)
  library(rstatix)
  
###########################################
### load & define Inversion and Ace ROI ###
###########################################

  output_results_window <- "/project/berglandlab/thermal_glm_dest/window_analysis_output.nested.qb.Rdata"
  inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"
  
    ### load suppl data
  inv.dt <- fread(inversion_map)
  setnames(inv.dt, "chrom", "chr.x")

  ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)


#### load windows
  #load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/3.windowAnalysis/window_analysis_output.Rdata")
  # scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.Rdata ~/.
  #load("~/window_analysis_output.nested.Rdata")

  #scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.qb.Rdata ~/.
  load(output_results_window)

### get sig thresholds from perms
 ## win.minp.ag <- win.out[pr==0.01 & nSNPs>100,
 ##         list(q05=quantile(wZa.p, 0.025, na.rm=T), q5=quantile(rbinom.p, .05, na.rm=T), .N),
 ##         list(perm, mod, locality)]$N %>% mean
 ## win.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp+year_factor"]

   ## win.minp.ag.ag <- win.minp.ag[perm!=0, list(q5.minp=quantile(minp, .05, na.rm=T), min.q5=min(q5, na.rm=T)), list(locality, mod, chr.x)]

### v2 -- stablish significance
###
win.minp.ag <- win.out[pr==0.05 & nSNPs>100 & perm!=0,
        list(lci=quantile(rbinom.p, 0.025, na.rm=T), uci=quantile(rbinom.p, .975, na.rm=T), .N),
        list(mod, chr.x=chr.x, locality, win.i, start, end)]

###win.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp+year_factor"]

###win.minp.ag.ag <- win.minp.ag[perm!=0, 
## list(q5.minp=quantile(minp, .05, na.rm=T), min.q5=min(q5, na.rm=T)), list(locality, mod, chr.x)]


### Part 1.  ---> the manhattan plot
  
### basic MH plot
  mh.plot.wza <-
  ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  geom_ribbon(data=win.minp.ag[mod=="aveTemp+year_factor"],
              aes(x=(start/2 +end/2)/1e6, ymin=-1*log10(uci), ymax=-1*log10(lci)),
              color="grey", fill="grey", alpha=0.75) +
  geom_point(data=win.out[mod=="aveTemp+year_factor"][pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rbinom.p),
                color=(rnp.pr)), size=.95) +
  geom_hline(yintercept=-log10(.01/1800)) +
  #geom_hline(data=win.minp.ag.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
  facet_grid(locality~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0),
         legend.box = "vertical",
         legend.key.size=unit(1/8, 'in'),
         legend.text=element_text(size=6),
         legend.title=element_text(size=8)) +
  labs(color="Prop. top 1%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab("-log10(Window p)")

  ggsave(mh.plot.wza, file="nested_qb.png", h=8.5, w=11)

  ### Part 2.  ---> the entire density plot for the supplement
  
### density plots using ggplot geom_density
### Using a probability cutoff of 0.01
  pr.densPlot <-
  ggplot() +
  geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_density(data=win.out[locality=="VA_ch"][perm!=0][nSNPs>=100][pr==0.01],
            aes(x=qlogis(rnp.pr),
                group=interaction(perm, (invName=="none")),
                linetype=as.factor(invName=="none")),
            color="grey") +
  geom_density(data=win.out[locality=="VA_ch"][perm==0][nSNPs>=100][pr==0.01],
              aes(x=qlogis(rnp.pr),
                  group=interaction(perm, (invName=="none")),
                  linetype=as.factor(invName=="none")),
              color="black", size=1) +
  facet_grid(chr.x~mod+as.factor(invName=="none")) +
  theme_bw() +
  theme(legend.position="none")


  wza.densPlot <-
  ggplot() +
  #geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_density(data=win.out[locality=="VA_ch"][perm!=0][nSNPs>=100][pr==0.01],
            aes(x=-wZa,
                group=interaction(perm, (invName=="none")),
                linetype=as.factor(invName=="none")),
            color="grey") +
  geom_density(data=win.out[locality=="VA_ch"][perm==0][nSNPs>=100][pr==0.01],
              aes(x=-wZa,
                  group=interaction(perm, (invName=="none")),
                  linetype=as.factor(invName=="none")),
              color="black", size=1) +
  facet_grid(chr.x~mod+as.factor(invName=="none")) +
  theme_bw() +
  theme(legend.position="none")


  mega <-
  pr.densPlot + wza.densPlot +
  plot_annotation(tag_levels="A")


  ggsave(mega, file="densPlot_nested.qb.split.pdf", h=10, w=10)


  ### Part 3.  ---> prepare figure 3
  ### 
  
  win.out[locality=="VA_ch"][nSNPs>=100][pr==0.01] -> data_va_p01
  #[perm!=0]

  data_va_p01$mod = gsub("year", "1.year", data_va_p01$mod)
  
  data_va_p01 %>%
    mutate(dat_type = case_when(perm == 0 ~ "real",
                                perm != 0 ~ "shuffle")) %>%
    group_by(dat_type, mod, invName, chr.x, i ) %>%
    summarize(mean_rnp.pr = mean(rnp.pr),
              sd_rnp.pr = sd(rnp.pr)) %>%
    mutate(inv_type = case_when(invName == "none" ~ "outside.inv",
                                invName != "none" ~ "inside.inv",
                                ))-> data_va_p01_sum
  
  ### add summarized Z test
  data_va_p01 %>%
    mutate(dat_type = case_when(perm == 0 ~ "real",
                                perm != 0 ~ "shuffle")) %>%
    group_by(dat_type, mod, invName, chr.x, i ) %>%
    summarize(mean_wza = mean(-wZa),
              sd_wza = sd(-wZa)) %>%
    mutate(inv_type = case_when(invName == "none" ~ "outside.inv",
                                invName != "none" ~ "inside.inv",
    ))-> data_va_p01_sum_wZa
  
    ## start plotting
    ggplot() +
    geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
    geom_density(data=data_va_p01_sum,
                 aes(x=qlogis(mean_rnp.pr),
                     #group=interaction(perm, (invName=="none")),
                     linetype=inv_type,
                 color=dat_type),
                 adjust = 1) +
      theme_classic() +
      scale_color_brewer(palette = "Set1") +
    facet_grid(mod~chr.x) -> part1
  
  
  ggplot() +
    geom_vline(xintercept=0, linetype="dashed", color="black") +
    geom_density(data=data_va_p01_sum_wZa,
                 aes(x=mean_wza,
                     #group=interaction(perm, (invName=="none")),
                     linetype=inv_type,
                     color=dat_type),
                 adjust = 1) +
    xlim(-5,40) +
    scale_colour_brewer(palette = "Set1") +
    theme_classic() +
    facet_grid(mod~chr.x) -> part2
  

  ggsave(part1/part2, file = "qlogis.wza.pdf", 
         width = 7,
         height = 5)

## Pretty wZa
  ggplot() +
    geom_density(data=data_va_p01_sum_wZa,
                 aes(x=mean_wza,
                     #group=interaction(perm, (invName=="none")),
                     fill=inv_type,
                     linetype=dat_type),
                 alpha = 0.5) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    xlim(-5,45) +
    facet_grid(mod~chr.x) -> pretty_wza
  
  ggsave(pretty_wza, file = "pretty_wza.pdf", 
         width = 8,
         height = 4)
  

  
# compare real data to simulations in a t test
  
 #####
# This makes the large plot of t-tests
 
 #data_va_p01$mod = gsub("year", "1.year", data_va_p01$mod)
  mod_list = list()
  mods_ops = c("aveTemp+1.year_factor", "1.year_factor")

  data_va_p01 %<>%
  mutate(inv_type = case_when(invName == "none" ~ "outside.inv",
                              invName != "none" ~ "inside.inv"))
                              
for(j in 1:2){
 data_va_p01 %>%
   filter(mod == mods_ops[j]) -> object_for_analysis_pre
 
 t_test_list = list()

 for(i in 1:100){
   
   if(i == 99){break}
   
   chrs = c("2L", "2R", "3L", "3R")
   inner_list = list()
   for(k in 1:4){
     print(k)
    
   inversion_list = list()
   invs_typ = c("inside.inv", "outside.inv")
   for(h in 1:2){  
   
    object_for_analysis_pre %>%
       filter(inv_type == invs_typ[h]) ->
       object_for_analysis
   
   object_for_analysis$rnp.pr[which(object_for_analysis$perm == 0 & 
                                        object_for_analysis$chr.x == chrs[k])] %>%
   qlogis() -> obs_pt
   obs_pt[!is.finite(obs_pt)] <- NA
   
   object_for_analysis$rnp.pr[which(object_for_analysis$perm == i & 
                                      object_for_analysis$chr.x == chrs[k])] %>%
     qlogis() -> shuff_pt
   shuff_pt[!is.finite(shuff_pt)] <- NA
   
   t.test(obs_pt, shuff_pt) -> t.test.out
   
   df_d <- data.frame(
     obs_d  = obs_pt,
     shuff_d = shuff_pt
   )
   
   df_d %>% 
     .[complete.cases(.),] %>%
     melt() %>% 
     cohens_d(value ~ variable, var.equal = TRUE) -> cohens_d
   
   outdf <- data.frame(matrix(ncol = 9, nrow = 1))
   colnames(outdf) <- c('mod','ith','chr', 'inv' , 'est' ,'stat', 'p.val', 'low_c', 'high_c')
   
   outdf$ith[1] = i
   outdf$mod[1] = mods_ops[j]
   outdf$inv[1] = invs_typ[h]
   outdf$chr[1] = chrs[k]
   outdf$est[1] = t.test.out$estimate[1]-t.test.out$estimate[2]
   outdf$stat[1] = t.test.out$statistic
   outdf$p.val[1] = t.test.out$p.value
   outdf$low_c[1] = t.test.out$conf.int[1]
   outdf$high_c[1] = t.test.out$conf.int[2]
   outdf$cohens_d[1] = cohens_d$effsize
   outdf$cohens_amg[1] = cohens_d$magnitude
   
   inversion_list[[h]] = outdf
   } # close h
   
   inner_list[[k]] = do.call(rbind, inversion_list)
   } # close k
   
   t_test_list[[i]] = do.call(rbind, inner_list)
 } # close i
 
 mod_list[[j]] =  do.call(rbind, t_test_list)
} ## close j 

  t_test_out = do.call(rbind, mod_list)
  
  save(t_test_out, file = "t_test_out.Rdata")
  
  load("t_test_out.Rdata")
  
  t_test_out %>%
    group_by(mod, inv, chr) %>%
    summarise(mean_D = mean((cohens_d)))
  
  t_test_out %>%
    ggplot(aes(
      x=(cohens_d),
      fill=inv
    )) +
    geom_vline(xintercept = 0.8, linetype = "dashed", size = 0.4) +
    geom_vline(xintercept = -0.8, linetype = "dashed", size = 0.4) +
    geom_histogram() +
    ggtitle("Effect size analysis: real vs permutation rnPvs") +
    xlab("Cohen's D") +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(mod~chr)-> 
    cohens_D_dist
  ggsave(cohens_D_dist, file = "cohens_D_dist.pdf", 
         width = 8,
         height = 4)
  
  ## T-test distribution
  t_test_out %>%
    ggplot(aes(
      x=ith,
      y=est,
      ymin=low_c,
      ymax=high_c,
      color=chr
    )) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(width = 0.5) +
    geom_point(size = 0.5) +
    ggtitle(paste("t-test estimates", sep = "")) +
    facet_grid(inv~mod) ->
    t_test_bars
  
  ggsave(t_test_bars, file = "t_test_bars.pdf", 
         width = 8,
         height = 4)

  ## distirbution of p-values
  t_test_out %>%
    ggplot(aes(
      x=ith,
      y=-log10(p.adjust(p.val, "bonferroni")),
      color=chr,
      shape=inv
    )) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gold") +
    ggtitle(paste("t-test p values", sep = "")) +
    facet_grid(chr~mod) ->
    t_test_pvals
  
  ggsave(t_test_pvals, file = "t_test_pvals.pdf", 
         width = 8,
         height = 4)
  