scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/environmental_ombibus_global_permV2/bestAIC/bestAIC.v2.Rdata ~/.

### libraries
  library(ggplot2)
  library(data.table)

### load
  load("~/bestAIC.v2.Rdata")
  o2.ag <- o2.ag[!cluster%in%c("2.North_America_I95", "2.North_America_Mid", "2.North_America_W")]

### remake modRank to order by old-2L-inv

  focal_order <- o2.ag[perm_strategy=="old"][chr=="2L"][inv=="Inside Inversion"][cluster=="5.Cville"][,c("var", "mod", "modRank"), with=F]
  setnames(focal_order, "modRank", "modRank_focal")
  setkey(focal_order, var, mod)
  setkey(o2.ag, var, mod)
  o2.ag <- merge(o2.ag, focal_order)

### plot
  o2.ag[,sig:=log2(prop.real/prop.perm.uci)>0 | log2(prop.real/prop.perm.lci)<0]
  o2.ag[var=="null", label:="null"]
  o2.ag[var=="pop_year", label:="pop_year"]
  o2.ag[is.na(label), label:="environment"]

  o2.ag[,perm_strategy:=factor(perm_strategy, levels=c("old", "new"))]
  o2.ag[,cluster:=factor(cluster, levels=c("5.Cville", "1.Europe_W", "3.Europe_E", "2.North_America_E"))]

  o2.ag[is.na(stage), stage:="stage0"]

  o2.ag.use <- o2.ag[stage=="stage0" | stage=="stage1" | (stage=="stage2" & !var%in%c("null", "pop_year"))]


  p1 <- ggplot(data=o2.ag.use[sig==T], aes(x=modRank_focal, y=log2(rr), color=perm_strategy, shape=label)) +
  facet_grid(cluster~chr+inv+perm_strategy) +
  geom_point(data=o2.ag.use[sig==F], aes(x=modRank_focal, y=log2(rr)), color="grey", alpha=.5) +
  geom_point() +
  geom_segment(data=o2.ag.use[var=="temp.max" & mod==2][stage=="stage2"],
              aes(x=modRank_focal+5, xend=modRank_focal, yend=log2(rr), y=-5), arrow = arrow(length = unit(0.1, "cm")), color="black")


  ggsave(p1, file="~/new_old.pdf", width=18, height=7)


### which environmentral variables behave nicely across the board?
  o3 <- o2.ag.use[,list(rr.mu=mean(log2(rr)), rr.sd=sd(log2(rr)),
                        nSig_pos=mean(sig==T & rr>1),
                        cville=any(sig[cluster=="5.Cville"]==T)),
                   list(perm_strategy, mod, var, stage)]
o3[mod>=1,group:=tstrsplit(var, "\\.")[[1]]]
o3[mod<1, group:="base"]

ggplot(data=o3, aes(x=nSig_pos, y=interaction(mod, var), color=cville)) +
geom_point() + facet_grid(group~perm_strategy, scales="free_y", space="free_y")

  ggplot(o2.ag.use, aes(x=interaction(mod, var), y=log2(rr))) + geom_boxplot() + facet_grid(~perm_strategy)

    ggsave(all_combos.plot, file="/Users/alanbergland/Documents/GitHub/Overwintering_18_19/EnvironmentalOmnibus_Global/redoPerm/plot_scripts/new_old_rank.pdf", width=12, height=7)

  inv2L.plot <- ggplot(data=o2.ag[chr=="2L"][inv=="Inside Inversion"][cluster=="5.Cville"],
                aes(x=modRank_focal, y=log2(rr), color=perm_strategy, group=interaction(mod, var))) +
                geom_point() + facet_grid(cluster~chr+inv)


### pairwise
  o2.wide <- dcast(o2.ag, mod+var+chr+inv+cluster~perm_strategy, value.var=c("rr", "modRank_focal", "modRank", "prop.real", "prop.perm.mu"))

  lim <- 5
  ggplot(data=o2.wide, aes(x=log2(rr_new), y=log2(rr_old))) +
  geom_point() +
  geom_point(data=o2.wide[var=="null"], aes(x=    log2(rr_new), y=log2(rr_old)), color="red") +
  geom_point(data=o2.wide[var=="pop_year"], aes(x=log2(rr_new), y=log2(rr_old)), color="blue") +
  ylim(-lim, lim) + xlim(-lim, lim) + facet_grid(cluster~chr+inv)



  lim <- 5
  ggplot(data=o2.wide, aes(x=log10(prop.real_new), y=log10(prop.real_old))) +
  geom_point() +
  geom_point(data=o2.wide[var=="null"], aes(x=    log10(prop.real_new), y=log10(prop.real_old)), color="red") +
  geom_point(data=o2.wide[var=="pop_year"], aes(x=log10(prop.real_new), y=log10(prop.real_old)), color="blue") +
  facet_grid(cluster~chr+inv)



  ggplot(data=o2.wide, aes(x=log10(prop.perm.mu_new), y=log10(prop.perm.mu_old))) +
  geom_point() +
  geom_point(data=o2.wide[var=="null"], aes(x=    log10(prop.perm.mu_new), y=log10(prop.perm.mu_old)), color="red") +
  geom_point(data=o2.wide[var=="pop_year"], aes(x=log10(prop.perm.mu_new), y=log10(prop.perm.mu_old)), color="blue") +
  facet_grid(cluster~chr+inv)



  ggplot(data=o2.wide, aes(x=log10(prop.real_new), y=log10(prop.perm.mu_new))) +
  geom_point() +
  geom_point(data=o2.wide[var=="null"], aes(x=    log10(prop.real_new), y=log10(prop.perm.mu_new)), color="red") +
  geom_point(data=o2.wide[var=="pop_year"], aes(x=log10(prop.real_new), y=log10(prop.perm.mu_new)), color="green") +
  facet_grid(cluster~chr+inv) +
  ylim(-6, 0) + xlim(-6, 0) + geom_abline(slope=1, intercept=0)



### how


### how many time is each mod-var significantly positive?
  o2.ag[,list(nSigPos=sum(sig==T & rr>1), nSigNeg=sum(sig==T & rr<1)), list(mod, var, perm_strategy)][order(nSigPos)]
  o2.ag[var=="temp.max"][mod==4]






        ggplot(data=o2.ag[chr=="2L"][inv=="Inside Inversion"][cluster=="5.Cville"],
                      aes(x=modRank_focal, y=prop.real, group=interaction(mod, var))) +
        geom_point() +
        geom_point(data=o2.ag[chr=="2L"][inv=="Inside Inversion"][cluster=="5.Cville"],
                      aes(x=modRank_focal, y=prop.perm.mu, group=interaction(mod, var)))+
        facet_grid(perm_strategy~chr+inv)



        ggplot(data=o2.ag[chr=="2L"][inv=="Inside Inversion"][cluster=="1.Europe_W"],
                      aes(x=qlogis(prop.perm.mu), y=qlogis(prop.real), group=interaction(mod, var), color=log2(rr))) +
        geom_point() + ylim(-12, -1) + xlim(-12, -1) + geom_abline(slope=1,intercept=0) + facet_grid(~perm_strategy)


o2.ag[cluster=="5.Cville"][perm_strategy=="new"][prop.real<.2][which.max(rr)]

log2(


  o2.wide[modRank_old==112]

  o2.wide <- dcast(o2.ag, mod+var+chr+inv+cluster+sig~perm_strategy, value.var="rr")
  ggplot(data=o2.wide, aes(x=new, y=old, color=sig)) + geom_point() + facet_grid(cluster~chr+inv) + ylim(0,40) + xlim(0,40) + geom_abline(intercept=0,slop=11)


  o2.wide <- dcast(o2.ag, mod+var+chr+inv+perm_strategy~cluster, value.var="sig")
  setnames(o2.wide, c("1.Europe_W", "3.Europe_E", "5.Cville"), c("Europe_W", "Europe_E", "Cville"))
  fisher.test(table(o2.wide[perm_strategy=="old"][chr=="2L"][inv=="Inside Inversion"]$"1.Europe_W",
                    o2.wide[perm_strategy=="old"][chr=="2L"][inv=="Inside Inversion"]$"3.Europe_E"))


o2.wide[perm_strategy=="new"][Europe_W==T][Europe_E==T][Cville==T]

  o2.ag[cluster=="3.Europe_E"][mod==-1][chr=="2L"]
