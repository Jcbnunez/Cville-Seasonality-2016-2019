### plot enrichment analyses
### 

### libraries
library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(4)
library(tidyr)
library(patchwork)

###
files_chr = system("ls ./out.enr/*.Rdata | grep CHR", intern = T)
files_win = system("ls ./out.enr/*.Rdata | grep win", intern = T)
files_Cvile = system("ls GLM_omnibus_window_analysis/*.Rdata | grep Cville", intern = T)


enrich_chr = 
foreach(i = 1:length(files_chr), .combine = "rbind")%do%{
  
  load(files_chr[i])
  return(o)
  
}

enrich_win = 
  foreach(i = 1:length(files_win), .combine = "rbind")%do%{
    
    message(i)
    load(files_win[i])
    return(o.win)
    
  }

cville_win = get(load(files_Cvile)) 

cville_win %>%
 mutate(data_type = case_when(perm == 0 ~ "real",
                              perm != 0 ~ "perm")) %>%
 group_by(data_type, pos_mean, chr) %>%
 summarise(uci.rnvp = quantile(rnp.binom.p, 0.005),
           uci.wZa = quantile(rnp.wZa.p, 0.005),
           ) -> cville_win.summarized

cville_win.summarized %>%
filter(uci.wZa < 1e-70 & data_type == "real" & chr == "2L" & pos_mean > 2e6)


#######
#######

cville_win.summarized %>%
  filter(chr == "2L") %>%
  ggplot(aes(
    x=pos_mean/1e6,
    #y=-log10(uci),
    color=data_type
  )) +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) +
  geom_vline(xintercept = 5.1-0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 9.6-0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 5.1+0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 9.6+0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 6.2-0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 6.2+0.2, color = "blue", linetype = "dashed") +
  geom_line(aes(y=-log10(uci.rnvp))) +
  geom_line(aes(y=log10(uci.wZa))) ->
  cville_plot

#enrich_win[p < 0.1][chr.x == "2L"][anchor.model == "humidity.ave;8;1.Europe_W"][win.start > 2e6]

#enrich_win %>%
#  dplyr::select(chr.x, win.start, win.end,p, anchor.model ) %>%
#  dcast(chr.x+ win.start+ win.end~anchor.model, value.var = "p") %>%
#  filter(chr.x == "2L") %>%
#  filter(`humidity.ave;8;1.Europe_W` < 0.01 & `humidity.ave;8;1.Europe_W` < 0.01)

### use perm to inf ucis
enrich_win %>%
  mutate(data_type = case_when( 
    perm == 0 ~ "real",
    perm != 0 ~ "perm"
    )) %>%
  group_by(chr.x, win.start, win.end, data_type, anchor.model) %>%
  summarise(uci.p = quantile(p, 0.05)) %>%
  dcast(anchor.model + chr.x+win.start+win.end~data_type, value.var = "uci.p") %>% 
  .[complete.cases(.),] %>%
  mutate(sig.vis.perm = case_when(real < perm ~ "sig",
                                  real >= perm ~ "not.sig")) ->
  enrich_win.processed

enrich_win %>%
  filter(perm == 0) %>%
  right_join( dplyr::select(enrich_win.processed, anchor.model, chr.x, win.start, win.end,sig.vis.perm) ) %>% 
  filter(chr.x == "2L") %>%
  ggplot(aes(
    x=((win.start+win.end)/2)/1e6,
    #y=-log10(p),
    y=log2(or),
    color=sig.vis.perm,
  )) +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) +
  geom_vline(xintercept = 5.1-0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 5.1+0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 9.6-0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 9.6+0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 6.2-0.2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 6.2+0.2, color = "blue", linetype = "dashed") +
  geom_point(size = 1) +
  geom_hline(yintercept = -log2(1), linetype = "dashed") +
  facet_grid(anchor.model~.) ->
  enrich_plot_plot

ggsave(cville_plot/enrich_plot_plot, file = "tester.enrh.pdf", w =9)

#######

enrich_chr %>%
  ggplot() +
  geom_violin(
    data= filter(enrich_chr, perm != 0),
    aes(x=chr.x,
        y=-log10(p),
        fill=inv,
        ), position=position_dodge(width=0.5)
    ) +
  geom_point(
    data= filter(enrich_chr, perm == 0),
    aes(x=chr.x,
        y=-log10(p),
        color=inv
        ), position=position_dodge(width=0.5)
  ) +
  facet_grid(anchor.model~.) -> p.vis.perm

ggsave(p.vis.perm, file = "p.vis.perm.pdf", w = 5, h = 3)


enrich_chr %>%  
  filter(perm == 0) %>%
  ggplot(aes(
    x=chr.x,
    y=log2(or),
    ymin=log2(lci),
    ymax=log2(uci),
    color=inv
  )) +
  geom_point(size = 2, position=position_dodge(width=0.5)) +
  geom_errorbar(
    width = 0.1,
    position=position_dodge(width=0.5))  +
  geom_hline(yintercept = 0) +
  facet_grid(anchor.model~.) ->
  chr_enric.plot

ggsave(chr_enric.plot, file = "chr_enric.plot.pdf", w = 5, h = 3)

