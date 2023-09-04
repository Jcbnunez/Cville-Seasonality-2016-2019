### Reshuffle FST overwintering

library(tidyverse)
library(magrittr)
library(foreach)
library(reshape2)

dat <- get(load("Fig1.panelC.Rdata"))

dat$pop1 %>% unique -> pops.vec

####
#pop="Charlottesville"


shuffle.fst = 
foreach(pop = pops.vec,
.combine = "rbind")%do%{
### loop
message(pop)

dat %>%
filter(pop1 == pop & pop2 == pop) %>%
filter(bin_date != "3.Multi-Year") -> tmp

tmp %>%
filter(bin_date == "1.within") -> with.tmp

tmp %>%
filter(bin_date == "2.Overwinter") -> over.tmp

t.test(over.tmp$FST,with.tmp$FST) -> real.obs

real.t = real.obs$statistic
real.p =real.obs$p.value

###
data.frame(
pop,
t=real.t,
p=real.p,
perm = 0
) -> o.real

### permutation
perm.o = 
foreach(p = 1:100,
.combine = "rbind")%do%{

n =  dim(tmp)[1]
tmp.shuff = tmp
tmp.shuff$bin_date = 
tmp$bin_date[sample(1:n, n, replace = FALSE)]

tmp.shuff %>%
filter(bin_date == "1.within") -> with.shuff

tmp.shuff %>%
filter(bin_date == "2.Overwinter") -> over.shuff

t.test(with.shuff$FST,over.shuff$FST) -> perm.obs

perm.t = perm.obs$statistic
perm.p = perm.obs$p.value

data.frame(
pop,
t=perm.t,
p=perm.p,
perm = p
) -> o.perm
return(o.perm)
} # close inner

rbind(
o.real,
perm.o
)

} # close outer

### loop
shuffle.fst %>%
group_by(type = perm == 0, pop) %>%
summarize(uci = quantile(t, 0.95)) %>%
dcast(pop~type, value.var = "uci") %>%
mutate(signifi = case_when(
TRUE > FALSE ~ "sig",
TRUE < FALSE ~ "ns"
))


#### plot

ggplot() +
geom_boxplot(
data = filter(shuffle.fst, perm != 0),
aes(
x=pop,
y=t), outlier.shape = NA, width = 0.6) +
geom_point(
data = filter(shuffle.fst, perm == 0),
aes(
x=pop,
y=t),
size = 4,
shape = 23,
fill = "red",
color = "black") +
theme_bw() +
xlab("Population") +
ylab("t statistic") +
coord_flip() ->
fst.perms

ggsave(fst.perms, file = "fst.perms.pdf", 
h = 4, w = 4)




