o.ag <- rbindlist(o.ag)

o.ag.ag <- o.ag[,list(N=sum(N)), list(perm, mod, var, chr, inv, cluster)]
o2 <- o.ag.ag[,list(N, prop=N/sum(N), totalN=sum(N), mod, var),
                    list(perm , chr, inv, cluster)]

o2.ag <- o2[,list(prop.real=prop[perm==0], totalN=totalN[perm==0],
                  prop.perm.mu=mean(prop[perm!=0]),
                  prop.perm.lci=quantile(prop[perm!=0], .01),
                  prop.perm.uci=quantile(prop[perm!=0], .99),
                  prop.perm.med=median(prop[perm!=0]),
                  prop.rr=mean(log2(prop[perm==0]/prop[perm!=0])),
                  prop.sd=sd(log2(prop[perm==0]/prop[perm!=0]))),
              list(chr, inv, mod, var, cluster)]
o2.ag[,rr:=prop.real/prop.perm.mu]

o2.ag[order(-prop.real)][]
o2.ag[prop.real>(prop.perm.uci)][order(rr)]
o2.ag[prop.rr-2*prop.sd>0][cluster=="2.North_America_E"][order(-prop.real)][]
o2.ag[,p:=pnorm(0, prop.rr, prop.sd)]


save(o.ag, o.ag.ag, file="~/environmentalOmnibus_global.summary.Rdata")
