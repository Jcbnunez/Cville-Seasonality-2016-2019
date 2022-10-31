### librarues
  library(data.table)
  library(ggplot2)
  library(patchwork)

### tranlation table

  a <- c(
    "AGA_14", "GA_at", 2014,
    "BA_12",  "ES_ba", 2012,
    "BHM_14", "MI_bh", 2014,
    "CUA_14", "VA_ch", 2014,
    "CUA_15", "VA_ch",2015,
    "CWI_14", "WI_cp",2014,
    "LMA_14", "MA_la",2014,
    "MA_12", "MA_la",2012,
    "NY_12", "NA_it",2012,
    "PA_12", "PA_li",2012,
    "PA_14", "PA_li",2014,
    "PA_9",  "PA_li",2009,
    "SCPA_14", "PA_st",2014,
    "SON_15", "ON_su",2015,
    "TKA_14", "KA_to",2014,
    "VI_12", "AT_gr",2012,
    "WI_12", "WI_cp",2012,
    "WI_13", "WI_cp",2013,
    "co_13", "CA_tu",2013,
    "rd_12", "CA_es",2012,
    "rd_13", "CA_es",2013)

### load new temperature data
  load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/Core20_enrichment/DEST_core20/core20_samps.Rdata")

  load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/global_temperature/weatherAve.Rdata")

  m <- merge(samps.use, weather.ag, by="sampleId")

### load weather data used in Machado et al 2020
  load("/Users/alanbergland/Documents/GitHub/dmel_seasonal_RTEC/predictability_model/weather_and_pbs.fix.Rdata")
  load("/Users/alanbergland/Documents/GitHub/dmel_seasonal_RTEC/predictability_model/figure3/Fig3_data/popStats.ag.Rdata")


  tdt <- data.table(pop=a[seq(1, length(a), 3)],
                    locality.x=a[seq(1, length(a), 3)+1],
                    year.x=a[seq(1, length(a), 3)+2])

  setkey(popStats.ag, pop)
  setkey(tdt, pop)
  popStats.ag <- merge(popStats.ag, tdt)

  tlm <- data.table(y=c(popStats.ag$tmax.max, popStats.ag$tmin.min), beta=sign(popStats.ag$beta.obs),
                      season=rep(c("spring", "fall"), each=21), locality.x=popStats.ag$locality.x, year.x=popStats.ag$year.x)

  m[,year.x:=as.numeric(as.character(year.x))]
  tlm[,year.x:=as.numeric(as.character(year.x))]

  setkey(tlm, locality.x, year.x, season)
  setkey(m, locality.x, year.x, season)
  m <- merge(m, tlm)


### plots
  temp.plot <-
  ggplot() +
  geom_point(data=m[beta==1],
              aes(x=season, y=aveTemp, group=interaction(locality.x, year.x))) +
  geom_point(data=m[beta==-1],
              aes(x=season, y=aveTemp, group=interaction(locality.x, year.x), color=interaction(locality.x, year.x))) +
  geom_line(data=m[beta==1],
              aes(x=season, y=aveTemp, group=interaction(locality.x, year.x))) +
  geom_line(data=m[beta==-1],
              aes(x=season, y=aveTemp, group=interaction(locality.x, year.x), color=interaction(locality.x, year.x)))


  tlm.plot <-
  ggplot() +
  geom_point(data=m[beta==1],
              aes(x=season, y=y, group=interaction(locality.x, year.x))) +
  geom_point(data=m[beta==-1],
              aes(x=season, y=y, group=interaction(locality.x, year.x), color=interaction(locality.x, year.x))) +
  geom_line(data=m[beta==1],
              aes(x=season, y=y, group=interaction(locality.x, year.x))) +
  geom_line(data=m[beta==-1],
              aes(x=season, y=y, group=interaction(locality.x, year.x), color=interaction(locality.x, year.x)))


temp.plot + tlm.plot




a[seq(1, length(a), 6)]
