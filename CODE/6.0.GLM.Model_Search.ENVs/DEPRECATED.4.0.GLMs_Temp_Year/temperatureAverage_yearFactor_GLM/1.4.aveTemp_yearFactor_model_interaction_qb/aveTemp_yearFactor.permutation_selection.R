### libraries
  library(data.table)
  library(gdata)
  library(lubridate)
  library(foreach)


### weather data
  load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/weatherAve.Rdata")
  setnames(weather.ave, "V1", "sampleId")

### samps
  #samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
  samps <- fread("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/DEST_10Mar2021_POP_metadata.csv")
  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]
  samps <- samps[set!="dgn"]
  samps[,Date:=date(paste(year, month, day, sep="-"))]

  samps[locality=="UA_od", locality:="UA_Ode"]

  samps <- merge(samps, weather.ave[,c("aveTemp", "sampleId")], by="sampleId")


  targetLocales <- list("FI_Aka"=list(locality="FI_Aka",   minYear=2000),
                        "DE_Bro"=list(locality="DE_Bro",   minYear=2000),
                        "VA_ch"=list(locality="VA_ch",     minYear=2014),
                        "PA_li"=list(locality="PA_li",     minYear=2000),
                        "DE_Mun"=list(locality="DE_Mun",   minYear=2000),
                        "UA_Ode"=list(locality=c("UA_Ode"), minYear=2010),
                        "TR_Yes"=list(locality="TR_Yes",   minYear=2000))

### iterate through localities
  perms <- foreach(i=1:5000, .combine="rbind")%do%{
    set.seed(i)
    ord.tmp <- samps[,list(aveTemp=aveTemp,
                         year=year,
                         sampleId=sample(sampleId),
                         sampleId.old=sampleId),
                    list(locality)]

    ord.tmp <- merge(ord.tmp, samps[,c("sampleId", "aveTemp", "year"), with=F], by.x="sampleId", by.y="sampleId")
    temp.cor <- cor.test(ord.tmp$aveTemp.x, ord.tmp$aveTemp.y)
    year.cor <- cor.test(ord.tmp$year.x, ord.tmp$year.y)

    ord.tmp[,list(cor_temp=cor(aveTemp.x, aveTemp.y),
                  cor_year=cor(year.x, year.y),
                  i=i),
             list(locality)]
  }

  perms <- perms[,list(i=i, rank_temp=rank(abs(cor_temp), ties="random"), rank_year=rank(abs(cor_year), ties="random"),
                      cor_temp, cor_year),
                  list(locality)]


perms[rank_temp<300 & rank_year<300][,list(i=sample(as.character(i), 100)), locality]


  perml <- melt(perms, measure.vars=c("cor_temp", "cor_year"), id.vars=c("i", "locality"))
  perml <- perml[,list(i=i, value=value, rank=rank(abs(value)-0, ties="random")), list(locality, variable)]

  ggplot(data=perml) +
  geom_histogram(aes(value, group=as.factor(rank<=100), fill=as.factor(rank<=100))) +
  facet_grid(variable~locality)

  perml[locality=="VA_ch"][rank<=100][variable=="cor_temp"]
