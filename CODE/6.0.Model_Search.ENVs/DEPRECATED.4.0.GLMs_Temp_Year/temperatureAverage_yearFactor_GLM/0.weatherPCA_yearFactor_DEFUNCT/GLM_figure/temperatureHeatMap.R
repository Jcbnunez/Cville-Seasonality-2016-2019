### libraries
  library(data.table)
  library(ggplot2)
  library(viridis)
  library(lubridate)
  library(foreach)
  library(patchwork)

### load samps
  samps <- fread("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/DEST_10Mar2021_POP_metadata.csv")
  samps <- samps[set!="dgn"]

  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]
  samps[,Date:=date(paste(year, month, day, sep="-"))]

  samps[locality=="UA_od", locality:="UA_Ode"]

  targetLocales <- list("FI_Aka"=list(locality="FI_Aka",   minYear=2000),
                        "DE_Bro"=list(locality="DE_Bro",   minYear=2000),
                        "VA_ch"=list(locality="VA_ch",     minYear=2014),
                        "PA_li"=list(locality="PA_li",     minYear=2000),
                        "DE_Mun"=list(locality="DE_Mun",   minYear=2000),
                        "UA_Ode"=list(locality=c("UA_Ode"), minYear=2010),
                        "TR_Yes"=list(locality="TR_Yes",   minYear=2000))
  samps[,id:=sampleId]

### load average weather data
  load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/weatherAve.Rdata")
  setnames(weather.ave, "V1", "sampleId")

  weather.ave <- merge(weather.ave, samps, by="sampleId")
  setnames(weather.ave, "locality.x", "locality")
  weather.ave[,ord:=as.factor(paste(year, yday, "_"))]

### load raw weather data
  load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/weather1.Rdata")

### add in delta

	weather30 <- foreach(i=1:length(targetLocales), .combine="rbind")%do%{
		#i <- j <- 4
		samps.tmp <- samps[locality==targetLocales[[i]]$locality & year>=targetLocales[[i]]$minYear]

		foreach(j=1:dim(samps.tmp)[1], .combine="rbind")%do%{

			weather.tmp <- weather1[J(data.table(locality=samps.tmp[j]$locality))]
			weather.tmp[,delta:=date-samps.tmp[j]$Date]
			data.table(sampleId=samps.tmp[j]$sampleId, year=samps.tmp[j]$year, yday=samps.tmp[j]$yday,
                 aveTemp=weather.tmp[delta>= -30 & delta<=0]$tempAvg,
                 delta=weather.tmp[delta>= -30 & delta<=0]$delta,
                 locality=weather.tmp[1]$locality)
		}
	}
  weather30[,ord:=as.factor(paste(year, yday, "_"))]
  ]
### delta_plot
  delta_plot <- ggplot(data=weather30,
        aes(x=delta*-1,
            y=ord,
            fill=aveTemp/10)) +
  geom_tile() +
  facet_grid(.~locality, scales="free_x", space="free_x") +
  scale_fill_viridis(limits=c(-5, 35), option="B", name="Temp °C") +
  #scale_fill_gradient2(limits=c(-5, 35), midpoint=15, name="Temp °C", low="steelblue", high="firebrick") +
  xlab("Days prior to collection") + ylab("") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x=element_blank())


### average plot

  avePlot <-
  ggplot(data=weather.ave) +
  geom_line(aes(y=aveTemp/10, x=ord, group=year)) +
  geom_point(aes(y=aveTemp/10, x=ord, color=aveTemp/10), size=2) +
  facet_grid(.~locality, scales="free_x", space="free_x") +
  scale_color_viridis(limits=c(-5, 35), option="B", name="Temp °C") +
  #scale_color_gradient2(limits=c(-5, 35), midpoint=15, name="Temp °C", low="steelblue", high="firebrick") +
  ylab("Average Temperature") + xlab("") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, size=8)) +
  scale_x_discrete(breaks=weather.ave$ord, labels=weather.ave$sampleId)

#scale_fill_gradient2(min="steelblue", max = "firebrick")
### mega plot


layout <- "
AAAAAA
AAAAAA
BBBBBB
"

mega <-
delta_plot + avePlot +
plot_layout(design=layout, guides = 'collect')

ggsave(mega, height=7, w=12,
      file="/Users/alanbergland/Documents/GitHub/Overwintering_18_19/GLM_figure/mega_environment.pdf")


ggsave(mega, height=7, w=12,
      file="/Users/alanbergland/Documents/GitHub/Overwintering_18_19/GLM_figure/mega_environment.png")
