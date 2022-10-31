### libraries
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(viridis)

### load samps data
  samps <- fread("~/DEST_10Mar2021_POP_metadata.csv")

  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]

  setkey(samps, sampleId)
  samps <- samps[J(unique(wm.impute$sampleId))]

### load weather data
  load(file="~/weather_impute.Rdata")
  wm.impute[,term:=paste(variable, dateDelta, sep="_")]

  ww.temps <-dcast(wm.impute[J(samps[locality=="VA_ch"][year>2012]$sampleId)][variable=="tmax"],
                sampleId~term, value.var="value")

  ww.mat <- as.matrix(ww.temps[,-"sampleId", with=F])
  row.names(ww.mat) <- as.matrix(ww.temps[,"sampleId", with=F])[,1]

  ww.temps.pca <- princomp(ww.mat)
  temps.pca <- data.table(sampleId=row.names(ww.temps.pca$scores),
                          pc1=ww.temps.pca$scores[,1])

  temps.pca <- merge(temps.pca, samps, by="sampleId")

  temps.summary <- wm.impute[J(samps[locality=="VA_ch"][year>2012]$sampleId)][variable=="tmax"][,list(mu=mean(value)), list(sampleId)]

  temps.pca <- merge(temps.pca, temps.summary, by="sampleId")


### plots
  wm.impute[,degreesC:=value/10]

  temp.plot <- ggplot(data=wm.impute[J(samps[locality=="VA_ch"][year>2012]$sampleId)],
          aes(x=delta, y=gsub("VA_ch_", "", sampleId), fill=degreesC)) +
          geom_tile() + theme_bw() +
          theme(axis.text.y=element_text(size=6)) +
          facet_grid(~variable) +
          xlab("Days before collection") +
          ylab("") +
          scale_fill_viridis(option="inferno") +
          labs(fill="Max Temp °C")


  pca.yday.plot <- ggplot(data=temps.pca, aes(x=yday, y=pc1, color=as.factor(year))) +
  geom_point() + theme_bw() +
  ylab("Principal\nComponent 1") +
  xlab("Julian Day") + labs(color="Year")

  pca.avaTemp.plot <- ggplot(data=temps.pca, aes(x=mu/10, y=pc1, color=as.factor(year))) +
  geom_point() + theme_bw() +
  ylab("Principal\nComponent 1") +
  xlab("Average Max Temperature, °C")+ labs(color="Year")

  layout <-"
  AAABB
  AAACC
  "

mega_weather <-
  temp.plot + pca.yday.plot + pca.avaTemp.plot +
  plot_layout(design =layout) +
  plot_annotation(tag_levels="A")


ggsave(mega_weather, file="~/mega_weather.png", h=4, w=8)
