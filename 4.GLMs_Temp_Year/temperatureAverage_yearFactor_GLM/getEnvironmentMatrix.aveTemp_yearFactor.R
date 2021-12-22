
### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
	library(data.table)
	library(gdata)
	library(foreach)
	library(ggplot2)
	library(ggmap)
	library(maps)
	library(mapdata)
	library(rnoaa)
	library(sp)
  library(lubridate)
  library(doMC)
  registerDoMC(20)

### load in samps
	samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
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

### get weather data
	weather <- foreach(i=1:length(targetLocales))%dopar%{
		message(i)
		stations <- meteo_nearby_stations(
				lat_lon_df=as.data.frame(samps[locality==targetLocales[[i]]$locality]),
				lat_colname = "lat",
				lon_colname = "long",
				station_data = ghcnd_stations(),
				var = c("TAVG", "TMIN", "TMAX"),
				year_min = min(samps[locality==targetLocales[[i]]$locality]$year),
				year_max = max(samps[locality==targetLocales[[i]]$locality]$year),
				radius = NULL,
				limit = 3
		)

		stations <- rbindlist(stations)
		stations.ag <- stations[,list(name=name[i], lat=latitude[1], long=longitude[1], dist=distance[1]), id]

		data <- meteo_pull_monitors(monitors=stations$id,
																date_min=min(samps[locality==targetLocales[[i]]$locality]$Date)-60,
																date_max=max(samps[locality==targetLocales[[i]]$locality]$Date),
																var="all")

		data <- as.data.table(data)
		data[,locality:=targetLocales[[i]]$locality]
		data <- merge(data, stations.ag, by="id")
		data
	}

	weather <- rbindlist(weather, fill=T)

	### quick summary and check
		weather[,list(tavg=mean(is.na(tavg)), tmin=mean(is.na(tmin)), tmax=mean(is.na(tmax)), dist=mean(dist)), list(id, locality)]
		weather[,temp:=tmin/2 + tmax/2]
		summary(lm(temp~tavg, weather))

	### simplify and make average temp
		weather <- weather[,c("id", "date", "tavg", "tmax", "tmin", "locality", "dist"), with=F]

		weather[,tempAvg:=tavg]
		weather[is.na(tempAvg), tempAvg:=tmax/2 + tmin/2]

	### rank stations by closeness
		weather.ag <- weather[,list(dist=dist[1], minYear=min(date), maxYear=max(date)), list(locality, id)]
		weather.ag <- weather.ag[,list(rank_dist=rank(dist), dist=dist, id=id, minYear=minYear, maxYear=maxYear), list(locality)]

		### fix PA site
			weather.ag[locality=="PA_li" & rank_dist==1, rank_dist:=4]
			weather.ag[locality=="PA_li" & rank_dist==2, rank_dist:=1]

		# weather <- weather[,-c("rank_dist"), with=F]
		weather <- merge(weather, weather.ag[,c("rank_dist", "id"), with=F], by="id")

		table(is.na(weather$tempAvg), weather$rank_dist)

	### extract out closest station and pad missing data with next closest station
		weather1 <- weather[rank_dist==1]
		setkey(weather, locality, rank_dist, date)

		### pad TR_yes
			weather1[is.na(tempAvg) & locality=="TR_Yes",
							rank_dist:=2]
			weather1[is.na(tempAvg) & locality=="TR_Yes",
							tempAvg:=weather[J(data.table(locality="TR_Yes", rank_dist=2, date=weather1[is.na(tempAvg)][locality=="TR_Yes"]$date))]$tempAvg]

		### pad VA_ch
			weather1[is.na(tempAvg) & locality=="VA_ch",
							rank_dist:=3]
			weather1[is.na(tempAvg) & locality=="VA_ch",
							tempAvg:=weather[J(data.table(locality="VA_ch", rank_dist=3, date=weather1[is.na(tempAvg)][locality=="VA_ch"]$date))]$tempAvg]

		#### sub in PA_li station rank 2 for station rank 1. r1 does not
		#	weather1[locality=="PA_li",
		#					rank_dist:=2]
		#	weather1[locality=="PA_li",
		#					tempAvg:=weather[J(data.table(locality="PA_li", rank_dist=2, date=weather1[locality=="PA_li"]$date))]$tempAvg]
#

		table(is.na(weather1$tempAvg), weather1$locality, weather1$rank_dist)

	### save
		save(weather1, weather, file="~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/weather1.Rdata")

### calculate average temperature 30 days prior to collection date
	setkey(weather1, locality)

	weather.ave <- foreach(i=1:length(targetLocales), .combine="rbind")%do%{
		#i <- j <- 4
		samps.tmp <- samps[locality==targetLocales[[i]]$locality & year>=targetLocales[[i]]$minYear]

		foreach(j=1:dim(samps.tmp)[1], .combine="rbind")%do%{

			weather.tmp <- weather1[J(data.table(locality=samps.tmp[j]$locality))]
			weather.tmp[,delta:=date-samps.tmp[j]$Date]
			data.table(samps.tmp[j]$sampleId, aveTemp=mean(weather.tmp[delta>= -30 & delta<=0]$tempAvg), locality=weather.tmp[1]$locality)
		}
	}

	save(weather.ave, file="~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/weatherAve.Rdata")
