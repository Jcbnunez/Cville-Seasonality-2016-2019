
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

### load in new Samps
	samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
	samps <- samps[set!="dgn"]

  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]
	samps[,Date:=date(paste(year, month, day, sep="-"))]

### load in MBE samps which has the "DateExact" column
  popInfo <- fread("/scratch/aob2x/dest/DEST_freeze1/populationInfo/samps_10Nov2020.csv")
	samps <- merge(samps, popInfo[,c("sampleId", "DateExact"), with=F], all.x=T, by="sampleId")
	samps[set=="CvilleSet", DateExact:=T]
	samps[is.na(DateExact), DateExact:=F]

	table(samps$set, samps$DateExact)
	samps <- samps[DateExact==T]

### get stations
	samps[,id:=sampleId]
	stations <- meteo_nearby_stations(
			lat_lon_df=as.data.frame(samps),
			lat_colname = "lat",
			lon_colname = "long",
			station_data = ghcnd_stations(),
			var = c("TMAX","TMIN"),
			year_min = min(samps$year),
			year_max = max(samps$year),
			radius = NULL,
			limit = 1
	)
	stations <- rbindlist(stations)
	stations <- as.data.table(stations)
	stations[,sampleId:=samps$sampleId]
	samps <- merge(samps, stations, by="sampleId")

	samps[locality=="UA_Ode",c("sampleId", "name", "lat", "latitude", "long", "longitude"), with=F]

### get weather data
  getWeatherData <- function(station, collectionDate, range, sampleId) {
    # station<- popInfo$stationId[1]; collectionDate=date(popInfo$collectionDate[1]); range=30; sampleId=popInfo$sampleId[1]

    tmp <- ghcnd_search(
      stationid=station,
      date_min = collectionDate-range,
      date_max = collectionDate,
      var = c("TMAX","TMIN"),
      refresh = F)

    tmax <- as.data.table(tmp$tmax)
    tmin <- as.data.table(tmp$tmin)

    setkey(tmax, id, date)
    setkey(tmin, id, date)

    w <- merge(tmax[,c("id", "tmax", "date"), with=F], tmin[,c("id", "tmin", "date"), with=F])

    w[,dateDelta:=paste("delta", as.numeric(collectionDate-date(date)), sep="")]
		w[,delta:=as.numeric(collectionDate-date(date))]
    w[,sampleId:=sampleId]
    w
    #dcast(w, sampleId ~ dateDelta, value.var=c("tmax", "tmin"))

  }

###
  weatherOut <- foreach(i=samps[locality=="UA_Ode"]$sampleId, .errorhandling="remove")%dopar%{
			print(paste(which(i==samps$sampleId), length(samps$sampleId), sep=" / "))
      #i<-samps$sampleId[2]
      station=samps[sampleId==i]$id.y
      collectionDate=samps[sampleId==i]$Date
      range=30
      sampleId=i
      message(sampleId)
      o <- getWeatherData(station, collectionDate, range, sampleId)
			o
    }
  weatherOut <- rbindlist(weatherOut, fill=T)
	weatherOut <- melt(weatherOut, id.vars=c("id", "date", "dateDelta", "delta", "sampleId"), measure.vars=c("tmax", "tmin"))
	setkey(weatherOut, sampleId, delta, variable)

### imputation
	wm.impute <- foreach(pop=unique(weatherOut$sampleId), .combine="rbind", .errorhandling="remove")%dopar%{
			message(pop)
			foreach(i=c("tmax", "tmin"), .combine="rbind", .errorhandling="remove")%do%{
				### pop=unique(weatherOut$sampleId)[155]; i="tmax"
				tmp <- weatherOut[sampleId==pop][variable==i]

				sp <- split(tmp[is.na(value)]$delta, cumsum(c(1, diff(tmp[is.na(value)]$delta) != 1)))

				if(sum(is.na(tmp$value))>0) {
					foreach(j=1:length(sp))%do%{
						mu <- mean(tmp[(first(sp[[j]]) - 0):(last(sp[[j]]) + 2)]$value, na.rm=T)
						tmp[sp[[j]]+1, value:=round(mu)]
						tmp[sp[[j]]+1, impute:=T]

					}
					tmp[is.na(impute), impute:=F]
				} else {
					tmp[, impute:=F]
				}
				tmp
			}
		}
		wm.impute.ag <- wm.impute[,list(n=sum(impute)), list(sampleId, variable)]
		table(wm.impute.ag$n)

		setkey(wm.impute, sampleId)
		wm.impute <- wm.impute[J(unique(wm.impute.ag[n<=5]$sampleId))]

		save(wm.impute, file="/project/berglandlab/alan/weather_impute.Rdata")
