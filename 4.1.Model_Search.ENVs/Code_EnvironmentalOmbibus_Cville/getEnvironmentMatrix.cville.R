# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### libraries
	library(data.table)
	library(gdata)
	library(foreach)
	library(rnoaa)
	library(sp)
  library(lubridate)
  library(doMC)
  registerDoMC(5)

### load in samps
	samps <- fread("Overwintering_18_19/temperatureAverage_yearFactor_GLM/DEST_10Mar2021_POP_metadata.csv")
	samps <- samps[set!="dgn"]

  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]
	samps[,Date:=date(paste(year, month, day, sep="-"))]

	samps[locality=="UA_od", locality:="UA_Ode"]

  targetLocales <- list("VA_ch"=list(locality="VA_ch",     minYear=2014))
	samps[,id:=sampleId]

### get weather data

	tmp <- nearest_stations_nooa(country = "UNITED STATES",
	  date = Sys.Date(),
	  add_map = F,
	  point = c(samps[locality=="VA_ch"]$long[1], samps[locality=="VA_ch"]$lat[1]),
	  no_of_stations = 2
	)


	lcd.dt <- foreach(year.i=2012:2018, .combine="rbind")%do%{
		message(year.i)
		lcd.tmp <- as.data.table(lcd(station="72401693736", year=year.i))
		lcd.dt <- lcd.tmp[,c("date", "hourlydrybulbtemperature", "hourlyrelativehumidity", "hourlyprecipitation", "hourlywindspeed"), with=F]
		lcd.dt[,date:=ymd_hms(gsub("T", " ", date))]
		lcd.dt[,temp:=as.numeric(as.character(hourlydrybulbtemperature))]
		lcd.dt[,humidity:=as.numeric(as.character(hourlyrelativehumidity))]
		lcd.dt[,precip:=as.numeric(as.character(hourlyprecipitation))]
		lcd.dt[,wind:=as.numeric(as.character(hourlywindspeed))]
		lcd.dt
	}


	sets <- data.table(mod=c(1:11),
										 start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
										 end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))
	i <- 1

	f2c <- function(x) (x - 32)/(9/5)

	weather.ave <- foreach(j=sets, .combine="rbind")%do%{
		samps.tmp <- samps[locality==targetLocales[[i]]$locality & year>=targetLocales[[i]]$minYear]

		foreach(j=1:dim(samps.tmp)[1], .combine="rbind")%do%{
			lcd.tmp <- lcd.dt
			lcd.tmp[,delta:=julian(date)-julian(samps.tmp[j]$Date)]

			lcd.tmp.ag <- lcd.tmp[,list(dailyMax=max(temp, na.rm=T), dailyMin=min(temp, na.rm=T)), list(delta=round(delta))]

			foreach(k=1:dim(sets)[1], .combine="rbind")%do%{
				lcd.mod <- lcd.tmp[delta>= -1*(sets[k]$end) & delta<=-1*(sets[k]$start)]
				lcd.mod.ag <- lcd.mod[,list(dailyMax=max(f2c(temp), na.rm=T), dailyMin=min(f2c(temp), na.rm=T)), list(delta=round(delta))]

				data.table(samps.tmp[j]$sampleId, locality=samps.tmp[1]$locality, mod=k,
									temp.ave=mean(lcd.mod$temp, na.rm=T),
									temp.var=var(lcd.mod$temp, na.rm=T),
									temp.max=max(lcd.mod$temp, na.rm=T),
									temp.min=min(lcd.mod$temp, na.rm=T),
									temp.propMax=mean(lcd.mod.ag$dailyMax>32),
									temp.propmin=mean(lcd.mod.ag$dailyMin<5),

									humidity.ave=mean(lcd.mod$humidity, na.rm=T),
									humidity.var=var(lcd.mod$humidity, na.rm=T),

									precip.ave=mean(lcd.mod$precip, na.rm=T),
									precip.var=var(lcd.mod$precip, na.rm=T),

									wind.ave=mean(lcd.mod$wind, na.rm=T),
									wind.var=var(lcd.mod$wind, na.rm=T))
			}

		}

	}

save(weather.ave, file="/Users/alanbergland/Documents/GitHub/Overwintering_18_19/EnvironmentalOmbibus_Cville/weather.Rdata")

p1 <- ggpairs(weather.ave[mod==4,c(4:15)],
    upper = list(continuous = wrap("cor", size = 2)),
		lower = list(continuous = wrap("points", alpha = 0.3,    size=0.5))) +
		theme(strip.text = element_text(size = 5), axis.text=element_text(size=4))

ggsave(p1, file="~/modWeather.png")






lcd.dt.ag <- lcd.dt[,list(temp.mu=mean(temp, na.rm=T),
													temp.var=var(temp, na.rm=T)/mean(temp, na.rm=T),
													temp.max=max(temp, na.rm=T),
													temp.min=min(temp, na.rm=T),
													hum.mu=mean(humidity, na.rm=T),
													hum.var=var(humidity, na.rm=T)/mean(humidity, na.rm=T),
													precip.mu=mean(precip, na.rm=T),
													precip.var=var(precip, na.rm=T)/mean(precip, na.rm=T),
													wind.mu=mean(wind, na.rm=T),
													wind.var=var(wind, na.rm=T)/mean(wind, na.rm=T)),
										list(date=date(date))]
lcd.dt.ag[,yday:=yday(date)]
lcd.dt.ag[,year:=year(date)]

















ggplot(data=lcd.dt.ag, aes(x=temp.var, y=temp.delta)) + geom_point()


p1 <-
ggplot(data=lcd.dt.ag, aes(x=yday, y=temp.mu, group=year, color=as.factor(year))) +
geom_line(size=.5, alpha=.5) +
geom_smooth()

p2 <-
ggplot(data=lcd.dt.ag, aes(x=yday, y=temp.var, group=year, color=as.factor(year))) +
geom_line(size=.5, alpha=.5)+
geom_smooth()

p3 <-
ggplot(data=lcd.dt.ag, aes(x=yday, y=hum.mu, group=year, color=as.factor(year))) +
geom_line(size=.5, alpha=.5) +
geom_smooth()

p4 <-
ggplot(data=lcd.dt.ag, aes(x=yday, y=hum.var, group=year, color=as.factor(year))) +
geom_line(size=.5, alpha=.5)+
geom_smooth()


p5 <-
ggplot(data=lcd.dt.ag, aes(x=yday, y=precip.mu, group=year, color=as.factor(year))) +
geom_line(size=.5, alpha=.5) +
geom_smooth()

p6 <-
ggplot(data=lcd.dt.ag, aes(x=yday, y=precip.var, group=year, color=as.factor(year))) +
geom_line(size=.5, alpha=.5)+
geom_smooth()


layout <- "
AB
CD
EF"

mega <- p1 + p2 + p3 + p4 + p5 + p6 +
plot_layout(guides = 'collect', design=layout)

ggsave(mega, file="~/yearly_weather.png")



ggp <- ggpairs(lcd.dt.ag[,c(2,3,4,5,6,7,8,9,10,11,12)])
ggsave(ggp, file="~/ggp.png")

ggplot(data=lcd.dt.ag, aes(x=temp.mu, y=temp.var)) + geom_point()

lcd.dt.ag <- na.omit(lcd.dt.ag)
pc <- prcomp((lcd.dt.ag[,c(2,3,4,5,6,7)]), center=T, scale=T)
lcd.dt.ag[,pc1:=pc$x[,1]]
lcd.dt.ag[,pc2:=pc$x[,2]]

ggplot(data=lcd.dt.ag, aes(x=yday, y=pc2, group=year)) +
geom_line()


	 noaa = meteo_noaa_hourly(station = "724016-93736",
                  year = 2014:2019) # poznan, poland
head(noaa)



tmp <- nearest_stations_nooa(country = "UNITED STATES",
	date = Sys.Date(),
	add_map = F,
	point = c(samps[locality=="VA_ch"]$long[1], samps[locality=="VA_ch"]$lat[1]),
	no_of_stations = 3
)














meteo_noaa_hourly = function(station = NULL, year, fm12 = TRUE){
	# station = "724016-93736"
	# year = 2014:2018
  stopifnot(is.character(station))
  #options(RCurlOptions = list(ssl.verifypeer = FALSE)) # required on windows for RCurl

  base_url = "https://www1.ncdc.noaa.gov/pub/data/noaa/"

  all_data = NULL

  for (i in seq_along(year)){

      address = paste0(base_url, year[i], "/", station, "-", year[i], ".gz")
      temp = tempfile()
      test_url(address, temp)

      # run only if downloaded file is valid
      dat = NULL
      if(!is.na(file.size(temp)) & (file.size(temp) > 0)) {

      dat = read.fwf(gzfile(temp,'rt'),header=F,
                   c(4, 6, 5, 4, 2, 2, 2, 2, 1, 6,
                     7, 5, 5, 5, 4, 3, 1, 1, 4, 1,
                     5, 1, 1, 1, 6, 1, 1, 1, 5, 1, 5, 1, 5, 1))
      unlink(temp)

      if(fm12){
      dat = dat[dat$V12 == "FM-12",] # take only FM-12 records
      }

      dat = dat[, c(4:7, 10:11, 13, 16, 19, 25, 29, 31, 33)]
      colnames(dat) = c("year", "month", "day", "hour", "lat", "lon", "alt",
                        "wd", "ws",  "visibility", "t2m", "dpt2m",
                        "slp")

      dat$date = ISOdatetime(year = dat$year,
                             month = dat$month,
                             day = dat$day,
                             hour = dat$hour, 0, 0, tz = "UTC")


      dat$t2m[dat$t2m == 9999] = NA
      dat$dpt2m[dat$dpt2m == 9999] = NA
      dat$ws[dat$ws == 9999] = NA
      dat$wd[dat$wd == 999] = NA
      dat$slp[dat$slp == 99999] = NA
      dat$visibility[dat$visibility == 999999] = NA

      dat$lon = dat$lon/1000
      dat$lat = dat$lat/1000
      dat$ws = dat$ws/10
      dat$t2m = dat$t2m/10
      dat$dpt2m = dat$dpt2m/10
      dat$slp = dat$slp/10

      } else {

       cat(paste0("  Check station name or year. The created link is not working properly:\n  ", address))

      }  # end of if statement for empty files


      all_data[[length(all_data) + 1]] = dat
      } # end of loop for years

  if(is.list(all_data)){
    all_data = do.call(rbind, all_data)
  }


  if(!is.null(all_data)){ # run only if there are some data downloaded:

    # order columns:
    all_data = all_data[, c("date","year", "month", "day", "hour", "lon", "lat", "alt",
                            "t2m", "dpt2m", "ws", "wd", "slp", "visibility") ]

    # sort data
    all_data = all_data[order(all_data$date), ]

    #closeAllConnections() # just in case sth left while testing...
  }

  return(all_data)
} # koniec funkcji meteo_terminowe
