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
	library(climate)
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
	samps[,id:=sampleId]

### Based on Machado, some samps do not have excat collection dates. We will remove those
	# https://cdn.elifesciences.org/articles/67577/elife-67577-supp1-v2.xlsx
	machado_samps <- read.xls("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/EnvironmentalOmnibus_Global/elife-67577-supp1.xlsx")
	machado_samps <- as.data.table(machado_samps)
	machado_samps[,exactDate:=!is.na(Day)]

	samps <- merge(samps, machado_samps[,c("Sample", "exactDate")], by.x="sampleId", by.y="Sample", all.x=T)
	samps[is.na(exactDate), exactDate:=T]
	table(samps$exactDate)

	samps <- samps[exactDate==T]

### test two weather stats
	getPower <- function(i) {
			daily_single_ag <- get_power(
			  community = "ag",
			  lonlat = c(samps[i]$long, samps[i]$lat),
			  pars = c("RH2M", "T2M", "PRECTOTCORR"),
			  dates = c(paste(samps[i]$year, "-01-01", sep=""), paste(samps[i]$year, "-12-31", sep="")),
			  temporal_api = "hourly",
				time_standard="UTC"
			)
			daily_single_ag <- as.data.table(daily_single_ag)
			daily_single_ag
	}

	getLCD <- function(i) {
		tmp <- nearest_stations_nooa(country = (samps[i]$COUNTRY),
		  date = Sys.Date(),
		  add_map = F,
		  point = c(samps[i]$long, samps[i]$lat),
		  no_of_stations = 5
		)
		tmp <- as.data.table(tmp)
		setnames(tmp, "distance [km]", "dist")
		tmp[,begin:=ymd(BEGIN)]
		tmp[,end:=ymd(END)]
		tmp <- tmp[begin<=samps[i]$Date][end >=samps[i]$Date][order(dist)]


		station <- paste(tmp$USAF, tmp$WBAN, sep="")
		message(paste(samps[i]$sampleId, samps[i]$year, i, dim(samps)[1], sep=" / "))
		lcd.dt <- foreach(j=1:length(station), .errorhandling="remove")%do%{
			#lcd.tmp <- as.data.table(lcd(station=station[j], year=samps[j]$year))

			path <- rnoaa:::lcd_get(station=station[j], year=samps[i]$year)
			lcd.tmp <- fread(path)
			setnames(lcd.tmp, names(lcd.tmp), tolower(names(lcd.tmp)))
			lcd.dt <- lcd.tmp[,c("date", "hourlydrybulbtemperature", "hourlyrelativehumidity", "hourlyprecipitation", "hourlywindspeed"), with=F]
			lcd.dt[,station:=station[j]]
			lcd.dt
		}
		lcd.dt <- rbindlist(lcd.dt)
		lcd.dt[,date:=ymd_hms(gsub("T", " ", date))]
		lcd.dt[,temp:=as.numeric(as.character(hourlydrybulbtemperature))]
		lcd.dt[,humidity:=as.numeric(as.character(hourlyrelativehumidity))]
		lcd.dt[,precip:=as.numeric(as.character(hourlyprecipitation))]
		lcd.dt[,wind:=as.numeric(as.character(hourlywindspeed))]
		lcd.dt[,sampleId:=samps[i]$sampleId]
		lcd.dt

	}

### use NASA power
	power.dt <- list()

	for(i in c(125:dim(samps)[1])) {
		message(paste(i, dim(samps)[1]))
		power.dt[[i]] <- getPower(i)
		Sys.sleep(5)
		power.dt[[i]][,sampleId:=samps[i]$sampleId]
		power.dt[[i]][,date:=ymd_hms(paste(paste(YEAR, MO, DY, sep="-"), paste(HR, ":00:00", sep=""), sep=" "))]
		power.dt[[i]][,queryDate:=Sys.time()]
	}

	power.dt <- rbindlist(power.dt)

	save(power.dt, file="Overwintering_18_19/EnvironmentalOmnibus_Global/nasa_power_weather.raw.Rdata")

### add in collection date
	power.dt <- merge(power.dt, samps[,c("sampleId", "Date"),with=F], by="sampleId")
	setnames(power.dt, "Date", "collectionDate")
	setnames(power.dt, c("T2M", "RH2M", "PRECTOTCORR"), c("temp", "humidity", "precip"))

### define the windows to summarize over
	sets <- data.table(mod=c(1:11),
										 start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
										 end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

### summarize
	setkey(power.dt, sampleId)
	weather.ave <- foreach(i=unique(power.dt$sampleId), .combine="rbind")%dopar%{
		# i <- unique(power.dt$sampleId)[100]
		power.tmp <- power.dt[J(i)]

		power.tmp[,delta:=julian(date)-julian(collectionDate)]

		foreach(k=1:dim(sets)[1], .combine="rbind")%do%{
			message(paste(i, k, sep=" / "))
			power.mod <- power.tmp[delta>= -1*(sets[k]$end) & delta<=-1*(sets[k]$start)]

			power.mod.ag <- power.mod[,list(dailyMax=max((temp), na.rm=T),
																	dailyMin=min((temp), na.rm=T)),
														 list(delta=round(delta))]

			data.table(sampleId=i, mod=k,
								temp.ave=mean(power.mod$temp),
								temp.var=var(power.mod$temp),
								temp.max=max(power.mod$temp),
								temp.min=min(power.mod$temp),
								temp.propMax=mean(power.mod.ag$dailyMax>32),
								temp.propmin=mean(power.mod.ag$dailyMin<5),

								humidity.ave=mean(power.mod$humidity),
								humidity.var=var(power.mod$humidity),

								precip.ave=mean(power.mod$precip),
								precip.var=var(power.mod$precip))
			}

		}

	summary(weather.ave)

### save
	save(weather.ave, file="/Users/alanbergland/Documents/GitHub/Overwintering_18_19/EnvironmentalOmnibus_Global/nasa_power.weather.mod.Rdata")

### plots
	load(file="/Users/alanbergland/Documents/GitHub/Overwintering_18_19/EnvironmentalOmnibus_Global/weather.Rdata")
	ws <- merge(weather.ave, samps[,c("sampleId", "season", "locality", "year", "lat"), with=F], by="sampleId")

	ggplot(data=ws[mod==4], aes(x=yday, y=temp.ave, group=interaction(locality, year))) +
	geom_point() +
	geom_line()

	p1 <- ggpairs(weather.ave[mod==4,c(3:12)],
	    upper = list(continuous = wrap("cor", size = 2)),
			lower = list(continuous = wrap("points", alpha = 0.3,    size=0.5))) +
			theme(strip.text = element_text(size = 5), axis.text=element_text(size=4))

	ggsave(p1, file="~/modWeather.global.png")

### check against Cville LCD dataset
	load(file="/Users/alanbergland/Documents/GitHub/Overwintering_18_19/EnvironmentalOmnibus_Global/nasa_power.weather.mod.Rdata")
	weather.nasapower <- weather.ave
	load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/EnvironmentalOmbibus_Cville/weather.Rdata")
	setnames(weather.ave, "V1", "sampleId")
	setkey(weather.nasapower, sampleId, mod)
	setkey(weather.ave, sampleId, mod)

	m <- merge(weather.nasapower, weather.ave)

	p1 <- ggplot(data=m, aes(x=temp.ave.x, y=temp.ave.y)) + geom_point()
	p2 <- ggplot(data=m, aes(x=humidity.ave.x, y=humidity.ave.y)) + geom_point()
	p3 <- ggplot(data=m, aes(x=precip.ave.x, y=precip.ave.y)) + geom_point()

	p1 + p2 + p3

	ggsave(p2, file="~/modWeather.nasapower_vs_lcd.png")
















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
