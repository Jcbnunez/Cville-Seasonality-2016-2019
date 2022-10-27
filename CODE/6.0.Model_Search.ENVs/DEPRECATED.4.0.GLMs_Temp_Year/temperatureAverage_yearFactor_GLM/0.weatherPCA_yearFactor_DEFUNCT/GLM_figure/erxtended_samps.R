
	### find sites with multiple time points
		samps.ag <- samps[,list(nSamps=length(locality),
							nSpring=sum(season=="spring"),
							nFall=sum(season=="fall"),
							nTime=length(unique(collectionDate)),
							maxDelta=max(yday) - min(yday),
							lat=mean(lat),
							long=mean(long)),
					list(locality, year, continent, set)]

		setkey(samps.ag, locality, year)
		setkey(samps, locality, year)




    		multi_sample <- ggplot() +
    		geom_line(data= samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality, linetype=continent)) +
    		geom_point(data=samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality, color=continent)) +
    		facet_grid(.~year) + scale_colour_manual(values = brewer.pal(8,"Set2")[c(1,2)]) +
    		scale_x_date(date_labels = "%b", limits = as.Date(c(110,355), origin = as.Date("2018-01-01"))) +
    		xlab("Collection Date") + ylab("Latitude") +
    		theme_bw() +
    		theme(axis.text = element_text(angle = 45, hjust = 1, size=12),
    					legend.text=element_text(size=12))

ggsave(multi_sample, file="~/multiSample_cville.png", h=6, w=12)
