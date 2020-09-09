setwd("~/IVLarkIsotopesJoO")

### Sampling map of points from Imperial Valley ###
install.packages(c("xlsx","rvertnet","ggmap","devtools","rlang","maptools","rgdal","rgeos","raster","png","viridis"))
install_github("oswaldosantos/ggsn")

require(rvertnet)
require("ggmap")
require("devtools")
require("rlang")
require("ggsn")
require("maptools")
require(rgdal)
require(rgeos)
require(raster)
require(png)
require(viridis)

feathers<-read.csv("./Appendices/IsotopeVoucherTable_v4.csv",stringsAsFactors=F)
feathers$histfac<-factor(c("historical","contemporary")[as.numeric(feathers$Year>1950)+1])

feathers_sampsize<-table(feathers$Latitude, feathers$Longitude)

feathers_sampsize_ind<-sapply(which(feathers_sampsize!=0), function(x) arrayInd(x,dim(feathers_sampsize)))

fss_lat<-rep(NA,ncol(feathers_sampsize_ind))
fss_lon<-rep(NA,ncol(feathers_sampsize_ind))
for(i in 1:ncol(feathers_sampsize_ind)){
	fss_lat[i]<-rownames(feathers_sampsize)[feathers_sampsize_ind[1,i]]
	fss_lon[i]<-colnames(feathers_sampsize)[feathers_sampsize_ind[2,i]]	
}

feathers_ss<-data.frame(samplesize=feathers_sampsize[which(feathers_sampsize!=0)], lat=fss_lat,lon= fss_lon)

feathers_ss_histfac<-rep(NA,nrow(feathers_ss))
for(i in 1:nrow(feathers_ss)){
	latfoo<-which(feathers$Latitude == feathers_ss$lat[i])
	lonfoo<-which(feathers$Longitude == feathers_ss$lon[i])
	feathers_ss_histfac[i]<-c("contemporary","historical")[as.numeric(feathers[intersect(latfoo,lonfoo)[1],]$histfac)]
}

feathers_ss<-cbind(feathers_ss_histfac,feathers_ss)
feathers_ss$lat<-as.numeric(as.character(feathers_ss$lat))
feathers_ss$lon<-as.numeric(as.character(feathers_ss$lon))

plants<-read.csv("./Appendices/PlantIsotopeData_v1.csv",stringsAsFactors=F)
plants$Species<-factor(plants$Species)

plant_ss<-data.frame(hab=c(rep("agriculture",5),rep("scrub",5)),samplesize=rep(3,10),lat=as.numeric(unique(plants$Latitude)),lon=as.numeric(unique(plants$Longitude)))
plant_ss$hab<-factor(plant_ss $hab)

soil<-read.csv("./Appendices/SoilIsotopeData_v1.csv",stringsAsFactors=F)
soil$Habitat.Type<-factor(soil$Habitat.Type)

soil_ss<-data.frame(hab=c(rep("agriculture",5),rep("scrub",5)),samplesize=rep(3,10),lat=as.numeric(unique(soil$Latitude)),lon=as.numeric(unique(soil$Longitude)))
soil_ss$hab<-factor(soil_ss$hab)

### Set API Key ###
### Users need to get their own API key i think ###
ggmap::register_google(key = "INSERT_YOUR_API_KEYHERE")

### Read in lark PNG plate ###
larkplate<-readPNG("./alpestris.png")
g <- rasterGrob(larkplate, interpolate=TRUE)

### Create main map ###
iv <-get_map(location = c(lon = -115.25, lat = 33),
                    zoom = 9, scale = 2,
                    maptype ='satellite',
                    color = 'color',source="google")
                    
attr(iv,"bb")

### read in GADM shapefiles ###
usa1<-readOGR("~/GIS/GADM/USA_adm/USA_adm1.shp")
mex1<-readOGR("~/GIS/GADM/MEX_adm/MEX_adm1.shp")
usamex1<-raster:::bind(usa1,mex1)

usamex1_crop<-crop(usamex1,extent(c(attr(iv, "bb")$ll.lon, attr(iv, "bb")$ur.lon, attr(iv, "bb")$ll.lat, attr(iv, "bb")$ur.lat)))
usamex1_fort<-fortify(usamex1_crop)

sample_size_df<-data.frame(size=c(1,3,5,10),stringsAsFactors=F,lon=seq(-116.05,-115.75,length.out=4),lat=rep(33.625,4))


### Start plotting ###
iv_map <- 	ggmap(iv,extent="normal") +
			geom_polygon(data= usamex1_crop,aes(x=long,y=lat,group=group),col="black",fill=NA,size=2.5) +
			geom_point(data=plant_ss,aes(x=lon,y=lat),shape=22,fill=viridis(2)[plant_ss$hab],size=plant_ss$samplesize,alpha=0.8) +
			geom_point(data=feathers_ss,aes(x=lon,y=lat),shape=21,fill=viridis(2)[feathers_ss$feathers_ss_histfac],size=feathers_ss$samplesize,alpha=0.8) +

			theme(legend.position = "none")+
			labs(x="Longitude",y="Latitude")+
			
			annotate("rect", xmin = -116.105, xmax = -115.815, ymin = 32.32, ymax = 32.4,fill="white") +
			
ggsn:::scalebar(dist=10,dist_unit="km",transform=T,x.min=-116.25,x.max=-115.875,y.min=32.365,y.max=33.7,st.size=4)+
			coord_cartesian(xlim=c(attr(iv, "bb")$ll.lon, attr(iv, "bb")$ur.lon),ylim=c(attr(iv, "bb")$ll.lat, attr(iv, "bb")$ur.lat),expand=F) + 
			annotate("rect", xmin = -116.1, xmax = -115.675, ymin = 33.525, ymax = 33.71,fill="white") +
		geom_point(data=sample_size_df,aes(x=lon,y=lat),shape=1,size= sample_size_df$size,alpha=0.6) +
annotate("text",x=sample_size_df$lon,y=sample_size_df$lat-0.055,label=sample_size_df$size) +
	annotate("text",x=-115.89,y=33.68,label="Sample Size")+
	   annotation_custom(g, xmin=-114.8, xmax=-114.35, ymin=32.27, ymax=32.575)

### INSET MAP ###
usamex1_inset<-crop(usamex1,extent(c(-125,-107, 28, 50)))

inset<-	ggplot()+
		coord_fixed()+
		geom_polygon(data= usamex1_inset,aes(x=long,y=lat,group=group),col="black",fill=NA,size=0.5)+
		geom_rect(data= usamex1_inset,aes(xmin=attr(iv, "bb")$ll.lon,xmax=attr(iv, "bb")$ur.lon,ymin=attr(iv, "bb")$ll.lat,ymax=attr(iv, "bb")$ur.lat),fill=NA,col="red",size=0.65)+
		theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background = element_rect(fill = "white"))+
		theme(legend.position = "none")

### Combine into single plot ###
g1<-ggplotGrob(iv_map)

#quartz(width=6.5,height=6.5)
tiff(file="./IVSampMap_v3.tiff",width=6.5,height=6.5,units="in",res=600)
#png()
grid.draw(g1)
pushViewport(viewport(x=0.825, y=0.825,w=.25, h=.4,clip="off"))

g2 <- ggplotGrob(inset)
grid.draw(g2)


dev.off()

