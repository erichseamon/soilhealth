#-agcensus_function
#-plots a map of a particular year of agcensus data - only works with 2012 right now
#-agcensus controls include:
#   
#"tile_farms"            "tile_acres"           "tile_acres_avgfarm"    "ditches_farms"         "ditches_acres"         
#"consease_farms"        "consease_acres"        "consease_avgfarm"      "notill_farms"          "notill_acres"         
#"notill_avgfarm"        "constill_farms"        "constill_acres"        "constill_avgfarm"      "convtill_farms"       
#"convtill_acres"        "convtill_acres.1"      "cc_farms"              "cc_acres"              "cc_avgfarm"  
#"ditches_acres_avgfarm"
#
#-example: agcensus_function(2012,tile_farms)

agcensus_function <- function(year, agcensuscontrols) {
year <- 2012
agcensuscontrols <- c("tile_farms")  

library(classInt)
library(leaflet)
library(dplyr)
library(Hmisc)
library(RColorBrewer)
library(maptools)
library(raster)

temp <- tempfile()
download.file("http://dmine.io/waf/SHEAF/spatial_data/states_conus/states_conus.zip",temp)
outDir<-"/tmp"
unzip(temp,exdir=outDir)

setwd("/tmp/states_conus")

states <- readShapePoly('states_conus.shp',
                        proj4string=CRS
                        ("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
projection = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

temp <- tempfile()
download.file("http://dmine.io/waf/SHEAF/spatial_data/counties_conus/counties_conus.zip",temp)
outDir<-"/tmp"
unzip(temp,exdir=outDir)

setwd("/tmp/counties_conus")

counties <- readShapePoly('UScounties_conus.shp',
                          proj4string=CRS
                          ("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
projection = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


n1 <- read.csv("https://dmine.io/waf/SHEAF/data/AgCensus2012.csv", header=TRUE)
n1 <- n1[,1:25]


xxx <- subset(n1, year == year)

xxx$county <- tolower(xxx$county)
xxx$county <- capitalize(xxx$county)

colnames(xxx)[4] <- "NAME"
colnames(xxx)[3] <- "STATE_NAME"

id = 5:25
xxx[id] = as.numeric(unlist(xxx[id]))


m <- merge(counties, xxx, by=c("STATE_NAME", "NAME"))


palz1 <- brewer.pal(9, "GnBu")

palz <- colorRampPalette(palz1)


palData <- classIntervals(eval(parse(text=paste("m$", agcensuscontrols, sep=""))), style="hclust")
colors <- findColours(palData, palz(100))


#pal <- colorNumeric(palette = c("white", "orange", "darkorange", "red", "darkred"),
#                    domain = eval(parse(text=paste("m$",input$agcensuscontrols , sep=""))))

pal2 <- colorNumeric(brewer.pal(11, "GnBu"), na.color = "#ffffff",
                     domain = eval(parse(text=paste("m$", agcensuscontrols, sep=""))))



exte <- as.vector(extent(states))

label <- paste(sep = "<br/>", m$STATE_NAME, round(eval(parse(text=paste("m$", agcensuscontrols, sep=""))), 0))
markers <- data.frame(label)
labs <- as.list(eval(parse(text=paste("m$", agcensuscontrols, sep=""))))


leaflet(data = m) %>% addProviderTiles("Stamen.TonerLite") %>% fitBounds(exte[1], exte[3], exte[2], exte[4]) %>% addPolygons(color = ~pal2(eval(parse(text=paste("m$", input$agcensuscontrols, sep="")))), popup = markers$label,  weight = 1) %>%
  addLegend(pal = pal2, values = ~eval(parse(text=paste("m$", agcensuscontrols, sep=""))), opacity = 1, title = NULL,
            position = "bottomright")

}