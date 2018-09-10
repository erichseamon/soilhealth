library(classInt)
library(maptools)

setwd("/dmine/data/states/")

states <- readShapePoly('states.shp',
                          proj4string=CRS
                          ("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
projection = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#counties <- counties[grep("Idaho|Washington|Oregon|Montana", counties@data$STATE_NAME),]
#counties <- counties[grep(input$state2, counties@data$STATE_NAME),]

setwd("/dmine/data/soilses/data/")
files2 <- list.files(pattern = "nass2*")
myfiles2 = lapply(files2, read.csv, strip.white = TRUE, header = TRUE)
y <- do.call(rbind, myfiles2)
x <- as.data.frame(y)


nass1 <- read.csv("nass2_1982.csv", header = TRUE)
nass2 <- read.csv("nass2_1987.csv", header = TRUE)
nass3 <- read.csv("nass2_1992.csv", header = TRUE)
nass4 <- read.csv("nass2_1997.csv", header = TRUE)
nass5 <- read.csv("nass2_2002.csv", header = TRUE)
nass6 <- read.csv("nass2_2007.csv", header = TRUE)
nass7 <- read.csv("nass2_2012.csv", header = TRUE)

n1 <- rbind(nass1, nass2, nass3, nass4, nass5, nass6, nass7)

subset(n1, Year == "1982")


