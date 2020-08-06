
######### Conifer species distribution ##########
rm(list = ls())

library(data.table)
library(maptools)
library(raster)
library(sp)
library(taxize)
library(rgeos)

#function to download distribution
kewDist <- function(taxa) {
  id <- get_ids(taxa, db ="pow",rows=1)$pow[[1]]
  idn <- id[which(grepl("urn:",id))]
  nat <- pow_lookup(idn, include = "distribution")$meta$distribution$natives
  return(nat)
}


fl <- list.files("/Volumes/PiyalKaru/UPPSALA/Data/gbif", full.names = T)
xx <- NULL
for (i in seq_along(fl)){
  tp <- fread(fl[i])
  tp <- tp[,c("species","decimalLongitude","decimalLatitude")]
  xx <- rbind(xx,tp)
}

#dt <- xx[,c("species","decimalLongitude","decimalLatitude",)]
#dt <- xx[,c(10,22,23)]
dm <- SpatialPointsDataFrame(coords = xx[,c(2,3)], data =xx,
      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))


#### Native range
abies <- kewDist("Picea abies")
obovata <- kewDist("Picea obovata")
glauca_laxa <- kewDist("Picea laxa") # P. glauca is a synonym of P. laxa according to Kew db
contorta <- kewDist("Pinus contorta")
sitchensis <- kewDist("Picea sitchensis")
rubens <- kewDist("Picea rubens")

tdwg <- shapefile("/Users/piyalkaru/Desktop/R files/TDWG/level4/level4.shp")

abi.poly <- gUnaryUnion(tdwg[tdwg$Level3_cod %in% abies$tdwgCode, ])
obo.poly <- gUnaryUnion(tdwg[tdwg$Level3_cod %in% obovata$tdwgCode, ])
glau.poly <- gUnaryUnion(tdwg[tdwg$Level3_cod %in% glauca_laxa$tdwgCode,])
cont.poly <- gUnaryUnion(tdwg[tdwg$Level3_cod %in% contorta$tdwgCode, ])
sitch.poly <- gUnaryUnion(tdwg[tdwg$Level3_cod %in% sitchensis$tdwgCode, ])
rub.poly <- gUnaryUnion(tdwg[tdwg$Level3_cod %in% rubens$tdwgCode, ])


crs(dm) <- crs(tdwg)

p.abi <- subset(dm, dm$species=="Picea abies")
p.abi.cr <- cbind(p.abi, pr=over(p.abi,abi.poly))
pp.abi <- p.abi.cr@data[complete.cases(p.abi.cr@data),-4]

p.obo <- subset(dm, dm$species=="Picea obovata")
p.obo.cr <- cbind(p.obo, pr=over(p.obo,obo.poly))
pp.obo <- p.obo.cr@data[complete.cases(p.obo.cr@data),-4]

p.glau <- subset(dm, dm$species=="Picea glauca")
p.glau.cr <- cbind(p.glau, pr=over(p.glau,glau.poly))
pp.glau <- p.glau.cr@data[complete.cases(p.glau.cr@data),-4]

p.cont <- subset(dm, dm$species=="Pinus contorta")
p.cont.cr <- cbind(p.cont, pr=over(p.cont,cont.poly))
pp.cont <- p.cont.cr@data[complete.cases(p.cont.cr@data),-4]

p.sitch <- subset(dm, dm$species=="Picea sitchensis")
p.sitch.cr <- cbind(p.sitch, pr=over(p.sitch,sitch.poly))
pp.sitch <- p.sitch.cr@data[complete.cases(p.sitch.cr@data),-4]

p.rub <- subset(dm, dm$species=="Picea rubens")
p.rub.cr <- cbind(p.rub, pr=over(p.rub,rub.poly))
pp.rub <- p.rub.cr@data[complete.cases(p.rub.cr@data),-4]


pp <- rbind(pp.abi,pp.obo,pp.glau,pp.rub,pp.sitch,pp.cont)
pp <- SpatialPointsDataFrame(coords = pp[,c(2,3)], data =pp,
                             proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
w.map <- shapefile("/Users/piyalkaru/Desktop/R files/world_simple/TM_WORLD_BORDERS_SIMPL-0.3.shp")
w.map <- w.map[!w.map$NAME=="Antarctica",]

#par(mar=c(0,0,0,0))
pdf("/Volumes/PiyalKaru/UPPSALA/outputs/conifer_dist.pdf", height = 8, width = 12)
plot(w.map, col="grey", border=NA)
points(pp[,c(2,3)], col=factor(pp$species), pch=19, cex=0.3)
#par(mar=c(50,0,0,0))
legend("topright",legend = unique(pp$species), col=factor(unique(pp$species)), pch=19,
       fill = F, border = F, horiz = T, bty="n", cex = 0.8)
dev.off()








