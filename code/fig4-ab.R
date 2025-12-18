rm(list = ls()) 
#------------哥伦布空间邻接矩阵------------------------------------------
library('spdep')
data("columbus")
col.gal.nb
Ws = nb2listw(col.gal.nb,style = 'W')
Ws
Wn = listw2mat(Ws) 
Wn
str(columbus)
class(columbus)
summary(columbus$HOVAL)
#------------哥伦布Spatal.data.frame格式底图-----------------------------
library(rgdal)
# setEPS()
# postscript('fig4-a.eps', width=4.5, height=3.1, paper="special", horizontal=FALSE)
# par(mar=c(0,0,0,0), oma=c(0,0,0,0))
chi.poly <- rgdal::readOGR('columbus.shp')
class(chi.poly)
str(slot(chi.poly,"data"))
summary(chi.poly@data$HOVAL)
plot(chi.poly)
# box(col="black", lwd=2)
# dev.off() 
#------------底图加上openstreet在线地图----------------------------------
library(leaflet)
leaflet(chi.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.5, smoothFactor = 0.5) %>%
  addTiles() 
#------------哥伦布空间矩阵连接地图-------------------------------------
library('spdep')
# setEPS()
# postscript('fig4-b.eps', width=4.5, height=3.1, paper="special", horizontal=FALSE)
# par(mar=c(0,0,0,0), oma=c(0,0,0,0))
list.queen=poly2nb(chi.poly,queen=TRUE)
list.queen
W=nb2listw(list.queen,style="W",zero.policy=TRUE)
W
chi.poly
plot(W,coordinates(chi.poly))
# box(col="black", lwd=2)
# dev.off() 

