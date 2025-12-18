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
chi.poly <- rgdal::readOGR('columbus.shp')
class(chi.poly)
str(slot(chi.poly,"data"))
summary(chi.poly@data$HOVAL)
plot(chi.poly)

#------------底图加上openstreet在线地图----------------------------------
library(leaflet)
leaflet(chi.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.5, smoothFactor = 0.5) %>%
  addTiles() 

# Fig5-a
#------------哥伦布房价HOVAL空间分布地图---------------------------------
qpal<-colorQuantile("OrRd", chi.poly@data$HOVAL, n=9) 
leaflet(chi.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.8, smoothFactor = 0.2, color = ~qpal(HOVAL)
  ) %>%
  addTiles()

# Fig5-b
#------------哥伦布收入INC空间分布地图---------------------------------
qpal<-colorQuantile("OrRd", chi.poly@data$INC, n=9) 
leaflet(chi.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.8, smoothFactor = 0.2, color = ~qpal(INC)
  ) %>%
  addTiles()

# Fig5-c
#------------哥伦布犯罪率crime分布地图----------------------------------
qpal<-colorQuantile("OrRd", chi.poly@data$CRIME, n=9) 
leaflet(chi.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.8, smoothFactor = 0.2, color = ~qpal(CRIME)
  ) %>%
  addTiles()

# Tip：注释 %>% addTiles()可隐藏背景即得论文中的图。
#------------哥伦布距离DI空间分布地图---------------------------------
qpal<-colorQuantile("OrRd", chi.poly@data$DISCBD, n=9) 
leaflet(chi.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.8, smoothFactor = 0.2, color = ~qpal(DISCBD)
  ) 
  # %>%
  # addTiles()

