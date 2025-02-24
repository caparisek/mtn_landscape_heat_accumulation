# -----------------------------------------------------------------------
##
## This script: 
##     1) loads in {nhdR} data for USA
##     2) crops NHD to obtain mountain ranges for USA
##     3) attaches elevation from {elevatr}, keeps only LakePond type
##     3) write_csv a DF that is mountain lakes
##
# -----------------------------------------------------------------------


# read in packages --------------------------------------------------------
#install.packages("USAboundariesData", repos = "https://ropensci.r-universe.dev", type = "source")
#remotes::install_version("USAboundaries", "0.4.0")
library(USAboundaries) #state labels
library(tidyverse)     #basic use
library(fBasics)       #basic stats
library(sf)            #spatial manipulations
library(mapview)       #spatial 
library(leaflet)       #spatial
library(rgdal)         #spatial
library(raster)        #spatial
options(scipen=999)    #removes sci-notation
library(LaCroixColoR)  #color
library(RColorBrewer)  #color
library(elevatr)       #data - elevation from usgs + tiles
library(nhdR)          #data - https://github.com/jsta/nhdR



# data – NHD polygons -----------------------------------------------------
## NHD-metadata: https://s3.amazonaws.com/edap-nhdplus/NHDPlusV21/Documentation/NHDPlusV2_User_Guide.pdf
## Info-NHDPlus'-VPU: https://www.epa.gov/waterdata/get-nhdplus-national-hydrography-dataset-plus-data

all_vpu_shp<-nhdR::vpu_shp 
#VPU 17 - Pacific Northwest (RPU = 17a, 17b, 17c, 17d)
#VPU 18 - California (RPU = 18a,18b,18c)

# get vpu export (downloads locally; run this on a computer just 1x)
?nhd_plus_get
nhd_plus_get(vpu = 18, "NHDSnapshot")
nhd_plus_get(vpu = 18, "NHDPlusAttributes")
nhd_plus_get(vpu = 18, "NHDPlusCatchment")
nhd_plus_get(vpu = 17, "NHDSnapshot")
nhd_plus_get(vpu = 17, "NHDPlusAttributes")
nhd_plus_get(vpu = 17, "NHDPlusCatchment")

# explore - list layers
?nhd_plus_list
nhd_plus_list(vpu = 18, "NHDSnapshot")
nhd_plus_list(vpu = 18, "NHDPlusAttributes")
nhd_plus_list(vpu = 18, "NHDPlusCatchment")
nhd_plus_list(vpu = 17, "NHDSnapshot")
nhd_plus_list(vpu = 17, "NHDPlusAttributes")
nhd_plus_list(vpu = 17, "NHDPlusCatchment")

# explore - layer info
?nhd_plus_info
nhd_plus_info(vpu = 18, component="NHDSnapshot", dsn="NHDWaterbody")
nhd_plus_info(vpu = 17, component="NHDSnapshot", dsn="NHDWaterbody")

# load vpu
California<- nhd_plus_load(vpu = 18, "NHDSnapshot", "NHDWaterbody")        #18 California
PacificNorthwest <- nhd_plus_load(vpu = 17, "NHDSnapshot", "NHDWaterbody") #17 Pacific Northwest

a<- nhd_plus_load(vpu = 1, "NHDSnapshot", "NHDWaterbody") #01 Northeast
b<- nhd_plus_load(vpu = 2, "NHDSnapshot", "NHDWaterbody") #02 Mid Atlantic

cN<- nhd_plus_load(vpu = "03N", "NHDSnapshot", "NHDWaterbody") #03N South Atlantic North
cS<- nhd_plus_load(vpu = "03S", "NHDSnapshot", "NHDWaterbody") #03S South Atlantic South
cW<- nhd_plus_load(vpu = "03W", "NHDSnapshot", "NHDWaterbody") #03W South Atlantic West

d<- nhd_plus_load(vpu = 04, "NHDSnapshot", "NHDWaterbody") #04 Great Lakes
e<- nhd_plus_load(vpu = 5, "NHDSnapshot", "NHDWaterbody")  #05 Ohio
f<- nhd_plus_load(vpu = 6, "NHDSnapshot", "NHDWaterbody")  #06 Tennessee
g<- nhd_plus_load(vpu = 7, "NHDSnapshot", "NHDWaterbody")  #07 Upper Mississippi
h<- nhd_plus_load(vpu = 8, "NHDSnapshot", "NHDWaterbody")  #08 Lower Mississippi
j<- nhd_plus_load(vpu = 9, "NHDSnapshot", "NHDWaterbody")  #09 Souris-Red-Rainy

kU<- nhd_plus_load(vpu = "10U", "NHDSnapshot", "NHDWaterbody") #10U Upper Missouri 
kL<- nhd_plus_load(vpu = "10L", "NHDSnapshot", "NHDWaterbody") #10L Lower Missouri

l<- nhd_plus_load(vpu = 11, "NHDSnapshot", "NHDWaterbody") #11 Ark-Red-White
m<- nhd_plus_load(vpu = 12, "NHDSnapshot", "NHDWaterbody") #12 Texas
n<- nhd_plus_load(vpu = 13, "NHDSnapshot", "NHDWaterbody") #13 Rio Grande
o<- nhd_plus_load(vpu = 14, "NHDSnapshot", "NHDWaterbody") #14 Upper Colorado
p<- nhd_plus_load(vpu = 15, "NHDSnapshot", "NHDWaterbody") #15 Lower Colorado
q<- nhd_plus_load(vpu = 16, "NHDSnapshot", "NHDWaterbody") #16 Great Basin
#r<- nhd_plus_load(vpu = 19, "NHDSnapshot", "NHDWaterbody") #there is no #19 
s<- nhd_plus_load(vpu = 20, "NHDSnapshot", "NHDWaterbody") #warning, #20 Hawaii
t<- nhd_plus_load(vpu = 21, "NHDSnapshot", "NHDWaterbody") #21 Puerto Rico/U.S. Virgin Islands

## uA<- nhd_plus_load(vpu = "22AS", "NHDSnapshot", "NHDWaterbody") #22A American Samoa
## uG<- nhd_plus_load(vpu = "22GU", "NHDSnapshot", "NHDWaterbody") #22G Guam
## uM<- nhd_plus_load(vpu = "22MP", "NHDSnapshot", "NHDWaterbody") #22M Northern Mariana Islands

class(California) #sf
class(PacificNorthwest) #sf

st_crs(California)
st_crs(PacificNorthwest)
st_crs(a)
st_crs(b)
st_crs(cN)
st_crs(cS)
st_crs(cW)
st_crs(d)
st_crs(e)
st_crs(f)
st_crs(g)
st_crs(h)
st_crs(j)
st_crs(kU)
st_crs(kL)
st_crs(l)
st_crs(m)
st_crs(n)
st_crs(o)
st_crs(p)
st_crs(q)
st_crs(s)
st_crs(t)

# Note: Col# and Names don't match - requires fix.
colnames(California)
colnames(PacificNorthwest)
colnames(a)
colnames(b)
colnames(cN) #fix
colnames(cS) #fix
colnames(cW) #fix
colnames(d)
colnames(e)
colnames(f)
colnames(g)
colnames(h)
colnames(j)
colnames(kU)
colnames(kL)
colnames(l)
colnames(m)
colnames(n) #fix
colnames(o) #fix
colnames(p) #fix
colnames(q) #fix
colnames(s) #fix - 12 col. No ShapeLength or ShapeArea, instead "permanent_"
colnames(t) #fix - 12 col. No ShapeLength or ShapeArea, instead "permanent_"

CN<-cN %>% #newname=oldname
  dplyr::rename(COMID=ComID,
                FDATE=FDate,
                GNIS_NAME=GNIS_Name,
                AREASQKM=AreaSqKm,
                ELEVATION=Elevation,
                REACHCODE=ReachCode,
                FCODE=FCode,
                SHAPE_LENG=Shape_Leng,
                SHAPE_AREA=Shape_Area)

CS<-cS %>%
  dplyr::rename(COMID=ComID,
                FDATE=FDate,
                GNIS_NAME=GNIS_Name,
                AREASQKM=AreaSqKm,
                ELEVATION=Elevation,
                REACHCODE=ReachCode,
                FCODE=FCode,
                SHAPE_LENG=Shape_Leng,
                SHAPE_AREA=Shape_Area)

CW<-cW %>% 
  dplyr::rename(COMID=ComID,
                FDATE=FDate,
                GNIS_NAME=GNIS_Name,
                AREASQKM=AreaSqKm,
                ELEVATION=Elevation,
                REACHCODE=ReachCode,
                FCODE=FCode,
                SHAPE_LENG=Shape_Leng,
                SHAPE_AREA=Shape_Area)

N<-n %>% 
  dplyr::rename(COMID=ComID,
                FDATE=FDate,
                GNIS_NAME=GNIS_Name,
                AREASQKM=AreaSqKm,
                ELEVATION=Elevation,
                REACHCODE=ReachCode,
                FCODE=FCode,
                SHAPE_LENG=Shape_Leng,
                SHAPE_AREA=Shape_Area)

O<-o %>% 
  dplyr::rename(COMID=ComID,
                FDATE=FDate,
                GNIS_NAME=GNIS_Name,
                AREASQKM=AreaSqKm,
                ELEVATION=Elevation,
                REACHCODE=ReachCode,
                FCODE=FCode,
                SHAPE_LENG=Shape_Leng,
                SHAPE_AREA=Shape_Area)

P<-p %>% 
  dplyr::rename(COMID=ComID,
                FDATE=FDate,
                GNIS_NAME=GNIS_Name,
                AREASQKM=AreaSqKm,
                ELEVATION=Elevation,
                REACHCODE=ReachCode,
                FCODE=FCode,
                SHAPE_LENG=Shape_Leng,
                SHAPE_AREA=Shape_Area)

Q<-q %>% 
  dplyr::rename(COMID=ComID,
                FDATE=FDate,
                GNIS_NAME=GNIS_Name,
                AREASQKM=AreaSqKm,
                ELEVATION=Elevation,
                REACHCODE=ReachCode,
                FCODE=FCode,
                SHAPE_LENG=Shape_Leng,
                SHAPE_AREA=Shape_Area)

SS<-s %>% 
  dplyr::rename(COMID=ComID,
                FDATE=FDate,
                RESOLUTION=Resolution,
                GNIS_NAME=GNIS_Name,
                AREASQKM=AreaSqKm,
                ELEVATION=Elevation,
                REACHCODE=ReachCode,
                FCODE=FCode) %>% 
  dplyr::mutate(SHAPE_LENG=NA) %>% #add col other df's have by same name
  dplyr::mutate(SHAPE_AREA=NA) %>% 
  dplyr::select(-Permanent_) #remove column other df's don't have

TT<-t %>% 
  dplyr::rename(COMID=ComID,
                FDATE=FDate,
                RESOLUTION=Resolution,
                GNIS_NAME=GNIS_Name,
                AREASQKM=AreaSqKm,
                ELEVATION=Elevation,
                REACHCODE=ReachCode,
                FCODE=FCode) %>% 
  dplyr::mutate(SHAPE_LENG=NA) %>% #add col other df's have by same name
  dplyr::mutate(SHAPE_AREA=NA) %>% 
  dplyr::select(-Permanent_) #remove column other df's don't have

# bind into one NHD 
NHD<-bind_rows(California,PacificNorthwest,a,b,CN,CS,CW,d,e,f,g,h,j,kU,kL,l,m,N,O,P,Q,SS,TT)

class(NHD)
colnames(NHD)
st_is_valid(NHD)              
valid_NHD<-st_make_valid(NHD)

# NHD write_shp 
write_csv(NHD,"data output/NHD_rbind_USA.csv")
st_write(NHD,"data output/NHD_rbind_USA.shp", append=FALSE)
st_write(valid_NHD,"data output/NHD_rbind_USA_VALID.shp", append=FALSE)

# NHD map polygons - exploratory maps for interest
mapview(California,
         col.regions="red",
         map.types = c("Stamen.TerrainBackground",
                       "Esri.NatGeoWorldMap", 
                       "CartoDB.Positron", 
                       "CartoDB.DarkMatter", 
                       "Esri.WorldImagery", 
                       "OpenTopoMap")) 
mapview(PacificNorthwest,
         col.regions="red",
         map.types = c("Stamen.TerrainBackground",
                       "Esri.NatGeoWorldMap", 
                       "CartoDB.Positron", 
                       "CartoDB.DarkMatter", 
                       "Esri.WorldImagery", 
                       "OpenTopoMap")) 




# data – MTN polygons -----------------------------------------------------
# EPA ecoregions - https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states 
OmernickEcoregion <- sf::st_read("data/mtn_region_databases/EPA_ecoregion_Level4_Omernik/us_eco_l4.shp")

class(OmernickEcoregion) #sf
st_crs(OmernickEcoregion)

OmernickEcoregion_sf <- st_transform(OmernickEcoregion, crs=4269) #give .shp's same CRS
st_crs(OmernickEcoregion_sf)


# Locate mountainous ecoregions 
mapview(OmernickEcoregion_sf,
        map.types = c("Stamen.TerrainBackground",
                      "Esri.NatGeoWorldMap",
                      "OpenTopoMap"),zcol = "US_L3NAME")

# Keep mountainous ecoregions 
mountain_ecoregion<-OmernickEcoregion_sf %>% 
  dplyr::filter(NA_L3NAME %in% c("Sierra Nevada",
                                 "Cascades",
                                 "North Cascades",
                                 "Klamath Mountains",
                                 "Southern Rockies",
                                 "Columbia Mountains/Northern Rockies",
                                 "Middle Rockies",
                                 "Canadian Rockies",
                                 "Arizona/New Mexico Mountains",
                                 "Southwestern Appalachians",
                                 "Central Appalachians",
                                 "Northern Appalachian and Atlantic Maritime Highlands",
                                 "North Central Appalachians",
                                 "Blue Mountains",
                                 "Blue Ridge",
                                 "Wasatch and Uinta Mountains",
                                 "Idaho Batholith"))  

# Modify name for mountainous ecoregions to plot better
renamed_mtn<-mountain_ecoregion %>% 
  #Simple Name
  mutate(MountainName= case_when(
    NA_L3NAME=="Canadian Rockies"                    ~ "Rockies", 
    NA_L3NAME=="Columbia Mountains/Northern Rockies" ~ "Rockies",
    NA_L3NAME=="Middle Rockies"                      ~ "Rockies",
    NA_L3NAME=="Southern Rockies"                    ~ "Rockies",
    
    NA_L3NAME=="North Central Appalachians"                           ~ "Appalachians",
    NA_L3NAME=="Northern Appalachian and Atlantic Maritime Highlands" ~ "Appalachians",
    NA_L3NAME=="Central Appalachians"                                 ~ "Appalachians",
    NA_L3NAME=="Southwestern Appalachians"                            ~ "Appalachians",
    
    NA_L3NAME=="Cascades"       ~ "Cascades",
    NA_L3NAME=="North Cascades" ~ "Cascades",
    
    NA_L3NAME=="Sierra Nevada"                ~ "Sierra Nevada",
    NA_L3NAME=="Klamath Mountains"            ~ "Klamath",
    NA_L3NAME=="Blue Mountains"               ~ "Blue Mountains",
    NA_L3NAME=="Idaho Batholith"              ~ "Idaho Batholith",
    NA_L3NAME=="Blue Ridge"                   ~ "Blue Ridge",
    NA_L3NAME=="Arizona/New Mexico Mountains" ~ "AZ/NM Mountains",
    NA_L3NAME=="Wasatch and Uinta Mountains"  ~ "Wasatch-Uinta")) %>% 
  #Complex Name 
  mutate(MountainName_complex= case_when(
    NA_L3NAME=="Canadian Rockies"                    ~ "Rockies - Canadian", 
    NA_L3NAME=="Columbia Mountains/Northern Rockies" ~ "Rockies - Columbia Mtn/Northern",
    NA_L3NAME=="Middle Rockies"                      ~ "Rockies - Middle",
    NA_L3NAME=="Southern Rockies"                    ~ "Rockies - Southern",
    NA_L3NAME=="North Central Appalachians"                           ~ "Appalachians - North Central",
    NA_L3NAME=="Northern Appalachian and Atlantic Maritime Highlands" ~ "Appalachians - North Highland & Atl. Maritime",
    NA_L3NAME=="Central Appalachians"                                 ~ "Appalachians - Central",
    NA_L3NAME=="Southwestern Appalachians"                            ~ "Appalachians - Southwestern",
    NA_L3NAME=="Cascades"                 ~ "Cascades",
    NA_L3NAME=="North Cascades"           ~ "Cascades - Northern",
    NA_L3NAME=="Sierra Nevada"                  ~ "Sierra Nevada",
    NA_L3NAME=="Klamath Mountains"              ~ "Klamath",
    NA_L3NAME=="Blue Mountains"                 ~ "Blue Mountains",
    NA_L3NAME=="Idaho Batholith"                ~ "Idaho Batholith",
    NA_L3NAME=="Blue Ridge"                     ~ "Blue Ridge",
    NA_L3NAME=="Arizona/New Mexico Mountains"   ~ "AZ/NM Mountains",
    NA_L3NAME=="Wasatch and Uinta Mountains"    ~ "Wasatch-Uinta"))

mapview(renamed_mtn,
        zcol = "MountainName",
        col.regions=(c("coral",
                       "slateblue4",
                       "chartreuse",
                       "cyan",
                       "steelblue4",
                       "gold",
                       "green4",
                       "darkcyan",
                       "cornflowerblue",
                       "aliceblue")))

mapview(renamed_mtn,
        zcol = "MountainName", 
        col.regions=brewer.pal(10,"Paired"))
class(renamed_mtn)

# save – mountain range polygons
# st_write(renamed_mtn,"data output/Omernik_Mountain_Ecoregion_Polygons_sf.shp", append=FALSE)



# ~ -----------------------------------------------------------------------
# crop – NHD to MTN range -------------------------------------------------

valid_MTN<-sf::st_make_valid(renamed_mtn)
sf::sf_use_s2(FALSE) #https://stackoverflow.com/questions/68478179/how-to-resolve-spherical-geometry-failures-when-joining-spatial-data

sf_crop<- st_join(x=valid_NHD, y=valid_MTN, left = FALSE) #lefttrue = keep all. leftfalse = keep where pt falls in mtn polygon



# crop – NHD to keep only LakePond FTYPE  ---------------------------------
class(sf_crop)#sf
st_crs(sf_crop)#EPSG 4269
unique(sf_crop$FTYPE)

NHD_LakePond<-sf_crop %>% 
  dplyr::filter(FTYPE== "LakePond")
unique(NHD_LakePond$FTYPE)



# assign elevation  ---------------------------------------------------------------
## Can't use NHD elevation column because, from NHD page: "Some fields in the NHD attribute tables contain special coded values as follows: The value “-9998” signifies that the applicable value for the field is missing or undetermined. The value “-9999” signifies that there is no applicable value and one will never be assigned"

ll_prj <- "EPSG:4269"
elev_NHD2 <- elevatr::get_elev_point(NHD_LakePond, prj = ll_prj, src = "epqs") 



# save – NID + Omernik Mtns -------------------------------------------------
write_csv(elev_NHD2,"data output/NHD_mtns_omernik_state_elev.csv")
st_write(elev_NHD2,"data output/NHD_mtns_omernik_state_elev.shp", append=FALSE) 







