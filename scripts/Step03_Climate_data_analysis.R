# -----------------------------------------------------------------------
##
## This script: 
##      Reads in data output from script "02_Climate_Data_download"  (i.e., climate data)
##      Reads in data output from script "01_tidy_NHD_datasets"      (i.e., cropped NHD data that only includes mountain ranges)
##      Works with mountain landscape climate data (GDD/KDD, mixed-effect models, cluster, DFA, maps, skew/kurtosis...)
##
# -----------------------------------------------------------------------

# Packages ----------------------------------------------------------------

set.seed(1) 

options(scipen=999) #removes sci-notation in ggplot axes & table outputs
library(tidyverse)  #basic use
library(fBasics)    #basic stats

library(lmerTest)   #for mixed-effect models
library(MuMIn)      #for mixed-effect model's R^2

library(sf)         #spatial manipulations
#library(sp)         #spatial manipulations
library(raster)     #spatial manipulations
library(maps)       #for creating map
library(ggspatial)  #for creating map

#install.packages("devtools")                           #for map color scheme
#devtools::install_github("johannesbjork/LaCroixColoR") #for map color scheme
library(LaCroixColoR)                                  #for map color scheme

library(lubridate)  #for working with dates
library(utils)      #for the read.table and for loop read-in 
library(MASS)       #for DFA

library(cowplot)    #for stitching plots together
library(grid)       #for creating common x and y labels on plots
library(gridExtra)  #for creating common x and y labels on plots
library(stringr)    #for ggplot facet_grid(labeller = stringr::str_wrap())

# Note, if {lmerTest} is giving error when running models: (i.e., https://stackoverflow.com/questions/77481539/error-in-initializeptr-function-cholmod-factor-ldeta-not-provided-by-pack)




# Read Data – Climate CHELSA – Historic Year-Month ------------------------
a_files_historic <- list.files(paste(getwd(), "data/climate_databases/CHELSA_1980-2018_global_MeanDailyAirTemp_MONTHLY/", sep="/")) #list files in Data directory

a_image.files_historic <- a_files_historic[grepl(".tif", a_files_historic) & grepl("tas", a_files_historic)] #get only the tas.tiff files
a_image.files_historic #check - should be seeing the file name "xxx.tif"

a_tas.raw.historic <- raster::stack(paste("data/climate_databases/CHELSA_1980-2018_global_MeanDailyAirTemp_MONTHLY//", a_image.files_historic, sep="/")) #create raster stack of images


### Explore - but increases code run-time.
#      plot(a_tas.raw.historic) # 12 plots of "prec" data; 1=Jan, 12=Dec.
#      ?brick # {raster::brick} works on single multi-layer file
#      ?stack # {raster::stack} works on multiple files (e.g., 1 .tif per month = 1 year)
#      proj4string(a_tas.raw.historic) 
#      crs(a_tas.raw.historic)         
#      extent(a_tas.raw.historic)
#      class(a_tas.raw.historic)
###



# Read Data – NHD for Omernik mountain ranges ---------------------------------
b_NHD_MTN_wdups<-sf::st_read("data_output/NHD_mtns_omernik_state_elev.shp") #From Script01: NHD w Omernik Mtn Rng's crop



# Convert COMID ( = Lake ID #) to "character", not "numeric", for the models 
sum(is.na(b_NHD_MTN_wdups$COMID))                            # 0 NAs
class(b_NHD_MTN_wdups$COMID)                                 # numeric
b_NHD_MTN_wdups$COMID <- as.character(b_NHD_MTN_wdups$COMID) # make character
class(b_NHD_MTN_wdups$COMID)                                 # character




# Relabel Mountain Ranges
b_NHDClimate_MtnsRenamed<-b_NHD_MTN_wdups %>%
  # "MtnRange_SIMPLE" condenses the range names that are parsed into multiple (e.g., Rockies, Apps., Cascades)
  mutate(MtnRange_SIMPLE= case_when(NA_L3NA=="Canadian Rockies"                    ~ "Rockies",
                                    NA_L3NA=="Columbia Mountains/Northern Rockies" ~ "Rockies",
                                    NA_L3NA=="Middle Rockies"                      ~ "Rockies",
                                    NA_L3NA=="Southern Rockies"                    ~ "Rockies",
                                    
                                    NA_L3NA=="North Central Appalachians"                           ~ "Appalachians",
                                    NA_L3NA=="Northern Appalachian and Atlantic Maritime Highlands" ~ "Appalachians",
                                    NA_L3NA=="Central Appalachians"                                 ~ "Appalachians",
                                    NA_L3NA=="Southwestern Appalachians"                            ~ "Appalachians",
                                    
                                    NA_L3NA=="Cascades"       ~ "Cascades",
                                    NA_L3NA=="North Cascades" ~ "Cascades",
                                    
                                    NA_L3NA=="Sierra Nevada"                ~ "Sierra Nevada",
                                    NA_L3NA=="Klamath Mountains"            ~ "Klamath Mountains",
                                    NA_L3NA=="Blue Mountains"               ~ "Blue Mountains",
                                    NA_L3NA=="Idaho Batholith"              ~ "Idaho Batholith",
                                    NA_L3NA=="Blue Ridge"                   ~ "Blue Ridge",
                                    NA_L3NA=="Arizona/New Mexico Mountains" ~ "AZ-NM Mountains",
                                    NA_L3NA=="Wasatch and Uinta Mountains"  ~ "Wasatch-Uinta Mountains")) %>%
  # "MtnRange_COMPLEX" rewrites the existing differently so when plotting (e.g., the Rockies, Apps., Cascades) the ranges sort better in a legend.
  mutate(MtnRange_COMPLEX= case_when(NA_L3NA=="Canadian Rockies"                    ~ "Rockies - Canadian", 
                                     NA_L3NA=="Columbia Mountains/Northern Rockies" ~ "Rockies - Columbia Mtn/Northern",
                                     NA_L3NA=="Middle Rockies"                      ~ "Rockies - Middle",
                                     NA_L3NA=="Southern Rockies"                    ~ "Rockies - Southern",
                                     
                                     NA_L3NA=="North Central Appalachians"                           ~ "Appalachians - North Central",
                                     NA_L3NA=="Northern Appalachian and Atlantic Maritime Highlands" ~ "Appalachians - North & Atlantic Maritime",
                                     NA_L3NA=="Central Appalachians"                                 ~ "Appalachians - Central",
                                     NA_L3NA=="Southwestern Appalachians"                            ~ "Appalachians - Southwestern",
                                     
                                     NA_L3NA=="Cascades"       ~ "Cascades - General",
                                     NA_L3NA=="North Cascades" ~ "Cascades - Northern",
                                     
                                     NA_L3NA=="Sierra Nevada"                ~ "Sierra Nevada",
                                     NA_L3NA=="Klamath Mountains"            ~ "Klamath Mountains",
                                     NA_L3NA=="Blue Mountains"               ~ "Blue Mountains",
                                     NA_L3NA=="Idaho Batholith"              ~ "Idaho Batholith",
                                     NA_L3NA=="Blue Ridge"                   ~ "Blue Ridge",
                                     NA_L3NA=="Arizona/New Mexico Mountains" ~ "AZ-NM Mountains",
                                     NA_L3NA=="Wasatch and Uinta Mountains"  ~ "Wasatch-Uinta Mountains"))







# _________________________________ ---------------------------------------------
# Remove duplicates in NHD-Omernik  ---------------------------------------------
# NHD-Omernik file has duplicates because some ranges straddle multiple states or L3/L4 ecoregions 

# Explore:
x <- b_NHDClimate_MtnsRenamed %>% dplyr::filter(GNIS_NA=="Lake Tahoe" & !US_L3NA=="Blue Ridge")
x1 <- x %>% dplyr::filter(STATE_N=="California" & L4_KEY=="5c  Northern Sierra Upper Montane Forests") 
x2 <- x %>% dplyr::filter(STATE_N=="California" & L4_KEY=="5f  Northeastern Sierra Mixed Conifer-Pine Forests") 
x3 <- x %>% dplyr::filter(STATE_N=="Nevada"     & L4_KEY=="5f  Northeastern Sierra Mixed Conifer-Pine Forests") 

library(leaflet) # For exploratory mapping

leaflet(data = x) %>% 
  addTiles() %>% 
  addPolygons()

leaflet(data = x1) %>% # Polygon from DF x1 vs x2 is in same size/place when plotted
  addTiles() %>% 
  addPolygons()

leaflet(data = x2) %>% # Polygon from DF x1 vs x2 is in same size/place when plotted
  addTiles() %>% 
  addPolygons()

leaflet(data = x3) %>% 
  addTiles() %>% 
  addPolygons()

# Remove duplicates 
c_NHD_MTN_unique <- b_NHDClimate_MtnsRenamed[!duplicated(b_NHDClimate_MtnsRenamed$COMID),]




# Fig - Map of mtn regions ------------------------------------------------
colnames(c_NHD_MTN_unique)
x_map_states <- USAboundaries::us_states()
unique(x_map_states)

x_map_NAM <- x_map_states %>% 
  dplyr::filter(!state_name %in% c("Alaska", "Hawaii", "Puerto Rico", "Guam", "American Samoa", "Northern Mariana Islands", "U.S. Virgin Islands"))

x_map_NAM <- x_map_states %>% 
  dplyr::filter(state_name %in% c("California"))

x_temp_CA<- x %>% 
  dplyr::filter(STATE_N == "California")
ggplot()+
  geom_sf(data=x_map_NAM, fill="white", color = "gray90")+
  geom_sf(data=x_temp_CA, aes(fill=MtnRange_SIMPLE),color=NA)+
  scale_fill_viridis_d(name = "Mountain Ranges")+
  theme_minimal()





# Raster manipulations ----------------------------------------------------
# Convert shapefile to SpatialPolygonsDataFrame using as(), to convert to a "Spatial" object. This returns the same polygons, but as a SpatialPolygonsDataFrame. See: https://www.jamiecmontgomery.com/post/cropping-rasters-down-to-size/

class(c_NHD_MTN_unique)
c_NHD_MTN_unique_sp <- as(c_NHD_MTN_unique, 'Spatial') # SpatialPolygonsDataFrame
class(c_NHD_MTN_unique_sp)                             

proj4string(c_NHD_MTN_unique_sp) # Need to change this to match rasterstack object "a_tas.raw.historic"
crs(c_NHD_MTN_unique_sp)         
extent(c_NHD_MTN_unique_sp)      

d_NHD.transformed<-sp::spTransform(c_NHD_MTN_unique_sp, CRS(proj4string(a_tas.raw.historic))) # Reprojected to albers equal area projection
proj4string(d_NHD.transformed)                                                                # This should match rasterstack projection now
crs(d_NHD.transformed)






# Get lat-long coordinates from polygons ----------------------------------
# See: https://rfunctions.blogspot.com/2017/08/extracting-data-from-rasters-using.html?view=magazine

e_historic1 <- sp::coordinates(d_NHD.transformed) #c_NHD_MTN_unique_sp must be SpatialPolygonsDataFrame
e_historic2 <- sp::SpatialPoints(coords=e_historic1)
e_historic3 <- a_tas.raw.historic[e_historic2]    #a_tas.raw.historic is a raster stack
#"e_historic3" is a data frame where rows = geographic points and columns are climate vars from the raster stack

# plot(e_historic3) # This takes a while to plot; just a visual

# Get lat, long, values in columns, as DF (hmm.)
e_cbind_NHD_Climate <- cbind.data.frame(coordinates(d_NHD.transformed), e_historic3) 

#Note - Climate data TAS values are *10, because making cells "integer" (whole) rather than "double" (decimal) is less onerous on the raster. Convert to Kelvin (divide by 10)
colnames(e_cbind_NHD_Climate)
head(e_cbind_NHD_Climate)

e_Climate_Historic_Converted <- e_cbind_NHD_Climate %>% 
  rename(longitude = 1) %>% # Relabel Lat-Long column names
  rename(latitude = 2) %>% 
  mutate_at(vars(CHELSA_tas_01_1980_V.2.1:CHELSA_tas_12_2019_V.2.1),list(~./10)) # Divide by 10 to convert to KELVIN

# Convert Wide to Long
e_historic_data_longways <- gather(data  = e_Climate_Historic_Converted, 
                                   key   = timeperiod, 
                                   value = Climate_TAS, CHELSA_tas_01_1980_V.2.1:CHELSA_tas_12_2019_V.2.1, factor_key=FALSE)
# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)


# timeperiod needs to be in "date" format
## Current: CHELSA_tas_01_1980_V.2.1
## Transition step: t = 01_1980
## Needs to be: 01-01-1980 (the "day" is a generic value)
class(e_historic_data_longways$timeperiod)
e_data_dated<-e_historic_data_longways %>% 
  dplyr::mutate(t0 = substr(e_historic_data_longways$timeperiod, 15, 18)) %>%   #read char 15-18
  dplyr::mutate(t1 = substr(e_historic_data_longways$timeperiod, 12, 18)) %>%   #read char 12-18
  dplyr::mutate(t2 = str_replace(t1, pattern = "_", replacement = "-01-")) %>%  #find+replace
  dplyr::filter(!t0=="1979" ) #remove 1979; incomplete year

x <- e_data_dated %>%
  group_by(t0) %>% 
  tally() %>% 
  ungroup()

range(e_data_dated$t0)
range(e_data_dated$t1)
range(e_data_dated$t2)

# Turn the date column into date "format" 
class(e_data_dated$t2)
e_data_dated$t3 <- as.POSIXct(e_data_dated$t2, format="%m-%d-%Y")
class(e_data_dated$t3)
range(e_data_dated$t3)
e_data_dated$t4 <- as.Date(e_data_dated$t3,format="%m-%d-%Y")
class(e_data_dated$t4)

# Check Step
colnames(e_data_dated)
sum(is.na(e_data_dated$Climate_TAS))
sum(is.na(e_data_dated$t4))
sum(is.na(e_data_dated$longitude))
sum(is.na(e_data_dated$timeperiod))

# Convert DF to spatial for spatial-join (note: need to assign CoorRefSystm (CRS), e.g., "WGS84[4326]")
class(e_data_dated)      # is not sf_object - use st_as_sf() to fix this
st_crs(e_data_dated)
class(c_NHD_MTN_unique)  # is sf_object
st_crs(c_NHD_MTN_unique) # make sure CRS matches e_data_dated's

e_data_dated_sf <- st_as_sf(e_data_dated, coords = c("longitude", "latitude"), crs=4269, remove = FALSE)
class(e_data_dated_sf) # check it's an sf_object 
st_crs(e_data_dated_sf)

# Remove unnecessary T1-3 columns
colnames(e_data_dated_sf)
head(e_data_dated_sf)
e_climate_by_date<-e_data_dated_sf %>% 
  dplyr::select(-c(t1,t2,t3))
head(e_climate_by_date)





# Join NHD_MTN + Climate --------------------------------------------------
class(e_climate_by_date)
class(c_NHD_MTN_unique)

f_NHD_Climate_Historic <- sf::st_join(x = e_climate_by_date, y = c_NHD_MTN_unique, left=FALSE) #lefttrue = keep all. leftfalse = keep where x links to y. 

sum(is.na(c_NHD_MTN_unique$COMID))       # 0
sum(is.na(f_NHD_Climate_Historic$COMID)) # 0


# Drop spatial awareness --------------------------------------------------
# because it's not needed and overwhelms the computer.

f_HistoricNHDClimate_NONSPATIAL <- st_drop_geometry(f_NHD_Climate_Historic)
class(f_HistoricNHDClimate_NONSPATIAL)


# Create dates column(s) --------------------------------------------------
# Make column that is just a month and a column that is just the year. 

g_NHDClimate_MonthYear<-f_HistoricNHDClimate_NONSPATIAL %>% 
  dplyr::mutate( Year       = substr(f_HistoricNHDClimate_NONSPATIAL$timeperiod, 15, 18)) %>% # read char 15-18
  dplyr::mutate( Month_num  = substr(f_HistoricNHDClimate_NONSPATIAL$timeperiod, 12, 13)) %>% # read char 12-13
  dplyr::mutate( Month_char = case_when(Month_num=="01" ~ "JAN", 
                                        Month_num=="02" ~ "FEB", 
                                        Month_num=="03" ~ "MAR", 
                                        Month_num=="04" ~ "APR", 
                                        Month_num=="05" ~ "MAY", 
                                        Month_num=="06" ~ "JUN", 
                                        Month_num=="07" ~ "JUL", 
                                        Month_num=="08" ~ "AUG", 
                                        Month_num=="09" ~ "SEP", 
                                        Month_num=="10" ~ "OCT", 
                                        Month_num=="11" ~ "NOV", 
                                        Month_num=="12" ~ "DEC")) %>% 
  dplyr::mutate( Temp_C =   (Climate_TAS-273.15),
                 Temp_F = (((Climate_TAS-273.15) * (9/5)) + 32)) 
  
x <- g_NHDClimate_MonthYear %>% group_by(Year) %>% tally() %>% ungroup()

fBasics::basicStats(g_NHDClimate_MonthYear$Temp_F)
fBasics::basicStats(g_NHDClimate_MonthYear$Temp_C)

class(g_NHDClimate_MonthYear$Temp_C)
class(g_NHDClimate_MonthYear$Temp_F)
class(g_NHDClimate_MonthYear$Month_num) # Character - convert this to - Numeric
class(g_NHDClimate_MonthYear$Year)      # Character - convert this to - Numeric 
class(g_NHDClimate_MonthYear$COMID)     # Integer   - convert this to - Character (COMID = LakeID#)

g_NHDClimate_MonthYear$Month_num<-as.numeric(g_NHDClimate_MonthYear$Month_num)
class(g_NHDClimate_MonthYear$Month_num) # Month_num is now numeric

g_NHDClimate_MonthYear$Year<-as.numeric(g_NHDClimate_MonthYear$Year)
class(g_NHDClimate_MonthYear$Year) # Year is now numeric

g_NHDClimate_MonthYear$COMID<-as.character(g_NHDClimate_MonthYear$COMID)
class(g_NHDClimate_MonthYear$COMID) # COMID is now character; must be character/factor for model

x<-g_NHDClimate_MonthYear %>% dplyr::filter(COMID=="3800743") # Check to see each unique Lake-Year is 1 row




# ____________________________ ---------------------------------------------
# TEMPERATURE --------------------------------------------------------------
# ~ Prep for model --------------------------------------------------------

class(g_NHDClimateMonthYear_MeanLogTempC$COMID)
unique(g_NHDClimate_MonthYear$MtnRange_SIMPLE)
sum(is.na(g_NHDClimate_MonthYear$MtnRange_SIMPLE))
colnames(g_NHDClimate_MonthYear)

g_NHDClimateMonthYear_MeanLogTempC <- g_NHDClimate_MonthYear %>% 
  group_by(COMID, Year) %>%
  dplyr::mutate(MeanTempC = mean(Temp_C)) %>% 
  ungroup() %>% 
  arrange(COMID, Year, Month_num) %>% 
  distinct(COMID, Year, .keep_all = TRUE) %>% 
  dplyr::select(-c(Month_num,Month_char, Temp_F)) %>% 
  dplyr::mutate(MeanTempC_log10_xplus10 = (log10( 10 + MeanTempC)))     # no NaNs
  #dplyr::mutate(MeanTempC_log10_1      = (log10( 1 + MeanTempC))) %>%  # NaNs produced 
  #dplyr::mutate(MeanTempC_lnx1         = (log1p(MeanTempC)))           # NaNs produced 

fBasics::basicStats(g_NHDClimateMonthYear_MeanLogTempC$MeanTempC)
fBasics::basicStats(g_NHDClimateMonthYear_MeanLogTempC$MeanTempC_log10_xplus10)
#fBasics::basicStats(g_NHDClimateMonthYear_MeanLogTempC$MeanTempC_log10_1)
#fBasics::basicStats(g_NHDClimateMonthYear_MeanLogTempC$MeanTempC_lnx1)







# ~ Each range gets its own DF   -------------------------------------------
h_KLAMATH<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "Klamath Mountains")

h_SIERRA<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")

h_CASCADES_ALL<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "Cascades") 

h_ROCKIES_ALL<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "Rockies")

h_APPALACHIAN_ALL<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "Appalachians")

h_IDAHO<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "Idaho Batholith")

h_BLUERIDGE<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Ridge")

h_BLUEMTN<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Mountains")

h_AZNM<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "AZ-NM Mountains")

h_WASATCH<-g_NHDClimateMonthYear_MeanLogTempC %>% 
  dplyr::filter(MtnRange_SIMPLE == "Wasatch-Uinta Mountains")




# ~ Normality checks ----------------------------------------------
qqnorm(h_KLAMATH$Temp_C)
qqline(h_KLAMATH$Temp_C)

h_KLAMATH %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("Klamath's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
h_KLAMATH %>% # weird it produced NaNs
  ggplot(aes(x=(log1p(Temp_C))))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("Klamath's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(h_SIERRA$Temp_C)
qqline(h_SIERRA$Temp_C)
h_SIERRA %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("h_SIERRA's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(h_CASCADES_ALL$Temp_C)
qqline(h_CASCADES_ALL$Temp_C)
h_CASCADES_ALL %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("h_CASCADES_ALL's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(h_ROCKIES_ALL$Temp_C)
qqline(h_ROCKIES_ALL$Temp_C)
h_ROCKIES_ALL %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("h_ROCKIES_ALL's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(h_IDAHO$Temp_C)
qqline(h_IDAHO$Temp_C)
h_IDAHO %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("h_IDAHO's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(h_BLUEMTN$Temp_C)
qqline(h_BLUEMTN$Temp_C)
h_BLUEMTN %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("h_BLUEMTN's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(h_APPALACHIAN_ALL$Temp_C)
qqline(h_APPALACHIAN_ALL$Temp_C)
h_APPALACHIAN_ALL %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("h_APPALACHIAN_ALL's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(h_BLUERIDGE$Temp_C)
qqline(h_BLUERIDGE$Temp_C)
h_BLUERIDGE %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("h_BLUERIDGE's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(h_AZNM$Temp_C)
qqline(h_AZNM$Temp_C)
h_AZNM %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("h_AZNM's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(h_WASATCH$Temp_C)
qqline(h_WASATCH$Temp_C)
h_WASATCH %>% 
  ggplot(aes(x=Temp_C))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("h_WASATCH's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())





# ~ Models: ----------------------------------------------------------------
### Model structure: TEMP as func of YEAR. YEAR is fixed effect. LAKE is random effect (i.e., lake effects).

# ~ KLAMATH -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_KLAMATH_RslopeRintercept<-lmerTest::lmer(data = h_KLAMATH, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_KLAMATH_RslopeRintercept)
MuMIn::r.squaredGLMM(j_KLAMATH_RslopeRintercept) # R2M: marginal / R2C: conditional

# Extract random effects from model (random slope, random intercept)
ranef(j_KLAMATH_RslopeRintercept)$COMID #View
coef( j_KLAMATH_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_KLAMATH_coef <- coef(j_KLAMATH_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# * both manual and predict-method deviate at this spot. 

# .------ manual  method -------------------------------------------------------
# Create a data frame with all combinations of LAKE and YEAR
j_KLAMATH_YearCombos <- expand.grid(COMID = as.character(unique(h_KLAMATH$COMID)), Year = unique(h_KLAMATH$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_KLAMATH_YearSlopes <- left_join(j_KLAMATH_YearCombos, j_KLAMATH_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_KLAMATH_ymxb <-j_KLAMATH_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))  

# Plot
xPLOT_KLAMATH<- j_KLAMATH_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = -0.00973564, slope = 0.00063006), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("Klamath Mountains")

#Random Slope histogram
#  j_KLAMATH_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("KLAMATH slope")+
#    scale_y_continuous(labels = scales::comma_format())


# .------ predict() -----------------------------------------------------------

# Create a data frame with all combinations of Year and COMID
j_KLAMATH_predict_all_combinations <- expand.grid(Year = unique(h_KLAMATH$Year), COMID = unique(h_KLAMATH$COMID))
class(j_KLAMATH_predict_all_combinations$COMID)

# Predict outcomes for all combinations
j_KLAMATH_predict_data <- data.frame(j_KLAMATH_predict_all_combinations, 
                             Predicted = predict(j_KLAMATH_RslopeRintercept, 
                                                 newdata = j_KLAMATH_predict_all_combinations))

# Merge with the predicted data
j_KLAMATH_predict_merged_data <- merge(j_KLAMATH_predict_data, j_KLAMATH_coef, by = "COMID")

# Create a ggplot with trendlines
j_KLAMATH_predict_merged_data %>% 
  ggplot(aes(x=Year.x, y=Predicted, group=factor(COMID)))+
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_point(size=3, color="coral", alpha=0.2)+
  geom_abline(aes(intercept = -0.00973564, slope = 0.00063006), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    axis.ticks.y = element_blank(),
    legend.position= "none")+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("KLAMATH predict()")

#Random Slope histogram
#  j_KLAMATH_predict_merged_data %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("Klamath Slope PREDICT()")+
#    scale_y_continuous(labels = scales::comma_format())





# ~ SIERRA -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_SIERRA_RslopeRintercept<-lmerTest::lmer(data = h_SIERRA, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_SIERRA_RslopeRintercept)
MuMIn::r.squaredGLMM(j_SIERRA_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(j_SIERRA_RslopeRintercept)$COMID #View
coef( j_SIERRA_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_SIERRA_coef <- coef(j_SIERRA_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a data frame with all combinations of LAKE and YEAR
j_SIERRA_YearCombos <- expand.grid(COMID = as.character(unique(h_SIERRA$COMID)), Year =  unique(h_SIERRA$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_SIERRA_YearSlopes <- left_join(j_SIERRA_YearCombos, j_SIERRA_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_SIERRA_ymxb <-j_SIERRA_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))

# Plot
xPLOT_SIERRA<-j_SIERRA_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = -1.488752858, slope = 0.001288583), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("Sierra Nevada")


#Random Slope histogram
#  j_SIERRA_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("SIERRA slope")+
#    scale_y_continuous(labels = scales::comma_format())




# ~ CASCADES_ALL -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_CASCADES_ALL_RslopeRintercept<-lmerTest::lmer(data = h_CASCADES_ALL, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_CASCADES_ALL_RslopeRintercept)
MuMIn::r.squaredGLMM(j_CASCADES_ALL_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(j_CASCADES_ALL_RslopeRintercept)$COMID #View
coef( j_CASCADES_ALL_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_CASCADES_ALL_coef <- coef(j_CASCADES_ALL_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a data frame with all combinations of LAKE and YEAR
j_CASCADES_ALL_YearCombos <- expand.grid(COMID = as.character(unique(h_CASCADES_ALL$COMID)), Year = unique(h_CASCADES_ALL$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_CASCADES_ALL_YearSlopes <- left_join(j_CASCADES_ALL_YearCombos, j_CASCADES_ALL_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_CASCADES_ALL_ymxb <-j_CASCADES_ALL_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))

# Plot
xPLOT_CASCADES<-j_CASCADES_ALL_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = -0.474161145, slope = 0.000813254), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.title = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("Cascades")


#Random Slope histogram
#  j_CASCADES_ALL_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("CASCADES_ALL slope")+
#    scale_y_continuous(labels = scales::comma_format())



# ~ ROCKIES_ALL -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_ROCKIES_ALL_RslopeRintercept<-lmerTest::lmer(data = h_ROCKIES_ALL, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_ROCKIES_ALL_RslopeRintercept)
MuMIn::r.squaredGLMM(j_ROCKIES_ALL_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(j_ROCKIES_ALL_RslopeRintercept)$COMID #View
coef( j_ROCKIES_ALL_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_ROCKIES_ALL_coef <- coef(j_ROCKIES_ALL_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a data frame with all combinations of LAKE and YEAR
j_ROCKIES_ALL_YearCombos <- expand.grid(COMID = as.character(unique(h_ROCKIES_ALL$COMID)), Year = unique(h_ROCKIES_ALL$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_ROCKIES_ALL_YearSlopes <- left_join(j_ROCKIES_ALL_YearCombos, j_ROCKIES_ALL_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_ROCKIES_ALL_ymxb <-j_ROCKIES_ALL_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))

# Plot
xPLOT_ROCKIES<-j_ROCKIES_ALL_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = -1.24978940, slope = 0.00115174), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("Rockies")


#Random Slope histogram
#  j_ROCKIES_ALL_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("ROCKIES_ALL slope")+
#    scale_y_continuous(labels = scales::comma_format())



# ~ APPALACHIAN_ALL -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_APPALACHIAN_ALL_RslopeRintercept<-lmerTest::lmer(data = h_APPALACHIAN_ALL, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_APPALACHIAN_ALL_RslopeRintercept)
MuMIn::r.squaredGLMM(j_APPALACHIAN_ALL_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(j_APPALACHIAN_ALL_RslopeRintercept)$COMID #View
coef( j_APPALACHIAN_ALL_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_APPALACHIAN_ALL_coef <- coef(j_APPALACHIAN_ALL_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)


# Create a data frame with all combinations of LAKE and YEAR
j_APPALACHIAN_ALL_YearCombos <- expand.grid(COMID = as.character(unique(h_APPALACHIAN_ALL$COMID)), Year = unique(h_APPALACHIAN_ALL$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_APPALACHIAN_ALL_YearSlopes <- left_join(j_APPALACHIAN_ALL_YearCombos, j_APPALACHIAN_ALL_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_APPALACHIAN_ALL_ymxb <-j_APPALACHIAN_ALL_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))

# Plot
xPLOT_APPS<-j_APPALACHIAN_ALL_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = -0.144433594, slope = 0.000696567), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=11),
    axis.ticks.x = element_blank(),

    axis.title = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("Appalachians")



#Random Slope histogram
#  j_APPALACHIAN_ALL_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("APPALACHIAN_ALL slope")+
#    scale_y_continuous(labels = scales::comma_format())



# ~ IDAHO -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_IDAHO_RslopeRintercept<-lmerTest::lmer(data = h_IDAHO, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_IDAHO_RslopeRintercept)
MuMIn::r.squaredGLMM(j_IDAHO_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(j_IDAHO_RslopeRintercept)$COMID #View
coef( j_IDAHO_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_IDAHO_coef <- coef(j_IDAHO_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a data frame with all combinations of LAKE and YEAR
j_IDAHO_YearCombos <- expand.grid(COMID = as.character(unique(h_IDAHO$COMID)), Year = unique(h_IDAHO$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_IDAHO_YearSlopes <- left_join(j_IDAHO_YearCombos, j_IDAHO_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_IDAHO_ymxb <-j_IDAHO_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))

# Plot
xPLOT_IDAHO<-j_IDAHO_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = -1.84668882, slope = 0.00143494), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=11),
    axis.title = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("Idaho Batholith")


#Random Slope histogram
#  j_IDAHO_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("IDAHO slope")+
#    scale_y_continuous(labels = scales::comma_format())




# ~ BLUERIDGE -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_BLUERIDGE_RslopeRintercept<-lmerTest::lmer(data = h_BLUERIDGE, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_BLUERIDGE_RslopeRintercept)
MuMIn::r.squaredGLMM(j_BLUERIDGE_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(j_BLUERIDGE_RslopeRintercept)$COMID #View
coef( j_BLUERIDGE_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_BLUERIDGE_coef <- coef(j_BLUERIDGE_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a data frame with all combinations of LAKE and YEAR
j_BLUERIDGE_YearCombos <- expand.grid(COMID = as.character(unique(h_BLUERIDGE$COMID)), Year = unique(h_BLUERIDGE$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_BLUERIDGE_YearSlopes <- left_join(j_BLUERIDGE_YearCombos, j_BLUERIDGE_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_BLUERIDGE_ymxb <-j_BLUERIDGE_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))

# Plot
xPLOT_BLUERIDGE<-j_BLUERIDGE_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = 0.122740150, slope = 0.000619010), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("Blue Ridge")


#Random Slope histogram
#  j_BLUERIDGE_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("BLUERIDGE slope")+
#    scale_y_continuous(labels = scales::comma_format())




# ~ BLUEMTN -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_BLUEMTN_RslopeRintercept<-lmerTest::lmer(data = h_BLUEMTN, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_BLUEMTN_RslopeRintercept)
MuMIn::r.squaredGLMM(j_BLUEMTN_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(j_BLUEMTN_RslopeRintercept)$COMID #View
coef( j_BLUEMTN_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_BLUEMTN_coef <- coef(j_BLUEMTN_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a data frame with all combinations of LAKE and YEAR
j_BLUEMTN_YearCombos <- expand.grid(COMID = as.character(unique(h_BLUEMTN$COMID)), Year = unique(h_BLUEMTN$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_BLUEMTN_YearSlopes <- left_join(j_BLUEMTN_YearCombos, j_BLUEMTN_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_BLUEMTN_ymxb <-j_BLUEMTN_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))

# Plot
xPLOT_BLUEMTN<-j_BLUEMTN_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = -1.08649330, slope = 0.00112120), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("Blue Mountains")


#Random Slope histogram
#  j_BLUEMTN_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("BLUEMTN slope")+
#    scale_y_continuous(labels = scales::comma_format())




# ~ AZNM -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_AZNM_RslopeRintercept<-lmerTest::lmer(data = h_AZNM, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_AZNM_RslopeRintercept)
MuMIn::r.squaredGLMM(j_AZNM_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(j_AZNM_RslopeRintercept)$COMID #View
coef( j_AZNM_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_AZNM_coef <- coef(j_AZNM_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a data frame with all combinations of LAKE and YEAR
j_AZNM_YearCombos <- expand.grid(COMID = as.character(unique(h_AZNM$COMID)), Year = unique(h_AZNM$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_AZNM_YearSlopes <- left_join(j_AZNM_YearCombos, j_AZNM_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_AZNM_ymxb <-j_AZNM_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))

# Plot
xPLOT_AZNM<-j_AZNM_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = -1.056307336, slope = 0.001176783), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("AZ-NM Mountains")


#Random Slope histogram
#  j_AZNM_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("AZNM slope")+
#    scale_y_continuous(labels = scales::comma_format())



# ~ WASATCH -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

j_WASATCH_RslopeRintercept<-lmerTest::lmer(data = h_WASATCH, MeanTempC_log10_xplus10 ~ Year + (Year|COMID))
summary(j_WASATCH_RslopeRintercept)
MuMIn::r.squaredGLMM(j_WASATCH_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(j_WASATCH_RslopeRintercept)$COMID #View
coef( j_WASATCH_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
j_WASATCH_coef <- coef(j_WASATCH_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a data frame with all combinations of LAKE and YEAR
j_WASATCH_YearCombos <- expand.grid(COMID = as.character(unique(h_WASATCH$COMID)), Year = unique(h_WASATCH$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
j_WASATCH_YearSlopes <- left_join(j_WASATCH_YearCombos, j_WASATCH_coef, by = "COMID")

# Tidy up & "predict" trend over time -
j_WASATCH_ymxb <-j_WASATCH_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(MeanTempC_mxb = ((Slope * Year) + Intercept))

# Plot
xPLOT_WASATCH<-j_WASATCH_ymxb %>%
  ggplot(aes(x=Year, y=MeanTempC_mxb, group = factor(COMID)))+ 
  geom_line(linewidth=0.7, alpha=0.3, color = "seagreen")+
  geom_abline(aes(intercept = -2.49637944, slope = 0.00177777), linewidth=1.0, alpha=0.9, color="black")+
  scale_y_continuous(limits = c(0.570,1.490), breaks = seq(0.570,1.490, 0.22))+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5))+
  ylab(label = "Mean Temperature (Celcius, modeled)")+
  ggtitle("Wasatch-Uinta Mountains")


#Random Slope histogram
# j_WASATCH_ymxb %>% 
#   ggplot(aes(x=Slope)) +
#   geom_histogram(bins=15)+
#   theme_classic()+
#   ggtitle("WASATCH slope")+
#   scale_y_continuous(labels = scales::comma_format())



# Plot mountain range plots from above as one gridded plot
dev.off()
plots_tempc <- cowplot::plot_grid(xPLOT_APPS,
                                  xPLOT_AZNM,
                                  xPLOT_BLUEMTN,
                                  xPLOT_BLUERIDGE,
                                  xPLOT_CASCADES,
                                  xPLOT_IDAHO,
                                  xPLOT_KLAMATH,
                                  xPLOT_ROCKIES,
                                  xPLOT_SIERRA,
                                  xPLOT_WASATCH,ncol=5)

# Create common x and y labels
y.grob <- grid::textGrob("Mean Temperature (log1p Celcius, modeled)", gp=gpar(fontsize=15), rot=90)
x.grob <- grid::textGrob("Year", gp=gpar(fontsize=15))

# Add these to plot
plots_tempc_grid <- gridExtra::grid.arrange(arrangeGrob(plots_tempc, left = y.grob, bottom = x.grob))

# --Fig S4 ----------------------------------------------------------------
ggsave(plot = plots_tempc_grid, "figures/12.18.2024 - TempC LMER Models.jpeg", height = 5, width = 14, dpi = 300)



# * RE Slope plotted ------------------------------------------------------

unique(g_NHDClimateMonthYear_MeanLogTempC$MtnRange_SIMPLE)
colnames(j_KLAMATH_ymxb)
colnames(j_SIERRA_ymxb)
colnames(j_CASCADES_ALL_ymxb)
colnames(j_ROCKIES_ALL_ymxb)
colnames(j_IDAHO_ymxb)
colnames(j_BLUEMTN_ymxb)
colnames(j_APPALACHIAN_ALL_ymxb)
colnames(j_BLUERIDGE_ymxb)
colnames(j_AZNM_ymxb)
colnames(j_WASATCH_ymxb)

# RBind the model outputs. 

k_bound_ranges_TEMP<-bind_rows(j_KLAMATH_ymxb,
                               j_SIERRA_ymxb,
                               j_CASCADES_ALL_ymxb,
                               j_ROCKIES_ALL_ymxb,
                               j_IDAHO_ymxb,
                               j_BLUEMTN_ymxb,
                               j_APPALACHIAN_ALL_ymxb,
                               j_BLUERIDGE_ymxb,
                               j_AZNM_ymxb,
                               j_WASATCH_ymxb)

# Join this to original DF w lake attributes.

x<-k_bound_ranges_TEMP %>% dplyr::filter(COMID=="3800743") # Each Lake-Year is 1 datapoint (TOO MANY POINTS)
x<-c_NHD_MTN_unique    %>% dplyr::filter(COMID=="3800743") # Check this is the same format as above. 

class(c_NHD_MTN_unique$COMID)
c_NHD_MTN_unique$COMID<-as.character(c_NHD_MTN_unique$COMID)

k_RangesWithAttributes            <- full_join(c_NHD_MTN_unique, k_bound_ranges_TEMP, by = join_by(COMID))
k_RangesWithAttributes_NONSPATIAL <- st_drop_geometry(k_RangesWithAttributes)

class(k_RangesWithAttributes_NONSPATIAL)
colnames(k_RangesWithAttributes_NONSPATIAL)
x<-k_RangesWithAttributes_NONSPATIAL %>% dplyr::filter(COMID=="3800743") # Check; since you're plotting slope, each lake has 1 slope, keep distinct COMIDs
fBasics::basicStats(k_RangesWithAttributes_NONSPATIAL$MeanTempC_mxb) # (helps to set yaxis in plot object "plots_tempc_grid")

# Keep only distinct COMID to plot slope 
k_RangesWithAttributes_distinct <- k_RangesWithAttributes_NONSPATIAL %>% 
  distinct(COMID, .keep_all = TRUE)

x<-k_RangesWithAttributes_distinct %>% dplyr::filter(COMID=="3800743") # Check; since you're plotting slope, each lake has 1 slope, keep distinct COMIDs




# ~> Pearson's R Correlation (for Fig S5) -----------------------------------------------
head(k_RangesWithAttributes_distinct$COMID)
sum(is.na(k_RangesWithAttributes_distinct$Slope))  #1007
sum(is.na(k_RangesWithAttributes_distinct$elevatn))#1

k_A_PearsonPrep <- k_RangesWithAttributes_distinct %>% 
  dplyr::filter(!is.na(elevatn) & ! is.na(Slope)) %>% 
  mutate(Log1pSlope     = log10(Slope),
         Log1pElevation = log10(1+elevatn))

# Normality Checks
k_A_PearsonPrep %>% 
  ggplot()+
  geom_histogram(aes(x=Log1pSlope))

qqnorm(k_A_PearsonPrep$Slope)
qqline(k_A_PearsonPrep$Slope, col="red")

qqnorm(k_A_PearsonPrep$elevatn)
qqline(k_A_PearsonPrep$elevatn, col="blue")

# Split into 10 ranges 
k_BT_Sierra <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")

k_BT_Cascades <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "Cascades")

k_BT_Rockies <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "Rockies")

k_BT_Idaho <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "Idaho Batholith")

k_BT_BlueMountain <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Mountains")

k_BT_BlueRidge <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Ridge")

k_BT_Appalachians <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "Appalachians")

k_BT_Klamath <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "Klamath Mountains")

k_BT_AZNM <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "AZ-NM Mountains")

k_BT_WSUN <- k_A_PearsonPrep %>% 
  dplyr::filter(MtnRange_SIMPLE == "Wasatch-Uinta Mountains")

# Pearson's R 
corr <- cor.test(k_BT_Sierra$Log1pSlope, k_BT_Sierra$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = - 0.9296865 

corr <- cor.test(k_BT_Cascades$Log1pSlope, k_BT_Cascades$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.6224428 

corr <- cor.test(k_BT_Rockies$Log1pSlope, k_BT_Rockies$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = - 0.7903184 

corr <- cor.test(k_BT_Idaho$Log1pSlope, k_BT_Idaho$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.9422413 

corr <- cor.test(k_BT_BlueMountain$Log1pSlope, k_BT_BlueMountain$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = - 0.9276575 

corr <- cor.test(k_BT_BlueRidge$Log1pSlope, k_BT_BlueRidge$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = - 0.4714338 

corr <- cor.test(k_BT_Appalachians$Log1pSlope, k_BT_Appalachians$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = - 0.4642048 

corr <- cor.test(k_BT_Klamath$Log1pSlope, k_BT_Klamath$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = - 0.8060945 

corr <- cor.test(k_BT_AZNM$Log1pSlope, k_BT_AZNM$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = - 0.9554049 

corr <- cor.test(k_BT_WSUN$Log1pSlope, k_BT_WSUN$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = - 0.9231936 

# From Pearson's R for TEMPERATURE model: 
# Appalachians  - cor = -0.5
# AZNM          - cor = -1.0 
# Blue Mountain - cor = -0.9 
# Blue Ridge    - cor = -0.5
# Cascades      - cor = -0.6 
# Idaho         - cor = -0.9 
# Klamath       - cor = -0.8 
# Rockies       - cor = -0.8 
# Sierra        - cor = -0.9 
# WU            - cor = -0.9 

k_PearsonRs_for_facet <- c("-0.5",
                           "-1.0",
                           "-0.9",
                           "-0.5",
                           "-0.6",
                           "-0.9",
                           "-0.8",
                           "-0.8",
                           "-0.9",
                           "-0.9" )

k_PearsonRs_annotation <- data.frame(MtnRange_SIMPLE = c("Appalachians" ,
                                                         "AZ-NM Mountains",
                                                         "Blue Mountains",
                                                         "Blue Ridge",
                                                         "Cascades",
                                                         "Idaho Batholith",
                                                         "Klamath Mountains",
                                                         "Rockies",              
                                                         "Sierra Nevada",      
                                                         "Wasatch-Uinta Mountains"), 
                                     label = paste("R =", k_PearsonRs_for_facet),
                                     x=Inf, 
                                     y=Inf)

k_RangesWithAttributes_distinct %>%
  ggplot()+
  geom_point( aes(x=Slope, y=elevatn), alpha=0.3, size=1, color="seagreen" )+
  scale_y_continuous(breaks = seq(0, 4500, 500), labels = scales::comma_format())+
  scale_x_continuous(breaks = seq(0, 0.0019, 0.00062), limits = c(0, 0.0019))+
  facet_wrap(.~MtnRange_SIMPLE, ncol=5)+
  geom_label(data=k_PearsonRs_annotation, 
             aes(x=x, y=y, label=label),
             inherit.aes = FALSE,
             hjust=1, vjust=1, size=4, fill=NA, color="black", label.size=0)+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title = element_text(size=13),
        strip.text = element_text(size=13),
        panel.grid.minor = element_blank())+
  ylab(label = "Elevation (m)")

# --Fig S5 ----------------------------------------------------------------
ggsave("figures/12.18.2024 - TempC (Mean, Log10) - Elevation ~ REslope.jpeg",  height = 5, width = 12,  dpi = 300)





# ___________________________ ---------------------------------------------
# Clear RAM  --------------------------------------------------------------

#sapply(ls(), function(x) object.size(get(x))) %>% sort %>% tail(20) # ID which objects take up most space
rm(a_tas.raw.historic)
rm(b_NHD_MTN_wdups)
rm(b_NHDClimate_MtnsRenamed)
rm(d_NHD.transformed)
rm(e_historic1)
rm(e_historic2)
rm(e_historic3)
rm(f_HistoricNHDClimate_NONSPATIAL)
rm(e_climate_by_date)
rm(e_data_dated_sf)
rm(f_NHD_Climate_Historic)
rm(e_Climate_Historic_Converted)
rm(e_historic_data_longways)
rm(e_data_dated)
rm(e_cbind_NHD_Climate )
rm(a_files_historic )
rm(a_image.files_historic)
#gc() # Free Unused Memory







# ___________________________ ---------------------------------------------
# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
# KKD / GDD ...  ----------------------------------------------------------
## Pick up from the object before temp models, "g_NHDClimate_MonthYear".

# Explore temp histograms: 
fBasics::basicStats(g_NHDClimate_MonthYear$Temp_C)
g_NHDClimate_MonthYear %>% 
  ggplot()+
  geom_histogram(aes(x=Temp_C),bins=170)+
  scale_x_continuous(breaks=seq(-30,35,2))+
  scale_y_continuous(breaks=seq(0,700000,100000), labels = scales::comma)+
  theme_bw()




# .--- Calc thresholds ----------------------------------------------------

# First establish that each lake currently has 12 months of data: 
x <- g_NHDClimate_MonthYear %>% 
  dplyr::filter(COMID=="3800747") %>% 
  group_by(Year) %>% 
  tally() %>% 
  ungroup()

x <- g_NHDClimate_MonthYear %>% 
  dplyr::filter(MtnRange_SIMPLE=="Blue Ridge")
head(x$COMID)

# Set desired quantiles
q = c(.25, .5, .70, .75, .80, .85, .90, .95)

m_kdd_gdd_criteria <- g_NHDClimate_MonthYear %>%
  #dplyr::filter(COMID=="3800747") %>% 
  group_by(MtnRange_SIMPLE) %>% 
  dplyr::mutate(quant25 = quantile(Temp_C, probs = q[1]),
                quant50 = quantile(Temp_C, probs = q[2]),
                quant70 = quantile(Temp_C, probs = q[3]),
                quant75 = quantile(Temp_C, probs = q[4]),
                quant80 = quantile(Temp_C, probs = q[5]),
                quant85 = quantile(Temp_C, probs = q[6]),
                quant90 = quantile(Temp_C, probs = q[7]),
                quant95 = quantile(Temp_C, probs = q[8])) %>% 
  ungroup() %>% 
  dplyr::group_by(COMID, Year, Month_char) %>%
  dplyr::mutate(#gdd_neg5_threshold     = ((((min(Temp_C))+(max(Temp_C)))/2) - (-5)) ,
                #gdd_5_threshold        = ((((min(Temp_C))+(max(Temp_C)))/2) - (5))  ,
                #gdd_10_threshold       = ((((min(Temp_C))+(max(Temp_C)))/2) - (10)) ,
                gdd_0_threshold       = (((((Temp_C))) - (0)))      ,
                
                #kdd_25Perc_threshold   = ((max(Temp_C))- (quant25)) ,
                #kdd_50Perc_threshold   = ((max(Temp_C))- (quant50)) ,
                #kdd_70Perc_threshold   = ((max(Temp_C))- (quant70)) ,
                #kdd_75Perc_threshold   = ((max(Temp_C))- (quant75)) ,
                #kdd_80Perc_threshold   = ((max(Temp_C))- (quant80)) ,
                #kdd_85Perc_threshold   = ((max(Temp_C))- (quant85)) ,
                kdd_90Perc_threshold   = (((Temp_C))- (quant90)) ,
                kdd_95Perc_threshold   = (((Temp_C))- (quant95)) )%>% 
  dplyr::ungroup()

# Check that each Lake-Year-Month has unique threshhold
x<-m_kdd_gdd_criteria %>% 
  dplyr::filter(COMID=="3800747" | COMID == "14981164" | COMID == "4710204") %>% 
  dplyr::select(c(COMID,Year,Month_char,Temp_C, gdd_0_threshold, kdd_90Perc_threshold, quant90)) %>% 
  arrange(Year, Month_char)

# Use [m_historic_quantiles] to calc projected KDD later in script:
m_historic_quantiles <- g_NHDClimate_MonthYear %>%
  #dplyr::filter(COMID=="3800747" | COMID == "14981164" | COMID == "4710204") %>% 
  group_by(MtnRange_SIMPLE) %>% 
  dplyr::summarise(quant25 = quantile(Temp_C, probs = q[1]),
                   quant50 = quantile(Temp_C, probs = q[2]),
                   quant70 = quantile(Temp_C, probs = q[3]),
                   quant75 = quantile(Temp_C, probs = q[4]),
                   quant80 = quantile(Temp_C, probs = q[5]),
                   quant85 = quantile(Temp_C, probs = q[6]),
                   quant90 = quantile(Temp_C, probs = q[7]),
                   quant95 = quantile(Temp_C, probs = q[8])) %>% 
  ungroup()

write_csv(m_historic_quantiles, "data_output/Historic_Quantile_Per_MountainRange.csv")




# .--- Replace neg. w/ 0 --------------------------------------------------
# Not logical to have a negative GDD or KDD; convert negatives to zeros.

x_count_negative <- sum(m_kdd_gdd_criteria$gdd_0_threshold < 0) # Checking zeros pre-step
x_count_negative                                                # Checking zeros pre-step
x_count_zero <- sum(m_kdd_gdd_criteria$gdd_0_threshold == 0)    # Checking zeros pre-step
x_count_zero                                                    # Checking zeros pre-step
x_count_negative+x_count_zero                                   # Checking zeros pre-step

m_kdd_gdd_criteria_negs_are_zero <- m_kdd_gdd_criteria %>% 
  mutate(across(c(#gdd_neg5_threshold,
                  gdd_0_threshold,  
                  #gdd_0_threshold_mean,
                  #gdd_5_threshold,       
                  #gdd_10_threshold,      
                  #kdd_25Perc_threshold,
                  #kdd_50Perc_threshold,
                  #kdd_70Perc_threshold,
                  #kdd_75Perc_threshold,
                  #kdd_80Perc_threshold,
                  #kdd_85Perc_threshold,
                  kdd_90Perc_threshold,
                  kdd_95Perc_threshold), ~ifelse( . < 0, 0, . )))
  
x_count_zero <- sum(m_kdd_gdd_criteria_negs_are_zero$gdd_0_threshold == 0) # Checking zeros post-step
x_count_zero                                                               # Checking zeros post-step

# Checking which have the zeros: 
x <- m_kdd_gdd_criteria_negs_are_zero %>% 
  dplyr::filter(gdd_0_threshold == 0)
fBasics::basicStats(x$Temp_C)
(unique(x$COMID))





# .--- Expand & Fill KDD & GDD days ---------------------------------------
# Each lake-year-month has 1 row of data. Expand to 30 days in order to get GDD/KDD sums for each waterbody. 
m_expansionprep_1 <- m_kdd_gdd_criteria_negs_are_zero %>% 
  mutate(Date = ymd(t4))

m_expansionprep_2 <- m_expansionprep_1 %>%
  mutate(DaysInMonth = lubridate::days_in_month(make_date(Year, Month_num, 1)))

# x <- m_kddgdd_summarized %>% # Check the below behavior before running "distinct()"
#   arrange(COMID, Year) %>% 
#   dplyr::select(c(COMID, Year, Month_num, Month_GDD_0,Month_KDD_90Perc,Sum_GDD_0,Sum_KDD_90Perc,
#                   TempC_Annual_Mean,TempC_Annual_Min,TempC_Annual_Max,TempC_Annual_SD,TempC_Annual_CV))

m_kddgdd_summarized <- m_expansionprep_2 %>% 
  #dplyr::filter(COMID =="3800743") %>% 
  ungroup() %>% 
  dplyr::mutate(#Month_GDD_neg5     = (DaysInMonth * gdd_neg5_threshold    ),
                Month_GDD_0        = (DaysInMonth * gdd_0_threshold       ),
                #Month_GDD_5        = (DaysInMonth * gdd_5_threshold       ),
                #Month_GDD_10       = (DaysInMonth * gdd_10_threshold      ),
                #Month_KDD_25Perc   = (DaysInMonth * kdd_25Perc_threshold  ),
                #Month_KDD_50Perc   = (DaysInMonth * kdd_50Perc_threshold  ),
                #Month_KDD_70Perc   = (DaysInMonth * kdd_70Perc_threshold  ),
                #Month_KDD_75Perc   = (DaysInMonth * kdd_75Perc_threshold  ),
                #Month_KDD_80Perc   = (DaysInMonth * kdd_80Perc_threshold  ),
                #Month_KDD_85Perc   = (DaysInMonth * kdd_85Perc_threshold  ),
                Month_KDD_90Perc   = (DaysInMonth * kdd_90Perc_threshold  ),
                Month_KDD_95Perc   = (DaysInMonth * kdd_95Perc_threshold  )) %>% 
  group_by(COMID, Year) %>% 
  dplyr::mutate(#Sum_GDD_neg5     = sum(Month_GDD_neg5    ),
                Sum_GDD_0         = sum(Month_GDD_0       ),
                #Sum_GDD_5        = sum(Month_GDD_5       ),
                #Sum_GDD_10       = sum(Month_GDD_10      ),
                #Sum_KDD_25Perc   = sum(Month_KDD_25Perc  ),
                #Sum_KDD_50Perc   = sum(Month_KDD_50Perc  ),
                #Sum_KDD_70Perc   = sum(Month_KDD_70Perc  ),
                #Sum_KDD_75Perc   = sum(Month_KDD_75Perc  ),
                #Sum_KDD_80Perc   = sum(Month_KDD_80Perc  ),
                #Sum_KDD_85Perc   = sum(Month_KDD_85Perc  ),
                Sum_KDD_90Perc    = sum(Month_KDD_90Perc  ),
                Sum_KDD_95Perc    = sum(Month_KDD_95Perc  ),
                TempC_Annual_Mean = mean(Temp_C),  # (By Lake, Year) based on the 12months of a Lake-Year
                TempC_Annual_Min  = min(Temp_C) ,
                TempC_Annual_Max  = max(Temp_C) ,
                TempC_Annual_SD   = sd(Temp_C)  ,
                TempC_Annual_CV   = ((TempC_Annual_SD/TempC_Annual_Mean) * 100 )) %>% 
  ungroup() %>% 
  dplyr::mutate(Log1p_SumGDD = log1p(Sum_GDD_0),
                Log1p_KDD90  = log1p(Sum_KDD_90Perc) ) %>% 
  dplyr::distinct(COMID, Year, .keep_all = TRUE) %>% # No more Year-Month uniques; 1 365-day sum per year.
  dplyr::select(-c(Date, timeperiod, t4, Month_num, Month_char,
                   Climate_TAS, Temp_C, Temp_F, 
                   quant25,quant50,quant70,quant75,quant80,quant85,quant90,quant95,
                   #gdd_neg5_threshold,
                   gdd_0_threshold,
                   #gdd_5_threshold,
                   #gdd_10_threshold,
                   #kdd_25Perc_threshold,
                   #kdd_50Perc_threshold,
                   #kdd_70Perc_threshold,
                   #kdd_75Perc_threshold,
                   #kdd_80Perc_threshold,
                   #kdd_85Perc_threshold,
                   kdd_90Perc_threshold,
                   kdd_95Perc_threshold,
                   #Month_GDD_neg5 ,
                   Month_GDD_0     ,
                   #Month_GDD_5     ,
                   #Month_GDD_10    ,
                   #Month_KDD_25Perc,
                   #Month_KDD_50Perc,
                   #Month_KDD_70Perc,
                   #Month_KDD_75Perc,
                   #Month_KDD_80Perc,
                   #Month_KDD_85Perc,
                   Month_KDD_90Perc,
                   Month_KDD_95Perc)) #removing columns based on the year-month (not applicable now)

fBasics::basicStats(m_kddgdd_summarized$Log1p_SumGDD)
range(m_kddgdd_summarized$Log1p_SumGDD)

# Check on the above
# x <- m_kddgdd_summarized %>% 
#   dplyr::filter(COMID == "7951572" |  COMID == "7951700" |  COMID == "8062947" ) %>% # check on arbitrary COMIDs
#   dplyr::select(c(Date, COMID,Year, Month_num, DaysInMonth, gdd_neg5_threshold, gdd_0_threshold,gdd_5_threshold, 
#                   Month_GDD_neg5,Month_GDD_0,Month_GDD_5,Sum_GDD_neg5,Sum_GDD_0)) %>% 
#   arrange(COMID, Date)
  




# .--- * Florida plot ----------------------------------------------------------

x<- m_kddgdd_summarized %>% group_by(COMID) %>% tally() %>% ungroup() # SumGDD is "per unique lake-year combinations"

x_plot_gddelev<-m_kddgdd_summarized %>% 
  ggplot()+
  geom_point(aes(x=Sum_GDD_0, y=elevatn), alpha=0.2, color = "navy")+
  scale_y_continuous(labels = scales::comma_format(),breaks = seq(0,4000,1000), limits=c(0,4000))+
  scale_x_continuous(labels = scales::comma_format())+
  facet_wrap(.~MtnRange_SIMPLE, ncol = 5)+
  theme_bw()+
  xlab(label = "Growing Degree Day")+
  ylab(label = "Elevation (m)")+
  theme(strip.background = element_rect(fill=NA, color=NA),
        legend.position  = "none")

m_kddgdd_summarized_forplotting<-m_kddgdd_summarized 

x_plot_kddelev<-m_kddgdd_summarized_forplotting %>% 
  #dplyr::filter(!Sum_KDD_90Perc==0) %>% 
  ggplot()+
  geom_point(aes(x=Sum_KDD_90Perc, y=elevatn), alpha=0.2, color = "navy")+
  scale_y_continuous(labels = scales::comma_format(),breaks = seq(0,4000,1000), limits=c(0,4000))+
  scale_x_continuous(labels = scales::comma_format())+
  facet_wrap(.~MtnRange_SIMPLE, ncol = 5)+
  theme_bw()+
  xlab(label = "Killing Degree Day")+
  ylab(label = "Elevation (m)")+
  theme(strip.background = element_rect(fill=NA, color=NA),
        legend.position  = "none")

x_plot_gddkdd_byelev <- cowplot::plot_grid(x_plot_gddelev,x_plot_kddelev,ncol=1)

# Fig S2 ------------------------------------------------------------------
ggsave(plot = x_plot_gddkdd_byelev, "figures/Figure S2 now includes KDD zeros - 06.17.2025 - Elev ~ SumGDD0 and SumKDD90 - with SumkDD0==0 removed.jpeg", height = 8, width = 9, dpi = 300)







# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------


# Models ------------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

# .---------- Each range gets its own df ----------------------------------
# & leave in GDD zero on purpose # fBasics::basicStats(m_kddgdd_summarized$Sum_GDD_0)

x<-m_kddgdd_summarized %>% dplyr::filter(COMID=="3800743")

m_KLAMATH<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "Klamath Mountains")

m_SIERRA<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")

m_CASCADES_ALL<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "Cascades") 

m_ROCKIES_ALL<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "Rockies")

m_APPALACHIAN_ALL<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "Appalachians")

m_IDAHO<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "Idaho Batholith")

m_BLUERIDGE<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Ridge")

m_BLUEMTN<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Mountains")

m_AZNM<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "AZ-NM Mountains")

m_WASATCH<-m_kddgdd_summarized %>% 
  dplyr::filter(MtnRange_SIMPLE == "Wasatch-Uinta Mountains")


# .---------- Normality Check ----------------------------------------------
qqnorm(m_KLAMATH$Sum_GDD_0)
qqline(m_KLAMATH$Sum_GDD_0)
m_KLAMATH %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("Klamath's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_KLAMATH %>% 
  ggplot(aes(x=(log1p(Sum_GDD_0))))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("Klamath's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
qqnorm(m_KLAMATH$Log1p_SumGDD)
qqline(m_KLAMATH$Log1p_SumGDD)
m_KLAMATH %>% 
  ggplot(aes(x=Log1p_SumGDD))+
  geom_histogram(bins=150,color="gray90", fill="purple")+
  ggtitle("Klamath's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

fBasics::basicStats(m_SIERRA$Sum_GDD_0)
fBasics::basicStats(m_SIERRA$Log1p_SumGDD)
qqnorm(m_SIERRA$Sum_GDD_0)
qqline(m_SIERRA$Sum_GDD_0)
m_SIERRA %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("m_SIERRA's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_SIERRA %>% 
  ggplot(aes(x=Log1p_SumGDD))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("m_SIERRA's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(m_CASCADES_ALL$Sum_GDD_0)
qqline(m_CASCADES_ALL$Sum_GDD_0)
m_CASCADES_ALL %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("m_CASCADES_ALL's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_CASCADES_ALL %>% 
  ggplot(aes(x=log1p(Sum_GDD_0)))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("m_CASCADES_ALL's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(m_ROCKIES_ALL$Sum_GDD_0)
qqline(m_ROCKIES_ALL$Sum_GDD_0)
m_ROCKIES_ALL %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("m_ROCKIES_ALL's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_ROCKIES_ALL %>% 
  ggplot(aes(x=log1p(Sum_GDD_0)))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("m_ROCKIES_ALL's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(m_IDAHO$Sum_GDD_0)
qqline(m_IDAHO$Sum_GDD_0)
m_IDAHO %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("m_IDAHO's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_IDAHO %>% 
  ggplot(aes(x=log1p(Sum_GDD_0)))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("m_IDAHO's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(m_BLUEMTN$Sum_GDD_0)
qqline(m_BLUEMTN$Sum_GDD_0)
m_BLUEMTN %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("m_BLUEMTN's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_BLUEMTN %>% 
  ggplot(aes(x=log1p(Sum_GDD_0)))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("m_BLUEMTN's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(m_APPALACHIAN_ALL$Sum_GDD_0)
qqline(m_APPALACHIAN_ALL$Sum_GDD_0)
m_APPALACHIAN_ALL %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("m_APPALACHIAN_ALL's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_APPALACHIAN_ALL %>% 
  ggplot(aes(x=log1p(Sum_GDD_0)))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("m_APPALACHIAN_ALL's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(m_BLUERIDGE$Sum_GDD_0)
qqline(m_BLUERIDGE$Sum_GDD_0)
m_BLUERIDGE %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("m_BLUERIDGE's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_BLUERIDGE %>% 
  ggplot(aes(x=log1p(Sum_GDD_0)))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("m_BLUERIDGE's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(m_AZNM$Sum_GDD_0)
qqline(m_AZNM$Sum_GDD_0)
m_AZNM %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("m_AZNM's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_AZNM %>% 
  ggplot(aes(x=log1p(Sum_GDD_0)))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("m_AZNM's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())

qqnorm(m_WASATCH$Sum_GDD_0)
qqline(m_WASATCH$Sum_GDD_0)
m_WASATCH %>% 
  ggplot(aes(x=Sum_GDD_0))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue")+
  ggtitle("m_WASATCH's distribution of TempC (nonlog)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())
m_WASATCH %>% 
  ggplot(aes(x=log1p(Sum_GDD_0)))+
  geom_histogram(bins=150,color="gray90", fill="dodgerblue4")+
  ggtitle("m_WASATCH's distribution of TempC (log1p)")+
  theme_classic()+
  scale_y_continuous(labels = scales::comma_format())


# ~ KLAMATH -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_KLAMATH$COMID)

m_KLAMATH_RslopeRintercept <-lmerTest::lmer(data = m_KLAMATH, Log1p_SumGDD ~ Year + (Year|COMID))
summary(m_KLAMATH_RslopeRintercept)
MuMIn::r.squaredGLMM(m_KLAMATH_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(m_KLAMATH_RslopeRintercept)$COMID #View
coef( m_KLAMATH_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_KLAMATH_coef <- coef(m_KLAMATH_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# * both manual- and predict- method deviate at this spot. 




# .------ manual method ---------------------------------------------------

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_KLAMATH$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_KLAMATH_YearCombos <- expand.grid(COMID = as.character(unique(m_KLAMATH$COMID)), 
                                    Year = unique(m_KLAMATH$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_KLAMATH_YearSlopes <- left_join(m_KLAMATH_YearCombos, m_KLAMATH_coef, by = "COMID")

# Tidy up & "predict" trend over time
m_KLAMATH_ymxb <-m_KLAMATH_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot: Random Effects
xPLOT_gdd_KLAMATH<- m_KLAMATH_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  3.07317691, slope =  0.00246226), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.title = element_blank(),
    legend.position= "none",
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Klamath Mountains")

# Plot: Random Slope histogram
m_KLAMATH_ymxb %>% 
  ggplot(aes(x=Slope)) +
  geom_histogram(bins=15)+
  theme_classic()+
  ggtitle("KLAMATH slope")+
  scale_y_continuous(labels = scales::comma_format())


# .------ predict() method ------------------------------------------------

# Create a data frame with all combinations of Year and COMID
m_KLAMATH_predict_all_combinations <- expand.grid(Year = unique(m_KLAMATH$Year), COMID = unique(m_KLAMATH$COMID))
class(m_KLAMATH_predict_all_combinations$COMID)

# Predict outcomes for all combinations
m_KLAMATH_predict_data <- data.frame(m_KLAMATH_predict_all_combinations, 
                             Predicted = predict(m_KLAMATH_RslopeRintercept, 
                                                 newdata = m_KLAMATH_predict_all_combinations))

# Merge with the predicted data
m_KLAMATH_predict_merged_data <- merge(m_KLAMATH_predict_data, m_KLAMATH_coef, by = "COMID")
colnames(m_KLAMATH_predict_merged_data)

# Plot: Random Effects from predict() method [same plot as manual method]
m_KLAMATH_predict_merged_data %>% 
  ggplot(aes(x=Year.x, y=Predicted, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  #geom_abline(aes(intercept =  xxx, slope =  xxx), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    axis.text = element_text(),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Klamath predict()")

# Plot: Random Slope histogram
m_KLAMATH_predict_merged_data %>% 
  ggplot(aes(x=Slope)) +
  geom_histogram(bins=15)+
  theme_classic()+
  ggtitle("Klamath Slope PREDICT()")+
  scale_y_continuous(labels = scales::comma_format())




# ~ SIERRA -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_SIERRA$COMID)

m_SIERRA_RslopeRintercept<-lmerTest::lmer(data = m_SIERRA, Log1p_SumGDD ~ Year + (Year|COMID))
summary(m_SIERRA_RslopeRintercept)
MuMIn::r.squaredGLMM(m_SIERRA_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(m_SIERRA_RslopeRintercept)$COMID #View
coef( m_SIERRA_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_SIERRA_coef <- coef(m_SIERRA_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_SIERRA$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_SIERRA_YearCombos <- expand.grid(COMID = as.character(unique(m_SIERRA$COMID)),
                                   Year = unique(m_SIERRA$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_SIERRA_YearSlopes <- left_join(m_SIERRA_YearCombos, m_SIERRA_coef, by = "COMID")

# Tidy up & "predict" trend over time
m_SIERRA_ymxb <-m_SIERRA_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot
xPLOT_gdd_SIERRA<-m_SIERRA_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  -2.62386708, slope =  0.00504712), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Sierra Nevada")

# Plot: Random Slope histogram
#  m_SIERRA_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("SIERRA slope")+
#    scale_y_continuous(labels = scales::comma_format())




# ~ CASCADES_ALL -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_CASCADES_ALL$COMID)

m_CASCADES_ALL_RslopeRintercept<-lmerTest::lmer(data = m_CASCADES_ALL, Log1p_SumGDD ~ Year + (Year|COMID))
summary(m_CASCADES_ALL_RslopeRintercept)
MuMIn::r.squaredGLMM(m_CASCADES_ALL_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(m_CASCADES_ALL_RslopeRintercept)$COMID #View
coef( m_CASCADES_ALL_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_CASCADES_ALL_coef <- coef(m_CASCADES_ALL_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_CASCADES_ALL$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_CASCADES_ALL_YearCombos <- expand.grid(COMID = as.character(unique(m_CASCADES_ALL$COMID)), 
                                         Year = unique(m_CASCADES_ALL$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_CASCADES_ALL_YearSlopes <- left_join(m_CASCADES_ALL_YearCombos, m_CASCADES_ALL_coef, by = "COMID")

# Tidy up & "predict" trend over time -
m_CASCADES_ALL_ymxb <-m_CASCADES_ALL_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot
xPLOT_gdd_CASCADES<-m_CASCADES_ALL_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  -0.78510562, slope =  0.00417306), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Cascades")


# Plot: Random Slope histogram
#  m_CASCADES_ALL_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("CASCADES_ALL slope")+
#    scale_y_continuous(labels = scales::comma_format())





# ~ ROCKIES_ALL -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_ROCKIES_ALL$COMID)

m_ROCKIES_ALL_RslopeRintercept<-lmerTest::lmer(data = m_ROCKIES_ALL, Log1p_SumGDD ~ Year + (Year|COMID))
summary(m_ROCKIES_ALL_RslopeRintercept)
MuMIn::r.squaredGLMM(m_ROCKIES_ALL_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(m_ROCKIES_ALL_RslopeRintercept)$COMID #View
coef( m_ROCKIES_ALL_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_ROCKIES_ALL_coef <- coef(m_ROCKIES_ALL_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_ROCKIES_ALL$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_ROCKIES_ALL_YearCombos <- expand.grid(COMID = as.character(unique(m_ROCKIES_ALL$COMID)), 
                                             Year = unique(m_ROCKIES_ALL$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_ROCKIES_ALL_YearSlopes <- left_join(m_ROCKIES_ALL_YearCombos, m_ROCKIES_ALL_coef, by = "COMID")

# Tidy up & "predict" trend over time -
m_ROCKIES_ALL_ymxb <-m_ROCKIES_ALL_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot
xPLOT_gdd_ROCKIES<-m_ROCKIES_ALL_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  -0.8208767, slope =  0.0041255), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Rockies")


# Plot: Random Slope histogram
#  m_ROCKIES_ALL_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("ROCKIES_ALL slope")+
#    scale_y_continuous(labels = scales::comma_format())



# ~ APPALACHIAN_ALL -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_APPALACHIAN_ALL$COMID)

m_APPALACHIAN_ALL_RslopeRintercept<-lmerTest::lmer(data = m_APPALACHIAN_ALL, Log1p_SumGDD ~ Year + (Year|COMID))

# Extract random effects from model (random slope, random intercept)
ranef(m_APPALACHIAN_ALL_RslopeRintercept)$COMID #View
coef( m_APPALACHIAN_ALL_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_APPALACHIAN_ALL_coef <- coef(m_APPALACHIAN_ALL_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_APPALACHIAN_ALL$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_APPALACHIAN_ALL_YearCombos <- expand.grid(COMID = as.character(unique(m_APPALACHIAN_ALL$COMID)), 
                                                 Year = unique(m_APPALACHIAN_ALL$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_APPALACHIAN_ALL_YearSlopes <- left_join(m_APPALACHIAN_ALL_YearCombos, m_APPALACHIAN_ALL_coef, by = "COMID")

# Tidy up & "predict" trend over time -
m_APPALACHIAN_ALL_ymxb <-m_APPALACHIAN_ALL_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot
xPLOT_gdd_APPS<-m_APPALACHIAN_ALL_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  3.482333164, slope =  0.002315684), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Appalachians")

# Plot: Random Slope histogram
#  m_APPALACHIAN_ALL_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("APPALACHIAN_ALL slope")+
#    scale_y_continuous(labels = scales::comma_format())




# ~ IDAHO -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_IDAHO$COMID)

m_IDAHO_RslopeRintercept<-lmerTest::lmer(data = m_IDAHO, Log1p_SumGDD ~ Year + (Year|COMID))
summary(m_IDAHO_RslopeRintercept)
MuMIn::r.squaredGLMM(m_IDAHO_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(m_IDAHO_RslopeRintercept)$COMID #View
coef( m_IDAHO_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_IDAHO_coef <- coef(m_IDAHO_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_IDAHO$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_IDAHO_YearCombos <- expand.grid(COMID = as.character(unique(m_IDAHO$COMID)), 
                                  Year = unique(m_IDAHO$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_IDAHO_YearSlopes <- left_join(m_IDAHO_YearCombos, m_IDAHO_coef, by = "COMID")

# Tidy up & "predict" trend over time -
m_IDAHO_ymxb <-m_IDAHO_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot
xPLOT_gdd_IDAHO<-m_IDAHO_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  -4.20546649, slope =  0.00574502), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=11),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Idaho Batholith")


# Plot: Random Slope histogram
#  m_IDAHO_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("IDAHO slope")+
#    scale_y_continuous(labels = scales::comma_format())




# ~ BLUERIDGE -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_BLUERIDGE$COMID)

m_BLUERIDGE_RslopeRintercept<-lmerTest::lmer(data = m_BLUERIDGE, Log1p_SumGDD ~ Year + (Year|COMID))
summary(m_BLUERIDGE_RslopeRintercept)
MuMIn::r.squaredGLMM(m_BLUERIDGE_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(m_BLUERIDGE_RslopeRintercept)$COMID #View
coef( m_BLUERIDGE_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_BLUERIDGE_coef <- coef(m_BLUERIDGE_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_BLUERIDGE$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_BLUERIDGE_YearCombos <- expand.grid(COMID = as.character(unique(m_BLUERIDGE$COMID)), 
                                      Year = unique(m_BLUERIDGE$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_BLUERIDGE_YearSlopes <- left_join(m_BLUERIDGE_YearCombos, m_BLUERIDGE_coef, by = "COMID")

# Tidy up & "predict" trend over time -
m_BLUERIDGE_ymxb <-m_BLUERIDGE_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot
xPLOT_gdd_BLUERIDGE<-m_BLUERIDGE_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  3.56317296, slope =  0.00245128), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Blue Ridge")


# Plot: Random Slope histogram
#  m_BLUERIDGE_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("BLUERIDGE slope")+
#    scale_y_continuous(labels = scales::comma_format())



# ~ BLUEMTN -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_BLUEMTN$COMID)

m_BLUEMTN_RslopeRintercept<-lmerTest::lmer(data = m_BLUEMTN, Log1p_SumGDD ~ Year + (Year|COMID))
summary(m_BLUEMTN_RslopeRintercept)
MuMIn::r.squaredGLMM(m_BLUEMTN_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(m_BLUEMTN_RslopeRintercept)$COMID #View
coef( m_BLUEMTN_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_BLUEMTN_coef <- coef(m_BLUEMTN_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_BLUEMTN$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_BLUEMTN_YearCombos <- expand.grid(COMID = as.character(unique(m_BLUEMTN$COMID)), 
                                    Year = unique(m_BLUEMTN$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_BLUEMTN_YearSlopes <- left_join(m_BLUEMTN_YearCombos, m_BLUEMTN_coef, by = "COMID")

# Tidy up & "predict" trend over time -
m_BLUEMTN_ymxb <-m_BLUEMTN_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot
xPLOT_gdd_BLUEMTN<-m_BLUEMTN_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  -1.52816339, slope =  0.00461952), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Blue Mountains")


# Plot: Random Slope histogram
#  m_BLUEMTN_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("BLUEMTN slope")+
#    scale_y_continuous(labels = scales::comma_format())




# ~ AZNM -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_AZNM$COMID)

m_AZNM_RslopeRintercept<-lmerTest::lmer(data = m_AZNM, Log1p_SumGDD ~ Year + (Year|COMID))
summary(m_AZNM_RslopeRintercept)
MuMIn::r.squaredGLMM(m_AZNM_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(m_AZNM_RslopeRintercept)$COMID #View
coef( m_AZNM_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_AZNM_coef <- coef(m_AZNM_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_AZNM$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_AZNM_YearCombos <- expand.grid(COMID = as.character(unique(m_AZNM$COMID)), 
                                 Year = unique(m_AZNM$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_AZNM_YearSlopes <- left_join(m_AZNM_YearCombos, m_AZNM_coef, by = "COMID")

# Tidy up & "predict" trend over time -
m_AZNM_ymxb <-m_AZNM_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot
xPLOT_gdd_AZNM<-m_AZNM_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  -1.86469727, slope =  0.00504423), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("AZ/NM Mountains")

# Plot: Random Slope histogram
#  m_AZNM_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("AZNM slope")+
#    scale_y_continuous(labels = scales::comma_format())





# ~ WASATCH -----------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(m_WASATCH$COMID)

m_WASATCH_RslopeRintercept<-lmerTest::lmer(data = m_WASATCH, Log1p_SumGDD ~ Year + (Year|COMID))
summary(m_WASATCH_RslopeRintercept)
MuMIn::r.squaredGLMM(m_WASATCH_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(m_WASATCH_RslopeRintercept)$COMID #View
coef( m_WASATCH_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
m_WASATCH_coef <- coef(m_WASATCH_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(m_WASATCH$Year)

# Create a data frame with all combinations of LAKE and YEAR
m_WASATCH_YearCombos <- expand.grid(COMID = as.character(unique(m_WASATCH$COMID)), 
                                    Year = unique(m_WASATCH$Year))

# Perform a left join to retain the Random Intercept and Random Slope values
m_WASATCH_YearSlopes <- left_join(m_WASATCH_YearCombos, m_WASATCH_coef, by = "COMID")

# Tidy up & "predict" trend over time -
m_WASATCH_ymxb <-m_WASATCH_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept)) 

# Plot
xPLOT_gdd_WASATCH<-m_WASATCH_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "navy")+ 
  geom_abline(aes(intercept =  -5.14105300, slope =  0.00629277), linewidth=2, alpha=0.99, color="black")+
  scale_y_continuous(limits = c(6.260, 9.065), breaks = seq(6.260,9.065, 0.66))+
  scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  theme_bw()+
  theme(
    plot.title = element_text(size=14, hjust=0.5),
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    legend.position= "none")+
  ylab(label = "Growing Degree Days (log1p, modeled)")+
  ggtitle("Wasatch-Uinta Mountains")


# Plot: Random Slope histogram
#  m_WASATCH_ymxb %>% 
#    ggplot(aes(x=Slope)) +
#    geom_histogram(bins=15)+
#    theme_classic()+
#    ggtitle("WASATCH slope")+
#    scale_y_continuous(labels = scales::comma_format())




# ~ plot all modeled ranges -----------------------------------------------

dev.off()
plots_gdd0<-cowplot::plot_grid(xPLOT_gdd_APPS,
                               xPLOT_gdd_AZNM,
                               xPLOT_gdd_BLUEMTN,
                               xPLOT_gdd_BLUERIDGE,
                               xPLOT_gdd_CASCADES,
                               xPLOT_gdd_IDAHO,
                               xPLOT_gdd_KLAMATH,
                               xPLOT_gdd_ROCKIES,
                               xPLOT_gdd_SIERRA,
                               xPLOT_gdd_WASATCH,
                               ncol=5)

y.grob <- grid::textGrob("Growing Degree Days (log1p, modeled)", gp=gpar(fontsize=15), rot=90)
x.grob <- grid::textGrob("Year", gp=gpar(fontsize=15))
plots_gdd0_grid <- gridExtra::grid.arrange(arrangeGrob(plots_gdd0, left = y.grob, bottom = x.grob))

# Fig 2 -------------------------------------------------------------------
ggsave(plot = plots_gdd0_grid, "figures/12.19.2024 - LMER Log1pGDD ~ Year.jpeg", height = 4, width = 12, dpi = 300)





# Clear RAM ---------------------------------------------------------------

#sapply(ls(), function(x) object.size(get(x))) %>% sort %>% tail(15)
rm(plots)
rm(plots_kdd0)
rm(x_plot_gddelev)
rm(x_plot_kddelev)
#rm(m_kdd_gdd_criteria_negs_are_zero)
rm(m_expansionprep_1)
rm(xPLOT_gdd_APPS)
rm(xPLOT_gdd_AZNM)
rm(xPLOT_gdd_BLUEMTN)
rm(xPLOT_gdd_BLUERIDGE)
rm(xPLOT_gdd_CASCADES)
rm(xPLOT_gdd_IDAHO)
rm(xPLOT_gdd_KLAMATH)
rm(xPLOT_gdd_ROCKIES)
rm(xPLOT_gdd_SIERRA)
rm(xPLOT_gdd_WASATCH)
#gc() # Free Unused Memory



# . - . - . - . - . - . - . - . - . - . - . - . - . ---------------------------
# .--- df prep; meanGDD; distinct COMID ---------------------------------------

# RBind the model outputs. 
n_BoundRanges_GDD0<-bind_rows(m_KLAMATH_ymxb,
                               m_SIERRA_ymxb,
                               m_CASCADES_ALL_ymxb,
                               m_ROCKIES_ALL_ymxb,
                               m_IDAHO_ymxb,
                               m_BLUEMTN_ymxb,
                               m_APPALACHIAN_ALL_ymxb,
                               m_BLUERIDGE_ymxb,
                               m_AZNM_ymxb,
                               m_WASATCH_ymxb)

fBasics::basicStats(n_BoundRanges_GDD0$SumGDD0_mxb)

# Check, w arbitrary COMID
x<-n_BoundRanges_GDD0 %>% dplyr::filter(COMID=="3800743") # Each Lake-Year is 1 datapoint 
x<-c_NHD_MTN_unique   %>% dplyr::filter(COMID=="3800743") 

# Join this to original DF w lake attributes.
class(n_BoundRanges_GDD0$COMID)
class(c_NHD_MTN_unique$COMID)

n_RangeAttributes <- left_join(c_NHD_MTN_unique, n_BoundRanges_GDD0, by = join_by(COMID))

# Drop spatial awareness (again, because c_ had it)
class(n_RangeAttributes)
n_RangeAttributes_NONSPATIAL <- st_drop_geometry(n_RangeAttributes)
class(n_RangeAttributes_NONSPATIAL)

n_Mtn_Attrib_Slopes <- n_RangeAttributes_NONSPATIAL %>% dplyr::filter(!is.na(Slope))

fBasics::basicStats(n_RangeAttributes_NONSPATIAL$Slope)
fBasics::basicStats(n_Mtn_Attrib_Slopes$Slope)



# .---calc meangdd --------------------------------------------------------

# Link "MEAN GDD" calc from pre-model data to Lake-Int-Slope
###  Remove SumGDDCol and give each lake unique line. (Each lake already had 1 unique int + slope)

# "m_" was df that was used to split ranges pre-model; so has SumGDD(GDD), MeanTemp, MinTemp for each Distinct Lake-Year
m_MEANgdd <- m_kddgdd_summarized %>% 
  ungroup() %>% 
  group_by(COMID) %>% 
  mutate(Mean_GDD                   = mean(Sum_GDD_0),
         Mean_KDD                   = mean(Sum_KDD_90Perc),
         TempC_DistinctMeanFor_Mean = mean(TempC_Annual_Mean),
         TempC_DistinctMeanFor_Min  = mean(TempC_Annual_Min),
         TempC_DistinctMeanFor_Max  = mean(TempC_Annual_Max),
         TempC_DistinctMeanFor_SD   = mean(TempC_Annual_SD),
         TempC_DistinctMeanFor_CV   = mean(TempC_Annual_CV)) %>% 
  ungroup() %>% 
  dplyr::distinct(COMID, .keep_all = TRUE) %>% # 1 row per lake
  dplyr::select(c(COMID, latitude, longitude, 
                  Mean_GDD, Mean_KDD,
                  TempC_DistinctMeanFor_Mean,
                  TempC_DistinctMeanFor_Min,
                  TempC_DistinctMeanFor_Max,
                  TempC_DistinctMeanFor_SD,
                  TempC_DistinctMeanFor_CV)) 

n_unique_LakeSlopeInt_noYear<-n_Mtn_Attrib_Slopes %>% 
  group_by(COMID) %>% 
  dplyr::distinct(COMID, .keep_all = TRUE) %>% # 1 row per lake
  dplyr::select(-c(Year, SumGDD0_mxb)) %>% 
  ungroup()

n_Means_with_LakeIntSlope <- left_join(n_unique_LakeSlopeInt_noYear,m_MEANgdd, by=join_by(COMID))

x<-n_Means_with_LakeIntSlope %>% dplyr::filter(COMID=="3800743")


### Is there a latitude effect? Yes.
# plot(n_Means_with_LakeIntSlope$latitude,n_Means_with_LakeIntSlope$Mean_GDD)
# n_Means_with_LakeIntSlope %>% 
#   ggplot(aes(x=latitude, y=Mean_GDD))+
#   geom_point()+
#   geom_smooth(method="lm", se=TRUE, color="blue")+
#   theme_minimal()
# model <- lm(Mean_GDD ~ latitude, data = n_Means_with_LakeIntSlope)
# summary(model)
# cor.test(n_Means_with_LakeIntSlope$latitude,n_Means_with_LakeIntSlope$Mean_GDD)



# * RE Slope plotted ------------------------------------------------------
# Plot - RE Slope x Elevation (slope from SumGDD0 model)

# ~> Pearson's R Correlation (for Fig 3) -----------------------------------------------
sum(is.na(n_Means_with_LakeIntSlope$Slope))  #0
sum(is.na(n_Means_with_LakeIntSlope$elevatn))#1

n_A_PearsonReady <- n_Means_with_LakeIntSlope %>% 
  dplyr::filter(!is.na(elevatn)) %>% 
  mutate(Log1pSlope = log10(Slope),
         Log1pElevation = log10(1+elevatn))

# Normality checks
n_A_PearsonReady %>% 
  ggplot()+
  geom_histogram(aes(x=Log1pSlope))

qqnorm(n_A_PearsonReady$Slope)
qqline(n_A_PearsonReady$Slope, col="red")

qqnorm(n_A_PearsonReady$elevatn)
qqline(n_A_PearsonReady$elevatn, col="blue")

# Split into 10 ranges 
n_B_Sierra <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")

n_B_Cascades <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "Cascades")

n_B_Rockies <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "Rockies")

n_B_Idaho <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "Idaho Batholith")

n_B_BlueMountain <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Mountains")

n_B_BlueRidge <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Ridge")

n_B_Appalachians <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "Appalachians")

n_B_Klamath <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "Klamath Mountains")

n_B_AZNM <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "AZ-NM Mountains")

n_B_WSUN <- n_A_PearsonReady %>% 
  dplyr::filter(MtnRange_SIMPLE == "Wasatch-Uinta Mountains")

# Pearson's R 
corr <- cor.test(n_B_Sierra$Log1pSlope, n_B_Sierra$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.9468623 

corr <- cor.test(n_B_Cascades$Log1pSlope, n_B_Cascades$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.6353264

corr <- cor.test(n_B_Rockies$Log1pSlope, n_B_Rockies$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.8101993 

corr <- cor.test(n_B_Idaho$Log1pSlope, n_B_Idaho$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.9318127

corr <- cor.test(n_B_BlueMountain$Log1pSlope, n_B_BlueMountain$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.9361022 

corr <- cor.test(n_B_BlueRidge$Log1pSlope, n_B_BlueRidge$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.5067348

corr <- cor.test(n_B_Appalachians$Log1pSlope, n_B_Appalachians$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.489865

corr <- cor.test(n_B_Klamath$Log1pSlope, n_B_Klamath$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.8125765 

corr <- cor.test(n_B_AZNM$Log1pSlope, n_B_AZNM$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.956019 

corr <- cor.test(n_B_WSUN$Log1pSlope, n_B_WSUN$Log1pElevation, method = "pearson", use = "complete.obs")
print(corr) # cor = -0.9062461 

n_PearsonRs_GDDnum <- c("-0.5",#App # cor    = -0.5
                           "-1.0",#AZ # cor     = -1.0  
                           "-0.9",#BM ) # cor   = -0.9 
                           "-0.5",#BR # cor     = -0.5
                           "-0.6",#CAS # cor    = -0.6
                           "-0.9",#id # cor     = -0.9
                           "-0.8",#K # cor      = -0.8
                           "-0.8",#Rocl # cor   = -0.8 
                           "-0.9",#Sierra # cor = -0.9 
                           "-0.9" )#Wu # cor     = -0.9 

n_PearsonRs_GDD_ann <- data.frame(MtnRange_SIMPLE = c("Appalachians" ,
                                                         "AZ-NM Mountains",
                                                         "Blue Mountains",
                                                         "Blue Ridge",
                                                         "Cascades",
                                                         "Idaho Batholith",
                                                         "Klamath Mountains",
                                                         "Rockies",              
                                                         "Sierra Nevada",      
                                                         "Wasatch-Uinta Mountains"), 
                                     label = paste("R =", n_PearsonRs_GDDnum),
                                     x=Inf, 
                                     y=Inf)


n_N_of_Lakes_for_Fig_2_and_3 <- c("10,467", #Appalachians
                                  "1,033",  #AZ-NM Mountains
                                  "284",   #Blue Mountains            
                                  "464",   #Blue Ridge                
                                  "2,165",  #Cascades                 
                                  "1,035",  #Idaho Batholith          
                                  "245",   #Klamath Mountains         
                                  "9,661",  #Rockies                  
                                  "2,358",  #Sierra Nevada
                                  "988")   #Wasatch-Uinta Mountains

n_N_of_Lakes_for_Fig_2_and_3_ann <- data.frame(MtnRange_SIMPLE = c("Appalachians" ,
                                                      "AZ-NM Mountains",
                                                      "Blue Mountains",
                                                      "Blue Ridge",
                                                      "Cascades",
                                                      "Idaho Batholith",
                                                      "Klamath Mountains",
                                                      "Rockies",              
                                                      "Sierra Nevada",      
                                                      "Wasatch-Uinta Mountains"), 
                                  label = paste("n =", n_N_of_Lakes_for_Fig_2_and_3),
                                  x=Inf, 
                                  y=Inf)







n_Means_with_LakeIntSlope %>%
  ggplot()+
  geom_point( aes(x=Slope, y=elevatn), alpha=0.1, color="navy")+
  scale_y_continuous(breaks = seq(0, 4500, 500), labels = scales::comma_format())+
  scale_x_continuous(breaks = seq(0.0019, 0.007, 0.0017), limits = c(0.0019,0.007))+
  facet_wrap(.~MtnRange_SIMPLE, ncol=5)+
  geom_label(data=n_PearsonRs_GDD_ann, 
             aes(x=x, y=y, label=label),
             inherit.aes = FALSE,
             hjust=1, vjust=1, size=4, fill=NA, color="black", label.size=0)+
  geom_label(data=n_N_of_Lakes_for_Fig_2_and_3_ann, 
             aes(x=x, y=y, label=label),
             inherit.aes = FALSE,
             hjust=1, vjust=1.8, size=4, fill=NA, color="black", label.size=0)+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title = element_text(size=13),
        strip.text = element_text(size=13),
        panel.grid.minor = element_blank())+
  ylab(label = "Elevation (m)")

x <- n_Means_with_LakeIntSlope %>% 
  group_by(MtnRange_SIMPLE) %>% 
  tally() %>% 
  ungroup()

x



# Fig 3 -------------------------------------------------------------------
ggsave("figures/12.19.2024 - GDD0 RE Slope ~ Elevation.jpeg", height = 5, width = 12, dpi = 300)

n_Means_with_LakeIntSlope %>%
  ggplot(aes(x=reorder(MtnRange_SIMPLE, -Slope, FUN = median), 
             y=Slope))+
  geom_boxplot( alpha=0.99, fill="navy")+
  #geom_boxplot(coef=0,outlier.shape=NA)+
  scale_y_continuous(breaks = seq(0.0019, 0.007, 0.0017), limits = c(0.0019,0.007))+
  theme_bw()+
  theme(axis.title.y = element_text(size=18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13, angle = 45, vjust = 1, hjust=1),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab(label = "Velocity of Change")

### 1-way anova: Is there a sig diff bw slopes of the 10 ranges? Yes. 
# aov <- aov(Slope ~ MtnRange_SIMPLE, data=n_Means_with_LakeIntSlope)
# summary(aov)

# Fig 4 -------------------------------------------------------------------
ggsave("figures/12.19.2024 - GDD0 RE Slope ~ Boxplot Range.jpeg", height = 10, width = 10, dpi = 300)






# Clear RAM  --------------------------------------------------------------

#sapply(ls(), function(x) object.size(get(x))) %>% sort %>% tail(20)
#rm(m_kdd_gdd_criteria_negs_are_zero)
#rm(m_expansionprep_2)
rm(m_APPALACHIAN_ALL)
rm(m_ROCKIES_ALL )
rm(m_CASCADES_ALL )
rm(m_KLAMATH )
rm(m_AZNM )
rm(m_BLUEMTN )
rm(m_BLUERIDGE )
rm(m_IDAHO )
rm(m_WASATCH )
rm(m_SIERRA )
rm(n_RangeAttributes_NONSPATIAL)
rm( m_ROCKIES_ALL_ymxb )
rm( m_CASCADES_ALL_ymxb )
rm( m_APPALACHIAN_ALL_ymxb )
rm( m_SIERRA_ymxb )
rm( m_KLAMATH_ymxb )
rm( m_BLUEMTN_ymxb )
rm( m_BLUERIDGE_ymxb )
rm( m_WASATCH_ymxb )
rm( m_IDAHO_ymxb )
rm( m_AZNM_ymxb )
rm(n_Plot_slopeelev)
rm(plot)
rm(m_APPALACHIAN_ALL_YearSlopes)
rm(m_ROCKIES_ALL_YearSlopes )
rm(m_CASCADES_ALL_YearSlopes)
rm(m_SIERRA_YearSlopes)
rm(m_KLAMATH_YearSlopes)
rm(m_BLUEMTN_YearSlopes)
rm(m_BLUERIDGE_YearSlopes)
rm(m_WASATCH_YearSlopes)
rm(m_IDAHO_YearSlopes)
rm(m_AZNM_YearSlopes)
rm(m_CASCADES_ALL_RslopeRintercept)
rm(m_APPALACHIAN_ALL_RslopeRintercept)
rm(m_ROCKIES_ALL_RslopeRintercept)
rm(m_SIERRA_RslopeRintercept)
rm(m_KLAMATH_RslopeRintercept)
rm(m_IDAHO_RslopeRintercept)
rm(m_WASATCH_RslopeRintercept)
rm(m_AZNM_RslopeRintercept)
rm(m_BLUERIDGE_RslopeRintercept)
rm(m_BLUEMTN_RslopeRintercept)
# gc() # Free Unused Memory







# .----------------------------- ------------------------------------------
# Projected Data ----------------------------------------------------------
# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
 

# .--- Read Data - ssp370 ------------------------------------------------------
p_ssp370 <- list.files(paste("data/climate_databases/CHELSA_climatologies_ssp370_tas__2011-2040__2041-2070__2071-2100/", sep="/")) #list files in Data directory

p_image.files.ssp370 <- p_ssp370[grepl(".tif", p_ssp370) & grepl("tas", p_ssp370)] #get only the tas.tiff files
p_image.files.ssp370 #check

p_tas.ssp370 <- raster::stack(paste("data/climate_databases/CHELSA_climatologies_ssp370_tas__2011-2040__2041-2070__2071-2100/", p_image.files.ssp370, sep="/")) #create raster stack of images





# .--- Link NHD polygons to Climatologies Data ---------------------------------
class(c_NHD_MTN_unique)
#already done# c_NHD_MTN_unique_sp <- as(c_NHD_MTN_unique, 'Spatial') # SpatialPolygonsDataFrame
proj4string(c_NHD_MTN_unique_sp) # Need to change to match rasterstack a_tas.raw.historic
crs(c_NHD_MTN_unique_sp)         # Check behavior
extent(c_NHD_MTN_unique_sp)      # Check behavior
class(c_NHD_MTN_unique_sp)       # Check behavior
p_ssp370.x0.NHD<-sp::spTransform(c_NHD_MTN_unique_sp, CRS(proj4string(p_tas.ssp370))) # Reprojected to albers equal area projection
proj4string(p_ssp370.x0.NHD)                                                          # This should match rasterstack projection now
crs(p_ssp370.x0.NHD)




# Get lat-long coordinates from polygons
p_ssp370.x1<-sp::coordinates(p_ssp370.x0.NHD) #NHD_MTN_sp must be SpatialPolygonsDataFrame
p_ssp370.x2<-sp::SpatialPoints(coords=p_ssp370.x1)
p_ssp370.x3<-p_tas.ssp370[p_ssp370.x2] #p_tas.ssp370 is a raster stack
# get lat, long, values in columns, as DF (hmm.)
p_ssp370.x4 <- cbind.data.frame(coordinates(p_ssp370.x0.NHD),p_ssp370.x3) #https://gis.stackexchange.com/questions/227370/using-r-to-extract-data-from-worldclim




#Climate data TAS is *10 because making cells "integer" (whole) rather than "double" (decimal) is less onerous on the raster. Convert to Kelvin (divide by 10)
colnames(p_ssp370.x4)

p_ssp370.x5<-p_ssp370.x4 %>% 
  rename(longitude = 1) %>% # Relabel Lat-Long column names
  rename(latitude = 2) %>% 
  mutate_at(vars(CHELSA_gfdl.esm4_r1i1p1f1_w5e5_ssp370_tas_01_2011_2040_norm:CHELSA_gfdl.esm4_r1i1p1f1_w5e5_ssp370_tas_12_2071_2100_norm),list(~./10))





# Convert Wide to Long
p_ssp370.x6.data_long <- gather(data=p_ssp370.x5, 
                              key = timeperiod, 
                              value = Climate_TAS, 
                              CHELSA_gfdl.esm4_r1i1p1f1_w5e5_ssp370_tas_01_2011_2040_norm:CHELSA_gfdl.esm4_r1i1p1f1_w5e5_ssp370_tas_12_2071_2100_norm, 
                              factor_key=FALSE)
# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)





# timeperiod needs to be in "date" format
## Current: CHELSA_gfdl.esm4_r1i1p1f1_w5e5_ssp370_tas_01_2011_2040_norm
## Transition step: t = 01_2011_2040
## Needs to be: 01-01-1980 (the "day" is a dummy; but needs this format for heavewaveR)
class(p_ssp370.x6.data_long$timeperiod)
colnames(p_ssp370.x6.data_long)
p_ssp370.x7.data_dated<-p_ssp370.x6.data_long %>% 
  dplyr::mutate(t1 = substr(p_ssp370.x6.data_long$timeperiod, 43, 49)) %>%   #read char 12-18
  dplyr::mutate(t2= str_replace(t1, pattern = "_", replacement = "-01-"))  %>%  #find+replace
  dplyr::mutate(climatology = "ssp370")

colnames(p_ssp370.x7.data_dated)
unique(p_ssp370.x7.data_dated$climatology)
unique(p_ssp370.x7.data_dated$t1)
unique(p_ssp370.x7.data_dated$t2)





#turn the date column into date "format" 
p_ssp370.x7.data_dated$t3 <- as.POSIXct(p_ssp370.x7.data_dated$t2, format="%m-%d-%Y")
class(p_ssp370.x7.data_dated$t3)
range(p_ssp370.x7.data_dated$t3)
p_ssp370.x7.data_dated$t4<-as.Date(p_ssp370.x7.data_dated$t3,format="%m-%d-%Y")
class(p_ssp370.x7.data_dated$t4)

unique(p_ssp370.x7.data_dated$t3)
unique(p_ssp370.x7.data_dated$t4)

colnames(p_ssp370.x7.data_dated)
sum(is.na(p_ssp370.x7.data_dated$Climate_TAS))
sum(is.na(p_ssp370.x7.data_dated$t4))
sum(is.na(p_ssp370.x7.data_dated$longitude))   
sum(is.na(p_ssp370.x7.data_dated$timeperiod)) 





# Convert df to spatial for spatial-join - assign CoorRefSystm (CRS). 
class(p_ssp370.x7.data_dated)  # .csv's dataclass is not an sf_object yet - fix this.
st_crs(p_ssp370.x7.data_dated)
colnames(p_ssp370.x7.data_dated)
class(c_NHD_MTN_unique_sp)
st_crs(c_NHD_MTN_unique_sp)

p_ssp370.x8.SF <- st_as_sf(p_ssp370.x7.data_dated, coords = c("longitude", "latitude"), crs=4269,remove = FALSE)
class(p_ssp370.x8.SF) #check it's an sf_object 
st_crs(p_ssp370.x8.SF)

#remove unnecessary T1-3 columns
p_ssp370.x9.climate_by_date<-p_ssp370.x8.SF %>% 
  dplyr::select(-c(t1,t2,t3))







# .--- Join NHD_MTN + Climate --------------------------------------------------
class(p_ssp370.x9.climate_by_date)
class(c_NHD_MTN_unique)
class(c_NHD_MTN_unique_sp)

p_NHDssp370<-sf::st_join(x = p_ssp370.x9.climate_by_date,
                         y = c_NHD_MTN_unique,
                         left=FALSE) #lefttrue = keep all. leftfalse = keep where x links to y. 

#Explore 
sum(is.na(p_NHDssp370$NA_L3NA))
sum(is.na(p_NHDssp370$Climate_TAS))
sum(is.na(p_NHDssp370$t4))
sum(is.na(p_NHDssp370$COMID))
sum(is.na(c_NHD_MTN_unique$COMID)) 







# .--- Create dates column(s) --------------------------------------------------
# Make column that is just a month and a column that is just the year. 

p_NHD_Climate_MonthYear<-p_NHDssp370 %>% 
  dplyr::mutate( Year       = substr(p_NHDssp370$timeperiod, 51, 54)) %>%  #read char 51-54
  dplyr::mutate( Month_num  = substr(p_NHDssp370$timeperiod, 43, 44)) %>%  
  dplyr::mutate( Month_char = case_when(Month_num=="01" ~ "JAN", 
                                        Month_num=="02" ~ "FEB", 
                                        Month_num=="03" ~ "MAR", 
                                        Month_num=="04" ~ "APR", 
                                        Month_num=="05" ~ "MAY", 
                                        Month_num=="06" ~ "JUN", 
                                        Month_num=="07" ~ "JUL", 
                                        Month_num=="08" ~ "AUG", 
                                        Month_num=="09" ~ "SEP", 
                                        Month_num=="10" ~ "OCT", 
                                        Month_num=="11" ~ "NOV", 
                                        Month_num=="12" ~ "DEC")) %>% 
  dplyr::mutate(Temp_C = (Climate_TAS-273.15)) %>% 
  dplyr::mutate(Temp_F = (((Climate_TAS-273.15) * (9/5)) + 32))

unique(p_NHD_Climate_MonthYear$Year)

fBasics::basicStats(p_NHD_Climate_MonthYear$Climate_TAS)
fBasics::basicStats(p_NHD_Climate_MonthYear$Temp_C)
fBasics::basicStats(p_NHD_Climate_MonthYear$Temp_F)

class(p_NHD_Climate_MonthYear$Temp_C)
class(p_NHD_Climate_MonthYear$Temp_F)
class(p_NHD_Climate_MonthYear$Month_num) # Character - convert this to - Numeric 
class(p_NHD_Climate_MonthYear$Year)      # Character - convert this to - Numeric 
class(p_NHD_Climate_MonthYear$COMID)     # Character - Good.

p_NHD_Climate_MonthYear$Month_num<-as.numeric(p_NHD_Climate_MonthYear$Month_num)
class(p_NHD_Climate_MonthYear$Month_num) # Month_num is numeric

p_NHD_Climate_MonthYear$Year<-as.numeric(p_NHD_Climate_MonthYear$Year)
class(p_NHD_Climate_MonthYear$Year) # Year is numeric

p_NHD_Climate_MonthYear$COMID<-as.character(p_NHD_Climate_MonthYear$COMID)
class(p_NHD_Climate_MonthYear$COMID) # COMID is character







# Clear RAM  --------------------------------------------------------------

#sapply(ls(), function(x) object.size(get(x))) %>% sort %>% tail(20)
rm(p_NHDssp370)
rm(p_ssp370.x8.SF)
rm(p_ssp370.x7.data_dated)
rm(p_ssp370.x6.data_long)
rm(p_ssp370.x5)
rm(p_ssp370.x4)
rm(p_ssp370.x3)
rm(p_ssp370.x2)
rm(p_ssp370.x1)
rm(m_ROCKIES_ALL_coef)
rm(m_CASCADES_ALL_coef)
rm(m_APPALACHIAN_ALL_coef)
rm(m_SIERRA_coef)
rm(m_KLAMATH_coef)
rm(m_IDAHO_coef)
rm(m_BLUEMTN_coef)
rm(m_BLUERIDGE_coef)
rm(m_AZNM_coef)
rm(m_WASATCH_coef)
rm(m_ROCKIES_ALL_YearCombos)
rm(m_CASCADES_ALL_YearCombos)
rm(m_APPALACHIAN_ALL_YearCombos)
rm(m_SIERRA_YearCombos)
rm(m_IDAHO_YearCombos)
rm(m_BLUEMTN_YearCombos)
rm(m_BLUERIDGE_YearCombos)
rm(m_AZNM_YearCombos)
rm(m_KLAMATH_YearCombos)
rm(m_WASATCH_YearCombos)
#gc() # Free Unused Memory






# .--- Calc GDD & KDD -----------------------------------------------------

# Plot: temp histogram  
p_NHD_Climate_MonthYear %>% 
  ggplot()+
  geom_histogram(aes(x=Temp_C),bins=170)+
  #scale_x_continuous(breaks=seq(-30,35,2))+
  #scale_y_continuous(breaks=seq(0,700000,100000), labels = scales::comma)+
  theme_bw()

# Join historic quantiles to the projected rows
p_NHDClim_MnthYr_withQuantiles <- left_join(p_NHD_Climate_MonthYear, m_historic_quantiles, by = join_by(MtnRange_SIMPLE))

# Check behavior
#   x <- p_NHDClim_MnthYr_withQuantiles %>%   
#     dplyr::filter(COMID=="3800747" | COMID == "14981164" | COMID == "4710204") %>% 
#     dplyr::select(COMID, MtnRange_SIMPLE, Year, Month_num, Temp_C, quant25, quant70, quant90) %>% 
#     arrange(COMID, Year, Month_num)
  
q_kddgdd_criteria <- p_NHDClim_MnthYr_withQuantiles %>%  
  ungroup() %>% 
  dplyr::group_by( COMID, Year, Month_char) %>%
  dplyr::mutate(#gdd_neg5_threshold     = ((((min(Temp_C))+(max(Temp_C)))/2) - (-5))      ,
                #gdd_0_threshold        = ((((min(Temp_C))+(max(Temp_C)))/2) - (0))       ,
                #gdd_5_threshold        = ((((min(Temp_C))+(max(Temp_C)))/2) - (5))       ,
                #gdd_10_threshold       = ((((min(Temp_C))+(max(Temp_C)))/2) - (10))      ,
                gdd_0_threshold        = (((((Temp_C))) - (0)))       ,
                
                #kdd_25Perc_threshold   = ((max(Temp_C)) - (m_historic_quantiles$quant25)) ,
                #kdd_50Perc_threshold   = ((max(Temp_C)) - (m_historic_quantiles$quant50)) ,
                #kdd_70Perc_threshold   = ((max(Temp_C)) - (m_historic_quantiles$quant70)) ,
                #kdd_75Perc_threshold   = ((max(Temp_C)) - (m_historic_quantiles$quant75)) ,
                #kdd_80Perc_threshold   = ((max(Temp_C)) - (m_historic_quantiles$quant80)) ,
                #kdd_85Perc_threshold   = ((max(Temp_C)) - (m_historic_quantiles$quant85)) ,
                kdd_90Perc_threshold   = (((Temp_C)) - (quant90)) ,
                kdd_95Perc_threshold   = (((Temp_C)) - (quant95)) )%>% 
  dplyr::ungroup()

## Check that each Lake-Year-Month has unique threshhold
#        x<-q_kddgdd_criteria %>% 
#          group_by(COMID,Year,Month_char,kdd_95Perc_threshold) %>%
#          tally() %>% 
#          ungroup()

x <- q_kddgdd_criteria %>%   
  dplyr::filter(COMID=="3800747" | COMID == "14981164" | COMID == "4710204") %>% 
  dplyr::select(COMID, MtnRange_SIMPLE, Year, Month_num, quant25, quant90, gdd_0_threshold, kdd_90Perc_threshold ) %>% 
  arrange(COMID, Year, Month_num)

# To explore projected temperature quantiles under spp370:
x <- p_NHDClim_MnthYr_withQuantiles %>%  
  group_by( Year) %>% 
  dplyr::summarise(quant25 = quantile(Temp_C, probs = q[1]),
                   quant50 = quantile(Temp_C, probs = q[2]),
                   quant70 = quantile(Temp_C, probs = q[3]),
                   quant75 = quantile(Temp_C, probs = q[4]),
                   quant80 = quantile(Temp_C, probs = q[5]),
                   quant85 = quantile(Temp_C, probs = q[6]),
                   quant90 = quantile(Temp_C, probs = q[7]),
                   quant95 = quantile(Temp_C, probs = q[8])) %>% 
  ungroup()

fBasics::basicStats(q_kddgdd_criteria$gdd_0_threshold)
fBasics::basicStats(q_kddgdd_criteria$kdd_90Perc_threshold)
fBasics::basicStats(q_kddgdd_criteria$kdd_95Perc_threshold)

x_count_negative <-sum(q_kddgdd_criteria$gdd_0_threshold <0)  # Checking zeros pre-step
x_count_negative                                              # Checking zeros pre-step
x_count_zero <-sum(q_kddgdd_criteria$gdd_0_threshold == 0)    # Checking zeros pre-step
x_count_zero                                                  # Checking zeros pre-step
x_count_negative+x_count_zero                                 # Checking zeros pre-step

q_kddgdd_criteria_negs_are_zero <- q_kddgdd_criteria %>% 
  mutate(across(c(#gdd_neg5_threshold,
                  gdd_0_threshold,        
                  #gdd_5_threshold,       
                  #gdd_10_threshold,      
                  #kdd_25Perc_threshold,
                  #kdd_50Perc_threshold,
                  #kdd_70Perc_threshold,
                  #kdd_75Perc_threshold,
                  #kdd_80Perc_threshold,
                  #kdd_85Perc_threshold,
                  kdd_90Perc_threshold,
                  kdd_95Perc_threshold), ~ifelse( . < 0, 0, . )))

x_count_zero <-sum(q_kddgdd_criteria_negs_are_zero$gdd_0_threshold == 0) # Checking zeros post-step
x_count_zero                                                             # Checking zeros post-step

fBasics::basicStats(q_kddgdd_criteria_negs_are_zero$gdd_0_threshold)
fBasics::basicStats(q_kddgdd_criteria_negs_are_zero$kdd_90Perc_threshold)
fBasics::basicStats(q_kddgdd_criteria_negs_are_zero$kdd_95Perc_threshold)







# .--- GDD/KDD expansion --------------------------------------------------
# Next, expand months so that each Lake-Year-Month has days data too; fill days with data from first day.
## Because the projected data are not continuous-years data, this would require (group_by(COMID,Year)), but that is too computationally intensive. Instead, splitting df into 3 (one per year range), then joining back afterward. 


# Split into 3 dataframes; 1 per future range
q_kddgdd_2040 <- q_kddgdd_criteria_negs_are_zero %>% 
  mutate(Date = ymd(t4)) %>%
  dplyr::filter(Year == 2040) # 2011-2040
q_kddgdd_2070 <- q_kddgdd_criteria_negs_are_zero %>% 
  mutate(Date = ymd(t4)) %>% 
  dplyr::filter(Year == 2070) # 2041-2070
q_kddgdd_2100 <- q_kddgdd_criteria_negs_are_zero %>% 
  mutate(Date = ymd(t4)) %>% 
  dplyr::filter(Year == 2100) #2071-2100

unique(q_kddgdd_2040$Year)
unique(q_kddgdd_2070$Year)
unique(q_kddgdd_2100$Year)







# Get Sum GDD & KDD for 2011-2040

q_expansionprep1_2040 <- q_kddgdd_2040 %>% 
  mutate(Date = ymd(t4))

q_expansionprep2_2040 <- q_expansionprep1_2040 %>%
  mutate(DaysInMonth = lubridate::days_in_month(make_date(Year, Month_num, 1)))

q_KDDGDDsummarized_2040 <- q_expansionprep2_2040 %>% 
  ungroup() %>% 
  dplyr::mutate(#Month_GDD_neg5     = (DaysInMonth * gdd_neg5_threshold    ),
                Month_GDD_0        = (DaysInMonth * gdd_0_threshold       ),
                #Month_GDD_5        = (DaysInMonth * gdd_5_threshold       ),
                #Month_GDD_10       = (DaysInMonth * gdd_10_threshold      ),
                #Month_KDD_25Perc   = (DaysInMonth * kdd_25Perc_threshold  ),
                #Month_KDD_50Perc   = (DaysInMonth * kdd_50Perc_threshold  ),
                #Month_KDD_70Perc   = (DaysInMonth * kdd_70Perc_threshold  ),
                #Month_KDD_75Perc   = (DaysInMonth * kdd_75Perc_threshold  ),
                #Month_KDD_80Perc   = (DaysInMonth * kdd_80Perc_threshold  ),
                #Month_KDD_85Perc   = (DaysInMonth * kdd_85Perc_threshold  ),
                Month_KDD_90Perc   = (DaysInMonth * kdd_90Perc_threshold  ),
                Month_KDD_95Perc   = (DaysInMonth * kdd_95Perc_threshold  )) %>% 
  group_by(COMID, Year) %>% 
  dplyr::mutate(#Sum_GDD_neg5     = sum(Month_GDD_neg5    ),
    Sum_GDD_0         = sum(Month_GDD_0       ),
    #Sum_GDD_5        = sum(Month_GDD_5       ),
    #Sum_GDD_10       = sum(Month_GDD_10      ),
    #Sum_KDD_25Perc   = sum(Month_KDD_25Perc  ),
    #Sum_KDD_50Perc   = sum(Month_KDD_50Perc  ),
    #Sum_KDD_70Perc   = sum(Month_KDD_70Perc  ),
    #Sum_KDD_75Perc   = sum(Month_KDD_75Perc  ),
    #Sum_KDD_80Perc   = sum(Month_KDD_80Perc  ),
    #Sum_KDD_85Perc   = sum(Month_KDD_85Perc  ),
    Sum_KDD_90Perc    = sum(Month_KDD_90Perc  ),
    Sum_KDD_95Perc    = sum(Month_KDD_95Perc  ),
    TempC_Annual_Mean      = mean(Temp_C)     , 
    TempC_Annual_Min       = min(Temp_C)      ,
    TempC_Annual_Max       = max(Temp_C)      ,
    TempC_Annual_SD        = sd(Temp_C)        ,
    TempC_Annual_CV        = ((TempC_Annual_SD/TempC_Annual_Mean) * 100 )   ) %>% 
  ungroup() %>% 
  dplyr::mutate(Log1p_SumGDD = log1p(Sum_GDD_0),
                Log1p_KDD90  = log1p(Sum_KDD_90Perc) ) %>% 
  dplyr::distinct(COMID, Year, .keep_all = TRUE) %>% # No more Year-Month uniques; 1 365-day sum per year.
  dplyr::select(-c(Date, timeperiod, t4, Month_num, Month_char,
                   Climate_TAS, Temp_C, Temp_F, 
                   #gdd_neg5_threshold,
                   gdd_0_threshold,
                   #gdd_5_threshold,
                   #gdd_10_threshold,
                   #kdd_25Perc_threshold,
                   #kdd_50Perc_threshold,
                   #kdd_70Perc_threshold,
                   #kdd_75Perc_threshold,
                   #kdd_80Perc_threshold,
                   #kdd_85Perc_threshold,
                   kdd_90Perc_threshold,
                   kdd_95Perc_threshold,
                   #Month_GDD_neg5 ,
                   Month_GDD_0     ,
                   #Month_GDD_5     ,
                   #Month_GDD_10    ,
                   #Month_KDD_25Perc,
                   #Month_KDD_50Perc,
                   #Month_KDD_70Perc,
                   #Month_KDD_75Perc,
                   #Month_KDD_80Perc,
                   #Month_KDD_85Perc,
                   Month_KDD_90Perc,
                   Month_KDD_95Perc )) #removing columns based on the year-month (not applicable now)




# Get Sum GDD & KDD for 2040-2070

q_expansionprep1_2070 <- q_kddgdd_2070 %>% 
  mutate(Date = ymd(t4))

q_expansionprep2_2070 <- q_expansionprep1_2070 %>%
  mutate(DaysInMonth = lubridate::days_in_month(make_date(Year, Month_num, 1)))

q_KDDGDDsummarized_2070 <- q_expansionprep2_2070 %>% 
  ungroup() %>% 
  dplyr::mutate(#Month_GDD_neg5     = (DaysInMonth * gdd_neg5_threshold    ),
                Month_GDD_0        = (DaysInMonth * gdd_0_threshold       ),
                #Month_GDD_5        = (DaysInMonth * gdd_5_threshold       ),
                #Month_GDD_10       = (DaysInMonth * gdd_10_threshold      ),
                #Month_KDD_25Perc   = (DaysInMonth * kdd_25Perc_threshold  ),
                #Month_KDD_50Perc   = (DaysInMonth * kdd_50Perc_threshold  ),
                #Month_KDD_70Perc   = (DaysInMonth * kdd_70Perc_threshold  ),
                #Month_KDD_75Perc   = (DaysInMonth * kdd_75Perc_threshold  ),
                #Month_KDD_80Perc   = (DaysInMonth * kdd_80Perc_threshold  ),
                #Month_KDD_85Perc   = (DaysInMonth * kdd_85Perc_threshold  ),
                Month_KDD_90Perc   = (DaysInMonth * kdd_90Perc_threshold  ),
                Month_KDD_95Perc   = (DaysInMonth * kdd_95Perc_threshold  )) %>% 
  group_by(COMID, Year) %>% 
  dplyr::mutate(#Sum_GDD_neg5     = sum(Month_GDD_neg5    ),
    Sum_GDD_0         = sum(Month_GDD_0       ),
    #Sum_GDD_5        = sum(Month_GDD_5       ),
    #Sum_GDD_10       = sum(Month_GDD_10      ),
    #Sum_KDD_25Perc   = sum(Month_KDD_25Perc  ),
    #Sum_KDD_50Perc   = sum(Month_KDD_50Perc  ),
    #Sum_KDD_70Perc   = sum(Month_KDD_70Perc  ),
    #Sum_KDD_75Perc   = sum(Month_KDD_75Perc  ),
    #Sum_KDD_80Perc   = sum(Month_KDD_80Perc  ),
    #Sum_KDD_85Perc   = sum(Month_KDD_85Perc  ),
    Sum_KDD_90Perc    = sum(Month_KDD_90Perc  ),
    Sum_KDD_95Perc    = sum(Month_KDD_95Perc  ),
    TempC_Annual_Mean      = mean(Temp_C)     , 
    TempC_Annual_Min       = min(Temp_C)      ,
    TempC_Annual_Max       = max(Temp_C)      ,
    TempC_Annual_SD        = sd(Temp_C)        ,
    TempC_Annual_CV        = ((TempC_Annual_SD/TempC_Annual_Mean) * 100 )   ) %>% 
  ungroup() %>% 
  dplyr::mutate(Log1p_SumGDD = log1p(Sum_GDD_0),
                Log1p_KDD90  = log1p(Sum_KDD_90Perc) ) %>%
  dplyr::distinct(COMID, Year, .keep_all = TRUE) %>% # No more Year-Month uniques; 1 365-day sum per year.
  dplyr::select(-c(Date, timeperiod, t4, Month_num, Month_char,
                   Climate_TAS, Temp_C, Temp_F, 
                   #gdd_neg5_threshold,
                   gdd_0_threshold,
                   #gdd_5_threshold,
                   #gdd_10_threshold,
                   #kdd_25Perc_threshold,
                   #kdd_50Perc_threshold,
                   #kdd_70Perc_threshold,
                   #kdd_75Perc_threshold,
                   #kdd_80Perc_threshold,
                   #kdd_85Perc_threshold,
                   kdd_90Perc_threshold,
                   kdd_95Perc_threshold,
                   #Month_GDD_neg5 ,
                   Month_GDD_0     ,
                   #Month_GDD_5     ,
                   #Month_GDD_10    ,
                   #Month_KDD_25Perc,
                   #Month_KDD_50Perc,
                   #Month_KDD_70Perc,
                   #Month_KDD_75Perc,
                   #Month_KDD_80Perc,
                   #Month_KDD_85Perc,
                   Month_KDD_90Perc,
                   Month_KDD_95Perc)) #removing columns based on the year-month (not applicable now)




# Get Sum GDD & KDD for 2070-2100

q_expansionprep1_2100 <- q_kddgdd_2100 %>% 
  mutate(Date = ymd(t4))

q_expansionprep2_2100 <- q_expansionprep1_2100 %>%
  mutate(DaysInMonth = lubridate::days_in_month(make_date(Year, Month_num, 1)))

q_KDDGDDsummarized_2100 <- q_expansionprep2_2100 %>% 
  ungroup() %>% 
  dplyr::mutate(#Month_GDD_neg5     = (DaysInMonth * gdd_neg5_threshold    ),
                Month_GDD_0        = (DaysInMonth * gdd_0_threshold       ),
                #Month_GDD_5        = (DaysInMonth * gdd_5_threshold       ),
                #Month_GDD_10       = (DaysInMonth * gdd_10_threshold      ),
                #Month_KDD_25Perc   = (DaysInMonth * kdd_25Perc_threshold  ),
                #Month_KDD_50Perc   = (DaysInMonth * kdd_50Perc_threshold  ),
                #Month_KDD_70Perc   = (DaysInMonth * kdd_70Perc_threshold  ),
                #Month_KDD_75Perc   = (DaysInMonth * kdd_75Perc_threshold  ),
                #Month_KDD_80Perc   = (DaysInMonth * kdd_80Perc_threshold  ),
                #Month_KDD_85Perc   = (DaysInMonth * kdd_85Perc_threshold  ),
                Month_KDD_90Perc   = (DaysInMonth * kdd_90Perc_threshold  ),
                Month_KDD_95Perc   = (DaysInMonth * kdd_95Perc_threshold  )) %>% 
  group_by(COMID, Year) %>% 
  dplyr::mutate(#Sum_GDD_neg5     = sum(Month_GDD_neg5    ),
    Sum_GDD_0         = sum(Month_GDD_0       ),
    #Sum_GDD_5        = sum(Month_GDD_5       ),
    #Sum_GDD_10       = sum(Month_GDD_10      ),
    #Sum_KDD_25Perc   = sum(Month_KDD_25Perc  ),
    #Sum_KDD_50Perc   = sum(Month_KDD_50Perc  ),
    #Sum_KDD_70Perc   = sum(Month_KDD_70Perc  ),
    #Sum_KDD_75Perc   = sum(Month_KDD_75Perc  ),
    #Sum_KDD_80Perc   = sum(Month_KDD_80Perc  ),
    #Sum_KDD_85Perc   = sum(Month_KDD_85Perc  ),
    Sum_KDD_90Perc    = sum(Month_KDD_90Perc  ),
    Sum_KDD_95Perc    = sum(Month_KDD_95Perc  ),
    TempC_Annual_Mean      = mean(Temp_C)     , 
    TempC_Annual_Min       = min(Temp_C)      ,
    TempC_Annual_Max       = max(Temp_C)      ,
    TempC_Annual_SD        = sd(Temp_C)        ,
    TempC_Annual_CV        = ((TempC_Annual_SD/TempC_Annual_Mean) * 100 )   ) %>% 
  ungroup() %>% 
  dplyr::mutate(Log1p_SumGDD = log1p(Sum_GDD_0),
                Log1p_KDD90  = log1p(Sum_KDD_90Perc) ) %>% 
  dplyr::distinct(COMID, Year, .keep_all = TRUE) %>% # No more Year-Month uniques; 1 365-day sum per year.
  dplyr::select(-c(Date, timeperiod, t4, Month_num, Month_char,
                   Climate_TAS, Temp_C, Temp_F, 
                   #gdd_neg5_threshold,
                   gdd_0_threshold,
                   #gdd_5_threshold,
                   #gdd_10_threshold,
                   #kdd_25Perc_threshold,
                   #kdd_50Perc_threshold,
                   #kdd_70Perc_threshold,
                   #kdd_75Perc_threshold,
                   #kdd_80Perc_threshold,
                   #kdd_85Perc_threshold,
                   kdd_90Perc_threshold,
                   kdd_95Perc_threshold,
                   #Month_GDD_neg5 ,
                   Month_GDD_0     ,
                   #Month_GDD_5     ,
                   #Month_GDD_10    ,
                   #Month_KDD_25Perc,
                   #Month_KDD_50Perc,
                   #Month_KDD_70Perc,
                   #Month_KDD_75Perc,
                   #Month_KDD_80Perc,
                   #Month_KDD_85Perc,
                   Month_KDD_90Perc,
                   Month_KDD_95Perc )) #removing columns based on the year-month (not applicable now)



# Rejoin these 3 dfs
q_kddgddSUM_forssp <- rbind(q_KDDGDDsummarized_2040,
                            q_KDDGDDsummarized_2070,
                            q_KDDGDDsummarized_2100) # = 3 rows per lake (1 for each year)


x_a<-q_KDDGDDsummarized_2100 %>% dplyr::filter(COMID=="3800743") # Check this is the same format as above. 
x_b<-q_kddgddSUM_forssp      %>% dplyr::filter(COMID=="3800743") # Check this is the same format as above. 







# .---* Florida plot ----------------------------------------------------------

x_plot_kddelev<-q_kddgddSUM_forssp %>% 
  #dplyr::filter(!Sum_KDD_90Perc==0) %>% 
  ggplot()+
  geom_point(aes(x=((Sum_KDD_90Perc)), y=elevatn, color=factor(Year)), alpha=0.2)+
  scale_colour_manual(values = c("slategray4", "slategray3","slategray2" ))+
  facet_wrap(~Year ~ MtnRange_SIMPLE, ncol = 10)+
  scale_x_continuous(labels = scales::comma_format(),breaks = seq(0,2400,800), limits=c(0,2400))+
  scale_y_continuous(labels = scales::comma_format(),breaks = seq(0,4100,1000), limits=c(0,4100))+
  theme_bw()+
  xlab(label = "Killing Degree Days")+
  ylab(label = "Elevation (m)")+
  theme(strip.background = element_rect(fill=NA, color=NA),
        legend.position  = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title = element_text(size=15))

x_plot_kddelev <- q_kddgddSUM_forssp %>%
  #dplyr::filter(!Sum_KDD_90Perc == 0) %>%
  ggplot() +
  geom_point(aes(x = Sum_KDD_90Perc, y = elevatn, color = factor(Year)), alpha = 0.2) +
  scale_colour_manual(values = c("slategray4", "slategray3", "slategray2")) +
  facet_grid(rows = vars(Year), 
             cols = vars(MtnRange_SIMPLE), 
             labeller = labeller(Year = label_value, MtnRange_SIMPLE = label_value)) +
  scale_x_continuous(labels = scales::comma_format(), breaks = seq(0, 2400, 800), limits = c(0, 2400)) +
  scale_y_continuous(labels = scales::comma_format(), breaks = seq(0, 4100, 1000), limits = c(0, 4100)) +
  theme_bw() +
  xlab(label = "Killing Degree Days") +
  ylab(label = "Elevation (m)") +
  theme(
    strip.background = element_rect(fill = NA, color = NA),
    strip.placement = "outside",
    strip.text.y = element_text(size = 15, face = "bold"), # Styling Year
    strip.text.x = element_text(size = 10), # Styling MtnRange_SIMPLE
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 15))




x_plot_kddelev <- q_kddgddSUM_forssp %>%
  #dplyr::filter(!Sum_KDD_90Perc == 0) %>%
  ggplot() +
  geom_point(aes(x = Sum_KDD_90Perc, y = elevatn, color = factor(Year)), alpha = 0.2) +
  scale_colour_manual(values = c("slategray4", "slategray3", "slategray2")) +
  facet_grid(
    rows = vars(Year),
    cols = vars(MtnRange_SIMPLE),
    labeller = labeller(MtnRange_SIMPLE = function(x) stringr::str_wrap(x, width = 10)) # Adjust the width as needed
  ) +
  scale_x_continuous(labels = scales::comma_format(), breaks = seq(0, 2400, 800), limits = c(0, 2400)) +
  scale_y_continuous(labels = scales::comma_format(), breaks = seq(0, 4100, 1000), limits = c(0, 4100)) +
  theme_bw() +
  xlab(label = "Killing Degree Days") +
  ylab(label = "Elevation (m)") +
  theme(
    strip.background = element_rect(fill = NA, color = NA),
    strip.placement = "outside",
    strip.text.y = element_text(size = 15, face = "bold"), # Styling Year
    strip.text.x = element_text(size = 12), # Styling wrapped MtnRange_SIMPLE
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 15))

# Fig x -------------------------------------------------------------
ggsave(plot = x_plot_kddelev, "figures/Fig S3 now includes KDD zeros - 06.17.2025 - Elevation ~ SumKDD90 PROJECTED - SumKDD90==0 is removed --log1p.jpeg",height = 5, width = 13, dpi = 300)






# .--- Join spp + historic ------------------------------------------------
# Will perform the LMER for each range in order to get updated slopes for the DFA. Start by adding historic timeseries to timeperiod 1 (LMER T1), timerperiod 1+2 (LMER T2), and timeperiod 1+2+3 (LMER T3).

# Historic Data 
colnames(m_kddgdd_summarized)
unique(m_kddgdd_summarized$Year)
q_historic_Sum_GDD_0 <- m_kddgdd_summarized
  
# Projected Data

# 2040
q_2040_ssp_Sum_GDD_0 <- q_kddgddSUM_forssp %>%
  dplyr::filter(Year == 2040)
unique(q_2040_ssp_Sum_GDD_0$Year)

# 2070
q_2070_2040_ssp_Sum_GDD_0 <- q_kddgddSUM_forssp %>%
  dplyr::filter(!Year == 2100)
unique(q_2070_2040_ssp_Sum_GDD_0$Year)

# 2100
q_2100_2070_2040_ssp_Sum_GDD_0 <- q_kddgddSUM_forssp 
unique(q_2100_2070_2040_ssp_Sum_GDD_0$Year)


# Join the spp-timeperiods to the historic ones
q_for_LMER_2040           <- bind_rows(q_2040_ssp_Sum_GDD_0              , q_historic_Sum_GDD_0)
q_for_LMER_2070_2040      <- bind_rows(q_2070_2040_ssp_Sum_GDD_0         , q_historic_Sum_GDD_0)
q_for_LMER_2100_2070_2040 <- bind_rows(q_2100_2070_2040_ssp_Sum_GDD_0    , q_historic_Sum_GDD_0)

unique(q_for_LMER_2040$Year)
unique(q_for_LMER_2070_2040$Year)
unique(q_for_LMER_2100_2070_2040$Year)




# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
# . LMER T1 [2011-2040] ---------------------------------------------------

# Each range gets its own DF

q_KLAMATH_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Klamath Mountains")

q_SIERRA_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")

q_CASCADES_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Cascades") 

q_ROCKIES_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Rockies")

q_APPS_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Appalachians")

q_IDAHO_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Idaho Batholith")

q_BLUERIDGE_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Ridge")

q_BLUEMTN_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Mountains")

q_AZNM_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "AZ-NM Mountains")

q_WASATCH_2040<-q_for_LMER_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Wasatch-Uinta Mountains")

class(q_for_LMER_2040$Year)
class(q_for_LMER_2040$COMID)
class(q_for_LMER_2040$Log1p_SumGDD)
unique(q_KLAMATH_2040$Year)



# ~ KLAMATH ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model
class(q_KLAMATH_2040$COMID)

q_KLAMATH_2040_RslopeRintercept<-lmerTest::lmer(data = q_KLAMATH_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_KLAMATH_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_KLAMATH_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_KLAMATH_2040_RslopeRintercept)$COMID #View
coef( q_KLAMATH_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_KLAMATH_2040_coef <- coef(q_KLAMATH_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_KLAMATH_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_KLAMATH_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_KLAMATH_2040$COMID)), 
                                    Year = unique(q_KLAMATH_2040$Year))

unique(q_KLAMATH_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_KLAMATH_2040_YearSlopes <- left_join(q_KLAMATH_2040_YearCombos, q_KLAMATH_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_KLAMATH_2040_ymxb <-q_KLAMATH_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_KLAMATH_ssp2011.2040<-q_KLAMATH_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Klamath, R2c = 0.2189906 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ SIERRA ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_SIERRA_2040_RslopeRintercept<-lmerTest::lmer(data = q_SIERRA_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_SIERRA_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_SIERRA_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_SIERRA_2040_RslopeRintercept)$COMID #View
coef( q_SIERRA_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_SIERRA_2040_coef <- coef(q_SIERRA_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_SIERRA_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_SIERRA_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_SIERRA_2040$COMID)), 
                                       Year = unique(q_SIERRA_2040$Year))

unique(q_SIERRA_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_SIERRA_2040_YearSlopes <- left_join(q_SIERRA_2040_YearCombos, q_SIERRA_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_SIERRA_2040_ymxb <-q_SIERRA_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_SIERRA_ssp2011.2040<-q_SIERRA_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray2")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Sierra, R2c = 0.2312897 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ CASCADES ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_CASCADES_2040_RslopeRintercept<-lmerTest::lmer(data = q_CASCADES_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_CASCADES_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_CASCADES_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_CASCADES_2040_RslopeRintercept)$COMID #View
coef( q_CASCADES_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_CASCADES_2040_coef <- coef(q_CASCADES_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_CASCADES_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_CASCADES_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_CASCADES_2040$COMID)), 
                                       Year = unique(q_CASCADES_2040$Year))

unique(q_CASCADES_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_CASCADES_2040_YearSlopes <- left_join(q_CASCADES_2040_YearCombos, q_CASCADES_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_CASCADES_2040_ymxb <-q_CASCADES_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_CASCADES_ssp2011.2040<-q_CASCADES_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray2")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Cascades, R2c = 0.219897 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")




# ~ ROCKIES ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_ROCKIES_2040_RslopeRintercept<-lmerTest::lmer(data = q_ROCKIES_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_ROCKIES_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_ROCKIES_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_ROCKIES_2040_RslopeRintercept)$COMID #View
coef( q_ROCKIES_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_ROCKIES_2040_coef <- coef(q_ROCKIES_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_ROCKIES_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_ROCKIES_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_ROCKIES_2040$COMID)), 
                                       Year = unique(q_ROCKIES_2040$Year))

unique(q_ROCKIES_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_ROCKIES_2040_YearSlopes <- left_join(q_ROCKIES_2040_YearCombos, q_ROCKIES_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_ROCKIES_2040_ymxb <-q_ROCKIES_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_ROCKIES_ssp2011.2040<-q_ROCKIES_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray2")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Rockies, R2c = 0.2210353 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")




# ~ APPALACHIAN ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_APPS_2040_RslopeRintercept<-lmerTest::lmer(data = q_APPS_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_APPS_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_APPS_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_APPS_2040_RslopeRintercept)$COMID #View
coef( q_APPS_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_APPS_2040_coef <- coef(q_APPS_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_APPS_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_APPS_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_APPS_2040$COMID)), 
                                       Year = unique(q_APPS_2040$Year))

unique(q_APPS_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_APPS_2040_YearSlopes <- left_join(q_APPS_2040_YearCombos, q_APPS_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_APPS_2040_ymxb <-q_APPS_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_APPALACHIAN_ssp2011.2040<-q_APPS_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray2")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Appalachians, R2c = 0.2368787 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")




# ~ IDAHO ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_IDAHO_2040_RslopeRintercept<-lmerTest::lmer(data = q_IDAHO_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_IDAHO_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_IDAHO_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_IDAHO_2040_RslopeRintercept)$COMID #View
coef( q_IDAHO_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_IDAHO_2040_coef <- coef(q_IDAHO_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_IDAHO_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_IDAHO_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_IDAHO_2040$COMID)), 
                                       Year = unique(q_IDAHO_2040$Year))

unique(q_IDAHO_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_IDAHO_2040_YearSlopes <- left_join(q_IDAHO_2040_YearCombos, q_IDAHO_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_IDAHO_2040_ymxb <-q_IDAHO_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_IDAHO_ssp2011.2040<-q_IDAHO_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray2")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Idaho, R2c = 0.2501359 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")




# ~ BLUERIDGE ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_BLUERIDGE_2040_RslopeRintercept<-lmerTest::lmer(data = q_BLUERIDGE_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_BLUERIDGE_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_BLUERIDGE_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_BLUERIDGE_2040_RslopeRintercept)$COMID #View
coef( q_BLUERIDGE_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_BLUERIDGE_2040_coef <- coef(q_BLUERIDGE_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_BLUERIDGE_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_BLUERIDGE_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_BLUERIDGE_2040$COMID)), 
                                       Year = unique(q_BLUERIDGE_2040$Year))

unique(q_BLUERIDGE_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_BLUERIDGE_2040_YearSlopes <- left_join(q_BLUERIDGE_2040_YearCombos, q_BLUERIDGE_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_BLUERIDGE_2040_ymxb <-q_BLUERIDGE_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_BLUERIDGE_ssp2011.2040<-q_BLUERIDGE_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray2")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("BlueRidge, R2c = 0.2569465 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")




# ~ BLUEMTN ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_BLUEMTN_2040_RslopeRintercept<-lmerTest::lmer(data = q_BLUEMTN_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_BLUEMTN_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_BLUEMTN_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_BLUEMTN_2040_RslopeRintercept)$COMID #View
coef( q_BLUEMTN_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_BLUEMTN_2040_coef <- coef(q_BLUEMTN_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_BLUEMTN_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_BLUEMTN_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_BLUEMTN_2040$COMID)), 
                                       Year = unique(q_BLUEMTN_2040$Year))

unique(q_BLUEMTN_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_BLUEMTN_2040_YearSlopes <- left_join(q_BLUEMTN_2040_YearCombos, q_BLUEMTN_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_BLUEMTN_2040_ymxb <-q_BLUEMTN_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_BLUEMTN_ssp2011.2040<-q_BLUEMTN_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray2")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Blue Mountains, R2c = 0.2220566 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")




# ~ AZNM ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_AZNM_2040_RslopeRintercept<-lmerTest::lmer(data = q_AZNM_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_AZNM_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_AZNM_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_AZNM_2040_RslopeRintercept)$COMID #View
coef( q_AZNM_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_AZNM_2040_coef <- coef(q_AZNM_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_AZNM_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_AZNM_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_AZNM_2040$COMID)), 
                                       Year = unique(q_AZNM_2040$Year))

unique(q_AZNM_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_AZNM_2040_YearSlopes <- left_join(q_AZNM_2040_YearCombos, q_AZNM_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_AZNM_2040_ymxb <-q_AZNM_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_AZNM_ssp2011.2040<-q_AZNM_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray2")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("AZ-NM Mountains, R2c = 0.2515449 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")



# ~ WASATCH ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_WASATCH_2040_RslopeRintercept<-lmerTest::lmer(data = q_WASATCH_2040, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_WASATCH_2040_RslopeRintercept)
MuMIn::r.squaredGLMM(q_WASATCH_2040_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_WASATCH_2040_RslopeRintercept)$COMID #View
coef( q_WASATCH_2040_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_WASATCH_2040_coef <- coef(q_WASATCH_2040_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_WASATCH_2040$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_WASATCH_2040_YearCombos <- expand.grid(COMID = as.character(unique(q_WASATCH_2040$COMID)), 
                                       Year = unique(q_WASATCH_2040$Year))

unique(q_WASATCH_2040_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_WASATCH_2040_YearSlopes <- left_join(q_WASATCH_2040_YearCombos, q_WASATCH_2040_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_WASATCH_2040_ymxb <-q_WASATCH_2040_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_WASATCH_ssp2011.2040<-q_WASATCH_2040_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray2")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Wasatch, R2c = 0.2313163 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
# . LMER T2 [2041-2070][2011-2040] ---------------------------------------------------

# .--- Each range gets its own DF

q_KLAMATHssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Klamath Mountains")

q_SIERRAssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")

q_CASCADESssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Cascades") 

q_ROCKIESssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Rockies")

q_APPALACHIANssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Appalachians")

q_IDAHOssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Idaho Batholith")

q_BLUERIDGEssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Ridge")

q_BLUEMTNssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Mountains")

q_AZNMssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "AZ-NM Mountains")

q_WASATCHssp.2040.2070<-q_for_LMER_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Wasatch-Uinta Mountains")



# ~ KLAMATH ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_KLAMATHssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_KLAMATHssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_KLAMATHssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_KLAMATHssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_KLAMATHssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_KLAMATHssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_KLAMATHssp.2040.2070_coef <- coef(q_KLAMATHssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_KLAMATHssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_KLAMATHssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_KLAMATHssp.2040.2070$COMID)), 
                                       Year = unique(q_KLAMATHssp.2040.2070$Year))

unique(q_KLAMATHssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_KLAMATHssp.2040.2070_YearSlopes <- left_join(q_KLAMATHssp.2040.2070_YearCombos, q_KLAMATHssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_KLAMATHssp.2040.2070_ymxb <-q_KLAMATHssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_KLAMATH_2040.2070<-q_KLAMATHssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Klamath, R2c = 0.4754713 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")




# ~ SIERRA ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_SIERRAssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_SIERRAssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_SIERRAssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_SIERRAssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_SIERRAssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_SIERRAssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_SIERRAssp.2040.2070_coef <- coef(q_SIERRAssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_SIERRAssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_SIERRAssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_SIERRAssp.2040.2070$COMID)), 
                                                 Year = unique(q_SIERRAssp.2040.2070$Year))

unique(q_SIERRAssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_SIERRAssp.2040.2070_YearSlopes <- left_join(q_SIERRAssp.2040.2070_YearCombos, q_SIERRAssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_SIERRAssp.2040.2070_ymxb <-q_SIERRAssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_SIERRA_2040.2070<-q_SIERRAssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Sierra, R2c = 0.4599227 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")

# ~ CASCADES ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_CASCADESssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_CASCADESssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_CASCADESssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_CASCADESssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_CASCADESssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_CASCADESssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_CASCADESssp.2040.2070_coef <- coef(q_CASCADESssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_CASCADESssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_CASCADESssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_CASCADESssp.2040.2070$COMID)), 
                                                 Year = unique(q_CASCADESssp.2040.2070$Year))

unique(q_CASCADESssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_CASCADESssp.2040.2070_YearSlopes <- left_join(q_CASCADESssp.2040.2070_YearCombos, q_CASCADESssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_CASCADESssp.2040.2070_ymxb <-q_CASCADESssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_CASCADES_2040.2070<-q_CASCADESssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Cascades, R2c = 0.4830665 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")

# ~ ROCKIES ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_ROCKIESssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_ROCKIESssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_ROCKIESssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_ROCKIESssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_ROCKIESssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_ROCKIESssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_ROCKIESssp.2040.2070_coef <- coef(q_ROCKIESssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_ROCKIESssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_ROCKIESssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_ROCKIESssp.2040.2070$COMID)), 
                                                 Year = unique(q_ROCKIESssp.2040.2070$Year))

unique(q_ROCKIESssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_ROCKIESssp.2040.2070_YearSlopes <- left_join(q_ROCKIESssp.2040.2070_YearCombos, q_ROCKIESssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_ROCKIESssp.2040.2070_ymxb <-q_ROCKIESssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_ROCKIES_2040.2070<-q_ROCKIESssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Rockies, R2c = 0.5388623 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")

# ~ APPALACHIAN ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_APPALACHIANssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_APPALACHIANssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_APPALACHIANssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_APPALACHIANssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_APPALACHIANssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_APPALACHIANssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_APPALACHIANssp.2040.2070_coef <- coef(q_APPALACHIANssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_APPALACHIANssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_APPALACHIANssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_APPALACHIANssp.2040.2070$COMID)), 
                                                 Year = unique(q_APPALACHIANssp.2040.2070$Year))

unique(q_APPALACHIANssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_APPALACHIANssp.2040.2070_YearSlopes <- left_join(q_APPALACHIANssp.2040.2070_YearCombos, q_APPALACHIANssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_APPALACHIANssp.2040.2070_ymxb <-q_APPALACHIANssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_APPALACHIAN_2040.2070<-q_APPALACHIANssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Appalachians, R2c = 0.5064469 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")

# ~ IDAHO ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_IDAHOssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_IDAHOssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_IDAHOssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_IDAHOssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_IDAHOssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_IDAHOssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_IDAHOssp.2040.2070_coef <- coef(q_IDAHOssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_IDAHOssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_IDAHOssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_IDAHOssp.2040.2070$COMID)), 
                                                 Year = unique(q_IDAHOssp.2040.2070$Year))

unique(q_IDAHOssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_IDAHOssp.2040.2070_YearSlopes <- left_join(q_IDAHOssp.2040.2070_YearCombos, q_IDAHOssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_IDAHOssp.2040.2070_ymxb <-q_IDAHOssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_IDAHO_2040.2070<-q_IDAHOssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Idaho Mountains, R2c = 0.5214504 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")

# ~ BLUERIDGE ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_BLUERIDGEssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_BLUERIDGEssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_BLUERIDGEssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_BLUERIDGEssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_BLUERIDGEssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_BLUERIDGEssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_BLUERIDGEssp.2040.2070_coef <- coef(q_BLUERIDGEssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_BLUERIDGEssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_BLUERIDGEssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_BLUERIDGEssp.2040.2070$COMID)), 
                                                 Year = unique(q_BLUERIDGEssp.2040.2070$Year))

unique(q_BLUERIDGEssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_BLUERIDGEssp.2040.2070_YearSlopes <- left_join(q_BLUERIDGEssp.2040.2070_YearCombos, q_BLUERIDGEssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_BLUERIDGEssp.2040.2070_ymxb <-q_BLUERIDGEssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_BLUERIDGE_2040.2070<-q_BLUERIDGEssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Blueridge, R2c = 0.5365829 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")

# ~ BLUEMTN ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_BLUEMTNssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_BLUEMTNssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_BLUEMTNssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_BLUEMTNssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_BLUEMTNssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_BLUEMTNssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_BLUEMTNssp.2040.2070_coef <- coef(q_BLUEMTNssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_BLUEMTNssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_BLUEMTNssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_BLUEMTNssp.2040.2070$COMID)), 
                                                 Year = unique(q_BLUEMTNssp.2040.2070$Year))

unique(q_BLUEMTNssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_BLUEMTNssp.2040.2070_YearSlopes <- left_join(q_BLUEMTNssp.2040.2070_YearCombos, q_BLUEMTNssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_BLUEMTNssp.2040.2070_ymxb <-q_BLUEMTNssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_BLUEMTN_2040.2070<-q_BLUEMTNssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Blue Mountains, R2c = 0.4855026 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")

# ~ AZNM ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_AZNMssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_AZNMssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_AZNMssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_AZNMssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_AZNMssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_AZNMssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_AZNMssp.2040.2070_coef <- coef(q_AZNMssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_AZNMssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_AZNMssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_AZNMssp.2040.2070$COMID)), 
                                                 Year = unique(q_AZNMssp.2040.2070$Year))

unique(q_AZNMssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_AZNMssp.2040.2070_YearSlopes <- left_join(q_AZNMssp.2040.2070_YearCombos, q_AZNMssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_AZNMssp.2040.2070_ymxb <-q_AZNMssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_AZNM_2040.2070<-q_AZNMssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("AZ-NM Mountains, R2c = 0.5183487 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")

# ~ WASATCH ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_WASATCHssp.2040.2070_RslopeRintercept<-lmerTest::lmer(data = q_WASATCHssp.2040.2070, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_WASATCHssp.2040.2070_RslopeRintercept)
MuMIn::r.squaredGLMM(q_WASATCHssp.2040.2070_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_WASATCHssp.2040.2070_RslopeRintercept)$COMID #View
coef( q_WASATCHssp.2040.2070_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_WASATCHssp.2040.2070_coef <- coef(q_WASATCHssp.2040.2070_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_WASATCHssp.2040.2070$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_WASATCHssp.2040.2070_YearCombos <- expand.grid(COMID = as.character(unique(q_WASATCHssp.2040.2070$COMID)), 
                                                 Year = unique(q_WASATCHssp.2040.2070$Year))

unique(q_WASATCHssp.2040.2070_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_WASATCHssp.2040.2070_YearSlopes <- left_join(q_WASATCHssp.2040.2070_YearCombos, q_WASATCHssp.2040.2070_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_WASATCHssp.2040.2070_ymxb <-q_WASATCHssp.2040.2070_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_WASATCH_2040.2070<-q_WASATCHssp.2040.2070_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray3")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Wasatch, R2c = 0.4974892 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")



# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
# . LMER T3 [2071-2100][2041-2070][2011-2040] ---------------------------------------------------

q_KLAMATHssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Klamath Mountains")

q_SIERRAssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")

q_CASCADESssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Cascades") 

q_ROCKIESssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Rockies")

q_APPALACHIANssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Appalachians")

q_IDAHOssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Idaho Batholith")

q_BLUERIDGEssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Ridge")

q_BLUEMTNssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Blue Mountains")

q_AZNMssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "AZ-NM Mountains")

q_WASATCHssp.AllYears<-q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Wasatch-Uinta Mountains")

q_for_LMER_2100_2070_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "AZ-NM Mountains") %>% 
  ggplot()+
  geom_point(aes(x=Year, y=Log1p_SumGDD), alpha=0.3)+
  facet_wrap(.~MtnRange_SIMPLE,ncol=5)




# ~ KLAMATH ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_KLAMATHssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_KLAMATHssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_KLAMATHssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_KLAMATHssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_KLAMATHssp.AllYears_RslopeRintercept)$COMID #View
coef( q_KLAMATHssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_KLAMATHssp.AllYears_coef <- coef(q_KLAMATHssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_KLAMATHssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_KLAMATHssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_KLAMATHssp.AllYears$COMID)), 
                                                 Year = unique(q_KLAMATHssp.AllYears$Year))

unique(q_KLAMATHssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_KLAMATHssp.AllYears_YearSlopes <- left_join(q_KLAMATHssp.AllYears_YearCombos, q_KLAMATHssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_KLAMATHssp.AllYears_ymxb <-q_KLAMATHssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_KLAMATH_AllYears<-q_KLAMATHssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Klamath, R2c = 0.6377102 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")



# ~ SIERRA ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_SIERRAssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_SIERRAssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_SIERRAssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_SIERRAssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_SIERRAssp.AllYears_RslopeRintercept)$COMID #View
coef( q_SIERRAssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_SIERRAssp.AllYears_coef <- coef(q_SIERRAssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_SIERRAssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_SIERRAssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_SIERRAssp.AllYears$COMID)), 
                                                Year = unique(q_SIERRAssp.AllYears$Year))

unique(q_SIERRAssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_SIERRAssp.AllYears_YearSlopes <- left_join(q_SIERRAssp.AllYears_YearCombos, q_SIERRAssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_SIERRAssp.AllYears_ymxb <-q_SIERRAssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_SIERRA_AllYears<-q_SIERRAssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Sierra, R2c = 0.6377891 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ CASCADES ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_CASCADESssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_CASCADESssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_CASCADESssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_CASCADESssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_CASCADESssp.AllYears_RslopeRintercept)$COMID #View
coef( q_CASCADESssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_CASCADESssp.AllYears_coef <- coef(q_CASCADESssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_CASCADESssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_CASCADESssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_CASCADESssp.AllYears$COMID)), 
                                                Year = unique(q_CASCADESssp.AllYears$Year))

unique(q_CASCADESssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_CASCADESssp.AllYears_YearSlopes <- left_join(q_CASCADESssp.AllYears_YearCombos, q_CASCADESssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_CASCADESssp.AllYears_ymxb <-q_CASCADESssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_CASCADES_AllYears<-q_CASCADESssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Cascades, R2c = 0.6450512 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ ROCKIES ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_ROCKIESssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_ROCKIESssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_ROCKIESssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_ROCKIESssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_ROCKIESssp.AllYears_RslopeRintercept)$COMID #View
coef( q_ROCKIESssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_ROCKIESssp.AllYears_coef <- coef(q_ROCKIESssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_ROCKIESssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_ROCKIESssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_ROCKIESssp.AllYears$COMID)), 
                                                Year = unique(q_ROCKIESssp.AllYears$Year))

unique(q_ROCKIESssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_ROCKIESssp.AllYears_YearSlopes <- left_join(q_ROCKIESssp.AllYears_YearCombos, q_ROCKIESssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_ROCKIESssp.AllYears_ymxb <-q_ROCKIESssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_ROCKIES_AllYears<-q_ROCKIESssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Rockies, R2c = 0.7007889 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ APPALACHIAN ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_APPALACHIANssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_APPALACHIANssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_APPALACHIANssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_APPALACHIANssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_APPALACHIANssp.AllYears_RslopeRintercept)$COMID #View
coef( q_APPALACHIANssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_APPALACHIANssp.AllYears_coef <- coef(q_APPALACHIANssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_APPALACHIANssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_APPALACHIANssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_APPALACHIANssp.AllYears$COMID)), 
                                                Year = unique(q_APPALACHIANssp.AllYears$Year))

unique(q_APPALACHIANssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_APPALACHIANssp.AllYears_YearSlopes <- left_join(q_APPALACHIANssp.AllYears_YearCombos, q_APPALACHIANssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_APPALACHIANssp.AllYears_ymxb <-q_APPALACHIANssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_APPALACHIAN_AllYears<-q_APPALACHIANssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Appalachians, R2c = 0.6776855 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ IDAHO ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_IDAHOssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_IDAHOssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_IDAHOssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_IDAHOssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_IDAHOssp.AllYears_RslopeRintercept)$COMID #View
coef( q_IDAHOssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_IDAHOssp.AllYears_coef <- coef(q_IDAHOssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_IDAHOssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_IDAHOssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_IDAHOssp.AllYears$COMID)), 
                                                Year = unique(q_IDAHOssp.AllYears$Year))

unique(q_IDAHOssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_IDAHOssp.AllYears_YearSlopes <- left_join(q_IDAHOssp.AllYears_YearCombos, q_IDAHOssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_IDAHOssp.AllYears_ymxb <-q_IDAHOssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_IDAHO_AllYears<-q_IDAHOssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Idaho, R2c = 0.698067 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ BLUERIDGE ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_BLUERIDGEssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_BLUERIDGEssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_BLUERIDGEssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_BLUERIDGEssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_BLUERIDGEssp.AllYears_RslopeRintercept)$COMID #View
coef( q_BLUERIDGEssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_BLUERIDGEssp.AllYears_coef <- coef(q_BLUERIDGEssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_BLUERIDGEssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_BLUERIDGEssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_BLUERIDGEssp.AllYears$COMID)), 
                                                Year = unique(q_BLUERIDGEssp.AllYears$Year))

unique(q_BLUERIDGEssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_BLUERIDGEssp.AllYears_YearSlopes <- left_join(q_BLUERIDGEssp.AllYears_YearCombos, q_BLUERIDGEssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_BLUERIDGEssp.AllYears_ymxb <-q_BLUERIDGEssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_BLUERIDGE_AllYears<-q_BLUERIDGEssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Blueridge, R2c = 0.7049202 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ BLUEMTN ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_BLUEMTNssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_BLUEMTNssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_BLUEMTNssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_BLUEMTNssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_BLUEMTNssp.AllYears_RslopeRintercept)$COMID #View
coef( q_BLUEMTNssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_BLUEMTNssp.AllYears_coef <- coef(q_BLUEMTNssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_BLUEMTNssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_BLUEMTNssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_BLUEMTNssp.AllYears$COMID)), 
                                                Year = unique(q_BLUEMTNssp.AllYears$Year))

unique(q_BLUEMTNssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_BLUEMTNssp.AllYears_YearSlopes <- left_join(q_BLUEMTNssp.AllYears_YearCombos, q_BLUEMTNssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_BLUEMTNssp.AllYears_ymxb <-q_BLUEMTNssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_BLUEMTN_AllYears<-q_BLUEMTNssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Blue Mountains, R2c = 0.6406844 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ AZNM ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_AZNMssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_AZNMssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_AZNMssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_AZNMssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_AZNMssp.AllYears_RslopeRintercept)$COMID #View
coef( q_AZNMssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_AZNMssp.AllYears_coef <- coef(q_AZNMssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_AZNMssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_AZNMssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_AZNMssp.AllYears$COMID)), 
                                                Year = unique(q_AZNMssp.AllYears$Year))

unique(q_AZNMssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_AZNMssp.AllYears_YearSlopes <- left_join(q_AZNMssp.AllYears_YearCombos, q_AZNMssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_AZNMssp.AllYears_ymxb <-q_AZNMssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_AZNM_AllYears<-q_AZNMssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("AZ-NM Mountains, R2c = 0.6806758 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")


# ~ WASATCH ------------------------------------------------------------
# LMER - Random Slope, Random Intercept Model

q_WASATCHssp.AllYears_RslopeRintercept<-lmerTest::lmer(data = q_WASATCHssp.AllYears, Log1p_SumGDD ~ Year + (Year|COMID))
summary(q_WASATCHssp.AllYears_RslopeRintercept)
MuMIn::r.squaredGLMM(q_WASATCHssp.AllYears_RslopeRintercept) 

# Extract random effects from model (random slope, random intercept)
ranef(q_WASATCHssp.AllYears_RslopeRintercept)$COMID #View
coef( q_WASATCHssp.AllYears_RslopeRintercept)$COMID #View

# Get DF with 1 slope & intercept per lake
q_WASATCHssp.AllYears_coef <- coef(q_WASATCHssp.AllYears_RslopeRintercept)$COMID %>% 
  rownames_to_column() %>% 
  rename(COMID = rowname) %>% 
  rename(Intercept = "(Intercept)") %>% 
  rename(Slope = "Year") %>% 
  mutate(Year = NA)

# Create a sequence of years from 1979 to 2019
fBasics::basicStats(q_WASATCHssp.AllYears$Year)

# Create a data frame with all combinations of LAKE and YEAR
q_WASATCHssp.AllYears_YearCombos <- expand.grid(COMID = as.character(unique(q_WASATCHssp.AllYears$COMID)), 
                                                Year = unique(q_WASATCHssp.AllYears$Year))

unique(q_WASATCHssp.AllYears_YearCombos$Year)

# Perform a left join to retain the Random Intercept and Random Slope values
q_WASATCHssp.AllYears_YearSlopes <- left_join(q_WASATCHssp.AllYears_YearCombos, q_WASATCHssp.AllYears_coef, by = "COMID")

# Tidy up & "predict" trend over time -
q_WASATCHssp.AllYears_ymxb <-q_WASATCHssp.AllYears_YearSlopes %>% 
  rename(Year = "Year.x") %>% 
  dplyr::select(-c("Year.y")) %>% 
  mutate(SumGDD0_mxb = ((Slope * Year) + Intercept))  

# Plot
x_PLOT_WASATCH_AllYears<-q_WASATCHssp.AllYears_ymxb %>%
  ggplot(aes(x=Year, y=SumGDD0_mxb, group = factor(COMID)))+
  geom_line(linewidth=0.7,alpha=0.2, color = "slategray4")+ 
  #geom_point(alpha=0.1, size=3)+
  #geom_abline(aes(intercept = -357.244941, slope = 0.231152), linewidth=2, alpha=0.99, color="black")+
  theme_bw()+
    scale_x_continuous(limits = c(1980, 2020),  breaks = seq(1980,2100, 10))+
  scale_y_continuous(limits = c(5.65,8.7), breaks = seq(5.65,8.7, 0.5))+
  ggtitle("Wasatch, R2c = 0.6549196 ")+
  ylab(label = "Sum GDD 0 (modeled)")+
  theme(
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    legend.position= "none")




# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
# .--- plot LMER trends (3x) -----------------------------------------------------

#dev.off()
#  plots_gdd0_2040<-cowplot::plot_grid(x_PLOT_KLAMATH_ssp2011.2040,
#                                              x_PLOT_SIERRA_ssp2011.2040,
#                                              x_PLOT_CASCADES_ssp2011.2040,
#                                              x_PLOT_ROCKIES_ssp2011.2040,
#                                              x_PLOT_APPALACHIAN_ssp2011.2040,
#                                              x_PLOT_IDAHO_ssp2011.2040,
#                                              x_PLOT_BLUERIDGE_ssp2011.2040,
#                                              x_PLOT_BLUEMTN_ssp2011.2040,
#                                              x_PLOT_AZNM_ssp2011.2040,
#                                              x_PLOT_WASATCH_ssp2011.2040,
#                                              ncol=10)
#  
#  plots_gdd0_2040_2070<-cowplot::plot_grid(x_PLOT_KLAMATH_2040.2070,
#                                              x_PLOT_SIERRA_2040.2070,
#                                              x_PLOT_CASCADES_2040.2070,
#                                              x_PLOT_ROCKIES_2040.2070,
#                                              x_PLOT_APPALACHIAN_2040.2070,
#                                              x_PLOT_IDAHO_2040.2070,
#                                              x_PLOT_BLUERIDGE_2040.2070,
#                                              x_PLOT_BLUEMTN_2040.2070,
#                                              x_PLOT_AZNM_2040.2070,
#                                              x_PLOT_WASATCH_2040.2070,
#                                              ncol=10)
#  
#  plots_gdd0_2040_2070_2100<-cowplot::plot_grid(x_PLOT_KLAMATH_AllYears,
#                                           x_PLOT_SIERRA_AllYears,
#                                           x_PLOT_CASCADES_AllYears,
#                                           x_PLOT_ROCKIES_AllYears,
#                                           x_PLOT_APPALACHIAN_AllYears,
#                                           x_PLOT_IDAHO_AllYears,
#                                           x_PLOT_BLUERIDGE_AllYears,
#                                           x_PLOT_BLUEMTN_AllYears,
#                                           x_PLOT_AZNM_AllYears,
#                                           x_PLOT_WASATCH_AllYears,
#                                           ncol=10)
#  
#  plots_gdd0<-cowplot::plot_grid(xPLOT_gdd_KLAMATH,
#                                 xPLOT_gdd_SIERRA,
#                                 xPLOT_gdd_CASCADES,
#                                 xPLOT_gdd_ROCKIES,
#                                 xPLOT_gdd_APPS,
#                                 xPLOT_gdd_IDAHO,
#                                 xPLOT_gdd_BLUERIDGE,
#                                 xPLOT_gdd_BLUEMTN,
#                                 xPLOT_gdd_AZNM,
#                                 xPLOT_gdd_WASATCH,
#                                 ncol=10)
#  
#
#plot_gdd_projected_ALLtrends <- cowplot::plot_grid(plots_gdd0,
#                                                   plots_gdd0_2040,
#                                                   plots_gdd0_2040_2070,
#                                                   plots_gdd0_2040_2070_2100, nrow=4)
#
# ggsave(plot = plot_gdd_projected_ALLtrends, "figures/06.13.2024 - 6-03-2024--1980-1999-2018--GDD0_Model_AllYears.jpeg", height = 6, width = 18, dpi = 300)
# 
# Plotting + save takes 37 minutes (April 16, 2024)








# .--- Join the three timeperiods -----------------------------------------------
# RBind the model outputs. 

r_GDD0_Bound_2040 <-bind_rows(q_KLAMATH_2040_ymxb,
                                  q_SIERRA_2040_ymxb,
                                  q_CASCADES_2040_ymxb,
                                  q_ROCKIES_2040_ymxb,
                                  q_IDAHO_2040_ymxb,
                                  q_BLUEMTN_2040_ymxb,
                                  q_APPS_2040_ymxb,
                                  q_BLUERIDGE_2040_ymxb,
                                  q_AZNM_2040_ymxb,
                                  q_WASATCH_2040_ymxb)

r_GDD0_Bound_2070.2040 <-bind_rows(q_KLAMATHssp.2040.2070_ymxb,
                                   q_SIERRAssp.2040.2070_ymxb,
                                   q_CASCADESssp.2040.2070_ymxb,
                                   q_ROCKIESssp.2040.2070_ymxb,
                                   q_IDAHOssp.2040.2070_ymxb,
                                   q_BLUEMTNssp.2040.2070_ymxb,
                                   q_APPALACHIANssp.2040.2070_ymxb,
                                   q_BLUERIDGEssp.2040.2070_ymxb,
                                   q_AZNMssp.2040.2070_ymxb,
                                   q_WASATCHssp.2040.2070_ymxb)

r_GDD0_Bound_2100.2070.2040 <-bind_rows(q_KLAMATHssp.AllYears_ymxb,
                                  q_SIERRAssp.AllYears_ymxb,
                                  q_CASCADESssp.AllYears_ymxb,
                                  q_ROCKIESssp.AllYears_ymxb,
                                  q_IDAHOssp.AllYears_ymxb,
                                  q_BLUEMTNssp.AllYears_ymxb,
                                  q_APPALACHIANssp.AllYears_ymxb,
                                  q_BLUERIDGEssp.AllYears_ymxb,
                                  q_AZNMssp.AllYears_ymxb,
                                  q_WASATCHssp.AllYears_ymxb)

fBasics::basicStats(r_GDD0_Bound_2040$SumGDD0_mxb)
fBasics::basicStats(r_GDD0_Bound_2070.2040$SumGDD0_mxb)
fBasics::basicStats(r_GDD0_Bound_2100.2070.2040$SumGDD0_mxb)




# Put 3 timeperiods into 1 df

colnames(r_GDD0_Bound_2040)
colnames(r_GDD0_Bound_2070.2040)
colnames(r_GDD0_Bound_2100.2070.2040)




# To bring these 3 into one df; rename "Intercept" "Slope" and "SumGDD0_mxb" to be unique 

r_modeled_slopes_2040_retitled<-r_GDD0_Bound_2040 %>% 
  rename(Intercept_2040 = Intercept) %>% 
  rename(Slope_2040     = Slope) %>% 
  rename(SumGDD0_mxb_2040   = SumGDD0_mxb) 

r_modeled_slopes_2070.2040_retitled<-r_GDD0_Bound_2070.2040 %>% 
  rename(Intercept_2070_2040 = Intercept) %>% 
  rename(Slope_2070_2040     = Slope) %>% 
  rename(SumGDD0_mxb_2070_2040   = SumGDD0_mxb) 

r_modeled_slopes_2100_2070_2040_retitled<-r_GDD0_Bound_2100.2070.2040 %>% 
  rename(Intercept_2100_2070_2040 = Intercept) %>% 
  rename(Slope_2100_2070_2040     = Slope) %>% 
  rename(SumGDD0_mxb_2100_2070_2040   = SumGDD0_mxb) 

colnames(r_modeled_slopes_2040_retitled)
colnames(r_modeled_slopes_2070.2040_retitled)
colnames(r_modeled_slopes_2100_2070_2040_retitled)


unique(r_modeled_slopes_2040_retitled$Year)
r_modeled_slopes_join1 <- full_join(r_modeled_slopes_2040_retitled,
                                    r_modeled_slopes_2070.2040_retitled, by = join_by(COMID,Year))
colnames(r_modeled_slopes_join1)
r_modeled_slopes_join2 <- full_join(r_modeled_slopes_join1,
                                    r_modeled_slopes_2100_2070_2040_retitled,  by = join_by(COMID,Year))
colnames(r_modeled_slopes_join2)



# .--- Join modeled outputs + Mtn Attributes ----------------------------------

# Check modeled-slopes join result, w arbitrary COMID
x<-r_modeled_slopes_join2 %>% 
  dplyr::filter(COMID=="3800743") %>%  # FYI - Each Lake-Year is 6 rows (NOT 1) [Addressed below]
  arrange(COMID,Year) 
x<-c_NHD_MTN_unique %>% 
  dplyr::filter(COMID=="3800743")  




# Join this to original DF w lake attributes.
r_NHD_LMER_spatial_ssp370 <- full_join(c_NHD_MTN_unique, r_modeled_slopes_join2, by = join_by(COMID))
x <- r_NHD_LMER_spatial_ssp370 %>% dplyr::filter(COMID=="3800743") # Check, w arbitrary COMID

# Drop spatial awareness (again, becSlope_2040# Drop spatial awareness (again, because c_ still had it)
class(r_NHD_LMER_spatial_ssp370)
r_NHD_LMER_spatial_ssp370_NONSPATIAL <- st_drop_geometry(r_NHD_LMER_spatial_ssp370)
class(r_NHD_LMER_spatial_ssp370_NONSPATIAL)

r_MtnSlopeAttribs_ssp370 <- r_NHD_LMER_spatial_ssp370_NONSPATIAL # No apparent outliers.






# .--- Mean GDD and distinct COMID ---------------------------------------------

# Link "MEAN GDD" calc from pre-model data to Lake-Int-Slope
### Use the split up version of "q_kddgddSUM_forssp"; which was the df right before model. 

colnames(q_kddgddSUM_forssp)
unique(q_2040_ssp_Sum_GDD_0$Year)
unique(q_2070_2040_ssp_Sum_GDD_0$Year)
unique(q_2100_2070_2040_ssp_Sum_GDD_0$Year)

q_MeanGDD_for_2040 <- q_2040_ssp_Sum_GDD_0 %>% # "q_" was df that was used to split ranges pre-model
  group_by(COMID) %>% 
  mutate(Mean_GDD_2040 = mean(Sum_GDD_0),
         Mean_KDD_2040 = mean(Sum_KDD_90Perc),
         TempC_2040_DistinctMeanFor_Mean = mean(TempC_Annual_Mean),
         TempC_2040_DistinctMeanFor_Min  = mean(TempC_Annual_Min),
         TempC_2040_DistinctMeanFor_Max  = mean(TempC_Annual_Max),
         TempC_2040_DistinctMeanFor_SD   = mean(TempC_Annual_SD),
         TempC_2040_DistinctMeanFor_CV   = mean(TempC_Annual_CV)) %>% 
  ungroup() %>% 
  dplyr::distinct(COMID, .keep_all = TRUE) %>% # 1 row per lake
  dplyr::select(c(latitude, longitude, 
                  COMID, 
                  Mean_GDD_2040, Mean_KDD_2040,
                  TempC_2040_DistinctMeanFor_Mean,
                  TempC_2040_DistinctMeanFor_Min,
                  TempC_2040_DistinctMeanFor_Max,
                  TempC_2040_DistinctMeanFor_SD,
                  TempC_2040_DistinctMeanFor_CV)) 

q_MeanGDD_for_2070_2040 <- q_2070_2040_ssp_Sum_GDD_0 %>% # "q_" was df that was used to split ranges pre-model
  group_by(COMID) %>% 
  mutate(Mean_GDD_2070_2040 = mean(Sum_GDD_0),
         Mean_KDD_2070_2040 = mean(Sum_KDD_90Perc),
         TempC_2070_2040_DistinctMeanFor_Mean = mean(TempC_Annual_Mean),
         TempC_2070_2040_DistinctMeanFor_Min  = mean(TempC_Annual_Min),
         TempC_2070_2040_DistinctMeanFor_Max  = mean(TempC_Annual_Max),
         TempC_2070_2040_DistinctMeanFor_SD   = mean(TempC_Annual_SD),
         TempC_2070_2040_DistinctMeanFor_CV   = mean(TempC_Annual_CV)) %>% 
  ungroup() %>% 
  dplyr::distinct(COMID, .keep_all = TRUE) %>% # 1 row per lake
  dplyr::select(c(COMID, 
                  Mean_GDD_2070_2040, Mean_KDD_2070_2040,
                  TempC_2070_2040_DistinctMeanFor_Mean,
                  TempC_2070_2040_DistinctMeanFor_Min,
                  TempC_2070_2040_DistinctMeanFor_Max,
                  TempC_2070_2040_DistinctMeanFor_SD,
                  TempC_2070_2040_DistinctMeanFor_CV)) 

q_MeanGDD_for_2100_2070_2040 <- q_2100_2070_2040_ssp_Sum_GDD_0 %>% # "q_" was df that was used to split ranges pre-model
  group_by(COMID) %>% 
  mutate(Mean_GDD_2100_2070_2040 = mean(Sum_GDD_0),
         Mean_KDD_2100_2070_2040 = mean(Sum_KDD_90Perc),
         TempC_2100_2070_2040_DistinctMeanFor_Mean = mean(TempC_Annual_Mean),
         TempC_2100_2070_2040_DistinctMeanFor_Min  = mean(TempC_Annual_Min),
         TempC_2100_2070_2040_DistinctMeanFor_Max  = mean(TempC_Annual_Max),
         TempC_2100_2070_2040_DistinctMeanFor_SD   = mean(TempC_Annual_SD),
         TempC_2100_2070_2040_DistinctMeanFor_CV   = mean(TempC_Annual_CV)) %>% 
  ungroup() %>% 
  dplyr::distinct(COMID, .keep_all = TRUE) %>% # 1 row per lake
  dplyr::select(c(COMID, 
                  Mean_GDD_2100_2070_2040, Mean_KDD_2100_2070_2040,
                  TempC_2100_2070_2040_DistinctMeanFor_Mean,
                  TempC_2100_2070_2040_DistinctMeanFor_Min,
                  TempC_2100_2070_2040_DistinctMeanFor_Max,
                  TempC_2100_2070_2040_DistinctMeanFor_SD,
                  TempC_2100_2070_2040_DistinctMeanFor_CV)) 


r_Distinct_LakeSlopeInt_ssp370<-r_MtnSlopeAttribs_ssp370 %>% 
  group_by(COMID) %>% 
  dplyr::distinct(COMID, .keep_all = TRUE) %>% # 1 row per lake
  dplyr::select(-c(Year, 
                   SumGDD0_mxb_2040,
                   SumGDD0_mxb_2070_2040,
                   SumGDD0_mxb_2100_2070_2040)) %>% 
  ungroup() %>% 
  dplyr::filter(!is.na(Slope_2040))

x<-r_Distinct_LakeSlopeInt_ssp370 %>% dplyr::filter(COMID=="3800743") #Each lake as 1 unique slope/intercept

s_MeanGDD_with_LakeIntSlope_2040 <- left_join(r_Distinct_LakeSlopeInt_ssp370,
                                              q_MeanGDD_for_2040, 
                                              by=join_by(COMID))

s_MeanGDD_with_LakeIntSlope_2070_2040 <- left_join(s_MeanGDD_with_LakeIntSlope_2040,
                                                   q_MeanGDD_for_2070_2040, 
                                                   by=join_by(COMID))

s_MeanGDD_with_LakeIntSlope_2100_2070_2040 <- left_join(s_MeanGDD_with_LakeIntSlope_2070_2040,
                                                q_MeanGDD_for_2100_2070_2040, 
                                                by=join_by(COMID))






# .---* RE Slope plotted (3x) ------------------------------------------------------
# Plot - RE Slope x Elevation (slope from SumGDD0 model)

r_Plot_slopeelev_2040<-s_MeanGDD_with_LakeIntSlope_2100_2070_2040 %>%
  ggplot()+
  geom_point( aes(x=Slope_2040, y=elevatn), color="slategray2", alpha=0.1)+
  facet_wrap(~MtnRange_SIMPLE, ncol = 10)+
  #scale_y_continuous(breaks = seq(0, 4500, 1000), labels = scales::comma_format())+
  #scale_x_continuous(breaks = seq(5, 18.5,  10), limits = c(0,18.5))+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA, color=NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab(label = "Elevation (m)")

r_Plot_slopeelev_2070<-s_MeanGDD_with_LakeIntSlope_2100_2070_2040 %>%
  ggplot()+
  geom_point( aes(x=Slope_2070_2040, y=elevatn), color="slategray3", alpha=0.1)+
  facet_wrap(~MtnRange_SIMPLE, ncol = 10)+
  #scale_y_continuous(breaks = seq(0, 4500, 1000), labels = scales::comma_format())+
  #scale_x_continuous(breaks = seq(5, 18.5,  10), limits = c(0,18.5))+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA, color=NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab(label = "Elevation (m)")

r_Plot_slopeelev_2100<-s_MeanGDD_with_LakeIntSlope_2100_2070_2040 %>%
  ggplot()+
  geom_point( aes(x=Slope_2100_2070_2040, y=elevatn), color="slategray4", alpha=0.1)+
  facet_wrap(~MtnRange_SIMPLE, ncol = 10)+
  #scale_y_continuous(breaks = seq(0, 4500, 1000), labels = scales::comma_format())+
  #scale_x_continuous(breaks = seq(5, 18.5,  10), limits = c(0,18.5))+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA, color=NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab(label = "Elevation (m)")

plots_slopeelev_projected<-cowplot::plot_grid(r_Plot_slopeelev_2040,
                               r_Plot_slopeelev_2070,
                               r_Plot_slopeelev_2100, 
                               nrow=3)
# Fig x
ggsave(plot = plots_slopeelev_projected, "figures/12.19.2024 - GDD0 REslope ~ Elevation - ssp370.jpeg", height = 5, width = 16, dpi = 300)







# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
# * Cluster Analysis ------------------------------------------------------
# Note: Cluster analysis is for historical data only. 


# Pre-Cluster Checks
sum(is.na(n_Means_with_LakeIntSlope$Slope))
sum(is.na(n_Means_with_LakeIntSlope$Mean_GDD))
sum(is.na(n_Means_with_LakeIntSlope$elevatn)) 

sum(n_Means_with_LakeIntSlope$Slope == 0)
sum(n_Means_with_LakeIntSlope$Mean_GDD == 0)
sum(n_Means_with_LakeIntSlope$elevatn == 0)

class((n_Means_with_LakeIntSlope$Slope))
class((n_Means_with_LakeIntSlope$Mean_GDD))
class((n_Means_with_LakeIntSlope$elevatn))

# LOG the column(s) being clustered upon
o_Cluster_LogScale<-n_Means_with_LakeIntSlope %>% 
  dplyr::mutate(Log1p_Slope              = log1p(Slope),
                Log1p_MeanGDD            = log1p(Mean_GDD),
                Log1p_Elev               = log1p(elevatn),
                Log10_TempC_MeanDistinct = log10( 10 + TempC_DistinctMeanFor_Mean))

o_clustered_data_list <- list()# Initialize an empty list to store clustered data

set.seed(3)

for ( region in unique(o_Cluster_LogScale$MtnRange_SIMPLE)) {
  
  #Subset data for the current region & columns desired
  region_data<-o_Cluster_LogScale %>% 
    dplyr::filter(MtnRange_SIMPLE==region)
  
  # Columns to cluster upon
  cluster_columns<-region_data %>% 
    dplyr::select(Log1p_MeanGDD)
  
  #Perform kmeans clustering
  clusters <- kmeans(cluster_columns, centers=3, iter.max = 1000)
  
  #Store clusterd data
  region_data$cluster<-clusters$cluster
  
  o_clustered_data_list[[region]] <- region_data
}

#Combine clusterd data for all regions into one dataset
o_combined_data<-do.call(rbind,o_clustered_data_list)


# Re-color clusters bc that assignment was randomized
unique(o_combined_data$MtnRange_SIMPLE)

o_recol_KLAMATH  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("Klamath Mountains"))
o_recol_KLAMATH2 <- o_recol_KLAMATH %>% mutate(cluster_fix = case_when(cluster==1~3, cluster==2~1, cluster==3~2)) %>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))
o_recol_KLAMATH2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)

o_recol_SIERRA  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("Sierra Nevada"))
o_recol_SIERRA2 <- o_recol_SIERRA %>% mutate(cluster_fix = case_when(cluster==1~2, cluster==2~3, cluster==3~1))%>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))  
o_recol_SIERRA2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)


o_recol_CASCADES  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("Cascades"))
o_recol_CASCADES2 <- o_recol_CASCADES %>% mutate(cluster_fix = case_when(cluster==1~1, cluster==2~2, cluster==3~3))%>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))  
o_recol_CASCADES2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)


o_recol_ROCKIES  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("Rockies"))
o_recol_ROCKIES2 <- o_recol_ROCKIES %>% mutate(cluster_fix = case_when(cluster==1~3, cluster==2~1, cluster==3~2)) %>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))
o_recol_ROCKIES2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)


o_recol_IDAHO  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("Idaho Batholith"))
o_recol_IDAHO2 <- o_recol_IDAHO %>% mutate(cluster_fix = case_when(cluster==1~3, cluster==2~1, cluster==3~2)) %>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))
o_recol_IDAHO2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)


o_recol_BLUEMTN  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("Blue Mountains"))
o_recol_BLUEMTN2 <- o_recol_BLUEMTN %>% mutate(cluster_fix = case_when(cluster==1~3, cluster==2~2, cluster==3~1)) %>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))
o_recol_BLUEMTN2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)


o_recol_APPS  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("Appalachians"))
o_recol_APPS2 <- o_recol_APPS %>% mutate(cluster_fix = case_when(cluster==1~1, cluster==2~3, cluster==3~2)) %>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))
o_recol_APPS2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)

o_recol_BLUERIDGE  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("Blue Ridge"))
o_recol_BLUERIDGE2 <- o_recol_BLUERIDGE %>% mutate(cluster_fix = case_when(cluster==1~3, cluster==2~1, cluster==3~2)) %>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))
o_recol_BLUERIDGE2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)

o_recol_AZNM  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("AZ-NM Mountains"))
o_recol_AZNM2 <- o_recol_AZNM %>% mutate(cluster_fix = case_when(cluster==1~2, cluster==2~3, cluster==3~1)) %>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))
o_recol_AZNM2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)

o_recol_WU  <- o_combined_data %>% subset(MtnRange_SIMPLE %in% c("Wasatch-Uinta Mountains"))
o_recol_WU2 <- o_recol_WU %>% mutate(cluster_fix = case_when(cluster==1~1, cluster==2~3, cluster==3~2)) %>% 
  mutate(cluster_label = case_when( cluster_fix==1 ~ "1 Hot",
                                    cluster_fix==2 ~ "2 Transitional",
                                    cluster_fix==3 ~ "3 Cold"))
o_recol_WU2 %>% ggplot(aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+ geom_point(alpha=0.6,size=5)

o_recoded <- rbind(o_recol_KLAMATH2,
                   o_recol_SIERRA2,
                   o_recol_CASCADES2,
                   o_recol_ROCKIES2,
                   o_recol_APPS2,
                   o_recol_IDAHO2,
                   o_recol_BLUERIDGE2,
                   o_recol_BLUEMTN2,
                   o_recol_AZNM2,
                   o_recol_WU2)




# .-----Clusterplots --------------------------------------------------------------
# Plot with clusters

x1 <- ggplot(o_recoded, aes(x=Mean_GDD, y=Slope, color=factor(cluster_label)))+
  geom_point(alpha=0.2, size=3)+
  scale_x_continuous(breaks = seq(570, 7830, 2420), limits = c(570,7830), labels = scales::comma_format())+
  scale_y_continuous(breaks = seq(0.0019, 0.007, 0.0017), limits = c(0.0019,0.007))+
  scale_color_manual( values = c("coral2", "slategray", "navy"))+
  facet_wrap(~MtnRange_SIMPLE,ncol=5)+
  labs(color= "Thermal Class")+
  xlab(label= "Mean GDD")+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA, color=NA),
        strip.text   = element_text(size=20),
        axis.title.y = element_text(size=28),
        axis.text    = element_text(size=15),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "none") 
        #plot.margin = unit(c(0, 1.4, 0, 0),  "inches"))

x2 <- ggplot(o_recoded, aes(x=Mean_GDD, y=elevatn, color=factor(cluster_label)))+
  geom_point(alpha=0.2, size=3)+
  scale_x_continuous(breaks = seq(570, 7830, 2420), limits = c(570,7830), labels = scales::comma_format())+
  scale_y_continuous(breaks = seq(0, 4500,1000), labels = scales::comma_format())+
  scale_color_manual( values = c("coral2", "slategray", "navy"))+
  facet_wrap(~MtnRange_SIMPLE,ncol=5)+
  labs(color= "Thermal Class")+
  xlab(label="Mean GDD")+
  ylab(label="Elevation (m)")+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA, color=NA),
        strip.text   = element_text(size=20),
        axis.title.y = element_text(size=28),
        axis.text    = element_text(size=15),
        legend.text  = element_text(size=20),
        legend.title = element_text(size=22),
        axis.title.x = element_blank())
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),)


legend <- get_legend(x2)

x2 <- x2 +  theme(legend.position = "none")

x3 <- cowplot::plot_grid(x1,x2,ncol=1)

x4 <- plot_grid(x3, legend, ncol=2, rel_widths = c(2,0.26))

#y.grob <- grid::textGrob("Mean GDD", gp=gpar(fontsize=15), rot=90)
x.grob <- grid::textGrob("Mean GDD", gp=gpar(fontsize=24))
x5 <- gridExtra::grid.arrange( arrangeGrob(x4, 
                                          #left = y.grob, 
                                          bottom = x.grob))


# Fig 5 -------------------------------------------------------------------
ggsave(plot = x5, "figures/12.19.2024 -- Fig5 5 -- RCP7.0 Cluster on Log1p MeanGDD -notSCALE- -- Clusterplot MeanGDD ~ Slope or Elevation.jpeg", height = 12, width = 19.5, dpi = 300)




# .-----Boxplots --------------------------------------------------------------
# Boxplot comparing KDD/GDD to clusters

x_meangdd <- ggplot(o_recoded, aes(x=cluster_label, y=Mean_GDD, fill=factor(cluster_label)))+
  geom_boxplot(coef=0,outlier.shape=NA)+
  scale_y_continuous(breaks = seq(570, 7830, 2420), limits = c(570,7830), labels = scales::comma_format())+
  scale_fill_manual( values = c("coral2", "slategray", "navy"))+
  facet_wrap(~MtnRange_SIMPLE, ncol=10)+
  labs(fill= "Thermal Class")+
  ylab(label="Mean GDD")+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA, color=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        strip.text = element_text(size=16),
        axis.text.y = element_text(size=13),
        legend.position = "none")


x_meakdd <- o_recoded %>% 
  #dplyr::filter(!Mean_KDD>1) %>% 
  ggplot(aes(x=cluster_label, y=Mean_KDD, fill=factor(cluster_label)))+
  scale_fill_manual( values = c("coral2", "slategray", "navy"))+
  geom_boxplot(coef=0,outlier.shape=NA)+
  scale_y_continuous(breaks = seq(0, 1500, 375), limits = c(0,1500), labels = scales::comma_format())+
  facet_wrap(~MtnRange_SIMPLE,ncol=10)+
  labs(fill= "Thermal Class")+
  ylab(label="Mean KDD")+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA, color=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        strip.text = element_text(size=16),
        axis.text.y = element_text(size=13),
        legend.position = "none")

x_elevation <- ggplot(o_recoded, aes(x=cluster_label, y=elevatn, fill=factor(cluster_label)))+
  geom_boxplot(coef=0,outlier.shape=NA)+
  scale_y_continuous(labels = scales::comma_format())+
  scale_fill_manual( values = c("coral2", "slategray", "navy"))+
  facet_wrap(~MtnRange_SIMPLE, ncol=10)+
  labs(fill= "Thermal Class")+
  ylab(label="Elevation (m)")+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA, color=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        strip.text = element_text(size=16),
        axis.text.y = element_text(size=13),
        legend.position = "none")

x_slope <- ggplot(o_recoded, aes(x=cluster_label, y=Slope, fill=factor(cluster_label)))+
  geom_boxplot(coef=0,outlier.shape=NA)+
  scale_y_continuous(breaks = seq(0.0019, 0.007, 0.0017), limits = c(0.0019,0.007))+
  scale_fill_manual( values = c("coral2", "slategray", "navy"))+
  facet_wrap(~MtnRange_SIMPLE,ncol=10)+
  labs(fill= "Thermal Class")+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA, color=NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=25),
        strip.text = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_text(size=25))

legend <- get_legend(x_slope)

x_slope <- x_slope +  theme(legend.position = "none")

x_boxplots <- cowplot::plot_grid(x_meangdd,
                                 x_meakdd,
                                 x_elevation,
                                 x_slope,
                                 ncol=1)

x_finalbox <- plot_grid(x_boxplots, legend, ncol=2, rel_widths = c(4,0.4))


# Fig S5 ------------------------------------------------------------------
ggsave(plot = x_finalbox, "figures/12.19.2024 -- Fig S5 -- RCP7.0 Cluster on Log1p MeanGDD -notSCALE- -- BoxPlot Slope.jpeg", height = 0, width = 25, dpi = 300)




# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------

# * DFA ---------------------------------------------------------------------
# Use this object for projected: s_MeanGDD_with_LakeIntSlope_2100_2070_2040 
# Use this object for historic:  o_recoded

# Prep Projected Data (log + scale)
t_almost_test_data  <- s_MeanGDD_with_LakeIntSlope_2100_2070_2040 %>% 
  
  mutate(LOG1P_Slope_2040           = log1p(Slope_2040),
         LOG1P_Slope_2070_2040      = log1p(Slope_2070_2040),
         LOG1P_Slope_2100_2070_2040 = log1p(Slope_2100_2070_2040),
         
         LOG1P_MeanGDD_2040           = log1p(Mean_GDD_2040),
         LOG1P_MeanGDD_2070_2040      = log1p(Mean_GDD_2070_2040),
         LOG1P_MeanGDD_2100_2070_2040 = log1p(Mean_GDD_2100_2070_2040),
         
         LOG1P_TEMPC_DISTINCT_MEAN_2040           = log10( 10 + TempC_2040_DistinctMeanFor_Mean),
         LOG1P_TEMPC_DISTINCT_MEAN_2070_2040      = log10( 10 + TempC_2070_2040_DistinctMeanFor_Mean),
         LOG1P_TEMPC_DISTINCT_MEAN_2100_2070_2040 = log10( 10 + TempC_2100_2070_2040_DistinctMeanFor_Mean),
         
         LOG1P_TEMPC_DISTINCT_MIN_2040           = log10(20+TempC_2040_DistinctMeanFor_Min),
         LOG1P_TEMPC_DISTINCT_MIN_2070_2040      = log10(20+TempC_2070_2040_DistinctMeanFor_Min),
         LOG1P_TEMPC_DISTINCT_MIN_2100_2070_2040 = log10(20+TempC_2100_2070_2040_DistinctMeanFor_Min),
         
         LOG1P_TEMPC_DISTINCT_MAX_2040           = log10(20+TempC_2040_DistinctMeanFor_Max),
         LOG1P_TEMPC_DISTINCT_MAX_2070_2040      = log10(20+TempC_2070_2040_DistinctMeanFor_Max),
         LOG1P_TEMPC_DISTINCT_MAX_2100_2070_2040 = log10(20+TempC_2100_2070_2040_DistinctMeanFor_Max),
         
         LOG1P_TEMPC_DISTINCT_SD_2040            = log10(20+TempC_2040_DistinctMeanFor_SD),
         LOG1P_TEMPC_DISTINCT_SD_2070_2040       = log10(20+TempC_2070_2040_DistinctMeanFor_SD),
         LOG1P_TEMPC_DISTINCT_SD_2100_2070_2040  = log10(20+TempC_2100_2070_2040_DistinctMeanFor_SD))





t_almost_test_data$LOGSCALE_Slope_2040           <- scale(t_almost_test_data$LOG1P_Slope_2040)
t_almost_test_data$LOGSCALE_Slope_2070_2040      <- scale(t_almost_test_data$LOG1P_Slope_2070_2040)
t_almost_test_data$LOGSCALE_Slope_2100_2070_2040 <- scale(t_almost_test_data$LOG1P_Slope_2100_2070_2040)

t_almost_test_data$LOGSCALE_MeanGDD_2040           <- scale(t_almost_test_data$LOG1P_MeanGDD_2040)
t_almost_test_data$LOGSCALE_MeanGDD_2070_2040      <- scale(t_almost_test_data$LOG1P_MeanGDD_2070_2040)
t_almost_test_data$LOGSCALE_MeanGDD_2100_2070_2040 <- scale(t_almost_test_data$LOG1P_MeanGDD_2100_2070_2040)

t_almost_test_data$LOGSCALE_TEMPC_distinct_MEAN_2040           <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_MEAN_2040)
t_almost_test_data$LOGSCALE_TEMPC_distinct_MEAN_2070_2040      <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_MEAN_2070_2040)
t_almost_test_data$LOGSCALE_TEMPC_distinct_MEAN_2100_2070_2040 <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_MEAN_2100_2070_2040)

t_almost_test_data$LOGSCALE_TEMPC_distinct_MIN_2040            <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_MIN_2040)
t_almost_test_data$LOGSCALE_TEMPC_distinct_MIN_2070_2040       <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_MIN_2070_2040)
t_almost_test_data$LOGSCALE_TEMPC_distinct_MIN_2100_2070_2040  <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_MIN_2100_2070_2040)

t_almost_test_data$LOGSCALE_TEMPC_distinct_MAX_2040            <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_MAX_2040)
t_almost_test_data$LOGSCALE_TEMPC_distinct_MAX_2070_2040       <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_MAX_2070_2040)
t_almost_test_data$LOGSCALE_TEMPC_distinct_MAX_2100_2070_2040  <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_MAX_2100_2070_2040)

t_almost_test_data$LOGSCALE_TEMPC_distinct_SD_2040             <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_SD_2040)
t_almost_test_data$LOGSCALE_TEMPC_distinct_SD_2070_2040        <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_SD_2070_2040)
t_almost_test_data$LOGSCALE_TEMPC_distinct_SD_2100_2070_2040   <- scale(t_almost_test_data$LOG1P_TEMPC_DISTINCT_SD_2100_2070_2040)





# Finish prepping Projected Data (for PCA and then DFA)                                            
colnames(t_almost_test_data)

t_test_data_2040 <- t_almost_test_data %>% 
  dplyr::select(c(COMID, MtnRange_SIMPLE, 
                  
                  Slope_2040,
                  LOG1P_Slope_2040,
                  LOGSCALE_Slope_2040,
                  
                  Mean_GDD_2040,
                  LOG1P_MeanGDD_2040,
                  LOGSCALE_MeanGDD_2040,
                  
                  TempC_2040_DistinctMeanFor_Mean,
                  LOG1P_TEMPC_DISTINCT_MEAN_2040,
                  LOGSCALE_TEMPC_distinct_MEAN_2040 )) %>% 

  rename(LogScaleSLOPE            = LOGSCALE_Slope_2040,
         Log1p_Slope              = LOG1P_Slope_2040,
         Log10_TempC_MeanDistinct = LOG1P_TEMPC_DISTINCT_MEAN_2040,
         Log1p_MeanGDD            = LOG1P_MeanGDD_2040)




t_test_data_2070_2040 <- t_almost_test_data %>% 
  dplyr::select(c(COMID, MtnRange_SIMPLE, 
                  
                  Slope_2070_2040,
                  LOG1P_Slope_2070_2040,
                  LOGSCALE_Slope_2070_2040,
                  
                  Mean_GDD_2070_2040,
                  LOG1P_MeanGDD_2070_2040,
                  LOGSCALE_MeanGDD_2070_2040,
                  
                  TempC_2070_2040_DistinctMeanFor_Mean,
                  LOG1P_TEMPC_DISTINCT_MEAN_2070_2040,
                  LOGSCALE_TEMPC_distinct_MEAN_2070_2040 )) %>% 
  
  rename(LogScaleSLOPE            = LOGSCALE_Slope_2070_2040, 
         Log1p_Slope              = LOG1P_Slope_2070_2040,
         Log10_TempC_MeanDistinct = LOG1P_TEMPC_DISTINCT_MEAN_2070_2040,
         Log1p_MeanGDD            = LOG1P_MeanGDD_2070_2040)





t_test_data_2100_2070_2040 <- t_almost_test_data %>% 
  dplyr::select(c(COMID, MtnRange_SIMPLE, 
                  
                  Slope_2100_2070_2040,
                  LOG1P_Slope_2100_2070_2040,
                  LOGSCALE_Slope_2100_2070_2040,
                  
                  Mean_GDD_2100_2070_2040,
                  LOG1P_MeanGDD_2100_2070_2040,
                  LOGSCALE_MeanGDD_2100_2070_2040,
                  
                  TempC_2100_2070_2040_DistinctMeanFor_Mean,
                  LOG1P_TEMPC_DISTINCT_MEAN_2100_2070_2040,
                  LOGSCALE_TEMPC_distinct_MEAN_2100_2070_2040 )) %>% 
  rename(LogScaleSLOPE            = LOGSCALE_Slope_2100_2070_2040,
         Log1p_Slope              = LOG1P_Slope_2100_2070_2040,
         Log10_TempC_MeanDistinct = LOG1P_TEMPC_DISTINCT_MEAN_2100_2070_2040,
         Log1p_MeanGDD            = LOG1P_MeanGDD_2100_2070_2040)





# Prep Historic Data (no prep) for DFA
t_hist_for_ConMat_loop <- o_recoded %>% 
  dplyr::select(COMID, MtnRange_SIMPLE, cluster_label, Log1p_MeanGDD)
colnames(o_recoded)







# .--- ConfusionMatrix FOR LOOP ---------------------------------------
# To test model accuracy ON HISTORIC data

models                    <- list() # List to store models
test_predictions_list     <- list() # List to store predictions
accuracies                <- list() # List to store accuracies
conf_matrices             <- list() # List to store confusion matrices
combined_predictions_list <- list() # List to store all-regions predicted clusters 

all_regions <- unique(t_hist_for_ConMat_loop$MtnRange_SIMPLE)
all_regions

# Loop through each region
for (region in all_regions) {
  
  # Subset for the current region
  region_data <- subset(t_hist_for_ConMat_loop, MtnRange_SIMPLE == region)
  
  # Generate indices for training and testing data
  train_index <- sample(1:nrow(region_data), 0.8 * nrow(region_data)) # 80% training (0.8)
  
  # Subset train & test data for the current region
  train_region_data <- region_data[ train_index, ] 
  test_region_data  <- region_data[-train_index, ]
  
  # Perform LDA model for the current region
  model <- MASS::lda(cluster_label ~ Log1p_MeanGDD, data = train_region_data)
  
  # Store the model with a name indicating the region
  model_name <- paste("model_region", region, sep = "_")
  assign(model_name, model, envir = .GlobalEnv)
  
  # Store the model in the list
  models[[region]] <- model
  
  # Predict the cluster for test data of the current region
  test_predictions_list <- predict(model, newdata = test_region_data)
  
  # Store predictions with attributes
  test_predictions_df <- data.frame(COMID             = test_region_data$COMID, 
                                    Cluster_Predicted = test_predictions_list$class)
  
  #Save Result
  assign(paste0(region, "test_predictions_df"), test_predictions_df, envir = .GlobalEnv)
  
  # Store predictions in a list
  test_predictions_list[[region]] <- test_predictions_df
  
  # Model Accuracy for the current region
  test_region_data$cluster_label <- as.factor(test_region_data$cluster_label)
  
  accuracy <- mean(test_predictions_list$class == test_region_data$cluster_label)
  accuracies[[region]] <- accuracy
  
  # Confusion Matrix for the current region
  conf_matrix <- caret::confusionMatrix(test_predictions_list$class, test_region_data$cluster_label)
  conf_matrices[[region]] <- conf_matrix
  
}

models
accuracies # Historic
conf_matrices




# .--- ConfusionMatrix MANUAL ---------------------------------------------
colnames(o_recoded)
t_hist_for_ConMat_manual <- o_recoded

unique(t_hist_for_ConMat_manual$MtnRange_SIMPLE)

o_recol_SIERRA2 <- t_hist_for_ConMat_manual%>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")


p_train_index <- sample(1:nrow(o_recol_SIERRA2), 0.8 * nrow(o_recol_SIERRA2)) # 80% training dataset
t_SIERRA_train <- o_recol_SIERRA2[ p_train_index, ] # Training dataset
t_SIERRA_test  <- o_recol_SIERRA2[-p_train_index, ] # Testing dataset

# L.DFA model
t_LDAmodel <- MASS::lda(cluster_label ~ Log1p_MeanGDD, data = t_SIERRA_train)

# Summary Stats/Plots
summary(t_LDAmodel)
t_LDAmodel # Coefficients are here. 
plot(t_LDAmodel)


#Predicting the cluster for new data
t_test_predictions_sierra    <- predict(t_LDAmodel, newdata = t_SIERRA_test)
class(t_test_predictions_sierra)            #List
head(t_test_predictions_sierra$class,6)     #Explore
head(t_test_predictions_sierra$posterior,6) #Explore
head(t_test_predictions_sierra$x,3)         #Explore

t_test_predictions_sierra_df <- data.frame(COMID=t_SIERRA_test$COMID, Cluster_Predicted = t_test_predictions_sierra)

# This is just for plotting; not used for anything else.
x_ploting <- cbind(t_SIERRA_train, predict(t_LDAmodel)$x)
ggplot(x_ploting, aes(cluster_label, LD1))+
  geom_point(aes(color=cluster_label))

# Model Accuracy
t_SIERRA_test$cluster_label<-as.factor(t_SIERRA_test$cluster_label)

caret::confusionMatrix( t_test_predictions_sierra_df$Cluster_Predicted.class, t_SIERRA_test$cluster_label) # Accuracy is here.
mean(t_test_predictions_sierra$class==t_SIERRA_test$cluster_label) # Accuracy is here. 






# .--- Test projected for loop: SIERRA ------------------------------------
t_hist_forloop_manualtest <- o_recoded

class(o_recoded$cluster_label)
colnames(o_recoded)

colnames(t_test_data_2040)

unique(t_hist_forloop_manualtest$MtnRange_SIMPLE)

t_SIERRA_train <- t_hist_forloop_manualtest%>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")
t_SIERRA_test <- t_test_data_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Sierra Nevada")

# L.DFA model
t_LDAmodel <- MASS::lda(cluster_label ~ Log1p_MeanGDD, data = t_SIERRA_train)

# Summary Stats/Plots
summary(t_LDAmodel)
t_LDAmodel # Coefficients are here. 
plot(t_LDAmodel)


#Predicting the cluster for new data
t_test_predictions_2040    <- predict(t_LDAmodel, newdata = t_SIERRA_test)
class(t_test_predictions_2040)            #List
head(t_test_predictions_2040$class,6)     #Explore
head(t_test_predictions_2040$posterior,6) #Explore
head(t_test_predictions_2040$x,3)         #Explore

t_test_predictions_2040_df <- data.frame(COMID=t_SIERRA_test$COMID, Cluster_Predicted = t_test_predictions_2040)
class(t_test_predictions_2040_df)

# This is just for plotting; not used for anything else.
x_ploting <- cbind(t_SIERRA_train, predict(t_LDAmodel)$x)
ggplot(x_ploting, aes(cluster_label, LD1))+
  geom_point(aes(color=cluster_label))

# Model Accuracy
class(t_SIERRA_train)
class(t_test_predictions_2040_df)

class(t_test_predictions_2040_df$Cluster_Predicted.class)
class(t_SIERRA_train$cluster_label)

t_SIERRA_train$cluster_label<-as.factor(t_SIERRA_train$cluster_label)

caret::confusionMatrix( t_test_predictions_2040_df$Cluster_Predicted.class, t_SIERRA_train$cluster_label) # Accuracy is here.
mean(t_test_predictions_2040$class==t_SIERRA_train$cluster_label) # Accuracy is here. 



# .--- Test projected for loop: CASCADES ------------------------------------

unique(t_hist_forloop_manualtest$MtnRange_SIMPLE)

t_CASCADES_train <- t_hist_forloop_manualtest%>% 
  dplyr::filter(MtnRange_SIMPLE == "Cascades")
t_CASCADES_test <- t_test_data_2040 %>% 
  dplyr::filter(MtnRange_SIMPLE == "Cascades")

# L.DFA model
t_LDAmodel <- MASS::lda(cluster_label ~ Log1p_MeanGDD, data = t_CASCADES_train)

# Summary Stats/Plots
summary(t_LDAmodel)
t_LDAmodel # Coefficients are here. 
plot(t_LDAmodel)


#Predicting the cluster for new data
t_test_predictions_2040    <- predict(t_LDAmodel, newdata = t_CASCADES_test)
class(t_test_predictions_2040)            #List
head(t_test_predictions_2040$class,6)     #Explore
head(t_test_predictions_2040$posterior,6) #Explore
head(t_test_predictions_2040$x,3)         #Explore

t_test_predictions_2040_df <- data.frame(COMID=t_CASCADES_test$COMID, Cluster_Predicted = t_test_predictions_2040)
class(t_test_predictions_2040_df)

# This is just for plotting; not used for anything else.
x_ploting <- cbind(t_CASCADES_train, predict(t_LDAmodel)$x)
ggplot(x_ploting, aes(cluster_label, LD1))+
  geom_point(aes(color=cluster_label))

# Model Accuracy
class(t_CASCADES_train)
class(t_test_predictions_2040_df)

class(t_test_predictions_2040_df$Cluster_Predicted.class)
class(t_CASCADES_train$cluster_label)

t_CASCADES_train$cluster_label<-as.factor(t_CASCADES_train$cluster_label)

caret::confusionMatrix( t_test_predictions_2040_df$Cluster_Predicted.class, t_CASCADES_train$cluster_label) # Accuracy is here. 
mean(t_test_predictions_2040$class==t_CASCADES_train$cluster_label) # Accuracy is here. 





# .--- 2040 predicted clusters --------------------------------------------

t_hist_forloop_LOOP2040 <- o_recoded

models                    <- list() # List to store models
test_predictions_list     <- list() # List to store predictions
accuracies                <- list() # List to store accuracies
conf_matrices             <- list() # List to store confusion matrices
combined_predictions_list <- list() # List to store all-regions predicted clusters 

all_regions <- unique(t_hist_forloop_LOOP2040$MtnRange_SIMPLE)
all_regions



# Loop through each region
for (region in all_regions) {
  
  # Subset train data for the current region
  train_region_data <- subset(t_hist_forloop_LOOP2040, MtnRange_SIMPLE == region)
  
  # Subset test data for the current region
  test_region_data <- subset(t_test_data_2040, MtnRange_SIMPLE == region)
  
  # Perform LDA model for the current region
  model <- MASS::lda(cluster_label ~ Log1p_MeanGDD, data = train_region_data)
  
  # Store the model with a name indicating the region
  model_name <- paste("model_region", region, sep = "_")
  assign(model_name, model, envir = .GlobalEnv)
  
  # Store the model in the list
  models[[region]] <- model
  
  # Predict the cluster for test data of the current region
  test_predictions_list <- predict(model, newdata = test_region_data)
  
  # Store predictions with attributes
  test_predictions_df <- data.frame(COMID_a             = test_region_data$COMID, 
                                    Cluster_Predicted = test_predictions_list$class)
  
  #Save Result
  assign(paste0(region, "test_predictions_df"), test_predictions_df, envir = .GlobalEnv)
  
  # Store predictions in a list
  test_predictions_list[[region]] <- test_predictions_df
  
  # Model Accuracy for the current region
  train_region_data$cluster_label <- as.factor(train_region_data$cluster_label)
  
  accuracy <- mean(test_predictions_list$class == train_region_data$cluster_label)
  accuracies[[region]] <- accuracy
  
  # Confusion Matrix for the current region
  conf_matrix <- caret::confusionMatrix(test_predictions_list$class, train_region_data$cluster_label)
  conf_matrices[[region]] <- conf_matrix
  
  
  # Combine predictions with attributes
  combined_predictions <- cbind(train_region_data, test_predictions_df)
  combined_predictions_list[[region]] <- combined_predictions
  
}

models
accuracies #2040
conf_matrices

u_DFA_Finished_2040 <- do.call(rbind,combined_predictions_list)
colnames(u_DFA_Finished_2040)





# .--- 2070 predicted clusters --------------------------------------------

t_hist_forloop_LOOP2070 <- o_recoded

models                    <- list() # List to store models
test_predictions_list     <- list() # List to store predictions
accuracies                <- list() # List to store accuracies
conf_matrices             <- list() # List to store confusion matrices
combined_predictions_list <- list() # List to store all-regions predicted clusters 

all_regions <- unique(t_hist_forloop_LOOP2070$MtnRange_SIMPLE)
all_regions



# Loop through each region
for (region in all_regions) {
  
  # Subset train data for the current region
  train_region_data <- subset(t_hist_forloop_LOOP2070, MtnRange_SIMPLE == region)
  
  # Subset test data for the current region
  test_region_data <- subset(t_test_data_2070_2040, MtnRange_SIMPLE == region)
  
  # Perform LDA model for the current region
  model <- MASS::lda(cluster_label ~ Log1p_MeanGDD, data = train_region_data)
  
  # Store the model with a name indicating the region
  model_name <- paste("model_region", region, sep = "_")
  assign(model_name, model, envir = .GlobalEnv)
  
  # Store the model in the list
  models[[region]] <- model
  
  # Predict the cluster for test data of the current region
  test_predictions_list <- predict(model, newdata = test_region_data)
  
  # Store predictions with attributes
  test_predictions_df <- data.frame(COMID_a             = test_region_data$COMID, 
                                    Cluster_Predicted = test_predictions_list$class)
  
  #Save Result
  assign(paste0(region, "test_predictions_df"), test_predictions_df, envir = .GlobalEnv)
  
  # Store predictions in a list
  test_predictions_list[[region]] <- test_predictions_df
  
  # Model Accuracy for the current region
  train_region_data$cluster_label <- as.factor(train_region_data$cluster_label)
  
  accuracy <- mean(test_predictions_list$class == train_region_data$cluster_label)
  accuracies[[region]] <- accuracy
  
  # Confusion Matrix for the current region
  conf_matrix <- caret::confusionMatrix(test_predictions_list$class, train_region_data$cluster_label)
  conf_matrices[[region]] <- conf_matrix
  
  
  # Combine predictions with attributes
  combined_predictions <- cbind(train_region_data, test_predictions_df)
  combined_predictions_list[[region]] <- combined_predictions
  
}


models
accuracies #2070
conf_matrices

u_DFA_Finished_2070 <- do.call(rbind,combined_predictions_list)





# .--- 2100 predicted clusters --------------------------------------------


t_hist_forloop_LOOP2100 <- o_recoded

models                    <- list() # List to store models
test_predictions_list     <- list() # List to store predictions
accuracies                <- list() # List to store accuracies
conf_matrices             <- list() # List to store confusion matrices
combined_predictions_list <- list() # List to store all-regions predicted clusters 

all_regions <- unique(t_hist_forloop_LOOP2100$MtnRange_SIMPLE)
all_regions


# Loop through each region
for (region in all_regions) {
  
  # Subset train data for the current region
  train_region_data <- subset(t_hist_forloop_LOOP2100, MtnRange_SIMPLE == region)
  
  # Subset test data for the current region
  test_region_data <- subset(t_test_data_2100_2070_2040, MtnRange_SIMPLE == region)
  
  # Perform LDA model for the current region
  model <- MASS::lda(cluster_label ~ Log1p_MeanGDD, data = train_region_data)
  
  # Store the model with a name indicating the region
  model_name <- paste("model_region", region, sep = "_")
  assign(model_name, model, envir = .GlobalEnv)
  
  # Store the model in the list
  models[[region]] <- model
  
  # Predict the cluster for test data of the current region
  test_predictions_list <- predict(model, newdata = test_region_data)
  
  # Store predictions with attributes
  test_predictions_df <- data.frame(COMID_a             = test_region_data$COMID, 
                                    Cluster_Predicted = test_predictions_list$class)
  
  #Save Result
  assign(paste0(region, "test_predictions_df"), test_predictions_df, envir = .GlobalEnv)
  
  # Store predictions in a list
  test_predictions_list[[region]] <- test_predictions_df
  
  # Model Accuracy for the current region
  train_region_data$cluster_label <- as.factor(train_region_data$cluster_label)
  
  accuracy <- mean(test_predictions_list$class == train_region_data$cluster_label)
  accuracies[[region]] <- accuracy
  
  # Confusion Matrix for the current region
  conf_matrix <- caret::confusionMatrix(test_predictions_list$class, train_region_data$cluster_label)
  conf_matrices[[region]] <- conf_matrix
  
  
  # Combine predictions with attributes
  combined_predictions <- cbind(train_region_data, test_predictions_df)
  combined_predictions_list[[region]] <- combined_predictions
  
}


models
accuracies #2100
conf_matrices


u_DFA_Finished_2100 <- do.call(rbind,combined_predictions_list)

colnames(u_DFA_Finished_2100)





# ~ DFA % Change ------------------------------------------------


# Check COMID column
x<- u_DFA_Finished_2040 %>% 
  mutate(identcols = ifelse(COMID==COMID_a, "Identical", "Different"))
x_2<- x %>% 
  group_by(identcols) %>% 
  tally() %>% 
  ungroup()

x<- u_DFA_Finished_2070 %>% 
  mutate(identcols = ifelse(COMID==COMID_a, "Identical", "Different"))
x_2<- x %>% 
  group_by(identcols) %>% 
  tally() %>% 
  ungroup()

x<- u_DFA_Finished_2100 %>% 
  mutate(identcols = ifelse(COMID==COMID_a, "Identical", "Different"))
x_2<- x %>% 
  group_by(identcols) %>% 
  tally() %>% 
  ungroup()

x<-u_DFA_Finished_2040 %>% 
  group_by(MtnRange_SIMPLE,cluster_fix,cluster_label) %>% 
  tally() %>% 
  ungroup()

x<-u_DFA_Finished_2040 %>% 
  group_by(MtnRange_SIMPLE,cluster_label,Cluster_Predicted) %>% 
  tally() %>% 
  ungroup()

x<-u_DFA_Finished_2070 %>% 
  group_by(MtnRange_SIMPLE,cluster_label,Cluster_Predicted) %>% 
  tally() %>% 
  ungroup()

x<-u_DFA_Finished_2100 %>% 
  group_by(MtnRange_SIMPLE,cluster_label,Cluster_Predicted) %>% 
  tally() %>% 
  ungroup()




# .--- For All Lakes ------------------------------------------------------

# CLUSTER COUNTS - HISTORIC
v_cluster_counts_historic <- o_recoded %>% 
  group_by(cluster_label) %>% 
  summarize( count_historic = n() ) %>% 
  mutate(  percent_historic = ((count_historic / (sum(count_historic))) *100)) %>% 
  ungroup() %>% 
  rename(Cluster = cluster_label)

# CLUSTER COUNTS - 2040
v_cluster_counts_2040 <- u_DFA_Finished_2040 %>% 
  group_by(Cluster_Predicted) %>% 
  summarize( count_2040 = n() ) %>% 
  mutate(  percent_2040 = ((count_2040 / (sum(count_2040))) *100)) %>% 
  ungroup() %>% 
  rename(Clusterb = Cluster_Predicted)

# CLUSTER COUNTS - 2070
v_cluster_counts_2070 <- u_DFA_Finished_2070 %>% 
  group_by(Cluster_Predicted) %>% 
  summarize( count_2070 = n() ) %>% 
  mutate(  percent_2070 = ((count_2070 / (sum(count_2070))) *100)) %>% 
  ungroup() %>% 
  rename(Clusterc = Cluster_Predicted)
  
# CLUSTER COUNTS - 2100
v_cluster_counts_2100 <- u_DFA_Finished_2100 %>% 
  group_by(Cluster_Predicted) %>% 
  summarize( count_2100 = n() ) %>% 
  mutate(  percent_2100 = ((count_2100 / (sum(count_2100))) *100)) %>% 
  ungroup() %>% 
  rename(Clusterd = Cluster_Predicted)

v_cluster_counts_all <- cbind(v_cluster_counts_historic,
                                  v_cluster_counts_2040,
                                  v_cluster_counts_2070,
                                  v_cluster_counts_2100)

v_cluster_counts_all <- v_cluster_counts_all %>% 
  dplyr::select(-c(Clusterb, Clusterc, Clusterd))


x1<-o_recoded %>% 
  tally() %>% 
  ungroup()

x2<-u_DFA_Finished_2040 %>% 
  tally() %>% 
  ungroup()


v_percentchange_all <- v_cluster_counts_all %>% 
  group_by(Cluster) %>% 
  mutate ( PChange_historic_to_2040 = (((count_2040 - count_historic) / count_historic) * 100),
           PChange_historic_to_2070 = (((count_2070 - count_historic) / count_historic) * 100),
           PChange_historic_to_2100 = (((count_2100 - count_historic) / count_historic) * 100),
           
           PDiff_historic_to_2040 = (percent_2040 - percent_historic),
           PDiff_historic_to_2070 = (percent_2070 - percent_historic),
           PDiff_historic_to_2100 = (percent_2100 - percent_historic)) %>% 
  ungroup()


write_csv(v_percentchange_all, "data_output/Percent_Change_Difference_2024-06-26_AllRegionsCombined.csv")



# .--- For Each Regoin ------------------------------------------------------
class(o_recoded)

# CLUSTER COUNTS - HISTORIC
v_cluster_counts_historic <- o_recoded %>% 
  group_by(MtnRange_SIMPLE, cluster_label) %>% 
  summarize( count_historic = n(), .groups = "keep") %>% 
  ungroup() %>% 
  group_by(MtnRange_SIMPLE) %>% 
  mutate(  percent_historic = ((count_historic / (sum(count_historic))) *100)) %>% 
  ungroup() %>% 
  rename(Cluster = cluster_label)

# CLUSTER COUNTS - 2040
v_cluster_counts_2040 <- u_DFA_Finished_2040 %>% 
  group_by(MtnRange_SIMPLE, Cluster_Predicted) %>% 
  summarize( count_2040 = n(), .groups = "keep" ) %>%
  ungroup() %>% 
  group_by(MtnRange_SIMPLE) %>%
  mutate(  percent_2040 = ((count_2040 / (sum(count_2040))) *100)) %>% 
  ungroup() %>% 
  mutate(Cluster = as.character(Cluster_Predicted)) %>% 
  dplyr::select(-Cluster_Predicted)

# CLUSTER COUNTS - 2070
v_cluster_counts_2070 <- u_DFA_Finished_2070 %>% 
  group_by(MtnRange_SIMPLE, Cluster_Predicted) %>% 
  summarize( count_2070 = n(), .groups = "keep" ) %>% 
  ungroup() %>% 
  group_by(MtnRange_SIMPLE) %>%
  mutate(  percent_2070 = ((count_2070 / (sum(count_2070))) *100)) %>% 
  ungroup() %>% 
  mutate(Cluster = as.character(Cluster_Predicted))%>% 
  dplyr::select(-Cluster_Predicted)

# CLUSTER COUNTS - 2100
v_cluster_counts_2100 <- u_DFA_Finished_2100 %>% 
  group_by(MtnRange_SIMPLE, Cluster_Predicted) %>% 
  summarize( count_2100 = n(), .groups = "keep" ) %>% 
  ungroup() %>% 
  group_by(MtnRange_SIMPLE) %>%
  mutate(  percent_2100 = ((count_2100 / (sum(count_2100))) *100)) %>% 
  ungroup() %>% 
  mutate(Cluster = as.character(Cluster_Predicted))%>% 
  dplyr::select(-Cluster_Predicted)

v_cluster_counts_hist_2040     <- left_join(v_cluster_counts_historic,v_cluster_counts_2040,      by = c("MtnRange_SIMPLE","Cluster"))
v_cluster_counts_hist2040_2070 <- left_join(v_cluster_counts_hist_2040,v_cluster_counts_2070,     by = c("MtnRange_SIMPLE","Cluster"))
v_cluster_counts_all           <- left_join(v_cluster_counts_hist2040_2070,v_cluster_counts_2100, by = c("MtnRange_SIMPLE","Cluster"))

# Fill NA with zero
v_cluster_counts_all[is.na(v_cluster_counts_all)] <- 0

class(v_cluster_counts_all)
v_percentchange_all <- v_cluster_counts_all %>% 
  group_by(MtnRange_SIMPLE, Cluster) %>% 
  mutate ( PChange_historic_to_2040 = (((count_2040 - count_historic) / count_historic) * 100),
           PChange_historic_to_2070 = (((count_2070 - count_historic) / count_historic) * 100),
           PChange_historic_to_2100 = (((count_2100 - count_historic) / count_historic) * 100),
           PDiff_historic_to_2040 = (percent_2040 - percent_historic),
           PDiff_historic_to_2070 = (percent_2070 - percent_historic),
           PDiff_historic_to_2100 = (percent_2100 - percent_historic)) %>% 
  ungroup()

write_csv(v_percentchange_all, "data_output/Percent_Change_Difference_2024-06-26--RCP7.0--Cluster.and.DFA.on--Log1p-MeanGDD--NotScaled.csv")





# * KDD % Change Over TIme ------------------------------------------------

# Mean KDD - per range - HISTORIC
w_meanKDDperRange_historic<- o_recoded %>% 
  group_by(MtnRange_SIMPLE) %>% 
  summarise( MeanKDDperRange_2018 = mean(Mean_KDD)) %>% 
  ungroup()

colnames(s_MeanGDD_with_LakeIntSlope_2100_2070_2040)
# Mean KDD - per range - 2040
w_meanKDDperRange_2040 <- s_MeanGDD_with_LakeIntSlope_2100_2070_2040 %>% 
  group_by(MtnRange_SIMPLE) %>% 
  summarise( MeanKDDperRange_2040 = mean(Mean_KDD_2040)) %>% 
  ungroup()

# Mean KDD - per range - 2070
w_meanKDDperRange_2070 <- s_MeanGDD_with_LakeIntSlope_2100_2070_2040 %>% 
  group_by(MtnRange_SIMPLE) %>% 
  summarise( MeanKDDperRange_2070 = mean(Mean_KDD_2070_2040)) %>% 
  ungroup()

# Mean KDD - per range - 2100
w_meanKDDperRange_2100 <- s_MeanGDD_with_LakeIntSlope_2100_2070_2040 %>% 
  group_by(MtnRange_SIMPLE) %>% 
  summarise( MeanKDDperRange_2100 = mean(Mean_KDD_2100_2070_2040)) %>% 
  ungroup() 

w_meanKDDperrange_hist_2040     <- left_join(w_meanKDDperRange_historic,      w_meanKDDperRange_2040, by = c("MtnRange_SIMPLE"))
w_meanKDDperrange_hist2040_2070 <- left_join(w_meanKDDperrange_hist_2040,     w_meanKDDperRange_2070, by = c("MtnRange_SIMPLE"))
w_meanKDDperrange_all           <- left_join(w_meanKDDperrange_hist2040_2070, w_meanKDDperRange_2100, by = c("MtnRange_SIMPLE"))

w_percentchange_KDDMEAN_all <- w_meanKDDperrange_all %>% 
  group_by(MtnRange_SIMPLE) %>% 
  mutate ( PChange_historic_to_2040 = (((MeanKDDperRange_2040 - MeanKDDperRange_2018) / MeanKDDperRange_2018) * 100),
           PChange_historic_to_2070 = (((MeanKDDperRange_2070 - MeanKDDperRange_2018) / MeanKDDperRange_2018) * 100),
           PChange_historic_to_2100 = (((MeanKDDperRange_2100 - MeanKDDperRange_2018) / MeanKDDperRange_2018) * 100)) %>% 
  ungroup()

write_csv(v_percentchange_all, "data_output/Percent_Change_Difference_2024-06-26--RCP7.0--KDD.csv")





# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
# . - . - . - . - . - . - . - . - . - . - . - . - . -----------------------
# * MAP ------------------------------------------------------------------
# Figure 1

# .--- map prep: projected ------------------------------------------------
mapdata_projected <- s_MeanGDD_with_LakeIntSlope_2100_2070_2040 %>% 
  dplyr::filter(!is.na(latitude))

fBasics::basicStats(mapdata_projected$Mean_KDD_2040)
fBasics::basicStats(mapdata_projected$Mean_KDD_2070_2040)
fBasics::basicStats(mapdata_projected$Mean_KDD_2100_2070_2040)

fBasics::basicStats(mapdata_projected$Mean_GDD_2040)
fBasics::basicStats(mapdata_projected$Mean_GDD_2070_2040)
fBasics::basicStats(mapdata_projected$Mean_GDD_2100_2070_2040)

x<- mapdata_projected %>% 
  dplyr::filter(Mean_KDD_2040>1500)
unique(x$US_L3NA)


fBasics::basicStats(o_recoded$Mean_GDD)
fBasics::basicStats(mapdata_projected$Mean_GDD_2040)
fBasics::basicStats(mapdata_projected$Mean_GDD_2070_2040)
fBasics::basicStats(mapdata_projected$Mean_GDD_2100_2070_2040)

sum(is.na(mapdata_projected$latitude))
sum(is.na(mapdata_projected$longitude))

crsnad83<-'+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0'
crsaea <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

nad83spatial<-with(mapdata_projected, SpatialPoints(coords=cbind(mapdata_projected$longitude, mapdata_projected$latitude), proj4string = CRS(crsnad83)))

sptransformnad83<-spTransform(nad83spatial, CRSobj=CRS(crsaea)) #reprojected to albers equal area projection

generic_nad83 <- as.data.frame(sptransformnad83@coords) #with just lat/long, back to df 

wgs1984.proj <- '+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'

names(generic_nad83)<-c("Longitude_albers", "Latitude_albers")
mapdata_projected_POINTS<-cbind(mapdata_projected, generic_nad83)

state<-map_data('state')
##  view(state)                              #explore
##  ca<-subset(states, region=="california") #explore; CA lat/long data

CAP_map_coords <- SpatialPoints(coords = with(state, data.frame(x = long, y = lat)), proj4string = CRS(wgs1984.proj)) #here are the coordinates in that frame, and this is the CRS they're in 
res_albers <- spTransform(CAP_map_coords, CRSobj = CRS(crsaea)) #put latlong to this CRS
state$long <- res_albers@coords[,1]
state$lat <- res_albers@coords[,2]

usa <- map_data("usa")
CAP_map_coords <- SpatialPoints(coords = with(usa, data.frame(x = long, y = lat)), proj4string = CRS(wgs1984.proj)) 
res_albers <- spTransform(CAP_map_coords, CRSobj = CRS(crsaea)) #put latlong to this CRS
usa$long <- res_albers@coords[,1]
usa$lat <- res_albers@coords[,2]



# .--- GDD projected ------------------------------------------------------
colnames(mapdata_projected_POINTS)

size <- 2
alpha <- 0.02 # For manuscript figure readability
alpha <- 0.5 # For thumbnail TOC Graphic Art
shape <- 16

plotmap_GDD_projected1 <- ggplot()+
  geom_polygon(data= state, aes(x=long, y=lat, group=group),fill="white", color="gray90")+
  geom_polygon(data= usa,   aes(x=long, y=lat, group=group),fill=NA,      color="black")+
  coord_fixed(1.1)+
  geom_point(data=mapdata_projected_POINTS, 
             aes(x=Longitude_albers, y=Latitude_albers, color=Mean_GDD_2040), 
             size=size, alpha=alpha, shape=shape)+
  scale_color_gradientn(colors=rev(lacroix_palette("Coconut", n = 6, type = "continuous")), 
                        breaks = seq(0,9338,2000),
                        limits = c(0,9338))+ # GOOD 6
  theme_bw()+
  labs(color= "Mean GDD 2011-2040")+
  theme(panel.border = element_rect(colour = "white", linewidth=0.4),
        #legend.position = c(1.0,0.98),
        #legend.margin   = margin(t=6, r=6, b=6, l=6),
        #theme(plot.margin = unit(c(1,1,1,1), "mm")),
        axis.ticks = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text  = element_blank(),
        panel.grid = element_blank())


ggsave(plot = plotmap_GDD_projected1, "figures/06.18.2025 - TOC Graphic Art.tiff", height = 1.5, width = 3, units="in", dpi = 900)



+plotmap_GDD_projected2 <- ggplot()+
  geom_polygon(data= state, aes(x=long, y=lat, group=group),fill="white", color="gray90")+
  geom_polygon(data= usa,   aes(x=long, y=lat, group=group),fill=NA,      color="black")+
  coord_fixed(1.1)+
  geom_point(data=mapdata_projected_POINTS, 
             aes(x=Longitude_albers, y=Latitude_albers, color=Mean_GDD_2070_2040), 
             size=size, alpha=alpha, shape=shape)+
  scale_color_gradientn(colors=rev(lacroix_palette("Coconut", n = 6, type = "continuous")), 
                        breaks = seq(0,9338,2000),
                        limits = c(0,9338))+ # GOOD 6
  theme_bw()+
  labs(color= "Mean GDD 2041-2070")+
  theme(panel.border = element_rect(colour = "white", linewidth=0.4),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text  = element_blank())



plotmap_GDD_projected3 <- ggplot()+
  geom_polygon(data= state, aes(x=long, y=lat, group=group),fill="white", color="gray90")+
  geom_polygon(data= usa,   aes(x=long, y=lat, group=group),fill=NA,      color="black")+
  coord_fixed(1.1)+
  geom_point(data=mapdata_projected_POINTS, 
             aes(x=Longitude_albers, y=Latitude_albers, color=Mean_GDD_2100_2070_2040), 
             size=size, alpha=alpha, shape=shape)+
  scale_color_gradientn(colors=rev(lacroix_palette("Coconut", n = 6, type = "continuous")), 
                        breaks = seq(0,9338,2000),
                        limits = c(0,9338))+ # GOOD 6
  theme_bw()+
  labs(color= "Mean GDD 2071-2100")+
  theme(panel.border = element_rect(colour = "white", linewidth=0.4),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text  = element_blank())+
  ggspatial::annotation_scale(location = "bl", 
                              style = "ticks",
                              line_width = 0.5,
                              height = unit(0.1,"cm"),
                              tick_height = 0.4,
                              pad_x = unit(0.05,"in"), 
                              pad_y=unit(0.2,"in"))


# .--- KDD projected ------------------------------------------------------

mapdata_kdd2040_nozero <- mapdata_projected_POINTS %>% 
  dplyr::filter(!Mean_KDD_2040 == 0)

plotmap_KDD_projected1 <- ggplot()+
  geom_polygon(data= state, aes(x=long, y=lat, group=group),fill="white", color="gray90")+
  geom_polygon(data= usa,   aes(x=long, y=lat, group=group),fill=NA,      color="black")+
  coord_fixed(1.1)+
  geom_point(data=mapdata_kdd2040_nozero, 
             aes(x=Longitude_albers, y=Latitude_albers, color=Mean_KDD_2040), 
             size=size, alpha=alpha, shape=shape)+
  scale_color_gradientn(colors=rev(lacroix_palette("PassionFruit", n = 6, type = "continuous")), 
                        breaks = seq(0,2018,500),
                        limits = c(0,2018)) + # GOOD: 6, 10
  theme_bw()+
  labs(color= "Mean KDD 2011-2040")+
  theme(panel.border = element_rect(colour = "white", linewidth=0.4),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text  = element_blank())



mapdata_kdd2070_nozero <- mapdata_projected_POINTS %>% 
  dplyr::filter(!Mean_KDD_2070_2040 == 0)

plotmap_KDD_projected2 <- ggplot()+
  geom_polygon(data= state, aes(x=long, y=lat, group=group),fill="white", color="gray90")+
  geom_polygon(data= usa,   aes(x=long, y=lat, group=group),fill=NA,      color="black")+
  coord_fixed(1.1)+
  geom_point(data=mapdata_kdd2070_nozero, 
             aes(x=Longitude_albers, y=Latitude_albers, color=Mean_KDD_2070_2040), 
             size=size, alpha=alpha, shape=shape)+
  scale_color_gradientn(colors=rev(lacroix_palette("PassionFruit", n = 6, type = "continuous")), 
                        breaks = seq(0,2018,500),
                        limits = c(0,2018)) + # GOOD: 6, 10
  theme_bw()+
  labs(color= "Mean KDD 2041-2070")+
  theme(panel.border = element_rect(colour = "white", linewidth=0.4),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text  = element_blank())


mapdata_kdd2100_nozero <- mapdata_projected_POINTS %>% 
  dplyr::filter(!Mean_KDD_2100_2070_2040 == 0)

plotmap_KDD_projected3 <- ggplot()+
  geom_polygon(data= state, aes(x=long, y=lat, group=group),fill="white", color="gray90")+
  geom_polygon(data= usa,   aes(x=long, y=lat, group=group),fill=NA,      color="black")+
  coord_fixed(1.1)+
  geom_point(data=mapdata_kdd2100_nozero, 
             aes(x=Longitude_albers, y=Latitude_albers, color=Mean_KDD_2100_2070_2040), 
             size=size, alpha=alpha, shape=shape)+
  scale_color_gradientn(colors=rev(lacroix_palette("PassionFruit", n = 6, type = "continuous")), 
                        breaks = seq(0,2018,500),
                        limits = c(0,2018)) + # GOOD: 6, 10
  theme_bw()+
  labs(color= "Mean KDD 2071-2100")+
  theme(panel.border = element_rect(colour = "white", linewidth=0.4),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text  = element_blank())+
  ggspatial::annotation_north_arrow(location = "br",
                                    width = unit(.30,"in"),
                                    height = unit(.49,"in"),
                                    style = north_arrow_minimal,
                                    pad_x = unit(0.2,"in"), 
                                    pad_y=unit(0.2,"in"))








# .--- map prep: historic -------------------------------------------------

mapdata_historic<-o_recoded %>% 
  dplyr::filter(!is.na(latitude))

fBasics::basicStats(mapdata_historic$Mean_GDD)
sum(is.na(mapdata_historic$latitude))
sum(is.na(mapdata_historic$longitude))

crsnad83<-'+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0'
crsaea  <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

nad83spatial<-with(mapdata_historic, SpatialPoints(coords=cbind(mapdata_historic$longitude, mapdata_historic$latitude), proj4string = CRS(crsnad83)))

sptransformnad83<-spTransform(nad83spatial, CRSobj=CRS(crsaea)) #reprojected to albers equal area projection

generic_nad83 <- as.data.frame(sptransformnad83@coords) #with just lat/long, back to df 

wgs1984.proj <- '+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'

names(generic_nad83)<-c("Longitude_albers", "Latitude_albers")
mapdata_historicPOINT<-cbind(mapdata_historic, generic_nad83)

state<-map_data('state')
##  view(state)                              #explore
##  ca<-subset(states, region=="california") #explore; CA lat/long data

CAP_map_coords <- SpatialPoints(coords = with(state, data.frame(x = long, y = lat)), proj4string = CRS(wgs1984.proj)) #here are the coordinates in that frame, and this is the CRS they're in 
res_albers <- spTransform(CAP_map_coords, CRSobj = CRS(crsaea)) #put latlong to this CRS
state$long <- res_albers@coords[,1]
state$lat <- res_albers@coords[,2]

usa <- map_data("usa")
CAP_map_coords <- SpatialPoints(coords = with(usa, data.frame(x = long, y = lat)), proj4string = CRS(wgs1984.proj)) 
res_albers <- spTransform(CAP_map_coords, CRSobj = CRS(crsaea)) #put latlong to this CRS
usa$long <- res_albers@coords[,1]
usa$lat <- res_albers@coords[,2]



# .--- GDD historic -------------------------------------------------------
colnames(mapdata_historicPOINT)



plotmap_GDD_historic <- ggplot()+
  geom_polygon(data= state, aes(x=long, y=lat, group=group),fill="white", color="gray90")+
  geom_polygon(data= usa,   aes(x=long, y=lat, group=group),fill=NA,      color="black")+
  coord_fixed(1.1)+
  geom_point(data=mapdata_historicPOINT, 
             aes(x=Longitude_albers, y=Latitude_albers, color=Mean_GDD), 
             size=size, alpha=alpha, shape=shape)+
  #scale_color_gradientn(colors=rev(lacroix_palette("Pamplemousse", n = 5, type = "continuous")), 
  #                      breaks = seq(0,300,50))+ # GOOD 5
  #scale_color_gradientn(colors=rev(lacroix_palette("PassionFruit", n = 6, type = "continuous")), 
  #                      breaks = seq(0,300,50))+ # GOOD: 6, 10
  scale_color_gradientn(colors=rev(lacroix_palette("Coconut", n = 6, type = "continuous")), 
                        breaks = seq(0,9338,2000),
                        limits = c(0,9338),
                        label = scales::comma)+ # GOOD 6
  theme_bw()+
  labs(color= "Mean GDD")+
  theme(panel.border = element_rect(colour = "white", linewidth=0.4),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text  = element_blank(),
        legend.title = element_text(size=13),
        legend.text = element_text(size=11))


# .--- KDD historic -------------------------------------------------------

# this object is created by the end of the prof.maps part
mapdata_kdd2018_nozero <- mapdata_historicPOINT %>% 
  dplyr::filter(!Mean_KDD == 0)

fBasics::basicStats(mapdata_kdd2018_nozero$Mean_KDD)                # 0 - 1472.151282
fBasics::basicStats(mapdata_kdd2040_nozero$Mean_KDD_2040)           # 0 - 1751.200000
fBasics::basicStats(mapdata_kdd2070_nozero$Mean_KDD_2070_2040)      # 0 - 1873.600000
fBasics::basicStats(mapdata_kdd2100_nozero$Mean_KDD_2100_2070_2040) # 0 - 2017.266667



plotmap_KDD_historic <- ggplot()+
  geom_polygon(data= state, aes(x=long, y=lat, group=group),fill="white", color="gray90")+
  geom_polygon(data= usa,   aes(x=long, y=lat, group=group),fill=NA,      color="black")+
  coord_fixed(1.1)+
  geom_point(data=mapdata_kdd2018_nozero, 
             aes(x=Longitude_albers, y=Latitude_albers, color=Mean_KDD), 
             size=size, alpha=alpha, shape=shape)+
  #scale_color_gradientn(colors=rev(lacroix_palette("Pamplemousse", n = 5, type = "continuous")), 
  #                      breaks = seq(0,300,50))+ # GOOD 5
  scale_color_gradientn(colors=rev(lacroix_palette("PassionFruit", n = 6, type = "continuous")), 
                        breaks = seq(0,2018,500),
                        limits = c(0,2018),
                        label = scales::comma) + # GOOD: 6, 10
  theme_bw()+
  labs(color= "Mean KDD")+
  theme(panel.border = element_rect(colour = "white", linewidth=0.4),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text  = element_blank(),
        legend.title = element_text(size=13),
        legend.text = element_text(size=11))




# .--- plot maps: GDD KDD -------------------------------------------------

x_gdd_legend <- get_legend(plotmap_GDD_historic)
x_kdd_legend <- get_legend(plotmap_KDD_historic)

plotmap_GDD_historic <- plotmap_GDD_historic +  theme(legend.position = "none")
plotmap_KDD_historic <- plotmap_KDD_historic +  theme(legend.position = "none")

y.grob_A <- grid::textGrob("1980 - 2019", gp=gpar(fontsize=15), rot=90)
y.grob_B <- grid::textGrob("2011 - 2040", gp=gpar(fontsize=15), rot=90)
y.grob_C <- grid::textGrob("2041 - 2070", gp=gpar(fontsize=15), rot=90)
y.grob_D <- grid::textGrob("2071 - 2100", gp=gpar(fontsize=15), rot=90)

plots_map_A <- plot_grid(plotmap_GDD_historic  , plotmap_KDD_historic  , ncol=2 )
plots_map_B <- plot_grid(plotmap_GDD_projected1, plotmap_KDD_projected1, ncol=2 )
plots_map_C <- plot_grid(plotmap_GDD_projected2, plotmap_KDD_projected2, ncol=2 )
plots_map_D <- plot_grid(plotmap_GDD_projected3, plotmap_KDD_projected3, ncol=2 )

plot_map_labeled_A <- gridExtra::grid.arrange(arrangeGrob(plots_map_A, left = y.grob_A))
plot_map_labeled_B <- gridExtra::grid.arrange(arrangeGrob(plots_map_B, left = y.grob_B))
plot_map_labeled_C <- gridExtra::grid.arrange(arrangeGrob(plots_map_C, left = y.grob_C))
plot_map_labeled_D <- gridExtra::grid.arrange(arrangeGrob(plots_map_D, left = y.grob_D))

x_combined_maps<-cowplot::plot_grid(plot_map_labeled_A,
                                    plot_map_labeled_B,
                                    plot_map_labeled_C,
                                    plot_map_labeled_D,
                                    ncol=1, align = "v")

x_spacer <- plot_grid(NULL, ncol=1, rel_heights = c(1))

x_legendscombined <- plot_grid(x_spacer, x_gdd_legend, x_kdd_legend, x_spacer, ncol=1)

# plots_all_with_legend <- plot_grid(x_combined_maps, 
#                                    x_legendscombined,
#                                    ncol=2, rel_widths = c(4,0.5)) # Original

plots_all_with_legend <- plot_grid(x_gdd_legend,
                                   x_combined_maps, 
                                   x_kdd_legend,
                                   ncol=3, rel_widths = c(0.7,4,0.65))

ggsave(plot = plots_all_with_legend, "figures/06.17.2025 - Map.jpeg", height = 10, width = 8, units="in", dpi = 300)


# . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ -------------------------------------
# . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ -------------------------------------
# . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ -------------------------------------
# * HISTOGRAMS & SKEW -------------------------------------------------
# Figure S1 & Table S1

rm(x1)
rm(x2)
rm(x3)

# . Get distinct lakes ----------------------------------------------------
# TAKE MEAN TEMP FOR A LAKE ACROSS YEARS. THEN GET DF WITH DISTINCT LAKES NOT LAKE-YEAR. 

x<-g_NHDClimate_MonthYear %>% dplyr::filter(COMID=="3800743") # Explore

z1_ClimateCountsData<-g_NHDClimate_MonthYear %>% 
  group_by(COMID) %>% 
  mutate( MeanTemp = (mean(Temp_C))) %>% 
  ungroup() %>% 
  dplyr::distinct(COMID, .keep_all = TRUE)  # 1 row per lake

x<-z1_ClimateCountsData %>% dplyr::filter(COMID=="3800743") # Explore

x<-z1_ClimateCountsData %>% 
  group_by(MtnRange_SIMPLE) %>% 
  tally() %>% 
  ungroup()

# . Get Skew & Kurtosis values --------------------------------------------
library(moments)

colnames(z1_ClimateCountsData)

x2_Skew_and_Kurtosis_Results <- z1_ClimateCountsData %>%
  group_by(MtnRange_SIMPLE) %>%
  summarize(Skewness_elevation_M      = skewness(elevatn, na.rm = TRUE),
            Kurtosis_elevation_M      = kurtosis(elevatn, na.rm = TRUE),
            Mean_elvation_M           = mean(elevatn, na.rm = TRUE),
            Skewness_SurfaceArea_SQKM = skewness(AREASQK, na.rm = TRUE),
            Kurtosis_SurfaceArea_SQKM = kurtosis(AREASQK, na.rm = TRUE),
            Mean_SurfaceArea_SQKM     = mean(AREASQK, na.rm = TRUE))

write_csv(x2_Skew_and_Kurtosis_Results,"data output/Skew_and_Kurtosis_Results.csv")

# fig: x=temp, y=elev -----------------------------------------------------

TempELEV <- z1_ClimateCountsData %>% 
  ggplot(aes(x=elevatn, y=MeanTemp, color=MtnRange_SIMPLE))+ 
  geom_line(size=0.5,alpha=0.5)+
  facet_wrap(~MtnRange_SIMPLE,ncol=10,
             labeller = labeller(MtnRange_SIMPLE = label_wrap_gen(width = 8)))+
  labs(x="Elevation (m)", y = expression(atop("Mean", "Temperature (" * degree * C * ")")))+
  scale_y_continuous(breaks = seq(-10,20,5))+
  scale_x_continuous(labels = scales::comma)+
  theme_bw()+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=16),
        #strip.text = element_blank(),
        #strip.placement = "outside",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=13,,angle = -70, hjust=0.1, vjust=0.5),
        axis.text.y = element_text(size=13),
        axis.ticks.y = element_blank(),
        plot.margin = margin(1,6,1,1), 
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none")

TempELEV



# fig: x=temp.mean, y=freq -----------------------------------------------------

Temp <- z1_ClimateCountsData %>%  
  ggplot(aes(x=MeanTemp,fill=MtnRange_SIMPLE))+ 
  geom_histogram(aes(y=(after_stat(density)/max(after_stat(density))), alpha=0.9), bins=200)+
  facet_wrap(~MtnRange_SIMPLE,ncol=10)+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks = seq(-10,20,5))+
  labs(y="Frequency (%)", x = expression("Mean Temperature ("*degree*C*")"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=13, hjust=0.5),
        axis.text.y = element_text(size=13),
        plot.margin = margin(1,6,1,1), 
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none")

Temp




# fig: x=area, y=freq -----------------------------------------------------

fBasics::basicStats(z1_ClimateCountsData$AREASQK)
log10(496.922000)
log10(0)

Area <- z1_ClimateCountsData %>% 
  ggplot(aes(x=log10(AREASQK), fill=MtnRange_SIMPLE)) +
  geom_histogram(aes(y=(after_stat(density)/max(after_stat(density))), alpha=0.9), bins=200)+
  #geom_freqpoly(aes(y=(after_stat(density)/max(after_stat(density))), 
  #                  alpha=0.5, color = MtnRange_SIMPLE ), bins=30)+
  facet_wrap(~MtnRange_SIMPLE,ncol=10)+
  scale_y_continuous(labels = scales::percent)+
  #scale_x_continuous(breaks = seq(0, 3, 0.4))+
  labs(y="Frequency (%)", x = expression("Surface area (log, "* km^2 *")"))+
  #xlab(expression("Surface Area"~('log,'~km^2))) + 
  theme_bw()+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=13, hjust=0.5),
        axis.text.y = element_text(size=13),
        axis.ticks.y = element_blank(),
        plot.margin = margin(1,6,1,1), 
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none")

Area


# fig: x=elev, y=freq -----------------------------------------------------

Elev <- z1_ClimateCountsData %>% 
  ggplot(aes(x=elevatn,fill=MtnRange_SIMPLE))+ 
  geom_histogram(aes(y=(after_stat(density)/max(after_stat(density))), alpha=0.9), bins=200)+
  facet_wrap(~MtnRange_SIMPLE,ncol=10)+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::comma)+
  labs(y="Frequency (%)", x = expression("Elevation ("* m *")"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=13,angle = -70,hjust=0.1),
        axis.text.y = element_text(size=13),
        axis.ticks.y = element_blank(),
        plot.margin = margin(1,6,1,1), 
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none")

Elev

z1_ClimateCountsData %>% 
  ggplot(aes(y=elevatn,fill=MtnRange_SIMPLE))+ 
  geom_histogram(aes(x=(after_stat(density)/max(after_stat(density))), alpha=0.9), bins=200)+
  facet_wrap(~MtnRange_SIMPLE,ncol=2)+
  scale_x_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::comma)+
  labs(y="Elevation (m)",x="Frequency (%)")+
  theme_bw()+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=13),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")




# fig: x=elev, y=SA-----------------------------------------------------

SAELEV <- z1_ClimateCountsData %>% 
  ggplot(aes(x=elevatn, y=(AREASQK), color=MtnRange_SIMPLE))+ 
  geom_line(size=0.5,alpha=0.5)+
  scale_x_continuous(labels = scales::comma)+
  facet_wrap(~MtnRange_SIMPLE,ncol=10,strip.position = "right",
             labeller = labeller(MtnRange_SIMPLE = label_wrap_gen(width = 10)))+
  labs(x="Elevation (m)")+
  ylab(expression("Surface Area"~(''~km^2))) + 
  #ylab(expression(atop("Surface Area", (km^2))))+
  #ylab(expression("Surface Area"~('log,'~km^2))) + 
  #scale_color_manual(values = colorRampPalette(sunset(1))(LEVELS)) +
  theme_bw()+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=13,angle = -70, hjust=0.1, vjust=0.5),
        axis.text.y = element_text(size=13),
        axis.ticks.y = element_blank(),
        plot.margin = margin(1,6,1,1), 
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none")

SAELEV

# bind figs ---------------------------------------------------------------
library(cowplot)
plots<-cowplot::plot_grid(TempELEV, Temp, Area, Elev, SAELEV, ncol=1)

plots <- cowplot::plot_grid(
  TempELEV, Temp, Area, Elev, SAELEV,
  ncol = 1, 
  align = "hv",       # Align horizontally and vertically
  axis = "tblr")       # Align axes on all sides (top, bottom, left, right)

#plots <- cowplot::plot_grid(
#  Temp, Area, Elev, SAELEV, 
#  ncol = 4, 
#  align = "hv", 
#  axis = "tblr", 
#  plot_spacing = 0.01)          # Adjust row heights if necessary

ggsave(plot = plots,  
       "figures/Histograms_Elev_Area_Temp_2025-01-09.jpeg", 
       units = "in", 
       height = 13.0, 
       width  = 17, 
       dpi = 300)






