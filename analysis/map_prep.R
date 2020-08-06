# Load packages
library(tidyverse)

library(rgdal)
library(maptools)
if (!require(gpclib)) install.packages("gpclib", type="source")
gpclibPermit()

# Load Myvatn map
myv = readOGR('data/myvatn_shape/Myvatn_WSGUTM28.shp',
              layer='Myvatn_WSGUTM28', p4s = "+proj=utm +zone=28")


# Data frame of coordinates for plotting
# Filter out islands ("holes")
myv_2 = myv
myv_2@data$id = rownames(myv_2@data)
myv_points = fortify(myv_2, region="id")
myv_df = as_tibble(left_join(myv_points, myv_2@data, by="id")) %>% 
  select(lat, long, piece) 
