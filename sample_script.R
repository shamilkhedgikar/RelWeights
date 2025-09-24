# Download data using - https://github.com/shamilkhedgikar/sample_shps.git

require(spatialreg)
require(sf)
require(mapview)
library(RColorBrewer)

rel_layer <- st_read("data/india_2012-17_AC/india_2012-17_AC.shp")
rel_layer <- rel_layer %>% st_make_valid()

inh_layer <- st_read("data/Sub-basins/Sub-basins.shp")
inh_layer <- inh_layer %>% st_transform(4326)
inh_layer <- inh_layer %>% st_make_valid()

rel_cols <- brewer.pal(max(3, nrow(rel_layer)), "Reds")
inh_cols <- brewer.pal(max(3, nrow(inh_layer)), "Blues")

rel_map <- mapview(
  rel_layer,
  layer.name    = "Rel Layer",
  alpha.regions = 0.45,
  color         = "firebrick3",
  lwd           = 1.2
)

inh_map <- mapview(
  inh_layer,
  layer.name    = "Inherited Layer",
  alpha.regions = 0.25,
  color         = "navy",
  lwd           = 1.2
)

rel_map + inh_map

# Load Functions tab2relweights_spdep or tab2relweights_sfdep
source("relweights_spdep.R")
source("relweights_sfdep.R")

tab2relweights_spdep(rel_layer = rel_layer, inh_layer = inh_layer, output_gal = "relweights_sample_new.gal")

tab2relweights_sfdep(rel_layer = rel_layer, inh_layer = inh_layer, style = "W", binary = TRUE, zero.policy = TRUE, output_gal = "relweights_sample_new_sfdep_v2.gal")
