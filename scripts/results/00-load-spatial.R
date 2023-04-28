
# read in elevation raster
base_rast <- rast(paste0(dd_env, "ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif"))
base_rast <- base_rast[[1]]

# all env rast
all_env_rast <- trim(rast(paste0(dd_env, "ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif")))

# set crs
target_crs <- "epsg:3035"

# read in lakes
lakes <- paste0(dd_env, "ch-lake-reference.shp")
lakes <- st_union(st_read(lakes))

# read in rivers
rivers <- paste0(dd_env, "CH_HYDRO.shp")
rivers <- st_read(rivers)
rivers <- rivers %>%
  filter(CUM_LEN > 500, STRAHLE > 1) %>%
  mutate(CUM_LEN_LOG = log(.$CUM_LEN))

# read in catchment object
subcatchment_file <- paste0(dd_ch, "swiss-2km-subcatchments/EZG_Gewaesser.gdb")
catchments_rhine <- c("Aare", "Reuss", "Limmat", "Rhein")

# read in subcatchments, transform, union
subcatchments_rhine_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(FLUSSGB %in% catchments_rhine) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm()

# mask over rhein
base_rast <- terra::trim(terra::mask(base_rast, vect(subcatchments_rhine_union)))
lakes <- st_intersection(lakes, subcatchments_rhine_union)
rivers <- st_intersection(rivers, subcatchments_rhine_union)
rivers_int <- st_union(rivers %>% dplyr::select(geometry))
river_intersect_lakes <- st_difference(rivers_int, lakes)

# make non-unioned object
subcatchments_rhine <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(FLUSSGB %in% catchments_rhine) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # remove z and m properties that can cause errors later
  st_zm()

# create croping area to inside of switzerland
cropping_sf <-
  read_sf(
    dsn = paste0(dd_ch, "hydrographische-gliederungderschweiz/Hydrografische+Gliederung/Hydrografische Gliederung_LV95"),
    layer = "basis04"
  ) %>%
  st_union() %>%
  st_transform(., crs = target_crs)

# get intersection between rhine subcatchments and switzerland borders
subcatchments_rhine_v2 <- st_intersection(subcatchments_rhine, cropping_sf)
subcatchments_rhine    <- subcatchments_rhine_v2

# remove lakes from objects
subcatchments_final <- st_difference(subcatchments_rhine, lakes)
river_intersect_lakes <- st_intersection(river_intersect_lakes, cropping_sf)

## extract environmental data
# extract polygon values over environmental data
all_env_subcatchments <- terra::extract(all_env_rast, 
                                        terra::vect(subcatchments_final), 
                                        fun = function(x) mean(x, na.rm = T))

all_env_subcatchments$TEILEZGNR <- subcatchments_final$TEILEZGNR
all_env_subcatchments$ID <- NULL



# read in swiss boundaries
ch <- st_read(paste0(dd, 'swissboundaries3d_2021-07_2056_5728.shp/SHAPEFILE_LV95_LN02'), layer = 'swissBOUNDARIES3D_1_3_TLM_LANDESGEBIET')
ch <- st_transform(st_union(ch), crs = target_crs)

