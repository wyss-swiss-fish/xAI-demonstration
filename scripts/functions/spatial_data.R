## script to perform spatial operations on swiss data that define the domains of the xAI analysis

#### Set standard crs ----

# read in elevation raster
base_rast <- rast(paste0(dd_env, "ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif"))
base_rast <- base_rast[[1]]

# set crs
target_crs <- 'epsg:3035'

# read in lakes
lakes <- paste0(dd_env, "ch-lake-reference.shp")
lakes <- st_read(lakes)

# read in rivers
rivers <- paste0(dd_env, "CH_HYDRO.shp")
rivers <- st_read(rivers)
rivers <- rivers %>%
  filter(CUM_LEN > 500, STRAHLE > 1) %>%
  mutate(CUM_LEN_LOG = log(.$CUM_LEN))

# read in catchment object
subcatchment_file <- paste0(dd_ch, 'swiss-2km-subcatchments/EZG_Gewaesser.gdb')
catchments_rhine <- c('Aare', 'Reuss', 'Limmat', 'Rhein')

# read in subcatchments, transform, union
subcatchments_rhine_union <- st_read(subcatchment_file, layer = 'Teileinzugsgebiet') %>% 
  filter(FLUSSGB %in% catchments_rhine) %>% 
  # convert to target crs
  st_transform(., crs = target_crs) %>% 
  # union together
  st_union() %>% 
  # remove z and m properties that can cause errors later
  st_zm()

# mask over rhein
base_rast  <- terra::trim(terra::mask(base_rast, vect(subcatchments_rhine_union)))
lakes      <- st_intersection(lakes, subcatchments_rhine_union)
rivers     <- st_intersection(rivers, subcatchments_rhine_union)
rivers_int <- st_union(rivers %>% dplyr::select(geometry))
river_intersect_lakes <- st_difference(rivers_int, st_union(lakes))

##### REMOVE AREAS OUTSIDE OF SWISS BORDERS

# get subcatchments but keep as seperate entities for aggregations
# read in subcatchments, transform, union
subcatchments_rhine <- st_read(subcatchment_file, layer = 'Teileinzugsgebiet') %>% 
  filter(FLUSSGB %in% catchments_rhine) %>% 
  # convert to target crs
  st_transform(., crs = target_crs) %>% 
  # remove z and m properties that can cause errors later
  st_zm() %>% 
  st_difference(., lakes)

cropping_sf <- 
  read_sf(dsn = paste0(dd_ch, 'hydrographische-gliederungderschweiz/Hydrografische+Gliederung/Hydrografische Gliederung_LV95'),
          layer = 'basis04') %>% 
  st_union() %>% 
  st_transform(., crs = target_crs)

# get intersection
subcatchments_rhine_v2 <- st_intersection(subcatchments_rhine, cropping_sf)
subcatchments_rhine <- subcatchments_rhine_v2

# remove lakes
subcatchments_final   <- st_difference(subcatchments_rhine, st_union(lakes))
river_intersect_lakes <- st_intersection(river_intersect_lakes, cropping_sf)


### #### Switzerland borders ----
### 
### # catchment for whole system
### cropping_sf <- 
###   read_sf(dsn = paste0(dd, 'hydrographische-gliederungderschweiz/Hydrografische+Gliederung/Hydrografische Gliederung_LV95'),
###           layer = 'basis04') %>% 
###   st_union()
### 
### # lake cropping
### inverse_cropping_sf <- 
###   read_sf(dsn = paste0(dd, "swiss-tlmregio-2021-2056/FGDB/swissTLMRegio_Product_LV95.gdb"),
###           layer = "TLMRegio_Lake") %>% 
###   filter(SHAPE_Area > 1000000)
### 
### # combine the lakes that border switzerland into the catchment area
### #inverse_cropping_sf[order(inverse_cropping_sf$SHAPE_Area, decreasing = T),][c(1,2,5),1] %>% plot
### lakes <- st_make_valid(st_union(inverse_cropping_sf[order(inverse_cropping_sf$SHAPE_Area, decreasing = T),][c(1,2,5),1]))
### cropping_sf <- st_union(c(cropping_sf, lakes))
### cropping_sf <- st_transform(cropping_sf, crs = target_crs)
### cropping_sf <- st_make_valid(cropping_sf)
### 
### # get lakes inside ch
### lakes_in_ch <- st_intersection(st_transform(inverse_cropping_sf, crs = target_crs), cropping_sf)
### 
### #### Subcatchments  ---- 
### 
### # load in 40km catchments
### subcatchment_file <- paste0(dd, 'swiss-2km-subcatchments/EZG_Gewaesser.gdb')
### subcatchments <- st_read(subcatchment_file, layer = 'Teileinzugsgebiet')
### subcatchments <- st_transform(subcatchments, crs = target_crs)
### subcatchments <- st_cast(subcatchments, 'POLYGON') %>% st_zm()
### 
### # union the subcatchments into 40km catchments
### subcatchments_40kms <- subcatchments %>%
###   mutate(TEZGNR40 = as.character(TEZGNR40)) %>% 
###   group_by(TEZGNR40) %>% 
###   dplyr::summarize(Shape = st_union(Shape))
### 
### # intersect with ch borders
### subcatchments_final <- st_intersection(subcatchments_40kms, cropping_sf)
### 
### 
### #### River network ----
### 
### streams_path <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/eu-hydro/CH_HYDRO.shp'
### ch_rivers <- st_read(streams_path)
### ch_rivers$CUM_LEN_LOG <- log(ch_rivers$CUM_LEN+1)
### 
### 
### ## remove lakes from rivers
### streams_in_lakes <- st_intersects(lakes_in_ch, ch_rivers)
### ch_rivers_2 <- ch_rivers[-unique(unlist(flatten(streams_in_lakes))),]
