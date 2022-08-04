## script to perform spatial operations on swiss data that define the domains of the xAI analysis

#### Set standard crs ----

target_crs <- 'epsg:3035'

#### Switzerland borders ----

# catchment for whole system
cropping_sf <- 
  read_sf(dsn = paste0(dd, 'hydrographische-gliederungderschweiz/Hydrografische+Gliederung/Hydrografische Gliederung_LV95'),
          layer = 'basis04') %>% 
  st_union()

# lake cropping
inverse_cropping_sf <- 
  read_sf(dsn = paste0(dd, "swiss-tlmregio-2021-2056/FGDB/swissTLMRegio_Product_LV95.gdb"),
          layer = "TLMRegio_Lake") %>% 
  filter(SHAPE_Area > 1000000)

# combine the lakes that border switzerland into the catchment area
#inverse_cropping_sf[order(inverse_cropping_sf$SHAPE_Area, decreasing = T),][c(1,2,5),1] %>% plot
lakes <- st_make_valid(st_union(inverse_cropping_sf[order(inverse_cropping_sf$SHAPE_Area, decreasing = T),][c(1,2,5),1]))
cropping_sf <- st_union(c(cropping_sf, lakes))
cropping_sf <- st_transform(cropping_sf, crs = target_crs)
cropping_sf <- st_make_valid(cropping_sf)

# get lakes inside ch
lakes_in_ch <- st_intersection(st_transform(inverse_cropping_sf, crs = target_crs), cropping_sf)

#### Subcatchments  ---- 

# load in 40km catchments
subcatchment_file <- paste0(dd, 'swiss-2km-subcatchments/EZG_Gewaesser.gdb')
subcatchments <- st_read(subcatchment_file, layer = 'Teileinzugsgebiet')
subcatchments <- st_transform(subcatchments, crs = target_crs)
subcatchments <- st_cast(subcatchments, 'POLYGON') %>% st_zm()

# union the subcatchments into 40km catchments
subcatchments_40kms <- subcatchments %>%
  mutate(TEZGNR40 = as.character(TEZGNR40)) %>% 
  group_by(TEZGNR40) %>% 
  dplyr::summarize(Shape = st_union(Shape))

# intersect with ch borders
subcatchments_final <- st_intersection(subcatchments_40kms, cropping_sf)


#### River network ----


streams_path <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/eu-hydro/CH_HYDRO.shp'
ch_rivers <- st_read(streams_path)


