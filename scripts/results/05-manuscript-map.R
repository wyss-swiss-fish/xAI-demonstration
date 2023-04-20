### load objects for creating better maps

# get paths to spatial objects
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"


# get the teilenzugsgebeit for the sense
sense <- data.frame(catchment = 'Sense', TEZGNR40 = c(100184, 101484, 105044, 105663, 102204, 102377))
emme <- data.frame(catchment = 'Emme', TEZGNR40 = c(106570, 106299, 101728, 101280, 106548, 103982, 108894))

#### swiss borders ----

# national borders
st_layers(paste0(dd, 'swissboundaries3d_2021-07_2056_5728.shp/SHAPEFILE_LV95_LN02'))
ch_poly <- st_read(paste0(dd, 'swissboundaries3d_2021-07_2056_5728.shp/SHAPEFILE_LV95_LN02'), 
        'swissBOUNDARIES3D_1_3_TLM_LANDESGEBIET')
ch_poly <- ch_poly %>% filter(NAME == 'Schweiz')
ch_poly <- st_union(ch_poly) %>%  st_transform(., crs = target_crs)

# subcatchment borders
subcatchment_file <- paste0(dd_ch, "swiss-2km-subcatchments/EZG_Gewaesser.gdb")
subcatchment <- st_read(subcatchment_file, layer = "Teileinzugsgebiet")
catchments_rhine <- c("Aare", "Reuss", "Limmat")

subcatchments_rhine_union <- subcatchment %>%
  filter(FLUSSGB %in% catchments_rhine) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm()

subcatchments_rhine_union <- st_intersection(ch_poly, subcatchments_rhine_union)

#### river networks ----

rivers <- paste0(dd_env, "CH_HYDRO.shp")
rivers <- st_read(rivers)
rivers <- rivers %>% st_transform(rivers, crs = target_crs)
rivers_rhein <- st_intersection(rivers, subcatchments_rhine_union)

#### lakes ----

lakes <- paste0(dd_env, "ch-lake-reference.shp")
lakes <- st_union(st_read(lakes))

st_union(lakes, ch_poly)

#### focal subcatchments ----

# get the teilenzugsgebeit for the sense and emme systems
sense <- data.frame(catchment = 'Sense', TEZGNR40 = c(100184, 101484, 105044, 105663, 102204, 102377))
emme <- data.frame(catchment = 'Emme', TEZGNR40 = c(106570, 106299, 101728, 101280, 106548, 103982, 108894))
catchments <- rbind(sense, emme)

subcatchment_file <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/swiss-2km-subcatchments/EZG_Gewaesser.gdb'

# read in subcatchments, transform, union
emme_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% emme$TEZGNR40) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm() %>% 
  st_as_sf() %>% 
  mutate(name = 'Emme')

sense_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% sense$TEZGNR40) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm() %>% 
  st_as_sf() %>% 
  mutate(name = 'Sense')

bbox_catch <- st_bbox(st_buffer(rbind(emme_union, sense_union), 10000))


#### urban areas ----

#### place names ----

#### elevation ----

tmap_mode('plot')

pdf('figures/map_of_catchments.pdf', width = 7, height = 5)
tm_shape(ch_poly) + 
  tm_borders(col='transparent') + 
  #tm_shape(subcatchments_rhine_union) + 
  #tm_polygons(border.col = 'transparent', fill = 'gray90')+
  tm_scale_bar(text.size = 1, position = c('left', 'bottom'), color.light = 'gray75') + 
  tm_shape(st_union(lakes, ch_poly)) + 
  tm_borders(col = 'black') + 
  tm_shape(emme_union) + 
  tm_polygons(border.col = 'red', col = 'gray80', legend.show = F, lwd = 2) + 
  tm_shape(sense_union) + 
  tm_polygons(border.col = 'red', col = 'gray80', legend.show = F, lwd = 2) + 
  tm_shape(rivers_rhein %>% arrange(STRAHLE)) + 
  tm_lines(lwd = 'STRAHLE', scale = 3, legend.lwd.show = F) + 
  tm_shape(lakes) + 
  tm_polygons(col = 'lightblue', border.col = 'black') + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()  





