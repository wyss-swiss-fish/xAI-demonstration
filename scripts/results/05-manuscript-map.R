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

# testing taking the 2km2 subcatchments containig the main stems
sense <- data.frame(catchment = 'Sense', 
                    TEILEZGNR = c(54206, 49238, 36039, 53576, 92800, 
                                  73564, 79104, 55745, 61040, 21140, 
                                  76158, 40046, 95372, 78251, 21268, 
                                  96970, 1492, 64847, 85482, 21499, 
                                  5407, 72778, 98155, 18709, 9298))
emme <- data.frame(catchment = 'Emme', 
                   TEILEZGNR = c(48935, 519, 76613, 38042, 23267,64482, 
                                 17698, 39457, 59083, 8682, 42744, 
                                 14206, 88728, 80592, 44084, 92441,
                                 62780, 47744, 19469, 80840, 66777, 
                                 54624, 54833, 26718, 81132, 31679, 
                                 4435, 94151))
catchments <- rbind(sense, emme)


subcatchment_file <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/swiss-2km-subcatchments/EZG_Gewaesser.gdb'

# read in subcatchments, transform, union
emme_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEILEZGNR %in% emme$TEILEZGNR) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm() %>% 
  st_as_sf() %>% 
  mutate(name = 'Emme')

sense_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEILEZGNR %in% sense$TEILEZGNR) %>%
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
  tm_polygons(col = 'red', legend.show = F, border.col = 'red', lwd = 5) + 
  tm_shape(sense_union) + 
  tm_polygons(col = 'red', legend.show = F, border.col = 'red', lwd = 5) + 
  tm_shape(rivers_rhein %>% arrange(STRAHLE)) + 
  tm_lines(lwd = 'STRAHLE', scale = 3, legend.lwd.show = F) + 
  tm_shape(lakes) + 
  tm_polygons(col = 'lightblue', border.col = 'black') + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()  




#### read in environmental data and extract values for emme and sense ----

# set vector of focal variables for assigning new names, this must be consistent with the models fitted
vars <- c('ecoF_discharge_max_log10', 
          'ecoF_slope_min_log10', 
          'ecoF_flow_velocity_mean', 
          'stars_t_mx_m_c', 
          'stars_t_mn_m_c',
          'local_asym_cl_log10',
          'local_dis2lake',
          'ecoF_eco_mean_ele_residual', 
          'local_imd_log10_ele_residual',
          'local_wet',
          'local_flood')

vars_renamed = c('discharge', 
                 'slope', 
                 'flow velocity', 
                 'temperature max', 
                 'temperature min', 
                 'connectivity', 
                 'distance to lake', 
                 'ecomorphology', 
                 'urbanisation',
                 'wetland', 
                 'floodplains')

# reset environmental names
vars_renamed_V2 <- cbind(vars, vars_renamed)

emme_env <- extract(env_data, vect(emme_union)) %>% 
  select(all_of(vars)) %>% 
  pivot_longer(colnames(.)) %>% 
  mutate(catchment = 'Emme')

sense_env <- extract(env_data, vect(sense_union)) %>% 
  select(all_of(vars)) %>% 
  pivot_longer(colnames(.)) %>% 
  mutate(catchment = 'Sense')

catchment_envs <- rbind(emme_env, sense_env)

# organise variable names
catchment_envs <- left_join(catchment_envs, data.frame(vars_renamed_V2), by = c('name' = 'vars'))

pdf(paste0(force_dir, '/emme_sense_env.pdf'))
ggplot(data = catchment_envs) + 
  geom_boxplot(aes(x = catchment, y = value, col = catchment)) + 
  facet_wrap(~vars_renamed, scale = 'free') + 
  theme_bw() + 
  theme(panel.grid = element_blank())
dev.off()




#### Make figure of raw environmental data ----

# ensure only catchments of interest
env_data_mask <- trim(mask(env_data, vect(subcatchments_rhine))) 

# take only variables of interest
env_data_mask <- env_data_mask[[vars]]

names(env_data_mask) <- vars_renamed_V2[,'vars_renamed'][which(names(env_data_mask) %in% vars_renamed_V2[,'vars'])]

pdf('figures/env_data.pdf', width = 10, height = 7)
tm_shape(env_data_mask) + 
  tm_raster(midpoint = NA, 
            style = 'cont', 
            palette = 'viridis', 
            title = '') + 
  tm_facets(ncol = 4, 
            free.scales.raster = T) + 
  tm_layout(frame = F)
dev.off()









