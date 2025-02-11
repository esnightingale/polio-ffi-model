# Calculate district coverage of active sites
get_es_coverage <- function(data, # current district dataset
                            es_linelist, # sample linelist for period of interest
                            shapes, # district shapefiles
                            poprast, # population raster
                            proj, # Projection of raster/shapes
                            buffer_km, # chosen buffer for site coverage
                            return_catchment_shapes # option to return catchment areas by district
                            ){

  # Define district centroid if coords missing
  es_linelist %>% 
    left_join(st_drop_geometry(shapes)) %>%
    mutate(x = coalesce(x, center_lon),
           y = coalesce(y, center_lat)) -> es_linelist
  
  # Aggregate samples to unique sites, define sf and Xkm buffer
  site_buff <- es_linelist %>%
    dplyr::select(guid, site_id, #site_class, 
                  coord_imp, x, y) %>%
    dplyr::distinct() %>% 
    sf::st_as_sf(coords = c("x","y"), crs = st_crs(4326), remove = F) %>% 
    # Project to match raster/shapes
    sf::st_transform(proj) %>% 
    sf::st_buffer(dist = set_units(buffer_km, "km"))
  
  print(paste0(nrow(site_buff), " unique sites"))
  
  # browser()
  
  if(nrow(site_buff) > 0){
    
  # Intersect site buffer and district polygons
  intersections <- sf::st_intersects(x = site_buff, y = shapes) 
  
  pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", 
                                   total = dim(site_buff)[1])
  
  site_dist <- purrr::map_dfr(1:dim(site_buff)[1], 
                              function(ix){
                                pb$tick()
                                sf::st_intersection(x = site_buff[ix,], 
                                                    y = shapes[intersections[[ix]],])
                                })

  print(paste0(nrow(site_dist), " segments by district"))
  print(paste0(n_distinct(site_dist$guid), " districts with coverage"))
  
  # site_dist %>%  group_by(guid) %>% tally() %>% View
  
  # Check validity of segments
  if(any(!st_is_valid(site_dist))){
    print("Validity of segments:")
    print(table(st_is_valid(site_dist)))
    print("Making valid...")
    site_dist <- st_make_valid(site_dist)
  }
  
  # Check geometry types
  # print("Geometry types:")
  # print(table(sf::st_geometry_type(site_dist)))
  # Remove any non-polygonal segments
  if(any(!sf::st_geometry_type(site_dist) %in% c("POLYGON","MULTIPOLYGON"))){
    print("Removing non-polygon geometries...")
    site_dist <- site_dist[grepl("POLYGON", sf::st_geometry_type(site_dist)),]
  }

  # browser()
  
  # Filter to segments where sample district equals coverage district, then union
  by_dist <- site_dist %>% 
    filter(guid == guid.1) %>% 
    dplyr::group_by(guid.1) %>%   
    dplyr::summarise(geometry = sf::st_union(geometry),  
                     total_pop = unique(total_pop)) %>% 
    dplyr::ungroup()
  
  print(paste0(nrow(by_dist), " districts with coverage"))
  
  # print(table(sf::st_geometry_type(by_dist)))
  # Remove any non-polygonal segments
  if(any(!sf::st_geometry_type(by_dist) %in% c("POLYGON","MULTIPOLYGON"))){
    print("Removing non-polygon geometries...")
    by_dist <- by_dist[grepl("POLYGON", sf::st_geometry_type(by_dist)),]
  }

  # Extract the population within each of these new district areas covered by sites
  by_dist$catchment_pop <- exactextractr::exact_extract(poprast,
                                                        by_dist, 
                                                        fun = "sum")  
  
  dist_covg <- by_dist %>% 
    rowwise() %>% 
    dplyr::mutate(es_coverage = min(catchment_pop/total_pop, 0.999)) %>% 
    ungroup() %>% 
    sf::st_drop_geometry() %>% 
    select(-total_pop)
  
  }else{
    # If no coverage at all, define empty datasets
    site_dist = NULL
    dist_covg <- data.frame(guid.1 = character(0), catchment_pop = numeric(0), es_coverage = numeric(0))
    print("Zero districts with coverage")
  }
  
  data <- dplyr::left_join(data, 
                           dist_covg, 
                           by = c("guid" = "guid.1")) %>%
    dplyr::mutate(across(catchment_pop:es_coverage, replace_na, 0))
  
  if (return_catchment_shapes){
    return(list(data = data, 
                n_sites = nrow(site_buff), 
                n_samples = nrow(es_linelist),
                catchment_shapes = site_dist))
  }else{
    return(list(data = data, 
                n_sites = nrow(site_buff),
                n_samples = nrow(es_linelist),
                catchment_shapes = NULL))
  }

}
