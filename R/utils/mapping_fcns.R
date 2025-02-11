################################################################################
# Mapping functions
################################################################################

# Define a mapping theme
map_theme = list(
  theme_void(),  # Empty theme without axis lines and texts
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )
)

# Function to add observed cVDPV2 cases and positive ES as labelled points
case_points = function(afp = afp, es = es, from,to=tdy,type = "WPV1", labels=T){
  pts = list(
    geom_point(aes(x=x,y=y),
               data = afp %>% filter(final_class == type,
                                       between(date,from,to)),
               fill = 'red',shape = 21,size=5,inherit.aes = FALSE),
    geom_point(aes(x=x,y=y),
               data = es %>% filter(final_class == type,
                                    between(collection_date,from,to)),
               fill = 'orange',shape = 22,size=5,inherit.aes = FALSE)
  )
  labs = list(
    geom_text_repel(aes(x=x,y=y,
                        label = paste0(lubridate::month(donset,label=T),day(donset))),
                    data = afp %>% filter(final_class == type,
                                            between(date,from,to)),
                    inherit.aes = FALSE),
    geom_text_repel(aes(x=x,y=y,
                        label = paste0(lubridate::month(collection_date,label=T),day(collection_date))),
                    data = es %>% filter(final_class == type,
                                         between(collection_date,from,to)) %>% group_by(guid) %>% arrange(desc(collection_date)) %>% slice(1),
                    inherit.aes = FALSE)
  )
  if(labels) return(c(pts,labs))
  return(pts)
}

# Function to add negative ES as labelled points
esneg_points = function(es = es, from,to=tdy,labels=T){

  pts = list(
    # if NA then assume negative
    # actually create new var that encompassess all the options
    geom_point(aes(x=x,y=y),
               data = es %>% filter(!flag,
                                    between(collection_date,from,to)),
               fill = 'grey99',shape = 22,size=2,inherit.aes = FALSE)
  )

  labs = list(
    geom_text_repel(aes(x=x,y=y,
                        label = paste0(lubridate::month(collection_date,label=T),day(collection_date))),
                    data = es %>% filter(flag,
                                         between(collection_date,from,to)) %>% group_by(guid) %>% arrange(desc(collection_date)) %>% slice(1),
                    inherit.aes = FALSE)
  )
  #browser()
  if(labels) return(c(pts,labs))
  return(pts)
}

# Plot relative risk from last 12m detections and estimated immunity
plot_risk <- function(out_1){

  risk_period <- ymd(out_1$risk$risk_period)
  es_1 <- out_1$risk$es
  afp_1 <- out_1$risk$afp

  out_1$sens$data %>%
    left_join(st_transform(shape2,4326) ,by = 'guid') %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill= risk_est_sc)) +
    scale_fill_viridis_c("Scaled Risk", limits = c(0,1)) + # ,trans = 'log'
    map_theme +
    # # Negative AFP
    # geom_point(aes(x=x,y=y),
    #            data = afp_1 %>% filter(!flag),
    #            col = 'white',shape = 4, size=0.5,inherit.aes = FALSE) +
    # Positive samples
    geom_jitter(aes(x=x,y=y, shape = "Env. sample"),
                data = es_1 %>% filter(flag),
                height = 0.01,
                width = 0.01,
                fill = "red", size=4,inherit.aes = FALSE) +

    # Positive AFP
    geom_jitter(aes(x=x,y=y, shape = "AFP case"),
                data = afp_1 %>% filter(flag),
                fill = "red",
                height = 0.01,
                width = 0.01,
                size=4,inherit.aes = FALSE) +
    # Negative samples
    geom_jitter(aes(x=x,y=y),
                data = es_1 %>% filter(!flag),
                height = 0.01,
                width = 0.01,
                fill = 'white',shape = 22,size=2,inherit.aes = FALSE) +
    scale_shape_manual(values = c(21,22),
                       labels = c("AFP case",
                                  "Env. sample")) +
    guides(fill = "none", col = "none", shape = "none") +
    labs(title = "WPV1 Circulation Risk",
         subtitle = paste(risk_period[1], "to", risk_period[2]),
         shape = NULL,
         caption = "Negative environmental samples shown in white") -> p_prior_risk
  p_prior_risk

  ggsave(here(figdir,"gif","risk",paste0("prior_risk_",risk_period[2],".png")), p_prior_risk, height=6,width=7)

}

# Map ENV samples collected
map_samps_env <- function(out_1, pts = F, legend = T){

  sens_period <- ymd(out_1$sens$sens_period)
  # es_1 <- out_1$risk$es
  # afp_1 <- out_1$risk$afp

  out_1$sens$data %>%
    left_join(shape2,by = 'guid') %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill= es_n)) +
    scale_fill_viridis_c(option = "rocket", "No. samples", limits = c(0,80), trans = "sqrt",
                         direction = -1) + # ,trans = 'log'
    map_theme +
    theme(legend.position = c(0.9,0.25)) +
    # guides(fill = "none", col = "none", shape = "none") +
    labs(title = "Environmental surveillance - Samples collected per district",
         subtitle = paste(sens_period[1], "to", sens_period[2]),
         shape = NULL) -> p_map_samps_env
  p_map_samps_env

  if(legend == F){
    p_map_samps_env <- p_map_samps_env +
      guides(fill = "none", col = "none")
  }

  ggsave(here(figdir,"gif","env",paste0("samps_env_",sens_period[2],"_nopts.png")), p_map_samps_env, height=6,width=7)

}

# Map estimated ENV sensitivity
map_sens_env <- function(out_1, pts = F, legend = T){

  sens_period <- ymd(out_1$sens$sens_period)
  # es_1 <- out_1$risk$es
  # afp_1 <- out_1$risk$afp

  out_1$sens$data %>%
    left_join(shape2,by = 'guid') %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill= SeR_ESi)) +
    scale_fill_viridis_c(option = "rocket", "Sensitivity", limits = c(0,1), trans = "sqrt",
                         direction = -1) + # ,trans = 'log'
    map_theme +
    theme(legend.position = c(0.9,0.25)) +
    scale_shape_manual(values = 22,
                       labels = "Env. sample") +
    # guides(fill = "none", col = "none", shape = "none") +
    labs(title = "Environmental surveillance sensitivity",
         subtitle = paste(sens_period[1], "to", sens_period[2]),
         shape = NULL) -> p_map_sens_env
  p_map_sens_env

  if(legend == F){
    p_map_sens_env <- p_map_sens_env +
      guides(fill = "none", col = "none")
  }

  if(pts == T){
    # p_map_sens_env <- p_map_sens_env +
    #     # All samples
    #     geom_jitter(aes(x=x,y=y, shape = "Env. sample"),
    #                 data = es_1,
    #                 height = 0.01,
    #                 width = 0.01,
    #                 fill = "white",
    #                 shape = 22,
    #                 size=2,inherit.aes = FALSE) +
    # ggsave(here(figdir,"gif","env",paste0("sens_env_",sens_period[2],".png")), p_map_sens_aenv, height=6,width=7)
  }else{
    ggsave(here(figdir,"gif","env",paste0("sens_env_",sens_period[2],"_nopts.png")), p_map_sens_env, height=6,width=7)
  }

}

# Map AFP reported
map_samps_afp <- function(out_1, pts = F, legend = T){

  sens_period <- ymd(out_1$sens$sens_period)
  # es_1 <- out_1$risk$es
  # afp_1 <- out_1$risk$afp

  out_1$sens$data %>%
    mutate(afp_r = afp_n*1e5/(0.45*total_pop)) %>%
    left_join(shape2,by = 'guid') %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = afp_r)) +
    scale_fill_viridis_c(option = "rocket", "Rate", limits = c(0,410), trans = "sqrt",
                         direction = -1) + # ,trans = 'log'
    map_theme +
    theme(legend.position = c(0.9,0.25)) +
    labs(title = "AFP surveillance - cases reported per district, per 100,000 U15s",
         subtitle = paste(sens_period[1], "to", sens_period[2]),
         shape = NULL) -> p_map_samps_afp
  # p_map_sens_afp

  if(legend == F){
    p_map_samps_afp <- p_map_samps_afp +
      guides(fill = "none", col = "none")
  }

  if(pts == T){
    p_map_samps_afp <- p_map_samps_afp +
      # Negative AFP
      # geom_point(aes(x=x,y=y),
      #            data = afp_1 %>% filter(!flag),
      #            col = 'white',shape = 4, size=0.5,inherit.aes = FALSE) +
      # All AFP
      geom_jitter(aes(x=x,y=y, shape = "AFP case"),
                  data = afp_1,
                  height = 0.01,
                  width = 0.01,
                  fill = "white",
                  col = "white",
                  shape = 4,
                  alpha = 0.3,
                  size=0.3,inherit.aes = FALSE)
    ggsave(here(figdir,"gif","afp",paste0("samps_afp_",sens_period[2],".png")), p_map_samps_afp, height=6,width=7)
  }else{
    ggsave(here(figdir,"gif","afp",paste0("samps_afp_",sens_period[2],"_nopts.png")), p_map_samps_afp, height=6,width=7)
  }

}

# Map estimated AFP sensitivity
map_sens_afp <- function(out_1, pts = F, legend = T){

  sens_period <- ymd(out_1$sens$sens_period)
  # es_1 <- out_1$risk$es
  # afp_1 <- out_1$risk$afp

  out_1$sens$data %>%
    left_join(shape2,by = 'guid') %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill= SeR_AFPi)) +
    scale_fill_viridis_c(option = "rocket", "Sensitivity", limits = c(0,1), trans = "sqrt",
                         direction = -1) + # ,trans = 'log'
    map_theme +
    theme(legend.position = c(0.9,0.25)) +
    labs(title = "AFP surveillance sensitivity",
         subtitle = paste(sens_period[1], "to", sens_period[2]),
         shape = NULL) -> p_map_sens_afp
  # p_map_sens_afp

  if(legend == F){
    p_map_sens_afp <- p_map_sens_afp +
      guides(fill = "none", col = "none")
  }

  if(pts == T){
    p_map_sens_afp <- p_map_sens_afp +
      # Negative AFP
      # geom_point(aes(x=x,y=y),
      #            data = afp_1 %>% filter(!flag),
      #            col = 'white',shape = 4, size=0.5,inherit.aes = FALSE) +
      # All AFP
      geom_jitter(aes(x=x,y=y, shape = "AFP case"),
                  data = afp_1,
                  height = 0.01,
                  width = 0.01,
                  fill = "white",
                  col = "white",
                  shape = 4,
                  alpha = 0.3,
                  size=0.3,inherit.aes = FALSE)
    ggsave(here(figdir,"gif","afp",paste0("sens_afp_",sens_period[2],".png")), p_map_sens_afp, height=6,width=7)
  }else{
    ggsave(here(figdir,"gif","afp",paste0("sens_afp_",sens_period[2],"_nopts.png")), p_map_sens_afp, height=6,width=7)
  }

}

