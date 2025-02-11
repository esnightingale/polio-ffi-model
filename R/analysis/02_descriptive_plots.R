################################################################################
################################################################################
# Descriptive plots of observed WPV1+ AFP cases and environmental samples
################################################################################

source("R/utils/setup_env.R")
source("R/utils/mapping_fcns.R")
cc_name <- tolower(paste0(cc$name, collapse = "_"))

# Virus type to analyse
type = cc$type

# Specify genetic cluster to analyse
clust = cc$cluster

# Specify emergence group to analyse
eg = cc$eg

# Set up directories
# Inputs (incl. all serotypes)
indir <- paste0("inputs/",cc_name)

# Outputs
if(!is.null(clust)){outdir <- here("results",cc_name,type,clust,substr(cc$start,1,4), "figures/descriptive")
}else if (!is.null(eg)){outdir <- here("results",cc_name,type,eg,substr(cc$start,1,4), "figures/descriptive")
}else {outdir <- here("results",cc_name,type,substr(cc$start,1,4), "figures/descriptive")}
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

################################################################################
# Load data

# District data
data <- readRDS(here(indir, "data.rds")) %>%
  # Retain only for specified serotype
  filter(serotype == gsub("[^0-9.-]", "", type))

# Radiation matrix
R <- readRDS(here(indir, "R.rds"))

# AFP cases
afp <- readRDS(here(indir, "afp.rds"))

# ES
es <- readRDS(here(indir, "es.rds"))

# Create extra field in ES for negative samples
es$final_neg <- 1
es$final_neg[es$final_class == type] <- 0
table(es$final_neg,es$final_class)

# District shapes
shape0 <- readRDS(here(indir, "shape0.rds"))
shape1 <- readRDS(here(indir, "shape1.rds"))
shape2 <- readRDS(here(indir, "shape2.rds"))

# ---------------------------------------------------------------------------- #

## Flag relevant samples/cases
if(!is.null(clust)){
  es$flag <- replace_na(grepl(clust, es$virus_cluster_s), FALSE)
  afp$flag <- replace_na(grepl(clust, afp$virus_clusters), FALSE)
}else if(!is.null(eg)){
  es$flag <- replace_na(grepl(eg, es$emergence_group_s), FALSE)
  afp$flag <- replace_na(grepl(eg, afp$emergence_group_s), FALSE)
}else{
  es$flag <- replace_na(es$final_class == type, FALSE)
  afp$flag <- replace_na(afp$final_class == type, FALSE)
}

# ---------------------------------------------------------------------------- #
# Check differences in population estimates

# if (cc_name =="nga"){
#   data %>%
#     ggplot(aes(total_pop, pop_hdx)) +
#     geom_abline(lty = "dashed") +
#     geom_point(alpha = 0.2) +
#     scale_x_continuous(trans = "log10") +
#     scale_y_continuous(trans = "log10") -> plot_pops_total
#
#   data %>%
#     ggplot(aes(population_u5, pop_u5)) +
#     geom_abline(lty = "dashed") +
#     geom_point(alpha = 0.2) +
#     scale_x_continuous(trans = "log10") +
#     scale_y_continuous(trans = "log10") -> plot_pops_u5
#
#   data %>%
#     ggplot(aes(pop_u15, pop_u15_polisdenom)) +
#     geom_abline(lty = "dashed") +
#     geom_point(alpha = 0.2) +
#     scale_x_continuous(trans = "log10") +
#     scale_y_continuous(trans = "log10") -> plot_pops_u15
#
#   plot_pops_total + plot_pops_u5 + plot_pops_u15
#   ggsave(here(outdir, "check_hdx_pops.png"),
#          height = 5, width = 15)
# }

# ---------------------------------------------------------------------------- #
# Time series

date_first <- date_to_period(cc$start)-1
date_last <- date_to_period(cc$end)

m = 3

tmp_afp <- afp %>%
  mutate(period = date_to_period(date)) %>%
  group_by(month, period) %>%
  summarise(afp_notifications = n(),
            cases = sum(flag, na.rm = T)) %>%
  ungroup() %>%
  mutate(date = period_to_date(period))

tmp_es <- es %>%
  mutate(period = date_to_period(collection_date)) %>%
  group_by(month, period) %>%
  summarise(es_samples = n(),
            es_positives = sum(flag, na.rm = T)) %>%
  ungroup() %>%
  mutate(date = period_to_date(period))

tmp_es$date <- factor(tmp_es$date, levels = unique(tmp_es$date))
tmp_afp$date <- factor(tmp_afp$date, levels = unique(tmp_afp$date))

msamples <- max(tmp_es$es_samples)
mafp <- max(tmp_afp$afp_notifications)

coeff1 <- round(mafp/max(tmp_afp$cases))
coeff2 <- round(msamples/max(tmp_es$es_positives))

desc <- case_when(!is.null(clust) ~ paste(type,clust,sep = ", "),
                  !is.null(eg) ~ paste(type,eg,sep = ", "),
                  is.null(clust) & is.null(eg) ~ type)

pal <- viridis::turbo(2, begin = 0.2, end = 0.8)
p1 <- ggplot(tmp_afp %>%
               filter(period >= date_first & period <= date_last),
             aes(x=date,cases)) +
  geom_bar(stat = "identity", fill = pal[2]) +
  geom_point(aes(period_to_date(period), afp_notifications/coeff1), pch = 22, col = pal[2], fill = "white") +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1,
                                 color=rep(c("black", rep("transparent", each =m-1)), n_distinct(tmp_afp$month)))) +
  scale_y_continuous(
    name = "Polio cases",
    sec.axis = sec_axis(~.*coeff1, name="Total AFP notifications")
  ) +
  labs(x = NULL)

p2 <- ggplot(tmp_es %>%
               filter(period >= date_first & period <= date_last),
             aes(x=date,y=es_positives)) +
  geom_bar(stat = "identity",fill=pal[1]) +
  geom_point(aes(x=period_to_date(period),y=es_samples/coeff2),pch = 22, col= pal[1], fill = "white") +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1,
                                 color=rep(c("black", rep("transparent", each =m-1)), n_distinct(tmp_es$month)))) +
  scale_y_continuous(
    name = "Positive ES samples",
    sec.axis = sec_axis(~.*coeff2, name="Total ES samples")
  ) +
  labs(x = NULL,
       caption = desc)

pdf(here(outdir,"time_series_plot.pdf"),height=5,width=7)
print(p1 / p2)
dev.off()

p <- p1 / p2
ggsave(here(outdir,"time_series_plot.png"),p,height=5,width=7)

# ---------------------------------------------------------------------------- #
# Samples collected per site per m

es %>%
  filter(date_to_period(collection_date) >= date_to_period(cc$start) &
           date_to_period(collection_date) <= date_to_period(cc$end)) %>%
  group_by(month, site_id) %>%
  summarise(n = n()) %>%
  ungroup() -> es_site_mth

es_site_mth %>%
  ggplot(aes(month, fill = as.factor(n))) +
  geom_bar() +
  scale_x_date(date_labels = "%Y-%b") +
  scale_fill_viridis_d(option = "turbo", begin = 0.1, end = 0.9) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = "Number of sites",fill = "Samples collected",
       title = "Samples collected per site per m")

ggsave(here(outdir,"sampling_frequency.png"), height=5,width=10)

# ---------------------------------------------------------------------------- #
# Construct maps

# Subset now to specific time period for mapping
data_sub <- data %>%
  filter(between(period, date_first, date_last))

tdy <- ymd(cc$end)

# Define plot area/dimensions
bbox = data_sub %>%
  left_join(shape2,by = 'guid') %>%
  st_as_sf() %>%
  st_bbox()
aspect_ratio = (bbox$xmax-bbox$xmin)/(bbox$ymax-bbox$ymin)
base_width = 6

# -------------------------------- #
# Plots

## Province label map
ggplot() +
  geom_text_repel(aes(label = adm1_name,x=center_lon,y=center_lat),
                  vjust =0,hjust=0,data = shape1) +
  map_theme
ggsave(here(outdir,
            paste0(cc_name,'_province_names.png')),
       width = 1+aspect_ratio*base_width,height=1+base_width,units='in',dpi=600)

## Clinical cases
ggplot() +
  geom_point(aes(x=x,y=y),
             data = afp %>% filter(flag,
                                   between(date_to_period(date),date_first,date_last)),
             col = 'red',shape = 16,size=3,inherit.aes = FALSE) +
  # geom_text_repel(aes(label = adm1_name,x=center_lon,y=center_lat),
  #                 vjust =0,hjust=0,data = filter(shape1, adm1_name == "Borno")) +
  # labs(caption = paste0(period_to_date(date_first),' to ',
                        # period_to_date(date_last))) +
  map_theme +
  labs(caption = desc) -> map_cases

map_cases
ggsave(here(outdir,
            paste(cc_name,type,'case.png', sep = "_")),
       map_cases,
       width = 5, height = 4)

## Cases/ES positives - past 12m
ggplot() +
  geom_sf(data = shape2, fill = '#00000000')+
  map_theme +
  case_points(from = tdy%m-%ms(12), labels=F)+
  ggtitle(paste0(type,' Cases and ES+ last 12 ms'),
          subtitle = paste0(format(tdy%m-%ms(12),'%d %b %Y'),
                            ' to ',
                            format(tdy,'%d %b %Y'))) +
  labs(caption = desc)

ggsave(here(outdir,
            paste(cc_name,
                   type,
                   'case_6mth.png', sep = "_")),
       width = 1+aspect_ratio*base_width,height=1+base_width,units='in',dpi=600)

## With ES negatives
ggplot() +
  geom_sf(data = shape2, fill = '#00000000')+
  map_theme +
  case_points(tdy%m-%ms(6), labels=F)+
  esneg_points(tdy%m-%ms(6),labels=F) +
  ggtitle(paste0(type,' Cases and ES+ last 6 ms'),
          subtitle = paste0(format(tdy%m-%ms(6),'%d %b %Y'),
                            ' to ',
                            format(tdy,'%d %b %Y'))) +
  labs(caption = desc)

ggsave(here(outdir,
            paste(cc_name,
                   type,
                   'case_6mth_wneg.png', sep = "_")),
       width = 1+aspect_ratio*base_width,height=1+base_width,units='in',dpi=600)

## Modelled immunity (proportion)
data_sub %>%
  mutate(period = round(period, 3)) %>%
  # Immunity one m prior to period of interest
  filter(period == round(date_to_period(cc[3])-(1/12),3)) %>%
  left_join(shape2,by = 'guid') %>% st_as_sf() %>%
  ggplot(aes(fill= immunity)) + geom_sf() +
  scale_fill_viridis_c('Immunity') + #limits = c(0,1)
  map_theme +
  ggtitle(paste0('Estimated ',type, ' Immunity ',
                 period_to_date(data_sub$period[1],'%b-%Y')))
ggsave(here(outdir,paste(cc_name,type,'immunity.png',sep = "_")),
       width = 1+aspect_ratio*base_width,height=base_width,units='in',dpi=600)

################################################################################
################################################################################
