################################################################################
# Plot estimated risks and surveillance sensitivities
################################################################################

# Set region and time period of interest
source("R/utils/setup_env.R")

# Other functions
source("R/utils/mapping_fcns.R")
source("R/utils/general_fcns.R")

################################################################################
# Read data / output

cc_name <- tolower(paste0(cc$name, collapse = "_"))

# Virus type to analyse
type = cc$type

# Specify genetic cluster to analyse
clust = cc$cluster

# Specify emergence group to analyse
eg = cc$eg

# Results directory
if(!is.na(clust)){indir <- here("results",cc_name,type,clust,substr(cc$start,1,4))
}else if (!is.na(eg)){indir <- here("results",cc_name,type,eg,substr(cc$start,1,4))
}else {indir <- here("results",cc_name,type,substr(cc$start,1,4))}

# Figure directory
figdir <- here(indir,"figures")

# Extracted outputs
data_all <- readRDS(here(indir, "data_all.rds"))
afp_all <- readRDS(here(indir, "afp_all.rds"))
afp_not_all <- readRDS(here(indir, "afp_not_all.rds"))
es_all <- readRDS(here(indir, "es_all.rds"))
es_npev_all <- readRDS(here(indir, "es_npev_all.rds"))
sens_all <- readRDS(here(indir, "sens_all.rds"))
sens_n_es <- readRDS(here(indir, "sens_n_es.rds"))
sens_env_wcovg <- readRDS(here(indir, "sens_env_wcovg.rds"))
sens_n_afp <- readRDS(here(indir, "sens_n_afp.rds"))
sens_n_afp_dist <- readRDS(here(indir, "sens_n_afp_dist.rds"))

# District shapes
shape2 <- readRDS(here("inputs",cc_name, "shape2.rds"))

tp <- unique(data_all$tdy)
tp_f <- format(ymd(tp), "%b-%Y")

################################################################################
# 12m risk

# Rolling 12m observed positives
data_all %>%
  pivot_longer(c("afp_pos","es_pos")) %>%
  ggplot(aes(tdy_f, value, fill = name)) +
  geom_col() +
  scale_x_discrete(drop = F) +
  theme_minimal() +
  theme(legend.position = c(0.8,0.8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = NULL,y = "Count",fill = NULL,
       title = "WPV1+ observations",
       subtitle = "Previous 12 months") -> p_risk_obs

# Rolling 12m risk
data_all %>%
  ggplot(aes(tdy_f, risk_est)) +
  # geom_boxplot() +
  geom_violin() +
  scale_x_discrete(drop = F) +
  scale_y_continuous(trans = "identity",
                     labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(legend.position = c(0.2,0.8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = NULL,y = "Estimated risk", fill = NULL,
       title = "Estimated risk",
       subtitle = "Based on observations during previous 12 months") -> p_risk_est

fig1 <- p_risk_obs + p_risk_est
fig1

ggsave(here(figdir,"risk_est.png"), fig1, height=5,width=10)

################################################################################
# Surveillance activity

## ES ##

# ES activity for each time point
sens_n_es %>%
  pivot_longer(-c(tdy,tdy_f,samples_per_site)) %>%
  ggplot(aes(ymd(tdy), value, group = name, col = name)) +
  geom_line() +
  scale_colour_viridis_d(option = "viridis", begin = 0.2, end = 0.8, direction = -1) +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  # theme_minimal() +
  theme(legend.position = c(0.2,0.8),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) +
  labs(x = NULL,y = "Count", col = NULL,
       title = "Surveillance activity per month",
       subtitle = "Number of active sites and number of collected samples") -> p_es_obs
p_es_obs

sens_n_es %>%
  ggplot(aes(ymd(tdy), samples_per_site, group = 1)) +
  geom_line() +
  geom_hline(aes(yintercept = mean(samples_per_site)), lty = "dashed", col = "red") +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  scale_y_continuous(limits = c(0,2.5)) +
  # theme_minimal() +
  theme(legend.position = c(0.2,0.8),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        ) +
  labs(x = NULL,y = "Mean per active site",
       title = NULL,
       subtitle = "Average samples collected per active site") -> p_avg_es_obs

fig2 <- p_es_obs + p_avg_es_obs
fig2

ggsave(here(figdir,"es_activity.png"), fig2, height=6,width=12)


## AFP ##

sens_n_afp %>%
  ggplot(aes(ymd(tdy), afp_r)) +
  geom_line() +
  geom_hline(aes(yintercept = 2), lty = "dashed", col = "red") +
  geom_hline(aes(yintercept = mean(afp_r)), lty = "dashed", col = "darkgrey") +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  scale_y_continuous(trans = "log10") +
  # theme_minimal() +
  theme(legend.position = c(0.2,0.8),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  labs(x = NULL,y = "AFP rate per 100,000",
       title = NULL,
       subtitle = "Observed rate of AFP reported\ncompared to expected rate of 2 per 100,000 (red)") -> p_afp_all

p_afp_all

ggsave(here(figdir,"afp_activity.png"), p_afp_all, height=5,width=10)

fig4 <- p_es_obs + p_afp_all
fig4

ggsave(here(figdir,"surveillance_activity.png"), fig4, height=5,width=10)

################################################################################
# Surveillance quality

afp_not_all %>%
  mutate(r_afp = r_afp*1e5) %>%
  group_by(guid) %>%
  mutate(cat_afp = case_when(all(r_afp >= 2) ~ "All",
                         sum(r_afp >= 2, na.rm = T) >= n()/2 ~ ">=50%",
                         sum(r_afp >= 2, na.rm = T) < n()/2 ~ "<50%",
                         !any(r_afp >= 2) ~ "None") %>% factor(levels = c("All",">=50%","<50%","None")),
         cat_adeq = case_when(all(r_adq >= 0.85) ~ "All",
                                  sum(r_adq >= 0.85, na.rm = T) >= n()/2 ~ ">=50%",
                                  sum(r_adq >= 0.85, na.rm = T) < n()/2 ~ "<50%",
                                  !any(r_adq >= 0.85) ~ "None") %>% factor(levels = c("All",">=50%","<50%","None"))) %>%
  ungroup() -> afp_not_wcat

afp_not_mth <- afp_not_all %>%
  group_by(tdy) %>%
  summarise(r_adq = mean(r_adq, na.rm = T)) %>%
  ungroup()

afp_not_wcat %>%
  group_by(guid) %>%
  summarise(cat_afp = unique(cat_afp),
            cat_adeq = unique(cat_adeq)) -> guid_afp_adeq

table(guid_afp_adeq$cat_afp, guid_afp_adeq$cat_adeq)
#       All >=50% <50% None
# All   630    68    1    0
# >=50%  63    10    0    0
# <50%    2     0    0    0
# None    0     0    0    0

afp_not_wcat %>%
  ggplot(aes(ymd(tdy), r_afp, col = cat_afp, group = interaction(guid, cat_afp))) +
  geom_line(lwd = 0.8) +
  geom_hline(aes(yintercept = 2), lty = "dashed") +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  scale_y_continuous(trans = "sqrt") +
  # scale_colour_manual(values = pal, breaks = c("high","mid","low"))+
  scale_colour_viridis_d(option = "turbo",begin = 0.1, end = 0.9) +
  theme(legend.position = "bottom") +
  labs(x = NULL,"Month of paralysis onset",
       y = "Rate per 100,000",
       col = "Months above\n2 per 100,000"
       # title = "AFP reporting rate, per 100,000 under-15s",
       # subtitle = "Rolling 12-month rate per district (by month of paralysis onset)"
       ) -> plot_afp
plot_afp
ggsave(here(figdir,"afp_not_rate.png"), plot_afp, height=5,width=10)

guid_afp_adeq %>%
  left_join(shape2,by = 'guid') %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill= cat_afp)) +
  scale_fill_viridis_d(name = NULL, # "Months above\n2 per 100,000",
                       option = "turbo", begin = 0.1, end = 0.9) +
  map_theme +
  theme(legend.position = "bottom") -> map_afp

comb_afp <- plot_afp + map_afp + plot_layout(widths = c(5,4))
comb_afp
ggsave(here(figdir,"plot_afp.png"), comb_afp, height=5,width=12)

afp_not_wcat %>%
  arrange(cat_adeq) %>%
  ggplot(aes(ymd(tdy), r_adq*100, col = cat_adeq, group = interaction(guid, cat_adeq))) +
  geom_line(lwd = 0.8, alpha = 0.7) +
  geom_hline(aes(yintercept = 85), lty = "dashed") +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  # scale_y_continuous(trans = "sqrt") +
  scale_colour_viridis_d(option = "turbo", begin = 0.1, end = 0.9) +
  theme(legend.position = "bottom") +
  labs(x = NULL, #"Month of paralysis onset",
       y = "% adequate",
       col = "Months above\n85% adequacy"
       # title = "Stool sample adequacy among reported AFP cases",
       # subtitle = "Rolling 12-month percentage per district (by month of paralysis onset)"
       ) -> plot_adeq
plot_adeq
ggsave(here(figdir,"afp_stool_adeq.png"), height=5,width=10)

# afp_not_wcat %>%
#   ggplot(aes(ymd(tdy), y = r_adq*100, ymin = r_adq_lwr*100, ymax = r_adq_upr*100, group = interaction(guid, cat_adeq), fill = cat_adeq, col = cat_adeq)) +
#   geom_ribbon(alpha = 0.005) +
#   geom_line(alpha = 0.5)+
#   # geom_smooth(aes(ymd(tdy), y = r_adq*100), se = T) +
#   geom_hline(aes(yintercept = 85), lty = "dashed") +
#   scale_x_date(labels = function(x) format(x, "%b-%Y")) +
#   scale_colour_viridis_d(option = "turbo", begin = 0.1, end = 0.9) +
#   scale_fill_viridis_d(option = "turbo", begin = 0.1, end = 0.9) +
#   theme(legend.position = "bottom") +
#   labs(x = NULL,
#        y = "% adequate",
#   ) -> plot_adeq2
# plot_adeq2

guid_afp_adeq %>%
  left_join(shape2,by = 'guid') %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill= cat_adeq)) +
  scale_fill_viridis_d(name = NULL, # "Months above\n85% adequacy",
                       option = "turbo", begin = 0.1, end = 0.9) +
  map_theme +
  theme(legend.position = "bottom") -> map_adeq

comb_adeq <- plot_adeq + map_adeq + plot_layout(widths = c(5,4))
comb_adeq
ggsave(here(figdir,"plot_adeq.png"), comb_adeq, height=5,width=12)

comb_afp / comb_adeq + plot_annotation(tag_levels = "A")
ggsave(here(figdir,"plot_afp_adeq.png"), height=10,width=12)

guid_afp_adeq %>%
  left_join(shape2,by = 'guid') %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill= (cat_afp == "<50%" & cat_adeq == "<50%"))) +
  scale_fill_viridis_d(NULL, option = "turbo", begin = 0.1, end = 0.9) +
  map_theme +
  ggtitle("Low AFP rate and stool adequacy")
ggsave(here(figdir,"map_low_afp_adeq.png"), height=5,width=7)

# ---------------------------------------------------------------------------- #
# ES #

data_all %>%
  group_by(tdy, tdy_f) %>%
  summarise(n_dist = n_distinct(guid),
            n_sites_tot = sum(n_sites),
            n_dist_wsites = sum(n_sites > 0),
            p_dist_wsites = n_dist_wsites/n_dist,
            es_covg_all = mean(es_coverage),
            es_covg_wsites = mean(es_coverage[n_sites > 0]),
            es_covg_wcovg = mean(es_coverage[es_coverage > 0])) %>%
  ungroup() -> dist_es_summ

dist_es_summ %>%
  ggplot(aes(ymd(tdy), p_dist_wsites*100)) +
  geom_col() +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  theme(plot.title = element_text(size = 14)) +
  labs(x = NULL, y = "Percentage",
       title = "% Districts with any active ES site(s)") -> p_any_site
p_any_site

dist_es_summ %>%
  pivot_longer(es_covg_all:es_covg_wcovg) %>%
  mutate(name = factor(name,
                       levels = c("es_covg_all","es_covg_wsites","es_covg_wcovg"),
                       labels = c("All districts",
                                  "Districts with any\nES site(s)",
                                  "Districts with any\nES coverage"))) %>%
  filter(name != "Districts with any\nES coverage") %>%
  ggplot(aes(ymd(tdy), value*100, lty = name)) +
  geom_line(lwd = 1) +
  # ylim(c(0,4)) +
  scale_y_continuous(trans = "log2") +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  theme(legend.position = c(0.8,0.5),
        plot.title = element_text(size = 14)) +
  labs(x = NULL, y = "Percentage", lty = NULL,
       title = "% Population within 5km of an ES site") -> p_any_covg
p_any_covg

comb_covg <- p_any_site + p_any_covg
comb_covg
ggsave(here(figdir,"es_dist_coverage.png"),comb_covg, height=5,width=10)

es_all %>%
  group_by(tdy, guid) %>%
  summarise(n = n(),
            n_npev = sum(npev == "Yes", na.rm = T),
            r_npev = n_npev/n) %>%
  group_by(guid) %>%
  mutate(cat_npev = case_when(all(r_npev >= 0.5) ~ "All",
                         sum(r_npev >= 0.5) >= n()/2 ~ ">=50%",
                         sum(r_npev >= 0.5) < n()/2 & any(r_npev >= 0.5) ~ "<50%",
                         !any(r_npev >= 0.5) ~ "None") %>% factor(levels = c("All",">=50%","<50%","None"))) %>%
  ungroup() -> es_wcat

es_wcat %>%
  group_by(guid) %>%
  summarise(cat_npev = unique(cat_npev)) -> guid_npev

guid_afp_adeq %>%
  full_join(guid_npev) -> guid_afp_es_qual

es_wcat %>%
  ggplot(aes(ymd(tdy), r_npev*100, group = interaction(guid, cat_npev))) +
  # geom_jitter(width = 2, alpha = 0.2) +
  geom_line(aes(colour = cat_npev), lwd = 1, alpha= 0.7) +
  # guides(col = "none") +
  geom_hline(yintercept = 50, lty = "dashed") +
  # geom_smooth() +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  scale_colour_viridis_d(option = "turbo", begin = 0.1, end = 0.9) +
  theme(legend.position = "bottom") +
  labs(x = NULL, #"Month of sample collection",
       y = "% NPEV detection",
       col = "Months above\n50% detection"
       # title = "NPEV detection rate among collected environmental samples",
       # subtitle = "Rolling 12-month percentage per district with any coverage (by month of sample collection))"
       ) -> plot_npev

plot_npev

ggsave(here(figdir,"es_npev_adeq.png"), height=5,width=10)

pal <- c(viridis::viridis(n = 4, option = "turbo", begin = 0.1, end = 0.9),"#FFFFFF")
guid_afp_es_qual %>%
  left_join(shape2,by = 'guid') %>%
  # mutate(cat_npev = replace_na(cat_npev,"No ES")) %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill= cat_npev)) +
  # scale_fill_manual(values = pal, "Months above\n50% detection") +
  scale_fill_viridis_d(NULL,#"Months above\n50% detection",
                       option = "turbo", begin = 0.1, end = 0.9,
                       na.value = "white", na.translate = F) +
  map_theme +
  theme(legend.position = "bottom") -> map_npev
map_npev

comb_npev <- plot_npev + map_npev + plot_layout(widths = c(5,4))
comb_npev
ggsave(here(figdir,"plot_npev.png"), comb_npev, height=5,width=12)

comb_covg / comb_npev + plot_annotation(tag_levels = "A", caption = "LGAs with no active ES are shown in white") #+ plot_layout(heights = c(1,2))
ggsave(here(figdir,"plot_es_adeq.png"), height=10,width=12)

# ES coverage
es_covg_mth <- data_all %>%
  group_by(tdy) %>%
  summarise(across(c(catchment_pop, es_coverage), mean, na.rm = T)) %>%
  ungroup()
es_covg_guid <- data_all %>%
  group_by(guid) %>%
  summarise(across(c(catchment_pop, es_coverage), mean, na.rm = T)) %>%
  ungroup()

sc <- max(es_covg_mth$catchment_pop)/max(es_covg_mth$es_coverage*100)
es_covg_mth %>%
  ggplot(aes(ymd(tdy))) +
  geom_line(aes(y = catchment_pop), col = "black") +
  geom_line(aes(y = es_coverage*100*sc), col = "steelblue", lty = "dashed") +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  scale_y_continuous(#limits = c(0,140),
                     name = "Catchment population count",
                     # Add a second axis and specify its features
                     sec.axis = sec_axis(trans = ~./sc,
                                         name = "LGA coverage (%)")) +
  labs(x = NULL,
       title = "Population catchment within 5km of active ES sites",
       subtitle = "Mean per LGA") -> p_es_catch
p_es_catch

ggsave(here(figdir,"es_pop_catch.png"), p_es_catch, height=5,width=10)

es_covg_guid %>%
  left_join(shape2,by = 'guid') %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill= es_coverage)) +
  scale_fill_viridis_c(trans = "sqrt",
                       "% ES coverage", option = "turbo", begin = 0.1, end = 0.9) +
  map_theme +
  theme(legend.position = "bottom") -> map_covg

map_covg + map_npev
ggsave(here(figdir,"map_covg_npev.png"), height=5,width=12)

data_all %>%
  group_by(tdy) %>%
  summarise(es_catch = mean(catchment_pop)) %>%
  ungroup() %>%
  ggplot(aes(ymd(tdy), es_catch, group = 1)) +
  # geom_jitter(alpha = 0.2) +
  # geom_line(aes(y = mean(es_coverage))) +
  geom_line() +
  scale_y_continuous() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  labs(x = NULL,y = "Count",title = "Mean population catchment of active ES sites") -> p_mean_catch
p_mean_catch

data_all %>%
  ggplot(aes(ymd(tdy), catchment_pop)) +
  geom_line(aes(group = guid), alpha = 0.2) +
  geom_smooth() +
  scale_y_continuous(trans = "sqrt") +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  labs(x = NULL,y = "Population count",
       title = "Population catchment within 5km of active ES sites",
       subtitle = "Per district") -> p_pop_catch
p_pop_catch

data_all %>%
  group_by(tdy) %>%
  summarise(es_coverage = mean(es_coverage)*100) %>%
  ungroup() %>%
  ggplot(aes(ymd(tdy), es_coverage, group = 1)) +
  # geom_jitter(alpha = 0.2) +
  # geom_line(aes(y = mean(es_coverage))) +
  geom_line() +
  scale_y_continuous() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  labs(x = NULL,y = "Percentage",title = "Mean district coverage") -> p_mean_covg
p_mean_covg

ggsave(here(figdir,"avg_es_covg.png"), p_mean_covg, height=5,width=6)

p_covg <- p_pop_catch + p_mean_covg + plot_layout(widths = c(2,1))
p_covg
ggsave(here(figdir,"es_pop_covg.png"), p_covg, height=5,width=12)

data_all %>%
  group_by(tdy) %>%
  summarise(n = n(),
            covd = sum(es_coverage > 0)) %>%
  ungroup() %>%
  ggplot(aes(ymd(tdy), covd)) +
  geom_col() +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  labs(x = NULL,y = "Count",
       title = "Number of districts with active ES") -> p_any_covg

data_all %>%
  group_by(tdy) %>%
  summarise(n = n(),
            covd = sum(es_coverage > 0),
            overall = mean(es_coverage)) %>%
  ungroup() %>%
  ggplot(aes(ymd(tdy), covd*100/n, group = 1)) +
  geom_line() +
  # scale_y_continuous(limits = c(0,70)) +
  scale_x_date(labels = function(x) format(x, "%b-%Y")) +
  labs(x = NULL,y = "Percentage",
       title = "Population coverage within districts with active ES") -> p_any_covg_perc

comb_es_covg <- p_any_covg + p_any_covg_perc
comb_es_covg

ggsave(here(figdir,"es_covg.png"), comb_es_covg, height=5,width=10)

# ES coverage relative to recent risk
data_all %>%
  ggplot(aes(risk_est_log, es_coverage)) +
  geom_jitter(alpha = 0.5) +
  # geom_smooth()+
  labs(x = "Log district risk estimate", y = "% ES coverage")
ggsave(here(figdir,"es_covg_vs_risk.png"), height=5,width=6)

################################################################################
# Sensitivity for each time point

## Main figure ##
m = 3
n = n_distinct(sens_all$tdy)
sens_all %>%
  pivot_longer(c("SSe_AFP","SSe_ENV")) %>%
  mutate(name = factor(name, labels = c("AFP","ENV"))) %>%
  ggplot(aes(tdy_f, value, colour = name)) +
  geom_boxplot() +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1,
                                   # Make every third month label transparent
                                   color=rep(c("black", rep("transparent", each = m-1)), n))) +
  scale_color_viridis_d(option = "turbo", begin = 0.2, end = 0.8, direction = -1) +
  theme(legend.position = "bottom") +
  labs(x = NULL,y = "Sensitivity",
       col = NULL
       # title = "Overall surveillance sensitivity",
       # subtitle = "By component"
       ) -> p_sens_bycomp
p_sens_bycomp
saveRDS(p_sens_bycomp, here(figdir,"p_sens_bycomp.rds"))

ggsave(here(figdir,"sens_bycomp.png"), p_sens_bycomp, height=5,width=10)

sens_all %>%
  ggplot(aes(tdy_f, SSe_ENV)) +
  geom_boxplot() +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1,
                                   # Make every third month label transparent
                                   color=rep(c("black", rep("transparent", each = m-1)), n))) +
  labs(x = NULL,y = "ENV Sensitivity",
       title = "Overall surveillance sensitivity",
       subtitle = "Environmental component only") -> p_sens_env
p_sens_env
ggsave(here(figdir,"sens_env.png"), p_sens_env, height=5,width=10)

sens_all %>%
  ggplot(aes(tdy_f, SSe_AFP)) +
  geom_boxplot() +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1,
                                   # Make every third month label transparent
                                   color=rep(c("black", rep("transparent", each = m-1)), n))) +
  labs(x = NULL,y = "AFP Sensitivity",
       title = "Overall surveillance sensitivity",
       subtitle = "AFP component only") -> p_sens_afp
p_sens_afp
ggsave(here(figdir,"sens_afp.png"), p_sens_afp, height=5,width=10)

sens_all %>%
  ggplot(aes(tdy_f, SSe_AE)) +
  geom_boxplot() +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1,
                                       # Make every third month label transparent
                                       color=rep(c("black", rep("transparent", each = m-1)), n))) +
  labs(x = NULL,y = "Sensitivity",
       title = "Overall surveillance sensitivity",
       subtitle = "AFP + ENV components") -> p_sens_ae
p_sens_ae
ggsave(here(figdir,"sens_ae.png"), p_sens_ae, height=5,width=10)

################################################################################
# Sensitivity within districts with ES

# Average catchment per LGA/per month
summary(data_all$es_coverage)
# 2014
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.00000 0.00000 0.01845 0.00000 0.99900

# 2016
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.00000 0.00000 0.03302 0.00000 0.99900

summary(data_all$es_coverage*data_all$total_pop)
# 2014
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0       0       0    7839       0  982321

# 2016
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0       0       0   12343       0  913787

sens_env_wcovg %>%
  group_by(guid, tdy) %>%
  summarise(across(es_coverage:SSe_ENV, mean)) %>%
  ungroup() %>%
  ggplot(aes(es_coverage, SSe_ENV)) +
  geom_jitter(alpha = 0.3) +
  scale_x_continuous(trans = "sqrt") +
  scale_y_continuous(trans = "sqrt") +
  labs(x = "ES coverage (% district population)",
       y = "Sensitivity")

ggsave(here(figdir,"SSeENV_vs_covg.png"), height=6,width=6)

perc_covg = 0.1
sens_env_grp <- sens_env_wcovg %>%
  filter(es_coverage > perc_covg) %>%
  group_by(guid, tdy) %>%
  summarise(SSe_ENV = mean(SSe_ENV)) %>%
  ungroup() %>%
  mutate(cat = paste0(">",perc_covg*100,"% coverage")) %>%
  bind_rows(sens_env_wcovg %>% group_by(guid, tdy) %>% summarise(SSe_ENV = mean(SSe_ENV)) %>% ungroup() %>% mutate(cat = "All"))

# Mean LGA/month sensitivity across (1) all LGAs/mths and (2) LGA/mths with > 10% ES coverage
sens_env_grp %>%
  group_by(cat) %>%
  summarise(mean_sens = mean(SSe_ENV),
            iqr_sens = paste(round(quantile(SSe_ENV, c(0.025,0.975)),4), collapse = ","))
# 2014-16
# 1 >10% coverage    0.959  0.6529,1
# 2 All              0.0255 0,0.3941

# 2016-20
# 1 >10% coverage    0.913  0.4287,1
# 2 All              0.0519 0,0.9859

sens_env_wcovg %>%
  mutate(group = es_coverage > 0.1) %>%
  ggplot(aes(x = group, y = SSe_ENV)) +
  geom_violin() +
  labs(y = "ES sensitivity", x = NULL) +
  theme(legend.position = c(0.7,0.7))

# sens_bycovg %>%
#   ggplot(aes(cat, SSe_ENV)) +
#   geom_violin() +
#   # geom_jitter(alpha = 0.01) +
#   # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   scale_y_continuous(trans = "sqrt") +
#   labs(x = NULL,y = "ENV Sensitivity",
#        col = NULL,
#        title = "Overall surveillance sensitivity",
#        subtitle = "Environmental component only") -> p_sens_env
# p_sens_env
# ggsave(here(figdir,"SSeENV_by_covg_sqrt.png"), p_sens_env, height=5,width=6)

hist(sens_env_wcovg$SSe_ENV, prob = T, xlab = "ES sensitivity", main = NULL)
summary(sens_env_wcovg$SSe_ENV)

################################################################################
# Relative sensitivity of AFP/ES

sens_all %>%
  mutate(ratio = SSe_ENV/SSe_AFP) %>%
  ggplot(aes(tdy_f, ratio)) +
  geom_hline(yintercept = 1, lty = "dashed", col = "grey40") +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_viridis_d(option = "turbo", begin = 0.2, end = 0.8, direction = -1) +
  labs(x = NULL,y = "Ratio",
       col = NULL,
       title = "Relative sensitivity of ENV versus AFP surveillance") -> p_sens_rel

p_sens_rel
ggsave(here(figdir,"sens_rel.png"), p_sens_rel, height=5,width=10)

################################################################################
# Total system sensitivity

# sens_all %>%
#   mutate(total = SSe_ENV+SSe_AFP) %>%
#   ggplot(aes(tdy_f, total)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,
#                                    # Make every third month label transparent
#                                    color=rep(c("black", rep("transparent", each = m-1)), n))) +
#   # scale_color_viridis_d(option = "turbo", begin = 0.2, end = 0.8, direction = -1) +
#   labs(x = NULL,y = "Sensitivity",
#        col = NULL,
#        title = "Total surveillance system sensitivity (ENV + AFP)") -> p_sens_tot
#
# p_sens_tot
# ggsave(here(figdir,"sens_tot.png"), p_sens_tot, height=5,width=10)

################################################################################
################################################################################
