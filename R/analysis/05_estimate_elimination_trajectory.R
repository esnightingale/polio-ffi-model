################################################################################
################################################################################
# Script to estimate the accumulating freedom from infection probability over a
# given period of zero detection, according to calculated system sensitivity.
#
# FFI calculated with respect to two assumptions about sensitivity:
# + Static, i.e. assuming sensitivity calculated as of Aug 2014 continues
# + Time-varying, i.e. updating each month's FFI probability with respect to the
#   surveillace sensitivity as of the previous month
################################################################################

# Set region and time period of interest
source("R/utils/setup_env.R")

cc_name <- tolower(paste0(cc$name, collapse = "_"))

# Local projection (for drawing buffers & estimating ES catchment)
proj_local = cc$proj

# Virus type to analyse
type = cc$type

# Specify genetic cluster to analyse
clust = cc$cluster

# Specify emergence group to analyse
eg = cc$eg

# Results directory
if(!is.null(clust)){outdir <- here("results",cc_name,type,clust,substr(cc$start,1,4),)
}else if (!is.null(eg)){outdir <- here("results",cc_name,type,eg,substr(cc$start,1,4))
}else {outdir <- here("results",cc_name,type,substr(cc$start,1,4))}

# Figure and table directories
figdir <- here(outdir,"figures")
tabdir <- here(outdir,"tables")
if(!dir.exists(tabdir)) dir.create(tabdir, recursive = T)

################################################################################

# Read sensitivity estimates for each time point
sens_all <- readRDS(here(outdir,"sens_all.rds")) %>%
  # subset to period after last detection (if sens calculated for longer period)
  filter(between(tdy, cc$start, cc$end))

# Defining a critical confidence value to illustrate on the plot
elim_criterion = 0.95

# Number of years ahead for projection
future = F
M = 5

################################################################################

source(here("R/utils/estimate_elimination_fcn.R"))

################################################################################
# Run FFI estimation

## Increase introduction risk for cVDPV2? By how much?

# Single run
prior = c(30,30)
out_ffi_tv <- elim_est_tv(sens_all,
                       intro = T,
                       future = future, M = M,
                       prior_parms = prior)

# Compare time-varying vs static
compare_tv <- lapply(c(T,F),
                     function(tv){elim_est_tv(sens_all,
                                              future = future, M = M,
                                              tv = tv)})

# Compare priors
prior_opts <- list(c(10,30),
                   # c(5,10),
                   c(30,30),
                   # c(10,5),
                   c(30,10))
compare_priors <- lapply(prior_opts,
                         function(prior){elim_est_tv(sens_all,
                                                     future = future, M = M,
                                                     prior_parms = prior)})

# Compare intro risk
intro_opts <- list("Low (1/1000)" = c(lower = 1/5000, best = 1/1000, upper = 1/500, p = 0.95),
                   "Mid (1/500)" = c(lower = 1/1000, best = 1/500, upper = 1/100, p = 0.95),
                   "High (1/100)" = c(lower = 1/1000, best = 1/100, upper = 1/10, p = 0.95))
compare_intro <- lapply(intro_opts,
                         function(intro){elim_est_tv(sens_all,
                                                     future = future, M = M,
                                                     prior_parms = prior,
                                                     intro_parms = intro)})

# ---------------------------------------------------------------------------- #
# Summary tables

# Summarise probability of FFI
bind_rows(
  mutate(out_ffi_tv$free_AFP_sims,grp = "AFP only"),
  mutate(out_ffi_tv$free_ENV_sims,grp = "ENV only"),
  mutate(out_ffi_tv$free_AE_sims, grp = "Both components")) %>%
  mutate(t = period_to_date(time),
         grp = factor(grp, levels = c("AFP only","ENV only","Both components"))) %>%
  select(grp, t, time, projection, future_flag, everything()) -> ffi_compare_surv

if(future){
  write.csv(ffi_compare_surv,
            here(tabdir,paste0(paste0("ffi_mth_tab_",year(cc[3]),"_ahead",M,".csv"))))

}else{
  write.csv(ffi_compare_surv,
            here(tabdir,paste0(paste0("ffi_mth_tab_",year(cc[3]),".csv"))))

}

# Estimated sensitivity per month:
ae_tmp <- out_ffi_tv$SSe_AE %>%
  as_tibble() %>%
  mutate(time = as.vector(out_ffi_tv$free_AE_sims$time)) %>%
  pivot_longer(-time) %>%
  group_by(time) %>%
  summarise(med = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>%
  mutate(grp = "Both components")
afp_tmp <- out_ffi_tv$SSe_AFP %>%
  as_tibble() %>%
  mutate(time = as.vector(out_ffi_tv$free_AFP_sims$time)) %>%
  pivot_longer(-time) %>%
  group_by(time) %>%
  summarise(med = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>%
  mutate(grp = "AFP only")
env_tmp <- out_ffi_tv$SSe_ENV %>%
  as_tibble() %>%
  mutate(time = as.vector(out_ffi_tv$free_ENV_sims$time)) %>%
  pivot_longer(-time) %>%
  group_by(time) %>%
  summarise(med = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>%
  mutate(grp = "ENV only")

ae_tmp %>%
  bind_rows(afp_tmp) %>%
  bind_rows(env_tmp) %>%
  mutate(t = period_to_date(time),
         grp = factor(grp, levels = c("AFP only","ENV only","Both components"))) %>%
  select(grp, t, time, everything()) -> sens_mth_tab

write.csv(sens_mth_tab,
          here(tabdir,paste0("sens_mth_tab_",year(cc[3]),".csv")))

# Average sensitivity overall
bind_rows(summary(as.vector(out_ffi_tv$SSe_AE)),
          summary(as.vector(out_ffi_tv$SSe_AFP)),
          summary(as.vector(out_ffi_tv$SSe_ENV))) %>%
  janitor::clean_names() %>%
  mutate(type = c("AE","AFP","ES"),
         across(-type, as.numeric)) %>%
  select(type, everything()) -> sens_tab
sens_tab

write.csv(sens_tab,
          here(tabdir,paste0("avg_sens_tab_",year(cc[3]),".csv")))

# ---------------------------------------------------------------------------- #
# Plots

# Setup
pd = position_dodge2(width = 0.05)
pal = viridis::turbo(3, begin = 0.2, end = 0.8)

# Highlight projected period if applicable
future <- min(ffi_compare_surv$time[ffi_compare_surv$future])

# Comparator time points for 2016-2020
detection = date_to_period("2016-07-01")
dietz = 2.9
declaration = date_to_period(ymd("2020-08-24"))

# ---------------------------------------------------------------------------- #
# Compare components

t1 = min(ffi_compare_surv$time)
tmax = max(ffi_compare_surv$time)

ffi_compare_surv %>%
  ggplot(aes(x=time,y=med)) +
  # Formatting
  geom_hline(aes(yintercept=0.95),linetype="dashed",col="grey50") +
  geom_hline(aes(yintercept=0.99),linetype="dashed",col="red") +
  annotate("text", x = t1, y = 0.93, label = "0.95", col="grey50", size = 2) +
  annotate("text", x = t1, y = 0.97, label = "0.99", col="red", size = 2) +
  # geom_vline(aes(xintercept = t1 + dietz)) +
  # annotate("text", x = eval(t1 + dietz - 0.45), y = 0.2,
  #          label = "Eichner & Dietz (1992) estimate\nfor risk of silent infections < 5%",
  #          size = 3) +
  # geom_vline(aes(xintercept = declaration)) +
  # annotate("text", x = eval(declaration - 0.35), y = 0.3,
  #          label = "Official declaration of\nWPV elimination in Nigeria",
  #          size = 3) +
  # geom_vline(aes(xintercept = detection), lty= "solid") +
  # annotate("text", x = eval(detection - 0.15), y = 0.2,
  #          label = "Detection of\nWPV1 in Borno state",
  #          size = 3) +
  # Monthly FFI probabilities
  geom_pointrange(aes(ymin=lwr, ymax=upr,
                      col = grp,
                      alpha = grp),
                  position = pd) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c(pal[c(3,1)],"black")) +
  scale_alpha_manual(values = c(0.7,0.7,1)) +
  scale_x_continuous(breaks = unique(round(ffi_compare_surv$time))) +
  guides(alpha = "none") +
  labs(x = "Month",
       y = "Pr(Infection free)",
       col = NULL,
       fill = NULL,
       title = NULL) + #"Freedom from infection (with introduction risk)") +
  ylim(0,1) -> p_ffi_comp_surv

if (!is.infinite(future)){
  p_ffi_comp_surv <- p_ffi_comp_surv +
  geom_rect(aes(xmin = future - 0.03,
                xmax = round(ifelse(is.infinite(future), Inf, tmax)),
                ymin = 0,
                ymax = 1),
            fill = "grey80",
            alpha = 0.01) +
    annotate("text",x = future + 0.2,y = 0.1, label = "Projection",col="grey", size = 4)
}

p_ffi_comp_surv
ggsave(here(figdir,
            paste0("compare_ffi_survtype_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),".png")),
       p_ffi_comp_surv, height=6,width=10)

# Combine with sensitivity estimates for MS
p_sens_bycomp <- readRDS(here(figdir,"p_sens_bycomp.rds")) #+ geom_vline(xintercept = "Jul-2016", lty = "dotted")
p_comb <- p_sens_bycomp / p_ffi_comp_surv +
  plot_layout(heights = c(2,3)) +
  plot_annotation(tag_levels = "A")
p_comb
ggsave(here(figdir,
            paste0("comb_compare_ffi_survtype_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),".png")),
       p_comb, height=9,width=10)

# ---------------------------------------------------------------------------- #
# Total system trajectory

ffi_compare_surv %>%
  filter(grp == "Both components") %>%
  ggplot(aes(x=time,y=med)) +
  # Formatting
  geom_hline(aes(yintercept=0.95),linetype="dashed",col="grey50") +
  geom_hline(aes(yintercept=0.99),linetype="dashed",col="red") +
  annotate("text", x = t1, y = 0.93, label = "0.95", col="grey50", size = 2) +
  annotate("text", x = t1, y = 0.97, label = "0.99", col="red", size = 2) +
  geom_rect(aes(xmin = future - 0.03,
                xmax = round(ifelse(is.infinite(future), Inf, tmax)),
                ymin = 0,
                ymax = 1),
            fill = "grey80",
            alpha = 0.01) +
  annotate("text",x = future + 0.5,y = 0.1, label = "Projection",col="grey", size = 4) +
  geom_vline(aes(xintercept = t1 + dietz)) +
  annotate("text", x = eval(t1 + dietz - 0.45), y = 0.2,
           label = "Eichner & Dietz (1992) estimate\nfor risk of silent infections < 5%",
           size = 3) +
# Monthly FFI probabilities
geom_pointrange(aes(ymin=lwr, ymax=upr),
                position = pd) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(round(ffi_compare_surv$time))) +
  guides(alpha = "none") +
  labs(x = "Month",
       y = "Pr(Infection free)",
       col = NULL,
       fill = NULL,
       title = NULL) +
  ylim(0,1) -> p_ffi_total

p_ffi_total

if(future){
  ggsave(here(figdir,
              paste0("ffi_total_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),"_ahead",M,".png")),
         p_ffi_total, height=6,width=10)
}else{
  ggsave(here(figdir,
              paste0("ffi_total_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),".png")),
         p_ffi_total, height=6,width=10)
}

# ---------------------------------------------------------------------------- #
# Probabilities at specific time points

# P[FFI] @ Jul-2016 when WPV1 re-emerged:
ffi_compare_surv %>%
  group_by(grp) %>%
  filter(t == period_to_date(detection)) %>%
  select(lwr, med, upr)
# grp               lwr   med   upr
# 1 AFP only        0.641 0.750 0.832
# 2 ENV only        0.534 0.652 0.755
# 3 Both components 0.772 0.853 0.902

# P[FFI] @ Aug-2020 when official declaration made:
ffi_compare_surv %>%
  group_by(grp) %>%
  filter(t == period_to_date(declaration)) %>%
  select(lwr, med, upr)
# 1 AFP only        0.833 0.886 0.923
# 2 ENV only        0.868 0.910 0.939
# 3 Both components 0.976 0.982 0.986

# ---------------------------------------------------------------------------- #
# Compare static vs time-varying sensitivity

lapply(compare_tv,
       function(x) x$free_AE_sims) %>%
  setNames(c("Time-varying","Static")) %>%
  bind_rows(.id = "sensitivity") %>%
  mutate(t = period_to_date(time)) -> ffi_compare_tv

# Plot
plot_ffi_compare <- function(ffi_compare,
                             groupvar,
                             grouplab,
                             pd = position_dodge2(width = 0.05)){
  ffi_compare %>%
    ggplot(aes(x=time,y=med)) +
    # Monthly FFI probabilities
    geom_pointrange(aes(ymin=lwr, ymax=upr,
                        col = !!sym(groupvar)),
                    position = pd) +
    # Formatting
    geom_hline(aes(yintercept=0.95),linetype="dashed",col="grey50") +
    geom_hline(aes(yintercept=0.99),linetype="dashed",col="red") +
    geom_text(aes(x = min(time),y = 0.93), label = "0.95",col="grey50", cex = 2) +
    geom_text(aes(x = min(time),y = 0.97), label = "0.99",col="red", cex = 2) +
    scale_colour_manual(values = c("grey70","black")) +
    scale_x_continuous(breaks = unique(round(ffi_compare$time))) +
    labs(x = "Month",
         y = "Pr(Infection free)",
         col = grouplab,
         fill = grouplab,
         title = NULL) + #"Freedom from infection (with introduction risk)") +
    theme(legend.position = "bottom") +
    ylim(0,1) -> p

  return(p)

}

p_ffi_comp_tv <- plot_ffi_compare(ffi_compare_tv, "sensitivity", "Sensitivity")
p_ffi_comp_tv

ggsave(here(figdir,
            paste0("compare_ffi_tv_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),".png")),
       p_ffi_comp_tv, height=6,width=10)

# Combine with plot of total sensitivity over time
m = 3
n = n_distinct(sens_all$tdy)
sens_static <- mutate(sens_all, grp = "static")
sens_static$SSe_AE <- sens_static$SSe_AE[sens_static$tdy == min(sens_static$tdy)]
sens_all %>%
  mutate(grp = "tv") %>%
  bind_rows(sens_static) %>%
  ggplot(aes(tdy_f, SSe_AE, col = grp)) +
  geom_boxplot() +
  scale_colour_manual(values = c("grey70","black")) +
  guides(col = "none") +
  theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1,
                                   # Make every third month label transparent
                                   color=rep(c("black", rep("transparent", each = m-1)), n))) +
  labs(x = NULL,y = "Total system sensitivity") -> p_sens
p_sens

p_comb <- p_sens / p_ffi_comp_tv +  plot_layout(heights = c(2,3)) + plot_annotation(tag_levels = "A")
p_comb
ggsave(here(figdir,
            paste0("comb_compare_ffi_tv_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),".png")),
       p_comb, height=10,width=10)

p_ffi_comp_tv$data %>%
  group_by(sensitivity) %>%
  filter(med > 0.95) %>%
  slice_min(time) -> ffi_comp_tv_tab
ffi_comp_tv_tab
#   sensitivity    lwr   med   upr  time t
# 1 Static       0.926 0.951 0.965 2020. 2019-Nov
# 2 Time-varying 0.929 0.953 0.968 2019. 2019-Jun

write.csv(ffi_comp_tv_tab,
          here(tabdir,paste0("ffi_comp_tv_tab_",year(cc[3]),".csv")))

# ---------------------------------------------------------------------------- #
# Compare between prior assumptions

lapply(compare_priors,
         function(x) x$free_AE_sims) %>%
    setNames(prior_opts) %>%
    bind_rows(.id = "prior") %>%
    mutate(t = period_to_date(time),
           prior = factor(prior, levels = c("c(10, 30)",
                                            "c(30, 30)",
                                            "c(30, 10)"))) -> ffi_compare_priors

  saveRDS(ffi_compare_priors, here::here(outdir, "ffi_compare_priors.rds"))

  # Identify time point (if any) at which confidence criterion for elimination is reached
  # time_AE_FFI95 <- min(ffi_all$time[ffi_all$med >= elim_criterion])

  # Plot
  ffi_compare_priors %>%
    ggplot(aes(x=time,y=med)) +
    # Monthly FFI probabilities
    geom_pointrange(aes(ymin=lwr, ymax=upr,
                        col = prior),
                    position = pd) +
    # Formatting
    geom_hline(aes(yintercept=0.95),linetype="dashed",col="grey50") +
    geom_hline(aes(yintercept=0.99),linetype="dashed",col="red") +
    geom_text(aes(x = min(time),y = 0.93), label = "0.95",col="grey50", cex = 2) +
    geom_text(aes(x = min(time),y = 0.97), label = "0.99",col="red", cex = 2) +
    # geom_vline(aes(xintercept = detection), lty= "solid") +
    # geom_text(aes(x = detection - 0.15, y = 0.2), label = "Detection of\nWPV1 in Borno state", size = 3) +
    geom_vline(aes(xintercept = min(time) + dietz)) +
    geom_text(aes(x = min(time) + dietz - 0.45, y = 0.2), label = "Eichner & Dietz (1992) estimate\nfor risk of silent infections < 95%", size = 3) +
    # geom_vline(aes(xintercept = declaration)) +
    # geom_text(aes(x = declaration - 0.35, y = 0.3), label = "Official declaration of\nWPV elimination in Nigeria", size = 3) +
    scale_colour_viridis_d(option = "mako", begin = 0.2, end = 0.8) +
    # scale_y_continuous(trans = "log") +
    scale_x_continuous(breaks = unique(round(ffi_compare_priors$time))) +
    theme(legend.position = "bottom") +
    ylim(0,1) +
    labs(x = "Month",
         y = "Pr(Infection free)",
         col = "Beta prior",
         fill = "Beta prior",
         title = NULL) -> p_ffi_comp_prior

  p_ffi_comp_prior

  ggsave(here(figdir,
              paste0("compare_ffi_priors_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),".png")),
         p_ffi_comp_prior, height=6,width=10)

# ---------------------------------------------------------------------------- #
# Compare between intro risk assumptions

lapply(compare_intro,
         function(x) x$free_AE_sims) %>%
    setNames(names(intro_opts)) %>%
    bind_rows(.id = "intro") %>%
    mutate(t = period_to_date(time),
           intro = factor(intro,
                          levels = names(intro_opts))) -> ffi_compare_intro

  saveRDS(ffi_compare_intro, here::here(outdir, "ffi_compare_intro.rds"))

  # Plot
  ffi_compare_intro %>%
    ggplot(aes(x=time,y=med)) +
    # Monthly FFI probabilities
    geom_pointrange(aes(ymin=lwr, ymax=upr,
                        col = intro),
                    position = pd) +
    # Formatting
    geom_hline(aes(yintercept=0.95),linetype="dashed",col="grey50") +
    geom_hline(aes(yintercept=0.99),linetype="dashed",col="red") +
    geom_text(aes(x = min(time),y = 0.93), label = "0.95",col="grey50", cex = 2) +
    geom_text(aes(x = min(time),y = 0.97), label = "0.99",col="red", cex = 2) +
    # geom_vline(aes(xintercept = detection), lty= "solid") +
    # geom_text(aes(x = detection - 0.15, y = 0.2), label = "Detection of\nWPV1 in Borno state", size = 3) +
    # geom_vline(aes(xintercept = min(time) + dietz)) +
    # geom_text(aes(x = min(time) + dietz - 0.45, y = 0.2), label = "Eichner & Dietz (1992) estimate\nfor risk of silent infections < 95%", size = 3) +
    # geom_vline(aes(xintercept = declaration)) +
    # geom_text(aes(x = declaration - 0.35, y = 0.3), label = "Official declaration of\nWPV elimination in Nigeria", size = 3) +
    scale_colour_viridis_d(option = "mako", begin = 0.2, end = 0.8) +
    # scale_y_continuous(trans = "log") +
    scale_x_continuous(breaks = unique(round(ffi_compare_intro$time))) +
    theme(legend.position = "bottom") +
    ylim(0,1) +
    labs(x = "Month",
         y = "Pr(Infection free)",
         col = "Monthly introduction risk",
         fill = "Introduction risk",
         title = NULL) -> p_ffi_comp_intro

  p_ffi_comp_intro

  ggsave(here(figdir,
              paste0("compare_ffi_intro_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),".png")),
         p_ffi_comp_intro, height=6,width=10)

################################################################################
################################################################################
