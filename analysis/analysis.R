#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lme4)
library(car)
library(AICcmodavg)
library(lemon)
library(lubridate)
library(cowplot)
library(ggsn)
source("analysis/bootstrap.R")
source("analysis/map_prep.R")


# import data
# z-score continuous variables
hobo <- read_csv("data/clean/hobo.csv") 
met_long <- read_csv("data/clean/met_long.csv") %>%
  mutate(midge_z = (midge_treat - mean(midge_treat, na.rm=T))/(2*sd(midge_treat, na.rm=T)),
         temp_z = (temperature - mean(temperature, na.rm=T))/(2*sd(temperature, na.rm=T)),
         par_z = (par - mean(par, na.rm=T))/(2*sd(par, na.rm=T)),
         date = as.factor(sampledate))
# note, metabolism data are in g m^-2 h^-1 (get converted below)
adults <- read_csv("data/clean/adults.csv") %>%
  mutate(midge_z = (midge_treat - mean(midge_treat, na.rm=T))/(2*sd(midge_treat, na.rm=T)))

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(1,1,1,1),
                  legend.margin = margin(0,0,0,-4),
                  legend.text = element_text(size = 8),
                  axis.text = element_text(size = 10, color="black", family = "sans"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,10,0,0)),
                  axis.title.x = element_text(margin = margin(10,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))

# define site colors
site_colors <- c("dodgerblue","firebrick","gray35")
names(site_colors) <- c("E2","E3","E5")




#==========
#========== GPP
#==========

# gpp_dat
gpp_dat <- met_long %>% filter(type == "gpp")

# fit model
gpp_m <- lmer(met ~ temp_z + par_z + (trial + midge_z + site)^2 
                + I(midge_z^2) + I(midge_z^3) + (1|rack/core), 
              data = gpp_dat)

# ANOVA
gpp_anova_iii <- Anova(gpp_m, type = 3, test.statistic = "F") %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "rowname") %>%
  set_names(c("term" ,"f","df_n","df_d","p")) %>%
  mutate(f  = signif(f, 2),
         df_n = signif(df_n, 3),
         df_d = signif(df_d, 3),
         p = ifelse(p < 0.001, 0, round(p, 3))) %>%
  filter(term != "(Intercept)")
gpp_anova_ii <- Anova(gpp_m, type = 2, test.statistic = "F") %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "rowname") %>%
  set_names(c("term" ,"f","df_n","df_d","p")) %>%
  mutate(f  = signif(f, 2),
         df_n = signif(df_n, 3),
         df_d = signif(df_d, 3),
         p = ifelse(p < 0.001, 0, round(p, 3))) 

# clean
gpp_anova <- tibble(Term = factor(gpp_anova_iii$term,
                                  levels = c("temp_z",
                                             "par_z",
                                             "trial",
                                             "midge_z",
                                             "site",
                                             "I(midge_z^2)", 
                                             "I(midge_z^3)",
                                             "trial:midge_z",
                                             "trial:site",
                                             "midge_z:site"),
                                  labels = c("temp_z",
                                             "par_z",
                                             "trial",
                                             "midge_z",
                                             "site",
                                             "I(midge_z^2)", 
                                             "I(midge_z^3)",
                                             "trial x midge",
                                             "trial x site",
                                             "midge x site")),
                    `F_iii` = paste0(gpp_anova_iii$f, "_{",
                                     gpp_anova_iii$df_n, ",",
                                     gpp_anova_iii$df_d,"}"),
                    `P_iii` = gpp_anova_iii$p,
                    `F_ii` = c(paste0(gpp_anova_ii$f, "_{",
                                      gpp_anova_ii$df_n, ",",
                                      gpp_anova_ii$df_d,"}")[1:5],
                               "","","","",""),
                    `P_ii` = c(gpp_anova_ii$p[1:5],
                               "","","","",""))
gpp_anova
gpp_anova %>% knitr::kable(format = "latex")

# random effects
summary(gpp_m)$varcor

# coefficients
coefs <- summary(gpp_m)$coefficients %>% round(3)
coefs

# standardize observations
gpp_dat = gpp_dat %>%
  mutate(temp_cor = coefs[which(rownames(coefs)=="temp_z"), 1]*temp_z,
         met_s = met - temp_cor + mean(temp_cor))

# generate predicted data
gpp_nd <- gpp_dat %>% 
  tidyr::expand(trial, site, midge_z, 
                temp_z = 0, par_z = 0, rack = 1, coried = 1) %>%
  left_join(gpp_dat %>% select(midge_z, midge_density) %>% unique)
gpp_preds <- predictSE.merMod(gpp_m, 
                              newdata = gpp_nd, 
                              REForm = NA, 
                              print.matrix = T)
gpp_nd$pred <- gpp_preds[,1]
gpp_nd$se <- gpp_preds[,2]


# plot
p_gpp <- gpp_nd %>%
  mutate(trial = factor(trial, levels = c(0,1), labels = c("Day 5 or 6", 
                                                           "Day 16 or 18"))) %>%
  ggplot(aes(midge_density/1000, pred, color = site))+
  facet_wrap(~trial, nrow = 1)+
  geom_line(aes(y = pred),
            size = 0.4)+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se, fill = site), 
                linetype = 0,
                alpha = 0.1)+
  geom_jitter(data = gpp_dat%>%
               mutate(trial = factor(trial, levels = c(0,1), labels = c("Day 5 or 6", 
                                                                        "Day 16 or 18"))),
             aes(y = met_s, fill = site),
             shape = 21,
             size = 1.5,
             stroke = 0,
             color = "black",
             height = 0,
             width = 0)+
  scale_color_manual("",values=site_colors)+
  scale_fill_manual(values=site_colors, guide = F)+
  scale_y_continuous(expression("GPP (mg "*O[2]~m^{-2}~h^{-1}*")"),
                     limits = c(0.0475, 0.205),
                     breaks = c(0.05, 0.10, 0.15, 0.20),
                     labels = 1000 * c(0.05, 0.10, 0.15, 0.20))+
  scale_x_continuous(expression("Larval density ("*1000~m^{-2}*")"),
                     breaks = c(0, 30, 60, 90),
                     limits = c(-5, 104))+
  theme(legend.position = c(0.2, 0.9),
        legend.direction = "horizontal",
        legend.text = element_text(margin = margin(l = -6)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.x = unit(0.4, "lines"),
        panel.border = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  guides(color = guide_legend(override.aes = list(size = 0.5, fill = NA)))
p_gpp





#==========
#========== Adults
#==========

# select final data of adults
adults_0_a <- adults  %>%
  mutate(duration = as.numeric(max(sampledate) - min(sampledate)))%>%
  filter(sampledate == max(sampledate))  %>%
  mutate(Y1 = cumem,
         Y2 = midge_treat - cumem)
cor.test(adults_0_a$midge_treat, 
         adults_0_a$cumem,
         method = "spearman")

adults_0 <- adults_0_a %>%
  filter(midge_treat >= 50)

adults_m <- glmer(cbind(Y1, Y2) ~ midge_z*site + (1|rack/core),
                  data = adults_0, family = "binomial")
adults_m_ii <- glmer(cbind(Y1, Y2) ~ midge_z + site + (1|rack/core),
                  data = adults_0, family = "binomial")

# p-values
Anova(adults_m, type = 3)
Anova(update(adults_m, .~.-midge_z:site))

# adults_m_boot <- bootLRT(adults_m, nboot = 2000, c = 2) 
# adults_m_boot_ii <- bootLRT(adults_m_ii, nboot = 2000, c = 2) 

# adults_m_boot %>%
#   write_csv("analysis/boot/adults_m_boot.csv")
# adults_m_boot_ii %>%
#   write_csv("analysis/boot/adults_m_boot_ii.csv")

adults_m_boot_in <- read_csv("analysis/boot/adults_m_boot.csv")
adults_m_boot_ii_in <- read_csv("analysis/boot/adults_m_boot_ii.csv")

adults_m_anova <- tibble(Term = factor(adults_m_boot_in$term,
                     levels = c("midge_z",
                                "site",
                                "midge_z:site"),
                     labels = c("density",
                                "site",
                                "density x site")),
       `X_iii` = round(adults_m_boot_in$chi, 1),
       `P_iii` = round(adults_m_boot_in$p, 3),
       `X_ii` = c(round(adults_m_boot_ii_in$chi, 1), ""),
       `P_ii` = c(round(adults_m_boot_ii_in$p, 3), ""))

adults_m_anova %>% knitr::kable(format = "latex")

# random effects
summary(adults_m)$varcor

# coefficients
coefs <- summary(adults_m)$coefficients %>% round(3)
coefs

# generate predicted values
adults_nd <- adults_0 %>%
  tidyr::expand(site, midge_z = seq(min(midge_z)-0.4, max(midge_z), length.out = 100), 
                rack = 1, coried = 1) %>%
  mutate(midge_treat = 2*midge_z*sd(adults_0$midge_treat) + mean(adults_0$midge_treat),
         midge_density = midge_treat*unique(adults_0$midge_density/adults_0$midge_treat))
adults_preds<- predictSE.merMod(adults_m, newdata = adults_nd, REForm = NA, print.matrix = T)
adults_nd$pred <- adults_preds[,1]
adults_nd$se <- adults_preds[,2]

# plot
adults_nd %>%
  ggplot(aes(midge_density/1000, pred, color = site))+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se, fill = site),
              alpha = 0.15, linetype = 0)+
  geom_line(aes(y = pred), size = 0.8)+
  geom_jitter(data = adults_0, aes(y = prop), size = 2, alpha = 0.7, width = 2)+
  scale_color_manual(values=site_colors)+
  scale_fill_manual(values=site_colors)+
  scale_y_continuous("Proportion Emerged")+
  scale_x_continuous(expression("Midge Density ("*1000~m^{-2}*")"))






#==========
#========== Adults ~ GPP
#==========

# calculate average GPP and average predicted GPP through time
adults_new <- adults_0 %>%
  full_join(gpp_dat %>%
              filter(midge_treat > 0)%>%
              group_by(core) %>%
              summarize(met_s = mean(met_s),
                        midge_density = unique(midge_density)) %>%
              mutate(met_per = 1000 * met_s / midge_density, 
                     # calculate production per midge
                     # convert production to mg O2 (x 1000)
                     met_per_z = (met_per - mean(met_per))/sd(met_per)))

# fit_model
adults_gpp <- glmer(cbind(Y1, Y2) ~ met_per_z + (1|rack/core),
                  data = adults_new, family = "binomial")

# p-values
Anova(adults_gpp)
# bootLRT(adults_gpp, nboot = 2000, c = 2) %>% 
#   write_csv("analysis/boot/adults_gpp_boot.csv")

# random effects
summary(adults_gpp)$varcor

# coefficients
coefs <- summary(adults_gpp)$coefficients %>% round(3)
coefs

# generate predicted values
adults_gpp_nd <- adults_new %>%
  tidyr::expand(nesting(met_per, met_per_z))
adults_gpp_preds<- predictSE.merMod(adults_gpp, 
                                    newdata = adults_gpp_nd, 
                                    REForm = NA, 
                                    print.matrix = T)
adults_gpp_nd$pred <- adults_gpp_preds[,1]
adults_gpp_nd$se <- adults_gpp_preds[,2]

# plot
p_adults <- adults_gpp_nd %>%
  ggplot(aes(met_per, pred))+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se),
              alpha = 0.1, linetype = 0)+
  geom_line(aes(y = pred), size = 0.5)+
  geom_point(data = adults_new, aes(y = prop, fill = site), 
             shape = 21,
             size = 2,
             stroke = 0,
             color = "black")+
  scale_fill_manual("",values=site_colors)+
  scale_y_continuous("Proportion emerged",
                     limits = c(0.09, 0.73),
                     breaks = c(0.1 ,0.3, 0.5, 0.7))+
  scale_x_continuous(expression("GPP"*~midge^{-1}*""),
                     limits = c(0.0009, 0.0056),
                     breaks = c(0.001, 0.003, 0.005))+
  theme(strip.text = element_blank(),
        legend.position = c(0.275, 0.925),
        legend.direction = "horizontal",
        legend.text = element_text(margin = margin(l = -6)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.x = unit(0.4, "lines"),
        panel.border = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(margin = margin(r = 13)),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
p_adults





#==========
#========== Combine panels
#==========

p_comb_a <- plot_grid(NULL, p_adults, NULL,
                    nrow = 1,
                    rel_widths = c(0.01, 1,0.3))

# combine
p_comb <- plot_grid(NULL, p_gpp, NULL, NULL, p_comb_a,
                    nrow = 5,
                    rel_heights = c(0.05,1,0.15,0.025, 0.9),
                    # align = "h",
                    labels = c("",
                               "a",
                               "",
                               "",
                               "b"),
                    label_size = 12,
                    label_fontface = "plain",
                    hjust = c(0, 0, 0, 0),
                    vjust = c(0, 0, 0, -0.75))
p_comb
# ggsave(plot = p_comb, file = "manuscript/p_comb.pdf", width = 3.5, height = 6)





#==========
#========== Feedback
#==========

# calculate predicted gpp for observed levels
gpp_mod <- gpp_dat %>% 
  filter(midge_treat >= 50) %>% # take non-zero densities
  tidyr::expand(trial, site, midge_z = seq(min(midge_z), max(midge_z)+0.3, length.out = 100), 
                temp_z = 0, par_z = 0, rack = 1, coried = 1) %>%
  # back calculate midge_treat and midge_density from midge_z
  mutate(midge_treat = 2*midge_z*sd(gpp_dat$midge_treat) + mean(gpp_dat$midge_treat),
         midge_density = midge_treat*unique(gpp_dat$midge_density/gpp_dat$midge_treat)) %>%
  na.omit()
gpp_mod$pred <- predict(gpp_m, newdata = gpp_mod, re.form=NA) # predictions from model fit

# GPP per midge
gpp_mod_sum <- gpp_mod %>%
  filter(midge_treat > 0) %>%
  group_by(site, midge_treat, midge_density) %>%
  summarize(pred = mean(pred)) %>%
  mutate(pred_per = 1000 * pred/midge_density, 
         # calculate GPP per midge
         # convert production to mg O2 (x 1000) 
         pred_per_z = (pred_per - mean(adults_new$met_per))/sd(adults_new$met_per)) 
# z-score predictd production per midge using mean and sd from adults_new

# predicted emergence from GPP per midge using GLMM for adults
gpp_mod_pred <-predict(adults_gpp, 
                       newdata = gpp_mod_sum %>% mutate(met_per_z = pred_per_z),
                       re.form = NA, type = "response")
gpp_mod_sum$pred <- gpp_mod_pred



# calculate predicted GPP with midges removed (midge_treat = 0)
gpp_null <- gpp_dat %>% 
  filter(midge_treat == 0) %>%
  tidyr::expand(trial, site, midge_z, temp_z = 0, par_z = 0, rack = 1, coried = 1)
gpp_null$pred <- predict(gpp_m, newdata = gpp_null, re.form=NA)  # predictions from model fit

# GPP per midge
gpp_null_sum <- gpp_null %>%
  group_by(site) %>%
  summarize(pred = mean(pred)) %>%
  tidyr::expand(nesting(site, pred), 
                midge_density = gpp_mod_sum$midge_density) %>%
  mutate(pred_per = 1000 * pred/midge_density,
         # calculate GPP per midge
         # convert production to mg O2 (x 1000) produced over 13 days (13 x 24 / day)
         pred_per_z = (pred_per - mean(adults_new$met_per))/sd(adults_new$met_per))
# z-score predictd production per midge using mean and sd from adults_new

# predicted emergence from GPP per midge using GLMM for adults
gpp_null_pred <-predict(adults_gpp, 
                        newdata = gpp_null_sum %>% mutate(met_per_z = pred_per_z), 
                        re.form = NA, 
                        type = "response")
gpp_null_sum$pred <- gpp_null_pred


# plot proportion
p_feed_a <- gpp_mod_sum %>%
  ggplot(aes(midge_density/1000, pred, color = site))+
  facet_wrap(~site, nrow = 1)+
  geom_line(aes(y = pred), 
            size = 0.5)+
  geom_line(data = gpp_null_sum,  aes(y = pred), 
            size = 0.5, linetype = 2)+
  geom_text(inherit.aes = F,
            data = tibble(site = "E3",
                          x = c(75, 85),
                          y = c(0.135, 0.32),
                          label = c("no midge effect",
                                    "with midge effect")),
            aes(x = x, y = y, label = label),
            size = 2.5)+
  scale_color_manual("", values=site_colors)+
  scale_y_continuous("Proportion emerged",
                     limits = c(0.1, 0.4),
                     breaks = c(0.1, 0.2, 0.3, 0.4))+
  scale_x_continuous(name = "",
                     breaks = c(40, 80, 120),
                     limits = c(25, 130),
                     labels = NULL)+
  theme(strip.text = element_blank(),
        legend.position = c(0.2, 0.95),
        legend.direction = "horizontal",
        legend.text = element_text(margin = margin(l = -6)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.x = unit(0.4, "lines"),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        plot.margin = margin(1,1,1,1))+
  guides(color = guide_legend(override.aes = list(size = 0.917 * 0.5, fill = NA)))
p_feed_a

# plot total
p_feed_b <- gpp_mod_sum %>%
  mutate(pred = pred * midge_density / 1000) %>%
  ggplot(aes(midge_density/1000, pred, color = site))+
  facet_wrap(~site, nrow = 1)+
  geom_line(aes(y = pred), 
            size = 0.5)+
  geom_line(data = gpp_null_sum %>%
              mutate(pred = pred * midge_density / 1000),  
            aes(y = pred), 
            size = 0.5, linetype = 2)+
  geom_text(inherit.aes = F,
            data = tibble(site = "E3",
                          x = c(75, 85),
                          y = c(0.135, 0.3),
                          label = c("no midge effect",
                                    "with midge effect")),
            aes(x = x, y = y, label = label),
            size = 2.5)+
  scale_color_manual("", values=site_colors, guide = F)+
  scale_y_continuous(expression("Total emerged ("*1000~m^{-2}*")"),
                     limits = c(7.8, 20),
                     breaks = c(8, 12, 16, 20))+
  scale_x_continuous(expression("Larval density ("*1000~m^{-2}*")"),
                     breaks = c(40, 80, 120),
                     limits = c(25, 130))+
  theme(strip.text = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        plot.margin = margin(1,1,1,1))

p_feed_b

# combine
p_feed <- plot_grid(NULL, p_feed_a, NULL, p_feed_b,
                nrow = 4,
                rel_heights = c(0.05,1, 0.025, 1),
                align = "v",
                labels = c("",
                           "a",
                           "",
                           "b"),
                label_size = 12,
                label_fontface = "plain",
                hjust = c(0, 0, 0, 0),
                vjust = c(0, 0, 0, -0.75))

p_feed
# ggsave(plot = p_feed, file = "manuscript/p_feed.pdf", width = 3.5, height = 4.5)





#==========
#========== HOBO
#==========

inc_dates <- met_long$sampledate %>% unique()

hobo %>%
  group_by(sampledate) %>%
  mutate(mi = mean(temperature),
         hi = max(temperature) - mi,
         lo = mi - min(temperature)) %>%
  ungroup() %>%
  select(sampledate, mi, hi, lo) %>%
  unique()

hobo_prep <- hobo %>%
  group_by(sampledate) %>%
  filter(sampledate < "2017-07-21",
         ! (sampledate %in% inc_dates & 
              abs((temperature - mean(temperature))) > 1.76)) %>%
  left_join(met_long %>%
              select(rack, site) %>%
              unique()) %>%
  group_by(site, sampledate) %>%
  summarize(temp = mean(temperature),
            par = mean(light_intensity/54)) %>%
  mutate(yday = yday(sampledate),
         day = yday - min(yday)) %>%
  ungroup()

# map
sites <- tibble(lat = c(7278009, 7276128, 7276396 - 100),
                long = c(0406276, 0407635, 0409775),
                site = c("E2","E3","E5"))
anchor_scale = c(0403100, 7281500)
anchor_north = c(0405500, 7281000)
names(anchor_scale) = c("x","y")
names(anchor_north) = c("x","y")
p_a <- myv_df %>%
  ggplot(aes(long,lat))+
  geom_polygon(aes(fill = piece), size = 0.3, color = "black")+
  coord_equal()+
  theme_void() + 
  theme(plot.margin = margin(0,0,0,0))+
  geom_text(inherit.aes = F,
            data = sites,
            fontface = "bold",
            aes(x = long, y = lat, label = site, color = site),
            size = 4)+
  geom_text(aes(x = 0407800,
                y = 7273500),
            label = "X",
            size = 5,
            color = "black")+
  scale_color_manual("",values = site_colors, guide = F)+
  scale_fill_manual("",values = c("gray95",rep("white", 18)), guide = F)+
  scalebar(myv_df, dist = 2, dist_unit = "km",st.bottom = F,
           location = "topleft",st.size = 2.5, st.dist = 0.03,
           border.size = 0.2,anchor = anchor_scale,
           transform = F, model = "WGS84")+
  north(myv_df,symbol = 3,
        anchor = anchor_north)
p_a

par_u <- mean(hobo_prep$par)
par_s <- sd(hobo_prep$par)
temp_u <- mean(hobo_prep$temp)
temp_s <- sd(hobo_prep$temp)
x <- (c(0, 50, 100, 150) - par_u ) / par_s

c <- 1


hobo_prep %>% filter(yday %in% c(190, 198)) %>% select(yday, sampledate) %>% unique

hobo_prep %>% group_by(site) %>% summarize(temp = mean(temp))

p_b <- hobo_prep %>%
  mutate(par = (par - par_u) / par_s,
         temp = (temp - temp_u) / temp_s + c) %>%
  gather(var, val, temp, par) %>%
  ggplot(aes(yday, val, color = site))+
  facet_wrap(~var, labeller = labeller(var = element_blank()))+
  geom_line(size = 0.4)+
  scale_color_manual("",values = site_colors)+
  scale_y_continuous("PAR ("*mu*mol~m^{-2}~s^{-1}*")",
                     breaks = x,
                     labels = round(x * par_s + par_u, 0),
                     limits = c(-2.5, 3),
                     sec.axis = sec_axis(name = "Temperature"~(degree~C),
                                         trans = ~ (. - c) * temp_s + temp_u,
                                         breaks = c(11, 12, 13, 14)))+
  scale_x_continuous("Day of year",
                     limits = c(187, 201),
                     breaks = c(190, 198),
                     labels = c("9 July", "17 July"))+
  theme(strip.text = element_blank(),
        legend.position = c(0.5, 1),
        legend.direction = "horizontal",
        legend.text = element_text(margin = margin(l = -6)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.y = unit(0, "lines"),
        panel.spacing = unit(0.5, "lines"),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        plot.margin = margin(10,0,10,0))
p_b



p_ab <- plot_grid(p_a, NULL,p_b,
                  nrow = 3,
                  rel_heights = c(1.3, 0, 1),
                  ncol = 1,
                  labels = c("a",
                             "",
                             "b"),
                  label_size = 12,
                  label_fontface = "plain",
                  hjust = c(-0.5, 0, -0.5),
                  vjust = c(2.5,0,-1)
)
p_ab
# ggsave(plot = p_ab, file = "manuscript/p_sites.pdf", width = 3.5, height = 5)
