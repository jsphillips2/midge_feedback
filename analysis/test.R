# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=20, margin=margin(0,0,12,0)),
                  legend.text = element_text(size = 20),
                  legend.title = element_text(size = 20),
                  axis.text=element_text(size=20, color="black"),
                  axis.title = element_text(size = 20),
                  axis.title.y=element_text(angle = 90 ,margin=margin(0,10,0,0)),
                  axis.title.x=element_text(margin=margin(10,0,0,0)),
                  panel.border = element_rect(size = 1.25, fill = NA)))

gpp_nd %>%
  mutate(trial = factor(trial, levels = c(0,1), labels = c("Day 7", "Day 20"))) %>%
  ggplot(aes(midge_density/1000,pred, color = site))+
  facet_wrap(~trial)+
  geom_line(aes(y = pred),
            size = 1.25)+
  geom_point(aes(y = pred),
             size = 4)+
  geom_errorbar(aes(ymin = pred - se, ymax = pred + se), width = 0, size = 1)+
  geom_point(data = gpp_dat%>%
               mutate(trial = factor(trial, levels = c(0,1), labels = c("Day 7", "Day 20"))), 
             aes(y = met_z), size = 3, alpha = 0.5)+
  scale_color_manual("Site",values=site_colors)+
  scale_fill_manual("Site",values=site_colors)+
  scale_y_continuous(expression("Gross Primary Production (g "*O[2]~m^{-2}~h^{-1}*")"))+
  scale_x_continuous(expression("Midge Density ("*1000~m^{-2}*")"),
                     breaks = c(15,50,85))+
  theme(legend.position = c(0.1,0.8))

# ggsave("analysis/figures/gpp.jpg", units = "in", width = 7.5, height = 5.5, device = "jpg", dpi = 900)


adults_gpp_nd %>%
  ggplot(aes(met_per, pred))+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se),
              alpha = 0.15, linetype = 0)+
  geom_line(aes(y = pred), size = 1)+
  geom_point(data = adults_new, aes(y = prop, color = site), size = 3.5, alpha = 0.7)+
  scale_color_manual("Site", values=site_colors)+
  scale_fill_manual("Site",values=site_colors)+
  scale_y_continuous("Proportion Emerged")+
  scale_x_continuous(expression("GPP"*~midge^{-1}*""))+
  theme(legend.position = c(0.1,0.85))

# ggsave("analysis/figures/emg.jpg", units = "in", width = 5.5, height = 5.5, device = "jpg", dpi = 900)

gpp_mod_sum %>%
  filter(site == "E2") %>%
  ggplot(aes(midge_density/1000, pred, color = site))+
  facet_wrap(~site, nrow = 1)+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se,fill = site),
              alpha = 0.15, linetype = 0)+
  geom_line(aes(y = pred), size = 1.25)+
  geom_line(data = gpp_null_sum %>%
              filter(site == "E2"),  aes(y = pred), size = 1.25, linetype = 2)+
  geom_ribbon(data = gpp_null_sum %>%
                filter(site == "E2"),aes(ymin = pred - se, ymax = pred + se,fill = site),
              alpha = 0.15, linetype = 0)+
  scale_color_manual("Site", values=site_colors, guide = F)+
  scale_fill_manual("Site", values=site_colors, guide = F)+
  scale_y_continuous("Proportion Emerged", limits = c(0.08,0.42), 
                     breaks = c(0.1,0.2,0.3,0.4))+
  scale_x_continuous(expression("Initial Midge Density ("*1000~m^{-2}*")"),
                     limits = c(25, 130),
                     breaks = c(40,80,120))

# ggsave("analysis/figures/feed_a.jpg", units = "in", width = 5, height = 5.5, device = "jpg", dpi = 900)


gpp_mod_sum %>%
  ggplot(aes(midge_density/1000, pred, color = site))+
  facet_wrap(~site, nrow = 1)+
  geom_ribbon(aes(ymin = pred - se, ymax = pred + se,fill = site),
              alpha = 0.15, linetype = 0)+
  geom_line(aes(y = pred), size = 1.25)+
  geom_line(data = gpp_null_sum,  aes(y = pred), size = 1.25, linetype = 2)+
  geom_ribbon(data = gpp_null_sum,aes(ymin = pred - se, ymax = pred + se,fill = site),
              alpha = 0.15, linetype = 0)+
  scale_color_manual("Site", values=site_colors, guide = F)+
  scale_fill_manual("Site", values=site_colors, guide = F)+
  scale_y_continuous("Proportion Emerged", limits = c(0.08,0.42), 
                     breaks = c(0.1,0.2,0.3,0.4))+
  scale_x_continuous(expression("Initial Midge Density ("*1000~m^{-2}*")"),
                     limits = c(25, 130),
                     breaks = c(40,80,120))

# ggsave("analysis/figures/feed.jpg", units = "in", width = 8.5, height = 5.5, device = "jpg", dpi = 900)



hobo %>%
  filter(sampledate < "2017-07-21") %>%
  left_join(met_long %>%
              select(rack, site) %>%
              unique()) %>%
  group_by(site, sampledate) %>%
  summarize(temp = mean(temperature),
            par = mean(light_intensity/54)) %>%
  ggplot(aes(sampledate, temp, color = site))+
  geom_line(size = 1.25)+
  scale_y_continuous("Temperature"~(degree~C),
                     breaks = c(11,12,13,14), limits = c(10.5,14.5))+
  scale_x_date("",breaks = 
                 lubridate::as_date(c("2017-07-08","2017-07-13","2017-07-18")),
               labels = c("08 June", "13 June", "18 June"))+
  scale_color_manual("Site", values=site_colors)+
  scale_fill_manual("Site", values=site_colors)+
  theme(legend.position = c(0.8,0.2))

# ggsave("analysis/figures/temp.jpg", units = "in", width = 7.5, height = 5.5, device = "jpg", dpi = 900)

tt <- read_csv("analysis/midges.csv")

tt %>% 
  gather(var, val, -year) %>%
  ggplot(aes(year, val + 1))+
  geom_line(size = 1.25)+
  scale_y_continuous("Tanytarsus catch",
                     trans = "log", breaks = c(10,1000,100000))+
  scale_x_continuous("Year")

# ggsave("analysis/figures/tanyt.jpg", units = "in", width = 7.5, height = 5.5, device = "jpg", dpi = 900)

