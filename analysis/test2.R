gpp_a <- lmer(met ~ temp_z + par_z + (trial + midge_z + site)^2 + I(midge_z^2) + I(midge_z^3)+
                (1|rack/core), 
              data = gpp_dat,
              REML = F)


gpp_b <- lmer(met ~ temp_z + par_z + (trial + midge_z + site)^2 + I(midge_z^2)+
                (1|rack/core), 
              data = gpp_dat,
              REML = F)

gpp_c <- lmer(met ~ temp_z + par_z + (trial + midge_z + site)^2+
                (1|rack/core), 
              data = gpp_dat,
              REML = F)

AIC(gpp_a, gpp_b, gpp_c)
list(gpp_a, gpp_b, gpp_c) %>% lapply(AICc) %>% unlist()


library(mgcv)

gpp_dat <- gpp_dat %>% mutate(Site = factor(site),
                              Trial = factor(trial),
                              Group = interaction(site, trial))

gpp_m <- gam(met ~ temp_z + par_z + Site + Trial + 
               s(midge_z, k = 3) + s(midge_z, k = 3, by = Group) ,
             data = gpp_dat)



# standardize observations
gpp_dat = gpp_dat %>%
  mutate(temp_cor = coef(gpp_m)[which(names(coef(gpp_m))=="temp_z")]*temp_z,
         met_z = met - temp_cor + mean(temp_cor))

# generate predicted data
gpp_nd <- gpp_dat %>% 
  tidyr::expand(nesting(Trial, Site, Group), midge_z, 
                temp_z = 0, par_z = 0) %>%
  left_join(gpp_dat %>% select(midge_z, midge_density) %>% unique)
gpp_preds <- predict(gpp_m, newdata = gpp_nd, se.fit = T)
gpp_nd$pred <- gpp_preds$fit
gpp_nd$se <- gpp_preds$se


# plot
gpp_nd %>%
  ggplot(aes(midge_density/1000,pred, color = Site))+
  facet_wrap(~Trial)+
  geom_line(aes(y = pred),
            size = 1)+
  geom_point(aes(y = pred),
             size = 3)+
  geom_errorbar(aes(ymin = pred - se, ymax = pred + se), width = 0, size = 1)+
  geom_point(data = gpp_dat,
             aes(y = met_z), size = 2, alpha = 0.5)+
  scale_color_manual("Site",values=site_colors)+
  scale_fill_manual("Site",values=site_colors)+
  scale_y_continuous(expression("Gross Primary Production (g "*O[2]~m^{-2}~h^{-1}*")"))+
  scale_x_continuous(expression("Midge Density ("*1000~m^{-2}*")"),
                     breaks = c(15,50,85))

