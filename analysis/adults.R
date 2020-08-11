#==========
#========== Adults
#==========

# select final data of adults
adults_0 <- adults  %>%
  mutate(duration = as.numeric(max(sampledate) - min(sampledate)))%>%
  filter(sampledate == max(sampledate)) %>%
  filter(midge_treat >= 50)

adults_m <- glmer(cbind(cumem, midge_treat - cumem) ~ midge_z*site + (1|rack/core),
                  data = adults_0, family = "binomial")

# p-values
Anova(adults_m, type = 3)
Anova(update(adults_m, .~.-midge_z:site))

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
              summarize(met_z = mean(met_z),
                        midge_density = unique(midge_density)) %>%
              mutate(met_per = 1000 * 24 * 13 * met_z / midge_density,
                     met_per_z = (met_per - mean(met_per))/sd(met_per))) 
# fit_model
adults_gpp <- glmer(cbind(cumem, midge_treat - cumem) ~ met_per_z + (1|rack/core),
                    data = adults_new, family = "binomial")

# p-values
Anova(adults_gpp)

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
  scale_x_continuous(expression("mg"~O[2]*~midge^{-1}*""),
                     limits = c(0.25, 1.75),
                     breaks = c(0.4, 1, 1.6))+
  theme(strip.text = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.direction = "horizontal",
        legend.text = element_text(margin = margin(l = -6)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.x = unit(0.4, "lines"),
        panel.border = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", bottom='both')
p_adults
# ggsave(file = "analysis/figures/p_adults.pdf",
#           width = 3.5, height = 3.5)