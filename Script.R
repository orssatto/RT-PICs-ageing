# Libraries
library(readxl)
library(naniar)
library(visdat)
library(dplyr)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(janitor)
library(cowplot)
library(emmeans)
library(lme4)
library(car)
library(rmcorr)
library(psycho)
library(sjstats)
library(pwr)
library(viridis)

#### 20% triangular contractions unmatched motor units
  ## Load dataset "d20". Variables: x = time (c1,c2, post) and y = deltaf, peak_dr, rt

d20 = read.csv("Delta F soleus 20pc_ALL_MUS_and_PAIRS_stats.csv") %>%
  clean_names() %>%
  mutate(
    time = as.factor(time),
    participant = as.factor(participant) 
  ) %>%

    # Model: deltaf
    fit20 <- lmer(deltaf ~ as.factor(time) + (1 | participant), data = d20)

    # Model diagnostics
    #plot(fit_20)
    qqnorm(residuals(fit20)); qqline(residuals(fit20))
    qqPlot(residuals(fit20))
    summary(fit20)
    
    #confint(fit_1way)
    anova_fit20 <- anova(fit20)
    omega_sq(anova_fit20)
    anova_fit20
    
    #emmeans pairwise d20
    fit20.emm.s <- emmeans(fit20, "time")
    pairs(fit20.emm.s, adjust = "bonferroni")
    
    # Fitted values
    (refgrid <- list (time=c("c1","c2","post")))
    mar_ef20 <- emmip(fit20, ~ time, at = refgrid, CIs = T, plotit = F)
    mar_ef20
    
        #plot deltaf 20
    ggplot(data = d20, aes(x = time, y = deltaf)) +
      geom_jitter(width = 0.15, alpha = 1, size = 3.5, aes(colour = participant)) +
      scale_color_viridis(discrete = TRUE, option = "D") +
      theme_bw(base_size = 14) +
      guides(fill = 'none', color = 'none') +
      geom_point(data = mar_ef20, aes(x = time, y = yvar), 
      position = position_nudge(x = -0.3), size = 3.5) +
      geom_errorbar(data = mar_ef20, aes(ymin = LCL, ymax = UCL, y = yvar),
      position = position_nudge(x = -0.3), width = 0, size = 2.0) +
      scale_x_discrete(breaks=c("c1", "c2", "post"), labels=c("-2 weeks", "0 weeks", "+6 weeks")) +
      ylim(-2.5,6.5) +
      theme(
        axis.text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      ) +
      labs(y = "Δ F (pps)", x = "Time")+
      facet_grid(~"Unmatched Motor Units") -> plot_delta_f
    plot_delta_f
    #ggsave(file = "deltaf.png", units="in", width = 7, height = 4.5, dpi = 300)
  
    # model: peak discharge rates (test units only)
    fit20dr <- lmer(peak_dr ~ as.factor(time) + (1 | participant), data = d20)  
    
    # Model diagnostics
    #plot(fit_20dr)
    qqnorm(residuals(fit20dr)); qqline(residuals(fit20dr))
    qqPlot(residuals(fit20dr))
    summary(fit20dr)
    
    #confint(fit_1way)
    anova_fit20dr <- anova(fit20dr)
    anova_fit20dr
    omega_sq(anova_fit20dr)
    
    #emmeans pairwise d20
    fit20dr.emm.s <- emmeans(fit20dr, "time")
    pairs(fit20dr.emm.s, adjust = "bonferroni")    
    
    ## Fitted values
    (refgrid <- list (time=c("c1","c2","post")))
    mar_ef20dr <- emmip(fit20dr, ~ time, at = refgrid, CIs = T, plotit = F)
    mar_ef20dr
    
    # model: rt (test units only)
    fit20rt <- lmer(rt ~ as.factor(time) + (1 | participant), data = d20)  
    
    # Model diagnostics
    #plot(fit_20rt)
    qqnorm(residuals(fit20rt)); qqline(residuals(fit20rt))
    qqPlot(residuals(fit20rt))
    summary(fit20rt)
    
    #confint(fit_1way)
    anova_fit20rt <- anova(fit20rt)
    omega_sq(anova_fit20rt)
    anova_fit20rt
    
    #emmeans pairwise d20
    fit20rt.emm.s <- emmeans(fit20rt, "time")
    pairs(fit20rt.emm.s, adjust = "bonferroni")   

    ## Fitted values
    (refgrid <- list (time=c("c1","c2","post")))
    mar_ef20rt <- emmip(fit20rt, ~ as.factor(time), at = refgrid, CIs = T, plotit = F)
    mar_ef20rt
    
    
    #### 20% triangular contractions matched motor units

# Load d20t data. Variables c1 vs c2: x = timecc (c1, c2) and y = deltafcc. 
                 #Variable control vs post: x = timecp (control, post) and y = deltafcp
    
    d20t = read_csv("Delta F soleus 20perc_TRACKED_MUs.csv") %>%
      clean_names() %>%
      mutate(
        timecc = as.factor(timecc),
        timecp = as.factor(timecp),
        participant = as.factor(participant))
      
  # Model: deltafcc, peak_dr_cc, rt_cc 20 (2-weeks Pre-training and Pre-training)
  fit20t_cc <- lmer(deltafcc ~ as.factor(timecc) + (1 | participant/mu_list), data = d20t)
  fit20t_dr_cc <- lmer(peak_dr_cc ~ as.factor(timecc) + (1 | participant/mu_list), data = d20t) 
  fit20t_rt_cc <- lmer(rt_cc ~ as.factor(timecc) + (1 | participant/mu_list), data = d20t)
  
  # Model diagnostics
  #plot(fit_20t_cc)
  qqnorm(residuals(fit20t_cc)); qqline(residuals(fit20t_cc))
  qqPlot(residuals(fit20t_cc))
  qqPlot(residuals(fit20t_dr_cc))
  qqPlot(residuals(fit20t_rt_cc))
  summary(fit20t_cc)
  summary(fit20t_dr_cc)
  summary(fit20t_rt_cc)
  
  #confint(fit_1way)
  anova(fit20t_cc)
  anova(fit20t_dr_cc)
  anova(fit20t_rt_cc)
  
  #emmeans pairwise d20t
  fit20t_cc.emm.s <- emmeans(fit20t_cc, "timecc")
  pairs(fit20t_cc.emm.s, adjust = "bonferroni")
  
  fit20t_dr_cc.emm.s <- emmeans(fit20t_dr_cc, "timecc")
  pairs(fit20t_dr_cc.emm.s, adjust = "bonferroni")
  
  fit20t_rt_cc.emm.s <- emmeans(fit20t_rt_cc, "timecc")
  pairs(fit20t_rt_cc.emm.s, adjust = "bonferroni")
  
  ## Fitted values
  (refgrid <- list (time=c("c1","c2")))
  mar_ef20t_cc <- emmip(fit20t_cc, ~ as.factor(timecc), at = refgrid, CIs = T, plotit = F)
  mar_ef20t_cc
  
  (refgrid <- list (time=c("c1","c2")))
  mar_ef20t_dr_cc <- emmip(fit20t_dr_cc, ~ as.factor(timecc), at = refgrid, CIs = T, plotit = F)
  mar_ef20t_dr_cc
  
  (refgrid <- list (time=c("c1","c2")))
  mar_ef20t_rt_cc <- emmip(fit20t_rt_cc, ~ as.factor(timecc), at = refgrid, CIs = T, plotit = F)
  mar_ef20t_rt_cc
  
  # Model: Delta F, peak discharge rates and recruitment threshold 20% triangular contractions 
  #(2-weeks Pre-training + Pre-training vs Post-training) 
  fit20t_cp <- lmer(deltafcp ~ as.factor(timecp) + (1 | participant/mu_list), data = d20t)
  fit20t_dr_cp <- lmer(peak_dr_cp ~ as.factor(timecp) + (1 | participant/mu_list), data = d20t) 
  fit20t_rt_cp <- lmer(rt_cp ~ as.factor(timecp) + (1 | participant/mu_list), data = d20t)
  
  
  # Model diagnostics
  #plot(fit_20t_cp)
  qqnorm(residuals(fit20t_cp)); qqline(residuals(fit20t_cp))
  qqPlot(residuals(fit20t_cp))
  summary(fit20t_cp)
  
  qqnorm(residuals(fit20t_dr_cp)); qqline(residuals(fit20t_dr_cp))
  qqPlot(residuals(fit20t_dr_cp))
  summary(fit20t_dr_cp)
  
  qqnorm(residuals(fit20t_rt_cp)); qqline(residuals(fit20t_rt_cp))
  qqPlot(residuals(fit20t_rt_cp))
  summary(fit20t_rt_cp)
  
  #confint(fit_1way)
  anova(fit20t_cp)
  anova(fit20t_dr_cp)
  anova(fit20t_rt_cp)
  
  
  #emmeans pairwise d20t
  fit20t_cp.emm.s <- emmeans(fit20t_cp, "timecp")
  pairs(fit20t_cp.emm.s, adjust = "bonferroni")
  
  fit20t_dr_cp.emm.s <- emmeans(fit20t_dr_cp, "timecp")
  pairs(fit20t_dr_cp.emm.s, adjust = "bonferroni")
  
  fit20t_rt_cp.emm.s <- emmeans(fit20t_rt_cp, "timecp")
  pairs(fit20t_rt_cp.emm.s, adjust = "bonferroni")
  
  ## Fitted values
  (refgrid <- list (time=c("control","post")))
  mar_ef20t_cp <- emmip(fit20t_cp, ~ as.factor(timecp), at = refgrid, CIs = T, plotit = F)
  mar_ef20t_cp
  
  (refgrid <- list (time=c("control","post")))
  mar_ef20t_dr_cp <- emmip(fit20t_dr_cp, ~ as.factor(timecp), at = refgrid, CIs = T, plotit = F)
  mar_ef20t_dr_cp
  
  (refgrid <- list (time=c("control","post")))
  mar_ef20t_rt_cp <- emmip(fit20t_rt_cp, ~ as.factor(timecp), at = refgrid, CIs = T, plotit = F)
  mar_ef20t_rt_cp
  
    #plot delta f 2-weeks Pre-training + Pre-training vs Post-training
  ggplot(data = d20t, aes(x = timecp, y = deltafcp)) +
    geom_jitter(width = 0.02, alpha = 1, size = 3.5, aes(colour = participant)) +
    scale_color_viridis(discrete = TRUE, option = "D") +
    theme_bw(base_size = 14) +
    #(~time) +
    guides(fill = F, color = F) +
    geom_line(data = d20t, aes(x = timecp, y = deltafcp, group = mu_list),color='grey50', size = .75, alpha = .25) +
    geom_point(data = mar_ef20t_cp, aes(x = timecp, y = yvar), 
               position = position_nudge(x = -0.15), size = 3.5) +
    geom_errorbar(data = mar_ef20t_cp, aes(ymin = LCL, ymax = UCL, y = yvar),
                  position = position_nudge(x = -0.15), width = 0.0, size = 2.0) +
    ylim(-2.5,6.5) +
    scale_x_discrete(breaks=c("control", "post"), labels=c("-2 and 0 weeks", "+6 weeks")) +
    theme(
      axis.text = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20),
      strip.text.x = element_text(size = 20),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(y = "Δ F (pps)", x = "Time") +
    facet_grid(~"Matched Motor Units") -> plot_deltaf_tcp
  plot_deltaf_tcp
  #ggsave(file = "peakdr.png", units="in", width = 7, height = 4.5, dpi = 300)
  
  
#########################################################################################################  

#### 20% trapezoidal and 40% triangular contractions
  # Load dataset d40 (2 factors): time(c1,c2,post) by contraction intensity(20%,40%)

  d40 = read_csv("delta_f_20_40.csv") %>%
  mutate(
    time = as.factor(time),
    participant = as.factor(participant),
    intensity = as.factor(intensity)
  )

# Model: deltaf 20 vs 40
fit2040 <- lmer(deltaf ~ as.factor(time)*as.factor(intensity) +  (1 | participant/mu_list), data = d40)

# Model diagnostics
#plot(fit2040)
qqnorm(residuals(fit20)); qqline(residuals(fit20))
qqPlot(residuals(fit2040))
summary(fit2040)

#confint(fit_2ways)
anova(fit2040)
avova_2040 <- anova(fit2040)
omega_sq(anova(fit2040))

#emmeans pairwise d2040
fit2040.emm.s <- emmeans(fit2040, "intensity", "time")
pairs(fit2040.emm.s, adjust = "bonferroni")

## Fitted values
(refgrid <- list(intensity=c("20","40"), time=c("c1","c2","post")))
mar_d40 <- emmip(fit2040.emm.s, ~ intensity|time, at = refgrid, CIs = T, plotit = F)
mar_d40

labels <- c(c1 = "-2 weeks", c2 = "0 weeks", post = "+6 weeks")


# plot d40
ggplot(data = d40, aes(x = as.factor(intensity), y = deltaf)) +
  geom_jitter(width = 0.05, alpha = 1, size = 3.5, aes(colour = participant)) +
  scale_color_viridis(discrete = TRUE, option = "D") +
  theme_bw(base_size = 14) +
  facet_grid(~time, labeller = labeller(time = labels)) +
  guides(fill = F, color = F) +
  geom_line(data = d40, aes(x = as.factor(intensity), y = deltaf, group = mu_list), size = .75, color='grey50', width = 0.05, alpha = .25) +
  geom_point(data = mar_d40, aes(x = as.factor(intensity), y = yvar), 
  position = position_nudge(x = -0.2), width = 0, size = 3.5) +
  geom_errorbar(data = mar_d40, aes(ymin = LCL, ymax = UCL, y = yvar),
  position = position_nudge(x = -0.2), width = 0, size = 2.0) +
  scale_x_discrete(breaks=c("20", "40"), labels=c("20%", "40%")) +
  theme(
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  scale_y_continuous(n.breaks = 8) +
  labs(y = "Δ F (pps)", x = "Contraction Intensity") -> plot_delta_f2040
plot_delta_f2040
ggsave(file = "plot_deltaf2040.png", units="in", width = 11.6, height = 9, dpi = 300)



#### Mean differences calculations

  #Mean difference d20 - 20% triangular contraction - Delta f unmatched motor units
emm <- emmeans(fit20, pairwise ~ time)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

md  <- ci_95$contrasts$estimate[1]
lower_90 <- ci_90$contrasts$lower.CL[1]
upper_90 <- ci_90$contrasts$upper.CL[1]
lower_95 <- ci_95$contrasts$lower.CL[1]
upper_95 <- ci_95$contrasts$upper.CL[1]
time <- 'c1-c2'

df_20c1c2 <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20c1c2

md  <- ci_95$contrasts$estimate[2]
lower_90 <- ci_90$contrasts$lower.CL[2]
upper_90 <- ci_90$contrasts$upper.CL[2]
lower_95 <- ci_95$contrasts$lower.CL[2]
upper_95 <- ci_95$contrasts$upper.CL[2]
time <- 'c1-post'

df_20c1post <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20c1post

md  <- ci_95$contrasts$estimate[3]
lower_90 <- ci_90$contrasts$lower.CL[3]
upper_90 <- ci_90$contrasts$upper.CL[3]
lower_95 <- ci_95$contrasts$lower.CL[3]
upper_95 <- ci_95$contrasts$upper.CL[3]
time <- 'c2-post'

timecontrast <- c('c1-c2','c1-post','c2-post')

df_20c2post <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20c2post

#Join d20
d20_md <- rbind(df_20c1c2, df_20c1post, df_20c2post) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
    time = recode_factor(time, 'c1-c2' = 'c1-c2', 'c1-post' = 'c1-post', 'c2-post' = 'c2-post')
        )

#plot md d20
ggplot(data = d20_md, aes(x = time, y = md)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), size = 0.8, width = 0) +
  geom_errorbar(aes(ymin = lower_90, ymax = upper_90), size = 2.0, width = 0) +
  theme_bw(base_size = 14) +
  labs(x = "Time", y = "ΔF Mean Difference (pps)") +
  geom_hline(yintercept=0, linetype="solid",color = "red", size=0.8) +
  #ylim(-1.0,0.5) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  facet_grid(~"UNMATCHED MOTOR UNITS") -> plot_md_20
plot_md_20

#Effect size d - Delta F unmatched motor units
conditional_effect <- emmeans(fit20, ~ time, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit20), edf = df.residual(fit20))

#Mean difference d20 - 20% triangular contraction - Peak discharge rates unmatched motor units
emm <- emmeans(fit20dr, pairwise ~ time)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

md  <- ci_95$contrasts$estimate[1]
lower_90 <- ci_90$contrasts$lower.CL[1]
upper_90 <- ci_90$contrasts$upper.CL[1]
lower_95 <- ci_95$contrasts$lower.CL[1]
upper_95 <- ci_95$contrasts$upper.CL[1]
time <- 'c1-c2'

df_20c1c2 <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20c1c2

md  <- ci_95$contrasts$estimate[2]
lower_90 <- ci_90$contrasts$lower.CL[2]
upper_90 <- ci_90$contrasts$upper.CL[2]
lower_95 <- ci_95$contrasts$lower.CL[2]
upper_95 <- ci_95$contrasts$upper.CL[2]
time <- 'c1-post'

df_20c1post <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20c1post

md  <- ci_95$contrasts$estimate[3]
lower_90 <- ci_90$contrasts$lower.CL[3]
upper_90 <- ci_90$contrasts$upper.CL[3]
lower_95 <- ci_95$contrasts$lower.CL[3]
upper_95 <- ci_95$contrasts$upper.CL[3]
time <- 'c2-post'

timecontrast <- c('c1-c2','c1-post','c2-post')

df_20c2post <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20c2post

#Join d20
d20_md <- rbind(df_20c1c2, df_20c1post, df_20c2post) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
    time = recode_factor(time, 'c1-c2' = 'c1-c2', 'c1-post' = 'c1-post', 'c2-post' = 'c2-post')
  )

#plot md d20 dr
ggplot(data = d20_md, aes(x = time, y = md)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), size = 0.8, width = 0) +
  geom_errorbar(aes(ymin = lower_90, ymax = upper_90), size = 2.0, width = 0) +
  theme_bw(base_size = 14) +
  labs(x = "Time", y = "Mean Difference\n(pps)") +
  geom_hline(yintercept=0, linetype="solid",color = "red", size=0.8) +
  #ylim(-1.0,0.5) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  facet_grid(~"Peak Discharge Rates") -> plot_md_20
plot_md_20

#Effect size d - Peak discharge rates unmatched motor units
conditional_effect <- emmeans(fit20dr, ~ time, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit20dr), edf = df.residual(fit20dr))

#Mean difference d20 - 20% triangular contraction - recruitment thresholds unmatched motor units
emm <- emmeans(fit20rt, pairwise ~ time)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

md  <- ci_95$contrasts$estimate[1]
lower_90 <- ci_90$contrasts$lower.CL[1]
upper_90 <- ci_90$contrasts$upper.CL[1]
lower_95 <- ci_95$contrasts$lower.CL[1]
upper_95 <- ci_95$contrasts$upper.CL[1]
time <- 'c1-c2'

df_20c1c2 <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20c1c2

md  <- ci_95$contrasts$estimate[2]
lower_90 <- ci_90$contrasts$lower.CL[2]
upper_90 <- ci_90$contrasts$upper.CL[2]
lower_95 <- ci_95$contrasts$lower.CL[2]
upper_95 <- ci_95$contrasts$upper.CL[2]
time <- 'c1-post'

df_20c1post <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20c1post

md  <- ci_95$contrasts$estimate[3]
lower_90 <- ci_90$contrasts$lower.CL[3]
upper_90 <- ci_90$contrasts$upper.CL[3]
lower_95 <- ci_95$contrasts$lower.CL[3]
upper_95 <- ci_95$contrasts$upper.CL[3]
time <- 'c2-post'

timecontrast <- c('c1-c2','c1-post','c2-post')

df_20c2post <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20c2post

#Join d20
d20_md <- rbind(df_20c1c2, df_20c1post, df_20c2post) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
    time = recode_factor(time, 'c1-c2' = 'c1-c2', 'c1-post' = 'c1-post', 'c2-post' = 'c2-post')
  )

#Effect size d - Recruitment threshold unmatched motor units
conditional_effect <- emmeans(fit20rt, ~ time, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit20rt), edf = df.residual(fit20rt))


#Mean difference d20t - 20% triangular contraction - Delta F matched motor units (2-weeks pre-training vs Pre-training)
emm <- emmeans(fit20t_cc, pairwise ~ timecc)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

md  <- ci_95$contrasts$estimate[1]
lower_90 <- ci_90$contrasts$lower.CL[1]
upper_90 <- ci_90$contrasts$upper.CL[1]
lower_95 <- ci_95$contrasts$lower.CL[1]
upper_95 <- ci_95$contrasts$upper.CL[1]
time <- 'c1-c2'

df_20tcc <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20tcc

#Join d20t control vs control
d20_tcc <- rbind(df_20tcc) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
  )

#Effect size d - Delta F tracked motor units
conditional_effect <- emmeans(fit20t_cc, ~ timecc, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit20t_cc), edf = df.residual(fit20t_cc))


#Mean difference d20t - 20% triangular contraction - Delta F matched motor units (2-weeks pre-training + Pre-training vs Post-training)
emm <- emmeans(fit20t_cp, pairwise ~ timecp)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

md  <- ci_95$contrasts$estimate[1]
lower_90 <- ci_90$contrasts$lower.CL[1]
upper_90 <- ci_90$contrasts$upper.CL[1]
lower_95 <- ci_95$contrasts$lower.CL[1]
upper_95 <- ci_95$contrasts$upper.CL[1]
time <- 'control-post'

df_20tcp <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20tcp

#Join d20t control vs post
d20_tcp <- rbind(df_20tcp) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
    )

#Effect size d
conditional_effect <- emmeans(fit20t_cp, ~ timecp, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit20t_cp), edf = df.residual(fit20t_cp))


#Mean difference d20t - 20% triangular contraction - Peak discharge rates matched motor units (2-weeks pre-training vs Pre-training)
emm <- emmeans(fit20t_dr_cc, pairwise ~ timecc)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

md  <- ci_95$contrasts$estimate[1]
lower_90 <- ci_90$contrasts$lower.CL[1]
upper_90 <- ci_90$contrasts$upper.CL[1]
lower_95 <- ci_95$contrasts$lower.CL[1]
upper_95 <- ci_95$contrasts$upper.CL[1]
time <- 'c1-c2'

df_20tcc_dr <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20tcc_dr

#Join d20t control vs control
d20_tcc <- rbind(df_20tcc_dr) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
  )

#Effect size d
conditional_effect <- emmeans(fit20t_dr_cc, ~ timecc, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit20t_dr_cc), edf = df.residual(fit20t_dr_cc))


#Mean difference d20t - 20% triangular contraction - Peak discharge rates matched motor units (2-weeks pre-training + Pre-training vs Post-training)
emm <- emmeans(fit20t_dr_cp, pairwise ~ timecp)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

md  <- ci_95$contrasts$estimate[1]
lower_90 <- ci_90$contrasts$lower.CL[1]
upper_90 <- ci_90$contrasts$upper.CL[1]
lower_95 <- ci_95$contrasts$lower.CL[1]
upper_95 <- ci_95$contrasts$upper.CL[1]
time <- 'control-post'

df_20tcp_dr <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20tcp_dr

#Join d20t control vs post
d20_tcp <- rbind(df_20tcp_dr) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
  )

#Effect size d
conditional_effect <- emmeans(fit20t_dr_cp, ~ timecp, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit20t_dr_cp), edf = df.residual(fit20t_dr_cp))

#Mean difference d20t - 20% triangular contraction - Recruitment Threshold matched motor units (2-weeks pre-training vs Pre-training)
emm <- emmeans(fit20t_rt_cc, pairwise ~ timecc)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

md  <- ci_95$contrasts$estimate[1]
lower_90 <- ci_90$contrasts$lower.CL[1]
upper_90 <- ci_90$contrasts$upper.CL[1]
lower_95 <- ci_95$contrasts$lower.CL[1]
upper_95 <- ci_95$contrasts$upper.CL[1]
time <- 'c1-c2'

df_20tcc_rt <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20tcc_rt

#Join d20t control vs control
d20_tcc <- rbind(df_20tcc_rt) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
  )

#Effect size d
conditional_effect <- emmeans(fit20t_rt_cc, ~ timecc, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit20t_rt_cc), edf = df.residual(fit20t_rt_cc))

#Mean difference d20t - 20% triangular contraction - Recruitment Threshold matched motor units (2-weeks pre-training + Pre-training vs Post)
emm <- emmeans(fit20t_rt_cp, pairwise ~ timecp)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

md  <- ci_95$contrasts$estimate[1]
lower_90 <- ci_90$contrasts$lower.CL[1]
upper_90 <- ci_90$contrasts$upper.CL[1]
lower_95 <- ci_95$contrasts$lower.CL[1]
upper_95 <- ci_95$contrasts$upper.CL[1]
time <- 'control-post'

df_20tcp_rt <- cbind(time,md,lower_90,upper_90,lower_95,upper_95)
df_20tcp_rt

#Join d20t control vs post
d20_tcp <- rbind(df_20tcp_rt) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
  )

#Effect size d
conditional_effect <- emmeans(fit20t_rt_cp, ~ timecp, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit20t_rt_cp), edf = df.residual(fit20t_rt_cp))

#mean difference d40
emm <- emmeans(fit2040, pairwise ~ intensity|time)
summary(emm)
ci_95 <- confint(emm, level = .95)
ci_90 <- confint(emm, level = .90)

#md c1
md  <- ci_95$contrasts$estimate[1]*-1
lower_90 <- ci_90$contrasts$lower.CL[1]*-1
upper_90 <- ci_90$contrasts$upper.CL[1]*-1
lower_95 <- ci_95$contrasts$lower.CL[1]*-1
upper_95 <- ci_95$contrasts$upper.CL[1]*-1
time <- 'c1'

df_c1 <- cbind(time, md,lower_90,upper_90,lower_95,upper_95)
df_c1

#md c2
md  <- ci_95$contrasts$estimate[2]*-1
lower_90 <- ci_90$contrasts$lower.CL[2]*-1
upper_90 <- ci_90$contrasts$upper.CL[2]*-1
lower_95 <- ci_95$contrasts$lower.CL[2]*-1
upper_95 <- ci_95$contrasts$upper.CL[2]*-1
time <- 'c2'

df_c2 <- cbind(time, md,lower_90,upper_90,lower_95,upper_95)
df_c2

#md post
md  <- ci_95$contrasts$estimate[3]*-1
lower_90 <- ci_90$contrasts$lower.CL[3]*-1
upper_90 <- ci_90$contrasts$upper.CL[3]*-1
lower_95 <- ci_95$contrasts$lower.CL[3]*-1
upper_95 <- ci_95$contrasts$upper.CL[3]*-1
time <- 'post'

df_post <- cbind(time, md,lower_90,upper_90,lower_95,upper_95)
df_post

#Join c1, c2, post md
d40_md <- rbind(df_c1, df_c2, df_post) %>%
  as_tibble() %>%
  mutate(
    md = as.numeric(md),
    lower_90 = as.numeric(lower_90),
    lower_95 = as.numeric(lower_95),
    upper_90 = as.numeric(upper_90),
    upper_95 = as.numeric(upper_95),
    time = recode_factor(time, 'c1' = 'c1', 'c2' = 'c2', 'post' = 'post')
  )

#plot md d40
ggplot(data = d40_md, aes(x = time, y = md)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), size = 2.0, width = 0) +
  theme_bw(base_size = 14) +
  labs(x = NULL, y = "Δ F 40%-20% (pps)")+
  geom_hline(yintercept=0, linetype="solid",color = "red", size=0.8) +
  scale_x_discrete(breaks=c("c1", "c2", "post"), labels=c("-2 weeks", "0 weeks", "+6 weeks")) +
  theme(
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  facet_grid(~"Matched Motor Units Within-Time and Unmatched Between-Time") -> plot_md_deltaf
plot_md_deltaf
#write.csv(df, file = "deltaf.csv", row.names = F)

#Effect size d
conditional_effect <- emmeans(fit2040, ~ time*intensity, adjust = "bonferroni")
eff_size(conditional_effect, sigma = sigma(fit2040), edf = df.residual(fit2040))

####Repeated measures correlation
#load rmcorr data
rmcorrdata <- read_excel("rmcorr_final.xlsx")
View(rmcorrdata)

#Delfa F vs Peak discharge rates - Unmatched motor units
rmc.dr <- rmcorr(as.factor(participant),
                  deltaf_untrack, dr_untrack, rmcorrdata, CIs = 
                    c("analytic", "bootstrap"), nreps = 1000,
                  bstrap.out = F)
rmc.dr

plot(rmc.dr, rmcorrdata, overall = T, 
     palette = NULL, xlab = "Delta F (pps)", ylab = "Peak Discharge Rate (pps)",
     overall.col = "gray60", overall.lwd =2, overall.lty = 3)

#Delfa F vs Peak discharge rates - Matched motor units
rmc.drt <- rmcorr(as.factor(mu_list),
                 deltafcp, peak_dr_cp, d20t, CIs = 
                   c("analytic", "bootstrap"), nreps = 1000,
                 bstrap.out = F)
rmc.drt

plot(rmc.drt, d20t, overall = T, 
     palette = NULL, xlab = "Delta F (pps)", ylab = "Peak Discharge Rate (pps)",
     overall.col = "gray60", overall.lwd =2, overall.lty = 3)

#Delfa F vs Peak torque - Unmatched motor units
rmc.pt <- rmcorr(as.factor(participant),
                  deltaf_untrack, pt, rmcorrdata, CIs = 
                    c("analytic", "bootstrap"), nreps = 1000,
                  bstrap.out = F)
rmc.pt

plot(rmc.pt, rmcorrdata, overall = T, 
     palette = NULL, xlab = "Delta F (pps)", ylab = "Peak Torque (N?m)",
     overall.col = "gray60", overall.lwd =2, overall.lty = 3)


#Delfa F vs TUG - Unmatched motor units
rmc.tug <- rmcorr(as.factor(participant),
                 deltaf_untrack, tug, rmcorrdata, CIs = 
                   c("analytic", "bootstrap"), nreps = 1000,
                 bstrap.out = F)
rmc.tug

rmc.1rm <- rmcorr(as.factor(participant),
                  deltaf_untrack, onerm, rmcorrdata, CIs = 
                    c("analytic", "bootstrap"), nreps = 1000,
                  bstrap.out = F)
rmc.1rm

rmc.1rm <- rmcorr(as.factor(participant),
                  deltaf_track, onerm, rmcorrdata, CIs = 
                    c("analytic", "bootstrap"), nreps = 1000,
                  bstrap.out = F)
rmc.1rm


rmc.five_sts <- rmcorr(as.factor(participant),
                  deltaf_untrack, five_sts, rmcorrdata, CIs = 
                    c("analytic", "bootstrap"), nreps = 100,
                  bstrap.out = F)
rmc.five_sts



rmc.thirty_sts <- rmcorr(as.factor(participant),
                       deltaf_untrack, thirty_sts, rmcorrdata, CIs = 
                         c("analytic", "bootstrap"), nreps = 100,
                       bstrap.out = F)
rmc.thirty_sts

rmc.cmj <- rmcorr(as.factor(participant),
                       deltaf_untrack, cmj_height, rmcorrdata, CIs = 
                         c("analytic", "bootstrap"), nreps = 100,
                       bstrap.out = F)
rmc.cmj


rmc.pt_t <- rmcorr(as.factor(participant),
                 deltaf_track, pt, rmcorrdata, CIs = 
                   c("analytic", "bootstrap"), nreps = 100,
                 bstrap.out = F)
rmc.pt_t

rmc.tug_t <- rmcorr(as.factor(participant),
                  deltaf_track, tug, rmcorrdata, CIs = 
                    c("analytic", "bootstrap"), nreps = 100,
                  bstrap.out = F)
rmc.tug_t

rmc.five_sts_t <- rmcorr(as.factor(participant),
                       deltaf_track, five_sts, rmcorrdata, CIs = 
                         c("analytic", "bootstrap"), nreps = 100,
                       bstrap.out = F)
rmc.five_sts_t

rmc.thirty_sts_t <- rmcorr(as.factor(participant),
                         deltaf_track, thirty_sts, rmcorrdata, CIs = 
                           c("analytic", "bootstrap"), nreps = 100,
                         bstrap.out = F)
rmc.thirty_sts_t

rmc.cmj_t <- rmcorr(as.factor(participant),
                  deltaf_track, cmj_height, rmcorrdata, CIs = 
                    c("analytic", "bootstrap"), nreps = 100,
                  bstrap.out = F)
rmc.cmj_t

#### Figures 3 and 4 layouts
#cowplot delta F 20% (unmatched + matched motor units)
deltaf_cowplot <- plot_grid(plot_delta_f, plot_deltaf_tcp, labels = "AUTO", label_size = 22, ncol = 2)
deltaf_cowplot

ggsave(file = "deltaf_cowplot.tiff", units="in", width = 11, height = 8.5, dpi = 600)

#cowplot delta F 20% vs 40%
deltaf2040_cowplot <- plot_grid(plot_delta_f2040, plot_md_deltaf, labels = "AUTO", rel_heights = c(2, 1), label_size = 22, ncol = 1)
deltaf2040_cowplot

ggsave(file = "deltaf2040_cowplot.tiff", units="in", width = 10, height = 12, dpi = 600)