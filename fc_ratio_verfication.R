### FLOW CYTOMETER RATIO VERIFICATION
## Alanna combined 4 isolates in ratios in all combinations and ratios of 0:100, 25:75, 50:50, 75:25, 100:0
## these were run in technical replicates of 3 on Flow Cytometer
## Massimo then applied his classification model 

setwd("~/Library/CloudStorage/OneDrive-epfl.ch/2025_isolates_alpine_2023_soils/ALANNA_EXPERIMENT/results/ratio_test_FC_classifier")

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr) 
library(vegan) #needed for permanova
library(ggrepel)
# library(wesanderson)  
library(viridis)
library(stringr)
library(growthcurver) ### N
library(lme4)
library(lmerTest) 
library(car)
library(ltc) # for pretty palettes :) 


#### setting up META ####
META <- readxl::read_xlsx(path = "FC_ratio_test_info.xlsx", sheet = "META") 

META <- META  %>% 
  filter(!Type %in% c("mono", "control"))  %>% 
  dplyr::select(-Type, -Round, -Day, -Fixed_dil, -Vol, -Vol_Total) %>%  #remove unneeded columns
  rename(Expected_A = Ratio_A, Expected_B = Ratio_B) %>% 
  mutate( Expected_A = as.numeric(Expected_A),
          Expected_B = as.numeric(Expected_B)
  )
# str(META)
META$model_dir <- paste(META$Isolate_A, "vs", META$Isolate_B, sep = "_")



### second test ####
#second test makes 6 models (for each model_dir), then applies to each replicate
df2_raw <- read.csv("model_validation_results.csv", header = TRUE) 
df2_raw <- df2_raw %>% dplyr::select(-count_A, -count_B, -N_events, -filename) #remove unneeded columns

df2 <- merge(df2_raw, META, by = "Well")

df2 <- df2 %>% 
  mutate(rel_ab_A = summed_proba_A / (summed_proba_A + summed_proba_B),
         rel_ab_B = summed_proba_B / (summed_proba_A + summed_proba_B),
         perc_diff = rel_ab_A - Expected_A) 

#### personal colour palettes####
# (https://github.com/loukesio/ltc-color-palettes)
hat_palette <- ltc("hat", n = length(unique(df2$model_dir)), type = "continuous")
expevo_palette <- ltc("expevo", n = length(unique(df2$model_dir)), type = "continuous")

casa_palette <- ltc("casa_natal") 
casa_palette <- casa_palette[-c(2, 8,9)] #remove 3/9 because we have 7 isolates
scales::show_col(casa_palette)


## plot percent difference for each mix
df2 %>% 
  filter(Isolate_A != "43") %>%
  filter(Isolate_B != "43") %>%
  # ggplot(aes(x = Community.x, y = perc_diff, colour = model_dir)) +
  ggplot(aes(x = Expected_A, y = perc_diff, colour = model_dir)) +
  geom_jitter(alpha = 0.7, size = 5, width = 0, height = 0) +
  scale_colour_manual(values = hat_palette) +
  # scale_colour_manual(values = casa_palette) +
  
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75)) +
  facet_wrap(~model_dir) +
  # scale_colour_manual(values = expevo_palette) +
  # labs(title = "second output") +
  theme_bw()
# ggsave(filename = "perc_diff_2.pdf", plot = last_plot(), height = 8, width = 10)
ggsave(filename = "perc_diff_no43.pdf", plot = last_plot(), height = 8, width = 10)


## plot real A ratio vs. calculated count 
df2 %>% 
  filter(Isolate_A != "43") %>%
  filter(Isolate_B != "43") %>%
  ggplot(aes(x = Expected_A, y = rel_ab_A, colour = model_dir, fill = model_dir)) +
  geom_point(alpha = 0.7, size = 4) + 
  scale_colour_manual(values = hat_palette) +
  scale_fill_manual(values = hat_palette) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  scale_x_continuous(limits = c(0.2, 0.9)) +
  scale_y_continuous(limits = c(0.2, 0.9)) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75)) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75)) +
  # labs(title = "second output") +
  theme_bw()
# ggsave(filename = "realA_vs_relA_2.pdf", plot = last_plot(), height = 8, width = 10)
ggsave(filename = "realA_vs_relA_no43.pdf", plot = last_plot(), height = 8, width = 10)

### linear regression
library(MASS) # for robust linear models (i.e., <20 samples for a regression)

model <- rlm(rel_ab_A ~ Expected_A * model_dir, data = df2)
summary(model)


rlm_results2 <- df2 %>% # from Claude , ## uses car pkg
  group_by(model_dir) %>%
  group_modify(~ {
    fit <- rlm(rel_ab_A ~ Expected_A, data = .x)
    test <- linearHypothesis(fit, c(
      "(Intercept) = 0", 
      "Expected_A = 1")) 
    tibble(
      intercept = coef(fit)[1],
      slope = coef(fit)[2],
      p_value = test$`Pr(>F)`[2]
    )
  })

rlm_results2
write.csv(rlm_results2, "rlm_results2.csv", row.names = FALSE)

detach("package:MASS", unload = TRUE) # but causes masking issues with dyplr, so detach once not needed









# ### first test ####
# #first test takes account replicate into model
## we decided to use second approach 
# df1_raw <- read.csv("coculture_pairwise_counts_round4.csv", header = TRUE) 
# df1_raw <- df1_raw %>% dplyr::select(-Round, -Temp, -Day, -Date, -count_A, -count_B, -isolate_A_mono_counts, -isolate_B_mono_counts, -N, -filename, -model_dir) #remove unneeded columns
# 
# df1 <- merge(df1_raw, META, by = "Well")
# 
# df1 <- df1 %>% 
#   mutate(rel_ab_A = summed_proba_A / (summed_proba_A + summed_proba_B),
#          rel_ab_B = summed_proba_B / (summed_proba_A + summed_proba_B),
#          perc_diff = rel_ab_A - Expected_A) 
# 
# ## plot percent difference for each mix
# df1 %>% ggplot(aes(x = Community.x, y = perc_diff, colour = model_dir )) +
#   geom_jitter(alpha = 0.7, size = 5, width = 0.2, height = 0) +
#   scale_colour_manual(values = hat_palette) +
#   labs(title = "first output") +
#   theme_bw()
# ggsave(filename = "perc_diff_1.pdf", plot = last_plot(), height = 8, width = 10)
# 
# 
# ## plot real A ratio vs. calculated count 
# df1 %>% 
#   # filter(isolate_A != "a43") %>%
#   # filter(isolate_B != "a43") %>%
#   ggplot(aes(x = Expected_A, y = rel_ab_A, colour = model_dir, fill = model_dir)) +
#   geom_point(alpha = 0.7, size = 4) + 
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
#   geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(title = "first output") +
#   theme_bw()
# ggsave(filename = "realA_vs_relA_1.pdf", plot = last_plot(), height = 8, width = 10)
# 
# ### linear regression
# model <- rlm(rel_ab_A ~ Expected_A * model_dir, data = df1)
# summary(model)
# 
# library(MASS) # for robust linear models (i.e., <20 samples for a regression)

# rlm_results1 <- df1 %>% # from Claude , ## uses car pkg
#   group_by(model_dir) %>%
#   group_modify(~ {
#     fit <- rlm(rel_ab_A ~ Expected_A, data = .x)
#     test <- linearHypothesis(fit, c("(Intercept) = 0", "Expected_A = 1")) 
#     tibble(
#       intercept = coef(fit)[1],
#       slope = coef(fit)[2],
#       p_value = test$`Pr(>F)`[2]
#     )
#   })
# 
# rlm_results1
# write.csv(rlm_results1, "rlm_results1.csv", row.names = FALSE)

# detach("package:MASS", unload = TRUE) # but causes masking issues with dyplr, so detach once not needed

