setwd("~/Library/CloudStorage/OneDrive-epfl.ch/2025_isolates_alpine_2023_soils/ALANNA_EXPERIMENT/results")

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr) 
library(vegan) #needed for permanova
library(ggrepel)
library(wesanderson)  
library(viridis)


##### import META info 
r01_d01 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r01_d01") 
r01_d02 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r01_d02") 
r01_d03 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r01_d03") 
r01_d07 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r01_d07") 
r02_d01 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r02_d01") 
r02_d02 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r02_d02") 
r02_d03 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r02_d03") 
r02_d07 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r02_d07") 
r03_d01 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r03_d01")
r03_d02 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r03_d02")
r03_d03 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r03_d03")
r03_d07 <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "r03_d07")
meta <- bind_rows(r01_d01, r01_d02, r01_d03, r01_d07, 
                  r02_d01, r02_d02, r02_d03, r02_d07,
                  r03_d01, r03_d02, r03_d03, r03_d07
                  ) 
pairs <- readxl::read_xlsx(path = "META_pairwise.xlsx", sheet = "pair_meta")
META <- left_join(meta, pairs, by = "Community")


#### Importing Massimo's model of pairwise ####
mass_r1 <- read.csv("FC_data/coculture_pairwise_counts_round1.csv", header = TRUE)
mass_r2 <- read.csv("FC_data/coculture_pairwise_counts_round2.csv", header = TRUE)
mass_r3 <- read.csv("FC_data/coculture_pairwise_counts_round3.csv", header = TRUE)

mass_combined <- bind_rows(mass_r1, mass_r2, mass_r3
                           ) %>%
  rename(Isolate_A = isolate_A, 
         Isolate_B = isolate_B)%>% #make upper case to match with META files 
  mutate(Isolate_A = as.character(Isolate_A), 
         Isolate_B = as.character(Isolate_B))

#### adjust for dilutions
mass_model <- left_join(mass_combined, META, by = c("Round", "Temp", "Day", "Well", "Community", "Isolate_A", "Isolate_B")) ### NEED NEW DATA FROM MASS ! 

MASS_MODEL <-mass_model %>% 
  mutate(count_A = count_A*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>% #calculate events per mL in original sample
  mutate(count_B = count_B*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>%
  mutate(summed_proba_A = summed_proba_A*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>% #calculate events per mL in original sample
  mutate(summed_proba_B = summed_proba_B*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>%
  mutate(isolate_A_mono_counts = isolate_A_mono_counts*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>% #calculate events per mL in original sample
  mutate(isolate_B_mono_counts = isolate_B_mono_counts*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>%
  mutate(N = N*(Vol_Total/Vol)*(1/Fixed_dil)*1000)
# 
# MASS_MODEL2 <- MASS_MODEL %>% 
#   bind_rows(
#     MASS_MODEL %>%
#       distinct(Community, Temp, Round, IsolateA, IsolateB, Type, model_dir) %>%  # add all relevant cols
#       mutate(Day = d00, N = 2*10^6)     
#   )
# 
#
# 
# #put fake d00 = 10^6 cells/mL (what we aimed for at d00 inoculation)
# mass_long2 <- mass_long %>%
#   bind_rows(
#     mass_long %>%
#       distinct(Community, Temp, Round, IsolateA, IsolateB, Type) %>%   # one row per unique Community/Temp combo
#       mutate(Day = 0, N = 2*10^6)     
#   )

### filter out bad data
MASS_MODEL <- MASS_MODEL %>%
  filter(Isolate_A != 6, Isolate_B != 6) # remove 6 because didnt grow 


#### Transformation magic  to get correct layout... :) 
mass_long <- MASS_MODEL %>%
  pivot_longer(
    cols = c(count_A, count_B),
    names_to = "count_type",
    values_to = "count") %>%
  mutate(
    summed_prob = ifelse(count_type == "count_A", summed_proba_A, summed_proba_B),
    mono_counts = ifelse(count_type == "count_A", isolate_A_mono_counts, isolate_B_mono_counts)) %>%
  mutate(Isolate_B = ifelse(count_type == "count_A", NA, Isolate_B))%>%
  mutate(Isolate_A = ifelse(count_type == "count_B", NA, Isolate_A)) %>%
  mutate(Isolate = ifelse(count_type == "count_A", Isolate_A, Isolate_B)) %>%
select(-summed_proba_A, -summed_proba_B, -isolate_A_mono_counts, -isolate_B_mono_counts, -Isolate_A, -Isolate_B) %>%
  filter()
# unique(mass_long$Isolate)

#### PLOT MONO CULTURE CURVES ####
mass_long %>%
  bind_rows( ## add "fake" day 00 = 10^6
    mass_long %>%
      distinct(Community, Temp, model_dir, Isolate) %>%   # one row per unique group
      mutate(Day = "0", mono_counts = 10^6)
  ) %>%
  mutate(Day = readr::parse_number(as.character(Day))) %>%   # "day03" -> 3 ##makes numeric
  ggplot(aes(x = Day, y = log10(mono_counts))) +
  scale_colour_manual(
    values = c("10C" = "skyblue3","20C" = "orange2")) + 
  scale_fill_manual(
    values = c("10C" = "skyblue3","20C" = "orange2")) + 
  stat_summary(                        # Shaded 95% CI ribbon
    aes(group = Temp, fill= Temp),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    colour = NA
  ) +
  stat_summary(                        # Mean lines
    aes(group = Temp, colour = Temp),  # colour, not fill
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_x_continuous(breaks = c(1, 2, 3, 7)) +
  facet_wrap(~ Isolate, ncol = 4, scales = "free_x")+
  # coord_cartesian(ylim = c(6, 12)) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/monos_FC_avg.pdf", plot = last_plot(), height = 8, width = 12)


### see each replicate round
mass_long %>%
  bind_rows( ## add "fake" day 00 = 10^6
    mass_long %>%
      distinct(Community, Temp, model_dir, Isolate, Round) %>%   # one row per unique group
      mutate(Day = "0", mono_counts = 10^6)
  ) %>%
  mutate(Day = readr::parse_number(as.character(Day))) %>%   # "day03" -> 3 ##makes numeric
  ggplot(aes(x = Day, y = log10(mono_counts))) +
  scale_colour_manual(
    values = c("10C" = "skyblue3","20C" = "orange2")) + 
  scale_fill_manual(
    values = c("10C" = "skyblue3","20C" = "orange2")) + 
  stat_summary(                        # Shaded 95% CI ribbon
    aes(group = interaction(Temp, Round), fill= Temp),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    colour = NA
  ) +
  stat_summary(                        # Mean lines
    aes(group = interaction(Temp, Round), colour = Temp, linetype = Round),  # colour, not fill
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_x_continuous(breaks = c(1, 2, 3, 7)) +
  scale_linetype_manual(values = c("round1" = "solid", "round2" = "dashed", "round3" = "dotted")) +  # adjust level names to match yours
  facet_wrap(~ Isolate, ncol = 4, scales = "free_x")+
  # coord_cartesian(ylim = c(6, 12)) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/monos_FC_3rounds.pdf", plot = last_plot(), height = 8, width = 12)




######## Testing summed prob. vs. binary counts ####
mass_long %>% 
  ggplot(aes(x = log(summed_prob), y = log(count), 
             # colour = Round, 
             colour = model_dir, 
             shape = Day, alpha = Temp)) +
  geom_point(size = 3) +
  # coord_cartesian(xlim = c(0, 30), ylim =) +
  scale_alpha_manual(values = c("10C" = 0.5, "20C" = 1)) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/count_vs_summed_prob.pdf", plot = last_plot(), height = 8, width = 12)


# model_A <- lm(log(count) ~ log(summed_prob), data = mass_long)
# summary(model_A)


#### plot each 3 replicate separate #####
# mass_long %>%
#   # filter(Temp == "10C") %>%
#   # filter(Temp == "20C") %>%
#   # ggplot(aes(group = interaction(model_dir, Day, Temp), x = model_dir, y = log10(count), fill = isolate)) +
#   ggplot(aes(group = interaction(model_dir, Day, Temp), x = model_dir, y = count, fill = isolate)) +
#   # geom_col(position = "dodge", color = "white", linewidth = 0.3) +
#   geom_col(position = "stack", color = "white", linewidth = 0.3) +
#   scale_fill_viridis_d(option = "turbo") +   # try also "viridis", "plasma", "magma", "turbo
#   scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
#   facet_grid(Temp + Day ~ model_dir, scales = "free_x") +
#   theme_bw(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )
# ggsave(filename = "figures/r1_model_sep_reps.pdf", plot = last_plot(), height = 8, width = 10)
# 
#   
#### plot averages ####
# mass_avg <- mass_long %>%
#   # group_by(model_dir, Temp, Day, isolate, count_type) %>%
#   group_by(model_dir, Round, Temp, Day, isolate, count_type) %>%
#   summarise(count = mean(count), .groups = "drop")
# 
# mass_avg %>% 
#   # filter(Temp == "10C") %>%
#   # filter(Temp == "20C") %>%
#   filter(Day %in% c("d03", "d07")) %>%
#   ggplot(aes(x = model_dir, y = log10(count), fill = isolate, alpha = Round)) +
#   # ggplot(aes(x = model_dir, y = count, fill = isolate)) +
#   # geom_col(position = "stack", color = "white", linewidth = 0.3)  +
#   geom_col(position = "dodge", color = "white", linewidth = 0.3)  +
#   scale_fill_viridis_d(option = "turbo") +   # try also "viridis", "plasma", "magma", "turbo
#   scale_alpha_manual(values = c(0.6, 1)) +  
#   # scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
#   facet_grid(Temp + Day ~ model_dir, scales = "free_x") +
#   theme_bw(base_size = 12) +
#   theme(
#     axis.text.x = element_blank(),
#     strip.text.x = element_text(angle = 90, hjust = 0, size = 12),
#     strip.text.y = element_text(angle = 0, hjust = 0, size = 12)
#   )
# ggsave(filename = "figures/r1_model_10C.pdf", plot = last_plot(), height = 8, width = 16)
# ggsave(filename = "figures/r1_model_20C.pdf", plot = last_plot(), height = 8, width = 16)
# ggsave(filename = "figures/r1_model_both_stacked.pdf", plot = last_plot(), height = 8, width = 16)
# ggsave(filename = "figures/r1_model_both.pdf", plot = last_plot(), height = 8, width = 16)
# ggsave(filename = "figures/avg_DAY3_7.pdf", plot = last_plot(), height = 8, width = 16)

# 

### relative abundances 
# mass_avg %>%
#   # filter(Temp == "10C") %>%
#   # filter(Temp == "20C") %>%
#   filter(Day %in% c("d03", "d07")) %>%
#   group_by(model_dir, Temp, Day, Round) %>%
#   mutate(prop = count / sum(count)) %>%
#   ungroup() %>% 
#   ggplot(aes(x = interaction(model_dir, Round), y = prop, fill = isolate, alpha = Round)) +
#   geom_col(position = "stack", color = "white", linewidth = 0.3) +
#   scale_fill_viridis_d(option = "turbo") +   # try also "viridis", "plasma", "magma", "turbo
#   scale_alpha_manual(values = c(0.6, 1)) +  
#   scale_y_continuous(labels = scales::percent_format()) +
#   facet_grid(Temp + Day ~ model_dir, scales = "free_x") +
#   labs(y = "Relative abundance") +
#   theme_bw(base_size = 12) +
# theme(
#   axis.text.x = element_blank(),
#   strip.text.x = element_text(angle = 90, hjust = 0, size = 12),
#   strip.text.y = element_text(angle = 0, hjust = 0, size = 12)
# )
# # ggsave(filename = "figures/r1_model_rel_ab.pdf", plot = last_plot(), height = 8, width = 16)
# ggsave(filename = "figures/rel_abu_DAY3_7.pdf", plot = last_plot(), height = 8, width = 16)



#### abundances overtime ######
## want to add "fake" d00 = 2*10^6 cells

mass_long %>%
  # bind_rows( ## add "fake" day 00 = 10^6
  #   mass_long %>%
  #     distinct(Temp, model_dir) %>%   # one row per unique group
  #     mutate(Day = "0", count = 10^6)
  #   ) %>% 
  # mutate(Day = 0, Events = 2*10^6)  %>%
  mutate(Day = readr::parse_number(as.character(Day))) %>%   # "day03" -> 3 ##makes numeric
  # ggplot(aes(x = Day, y = log10(count))) +
  ggplot(aes(x = Day, y = log10(summed_prob))) +
  scale_fill_viridis_d(option = "turbo") +   # try also "viridis", "plasma", "magma", "turbo
  scale_colour_viridis_d(option = "turbo") +   # try also "viridis", "plasma", "magma", "turbo
  # ggplot(aes(x = Day, y = count)) +
  stat_summary(                        # Shaded 95% CI ribbon
    aes(group = Isolate, fill = Isolate),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    colour = NA
  ) +
  stat_summary(                        # Mean lines
    aes(group = Isolate, colour = Isolate),  # colour, not fill
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_x_continuous(breaks = c(1, 2, 3, 7)) +
  # facet_grid(Temp ~ model_dir, scales = "free_x") +
  facet_wrap(~ model_dir + Temp, ncol = 10, scales = "free_x")+
  coord_cartesian(ylim = c(6, 12)) +
  theme_bw(base_size = 16) 
ggsave(filename = "figures/A_B_time_avg.pdf", plot = last_plot(), height = 12, width = 20)



### see each round separate
mass_long %>%
  # filter(Round == "round1") %>%
  filter(Temp == "10C") %>%
  mutate(Day = readr::parse_number(as.character(Day))) %>%
  # ggplot(aes(x = Day, y = log10(count))) +
  ggplot(aes(x = Day, y = log10(summed_prob))) +
  scale_fill_viridis_d(option = "turbo") +
  scale_colour_viridis_d(option = "turbo") +
  stat_summary(
    aes(group = interaction(Isolate, Round), fill = Isolate),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    colour = NA
  ) +
  stat_summary(
    aes(group = interaction(Isolate, Round), colour = Isolate, linetype = Round),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_linetype_manual(values = c("round1" = "solid", "round2" = "dashed", "round3" = "dotted")) +  # adjust level names to match yours
  scale_x_continuous(breaks = c(1, 2, 3, 7)) +
  coord_cartesian(ylim = c(6, 12)) +
  # facet_grid(Temp ~ model_dir, scales = "free_x") +
  # facet_wrap(~ model_dir + Temp, ncol = 10, scales = "free_x")+
  facet_wrap(~ model_dir + Temp, ncol = 5, scales = "free_x")+
  
  theme_bw(base_size = 16) +
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 0, size = 12),
    # strip.text.x = element_text(angle = 90, hjust = 0, size = 12),
    # strip.text.y = element_text(angle = 0, hjust = 0, size = 12)
  )
# ggsave(filename = "figures/A_B_time_3rounds.pdf", plot = last_plot(), height = 12, width = 20)
ggsave(filename = "figures/A_B_time_3rounds_10C.pdf", plot = last_plot(), height = 12, width = 20)

# ggsave(filename = "figures/A_B_time_round1.pdf", plot = last_plot(), height = 12, width = 20)


#### total N. Use Mass' data 
## TOTAL N , averages across all 
mass_long %>%
  filter(count_type == "count_A") %>% ### need to do this because otherwise CI for N is skewed b/c double data
  bind_rows( ## add "fake" day 00 = 10^6
    mass_long %>%
      distinct(Community, Temp, model_dir) %>%   # one row per unique group
      mutate(Day = "0", N = 2*10^6)
  ) %>%
  mutate(Day = readr::parse_number(as.character(Day))) %>%   # "day03" -> 3 ##makes numeric
  ggplot(aes(x = Day, y = log10(N))) +
  stat_summary(                        # Shaded 95% CI ribbon
    aes(group = Temp, fill = Temp),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    colour = NA
  ) +
  stat_summary(                        # Mean lines
    aes(group = Temp, colour = Temp),  # colour, not fill
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_colour_manual(
    values = c("10C" = "skyblue3","20C" = "orange2")) + 
  scale_fill_manual(
    values = c("10C" = "skyblue3","20C" = "orange2")) + 
  scale_x_continuous(breaks = c(1, 2, 3, 7)) +
  # facet_grid(Temp ~ model_dir, scales = "free_x") +
  facet_wrap(~ model_dir, nrow = 4, scales = "free_x")+
  coord_cartesian(ylim = c(6, 12)) +
  # labs(y = "log10(total N events/mL)") +
  theme_bw(base_size = 16) 
ggsave(filename = "figures/totalN_avg.pdf", plot = last_plot(), height = 12, width = 20)



## TOTAL N , per round
mass_long %>%
  filter(count_type == "count_A") %>% ### need to do this because otherwise CI for N is scewed b/c double data
  bind_rows( ## add "fake" day 00 = 10^6
    mass_long %>%
      distinct(Temp, model_dir, Round) %>%   # one row per unique group
      mutate(Day = "0", N = 2*10^6) 
  ) %>%
  mutate(Day = readr::parse_number(as.character(Day))) %>%   # "day03" -> 3 ##makes numeric
  ggplot(aes(x = Day, y = log10(N))) +
  stat_summary(  # Shaded 95% CI ribbon
    aes(group = interaction(Temp, Round), fill = Temp),
    fun.data = mean_cl_normal, #bootstrap CI (no normality assumption)
    geom = "ribbon",
    alpha = 0.2,
    colour = NA
  ) +
  stat_summary(  # Mean lines
    aes(group = interaction(Temp, Round), colour = Temp, linetype = Round),  # colour, not fill
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_colour_manual(
    values = c("10C" = "skyblue3","20C" = "orange2")) + 
  scale_fill_manual(
    values = c("10C" = "skyblue3","20C" = "orange2")) + 
  scale_x_continuous(breaks = c(1, 2, 3, 7)) +
  # facet_grid(Temp ~ model_dir, scales = "free_x") +
  facet_wrap(~ model_dir, nrow = 4, scales = "free_x")+
  coord_cartesian(ylim = c(6, 12)) +
  # labs(y = "log10(total N events/mL)") +
  theme_bw(base_size = 16) 
ggsave(filename = "figures/totalN_per_round.pdf", plot = last_plot(), height = 12, width = 20)




### Find which have alot of noise (i.e., alot of events outside P1 gate) #####
# noisy <- COUNTS %>%  #majority in day1 10 C (no growth yet), and day7 20C (cells breaking/consumed? )
#   mutate(P1_ratio = P1.Count/All.Count) %>%
#   filter(P1_ratio < 0.7) #141 samples total. 96 = day1 10C, 8 = day2 10C, 7 = day 3 20C, 6 = day 7 10 C, 24 = day 7 20C


# ggsave(filename = "raw cells per mL.png", plot = last_plot(), height = 8, width = 10)

# 
# # import calculated R and K values from growth curves off plate reader (done in "pilot_growth_curves_v3.R")
# curve_values <- read.csv("curve_values.csv", header = TRUE)
# curve_values <- curve_values %>%
#   rename(Temp = Temperature) %>%#make match to COUNTS df
#   mutate(Isolate = sprintf("%02d", Isolate)) #puts a 0 in front of 06
# 
# ## look at just neg 
# negative <- COUNTS %>% filter(Community =="neg")




##### trying INTERACTIONS each time point ####

### CALCULATION (assume resource competition)

### for interactions
INTERACT <- MASS_MODEL %>%  # MASS_MODEL already in wide format
  mutate(Interaction = N/ ((isolate_A_mono_counts + isolate_B_mono_counts)/2)) %>% #if assuming resource competition (yes)
  mutate(Sign = case_when(
    Interaction > 1.1  ~ "positive",
    Interaction < 0.9  ~ "negative",
    TRUE         ~ "neutral"))
INTERACT$Sign <- factor(INTERACT$Sign, levels = c("positive", "neutral", "negative")) # just keeps things in order

# ## calculation2
# INTERACT <- interact %>% 
#   mutate(Interaction = (Events_pair - Events_monoA - Events_monoB) / (Events_monoA + Events_monoB))

# ## calculation3
# INTERACT <- interact %>% 
#   mutate(Interaction = Events_pair / (Events_monoA + Events_monoB))

# 
# ## Kamrad et al. 2026 used >10% as determining neg or pos 
# pos <- INTERACT %>% filter(Interaction >= 1.1)          #   335 interactions
# neg <- INTERACT %>% filter(Interaction <= 0.9)         #  770 interactions
# neu <- INTERACT %>% filter(Interaction > 0.9 & Interaction < 1.1)  # 71 interactions


### overall just 10 vs 20 C
INTERACT %>%
  # filter(Day %in% c("d03", "d07")) %>%
  # ggplot(aes(x = Temp, y = log(Interaction), colour = Temp, fill = Temp, shape = Sign)) +
  ggplot(aes(x = Temp, y = log(Interaction), colour = Sign, fill = Sign)) +
  # ggplot(aes(x = Temp, y = Interaction, colour = Temp, fill = Temp)) +
  # geom_violin(alpha = 0.5) + 
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.5) +
  # scale_colour_manual(
  #   values = c("10C" = "skyblue3","20C" = "orange2"
  #   )) +
  # scale_fill_manual(
  #   values = c("10C" = "skyblue3","20C" = "orange2"
  #   )) +
  
  scale_colour_manual(
    values = c("negative" = "orchid1","neutral" = "gray50","positive" = "seagreen1"
    )) +
  scale_fill_manual(
    values = c("negative" = "orchid1","neutral" = "gray50","positive" = "seagreen1"
    )) +
  
  # geom_hline(yintercept = log(1.1), color = "green", linetype  = "dashed", linewidth = 1) +
  # geom_hline(yintercept = log(0.9), color = "red", linetype  = "dashed", linewidth = 1) +
  facet_wrap(~Day) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/interaction_day_temp_log.pdf", plot = last_plot(), height = 8, width = 10)
# ggsave(filename = "interaction_day_comm_out.png", plot = last_plot(), height = 8, width = 10)


INTERACT %>%
  count(Day, Round, Temp, Sign) %>%
  group_by(Day, Round, Temp) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = Day, y = pct, fill = Sign)) +
  geom_col() +
  # scale_colour_manual(
  #   values = c("negative" = "orchid1","neutral" = "gray50","positive" = "seagreen1"
  #   )) +
  scale_fill_manual(
    values = c("negative" = "orchid1","neutral" = "gray50","positive" = "seagreen1"
    )) +
  facet_grid(Round ~ Temp) +
  labs(x = "Day", y = "Proportion (%)", fill = "Interaction") +
  theme_bw(base_size = 16)
# ggsave(filename = "figures/interaction_Sign.pdf", plot = last_plot(), height = 8, width = 10)
ggsave(filename = "figures/interaction_Sign_perROUND.pdf", plot = last_plot(), height = 8, width = 10)


#### extra plots playing #####

# # separate by day
# INTERACT %>% 
#   # filter(Interaction > -2000)  %>%
#   ggplot(aes(x = Community, y = Interaction, colour = Temp)) +
#   geom_jitter(width = 0.2, size = 2.5, alpha = 0.5) + 
#   scale_colour_manual(
#     values = c("10C" = "skyblue3","20C" = "orange2"
#   )) +
#   facet_wrap(~Day) +
#   theme_bw()
# ggsave(filename = "figures/interaction_day_comm.pdf", plot = last_plot(), height = 8, width = 10)
# # ggsave(filename = "interaction_day_comm_out.png", plot = last_plot(), height = 8, width = 10)
# 
# 
# INTERACT %>% 
#   # filter(Interaction > -2000)  %>%
#   ggplot(aes(x = Community, y = log(Interaction), colour = Temp)) +
#   geom_jitter(width = 0.2, size = 2.5, alpha = 0.5) + 
#   scale_colour_manual(
#     values = c("10C" = "skyblue3","20C" = "orange2"
#     )) +
#   facet_wrap(~Day) +
#   theme_bw()
# ggsave(filename = "figures/interaction_day_comm_log.pdf", plot = last_plot(), height = 8, width = 10)
# # ggsave(filename = "interaction_day_comm_out.png", plot = last_plot(), height = 8, width = 10)


### separate by Community 
# INTERACT %>% 
#   ggplot(aes(x = Day, y = Interaction, colour = Temp)) +
# geom_jitter(width = 0.2) + 
# scale_colour_manual(
#     values = c("10C" = "skyblue3","20C" = "orange2"
#     )) +
#   facet_wrap(~Community) +
#   theme_bw()


### trying anova  ####
### checking assumptions https://statsandr.com/blog/anova-in-r/

res_aov <- aov(log(Interaction) ~ Temp * Day + Round, data = INTERACT)
summary(res_aov)

# Marginal means for each effect
library(emmeans)
emmeans(res_aov, ~ Temp)
emmeans(res_aov, ~ Day)
emmeans(res_aov, ~ Round)
emmeans(res_aov, ~ Temp * Day, type = "response") 

### testing normality, etc. 
hist(res_aov$residuals) # NOT normally distributed...even when remove big outliers 
shapiro.test(res_aov$residuals) 

# QQ-plot
# qqPlot(res_aov$residuals,
#        id = FALSE # id = FALSE to remove point identification
# )
 
# Kruskall Wallis test for non-normal data https://statsandr.com/blog/kruskal-wallis-test-nonparametric-version-anova/
kruskal.test(log(Interaction) ~ Temp , data = INTERACT)
kruskal.test(log(Interaction) ~ Day , data = INTERACT)
kruskal.test(log(Interaction) ~ Round , data = INTERACT)

library(ggstatsplot)

INTERACT %>%
  mutate(log_Interaction = log(Interaction)) %>%
  ggbetweenstats(
    x = Temp,
    # y = Interaction,
    y = log_Interaction,
    type = "nonparametric",
    plot.type = "box",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    bf.message = FALSE
  ) +
  scale_colour_manual(
    values = c("10C" = "skyblue3","20C" = "orange2"
    )) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/KW_inter_temp.pdf", plot = last_plot(), height = 8, width = 10)



INTERACT %>%
  mutate(log_Interaction = log(Interaction)) %>%
  ggbetweenstats(
    x = Day,
    y = log_Interaction,
    type = "nonparametric",
    plot.type = "box",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    bf.message = FALSE
  ) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/KW_inter_Day.pdf", plot = last_plot(), height = 8, width = 10)


# +
#   scale_colour_manual(
#     values = c("10C" = "skyblue3","20C" = "orange2"
#     ))


INTERACT %>%
  mutate(log_Interaction = log(Interaction)) %>%
  ggbetweenstats(
    x = Round,
    y = log_Interaction,
    type = "nonparametric",
    plot.type = "box",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    bf.message = FALSE
  ) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/KW_inter_Round.pdf", plot = last_plot(), height = 8, width = 10)





##### Standardise to monos ######
### nope, this also doesnt make sense.... 
mass_long <- mass_long %>%
  mutate(mono_stand = summed_prob/mono_counts) 

mass_long %>%
  # filter(Round == "round1") %>%
  filter(Temp == "10C") %>%
  mutate(Day = readr::parse_number(as.character(Day))) %>%
  ggplot(aes(x = Day, y = log10(mono_stand))) +
  scale_fill_viridis_d(option = "turbo") +
  scale_colour_viridis_d(option = "turbo") +
  stat_summary(
    aes(group = interaction(Isolate, Round), fill = Isolate),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    colour = NA
  ) +
  stat_summary(
    aes(group = interaction(Isolate, Round), colour = Isolate, linetype = Round),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_linetype_manual(values = c("round1" = "solid", "round2" = "dashed", "round3" = "dotted")) +  # adjust level names to match yours
  scale_x_continuous(breaks = c(1, 2, 3, 7)) +
  coord_cartesian(ylim = c(-2.5, 2.5)) +
  # facet_grid(Temp ~ model_dir, scales = "free_x") +
  facet_wrap(~ model_dir + Temp, ncol = 6, scales = "free_x")+
  theme_bw(base_size = 16) +
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 0, size = 12),
    # strip.text.x = element_text(angle = 90, hjust = 0, size = 12),
    # strip.text.y = element_text(angle = 0, hjust = 0, size = 12)
  )
ggsave(filename = "figures/mono_standardised.pdf", plot = last_plot(), height = 12, width = 20)



# ## looking at differences between A and B.... not working & doesn't make sense ######
# MASS_MODEL <- MASS_MODEL %>%
#   mutate(AB_diff = summed_proba_A - summed_proba_B)
# 
# MASS_MODEL %>%
#   # filter(Round == "round1") %>%
#   filter(Temp == "10C") %>%
#   mutate(Day = readr::parse_number(as.character(Day))) %>%
#   ggplot(aes(x = Day, y = log10(AB_diff))) +
#   # scale_fill_viridis_d(option = "turbo") +
#   # scale_colour_viridis_d(option = "turbo") +
#   stat_summary(
#     aes(group = Round),
#     fun.data = mean_cl_normal,
#     geom = "ribbon",
#     alpha = 0.2,
#     colour = NA
#   ) +
#   stat_summary(
#     aes(group = Round, linetype = Round),
#     fun = mean,
#     geom = "line",
#     linewidth = 1.2
#   ) +
#   scale_linetype_manual(values = c("round1" = "solid", "round2" = "dashed", "round3" = "dotted")) +  # adjust level names to match yours
#   scale_x_continuous(breaks = c(1, 2, 3, 7)) +
#   geom_hline(yintercept = log10(0.00000000001), colour = "blue")+
#   # coord_cartesian(ylim = c(-2.5, 2.5)) +
#   # facet_grid(Temp ~ model_dir, scales = "free_x") +
#   facet_wrap(~ model_dir + Temp, ncol = 6, scales = "free_x")+
#   theme_bw(base_size = 16) +
#   theme(
#     # axis.text.x = element_text(angle = 90, hjust = 0, size = 12),
#     # strip.text.x = element_text(angle = 90, hjust = 0, size = 12),
#     # strip.text.y = element_text(angle = 0, hjust = 0, size = 12)
#   )




#### testing if correlation between later outcome to d01 ####
library(ggpmisc) ### for adding p-values

MASS_MODEL <- MASS_MODEL %>%
  mutate(AB_ratio = summed_proba_A / summed_proba_B)

# need to set up for a regression 
MASS_wide <- MASS_MODEL %>%
  select(Community, model_dir, Round, Temp, Well, Day, AB_ratio) %>%
  pivot_wider(names_from = Day,
              values_from = AB_ratio,
              names_prefix = "day_")

## day 3 ~ day 1
lm_d1_d3 <- lm(day_d03 ~ day_d01 +Community + Round, data = MASS_wide) # want to control for each community , not sure if this is right 
# lm_d1_d3 <- lm(day_d03 ~ day_d01, data = MASS_wide) 
summary(lm_d1_d3) # day 1 significant... 

MASS_wide %>% 
  ggplot(aes(log(day_d01), log (day_d03), colour = Community, shape = Round)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "black", linewidth = 1,
              inherit.aes = FALSE,    # want just one big grouping mes                      
              aes(x = log(day_d01), y = log(day_d03))) +
  theme_bw(base_size = 16) +
  labs(x = "log(AB_ratio) Day 1",
       y = "log(AB_ratio) Day 3")
ggsave(filename = "figures/AB_ratio_d01vsd03.pdf", plot = last_plot(), height = 8, width = 10)


## day 7 ~ day 1
lm_d1_d7<- lm(day_d07 ~ day_d01 +Community + Round, data = MASS_wide)
# lm_d1_d7<- lm(day_d07 ~ day_d01, data = MASS_wide)
summary(lm_d1_d7) # day 1 NOT significant... 

MASS_wide %>% 
  ggplot(aes(log(day_d01), log (day_d07), colour = Community, shape = Round)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "black", linewidth = 1,
              inherit.aes = FALSE, # want just one big grouping 
              aes(x = log(day_d01), y = log(day_d07))) +
# +
#   stat_poly_eq(inherit.aes = FALSE,
#                aes(x = log(day_d01), y = log(day_d07),
#                    label = paste(after_stat(rr.label), 
#                                  after_stat(p.value.label), 
#                                  sep = "~~")),
#                label.x = "left", label.y = "top",
#                size = 5) +
  theme_bw(base_size = 16) +
  labs(x = "log(AB_ratio) Day 1",
       y = "log(AB_ratio) Day 7")
ggsave(filename = "figures/AB_ratio_d01vsd07.pdf", plot = last_plot(), height = 8, width = 10)
  



