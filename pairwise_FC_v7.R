### v7: 14 august 2026. 
# check / fix AUC of MONOS. something wrong 
# add notes to draft figure numbers (noted in ggsave name) and script notes elsewhere. 

setwd("~/Library/CloudStorage/OneDrive-epfl.ch/2025_isolates_alpine_2023_soils/ALANNA_EXPERIMENT/results")

### load pkgs
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr) 
library(vegan) #needed for permanova
library(ggrepel)
library(wesanderson)  
library(viridis)
library(stringr)
library(growthcurver) ### N
library(lme4)
library(lmerTest) 
# detach("package:MASS", unload = TRUE) # MASS can cause problems with other pkgs, so sometimes need to detatch


#### COLOUR PALETTES ####
#for very pretty colour palettes :) 
install.packages("ltc")
library(ltc)

# (https://github.com/loukesio/ltc-color-palettes) 
# hat_palette <- ltc("hat", n = 7, type = "continuous") # we have 7 isolates
# scales::show_col(hat_palette)

# hat_palette <- ltc("hat") 
# hat_palette <- hat_palette[-c(1, 9, 10)] #remove some to have 7 isolates
# scales::show_col(hat_palette)

expevo_palette <- ltc("expevo", n = 6)          # original 6 colours
expevo_palette <- c(expevo_palette, "#4DB27C")  # add my 7th colour
scales::show_col(expevo_palette)

# casa_palette <- ltc("casa_natal") 
# casa_palette <- casa_palette[-c(8,9)] #remove 2/9 because we have 7 isolates
# scales::show_col(casa_palette)


#### IMPORT DATA ####
## each df is META data for day*round of flow cytometry (FC) data
## "Community" should have been better named "Pair" or "Combo"
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

## Massimo's machine learning model output data of FC data (total counts and counts of each individual isolate)
mass_r1 <- read.csv("FC_data/coculture_pairwise_counts_round1.csv", header = TRUE)
mass_r2 <- read.csv("FC_data/coculture_pairwise_counts_round2.csv", header = TRUE)
mass_r3 <- read.csv("FC_data/coculture_pairwise_counts_round3.csv", header = TRUE)

mass_combined <- bind_rows(mass_r1, mass_r2, mass_r3) %>%
  rename(Isolate_A = isolate_A, 
         Isolate_B = isolate_B) %>% #make upper case to match with META files 
  mutate(Isolate_A = as.character(Isolate_A), 
         Isolate_B = as.character(Isolate_B))

## combine with META data
mass_model <- left_join(mass_combined, META, by = c("Round", "Temp", "Day", "Well", "Community", "Isolate_A", "Isolate_B"))  

## calculate events/mL based on dilutions & volumes used on FC plate 
MASS_MODEL <-mass_model %>% 
  mutate(count_A = count_A*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>% #calculate events per mL in original sample
  mutate(count_B = count_B*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>%
  mutate(summed_proba_A = summed_proba_A*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>% #calculate events per mL in original sample
  mutate(summed_proba_B = summed_proba_B*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>%
  mutate(isolate_A_mono_counts = isolate_A_mono_counts*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>% #calculate events per mL in original sample
  mutate(isolate_B_mono_counts = isolate_B_mono_counts*(Vol_Total/Vol)*(1/Fixed_dil)*1000) %>%
  mutate(N = N*(Vol_Total/Vol)*(1/Fixed_dil)*1000)

## filter out bad data
MASS_MODEL <- MASS_MODEL %>%
  filter(Isolate_A != 6, Isolate_B != 6) # remove 6 because didnt grow 

## Transformation/pivoting magic for correct layout :) 
MASS_LONG <- MASS_MODEL %>%
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

## add fake day 0 cells/mL, make numerical 
fake_d0 <- MASS_LONG %>%
  distinct(Community, Well, Round, Temp, model_dir, Isolate) %>%
  mutate(Day = "0", mono_counts = 2*1e6, summed_prob = 1e6, count = 1e6, N = 2*1e6)

mass_long <- bind_rows(MASS_LONG, fake_d0) %>%
  mutate(Day = as.numeric(str_remove(Day, "^d"))) %>%
  select(-c(Date, Sample, filename, Fixed_dil, Vol, Vol_Total, Type)) # remove extra META info no longer needed

## isolate numbers to also say genera
mass_long <- mass_long %>% # do we also want to make Italics always here?? requires R markdown & ggtext in ggplot.
  mutate(Isolate = case_match(
    Isolate,
    "12" ~ "12_Sphingomonas",
    "16" ~ "16_Pseudoarthrobacter",
    "20" ~ "20_Paenibacillus",
    "23" ~ "23_Variovorax",
    "32" ~ "32_Luteibacter",
    "43" ~ "43_Psychrobacillus",
    "56" ~ "56_Flavobacterium",
    .default = Isolate
  ))

#### test summed prob. vs. binary counts ####
## two ways to calculate number events of each isolate in pair. A) assigning 0 or 1, or B) 
# MASS_LONG %>% 
#   ggplot(aes(x = log(summed_prob), y = log(count), 
#              # colour = Round, 
#              colour = model_dir, 
#              shape = Day, alpha = Temp)) +
#   geom_point(size = 3) +
#   # coord_cartesian(xlim = c(0, 30), ylim =) +
#   scale_alpha_manual(values = c("10C" = 0.5, "20C" = 1)) +
#   theme_bw(base_size = 16)
# ggsave(filename = "figures/count_vs_summed_prob.pdf", plot = last_plot(), height = 8, width = 12)

## inspect for NAs (all fine)
# mass_long %>% filter(if_any(-count_type, is.na)) 






#### MONOS plotting ####
## MONOs average
mass_long %>%
  ggplot(aes(x = Day, y = log10(mono_counts))) +
  scale_colour_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07"),
    labels = c("10C" = "10°C", "20C" = "20°C"),
    name = "Temperature") +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07"),
    labels = c("10C" = "10°C", "20C" = "20°C"),
    name = "Temperature") +
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
  scale_x_continuous(breaks = c(0, 1, 2, 3, 7)) +
  facet_wrap(~ Isolate, ncol = 4, scales = "free_x")+
  labs(y = expression(log[10]*"(cells/mL)")) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/Fig3A_monos_mean.pdf", plot = last_plot(), height = 8, width = 12)



## MONOs replicate rounds
mass_long %>%
  ggplot(aes(x = Day, y = log10(mono_counts))) +
  scale_colour_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07"),
    labels = c("10C" = "10°C", "20C" = "20°C"),
    name = "Temperature") +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07"),
    labels = c("10C" = "10°C", "20C" = "20°C"),
    name = "Temperature") +
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
  scale_x_continuous(breaks = c(0, 1, 2, 3, 7)) +
  scale_linetype_manual(values = c("round1" = "solid", "round2" = "dashed", "round3" = "dotted")) +  # adjust level names to match yours
  facet_wrap(~ Isolate, ncol = 4, scales = "free_x")+
  labs(y = expression(log[10]*"(cells/mL)")) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/FigS2_monos_3rounds.pdf", plot = last_plot(), height = 8, width = 12)







### MONOS calculate r and K ##### 
MONOS <- mass_long %>%
  select(-c(count_type)) %>% # need to remove, otherwise d00 fake counts don't convert to wide format
  pivot_wider(
    id_cols      = c(Community, Well, Temp, Isolate, model_dir,  Round),
    names_from   = Day,
    values_from  = mono_counts, #only interested in mono_counts here
    names_prefix = "mono_d0"
  )
MONOS %>% filter(if_any(everything(), is.na)) ### PROBLEM... NAs for mono_d07. all NAs are from round1*day7. Isolate 32 (x3), 56 (x1), 16 (x1), 20 (x1), 12 (x1), 23 (x1)

## calculate a) maximum cells/mL reached, and b) growth rate between day 00 to day 03
MONOS <- MONOS %>%
  mutate(
    mono_max = pmax(mono_d01, mono_d02, mono_d03, mono_d07, na.rm = TRUE), # max density reached
    # mono_r_d07 = (mono_d00 + mono_d07) / 7, # growth rate to day 07
    mono_r_d03 = (mono_d00 + mono_d03) / 3) # growth rate to day 03

# collapse / summarise (i.e., mono_max, mono_r_d07, mono_r_d03 is same for each Isolate, so can group for Isolate, Temp, & Round)
MONOS_summary <- MONOS %>%
  select(-mono_d02, -mono_d07, -mono_d03, -mono_d01, -mono_d00) %>%
  group_by(Isolate, Temp, Round) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

## plot "K" 
MONOS_summary %>% 
  ggplot(aes(x = Isolate, y = mono_max, colour = Temp, fill = Temp)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(
    size = 5,
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15, jitter.height = 0)
  ) +
  scale_colour_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  theme_bw(base_size = 16)

## plot "r" up to day3
MONOS_summary %>% 
  ggplot(aes(x = Isolate, y = log(mono_r_d03), colour = Temp, fill = Temp)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(
    size = 5,
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15, jitter.height = 0)
  ) +
    scale_colour_manual(
      values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
      scale_fill_manual(
        values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  theme_bw(base_size = 16)


## PLOT RATE vs. MAX DENSITY ("r" vs "K")
## determine centroids of each 3 replicates 
centroids <- MONOS_summary %>%
  group_by(Isolate, Temp) %>%
  summarise(mono_r_d03 = mean(mono_r_d03, na.rm = TRUE),
            mono_max   = mean(mono_max, na.rm = TRUE),
            .groups = "drop")

MONOS_summary %>%
  ggplot(aes(x = log(mono_r_d03), y = log(mono_max), colour = Isolate, fill = Isolate, alpha = factor(Temp))) +
  geom_point(size = 4) +
  geom_point(data = centroids,
             aes(x = log(mono_r_d03), y = log(mono_max), colour = Isolate, alpha = factor(Temp)),
             shape = 17,      # fixed triangle shape
             size = 8,
             stroke = 2,
             show.legend = FALSE) +
  geom_line(data = centroids,
            aes(x = log(mono_r_d03), y = log(mono_max), colour = Isolate, group = Isolate),
            linewidth = 1,
            alpha = 1,
            inherit.aes = FALSE,
            show.legend = FALSE) +
  geom_smooth(data = MONOS_summary,
              aes(x = log(mono_r_d03), y = log(mono_max), group = Temp, linetype = Temp),
              method = "lm", inherit.aes = FALSE, color = "black", se = TRUE, alpha = 0.15) +
  scale_colour_manual(values = expevo_palette) +
  scale_alpha_manual(
    values = c("10C" = 0.6, "20C" = 1),
    labels = c("10C" = "10°C", "20C" = "20°C"),
    name = "Temperature") +
  scale_linetype_manual(
    values = c("10C" = "solid", "20C" = "dotted"),
    labels = c("10C" = "10°C", "20C" = "20°C"),
    name = "Temperature") +
  guides( # extends line in legend icon 
    alpha = guide_legend(keywidth = unit(2.2, "cm"), override.aes = list(linewidth = 1.2)),
    linetype = guide_legend(keywidth = unit(2.2, "cm"), override.aes = list(linewidth = 1.2))
  ) +
  theme_bw(base_size = 16)
ggsave(filename = "figures/Fig3B_monos_r_K_lines.pdf", plot = last_plot(), height = 8, width = 10)






#### MONOS r & K stats ####
## mean rate between d00 and d03
aov_mono_r <- aov(mono_r_d03 ~ Temp*Isolate, data = MONOS_summary)
summary(aov_mono_r)

tukey_mono_r <- TukeyHSD(aov_mono_r)
tukey_r_df <- as.data.frame(tukey_mono_r$`Temp:Isolate`)
write.csv(tukey_r_df, "tukey_temp_isolate_mono_r.csv", row.names = TRUE)

## max density 
aov_mono_max <- aov(mono_max ~ Temp*Isolate, data = MONOS_summary)
summary(aov_mono_max)

tukey_mono_max <- TukeyHSD(aov_mono_max)
tukey_max_df <- as.data.frame(tukey_mono_max$`Temp:Isolate`)
write.csv(tukey_max_df, "tukey_temp_isolate_mono_max.csv", row.names = TRUE)

## t.tests for Flavobacterium
t_test_df <- MONOS_summary %>% filter(Isolate == "56_Flavobacterium")
t.test(mono_r_d03 ~ Temp, data = t_test_df)
t.test(mono_max ~ Temp, data = t_test_df)

## test for overall growth rate trade-off 
# cor.test(MONOS_summary$mono_r_d03, MONOS_summary$mono_max, 
#          method = "pearson")
cor.test(MONOS_summary$mono_r_d03, MONOS_summary$mono_max, 
         method = "spearman")

mono_tradeoff_model <- lmer(mono_max ~ mono_r_d03 + Temp + (1 | Isolate), data = MONOS_summary)
summary(mono_tradeoff_model)







#### PAIRS plotting ######

## PAIRS average 3 rounds
mass_long %>%
  ggplot(aes(x = Day, y = log10(summed_prob))) +
  scale_colour_manual(values = expevo_palette) +
  scale_fill_manual(values = expevo_palette) +
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
  scale_x_continuous(breaks = c(0, 1, 2, 3, 7)) +
  facet_wrap(~ model_dir + Temp, ncol = 10, scales = "free_x")+
  coord_cartesian(ylim = c(6, 12)) +
  theme_bw(base_size = 16) 
ggsave(filename = "figures/Fig4_A_B_mean.pdf", plot = last_plot(), height = 12, width = 20)


## PAIRS 3 rounds separately 
mass_long %>%
  filter(Temp == "10C") %>%
  mutate(Day = readr::parse_number(as.character(Day))) %>%
  ggplot(aes(x = Day, y = log10(summed_prob))) +
  scale_colour_manual(values = expevo_palette) +
  scale_fill_manual(values = expevo_palette) +
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
  scale_linetype_manual(values = c("round1" = "solid", "round2" = "dashed", "round3" = "dotted")) +  
  scale_x_continuous(breaks = c(0, 1, 2, 3, 7)) +
  coord_cartesian(ylim = c(6, 12)) +
  facet_wrap(~ model_dir + Temp, ncol = 5, scales = "free_x")+
  theme_bw(base_size = 16)
ggsave(filename = "figures/pairs_A_B_3rounds_10C.pdf", plot = last_plot(), height = 12, width = 20)


## PAIRS - TOTAL N, average of 3 rounds
mass_long %>% 
  filter(count_type == "count_A" | is.na(count_type)) %>% ### need to do this because otherwise CI for N is skewed b/c double data for A, B (fake d00 is NAs)
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
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 7)) +
  # facet_grid(Temp ~ model_dir, scales = "free_x") +
  facet_wrap(~ model_dir, nrow = 4, scales = "free_x")+
  coord_cartesian(ylim = c(6, 12)) +
  # labs(y = "log10(total N events/mL)") +
  theme_bw(base_size = 16) 
ggsave(filename = "figures/pair_totalN_mean.pdf", plot = last_plot(), height = 12, width = 20)


## PAIRS - TOTAL N, rounds separate
mass_long %>%
  filter(count_type == "count_A" | is.na(count_type)) %>% ### need to do this because otherwise CI for N is scewed b/c double data
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
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  scale_x_continuous(breaks = c(1, 2, 3, 7)) +
  facet_wrap(~ model_dir, nrow = 4, scales = "free_x")+
  coord_cartesian(ylim = c(6, 12)) +
  theme_bw(base_size = 16) 
ggsave(filename = "figures/pair_totalN_per_round.pdf", plot = last_plot(), height = 12, width = 20)



#### PAIRS N & respective MONOS plotting ####
## see how total cell density of pairs compared to the respective monocultures

## some complicated to have correct colours/matching. 7 isolate monocultures & also N of the pair , so 8 unique colours total
isolate_levels <- sort(unique(mass_long$Isolate))
expevo_palette <- setNames(c(ltc("expevo", n = 6), "#4DB27C"), isolate_levels)  # 7 isolates, 7 colours # want to use same colour palette as other plots 
line_cols <- c(expevo_palette, "N" = "black")  # add N, total cell density only for colour scale, not fill scale

## plot for 10 C
mass_long %>%
  filter(Temp == "10C") %>%
  ggplot(aes(x = Day)) +
  stat_summary(   # Ribbon + line for mono_counts, split by Isolate
    aes(y = log10(mono_counts), group = Isolate, fill = Isolate),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    colour = NA
  ) +
  stat_summary(
    aes(y = log10(mono_counts), group = Isolate, colour = Isolate),
    fun = mean,
    geom = "line",
    linewidth = 1.2
    # solid is the default, so no linetype needed here
  ) +
  stat_summary(   # Ribbon + line for N (grey, single line per facet)
    aes(y = log10(N), group = 1, colour = "N"),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    fill = "grey40",
    colour = NA
  ) +
  stat_summary(
    aes(y = log10(N), group = 1, colour = "N"),
    fun = mean,
    geom = "line",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  scale_fill_manual(values = expevo_palette, name = "Isolate") +
  scale_colour_manual(values = line_cols, name = NULL) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 7)) +
  facet_wrap(~ model_dir, ncol = 7, scales = "free_x") +
  coord_cartesian(ylim = c(6, 12)) +
  labs(y = expression(log[10]*"(cells/mL)"),
       title = "A. 10°C") +
  theme_bw(base_size = 16)
ggsave(filename = "figures/FigS4Amono_counts_N_10C.pdf", plot = last_plot(), height = 12, width = 20)


## plot for 20 C
mass_long %>%
  filter(Temp == "20C") %>%
  ggplot(aes(x = Day)) +
  stat_summary(   # Ribbon + line for mono_counts, split by Isolate
    aes(y = log10(mono_counts), group = Isolate, fill = Isolate),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    colour = NA
  ) +
  stat_summary(
    aes(y = log10(mono_counts), group = Isolate, colour = Isolate),
    fun = mean,
    geom = "line",
    linewidth = 1.2
    # solid is the default, so no linetype needed here
  ) +
  stat_summary(   # Ribbon + line for N (grey, single line per facet)
    aes(y = log10(N), group = 1, colour = "N"),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.2,
    fill = "grey40",
    colour = NA
  ) +
  stat_summary(
    aes(y = log10(N), group = 1, colour = "N"),
    fun = mean,
    geom = "line",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  scale_fill_manual(values = expevo_palette, name = "Isolate") +
  scale_colour_manual(values = line_cols, name = NULL) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 7)) +
  facet_wrap(~ model_dir, ncol = 7, scales = "free_x") +
  coord_cartesian(ylim = c(6, 12)) +
  labs(y = expression(log[10]*"(cells/mL)"),
       title = "B. 20°C") +
  theme_bw(base_size = 16)
ggsave(filename = "figures/FigS4Bmono_counts_N_20C.pdf", plot = last_plot(), height = 12, width = 20)

  
### plot monos vs. them in pairwise (doesn't provide much to interpret for a visual)
# mass_long %>%
#   filter(Temp == "10C") %>%
##  filter(Temp == "20C") %>%
#   ggplot(aes(x = Day)) +
#   # summed_prob — solid lines
#   stat_summary(
#     aes(y = log10(summed_prob), group = Isolate, fill = Isolate),
#     fun.data = mean_cl_normal, geom = "ribbon", alpha = 0.15, colour = NA
#   ) +
#   stat_summary(
#     aes(y = log10(summed_prob), group = Isolate, colour = Isolate, linetype = "summed_prob"),
#     fun = mean, geom = "line", linewidth = 1.2
#   ) +
#   # mono_counts — dashed lines
#   stat_summary(
#     aes(y = log10(mono_counts), group = Isolate, fill = Isolate),
#     fun.data = mean_cl_normal, geom = "ribbon", alpha = 0.15, colour = NA
#   ) +
#   stat_summary(
#     aes(y = log10(mono_counts), group = Isolate, colour = Isolate, linetype = "mono_counts"),
#     fun = mean, geom = "line", linewidth = 1.2
#   ) +
#   scale_fill_manual(values = isolate_cols, name = "Isolate") +
#   scale_colour_manual(values = isolate_cols, name = "Isolate") +
#   scale_linetype_manual(values = c("summed_prob" = "solid", "mono_counts" = "dashed"), name = "Metric") +
#   scale_x_continuous(breaks = c(0, 1, 2, 3, 7)) +
#   facet_wrap(~ model_dir, ncol = 7, scales = "free_x") +
#   coord_cartesian(ylim = c(6, 12)) +
#   theme_bw(base_size = 16)
# ggsave(filename = "figures/summed_prob_mono_counts_10C.pdf", plot = last_plot(), height = 12, width = 20)
## ggsave(filename = "figures/summed_prob_mono_counts_20C.pdf", plot = last_plot(), height = 12, width = 20)






#### AUC - PAIRS ####
## calculate area under curve AUC 
mass_long %>% 
  filter(Day != 0) %>% 
  filter(if_any(everything(), is.na))

auc_pairs <- mass_long %>% 
  select(-c(count_type)) %>% 
  pivot_wider(
    id_cols      = c(Community, Well, Temp, Isolate, model_dir,  Round),
    names_from   = Day,
    values_from  = summed_prob, #this is then calculating for each isolate separate 
    names_prefix = "summ_prob_d0"
  )

auc_pairs %>% filter(if_any(everything(), is.na)) ## PROBLEM...8 NA values for summed_prob_d07 (4 wells)

## manually calculate each section's rectangle (trapezoidal rule), then add together 
auc_pairs <- auc_pairs %>%
  mutate(
    auc_1 = summ_prob_d00 + summ_prob_d01 / 2 ,
    auc_2 = (summ_prob_d01 + summ_prob_d02) / 2,
    auc_3 = (summ_prob_d02 + summ_prob_d03) / 2,
    auc_7 = ((summ_prob_d03 + summ_prob_d07)*4) / 2,
    total_pair_auc = auc_1 + auc_2+ auc_3 + auc_7 
  )

auc_pairs %>%
  ggplot(aes(x = Temp, y = log(total_pair_auc), colour = Temp, fill = Temp)) +
  scale_colour_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  geom_jitter(width = 0.2) +
  theme_bw(base_size = 16)
  
auc_pairs_wide <- auc_pairs %>% #calculate difference between 20 vs. 10 
  select(Community, Well, Temp, model_dir, Round, Isolate, total_pair_auc) %>%
  pivot_wider(
    names_from  = Temp,
    values_from = total_pair_auc,
    names_prefix = "auc_"
  ) %>%
  mutate(auc_pair_diff = auc_20C - auc_10C)

auc_pairs_wide%>%
  ggplot(aes(x = Isolate, y = auc_pair_diff, colour = model_dir)) +
  geom_jitter(width = 0.2) +
  geom_hline(yintercept =  0) +
  labs(y = "20C-10C area under curve PAIRS")+
  theme_bw(base_size = 16)
  




#### AUC - MONO ?? NEEDS FIXING ??  ####
auc_monos <- mass_long %>% ### maybe here 
  select(-c(count_type)) %>% # need to remove, otherwise d00 fake counts don't covnert to wide format
  pivot_wider(
    id_cols      = c(Community, Well, Temp, Isolate, model_dir, Round),
    names_from   = Day,
    values_from  = mono_counts,
    names_prefix = "mono_d0"
  )

# manually calculate each section's rectangle (trapezoidal rule), then add together
auc_monos<- auc_monos %>%
  mutate(
    auc_1 = (mono_d00 + mono_d01) / 2,
    auc_2 = (mono_d01 + mono_d02) / 2, 
    auc_3 = (mono_d02 + mono_d03) / 2,
    auc_7 = ((mono_d03 + mono_d07)*4)/2,
    total_mono_auc = auc_1 + auc_2+ auc_3 + auc_7 
  )

# this df for statistics need to collapse for monoculture AUC stats
auc_mono_collapsed <- auc_monos %>% 
  group_by(Isolate, Temp, Round) %>%
  summarise(
    total_mono_auc = mean(total_mono_auc),
    .groups = "drop"
  )

auc_mono_collapsed %>%
  ggplot(aes(x = Isolate, y = total_mono_auc, colour = Temp, fill = Temp)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
    size = 5, alpha = 0.7
  ) +
  stat_summary(
    aes(group = Temp),
    fun = median,
    geom = "crossbar",
    width = 0.5,
    fill = NA,
    position = position_dodge(width = 0.75)
  ) +
  scale_colour_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  theme_bw(base_size = 16)

# auc_mono_wide <- auc_monos %>% #calculate difference between 20 vs. 10 
#   select(Community, Well, Temp, model_dir, Round, Isolate, total_mono_auc) %>%
#   pivot_wider(
#     names_from  = Temp,
#     values_from = total_mono_auc,
#     names_prefix = "auc_"
#   ) %>%
#   mutate(auc_mono_diff = auc_20C - auc_10C)

### need to collapse because there's multiple repeated lines for unique Isolate*model_dir 
### SHOULD ONLY BE 3 POINTS 
# auc_mono_wide_collapsed <- auc_mono_wide %>% 
#   group_by(Isolate, model_dir) %>%
#   summarise(
#     auc_20C      = mean(auc_20C),
#     auc_10C      = mean(auc_10C),
#     auc_mono_diff = mean(auc_mono_diff),
#     .groups = "drop"
#   )

# auc_mono_wide_collapsed%>% ### SHOULD ONLY BE 3 POINTS 
#   ggplot(aes(x = Isolate, y = auc_mono_diff)) +
#   geom_jitter(width = 0.2) +
#   geom_hline(yintercept =  0) +
#   labs(y = "20C-10C area under curve MONOS")+
#   theme_bw(base_size = 16)




#### AUC - total N pairs ####
auc_N <- mass_long %>%
  select(-c(count_type)) %>%
  pivot_wider(
    id_cols      = c(Community, Well, Temp, Isolate, model_dir, Round),
    names_from   = Day,
    values_from  = N,
    names_prefix = "N_d0"
  )

# manually calculate each section's rectangle (trapezoidal rule), then add together
auc_N <- auc_N %>%
  mutate(
    auc_1 = (N_d00 + N_d01) / 2,
    auc_2 = (N_d01 + N_d02) / 2,
    auc_3 = (N_d02 + N_d03) / 2,
    auc_7 = ((N_d03 + N_d07) * 4) / 2,
    total_auc_N = auc_1 + auc_2 + auc_3 + auc_7
  )

auc_N %>%
  ggplot(aes(x = Temp, y = log(total_auc_N), colour = Temp, fill = Temp)) +
  scale_colour_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  geom_jitter(width = 0.2) +
  theme_bw(base_size = 16)

auc_N_wide <- auc_N %>% # calculate difference between 20 vs. 10
  select(Community, Well, Temp, model_dir, Round, Isolate, total_auc_N) %>%
  pivot_wider(
    names_from   = Temp,
    values_from  = total_auc_N,
    names_prefix = "auc_"
  ) %>%
  mutate(auc_diff_N = auc_20C - auc_10C)

auc_N_wide %>%
  ggplot(aes(x = Isolate, y = auc_diff_N, colour = model_dir)) +
  geom_jitter(width = 0.2) +
  geom_hline(yintercept = 0) +
  labs(y = "20C-10C area under curve PAIRS (N)") +
  theme_bw(base_size = 16)


#### AUC - STATS ####
## MONOS 
auc_mono_collapsed <- auc_monos %>% ### need to collapse for monoculture AUC stats
  group_by(Isolate, Temp, Round) %>%
  summarise(
    total_mono_auc = mean(total_mono_auc),
    .groups = "drop"
  )

aov_mono_auc <- aov(total_mono_auc ~ Temp*Isolate, data = auc_mono_collapsed) ### also include effect of model_dir ?? (i.e., which pair )
summary(aov_mono_auc)

## Tukey's HSD test
# tukey_mono_auc <- TukeyHSD(aov_mono_auc)
# tukey_mono_auc
# # plot(tukey_mono_auc)
# #export
# tukey_auc_mono <- as.data.frame(tukey_mono_auc$`Temp:Isolate`)
# write.csv(tukey_auc_mono, "tukey_temp_isolate_mono_auc.csv", row.names = TRUE)

## t.tests of Flavobacterium AUC
t_test_df2 <- auc_mono_collapsed %>% filter(Isolate == "56_Flavobacterium")
t.test(total_mono_auc ~ Temp, data = t_test_df2) 

## pairs total cell density N 
aov_pair_N_auc <- aov(total_auc_N ~ Temp*model_dir, data = auc_N) ### also include effect of model_dir (i.e., which pair )
summary(aov_pair_N_auc)




#### COMPETITIVE RANKS ####

# ## more basic approach, just as mean AUC across everything (not used currently)
# ## calculate mean and sd for 
# auc_pair_summary <- auc_pairs_wide %>%
#   group_by(Isolate) %>%
#   summarise(
#     mean_auc_20C = mean(auc_20C, na.rm = TRUE),
#     sd_auc_20C   = sd(auc_20C, na.rm = TRUE),
#     mean_auc_10C = mean(auc_10C, na.rm = TRUE),
#     sd_auc_10C   = sd(auc_10C, na.rm = TRUE),
#     n = n(),
#     .groups = "drop"
#   )
# auc_pair_summary
# 
# # Rank orders
# auc_pair_summary %>% arrange(desc(mean_auc_20C))
# auc_pair_summary %>% arrange(desc(mean_auc_10C))


## rank as number of pairs where AUC is significantly greater than its partner/competitor isolate 
## code from CLAUDE
auc_pairs_wide <- auc_pairs_wide %>%
  group_by(Community, Round, model_dir, Well) %>%
  mutate(
    auc_10C_is_max = auc_10C == max(auc_10C, na.rm = TRUE),
    auc_20C_is_max = auc_20C == max(auc_20C, na.rm = TRUE)
  ) %>%
  ungroup()

# check if always a TRUE & FALSE for each model_dir (currently problem here...)
# I need to remember/understand what this is testing here... 
auc_pairs_wide %>%
  group_by(Community, Round, model_dir, Well) %>%
  summarise(n = n(), 
            n_true_10C = sum(auc_10C_is_max, na.rm = TRUE), 
            n_true_20C = sum(auc_20C_is_max, na.rm = TRUE), 
            .groups = "drop") %>%
  filter(n != 2 | n_true_10C != 1 | n_true_20C != 1)

# count number of wins for each Isolate
win_counts <- auc_pairs_wide %>%
  group_by(Isolate) %>%
  summarise(
    n_wins_10C = sum(auc_10C_is_max, na.rm = TRUE),
    n_wins_20C = sum(auc_20C_is_max, na.rm = TRUE),
    n_total = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_wins_20C))
win_counts

# make long format (needed to plot and later for correlation)
wins_long <- win_counts %>%
  pivot_longer(
    cols = c(n_wins_10C, n_wins_20C),
    names_to = "Temp",
    values_to = "n_wins"
  ) %>%
  mutate(Temp = gsub("n_wins_", "", Temp))  # -> "10C", "20C"

## try plotting
wins_long %>%
  ggplot(aes(x = n_wins, y = reorder(Isolate, n_wins), fill = Temp)) +
  scale_fill_manual(
    values = c("10C" = "#00AFBB", "20C" = "#FC4E07")) +
  geom_col(position = "dodge", alpha = 0.7) +
  labs(title = "Isolate wins by temperature", x = "Number of wins", y = NULL) +
  theme_bw(base_size = 16)





#### CORRELATE monos ~ rank ####
## Get mean r and K for each monos ... but maybe stats better to not reduce to mean ???
MONOS_mean <- MONOS_summary %>%
  group_by(Isolate, Temp) %>%
  summarise(mean_mono_max = mean(mono_max, na.rm = TRUE), 
            mean_mono_r = mean(mono_r_d03, na.rm = TRUE),
            .groups = "drop")

# Join to win counts/ranks data
combined <- MONOS_mean %>%
  left_join(wins_long, by = c("Isolate", "Temp"))

## correlate RATE for 10C and 20C
cor.test(combined$mean_mono_r[combined$Temp == "10C"],
         combined$n_wins[combined$Temp == "10C"],
         # method = "spearman" #S = 8, p-value = 0.02381, rho 0.8571429
         method = "kendall" #T = 18, p-value = 0.03016, tau 0.7142857 
)

cor.test(combined$mean_mono_r[combined$Temp == "20C"],
         combined$n_wins[combined$Temp == "20C"],
         # method = "spearman" #S = 8.5747, p-value = 0.0162, rho 0.8468812 
         method = "kendall" #z = 2.1268, p-value = 0.03344, tau 0.6831301 
)

## correlate maximum reached "K" for 10C and 20C
cor.test(combined$mean_mono_max[combined$Temp == "10C"],
         combined$n_wins[combined$Temp == "10C"],
         # method = "spearman" #S = 10, p-value = 0.03413, rho 0.8214286
         method = "kendall" # T = 18, p-value = 0.03016, tau 0.7142857
         )

cor.test(combined$mean_mono_max[combined$Temp == "20C"],
         combined$n_wins[combined$Temp == "20C"],
         # method = "spearman" #S = 8.5747, p-value = 0.0162,  rho  0.8468812 
         method = "kendall" #z = 2.1268, p-value = 0.03344,  tau 0.6831301
)


### Spearman and Kendall only depend on the rank order of the values, not their magnitudes. 
## So if mean_mono_max and mean_mono_r happen to put the isolates in the exact same order within a given Temp group,
## you'll get identical rho/tau/p-values even though the underlying numbers are completely different.

# ## from Claude: permutation test is another option
# set.seed(123)
# obs_rho <- cor(combined$mean_mono_max[combined$Temp == "20C"],
#                combined$n_wins[combined$Temp == "20C"],
#                method = "spearman")
# 
# n_perm <- 10000
# perm_rhos <- replicate(n_perm, {
#   cor(combined$mean_mono_max[combined$Temp == "20C"],
#       sample(combined$n_wins[combined$Temp == "20C"]),
#       method = "spearman")
# })
# 
# perm_p <- mean(abs(perm_rhos) >= abs(obs_rho))
# perm_p





#### INTERACTIONS from AUC #####
## combine all the AUC data
AUC_INTERACT <- auc_monos %>%
  left_join(auc_pairs, by = c("Community", "Well", "Temp", "Isolate", "model_dir", "Round")) %>%
  left_join(auc_N, by = c("Community", "Well", "Temp", "Isolate", "model_dir", "Round")) %>%
  select(Community, Well, Temp, Isolate, model_dir, Round, total_mono_auc, total_pair_auc, total_auc_N)

## INTERACTION = N vs. expected from monos (i.e., mean of 2 respective monos)
AUC_INTERACT <- AUC_INTERACT %>%
  group_by(Community, Well, Temp, model_dir, Round) %>%
  mutate(mean_auc = mean(total_mono_auc, na.rm = TRUE)) %>% #calculate average auc of the 2 monocultures (i.e., isolate A and B). Equivalent of (monoA + mono)/2
  ungroup()

AUC_INTERACT <- AUC_INTERACT%>%
  mutate(
    pair_interaction = total_auc_N / mean_auc, #mean_auc is expected, assumes complete resource sharing. 
    isolate_interaction = total_pair_auc / mean_auc  )

### Martina/RIVER interaction calculation (relative to highest AUC of either monoculture). Done on replicate or "Well" level. 
### When I plot same heat map, all positive pair level, but most not significant, so doesnt change conclusions
# AUC_INTERACT <- AUC_INTERACT %>%
#   group_by(Community, Well, Temp, model_dir, Round) %>%
#   mutate(max_mono_auc = max(total_mono_auc, na.rm = TRUE)) %>%
#   ungroup()
# 
# AUC_INTERACT <- AUC_INTERACT%>%
#   mutate( 
#     pair_interaction = total_auc_N / max_mono_auc,
#     isolate_interaction = total_pair_auc / max_mono_auc  )


#### Plotting INTERACTIONS ####
## plot mean of each isolate
AUC_INTERACT %>% 
  ggplot(aes(x = Isolate, y = log(isolate_interaction), colour = Temp, fill = Temp)) +
  geom_point(
    size = 2, alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6)
  ) +  geom_boxplot(alpha = 0.2) +
  scale_colour_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07"),
    labels = c("10C" = "10°C", "20C" = "20°C"),
    name = "Temperature") +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07"),
    labels = c("10C" = "10°C", "20C" = "20°C"),
    name = "Temperature") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = expression(log[10]*"(Interaction)"),
       title = "B. Isolate-level Interaction") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename = "figures/Fig6A_Isolate_interaction.pdf", plot = last_plot(), height = 7, width = 8)
ggsave(filename = "figures/Fig6A_Isolate_MARTINA.pdf", plot = last_plot(), height = 7, width = 8)



## plot mean of pair-level interaction
AUC_INTERACT %>% 
  ggplot(aes(x = Temp, y = log(pair_interaction), colour = Temp, fill = Temp)) +
  geom_point(
    size = 2, alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.6)
  ) +  geom_boxplot(alpha = 0.2) +
  scale_colour_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  scale_fill_manual(
    values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename = "figures/pairlevel_mean_interaction.pdf", plot = last_plot(), height = 7, width = 8)


## pair-level interaction
AUC_INTERACT %>%
  distinct(Community, Well, Temp, model_dir, Round, pair_interaction) %>%
  separate(model_dir, into = c("Isolate1", "Isolate2"), sep = "_vs_", remove = FALSE) %>%
  mutate(log_comm_interaction = log(pair_interaction)) %>%
  group_by(Isolate1, Isolate2, Temp) %>%
  summarise(
    mean_log_comm_interaction = mean(log_comm_interaction, na.rm = TRUE),
    n = sum(!is.na(log_comm_interaction)),
    p_value = if (n >= 2) t.test(log_comm_interaction, mu = 0)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(sig_label = case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ ""
  )) %>%
  ggplot(aes(x = Isolate1, y = Isolate2, fill = mean_log_comm_interaction)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sig_label), colour = "black", size = 4, vjust = 0.7) +
  # scale_fill_gradient2(midpoint = 0, low = "#1D457F", mid = "white", high = "#E7B800") +
  scale_fill_gradient2(
    low = "#1D457F",
    mid = "white",
    high = "#E7B800",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    na.value = "grey90"
  ) +
  facet_wrap(~Temp) +
  labs(
    x = "Isolate",
    y = "Isolate",
    fill = "log(Interaction)",
    title = "A. Pair-level"
  ) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "figures/Fig6A_pairlevel_interaction.pdf", plot = last_plot(), height = 7, width = 14)
ggsave(filename = "figures/Fig6A_pairlevel_interaction.pdf", plot = last_plot(), height = 7, width = 14)


### ISOLATE SPECIFIC INTERACTIONS
AUC_INTERACT %>%
  group_by(model_dir, Isolate, Temp) %>%
  summarise(mean_log_isolate_interaction = mean(log(isolate_interaction), na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = model_dir, y = Isolate, fill = mean_log_isolate_interaction)) +
  geom_tile() +
  # scale_fill_viridis_c(option = "viridis") 
  scale_fill_gradient2(midpoint = 0, low = "#1D457F", mid = "white", high = "#E7B800", na.value = "grey90") +
  facet_wrap(~Temp) +
  labs(fill = "Mean log(isolate\ninteraction)") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### PLOT EVERYTHING TOGETHER  :)
### complicated code because now "Isolate" has num_genus naming. I.e., 12_Sphingomonas . ClaudeAI generated... but seems to work. 
plot_data <- bind_rows(
  AUC_INTERACT %>%
    # pair interaction
    distinct(model_dir, Temp, Community, Well, Round, pair_interaction) %>%
    mutate(log_value = log(pair_interaction)) %>%
    group_by(model_dir, Temp) %>%
    summarise(
      value    = mean(log_value, na.rm = TRUE),
      n        = sum(!is.na(log_value)),
      p_value  = if (n >= 2) t.test(log_value, mu = 0)$p.value else NA_real_,
      .groups  = "drop"
    ) %>%
    mutate(y_label = "Pair-level", isolate_text = ""),
  
# individual isolates' interaction 
  AUC_INTERACT %>%
    separate(model_dir, into = c("Isolate1", "Isolate2"), sep = "_vs_", remove = FALSE) %>%
    mutate(isolate_num      = str_extract(Isolate, "^[0-9]+"),
           isolate_position = if_else(isolate_num == Isolate1, "Isolate A", "Isolate B"),
           log_value        = log(isolate_interaction)) %>%
    group_by(model_dir, Temp, isolate_position) %>%
    summarise(
      value        = mean(log_value, na.rm = TRUE),
      n            = sum(!is.na(log_value)),
      p_value      = if (n >= 2) t.test(log_value, mu = 0)$p.value else NA_real_,
      isolate_text = paste0(str_extract(first(Isolate), "^[0-9]+"),
                            str_extract(first(Isolate), "(?<=_)[A-Za-z]{1,3}")),
      .groups      = "drop"
    ) %>%
    rename(y_label = isolate_position)
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig_label = case_when(
      is.na(p_adj)   ~ "",
      p_adj < 0.001  ~ "***",
      p_adj < 0.01   ~ "**",
      p_adj < 0.05   ~ "*",
      TRUE           ~ ""
    )
  )

#### plot interactions as heatplot
## non-significant shown as grey
plot_data %>%
  mutate(y_label = factor(y_label, levels = c("Isolate A", "Isolate B", "Pair-level")),
         fill_value = if_else(sig_label != "", value, NA_real_)) %>%
  ggplot(aes(x = model_dir, y = y_label, fill = fill_value)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(midpoint = 0, low = "#1D457F", mid = "white", high = "#E7B800",
                       na.value = "grey90") +
  facet_wrap(~Temp, nrow = 2) +
  labs(x = "Pair", y = NULL, fill = "Mean log(interaction)") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "figures/FigS5_AUC_interactions_SIGNIF.pdf", plot = last_plot(), height = 5, width = 10)



#### pairlevel INTERACTION stats ####
### checking assumptions https://statsandr.com/blog/anova-in-r/

## t-tests & anova - pairlevel
## t.test if pairlevel interactions differ from 0
# 10C
t.test(log(pair_interaction) ~ 1, 
       data = subset(AUC_INTERACT, Temp == "10C"), 
       mu = 0)
# 20C
t.test(log(pair_interaction) ~ 1, 
       data = subset(AUC_INTERACT, Temp == "20C"), 
       mu = 0)

res_aov <- aov(log(pair_interaction) ~ Temp*model_dir + Round , data = AUC_INTERACT)
summary(res_aov)

# Marginal means for each effect
library(emmeans)
emmeans(res_aov, ~ Temp)
emmeans(res_aov, ~ model_dir) #doesn't tell alot? because 10C and 20C grouped
emmeans(res_aov, ~ Round)
emmeans(res_aov, ~ Temp*model_dir, type = "response")

### testing normality, etc. 
hist(res_aov$residuals) # looks normal enough?
shapiro.test(res_aov$residuals) # p = 0.00137

## QQ-plot
# qqPlot(res_aov$residuals,
#        id = FALSE # id = FALSE to remove point identification
# )
 
#### Kruskall Wallis - pairlevel interaction
# used for non-normal data https://statsandr.com/blog/kruskal-wallis-test-nonparametric-version-anova/
kruskal.test(log(pair_interaction) ~ Temp , data = AUC_INTERACT) #p-value = 5.958e-06
kruskal.test(log(pair_interaction) ~ model_dir , data = AUC_INTERACT) #p-value = 1.769e-10
kruskal.test(log(pair_interaction) ~ Round , data = AUC_INTERACT) #p-value = 0.03315



#### Isolate INTERACTION stats #####
##  ANOVA - isolate interaction
res_aov <- aov(log(isolate_interaction) ~ Temp + Isolate + model_dir + Round + Isolate*Temp + Isolate*model_dir, data = AUC_INTERACT)
summary(res_aov)

# Marginal means for each effect
library(emmeans)
emmeans(res_aov, ~ Temp)
emmeans(res_aov, ~ Isolate) 
emmeans(res_aov, ~ Isolate*Temp, type = "response") 
emmeans(res_aov, ~ Isolate*model_dir, type = "response") 
emmeans(res_aov, ~ Round)

### testing normality, etc. 
hist(res_aov$residuals) # looks normal enough
shapiro.test(res_aov$residuals)  #p-value = 0.03183


#### Kruskall Wallis - isolate interaction
kruskal.test(log(isolate_interaction) ~ Temp , data = AUC_INTERACT) #p-value = 0.002339
kruskal.test(log(isolate_interaction) ~ Isolate , data = AUC_INTERACT) #p-value = 2.145e-07
kruskal.test(log(isolate_interaction) ~ model_dir , data = AUC_INTERACT) #p-value = 0.0001167
kruskal.test(log(isolate_interaction) ~ Round , data = AUC_INTERACT) #p-value = 0.1453


####  Wilcoxon Rank tests - isolate interaction
wilcoxon <- AUC_INTERACT %>%
  filter(!is.na(isolate_interaction)) %>%
  mutate(log_interaction = log(isolate_interaction)) %>%
  group_by(Isolate, Temp) %>%
  summarise(
    n = n(),
    median_log = median(log_interaction),
    p_raw = wilcox.test(log_interaction, mu = 0, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_raw, method = "BH"), #Benjamini-Hochberg FDR-corrected p-value (recommended for interpreting significance across the 14 tests)
    p_adj_bonferroni = p.adjust(p_raw, method = "bonferroni"), #stricter p-value adjustment alternative
    sig = p_adj_BH < 0.05
  )
print(wilcoxon)
write.csv(wilcoxon, "wilcoxon_results_isolate_interactions.csv", row.names = FALSE)


















#### OLD CODE NOT BEING USED ####

# 
# #### #### plotting kruskall wallis tests 
# library(ggstatsplot)
# AUC_INTERACT %>%
#   mutate(log_Interaction = log(pair_interaction)) %>%
#   ggbetweenstats(
#     x = Temp,
#     # y = Interaction,
#     y = log_Interaction,
#     type = "nonparametric",
#     plot.type = "box",
#     pairwise.comparisons = TRUE,
#     pairwise.display = "significant",
#     centrality.plotting = FALSE,
#     bf.message = FALSE
#   ) +
#   scale_colour_manual(
#     values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
#   scale_fill_manual(
#     values = c("10C" = "#00AFBB","20C" = "#FC4E07")) +
#   theme_bw(base_size = 16) +
#   labs(title = "pair Interaction")
# # ggsave(filename = "figures/KW_comm_inter_temp.pdf", plot = last_plot(), height = 8, width = 10)
# 
# AUC_INTERACT %>%
#   mutate(log_Interaction = log(pair_interaction)) %>%
#   ggbetweenstats(
#     x = model_dir,
#     y = log_Interaction,
#     type = "nonparametric",
#     plot.type = "box",
#     pairwise.comparisons = TRUE,
#     pairwise.display = "significant",
#     centrality.plotting = FALSE,
#     bf.message = FALSE
#   ) +
#   theme_bw(base_size = 16) +
#   labs(title = "Pair Interaction")
# # ggsave(filename = "figures/KW_comm_inter_pairs.pdf", plot = last_plot(), height = 8, width = 10)
# 
# AUC_INTERACT %>%
#   mutate(log_Interaction = log(pair_interaction)) %>%
#   ggbetweenstats(
#     x = Round,
#     y = log_Interaction,
#     type = "nonparametric",
#     plot.type = "box",
#     pairwise.comparisons = TRUE,
#     pairwise.display = "significant",
#     centrality.plotting = FALSE,
#     bf.message = FALSE
#   ) +
#   theme_bw(base_size = 16)+
#   labs(title = "Pair Interaction")
# # ggsave(filename = "figures/KW_comm_inter_Round.pdf", plot = last_plot(), height = 8, width = 10)
# 
# 
# 




# ### RELATIVE ABUNDANCES - 2023 SOILS ####
# # this from 16S data of Laureen. Created in separate R script "2023_alpine_soils_16S.R"
# genus <- read.table('AlpSoils23_16S_genus_relab.csv', sep = ',', header = TRUE)
# META_soils <- read.table('Fig2_AB_AlpineSoil23_metadata.tsv', sep = '\t', header = TRUE)
# ## list Genuses of the 7 isolates (Pseudarthrobacter and Arthrobacter both for isolate 16)
# names_to_keep <- c("Sphingomonas", "Pseudarthrobacter", "Arthrobacter", "Paenibacillus", "Variovorax", "Luteibacter", "Psychrobacillus", "Flavobacterium")
# 
# genus_filter <- genus %>%
#   filter(str_detect(Genus, paste(names_to_keep, collapse = "|"))) 
# 
# GENUS <- t(genus_filter)  #transpose dataset            
# colnames(GENUS) <- GENUS[1, ]   # make first row column names
# GENUS <- GENUS[-1, ]            # then remove it
# GENUS <- as.data.frame(GENUS) %>% #transposing made it not a df
#   rownames_to_column("ID")%>%     # make row names 
# mutate(across(-ID, as.numeric))
# 
# # nrow(GENUS) #86 rows 
# # nrow(META_soils) # 70 rows...so META not for all soils/sites 
# 
# # try filtering just sites of interest (Moiry = C, Bettmerkorn = K, Plaine Morte = L)
# GENUS <- GENUS %>% filter(str_detect(ID, "^(C_|K_|L_)"))
# META_soils <- META_soils %>% filter(str_detect(ID, "^(C_|K_|L_)"))
# # 
# # nrow(GENUS) #29 rows 
# # nrow(META_soils) # 26 rows...so META not for all soils/sites 
# 
# 
# GENUS_2023 <- left_join(GENUS, META_soils, by = "ID")
# GENUS_2023 <- GENUS_2023[, 2:10] # remove anything after row 11 (just want site code), remove ID
# 
# GENUS_2023_mean <- GENUS_2023 %>%#take average for all of site (regardless of subsite, depth, etc.)
#   group_by(site) %>%
#   summarise(across(where(is.numeric), mean))











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
# mass_mean <- mass_long %>%
#   # group_by(model_dir, Temp, Day, isolate, count_type) %>%
#   group_by(model_dir, Round, Temp, Day, isolate, count_type) %>%
#   summarise(count = mean(count), .groups = "drop")
# 
# mass_mean %>% 
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
# ggsave(filename = "figures/mean_DAY3_7.pdf", plot = last_plot(), height = 8, width = 16)

# 

### relative abundances 
# mass_mean %>%
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




#### extra interaction plots playing #####

# # separate by day
# INTERACT %>% 
# #   # filter(Interaction > -2000)  %>%
# #   ggplot(aes(x = pair, y = Interaction, colour = Temp)) +
# #   geom_jitter(width = 0.2, size = 2.5, alpha = 0.5) + 
# #   scale_colour_manual(
# #     values = c("10C" = "skyblue3","20C" = "orange2"
# #   )) +
# #   facet_wrap(~Day) +
# #   theme_bw(base_size = 16)
# # ggsave(filename = "figures/interaction_day_comm.pdf", plot = last_plot(), height = 8, width = 10)
# # # ggsave(filename = "interaction_day_comm_out.png", plot = last_plot(), height = 8, width = 10)
# # 
# # 
# # INTERACT %>% 
# #   # filter(Interaction > -2000)  %>%
# #   ggplot(aes(x = pair, y = log(Interaction), colour = Temp)) +
# #   geom_jitter(width = 0.2, size = 2.5, alpha = 0.5) + 
# #   scale_colour_manual(
# #     values = c("10C" = "skyblue3","20C" = "orange2"
# #     )) +
# #   facet_wrap(~Day) +
# #   theme_bw(base_size = 16)
# # ggsave(filename = "figures/interaction_day_comm_log.pdf", plot = last_plot(), height = 8, width = 10)
# # # ggsave(filename = "interaction_day_comm_out.png", plot = last_plot(), height = 8, width = 10)
# 
# 
# ### separate by pair 
# # INTERACT %>% 
# #   ggplot(aes(x = Day, y = Interaction, colour = Temp)) +
# # geom_jitter(width = 0.2) + 
# # scale_colour_manual(
# #     values = c("10C" = "skyblue3","20C" = "orange2"
# #     )) +
# #   facet_wrap(~pair) +
# #   theme_bw(base_size = 16)
# 
# 
# 
# 
# 
# 
# ##### Standardise to monos ######
# ### nope, this also doesnt make sense.... 
# mass_long <- mass_long %>%
#   mutate(mono_stand = summed_prob/mono_counts) 
# 
# mass_long %>%
#   # filter(Round == "round1") %>%
#   filter(Temp == "10C") %>%
#   mutate(Day = readr::parse_number(as.character(Day))) %>%
#   ggplot(aes(x = Day, y = log10(mono_stand))) +
#   scale_fill_viridis_d(option = "turbo") +
#   scale_colour_viridis_d(option = "turbo") +
#   stat_summary(
#     aes(group = interaction(Isolate, Round), fill = Isolate),
#     fun.data = mean_cl_normal,
#     geom = "ribbon",
#     alpha = 0.2,
#     colour = NA
#   ) +
#   stat_summary(
#     aes(group = interaction(Isolate, Round), colour = Isolate, linetype = Round),
#     fun = mean,
#     geom = "line",
#     linewidth = 1.2
#   ) +
#   scale_linetype_manual(values = c("round1" = "solid", "round2" = "dashed", "round3" = "dotted")) +  # adjust level names to match yours
#   scale_x_continuous(breaks = c(1, 2, 3, 7)) +
#   coord_cartesian(ylim = c(-2.5, 2.5)) +
#   # facet_grid(Temp ~ model_dir, scales = "free_x") +
#   facet_wrap(~ model_dir + Temp, ncol = 6, scales = "free_x")+
#   theme_bw(base_size = 16) +
#   theme(
#     # axis.text.x = element_text(angle = 90, hjust = 0, size = 12),
#     # strip.text.x = element_text(angle = 90, hjust = 0, size = 12),
#     # strip.text.y = element_text(angle = 0, hjust = 0, size = 12)
#   )
# ggsave(filename = "figures/mono_standardised.pdf", plot = last_plot(), height = 12, width = 20)



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


# #### testing if correlation between later outcome to d01 ####
# library(ggpmisc) ### for adding p-values
# 
# MASS_MODEL <- MASS_MODEL %>%
#   mutate(AB_ratio = summed_proba_A / summed_proba_B)
# 
# # need to set up for a regression 
# MASS_wide <- MASS_MODEL %>%
#   select(Community, model_dir, Round, Temp, Well, Day, AB_ratio) %>%
#   pivot_wider(names_from = Day,
#               values_from = AB_ratio,
#               names_prefix = "day_")
# 
# ## day 3 ~ day 1
# lm_d1_d3 <- lm(day_d03 ~ day_d01 +Community + Round, data = MASS_wide) # want to control for each community , not sure if this is right 
# # lm_d1_d3 <- lm(day_d03 ~ day_d01, data = MASS_wide) 
# summary(lm_d1_d3) # day 1 significant... 
# 
# MASS_wide %>% 
#   ggplot(aes(log(day_d01), log (day_d03), colour = Community, shape = Round)) +
#   geom_point(size = 4, alpha = 0.7) +
#   geom_smooth(method = "lm", se = TRUE,
#               colour = "black", linewidth = 1,
#               inherit.aes = FALSE,    # want just one big grouping mes                      
#               aes(x = log(day_d01), y = log(day_d03))) +
#   theme_bw(base_size = 16) +
#   labs(x = "log(AB_ratio) Day 1",
#        y = "log(AB_ratio) Day 3")
# ggsave(filename = "figures/AB_ratio_d01vsd03.pdf", plot = last_plot(), height = 8, width = 10)
# 
# 
# ## day 7 ~ day 1
# lm_d1_d7<- lm(day_d07 ~ day_d01 +Community + Round, data = MASS_wide)
# # lm_d1_d7<- lm(day_d07 ~ day_d01, data = MASS_wide)
# summary(lm_d1_d7) # day 1 NOT significant... 
# 
# MASS_wide %>% 
#   ggplot(aes(log(day_d01), log (day_d07), colour = Community, shape = Round)) +
#   geom_point(size = 4, alpha = 0.7) +
#   geom_smooth(method = "lm", se = TRUE,
#               colour = "black", linewidth = 1,
#               inherit.aes = FALSE, # want just one big grouping 
#               aes(x = log(day_d01), y = log(day_d07))) +
#   # +
#   #   stat_poly_eq(inherit.aes = FALSE,
#   #                aes(x = log(day_d01), y = log(day_d07),
#   #                    label = paste(after_stat(rr.label), 
#   #                                  after_stat(p.value.label), 
#   #                                  sep = "~~")),
#   #                label.x = "left", label.y = "top",
#   #                size = 5) +
#   theme_bw(base_size = 16) +
#   labs(x = "log(AB_ratio) Day 1",
#        y = "log(AB_ratio) Day 7")
# ggsave(filename = "figures/AB_ratio_d01vsd07.pdf", plot = last_plot(), height = 8, width = 10)


