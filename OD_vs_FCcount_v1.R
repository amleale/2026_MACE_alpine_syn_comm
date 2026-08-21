setwd("~/Library/CloudStorage/OneDrive-epfl.ch/2025_isolates_alpine_2023_soils/nov2025_jan2026_pilots/flow_cytometer/OD_vs_FC_pilot")

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr) 
library(vegan) #needed for permanova
library(ggrepel)
library(wesanderson)  
library(viridis)


#### importing FC counts #### 
counts1 <- read.csv("16012026_alanna_isolate_counts.csv", header = TRUE)
counts2 <- read.csv("29012026_alanna_isolate_counts.csv", header = TRUE)

meta1 <- readxl::read_xlsx(path = "OD_vs_FC_layout.xlsx", sheet = "dilutions_16jan2026") # these are all undiluted OD (so 100% OD)
meta2 <- readxl::read_xlsx(path = "OD_vs_FC_layout.xlsx", sheet = "dilutions_23jan2026") # these have 25% OD diluted and 100% 

counts1 <- left_join(counts1, meta1, by = "Well.ID")
counts2 <- left_join(counts2, meta2, by = "Well.ID")


COUNTS <- counts1 %>% full_join(counts2)
COUNTS$Date <- gsub(" .*", "", COUNTS$Run.Time)
COUNTS <- COUNTS %>% select(c(Well.ID, P1.Count, P1.Abs..Count, Date, Isolate, Sample_vol, Total_vol, OD_per))

# just test to see densities... what would be possible with FC and volumes 
max_min_abs.counts <- COUNTS %>%
  filter(!Isolate %in% c("25", "58", "54", "62", "MQ", "mq+", "PBS", "pbs+", "pbs+, F")) %>%
  filter(Date == "1/16/2026") %>%
  group_by(Isolate)  %>%
  summarise(
    max_count = max(P1.Abs..Count, na.rm = TRUE),
    mean_count = mean(P1.Abs..Count, na.rm = TRUE),
    min_count = min(P1.Abs..Count, na.rm = TRUE)
  )
max_min_abs.counts



## plot just undiluted samples 
COUNTS %>% 
  filter(!Isolate %in% c("25", "58", "MQ", "mq+", "PBS", "pbs+", "pbs+, F")) %>% # remove isolates that won't be used, and remove controls
  filter(OD_per == "p100") %>%
  ggplot(aes(
    x = Isolate,
    y = log ((P1.Abs..Count * (Total_vol / Sample_vol))*1000) # accounts for dilution onto plate i.e., 0.5 - 2 uL / 200uL. convert to be per mL
    # y = P1.Abs..Count * (Total_vol / Sample_vol) * 1000
  )) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    aes(colour = Date),
    width = 0.2,
    alpha = 0.7,
    size = 3
  ) +
  # geom_hline(yintercept = log(10^7)) + # add reference label
  labs(y = "log (cells/mL)") +
  # labs(y = "cells/mL") +
  theme_bw()
# ggsave(filename = "raw cells per mL.png", plot = last_plot(), height = 8, width = 10)



### now associate to OD measurements
ODs <- readxl::read_xlsx(path = "OD_vs_FC_AML_OD600_96well_20260128_095233.xlsx", sheet = "forR") 
# why is 56 missing ???

COUNTS_OD <- left_join(COUNTS, ODs, by = c("Isolate", "OD_per"))

COUNTS_OD %>% 
  filter(!Isolate %in% c("25", "58", "MQ", "mq+", "PBS", "pbs+", "pbs+, F")) %>% # remove isolates that won't be used, and remove controls
  ggplot(aes(
    colour = Isolate,
    alpha = OD_per,
    x = OD,
    # y = log ((P1.Abs..Count * (Total_vol / Sample_vol))*1000) # accounts for dilution onto plate i.e., 0.5 - 2 uL / 200uL. convert to be per mL
    y = log1p ((P1.Abs..Count * (Total_vol / Sample_vol))*1000) # deals with zeros
    # y = P1.Abs..Count * (Total_vol / Sample_vol) * 1000
  )) +
  geom_point(  size = 5) +
  geom_smooth(
    aes(group = Isolate),  # explicitly group by Isolate only
    method = "lm",
    se = FALSE,      # remove confidence interval shading
    linewidth = 0.8 ,
    alpha = 1         # override inherited alpha
  ) +
  scale_alpha_manual(values = c(1, 0.4)) +
  # geom_hline(yintercept = log(10^7)) + # add reference label
  labs(y = "log (cells/mL)",
       x = "OD600", 
       alpha = "culture dilution") +
  # labs(y = "cells/mL") +
  theme_bw(  base_size = 16) 
ggsave(filename = "FigS6_OD_vs_FC.png", plot = last_plot(), height = 8, width = 10)




