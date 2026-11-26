# Synthetic Alpine bacterial community

# 2026 - experiments performed
Experiment with 7 bacterial isolates grown in monocultures and pairwise combos in R2A broth at 10 C and 20 C, with abundances determined by flow cytometry (FC). Experiment was repeated three independent times, each termed a "Round" (r01, r02, r03). 

## flow cytometer data
From Massimo https://github.com/MACE-laboratory/MACE_flow_cytometry_classifiers/
- all_pairwise_stats_round01.csv (accuracy calculations)
- all_pairwise_stats_round02.csv
- all_pairwise_stats_round03.csv
- coculture_pairwise_counts_round1.csv (output from Massimo's model. Used for R script)
- coculture_pairwise_counts_round2.csv
- coculture_pairwise_counts_round3.csv
- d00_r01 ?? (FC of monocultures used to mix pairs in 1:1 ratio at day00)
- d00_r02 ??
- d00_r03 ??

### coculture_pairwise_counts_round explained
- Round, Temp, Day, Date, Well, Community: all these columns match META file (see next section)
- isolate_A, isolate_B: 06, 12, 16, 20, 23, 32, 43, 56 corresponding to one of the eight unique bacterial isolates
- count_A, count_B: total "events" determined by FC model, using binary assignment 
- isolate_A_mono_counts, isolate_B_mono_counts: total "events" of monoculture 
- N: total number of events in pairwise wells 
- model_dir: "06_vs_12" matches as combo of isolate A and B, matches also to "Community"
- summed_proba_A, summed_proba_B: total "events" determined by FC model using summed probabilities 
- filename: original FC raw data file for that unique Round * Temp * Day * Well * Community
- Sample: 1-96 for 96 wells on plate 

### META excel file 
- Well      :	well on 96-well black flow cytometry plate 
- Community : numbers s1-s28 are each pair, or monos (ixx), syn comm with all 8, or negative control
- Day	      : days since inoculation (1, 2, 3, 7)
- Fixed_dil :	dilution factor for formaldehyde fix. Always 20:80 = 0.2
- Vol_10C   :	volume of fixed cells from 10C cultures (needed to calculate dilution factor in 200uL)
- Vol_20C	  : volume of fixed cells from 20C cultures (needed to calculate dilution factor in 200uL)
- Vol_Total :	total volume in FC well. Always 200uL.
- Temp	    : temperature treatment. 10 C or 20 C.
- r01_d01	11.03.2026 (round 01, day 01, date) 
- r01_d02	12.03.2026
- r01_d03	13.03.2026
- r01_d07	16.03.2026
- r02_d01	07.04.2026
- r02_d02	08.04.2026
- r02_d03	09.04.2026
- r02_d07	14.04.2026
- r03_d01	19.05.2026
- r03_d02	20.05.2026
- r03_d03	21.05.2026
- r03_d07	25.05.2026
- pair_meta tab =	monocultures used to train model, only 1 rep per temp*day

## FC validation tests 
Three Isolates were mixed in ratios of 0.25, 0.5, 0.75 (expected ratios) based from monoculture FC readings, then run on FC and model applied. 
Ratios calculated from model were then compared to expected
- fc_ratio_verification.R
- FC_ratio_classifer_overnight_counts_22.07.2026.csv (FC counts used to make 25, 50, 75% mixes)
- FC_ratio_test_info.xlsx (used to calculate volumes from overnight counts to make 25, 50, 75% ratios, also includes a META data) 
- META_FC_ratio_test.csv (just extracted from FC_ratio_test_info.xlsx)

### DNA extractions
  - tabs "DNA_extract_layout" & "DNA_concentrations" in excel "file FC_ratio_test_info.xlsx"
  - "DNA_extract_layout" : number refers to mixture number (i.e., 1 to 17 mixtures, listed in META)


## FC vs. OD600 pilots
Measured OD600 of cultures, then diluted to 25% concentration and OD600 measured again. Both 100% and 25% cultures then run on FC. FC counts vs. OD600 compared.
- OD_vs_FCcount_v1.R 
- OD_vs_FC_AML_OD600_96well_20260128_095233.xlsx (OD600 values)
- 16012026_alanna_isolate_counts.csv (FC counts)
- 29012026_alanna_isolate_counts.csv (FC counts)
- OD_vs_FC_layout.xlsx (META file)
