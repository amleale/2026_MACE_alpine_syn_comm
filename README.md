# Synthetic Alpine bacterial community

# 2026 - experiments performed
Experiment with 7 bacterial isolates in laboratory. Where all monocultrues and pairs grown at 10 C and 20 C, with abdundances determined by flow cytometry (FC)

## flow cytometer data
From Massimo https://github.com/MACE-laboratory/MACE_flow_cytometry_classifiers/
- all_pairwise_stats_round01.csv (accuracy calculations)
- all_pairwise_stats_round02.csv
- all_pairwise_stats_round03.csv
- coculture_pairwise_counts_round1.csv (used for R script)
- coculture_pairwise_counts_round2.csv
- coculture_pairwise_counts_round3.csv

## coculture_pairwise_counts_round explained
- Round, Temp, Day, Date, Well, Community: all these columns match META file (see next section)
- isolate_A, isolate_B: 06, 12, 16, 20, 23, 32, 43, 56 corresponding to one of the eight unique bacterial isolates
- count_A, count_B: total "events" determined by FC model, using binary assignment 
- isolate_A_mono_counts, isolate_B_mono_counts: total "events" of monoculture 
- N: total number of events in pairwise wells 
- model_dir: "06_vs_12" matches as combo of isolate A and B, matches also to "Community"
- summed_proba_A, summed_proba_B: total "events" determined by FC model using summed probabilities 
- filename: original FC raw data file for that unique Round * Temp * Day * Well * Community
- Sample: 1-96 for 96 wells on plate 

## META excel file 
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
- 
