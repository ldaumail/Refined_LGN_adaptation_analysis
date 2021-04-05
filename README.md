# Refined_LGN_adaptation_analysis
### This repository includes scripts used to process neurophysiological data and plot the data. The paradigm used to collect this data is a drifting grating paradigm.  Most of the statistical analyses were performed on R (repository entitled: R_LGN_adaptation_analysis), 
### except for the analysis of the power for which a ROC analysis was performed on Matlab (not included in the manuscript)
### Most of the repository includes data processing/analysis of the monocular stimulation condition.
### One file however includes the comparison between monocular and binocular stimulation

## Regarding the location of all following files, one might need to add "OneDrive - Vanderbilt\" between 'daumail\" and "Documents\"
## Monocular stimulation condition scripts:
    #### "new_data_set.m": takes "good_single_units_data_4bmpmore.mat" located in:  "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\" and transforms it into "refined_dataset.mat" stored in "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\".
    #### This script removes some channels that were showing mistriggered trials.
    
    #### "get_clean_peaks_and_data.m": isolation of trials peak response values of the low pass filtered data and trials data from "refined_dataset.mat" for the data analysis. Trials and peak values are selected based on multiple criteria, including: at least 4 peaks per trial, no outlier in baseline activity, etc. Saves the data under the form of 5 different file types: 
		   - one file only including peak response values of the low pass filtered data of an individual single unit: Files saved in the following directory: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\individual_units\"

		   - one file gathering all peak response values of the low pass filtered data of all single units together: Saved in: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\all_data_peaks.mat"
		   - One file with all peak response locations of low-pass filtered data: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\clean_SUA_locs.mat"
		   - one file with all the trials time series data of all single units together: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\clean_origin_sup_50.mat" 
		   - one file with all low pass filtered time series data of all single units together: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\clean_SUA_sup_50.mat"

    "align_clean_data.m": retriggering all trials isolated in the script "get_clean_peaks_and_data" to peak 1 location.
   


## Comparison of adaptation between monocular and binocular stimulation: (BinocularAdaptationTrialSelection.m)
       -Data Location:- (some directories might need to be updated through inclusion of -daumail\ "OneDrive - Vanderbilt" \Documents-
       #### 1) Data used for pre-processing: (trial selection)
       #### "C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\single_units_ns6_metadata.mat"
       #### "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\refined_dataset.mat"
       #### 2) Data saved after preprocessing:
       #### "C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\NoFiltMultiContSUA.mat"
       #### "C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\all_unfiltered_data_peaks.mat"
       #### 4) Monocular adaptation pvalues (from R analysis): "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\lmer_results_orig_03032020_corrected.csv"
       #### 5) Monocular condition peak values: "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\all_data_peaks"
       #### 6) Monocular and Binocular neural data obtained through a pre-processing of the data in 1): 
       "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\all_orig_bs_zscore_trials"
       #### Interaction condition pvalues (from R analysis):
       "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\mixedmodel_pvals_anova_linearTrend.csv"


