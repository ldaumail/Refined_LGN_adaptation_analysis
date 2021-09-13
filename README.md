# Refined_LGN_adaptation_analysis
### This repository includes scripts used to process neurophysiological data and plot the data. The paradigm used to collect this data is a drifting grating paradigm.  Most of the statistical analyses were performed on R (repository entitled: R_LGN_adaptation_analysis), 
### except for the analysis of the power for which a ROC analysis was performed on Matlab (not included in the manuscript)
### Most of the repository includes data processing/analysis of the monocular stimulation condition.
### One file however includes the comparison between monocular and binocular stimulation

## Regarding the location of all following files, one might need to add "OneDrive - Vanderbilt\" between 'daumail\" and "Documents\"
## Monocular stimulation condition scripts:
### Preprocessing of spiking activity data
    #### "new_data_set.m": takes "good_single_units_data_4bmpmore.mat" located in:  "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\" and transforms it into "refined_dataset.mat" stored in "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\".
    #### This script removes some channels that were showing mistriggered trials.
    
    #### "get_clean_peaks_and_data.m": isolation of trials peak response values of the low pass filtered data and trials data from "refined_dataset.mat" for the next steps of data analysis. Trials and peak values are selected based on multiple criteria, including: at least 4 peaks per trial, no outlier in baseline activity, etc. Saves the data under the form of 5 different file types: 
		   - one file only including peak response values of the low pass filtered data of an individual single unit: Files saved in the following directory: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\individual_units\"
		   - one file gathering all peak response values of the low pass filtered data of all single units together: Saved in: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\all_data_peaks.mat"
		   - One file with all peak response locations of low-pass filtered data: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\clean_SUA_locs.mat"
		   - one file with all the trials time series data of all single units together: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\clean_origin_sup_50.mat" 
		   - one file with all low pass filtered time series data of all single units together: "C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\clean_SUA_sup_50.mat"
     #### (Not very useful:) "align_clean_data.m": retriggering all trials isolated in the script "get_clean_peaks_and_data.m" to peak 1 location isolated in "get_clean_peaks_and_data.m". The former version of this script "former_align_clean_data.m" accomplishes the same goal, but preprocesses the data in a less sofisticated way than "get_clean_peaks_and_data.m", thus, it is better to use "get_clean_peaks_and_data.m"+"align_clean_data.m" instead.
     #### "get_origin_peaks.m": following "get_clean_peaks_and_data.m", allows isolation of peak response values of origin, using peak locations of low pass filtered data (obtained in "get_clean_peaks_and_data.m") + time windows arround those locations to obtain the max peak value (Low pass peak time location+-125ms)
                   - Data saved at: "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\" where individual single unit peak values are stored in separate files in "\su_peaks_03032020_corrected\orig_peak_values\", all peak values saved in "\su_peaks_03032020_corrected\orig_peak_values\all_units\all_raw_data_peaks.m", mean peak values used for R plots: "\su_peaks_03032020_corrected\orig_peak_values\all_units\all_raw_mean_data_peaks.m", corresponding cell classes and single unit filenames: "\su_peaks_03032020_corrected\orig_peak_values\all_units\filenames_layers.m"
     #### "get_origin_troughs.m": similar to previous script for response trough values.
     
### Statistical analysis results of monocular adaptation of spiking activity data:
     #### "lmer_results_anal.m": follows "get_origin_peaks.m" and "individual_channels_analysis_origpeaks_KRpvalcorr.R" of the statistical analysis on R.This script ("lmer_results_anal.m"), allows to get proportion values of single units that show significant adapatation effects.
     #### "cohensd.m": effect sizes of both spiking activity data and analysis of the F1 component (power at 4Hz).
     #### "both_monkeys_or_not.m": determines which single unit showed adaptation, according to the session, allows to see which monkey it was recorded from. The goal was to see if adaptation is seen in both monkeys. Data saved in: "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\selected_units_filenames.m", "C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\selected_units_sessions"
     
     
### Analysis of F1 component (power, spectrograms, 4Hz response):
     #### "clean_data_ROCAnalysis.m": computes spectrogram, isolates response at 4Hz frequency, performs ROC analysis between two time windows (beginning of the trial versus end of the trial power responses)
### Plotting scripts of adaptation analysis in monocular stimulation condition:
     #### "paper_plots.m": summarizes most of the plots of adaptation manuscript, especially for monocular stimulation condition.
     #### "peak_alignment_figures.m": makes figures that show how peak responses were retriggered for plotting
     #### other plotting scripts (most of them are in "paper_plots.m" in a better version: "plot_clean_origin_data.m", "plot_clean_origin_peaks.m", "plot_clean_origin_troughs.m", "plot_spike_waveform.m", "plot_units_peak_aligned.m".

### Spike waveform shapes: 
    ##### "NarrowVsBroadWf.m"

### Multiple contrasts analysis:
     #### "multiple_contrasts_peak_isolation.m": allows to isolate monocular responses under various stimulation contrasts and to plot contrast response curves of LGN neurons.

### S-potential isolation attempt scripts: (another repository is dedicated to this analysis, most of s-potential analysis is in s-potential repository, not here, these are just very initial scripts, e.g. not so useful)
    #### "get_ns6_files.m": get the files containing data sampled at 30kHz and store them in specific folder architecture
    #### "s_potentials_isolation.m": not interesting, but preprocesses LFP and MUA 
### Microsaccades analysis: (another repository is dedicated to this analysis)
    #### "microsaccades_analysis.m": very initial script to look at electroocculograms. The main microsaccade isolation analysis is in the corresponding repository.
    

## Comparison of adaptation between monocular and binocular stimulation:
1) "BinocularAdaptationTrialSelection.m"
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
       
2) "adaptation_index.m": 
       - calculates an adaptation index, in order to compare the binocular and monocular adaptation conditions.
       - plots the data using gramm toolbox (https://github.com/piermorel/gramm)
       - Perform ROC analysis to test on difference of distributions
       - Performs other tests for difference between distributions (Qi square goodness of fit test...)
       
3) "adaptation_and_binocular_interaction_plots.m"
       - Allows to plot all results (horizontal histograms comparing each peak between monocular and binocular condition, jitter..)
 
 ## Variability quenching:
 1) Please refer to "fano_factor.m", the whole analysis pipeline is coded in this script, and refers to functions needed to run it, such as "peakTrigFano"
 



