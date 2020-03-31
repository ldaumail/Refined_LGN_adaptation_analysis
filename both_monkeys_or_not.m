
%this script intends to determine wether adaptation occurred in both monkeys or not. we look both
%at the spiking activity and the power

gendatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [gendatadir 'refined_dataset']; 
gen_data_file = load(channelfilename);

spikpvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
spikpvalfilename = [spikpvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
spikpvalues = dlmread(spikpvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

partsdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\part1_part2_norm_power.mat';
parts = load(partsdir);
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
sig95_idx = load( strcat(newdatadir,'roc_results95_stimonset_to1150ms.mat'));



layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];


filenames = cell(length(layer),4);
all_mean_data = nan(4, length(layer));
for i = 1:length(layer)
    if ~isempty(peakvals.peak_vals(i).peak)
filename = gen_data_file.new_data(i).channel_data.filename;
 filename = erase(filename, '.mat');
 filenames(i,1) = cellstr(filename);
 filenames(i,2) = cellstr(layer(i));

 mean_data = nanmean(peakvals.peak_vals(i).peak,2);
 all_mean_data(:,i) = mean_data;
  if all_mean_data(4,i) < all_mean_data(1,i) && spikpvalues(i,4) < .05
      
     filenames(i,3) = cellstr('spikadapt');
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
  if parts.parts(i,1) > parts.parts(i,2) && sig95_idx.all_sigs95(i) == 1
     filenames(i,4) = cellstr('powadapt');
  end
   end
end