%% Loic Daumail 12/1/2021
%%see how variable is the stimulus phase
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;

unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);


 allphases =0;
 for i =1:71
     if ~isempty(filenames{i})
    phases = unique(unitsData.new_data(i).channel_data.phase);
    allphases = unique([allphases; phases]);
     end
 end