%This script was developped to analyze the noise present in single units data

%get filenames where the data is located
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\s_potentials_analysis\analysis\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;

unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);

%get binary data
[binSpkTrials,NoFiltMultiContSUA, peakLocs] = binTrialSelection(unitsData, filenames);

allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\binary_trials_data_06022021';
save(strcat(allfilename, '.mat'), 'binSpkTrials');
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\origin_trials_data_06022021';
save(strcat(allfilename, '.mat'), 'NoFiltMultiContSUA');
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\peak_locs_data_06022021';
save(strcat(allfilename, '.mat'), 'peakLocs');

%{
%to load older data
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
peakLocs = load(strcat(datadir, 'all_locs_data_95CI_05022021'));
peakLocs = peakLocs.peakLocs;
NoFiltMultiContSUA = load(strcat(datadir,'NoFiltMultiContSUA_05022021'));
NoFiltMultiContSUA = NoFiltMultiContSUA.NoFiltMultiContSUA;
%}


%Load peakLocs, NoFiltMultiContSUA, binSpkTrials
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
peakLocs = load(strcat(datadir, 'peak_locs_data_06022021'));
peakLocs = peakLocs.peakLocs;
NoFiltMultiContSUA = load(strcat(datadir,'origin_trials_data_06022021'));
NoFiltMultiContSUA = NoFiltMultiContSUA.NoFiltMultiContSUA;
binSpkTrials = load(strcat(datadir,'binary_trials_data_06022021'));
binSpkTrials = binSpkTrials.binSpkTrials;

%align binary trials to peak locations for each peak
[wind_peak_vals, peak_aligned_trials] = binPeakTrigResps(peakLocs, NoFiltMultiContSUA, binSpkTrials);
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\peak_aligned_binary_trials_data_06022021';
save(strcat(allfilename, '.mat'), 'peak_aligned_trials');
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\binary_peak_winds30_06022021';
save(strcat(allfilename, '.mat'), 'wind_peak_vals');

%load the binary peak values data
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
wind_peak_vals = load(strcat(datadir, 'binary_peak_winds30_06022021'));
wind_peak_vals = wind_peak_vals.wind_peak_vals;


filenames = fieldnames(wind_peak_vals);
bins = [1,6];

mean_winds = nan(4, 2,length(filenames) ); %Mean spike # per window for each peak in each unit

for i = 1: length(filenames)
    filename = filenames{i};
    if length(fieldnames(wind_peak_vals.(filename))) == 2
        for b = 1:2
            binN = sprintf('bin%d',bins(b));
             %compute mean peak responses
             for p = 1:4
                 
                 count = sum(squeeze(wind_peak_vals.(filename).(binN)(p,:,:)),1);
                 mean_winds(p,b,i) = sum(count)/size(squeeze(wind_peak_vals.(filename).(binN)(p,:,:)),2);
             end

        end
    end
end

%%find a way to melt matrix in
%%p1p1p1p1....p2p2p2p2....p3p3p3p3....p4p4p4p4p4 fashion
peakvals = [squeeze(mean_winds(1,:,:))';squeeze(mean_winds(2,:,:))';squeeze(mean_winds(3,:,:))';squeeze(mean_winds(4,:,:))'];
%peakvals = reshape(peaks, [length(peaks(:,1,1))*length(peaks(1,1,:)), length(peaks(1,:,1))]);
linPeakVals = reshape(peakvals, [2*length(peakvals(:,1)),1]); 
condition = [repmat({'Monocular'},length(peakvals(:,1)),1); repmat({'Binocular'},length(peakvals(:,1)),1)];
unit = repmat(1:length(mean_winds(1,1,:)),1,8)';
peakLabel = repmat([repmat({'Pk1'}, length(mean_winds(1,1,:)),1);repmat({'Pk2'}, length(mean_winds(1,1,:)),1);repmat({'Pk3'}, length(mean_winds(1,1,:)),1);repmat({'Pk4'}, length(mean_winds(1,1,:)),1)],2,1);

%compute fano factor for each peak
%Fano factor = variance/ mean



