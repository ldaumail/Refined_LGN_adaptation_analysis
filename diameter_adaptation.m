%This script was developped to assess the influence of grating diameter on
%adaptation
%Loic Daumail 05/06/2022
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;


unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);


cellClass = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
cellClass([1,46,55]) = [];

%select trials with at least 4 peak values, of convenient quality. Keep peak locations and trial responses: 
[peakLocs, NoFiltMultiContSUA] = peakLocsTrialSelection(unitsData, filenames);

%Store Peaks and peak-triggered trials
[peak_vals, peak_aligned_trials] = peaksAndPeakTrigResps(peakLocs, NoFiltMultiContSUA);

%Get grating diameter
grating_diameters =[];
 for i =1:length(filenames)
     if ~isempty(filenames{i})
    diameter = unique(unitsData.new_data(i).channel_data.diameter);
    grating_diameters = [grating_diameters; diameter];
     end
 end

 %assess adaptation effect in monocular condition
 dp1p4 = struct('bin',[]);
 normdp1p4 = struct('bin',[]);
 midx = struct('bin',[]);
 meanpk = struct('bin',[]);
 filenames = fieldnames(peak_vals);
 contLims = [0,0.1,0.3,0.5,0.7,1];
 clear i
 for i = 1:length(fieldnames(peak_vals))
     for  n =1:length(contLims)
         if n ==1 || n ==6
         binNb = sprintf('bin%d', n);
          if isfield(peak_vals.(filenames{i}), binNb)
         dp1p4.bin.(binNb)(i) = -mean(peak_vals.(filenames{i}).(binNb)(1,:) - peak_vals.(filenames{i}).(binNb)(4,:));
         normdp1p4.bin.(binNb)(i) = -mean(peak_vals.(filenames{i}).(binNb)(1,:) - peak_vals.(filenames{i}).(binNb)(4,:))/mean(peak_vals.(filenames{i}).(binNb)(1,:));
         midx.bin.(binNb)(i) = -2*(mean(peak_vals.(filenames{i}).(binNb)(1,:) - peak_vals.(filenames{i}).(binNb)(4,:)))/(mean(peak_vals.(filenames{i}).(binNb)(1,:) + peak_vals.(filenames{i}).(binNb)(4,:)));
         meanpk.(binNb) = mean(peak_vals.(filenames{i}).(binNb),2);
          end
         end
     end
 end
 
 
 figure()
 for  n =1:length(contLims)
     if n ==1 || n ==6
         binNb = sprintf('bin%d', n);
         if isfield(midx.bin, binNb)
             plot(grating_diameters, midx.bin.(binNb), 'o')
             hold on
         end
     end
 end
  
 
 figure();
 plot()