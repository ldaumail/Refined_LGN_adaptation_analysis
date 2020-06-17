
%following the script get_clean_peaks_and_data.m (data cleaning and
%selection pipeline)
%this script was written to isolate peaks of the origin data in order to
%perform the statistical analysis of the peaks 
%some lines were commented out and replaced in order to also isolate
%normalize peak values, averaged across trials in order to plot the
%normalized average peak values of each unit on R


newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_origin_sup_50']; 
data_file = load(channelfilename);
channelfilename = [newdatadir 'clean_SUA_sup_50']; 
filt_data_file = load(channelfilename);
locsfilename = [newdatadir 'clean_SUA_locs'];
all_locsdSUA = load(locsfilename);
gendatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [gendatadir 'refined_dataset']; 
gen_data_file = load(channelfilename);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};

xabs = -199:1300;
nyq = 500;

channum = 1: length(data_file.clean_origin_data);
mean_origin_dSUA = struct();
mean_filtered_dSUA = struct();
suas_trials = struct();
peak_vals = struct();
mean_peak_vals = struct();
mean_peaks = nan(4,length(channum));
up_dist = nan(1, length(channum),4);
max_low_dist = nan(1, length(channum));
all_locsdSUA_filtered = nan(1,length(channum),4);
filenames = cell(length(channum),2);

for i = channum  
    if ~isempty(data_file.clean_origin_data(i).unit)
    trialidx = 1:length(data_file.clean_origin_data(i).unit(1,:));
    origin_dSUA = data_file.clean_origin_data(i).unit(401:1900,:) %- mean(data_file.clean_origin_data(i).unit(401:600,:),1);
    
    %create normalized origin trials data to plot average peaks for each unit with R
    
norm_unit = nan(size(origin_dSUA));
    clear tr
    for tr =trialidx
            min_unit =min(origin_dSUA(:,tr),[],1);
            max_unit = max(origin_dSUA(:,tr),[],1);
            norm_unit(:,tr) = (origin_dSUA(:,tr)-min_unit)./(max_unit - min_unit);
    end
    
    filtered_dSUA = filt_data_file.clean_high_SUA(i).namelist;
    
  
    %determine the peak location of interest for each trial of a given single
    %unit
    all_locsdSUA_trials = all_locsdSUA.peaks_locs(i).locs;
    
    up_dist_trials = nan(4,length(trialidx));
    clear pn
    for pn = 1:4
    locs_peak = all_locsdSUA_trials(pn, :);
    up_dist_trials(pn,:)= length(xabs)- locs_peak;
    end
    %get the max distance between the peakalign and the stimulus onset
    max_low_dist_unit = max(all_locsdSUA_trials,[],'all');
    %create new matrix with the length(max(d)+max(xabs - d))
    new_dist_unit = max_low_dist_unit + max(up_dist_trials,[],'all'); 
    fp_locked_trials = nan(new_dist_unit,length(origin_dSUA(1,:)),4);
    filtered_fp_locked_trials = nan(new_dist_unit,length(filtered_dSUA(1,:)),4);
     clear n pn
     for pn =1:4
           for n =trialidx
                  lower_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+1;
                  upper_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+length(xabs);
                  
                  %origin data of the statistical analysis
                  fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = origin_dSUA(:,n);
                  %normalized data for the plotting
                 % fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = norm_unit(:,n);
                  
                  filtered_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = filtered_dSUA(:,n);
           end
 
     end
    %get the aligned data if it exists for the unit 
    suas_trials(i).aligned= fp_locked_trials;
    max_low_dist(i) = max_low_dist_unit;
    
    
    clear pn
       for pn = 1:4
           %peak data for the stats
      peak_vals(i).peak(pn,:)= max(suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,pn), [],1);
       end
       %mean peaks for the R plots 
       mean_peaks(:,i) = mean(peak_vals(i).peak,2);
    else
        %peak data for the stats
      peak_vals(i).peak = [];
      %peak data for the R plots
      mean_peaks(:,i) = nan(4,1);
       
    end
 filename = [gen_data_file.new_data(i).channel_data.filename, f{2}];
 filename = erase(filename, '.mat');
 filenames(i,1) = cellstr(filename);
 filenames(i,2) = cellstr(layer(i));
 peaks = peak_vals(i).peak;
channelfilename = [gendatadir 'su_peaks_03032020_corrected\orig_peak_values\' filename];
%save(strcat(channelfilename, '.mat'), 'peaks');
end  
 %mean_peak_vals.peak = mean_peaks;
 allfilename = [gendatadir 'su_peaks_03032020_corrected\orig_peak_values\all_units\all_raw_data_peaks'];
 save(strcat(allfilename, '.mat'), 'peak_vals');
 allfilename = [gendatadir 'su_peaks_03032020_corrected\orig_peak_values\all_units\all_raw_mean_data_peaks'];
 %save(strcat(allfilename, '.mat'), 'mean_peaks');
 savefilename = [gendatadir 'su_peaks_03032020_corrected\orig_peak_values\all_units\filenames_layers'];
 %save(strcat(savefilename, '.csv'), 'filenames');
