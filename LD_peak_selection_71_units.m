%% Code below adapted from get_clean_peaks_and_data.m
%after saving the data with new_data_set.m, we isolate the peaks and trials for the
%new analysis of the refined data
%we also save the data we want to plot (only the clean data)
%this script is the data cleaning and selection pipeline
%written by Loic Daumail 
%edited on 06-05-2020
%edited on 10/30/2022 to remove all peak selection criteria, eccept the
%minimum of 10 trials condition.

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [newdatadir 'refined_dataset']; 
data_file = load(channelfilename);

%exclude 160517, (first unit, left empty, it is a K neuron)
%Reject 180806 p1 uclust17, M cell, as doesn't seem well triggered (46)
%Reject 181207 (B) uclust22, M cell, as doesn't seem well triggered (55)
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 %% find peak locations of smoothed data, to further allow us to isolate peak values of unfiltered data in order to analyze them on R and fit a LMER
 
 channum = 1: length(data_file.new_data);
 xabs = -199:1300;
 nyq = 500;
 mean_filtered_dSUA = struct();
 
 clean_high_SUA = struct();
 clean_origin_data = struct();
 data_peaks = struct();
 peaks_locs = struct();
 
 clear i
 for i = channum
     
     filename = [data_file.new_data(i).channel_data.filename, f{2}];
     filename = erase(filename, '.mat');
     
     blankcontrast = data_file.new_data(i).channel_data.contrast ==  0 & data_file.new_data(i).channel_data.fixedc ==  0;
     highcontrast = data_file.new_data(i).channel_data.contrast >=  0.5 & data_file.new_data(i).channel_data.fixedc ==  0;
     
     
     trialidx = 1:length(data_file.new_data(i).channel_data.sdftr_chan(1,:));
     raw_bs = nan(length(xabs), length(trialidx));
     filtered_dSUA = nan(length(xabs), length(trialidx));
     origin_data = nan(length(xabs)+401, length(trialidx));
     all_norm_lpdSUA= nan(length(xabs),length(trialidx));
     
     
%      powerstim = nan(length(trialidx),1025);
%      freqstim = nan(length(trialidx),1025);
%      fourhzpowerstim =nan(length(trialidx),1);
%      bsl = nan(1, length(trialidx));
%      mean_wnd1 = nan(1,length(trialidx));
     
     all_pks = nan(4,length(data_file.new_data(i).channel_data.sdftr_chan(1,highcontrast)));
     
     for tridx = trialidx
%          
         all_data = data_file.new_data(i).channel_data.sdftr_chan(401:1900,tridx);
         origin_data(:,tridx) = data_file.new_data(i).channel_data.sdftr_chan(:,tridx);
         raw_bs(:,tridx) = all_data(1:end)- mean(all_data(1:200));
         
         lpc       = 4.5; %low pass cutoff
         lWn       = lpc/nyq;
         [bwb,bwa] = butter(4,lWn,'low');
         lpdSUA      = filtfilt(bwb,bwa, raw_bs(:,tridx));
         
%          
         filtered_dSUA(:,tridx) = lpdSUA;
%          %all_norm_lpdSUA(:,tridx) = (lpdSUA - min(lpdSUA))/(max(lpdSUA)- min(lpdSUA));
%          mean_wnd1(tridx) = mean(lpdSUA(201:480)); %mean spiking response of a trial over 280ms following stim onset
%          
%          %%% power
%          
%          
%          [powerstim(tridx,:), freqstim(tridx,:)] = calcFFT(all_data(200:1350)); %fourrier transform
%          
%          %find the index of the frequency vector closest to 4hz and point to the
%          %power value of this index for every trial, and store the value in
%          %fourhzpower
%          [val,index] = min(abs(4-freqstim(tridx,:))); %index of value closest to 4Hz
%          fourhzpowerstim(tridx,1) = powerstim(tridx,index); %get power at that index, assuming this is 4Hz
%          
     end
     
        %%%%%%%%%%% %reject trials below Mean + 1.96*STD in the blank condition %%%%%%
    %power related variables
%      power0 = fourhzpowerstim(blankcontrast);
%      power5 = fourhzpowerstim(highcontrast);
%      
     %spiking activity related variables
%      mean_wnd1_5 =mean_wnd1(highcontrast);
     filtered_dSUA_high = filtered_dSUA(:, highcontrast);
     origin_data_high = origin_data(:, highcontrast);
     %first peak location related variables
%      sua_bsl =  mean(filtered_dSUA_high(1:200,:),1);
     
%      for tr = 1:length(power5)
%          if mean_wnd1_5(tr) > mean(sua_bsl)+1.96*std(sua_bsl)  && power5(tr) > mean(power0)+1.96*std(power0) %
%              
%              filtered_dSUA_high(:,tr) = filtered_dSUA_high(:,tr);
%              origin_data_high(:,tr) = origin_data_high(:,tr);
%          else
%              
%              filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
%              origin_data_high(:,tr) =  nan(length(origin_data_high(:,tr)),1);
%          end
%      end
     
     %determine the first peak location for each trial of a given single
     %unit
     all_locsdSUA_trials = nan(6,length(filtered_dSUA_high(1,:)));
     clear trial
     for trial = 1:length(filtered_dSUA_high(1,:))
         
         for ln = 1:550
             if filtered_dSUA_high(200+ln,trial) < filtered_dSUA_high(200+ln+1,trial) && ~all(isnan(filtered_dSUA_high(:,trial)))
                 
                 locsdSUA_trial_struct = findpeaks_Loic(filtered_dSUA_high(200+ln:1499,trial));
                 locsdSUA_trial = locsdSUA_trial_struct.loc;
                 %if peak1 is too small, peak2 becomes peak1
               %  if filtered_dSUA_high(locsdSUA_trial(1)+200+ln,trial) >= 0.4*filtered_dSUA_high(locsdSUA_trial(2)+200+ln)
                     %store first peak location
                     all_locsdSUA_trials(1:length(locsdSUA_trial),trial) = locsdSUA_trial(1:end)+200+ln;
                 %else
                  %   all_locsdSUA_trials(1:length(locsdSUA_trial(2:end)),trial) = locsdSUA_trial(2:end)+200+ln;
                     
               %  end
                 
                 break
             end
         end
         
         if nnz(~isnan(all_locsdSUA_trials(:,trial))) >= 4 && ~all(isnan(all_locsdSUA_trials(:,trial)))
             %adjust location to the first data point of lpsu (+ln),
             
             all_pks(:,trial) = filtered_dSUA_high(all_locsdSUA_trials(1:4,trial), trial);
             filtered_dSUA_high(:,trial) = filtered_dSUA_high(:,trial);
             all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
             origin_data_high(:,trial) = origin_data_high(:,trial);
         else
             filtered_dSUA_high(:,trial) = nan(length(filtered_dSUA_high(:,trial)),1);
             all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
             origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);
         end
         
         if ~all(isnan(all_locsdSUA_trials(:,trial))) && (all_locsdSUA_trials(4,trial) ~= 1500)
             %adjust location to the first data point of lpsu (+ln),
             
             all_pks(:,trial) = filtered_dSUA_high(all_locsdSUA_trials(1:4,trial), trial);
             filtered_dSUA_high(:,trial) = filtered_dSUA_high(:,trial);
             all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
             origin_data_high(:,trial) = origin_data_high(:,trial);
         else
             all_pks(:,trial) = nan(length(all_pks(:,trial)),1);
             filtered_dSUA_high(:,trial) = nan(length(filtered_dSUA_high(:,trial)),1);
             all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
             origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);
             
         end
     end
     %{
    figure(); plot(-199:1300, filtered_dSUA_high(1:1500,1))
    hold on
    plot(all_locsdSUA_trials(1:4,1)-200, all_pks(:,1))
    set(gca,'box','off')
     %}
     %%% reject outlier peaks and the corresponding trials in
     %%% filtered_dSUA_high
     
     
     %reject if there is a peak 1 outlier, if the max peak value in the
     %baseline is an outlier
     
     % First find peaks before stimulus onset
     
%      bsl_peaks = nan(1, length(filtered_dSUA_high(1,:)));
%      clear tr
%      for tr = 1:length(filtered_dSUA_high(1,:))
%          
%          for loc = 1:200
%              if filtered_dSUA_high(loc,tr) < filtered_dSUA_high(loc+1,tr) && ~all(isnan(filtered_dSUA_high(:,tr)))
%                  if length(filtered_dSUA_high(loc:200,tr)) >= 3
%                      if ~isempty(findpeaks_Loic(filtered_dSUA_high(loc:200,tr)))
%                          bsl_peak_locs_struct = findpeaks_Loic(filtered_dSUA_high(loc:200,tr));
%                          bsl_peak_locs = bsl_peak_locs_struct.loc;
%                          
%                          bsl_peaks(1,tr) = max(filtered_dSUA_high(bsl_peak_locs+loc,tr));
%                      else
%                          bsl_peaks(1,tr) = NaN;   
%                      end
%                  end
%                  break
%              end
%          end
%      end
%      
%      out_bsl_peaks = isoutlier(bsl_peaks);
     
%      p1outliers = isoutlier(all_pks(1,:));
%      clear tr
%      for tr = 1:length(filtered_dSUA_high(1,:))
%          %exclude trials
%          if p1outliers(tr) == 0 && ~all(isnan(all_pks(:,tr))) && out_bsl_peaks(tr) ==0
%              
%              filtered_dSUA_high(:,tr) = filtered_dSUA_high(:, tr);
%              all_pks(:, tr) = all_pks(:,tr);
%              all_locsdSUA_trials(:,tr) = all_locsdSUA_trials(:,tr);
%              origin_data_high(:,tr) = origin_data_high(:, tr);
%              
%          else
%              filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
%              all_pks(:,tr) = nan(length(all_pks(:,tr)),1);
%              all_locsdSUA_trials(:,tr) = nan(size(all_locsdSUA_trials(:,tr)));
%              origin_data_high(:,tr) = nan(length(origin_data_high(:,tr)),1);
%          end
%      end
     filtered_dSUA_high = filtered_dSUA_high(:,~all(isnan(filtered_dSUA_high))); % for nan - cols
     all_locsdSUA_trials =  all_locsdSUA_trials(:,~all(isnan(all_locsdSUA_trials)));
     all_pks = all_pks(:, ~all(isnan(all_pks)));
     origin_data_high = origin_data_high(:,~all(isnan(origin_data_high)));
     
     
     if length(filtered_dSUA_high(1,:)) >=10
         clean_high_SUA(i).namelist =  filtered_dSUA_high;
         clean_origin_data(i).unit = origin_data_high;
         peaks_locs(i).locs = all_locsdSUA_trials;
     elseif length(filtered_dSUA_high(1,:)) <10
         %all_pks(:,:) = nan(length(all_pks(:,1)),length(all_pks(1,:)));
         all_pks(:,:) = [];
         clean_high_SUA(i).namelist =  [];
         clean_origin_data(i).unit = [];
         peaks_locs(i).locs = [];
         
     end
     
     data_peaks(i).namelist = all_pks(:,~all(isnan(all_pks)));
     all_pks = all_pks(:,~all(isnan(all_pks)));
     %channelfilename = [newdatadir 'su_peaks_03032020_corrected\individual_units\' filename];
     %save(strcat(channelfilename, '.mat'), 'all_pks');
 end
 %allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\all_data_peaks'];
 %save(strcat(allfilename, '.mat'), 'data_peaks');
 
%  allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\clean_SUA_sup_50_03052020'];
%  save(strcat(allfilename, '.mat'), 'clean_high_SUA');
%  allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\clean_SUA_locs_03052020'];
%  save(strcat(allfilename, '.mat'), 'peaks_locs');
 allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\clean_origin_sup_50_10302022'];
 save(strcat(allfilename, '.mat'), 'clean_origin_data');
 
 allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\filt_SUA_sup_50_10302022'];
 save(strcat(allfilename, '.mat'), 'clean_high_SUA');
 allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\filt_SUA_locs_10302022'];
 save(strcat(allfilename, '.mat'), 'peaks_locs');
 
 
 %count number of remaining units after preproc
   cnt =0;
   for i =1:length(data_peaks)
       if ~isempty(data_peaks(i).namelist)
           cnt = cnt+1;
       end
        
   end
   
 %% code below adapted from get_origin_peaks.m
 
%following the script get_clean_peaks_and_data.m (data cleaning and
%selection pipeline)
%this script was written to isolate peaks of the origin data in order to
%perform the statistical analysis of the peaks 
%some lines were commented out and replaced in order to also isolate
%normalized peak values, averaged across trials in order to plot the
%normalized average peak values of each unit on R
%Written by Loic Daumail, last edited on 6/29/2020, then 10/23/2022

%newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_origin_sup_50_10302022']; 
data_file = load(channelfilename);
channelfilename = [newdatadir 'filt_SUA_sup_50_10302022']; 
filt_data_file = load(channelfilename);
locsfilename = [newdatadir 'filt_SUA_locs_10302022'];
all_locsdSUA = load(locsfilename);
gendatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
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
    origin_dSUA = data_file.clean_origin_data(i).unit(401:1900,:); %- mean(data_file.clean_origin_data(i).unit(401:600,:),1);
    
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
channelfilename = [gendatadir 'su_peaks_10302022\orig_peak_values\' filename];
%save(strcat(channelfilename, '.mat'), 'peaks');
end  
 mean_peak_vals.peak = mean_peaks;
 allfilename = [gendatadir 'su_peaks_10302022\orig_peak_values\all_units\all_raw_data_peaks'];
 save(strcat(allfilename, '.mat'), 'peak_vals');
 allfilename = [gendatadir 'su_peaks_10302022\orig_peak_values\all_units\all_raw_mean_data_peaks'];
 save(strcat(allfilename, '.mat'), 'mean_peaks');
 % Convert cell to a table and use first row as variable names
T = cell2table(filenames);
 savefilename = [gendatadir 'su_peaks_10302022\orig_peak_values\all_units\filenames_layers'];
% Write the table to a CSV file
writetable(T,strcat(savefilename, '.csv'))
 
 save(strcat(savefilename, '.csv'), 'filenames');
