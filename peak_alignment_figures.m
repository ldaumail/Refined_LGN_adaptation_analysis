%this script was developped by Loic Daumail. Last updated: 5/21/2021
% The purpose of this script is to plot an example of a unit with peaks
% triggered to stimulus onset, and another version with the peak responses
% triggered to each peak, to assess the change of the profile of the mean
% response between these two approaches

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_origin_sup_50']; 
data_file = load(channelfilename);
channelfilename = [newdatadir 'clean_SUA_sup_50']; 
filt_data_file = load(channelfilename);
locsfilename = [newdatadir 'clean_SUA_locs'];
all_locsdSUA = load(locsfilename);
xabs = -199:1300;
nyq = 500;

channum = 1: length(data_file.clean_origin_data);
mean_origin_dSUA = struct();
mean_filtered_dSUA = struct();
suas_trials = struct();
up_dist = nan(1, length(channum),4);
max_low_dist = nan(1, length(channum));
all_locsdSUA_filtered = nan(1,length(channum),4);

for i = channum  
    if ~isempty(data_file.clean_origin_data(i).unit)
   i =22;
    trialidx = 1:length(data_file.clean_origin_data(i).unit(1,:));
    origin_dSUA = data_file.clean_origin_data(i).unit(401:1900,:)%-mean(data_file.clean_origin_data(i).unit(401:600,:),1);
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
        for n = 1:length(origin_dSUA(1,:))
            lower_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+1;
            upper_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+length(xabs);
            fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = origin_dSUA(:,n);
            filtered_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = filtered_dSUA(:,n);
        end
        
        eval(['mean_origin_dSUA(i).mean_peakaligned' num2str(pn) '=  nanmean(fp_locked_trials(:,:,pn),2);'])
        eval(['mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '= nanmean(filtered_fp_locked_trials(:,:,pn),2);']) % for nan - cols
        
    end
    %get the aligned data if it exists for the unit
    suas_trials(i).aligned= fp_locked_trials;
    max_low_dist(i) = max_low_dist_unit;
    end
end
    
%% Plot unaligned individual trials of example single unit 
   
   figure(); 
   plot(xabs, origin_dSUA(:,:), 'linewidth', 1)
   set(gca, 'Box', 'off')
   set(gca, 'linewidth', 2)
   xlim([-50 1100])
   filename = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\unaligned_trials_unit22');
       saveas(gcf, strcat(filename, '.png'));
       saveas(gcf, strcat(filename, '.svg')); 
 
       
  %% Example trial
  figure(); 
   plot(xabs, origin_dSUA(:,1), 'linewidth', 1)
   set(gca, 'Box', 'off')
   set(gca, 'linewidth', 2)
   xlim([-50 1100])
   filename = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\example_trial_unit22');
       saveas(gcf, strcat(filename, '.png'));
       saveas(gcf, strcat(filename, '.svg')); 
       
 %% Low pass filtered example Trial
     lpc       = 4.5; %low pass cutoff
                lWn       = lpc/nyq;
                [bwb,bwa] = butter(4,lWn,'low');
                lpdSUA      = filtfilt(bwb,bwa, origin_dSUA(:,1));
    figure(); 
   plot(xabs, lpdSUA(:,1), 'linewidth', 1)
   set(gca, 'Box', 'off')
   set(gca, 'linewidth', 2)
   xlim([-50 1100])
   filename = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\example_lowpass_trial_unit22');
       saveas(gcf, strcat(filename, '.png'));
       saveas(gcf, strcat(filename, '.svg'));         
                
%% Mean of unaligned trials


mean_unit = nanmean(origin_dSUA,2);
stdev = std(origin_dSUA,[],2);

h =figure();
plot(xabs, mean_unit, 'linewidth', 1)
hold on
h1= ciplot( mean_unit+ stdev, mean_unit-stdev,[-200:1299],[40/255 40/255 40/255],0.1);
set(h1, 'edgecolor','none')
set(h,'position',get(h,'position').*[1 1 1.15 1])
set(gca, 'Box', 'off')
set(gca, 'linewidth', 2)
xlim([-50 1100])
ylim([0 180])

filename = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\unaligned_trials_unit22_std');
saveas(gcf, strcat(filename, '.png'));
saveas(gcf, strcat(filename, '.svg')); 

%%
    %{
   figure();
   x = 1:length(fp_locked_trials(:,1,4));
   plot(x,fp_locked_trials(:,:,4))
   set(gca, 'Box', 'off')
   set(gca, 'linewidth', 2)
   %hold on
   %plot(x, mean(fp_locked_trials(:,:,4),2),'LineWidth',1, 'Color', 'black')
   
   %mean_aligned = mean(fp_locked_trials_out,2);
   %nanmean_aligned = nanmean(fp_locked_trials_out,2);
   %}
   
%% Plot of aligned individual units
     
     figure();
    % mean_unit = squeeze(nanmean(suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,:),2));
    % stdev = squeeze(std(suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,:),[],2));
     for pn =1:4
         h = subplot(1,4,pn);
         plot(-125:124, suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,pn));
         hold on
        % h1= ciplot( mean_unit(:,pn)+ 1.96*stdev(:,pn)/sqrt(14), mean_unit(:,pn)-1.96*stdev(:,pn)/sqrt(14),[-125:124],[40/255 40/255 40/255],0.1);
         %set(h1, 'edgecolor','none')
         set(h,'position',get(h,'position').*[1 1 1.15 1])
         ylim([0 275])
         xlim([-125 125])
         set(gca,'box','off')
         set(gca, 'linewidth',2)
         ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
         if pn > 1
             ax1 = gca;
             ax1.YAxis.Visible = 'off';
         end
     end
     
     
     sgtitle({'M cell trials'}, 'Interpreter', 'none')
     xlabel('Resolution (ms)')
     set(gcf,'Units','inches')
     
filename = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\aligned_trials_unit22');
%saveas(gcf, strcat(filename, '.png'));
saveas(gcf, strcat(filename, '.svg')); 
 


%% Mean Example Unit plot

%%% plot aligned units to peak of interest
%layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
%'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
%'','','M','M','M','P','M','M','M','M','P','P'};
%layer([1,46,55]) = [];
%for i = 1:length(suas_trials)
 if ~isnan(max_low_dist(i))
     figure();
     mean_unit = squeeze(nanmean(suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,:),2));
     stdev = squeeze(std(suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,:),[],2));
     for pn =1:4
         h = subplot(1,4,pn);
         plot(-125:124, mean_unit(:,pn),'linewidth', 1);
         hold on
         h1= ciplot( mean_unit(:,pn)+ stdev(:,pn), mean_unit(:,pn)-stdev(:,pn),[-125:124],[40/255 40/255 40/255],0.1);
         set(h1, 'edgecolor','none')
         set(h,'position',get(h,'position').*[1 1 1.15 1])
         ylim([0 225])
         xlim([-125 125])
         set(gca,'box','off')
         set(gca, 'linewidth',2)
         ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
         if pn > 1
             ax1 = gca;
             ax1.YAxis.Visible = 'off';
         end
     end
     
     
     sgtitle({'M cell mean activity'}, 'Interpreter', 'none')
     xlabel('Resolution (ms)')
     set(gcf,'Units','inches')
     %set(gcf,'position',[1 1 8.5 11])
     filename = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\aligned_mean_unit1_std');
     % saveas(gcf, strcat(filename, '.png'));
     saveas(gcf, strcat(filename, '.svg'));
 end

