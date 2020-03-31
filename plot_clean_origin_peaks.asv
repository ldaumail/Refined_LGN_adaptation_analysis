newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
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
    trialidx = 1:length(data_file.clean_origin_data(i).unit(1,:));
    origin_dSUA = data_file.clean_origin_data(i).unit(401:1900,:)- mean(data_file.clean_origin_data(i).unit(401:600,:),1);
    filtered_dSUA = filt_data_file.clean_high_SUA(i).namelist;
    
    %{
       for  tr = trialidx
           lpc       = 4.5; %low pass cutoff
    lWn       = lpc/nyq;
    [bwb,bwa] = butter(4,lWn,'low');
    origin_dSUA(:,tr) = filtfilt(bwb,bwa, origin_dSUA(:,tr));
       end
    %}
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
    %{
        %compute the mean single unit activity if more than 10 trials
    filtered_mean = nanmean(filtered_fp_locked_trials(:,:,pn),2);
    lpc       = 4.5; %low pass cutoff
    lWn       = lpc/nyq;
    [bwb,bwa] = butter(4,lWn,'low');
    lpdSUA      = filtfilt(bwb,bwa, filtered_mean(~isnan(filtered_mean)));
           %}
     eval(['mean_origin_dSUA(i).mean_peakaligned' num2str(pn) '=  nanmean(fp_locked_trials(:,:,pn),2);']) 
     eval(['mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '= nanmean(filtered_fp_locked_trials(:,:,pn),2);']) % for nan - cols
    
     end
    %get the aligned data if it exists for the unit 
    suas_trials(i).aligned= fp_locked_trials;
    max_low_dist(i) = max_low_dist_unit;
     
    %{
    
   
   figure();
   x = 1:length(fp_locked_trials(:,1,4));
   plot(x,fp_locked_trials(:,:,4))
   hold on
   plot(x, mean(fp_locked_trials(:,:,4),2),'LineWidth',1, 'Color', 'black')
   
   %mean_aligned = mean(fp_locked_trials_out,2);
   %nanmean_aligned = nanmean(fp_locked_trials_out,2);
   %}
   
 
    end
end  


%%% plot aligned units to peak of interest
layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
for i = 1:length(suas_trials)
     if ~isnan(max_low_dist(i))
        figure();
        %plot(1:250, mean(suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,1),2))
        plot(1:250, suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,1))
        hold on
        plot(251:500, suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,2))
        hold on
        plot(501:750, suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,3))
        hold on
        plot(751:1000, suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,4))
      set(gca, 'linewidth',2)
      set(gca,'box','off')

    title({sprintf('%s | %s', num2str(i), char(layer(i))), 'responses, p<0.05, associated to adaptation pvalues'}, 'Interpreter', 'none')
    xlabel('Resolution (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
        filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020\plots\',strcat( sprintf('4_5hzfilt_trials_bscorr_origin_data_aligned_unit_%d', i)));
        saveas(gcf, strcat(filename, '.png'));
     end

end
    
nnz(~isnan(max_low_dist))

%% plot mean cell class activity with data aligned to peak of interest
layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

layer_idx = find(strcmp(layer, 'K'));
%%store all trials of a given peak, in a matrix, across all units
aligned_trials = nan(1200, length(layer_idx));
clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:4
           aligned_trials(250*(pn-1)+1:250*pn,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-1-124:max_low_dist(layer_idx(i))-1+125,:,pn),2);
       end
    end
end
aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

      figure();
      mean_origin =nanmean(aligned_trials,2);
      plot(1:250, mean_origin(1:250),'LineWidth',1, 'Color',[165/255 42/255 42/255])
      hold on
      plot(301:550, mean_origin(251:500),'LineWidth',1, 'Color',[165/255 42/255 42/255])
      hold on
      plot(601:850, mean_origin(501:750),'LineWidth',1, 'Color',[165/255 42/255 42/255])
      hold on
      plot(901:1150, mean_origin(751:1000),'LineWidth',1, 'Color',[165/255 42/255 42/255])
      hold on
      ci_low = mean_origin - 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(aligned_trials(1,:)));
      plot(1:250, ci_low(1:250),':', 'LineWidth',.7,'Color', [165/255 42/255 42/255])
      hold on
      plot(301:550, ci_low(251:500),':', 'LineWidth',.7,'Color', [165/255 42/255 42/255])
      hold on
      plot(601:850, ci_low(501:750),':', 'LineWidth',.7,'Color', [165/255 42/255 42/255])
      hold on
       plot(901:1150, ci_low(751:1000),':', 'LineWidth',.7,'Color', [165/255 42/255 42/255])
      hold on
      ci_high = mean_origin + 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(aligned_trials(1,:)));
      
      plot(1:250, ci_high(1:250),':', 'LineWidth',.7,'Color', [165/255 42/255 42/255])
      hold on
      plot(301:550, ci_high(251:500),':', 'LineWidth',.7,'Color', [165/255 42/255 42/255])
      hold on
      plot(601:850, ci_high(501:750),':', 'LineWidth',.7,'Color', [165/255 42/255 42/255])
      hold on
      plot(901:1150, ci_high(751:1000),':', 'LineWidth',.7,'Color', [165/255 42/255 42/255])
   
 %blue [49/255 130/255 189/255]
 %black [222/255 45/255 38/255]
 %orange [253/255 174/255 107/255]
 %brown [165/255 42/255 42/255]
 %red [215/255 25/255 28/255]
 %pink = [229/255, 49/255, 90/255])
 
      set(gca, 'linewidth',2)
      set(gca,'box','off')

    title({'K mean responses'}, 'Interpreter', 'none', 'FontSize', 20)
    xlabel('\fontsize{14}Resolution (ms)')
    ylh = ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
   
   set(gcf,'Units','inches') 
   %set(gcf,'position',[1 1 8.5 11])
   set(gcf,'position',[1 1 15 11])
   
   filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020\plots\',strcat('spaced_95ci_nanmean_bscorr_origin_data_aligned_Kcells_brown'));
   saveas(gcf, strcat(filename, '.png'));
   
   
   
   % h = subplot(1,4,1). then, set(h,'position',get(h,'position').*[1 1 1.25 1]); 
   %then for sublots 2,3,4, ... set(gca,'ycolor','w','yticklabel',[]); also for every subplot set ylims to the same values
   
   %% plot mean cell class activity with data aligned to peak of interest using subplot
layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

layer_idx = find(strcmp(layer, 'M'));
%%store all trials of a given peak, in a matrix, across all units
aligned_trials = nan(1204, length(layer_idx));
norm_aligned_trials = nan(1204, length(layer_idx));
clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:4
           aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
           
       end
       norm_aligned_trials(:,i) = (aligned_trials(:,i) - min(aligned_trials(:,i)))/(max(aligned_trials(:,i))-min(aligned_trials(:,i)));
    end
end
aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

      figure();
      mean_origin =nanmean(norm_aligned_trials,2);
      ci_low = mean_origin - 1.96*std(norm_aligned_trials,0,2, 'omitnan')./sqrt(length(aligned_trials(1,:)));
      ci_high = mean_origin + 1.96*std(norm_aligned_trials,0,2, 'omitnan')./sqrt(length(aligned_trials(1,:)));
      
      for pn = 1:4
      h =subplot(1,4,pn);
      plot(-125:125, mean_origin(250*(pn-1)+1:250*pn+1),'LineWidth',1, 'Color',[215/255 25/255 28/255])
      hold on
      plot(-125:125, ci_low(250*(pn-1)+1:250*pn+1),':', 'LineWidth',.7,'Color', [215/255 25/255 28/255])
      hold on
      plot(-125:125, ci_high(250*(pn-1)+1:250*pn+1),':', 'LineWidth',.7,'Color', [215/255 25/255 28/255])
     set(h,'position',get(h,'position').*[1 1 1.15 1])
      % P ylim([-20 85])
      % M ylim([-35 135])
      % K ylim([-20 200])
      ylim([0 1])
      xlim([-125 125])
      set(gca,'box','off')
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
      if pn >1 
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off';   
      end
      
      end
 
  
 %red [215/255 25/255 28/255] 
 %orange [253/255 174/255 97/255]
 %blue [44/255 123/255 182/255]

 %brown [165/255 42/255 42/255]
    currfig = gcf;
    
    title(currfig.Children(end),{'M mean responses'}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    ylh = ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020\plots\',strcat('subplots_norm_95ci_nanmean_bscorr_origin_data_aligned_Mcells_red'));
   saveas(gcf, strcat(filename, '.png'));
   
   %% plot with significant adapting units
   
pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
cellclass = [ 'M', 'P', 'K'];
red = [215/255 25/255 28/255];
 orange = [253/255 174/255 97/255];
 blue =[44/255 123/255 182/255];
colors = [red; orange; blue];
% M ylim([0.05 .75])
% P ylim([.1 .8])
% K ylim([-0.01 1])
normylims = [[0.05 .75];[.1 .8];[-0.01 1]];
 % P ylim([-20 85])
 % M ylim([-35 135])
 % K ylim([-20 200])
ylims = [[-30 135];[-20 85];[-20 200]];

for nc = 1:3
clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
%%store all trials of a given peak, in a matrix, across all units
%aligned_trials = nan(1204, length(layer_idx));
clear aligned_trials
aligned_trials = nan(1204, length(layer_idx));

%norm_aligned_trials = nan(1204, length(layer_idx));
%norm_aligned_trials = nan(1204, length(layer_idx));

clear i 
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:4
          % aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
           aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-1-124:max_low_dist(layer_idx(i))+125,:,pn),2);
       end
       %normalizing with max and min of each unit
       %norm_aligned_trials(:,i) = (aligned_trials(:,i) - min(aligned_trials(:,i)))/(max(aligned_trials(:,i))-min(aligned_trials(:,i)));
     
    end
end
  %normalizing with max and min across units
       maximum = max(aligned_trials,[],'all');
       minimum = min(aligned_trials,[],'all');
       norm_aligned_trials = (aligned_trials - minimum)/(maximum-minimum);

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))

 clear sig_su mean_sig_su
  cnt = 0;
 all_mean_data = nan(4, length(layer_idx));
 sig_adapsu = nan(length(aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
 mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cnt= cnt+1;
      sig_adapsu(:,cnt) = aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
      end
  end
  mean_sig_su = nanmean(sig_adapsu,2);

  
 
      figure();
      mean_origin =nanmean(aligned_trials,2);
      ci_low = mean_origin - 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(aligned_trials(1,:)));
      ci_high = mean_origin + 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(aligned_trials(1,:)));
      
      for pn = 1:4
      h =subplot(1,4,pn);
      plot(-125:125, mean_origin(250*(pn-1)+1:250*pn+1),'LineWidth',1, 'Color',colors(nc,:) )
      hold on
      plot(-125:125, ci_low(250*(pn-1)+1:250*pn+1),':', 'LineWidth',.7,'Color',colors(nc,:) )
      hold on
      plot(-125:125, ci_high(250*(pn-1)+1:250*pn+1),':', 'LineWidth',.7,'Color',  colors(nc,:))
      hold on
      plot(-125:125, mean_sig_su(250*(pn-1)+1:250*pn+1),  'LineWidth',1, 'Color',[141/255 140/255 140/255] )
     set(h,'position',get(h,'position').*[1 1 1.15 1])
     
        ylim([-35 200])
      
      
       % K norm ylim([-0.02 1.1])
      xlim([-125 125])
      set(gca,'box','off')
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
      if pn >1 
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off';   
      end
      
      end

 %brown [165/255 42/255 42/255]
    currfig = gcf;
    
    title(currfig.Children(end),{sprintf('%s mean responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    ylh = ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])

   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('sigadapt_subplots_common_95ci_nanmean_bscorr_origin_data_aligned_%s', cellclass(nc),'cells')));
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
end

   
   %% plot with significant adapting units and overall mean in two separate sublots, stacking all peaks together
   
 pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];
 %red [215/255 25/255 28/255] 
 %orange [253/255 174/255 97/255] %brown [121/255 96/255 76/255]; 
 %blue [44/255 123/255 182/255]
brown = [121/255 96/255 76/255]; 
red = [215/255 25/255 28/255];
orange = [253/255 174/255 97/255];
blue =[44/255 123/255 182/255];
colors = [red; orange; blue];
% M ylim([0.05 .75])
% P ylim([.1 .8])
% K ylim([-0.01 1])
normylims = [[0.05 .75];[.1 .8];[-0.01 1]];
 % P ylim([-20 85])
 % M ylim([-35 135])
 % K ylim([-20 200])
ylims = [[-30 135];[-20 85];[-20 200]];

for nc = 1:3
 clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
%%store all trials of a given peak, in a matrix, across all units
clear aligned_trials
aligned_trials = nan(1204, length(layer_idx));

clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:4
          % aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
           aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-1-124:max_low_dist(layer_idx(i))+125,:,pn),2);
       end
       %normalizing with max and min of each unit
       %norm_aligned_trials(:,i) = (aligned_trials(:,i) - min(aligned_trials(:,i)))/(max(aligned_trials(:,i))-min(aligned_trials(:,i)));
     
    end
end
  %normalizing with max and min across units
       maximum = max(aligned_trials,[],'all');
       minimum = min(aligned_trials,[],'all');
       norm_aligned_trials = (aligned_trials - minimum)/(maximum-minimum);

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))

 clear sig_su mean_sig_su
  cnt = 0;
 all_mean_data = nan(4, length(layer_idx));
 sig_adapsu = nan(length(norm_aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
 mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cnt= cnt+1;
      sig_adapsu(:,cnt) = aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
      end
  end
  mean_sig_su = nanmean(sig_adapsu,2);



     figure();
     mean_origin =nanmean(aligned_trials,2);
 
     col1 = [121/255 96/255 76/255]; 
     col2 = colors(nc,:);
     
      R = linspace(col1(1), col2(1),4);
      G = linspace(col1(2), col2(2),4);
      B = linspace(col1(3), col2(3),4);
      cmap = [R' G' B'];
      colormap(cmap); 
      
      h = subplot(1,2,1);
      
      nlines = 4;
      for nl = 1:nlines
      %cmap = jet(4); 
   
      plot(-125:125, mean_origin(250*(nl-1)+1:250*nl+1), 'LineWidth',1,'color',cmap(nl,:));
      hold on
      end
      set(h,'position',get(h,'position').*[1 1 1.15 1])
      ylim(ylims(nc,:))
      xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('Overall Mean')
      
      h1 = subplot(1,2,2);
      for nl = 1:nlines 
      plot(-125:125, mean_sig_su(250*(nl-1)+1:250*nl+1),  'LineWidth',1, 'Color',cmap(nl,:) )
      hold on
      end
      set(h1,'position',get(h1,'position').*[1 1 1.15 1])
      ylim(ylims(nc,:))
      xlim([-125 125])
      title('Adapting units mean')
      set(gca,'box','off')
      legend('peak 1', 'peak 2', 'peak 3', 'peak 4')
    
      
       % K norm ylim([-0.02 1.1])
      
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
     % if pn >1 
      %ax1 = gca;                   
      %ax1.YAxis.Visible = 'off';   
     % end

    currfig = gcf;
    
    title(currfig.Children(end),{sprintf('%s mean responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('sigadapt_stacked_common_nanmean_bscorr_origin_data_aligned_%s',cellclass(nc),'cells_cmap_corrected_scaled')));
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
end 


   %% plot with significant adapting units, non significant units and overall mean in three separate subplots, stacking all peaks together
   
 pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];

%{
brown = [121/255 96/255 76/255]; 
red = [215/255 25/255 28/255];
orange = [253/255 174/255 97/255];
blue =[44/255 123/255 182/255];
colors = [red; orange; blue];
%}
%maps
nlines = 4;
cmaps = struct();
cmaps(1).map =cbrewer2('Reds', nlines);
cmaps(2).map =cbrewer2('Oranges', nlines);
cmaps(3).map =cbrewer2('Blues', nlines);

% M ylim([0.05 .75])
% P ylim([.1 .8])
% K ylim([-0.01 1])
normylims = [[0.05 .75];[.1 .8];[-0.01 1]];
 % P ylim([-20 85])
 % M ylim([-30 135])
 % K ylim([-20 200])
ylims = [[-30 135];[-20 85];[-20 200]];

for nc = 1:3
 clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
%%store all trials of a given peak, in a matrix, across all units
clear aligned_trials
aligned_trials = nan(1204, length(layer_idx));

clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:4
          % aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
           aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-1-124:max_low_dist(layer_idx(i))+125,:,pn),2);
       end
       %normalizing with max and min of each unit
       %norm_aligned_trials(:,i) = (aligned_trials(:,i) - min(aligned_trials(:,i)))/(max(aligned_trials(:,i))-min(aligned_trials(:,i)));
     
    end
end
  %normalizing with max and min across units
       maximum = max(aligned_trials,[],'all');
       minimum = min(aligned_trials,[],'all');
       norm_aligned_trials = (aligned_trials - minimum)/(maximum-minimum);

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))

 clear sig_su non_sig_su mean_sig_su
 cntsig = 0;
 cntnotsig = 0;
 all_mean_data = nan(4, length(layer_idx));
 sig_adapsu = nan(length(norm_aligned_trials(:,1)),length(layer_idx));
 non_sig_su =nan(length(norm_aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
 mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cntsig= cntsig+1;
      sig_adapsu(:,cntsig) = norm_aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
  if  pvalues(layer_idx(nunit),4) > .05
      
      cntnotsig= cntnotsig+1;
      non_sig_su(:,cntnotsig) = norm_aligned_trials(:,nunit);
      
  end
      end
  end
  mean_sig_su = nanmean(sig_adapsu,2);
  mean_not_sig_su = nanmean(non_sig_su,2);
  mean_origin =nanmean(norm_aligned_trials,2);

     
     figure();
     %{
     col1 = [121/255 96/255 76/255]; 
     col2 = colors(nc,:);
     
      R = linspace(col1(1), col2(1),4);
      G = linspace(col1(2), col2(2),4);
      B = linspace(col1(3), col2(3),4);
      cmap = [R' G' B'];
      colormap(cmap); 
     %}
       
       cmap = flip(cmaps(nc).map) ;
      % axes('ColorOrder', cmap, 'NextPlot', 'ReplaceChildren');
       colormap(cmap); 
     
      
      h1 = subplot(1,3,1); 
      for nl = 1:nlines
      plot(-125:125, mean_not_sig_su(250*(nl-1)+1:250*nl+1), 'LineWidth',2,'color',cmap(nl,:));
      hold on
      end
      set(h1,'position',get(h1,'position').*[1 1 1.15 1])
     ylim(normylims(nc,:))
      % ylim([-0.01 1])
     xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('non sig') 
      
      h2 = subplot(1,3,2); 
      for nl = 1:nlines
      plot(-125:125, mean_origin(250*(nl-1)+1:250*nl+1), 'LineWidth',2,'color',cmap(nl,:));
      hold on
      end
      set(h2,'position',get(h2,'position').*[1 1 1.15 1])
      ylim(normylims(nc,:))
      %ylim([-0.01 1])
      xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('Overall Mean')
      
      h3 = subplot(1,3,3);
      for nl = 1:nlines 
      plot(-125:125, mean_sig_su(250*(nl-1)+1:250*nl+1),  'LineWidth',2, 'Color',cmap(nl,:) )
      hold on
      end
      set(h3,'position',get(h3,'position').*[1 1 1.15 1])
     ylim(normylims(nc,:))
  %   ylim([-0.01 1])
      xlim([-125 125])
      title('Adapting units mean')
      set(gca,'box','off')
      legend('peak 1', 'peak 2', 'peak 3', 'peak 4')
    
      
       % K norm ylim([-0.02 1.1])
      
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
     % if pn >1 
      %ax1 = gca;                   
      %ax1.YAxis.Visible = 'off';   
     % end

    currfig = gcf;
    %currfig.Children(end),
    sgtitle({sprintf('%s mean responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('signonsigadapt_stacked_commonnorm_nanmean_bscorr_origin_data_aligned_%s',cellclass(nc),'cells_cbrewer2_scaled')));
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
end 
%% Plot three subplots: nonsig| sig increase | sig adapt

pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];

%maps
nlines = 4;
cmaps = struct();
cmaps(1).map =cbrewer2('Reds', nlines);
cmaps(2).map =cbrewer2('Oranges', nlines);
cmaps(3).map =cbrewer2('Blues', nlines);

% M ylim([0.05 .75])
% P ylim([.1 .8])
% K ylim([-0.01 1])
normylims = [[0.05 .75];[.1 .8];[-0.01 1]];
 % P ylim([-20 85])
 % M ylim([-30 135])
 % K ylim([-20 200])
ylims = [[-30 135];[-20 85];[-20 200]];

for nc = 1:3
 clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
%%store all trials of a given peak, in a matrix, across all units
clear aligned_trials
aligned_trials = nan(1204, length(layer_idx));

clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:4
          % aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
           aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-1-124:max_low_dist(layer_idx(i))+125,:,pn),2);
       end
       %normalizing with max and min of each unit
       %norm_aligned_trials(:,i) = (aligned_trials(:,i) - min(aligned_trials(:,i)))/(max(aligned_trials(:,i))-min(aligned_trials(:,i)));
     
    end
end
  %normalizing with max and min across units
       maximum = max(aligned_trials,[],'all');
       minimum = min(aligned_trials,[],'all');
       norm_aligned_trials = (aligned_trials - minimum)/(maximum-minimum);

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))

 clear sig_adapsu sig_incsu non_sig_su mean_sig_su
 cntsigadapt = 0;
 cntsiginc = 0;
 cntnotsig = 0;
 all_mean_data = nan(4, length(layer_idx));
 sig_adapsu = nan(length(norm_aligned_trials(:,1)),length(layer_idx));
 sig_incsu = nan(length(norm_aligned_trials(:,1)),length(layer_idx));
 non_sig_su =nan(length(norm_aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
 mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cntsigadapt= cntsigadapt+1;
      sig_adapsu(:,cntsigadapt) = norm_aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
   if all_mean_data(4,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cntsiginc= cntsiginc+1;
      sig_incsu(:,cntsiginc) = norm_aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
  if  pvalues(layer_idx(nunit),4) > .05
      
      cntnotsig= cntnotsig+1;
      non_sig_su(:,cntnotsig) = norm_aligned_trials(:,nunit);
      
  end
      end
  end
  mean_not_sig_su = nanmean(non_sig_su,2);
  mean_sig_incsu = nanmean(sig_incsu,2);
  mean_sig_adapsu = nanmean(sig_adapsu,2);
 
 
  mean_origin =nanmean(norm_aligned_trials,2);

     
     figure();
       
       cmap = flip(cmaps(nc).map) ;
      % axes('ColorOrder', cmap, 'NextPlot', 'ReplaceChildren');
       colormap(cmap); 
     
      
      h1 = subplot(1,3,1); 
      for nl = 1:nlines
      plot(-125:125, mean_not_sig_su(250*(nl-1)+1:250*nl+1), 'LineWidth',2,'color',cmap(nl,:));
      hold on
      end
      set(h1,'position',get(h1,'position').*[1 1 1.15 1])
      ylim(normylims(nc,:))
     %  ylim([-0.01 1])
     xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('non sig') 
      
      h2 = subplot(1,3,2); 
      for nl = 1:nlines
      plot(-125:125, mean_sig_incsu(250*(nl-1)+1:250*nl+1), 'LineWidth',2,'color',cmap(nl,:));
      hold on
      end
      set(h2,'position',get(h2,'position').*[1 1 1.15 1])
      ylim(normylims(nc,:))
      %ylim([-0.01 1])
      xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('sig increase')
      
      h3 = subplot(1,3,3);
      for nl = 1:nlines 
      plot(-125:125, mean_sig_adapsu(250*(nl-1)+1:250*nl+1),  'LineWidth',2, 'Color',cmap(nl,:) )
      hold on
      end
      set(h3,'position',get(h3,'position').*[1 1 1.15 1])
     ylim(normylims(nc,:))
    % ylim([-0.01 1])
      xlim([-125 125])
      title('Adapting units mean')
      set(gca,'box','off')
      legend('peak 1', 'peak 2', 'peak 3', 'peak 4')
    
      
       % K norm ylim([-0.02 1.1])
      
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
     % if pn >1 
      %ax1 = gca;                   
      %ax1.YAxis.Visible = 'off';   
     % end

    currfig = gcf;
    %currfig.Children(end),
    sgtitle({sprintf('%s mean responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('signonsigadaptinc_stacked_commonnorm_nanmean_bscorr_origin_data_aligned_%s',cellclass(nc),'cells_cbrewer2_scaled')));
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
end

%% Plot two subplots: sig increase | sig adapt

pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];

%maps
nlines = 4;
cmaps = struct();
cmaps(1).map =cbrewer2('Reds', nlines);
cmaps(2).map =cbrewer2('Oranges', nlines);
cmaps(3).map =cbrewer2('Blues', nlines);

% M ylim([0.05 .75])
% P ylim([.1 .8])
% K ylim([-0.01 1])
normylims = [[0.05 .75];[.1 .8];[-0.01 1]];
 % P ylim([-20 85])
 % M ylim([-30 135])
 % K ylim([-20 200])
ylims = [[-30 135];[-20 85];[-20 200]];

for nc = 1:3
 clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
%%store all trials of a given peak, in a matrix, across all units
clear aligned_trials
aligned_trials = nan(1204, length(layer_idx));

clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:4
          % aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
           aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-1-124:max_low_dist(layer_idx(i))+125,:,pn),2);
       end
       %normalizing with max and min of each unit
       %norm_aligned_trials(:,i) = (aligned_trials(:,i) - min(aligned_trials(:,i)))/(max(aligned_trials(:,i))-min(aligned_trials(:,i)));
     
    end
end
  %normalizing with max and min across units
       maximum = max(aligned_trials,[],'all');
       minimum = min(aligned_trials,[],'all');
       norm_aligned_trials = (aligned_trials - minimum)/(maximum-minimum);

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))

 clear sig_adapsu sig_incsu non_sig_su mean_sig_su
 cntsigadapt = 0;
 cntsiginc = 0;
 cntnotsig = 0;
 all_mean_data = nan(4, length(layer_idx));
 sig_adapsu = nan(length(norm_aligned_trials(:,1)),length(layer_idx));
 sig_incsu = nan(length(norm_aligned_trials(:,1)),length(layer_idx));

  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
 mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cntsigadapt= cntsigadapt+1;
      sig_adapsu(:,cntsigadapt) = norm_aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
   if all_mean_data(4,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cntsiginc= cntsiginc+1;
      sig_incsu(:,cntsiginc) = norm_aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
 
      end
  end
 
  mean_sig_incsu = nanmean(sig_incsu,2);
  mean_sig_adapsu = nanmean(sig_adapsu,2);
 
 
  mean_origin =nanmean(norm_aligned_trials,2);

     
     figure();
       
       cmap = flip(cmaps(nc).map) ;
      % axes('ColorOrder', cmap, 'NextPlot', 'ReplaceChildren');
       colormap(cmap); 
     
      
      h2 = subplot(1,2,1); 
      for nl = 1:nlines
      plot(-125:125, mean_sig_incsu(250*(nl-1)+1:250*nl+1), 'LineWidth',2,'color',cmap(nl,:));
      hold on
      end
      set(h2,'position',get(h2,'position').*[1 1 1.15 1])
      ylim(normylims(nc,:))
      %ylim([-0.01 1])
      xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('sig increase')
      
      h3 = subplot(1,2,2);
      for nl = 1:nlines 
      plot(-125:125, mean_sig_adapsu(250*(nl-1)+1:250*nl+1),  'LineWidth',2, 'Color',cmap(nl,:) )
      hold on
      end
      set(h3,'position',get(h3,'position').*[1 1 1.15 1])
     ylim(normylims(nc,:))
    % ylim([-0.01 1])
      xlim([-125 125])
      title('Adapting units mean')
      set(gca,'box','off')
      legend('peak 1', 'peak 2', 'peak 3', 'peak 4')
    
      
       % K norm ylim([-0.02 1.1])
      
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
     % if pn >1 
      %ax1 = gca;                   
      %ax1.YAxis.Visible = 'off';   
     % end

    currfig = gcf;
    %currfig.Children(end),
    sgtitle({sprintf('%s mean responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('sigadaptinc_stacked_commonnorm_nanmean_bscorr_origin_data_aligned_%s',cellclass(nc),'cells_cbrewer2_scaled')));
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
end

   %{
   
      %% plot with significant addapting units normalizing adapting units and non adapting units
   
pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

layer_idx = find(strcmp(layer, 'M'));
%%store all trials of a given peak, in a matrix, across all units
%aligned_trials = nan(1204, length(layer_idx));
aligned_trials = nan(1204, length(layer_idx));

%norm_aligned_trials = nan(1204, length(layer_idx));
norm_aligned_trials = nan(1204, length(layer_idx));

clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:4
          % aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
           aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-1-124:max_low_dist(layer_idx(i))+125,:,pn),2);
       end
            
    end
end

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))

 clear sig_su mean_sig_su
  cnt = 0;
  cnt2 =0;
 all_mean_data = nan(4, length(layer_idx));
 sig_su = nan(length(aligned_trials(:,1)),length(layer_idx));
 non_sig_su = nan(length(aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
             mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
             all_mean_data(:,nunit) = mean_data;
                   if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
                           cnt= cnt+1;
                           sig_su(:,cnt) = aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
                   else
                       cnt2 = cnt2+1;
                       non_sig_su(:,cnt2) = aligned_trials(:,nunit);
                   end
      end
  end
    %normalizing with max and min across units of each group
       maxsig = max(sig_su,[],'all');
       minsig = min(sig_su,[],'all');
      % norm_sig_su = (sig_su - minsig)/(maxsig-minsig);
       
       maxnonsig = max(non_sig_su,[],'all');
       minnonsig = min(non_sig_su,[],'all');
       %norm_non_sig_su = (non_sig_su - minnonsig)/(maxnonsig-minnonsig);
       
       maximum = max(aligned_trials,[],'all');
       minimum = min(aligned_trials,[],'all');
       norm_sig_su = (sig_su - minimum)/(maximum-minimum);
       norm_non_sig_su = (non_sig_su - minimum)/(maximum-minimum);
  
      mean_sig_su = nanmean(norm_sig_su,2);
     
      mean_nonsig_su = nanmean(norm_non_sig_su,2);


      figure();
     
      ci_low = mean_nonsig_su - 1.96*std(norm_non_sig_su,0,2, 'omitnan')./sqrt(length(norm_non_sig_su(1,:)));
      ci_high = mean_nonsig_su + 1.96*std(norm_non_sig_su,0,2, 'omitnan')./sqrt(length(norm_non_sig_su(1,:)));
      
      for pn = 1:4
      h =subplot(1,4,pn);
      plot(-125:125, mean_nonsig_su(250*(pn-1)+1:250*pn+1),'LineWidth',1, 'Color', [215/255 25/255 28/255] )
      hold on
      plot(-125:125, ci_low(250*(pn-1)+1:250*pn+1),':', 'LineWidth',.7,'Color', [215/255 25/255 28/255] )
      hold on
      plot(-125:125, ci_high(250*(pn-1)+1:250*pn+1),':', 'LineWidth',.7,'Color',[215/255 25/255 28/255]  )
      hold on
      plot(-125:125, mean_sig_su(250*(pn-1)+1:250*pn+1),  'LineWidth',1, 'Color',[141/255 140/255 140/255] )
     set(h,'position',get(h,'position').*[1 1 1.15 1])
      % P ylim([-20 85])
      % M ylim([-35 135])
      % K ylim([-20 200])
      
      ylim([0.07 .75])
      
      %P ylim([0.1 .8])
      
      
     % K ylim([-0.02 1.07])
      %M ylim([0.07 .75])
       % K norm ylim([-0.02 1.1])
      xlim([-125 125])
      set(gca,'box','off')
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
      if pn >1 
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off';   
      end
      
      end
 
  
 %red [215/255 25/255 28/255] 
 %orange [253/255 174/255 97/255]
 %blue [44/255 123/255 182/255]

 %brown [165/255 42/255 42/255]
    currfig = gcf;
    
    title(currfig.Children(end),{'M mean responses'}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    ylh = ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
  
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020\plots\',strcat('sigadapt_subplots_sig_nonsig_allnorm_95ci_nanmean_bscorr_origin_data_aligned_Mcells_red'));
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
   %}




 x = 0:0.01:2;
       sc = [0.5; 1; 2];

       t0 = [0; 0.2; 0.4];

       t = bsxfun(@rdivide, bsxfun(@plus, x, t0), sc);

       y = sin(t * 2 * pi);

       cmap = cbrewer2('Reds', numel(sc));

       axes('ColorOrder', cmap, 'NextPlot', 'ReplaceChildren');

       plot(x, y);