%this script gathers most plots used for the adaptation manuscript
%Written and edited by Loic Daumail 26/6/2020

%%Figure 1: plot an example of a response (use i=1)
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
    origin_dSUA = data_file.clean_origin_data(i).unit(401:1900,:)%-mean(data_file.clean_origin_data(i).unit(401:600,:),1);
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

%% Figure 1: Example Unit plot

%%% plot aligned units to peak of interest
layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
for i = 1:length(suas_trials)
     if ~isnan(max_low_dist(i))
        figure();
      mean_unit = squeeze(nanmean(suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,:),2));
      stdev = squeeze(std(suas_trials(i).aligned(max_low_dist(i)-1-124:max_low_dist(i)-1+125,:,:),[],2));
      for pn =1:4
      h = subplot(1,4,pn);
      plot(-125:124, mean_unit(:,pn));
      hold on
      h1= ciplot( mean_unit(:,pn)+ 1.96*stdev(:,pn)/sqrt(14), mean_unit(:,pn)-1.96*stdev(:,pn)/sqrt(14),[-125:124],[40/255 40/255 40/255],0.1);
      set(h1, 'edgecolor','none') 
     
      set(h,'position',get(h,'position').*[1 1 1.15 1])
      ylim([0 190])
      xlim([-125 125])
      
      set(gca,'box','off')
      set(gca, 'linewidth',2)
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      
      if pn > 1
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off'; 
      end
      end
    

    sgtitle({'Baseline corrected M cell trials'}, 'Interpreter', 'none')
    xlabel('Resolution (ms)')
   set(gcf,'Units','inches') 
   %set(gcf,'position',[1 1 8.5 11])
        filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\aligned_trials_unit1');
        saveas(gcf, strcat(filename, '.png'));
        saveas(gcf, strcat(filename, '.svg')); 
     end

end

%% Figure 1: plot number of cells in a bar plot

y = [18 15 3];
x = categorical({'M' 'P' 'K'});
x = reordercats(x,{'M' 'P' 'K'});
%{
col(1,:) =[102/255 194/255 165/255] ; %
col(2,:) = [252/255 141/255 98/255]; % -- 
col(3,:) = [141/255 160/255 203/255]; % -- 

col(1,:) =[146/255 197/255 222/255] ; %--blue 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [166/255 219/255 160/255]; % -- green

col(1,:) =[194/255 165/255 207/255] ; %--purple
col(2,:) = [253/255 174/255 97/255]; % -- orange
col(3,:) = [166/255 219/255 160/255]; % -- green
%}
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue

figure();
b=bar(x,y,'FaceColor','flat', 'BarWidth', 0.5);
%b.FaceColor = 'flat';
b.CData = col;
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
set(gca,'box','off')
set(gca, 'linewidth',2)

 filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\q_cel1_class_sample_sizes_stimcol');
        saveas(gcf, strcat(filename, '.png'));
        saveas(gcf, strcat(filename, '.svg')); 

%% Figure 1: Plot the spike waveforms

gendatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [gendatadir 'refined_dataset']; 
gen_data_file = load(channelfilename);

%plot spike waveform of the retained channels for every retained signle unit
%use the peakvals to see which unit was retained in the data cleaning
%process
channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);


layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];
spikes = nan(81,25,3);

  
for nc = 1:length(cellclass)
 clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));

    for un = 1:length(layer_idx)
         if ~isempty(peakvals.peak_vals(layer_idx(un)).peak)
         spike_dat = gen_data_file.new_data(layer_idx(un)).channel_data.wf.waveForms;
%mean of just one unit

        mean_data = squeeze(mean(spike_dat,2));
        chan_idx = str2double(gen_data_file.new_data(layer_idx(un)).channel_data.chan(3:4)); 
        spikes(:,un,nc) = mean_data(chan_idx,:);
         end
    end 
    
end

%spikes = reshape(spikes,81,75);

%index = reshape(1:75, 3, 25).';


 figure();
   for nc = 1:3
        %for un =1:length(1:75)
     % h =subplot(length(1:25),3,index(un));
     h =subplot(length(1),3,nc);
    plot(-40:40, spikes(:,:,nc)-mean(spikes(:,:,nc),1),  'LineWidth',2)
   % hold on
   % plot(xlim, [0 0],'k')

      
       set(h,'position',get(h,'position').*[1 1 1.15 1])
       
      set(gca,'box','off')
      
       
      ax1 = gca; 
      if nc > 1
      ax1.YAxis.Visible = 'off';   
      ax1.XAxis.Visible = 'off';
      end
    end
    %end
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\spike_waveform\all_units_spike_waveforms_stacked');
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
   
   %% Figure 2: cell class response plots
   
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
    origin_dSUA = data_file.clean_origin_data(i).unit(401:1900,:); % -mean(data_file.clean_origin_data(i).unit(401:600,:),1);
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

 % plot overall mean
   
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
red = [215/255 25/255 28/255];
 orange = [253/255 174/255 97/255];
 blue =[44/255 123/255 182/255];
colors = [red; orange; blue];
%}
%{
col(1,:) =[102/255 194/255 165/255] ; %-- green
col(2,:) = [252/255 141/255 98/255]; % -- orange
col(3,:) = [141/255 160/255 203/255]; % -- blue

col(1,:) =[194/255 165/255 207/255] ; %--purple 
col(2,:) = [166/255 219/255 160/255]; % -- green
col(3,:) = [253/255 174/255 97/255]; % -- orange

col(1,:) =[146/255 197/255 222/255] ; %--blue 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [166/255 219/255 160/255]; % -- green


col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [238/255 58/255 104/255]; % -- pink
col(3,:) = [156/255 203/255 59/255]; % -- green


col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue


col(1,:) =[194/255 165/255 207/255] ; %--purple
col(2,:) = [253/255 174/255 97/255]; % -- orange
col(3,:) = [166/255 219/255 160/255]; % -- green

col(1,:) =[146/255 197/255 222/255] ; %--blue 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [166/255 219/255 160/255]; % -- green
%}
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue

% M ylim([0.05 .75])
% P ylim([.1 .8])
% K ylim([-0.01 1])
%normylims = [[0.05 .75];[.1 .8];[-0.01 1]];
 % P ylim([-20 85])
 % M ylim([-35 135])
 % K ylim([-20 200])
ylims = [[0 160];[0 105];[0 210]];
%ylims = [[-35 135];[-20 85];[-20 200]];

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
    %   maximum = max(aligned_trials,[],'all');
     %  minimum = min(aligned_trials,[],'all');
      % norm_aligned_trials = (aligned_trials - minimum)/(maximum-minimum);

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))
%
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
      sig_ci_low = mean_sig_su - 1.96*std(sig_adapsu,0,2, 'omitnan')./sqrt(length(find(~isnan(sig_adapsu(1,:)))));
      sig_ci_high = mean_sig_su + 1.96*std(sig_adapsu,0,2, 'omitnan')./sqrt(length(find(~isnan(sig_adapsu(1,:)))));
      
  
 
      figure();
      mean_origin =nanmean(aligned_trials,2);
     % ci_low = mean_origin - 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(find(~isnan(aligned_trials(1,:)))));
     % ci_high = mean_origin + 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(find(~isnan(aligned_trials(1,:)))));
     ci_low = mean_origin - 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(find(~isnan(aligned_trials(1,:)))));
     ci_high = mean_origin + 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(find(~isnan(aligned_trials(1,:)))));
      
      for pn = 1:4
      h =subplot(1,4,pn);
      %{
      h2= ciplot(sig_ci_low(250*(pn-1)+1:250*pn+1), sig_ci_high(250*(pn-1)+1:250*pn+1),[-125:125],[40/255 40/255 40/255],0.2);
      set(h2, 'edgecolor','none')  
      hold on
      %}
      
      plot(-125:125, mean_origin(250*(pn-1)+1:250*pn+1),'LineWidth',2, 'Color',[40/255 40/255 40/255] )
      hold on
      h1= ciplot(ci_low(250*(pn-1)+1:250*pn+1), ci_high(250*(pn-1)+1:250*pn+1),[-125:125],col(nc,:),0.5);
      set(h1, 'edgecolor','none') 
      hold on 
     
      plot(-125:125, mean_sig_su(250*(pn-1)+1:250*pn+1),'LineWidth',2, 'Color',[140/255 140/255 140/255] ) 
      set(h,'position',get(h,'position').*[1 1 1.15 1])
     
        ylim(ylims(nc,:))
      
      
       % K norm ylim([-0.02 1.1])
      xlim([-125 125])
      set(gca, 'linewidth',2)
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

   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('q_origin_data_aligned_%s_cells_95ci_sigadap', cellclass(nc))));
  % saveas(gcf, strcat(filename, '.svg')); 
   %saveas(gcf, strcat(filename, '.png')); 
end
   

%% (part of Figure 2) Plot overall across all cell classes and mean response of one cell for each cell class

%% find spike waveforme that correspond to single unit example plot (Figure 2)

    
   %% example single unit

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
    origin_dSUA = data_file.clean_origin_data(i).unit(401:1900,:); % -mean(data_file.clean_origin_data(i).unit(401:600,:),1);
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

 % plot example
   
pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

%waveform data
gendatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [gendatadir 'refined_dataset']; 
gen_data_file = load(channelfilename);


layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
cellclass = [ 'M', 'P', 'K'];

col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue

% M ylim([0.05 .75])
% P ylim([.1 .8])
% K ylim([-0.01 1])
%normylims = [[0.05 .75];[.1 .8];[-0.01 1]];
 % P ylim([-20 85])
 % M ylim([-35 135])
 % K ylim([-20 200])
ylims = [[0 240];[0 120];[0 230]];
%ylims = [[-35 135];[-20 85];[-20 200]];
unitNum = [5,4,1];
for nc = 1:3
    clear layer_idx
    layer_idx = find(strcmp(layer, cellclass(nc)));
    %%store all trials of a given peak, in a matrix, across all units
    %aligned_trials = nan(1204, length(layer_idx));
    clear numTrials
    numTrials = nan(length(layer_idx),1);

    clear i 
    for i = 1:length(layer_idx)
        if ~isnan(max_low_dist(layer_idx(i)))
             %store number of trials for each unit to determine the max length
             numTrials(i) = length(suas_trials(layer_idx(i)).aligned(1,:,1));  
        else
            numTrials(i) = 0;

        end
         
    end

    clear aligned_trials
    aligned_trials = nan(1204, max(numTrials), length(layer_idx));
    clear i
    for i = 1:length(layer_idx)
        if ~isnan(max_low_dist(layer_idx(i)))
           for pn = 1:4
               aligned_trials(250*(pn-1)+1:250*pn+1,1:numTrials(i),i)= suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-1-124:max_low_dist(layer_idx(i))+125,:,pn);

           end
         end
    end

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))
%
     clear sig_adapsu spikes spike_dat
      cnt = 0;
     all_mean_data = nan(4, length(layer_idx));
     sig_adapsu = nan(length(aligned_trials(:,1)), max(numTrials),length(layer_idx));
     spikes = nan(3000,81,length(layer_idx));
      for nunit = 1:length(layer_idx)
          if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
              mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
              all_mean_data(:,nunit) = mean_data;
              if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
                  cnt= cnt+1;
                  sig_adapsu(:,find(~all(isnan(aligned_trials(:,:,nunit)))),cnt) = aligned_trials(:, ~all(isnan(aligned_trials(:,:,nunit))),nunit); 
                  
                  %isolate corresponding spike waveform
                    spike_dat = gen_data_file.new_data(layer_idx(nunit)).channel_data.wf.waveForms;
                    %mean_data = squeeze(mean(spike_dat,2));%mean of just one unit
                    chan_idx = str2double(gen_data_file.new_data(layer_idx(nunit)).channel_data.chan(3:4)); 
                    spikes(:,:,cnt) = squeeze(spike_dat(1,:,chan_idx,:));
         % plot(x_stim,norm_chan(:, nunit)')
         %hold on
              end
          end
      end

      
      %plot example single unit
      idxs = find(~isnan(sig_adapsu(1,1,:)));
      clear  meanSigSu unit_trials
     % for n =1:length(find(~isnan(sig_adapsu(1,1,:))))
          %take trials data of only one significant unit to plot it
          unit_trials = squeeze(sig_adapsu(:,:,idxs(unitNum(nc))));
          meanSigSu = nanmean(unit_trials,2);

          %mean_sig_su = nanmean(sig_adapsu,2);
          sig_ci_low = meanSigSu - 1.96*std(unit_trials(:,~all(isnan(unit_trials))),0,2, 'omitnan')./sqrt(length(find(~isnan(unit_trials(1,:)))));
          sig_ci_high = meanSigSu + 1.96*std(unit_trials(:,~all(isnan(unit_trials))),0,2, 'omitnan')./sqrt(length(find(~isnan(unit_trials(1,:)))));



          figure();

          for pn = 1:4
              h =subplot(1,4,pn);

              plot(-125:125, meanSigSu(250*(pn-1)+1:250*pn+1),'LineWidth',2, 'Color',[40/255 40/255 40/255] )
              hold on
              h1= ciplot(sig_ci_low(250*(pn-1)+1:250*pn+1), sig_ci_high(250*(pn-1)+1:250*pn+1),[-125:125],col(nc,:),0.5);
              set(h1, 'edgecolor','none') 
              set(h,'position',get(h,'position').*[1 1 1.15 1])

                ylim(ylims(nc,:))


               % K norm ylim([-0.02 1.1])
              xlim([-125 125])
              set(gca, 'linewidth',2)
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

        title(currfig.Children(end),{sprintf('%s single unit mean responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)

       % xlabel('\fontsize{14}Resolution (ms)')
        ylh = ylabel({'\fontsize{14}Spike Rate (spikes/s)'});

       set(gcf,'Units','inches') 
       set(gcf,'position',[1 1 15 11])


       filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('example_%s_unit_95ci_sigadap', cellclass(nc))));
       %saveas(gcf, strcat(filename, '.svg')); 
       %saveas(gcf, strcat(filename, '.png')); 
       
       
       %plot corresponding spike waveform
       meanWf = mean(spikes(:,:,unitNum(nc)),1); %-mean(mean(spikes(:,:,unitNum(nc)),1));
       figure();
       plot(-40:40,meanWf,  'LineWidth',2,'Color',[40/255 40/255 40/255] )
       hold on 
       wf_ci_low = meanWf - 1.96*std(spikes(:,:,unitNum(nc)),0,1, 'omitnan')./sqrt(length(find(~isnan(spikes(:,1,unitNum(nc))))));
       wf_ci_high =  meanWf + 1.96*std(spikes(:,:,unitNum(nc)),0,1, 'omitnan')./sqrt(length(find(~isnan(spikes(:,1,unitNum(nc))))));
       h2= ciplot(wf_ci_low, wf_ci_high,[-40:40],col(nc,:),0.5);
       set(h2, 'edgecolor','none') 
            
       set(gcf,'position',get(gcf,'position').*[1 1 1.15 1])
       set(gca,'box','off')
      
       
      ax1 = gca; 
      
      ax1.YAxis.Visible = 'off';   
      ax1.XAxis.Visible = 'off';
       
      
      filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('example_%s_unit_spike_wf_95CI', cellclass(nc))));
       saveas(gcf, strcat(filename, '.svg')); 
      
end

 %% Figure 4: Power plots
 
%% Plots for the method process explaination
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_origin_sup_50']; 
data_file = load(channelfilename);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

%layer_idx = find(strcmp(layer, 'M'));
%log_p_layer = zeros(length(layer),1);
%log_p_layer(layer_idx) = logical(layer_idx);


%% compute the power spectrum (spectrogram) using mtspecgramc function

Ses = struct();
bs_data = struct();
channum = 1: length(data_file.clean_origin_data);
mean_S = nan(1174,38, length(channum));

xabs = -100:1301;

%filtered_dMUA = nan(length(xabs), length(channum));
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];

clear i ;
i=19;
 %for i = 1:length(channum)
     if ~isempty(data_file.clean_origin_data(i).unit)
data = squeeze(data_file.clean_origin_data(i).unit(401:1900,:,:));
   bsl = mean(data(1:200,:));
   
   norm_mean_bs = data(72:end, :) - bsl;
   namelist1(1,1:length(sprintf('chan_%d',i))) = sprintf('chan_%d',i);
   bs_data(i).namelist1 = norm_mean_bs;
clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,:) ,movingwin, params); 
 
namelist2(1,1:length(sprintf('S_%d',i))) = sprintf('S_%d',i);
Ses(i).namelist2 = S;
mean_S(:,:,i) = nanmean(S,3);

tvec = t*1000 -129;

%time_adj = 1:128;
%tvec = cat(2, time_adj , t*1000) ;
%we can also store tvec and f in a struct, but they are all identical
     end
 %end



%% Power spectrum plot
figure, 

imagesc(tvec,sort(f),mean_S(:,:,i)')
ylim([2 20]); 
set(gca,'ydir','normal')
title({'Mean spectrogram', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus onset(ms)')
    ylabel('Frequency band (Hz)')
filename = strcat('C:\Users\daumail\Documents\first_year_committee_meeting\mean_powersepc_imagesc');
saveas(gcf, strcat(filename, '.png')); 
saveas(gcf, strcat(filename, '.svg')); 
%% plot the mean power only in the 5Hz range for 1 unit across trials 

figure, 
%normspec = (mean_S(:,:,i) - min(mean_S(:,:,i)))./(max(mean_S(:,:,i)) - min(mean_S(:,:,i)));
%normS = (S(:,1,:)-min(S(:,1,:)))./(max(S(:,1,:)) - min(S(:,1,:)));
%cihigh = normspec(:,1) + 1.96*std(normS(:,1,:),[],3)/sqrt(length(S(1,1,:)));
%cilow = normspec(:,1) - 1.96*std(normS(:,1,:),[],3)/sqrt(length(S(1,1,:)));
cihigh = mean_S(:,1,i) + 1.96*std(S(:,1,:),[],3)/sqrt(length(S(1,1,:)));
cilow = mean_S(:,1,i) - 1.96*std(S(:,1,:),[],3)/sqrt(length(S(1,1,:)));
plot(tvec,mean_S(:,1,i))
hold on
h1= ciplot( cihigh,cilow,tvec,[40/255 40/255 40/255],0.1);
set(h1, 'edgecolor','none') 
set(gca, 'box', 'off')
      
title({'Mean power at 4 Hz', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus onset(ms)')
    ylabel('Power (spikes^2/sec^2)')
filename = strcat('C:\Users\daumail\Documents\first_year_committee_meeting\mean_power_4hz_1unit_origin');
saveas(gcf, strcat(filename, '.png')); 
saveas(gcf, strcat(filename, '.svg')); 


reps   = 10000;
all_sigs95 = nan(length(Ses),1);
all_sigs90 = nan(length(Ses),1);

part1 = nanmean(squeeze(Ses(i).namelist2(1:575,1,:)), 1);
part2 = nanmean(squeeze(Ses(i).namelist2(576:1150,1,:)), 1);

    if nanmean(part1) > nanmean(part2)
    cond1               = part1;
    cond2               = part2;
    else
    cond1               = part2;
    cond2               = part1;
    end
     %% Plot T1 and T2 power distributions of one unit across trials 
 figure(); boxplot([cond1' cond2'],'notch','off','labels',{'cond1','cond2'}); hold on    
x=ones(length(cond1),1).*(1+(rand(length(cond1),1)-0.5)/5);
x1=ones(length(cond1),1).*(1+(rand(length(cond1),1)-0.5)/10);
f1=scatter(x,cond1,'k','filled');f1.MarkerFaceAlpha = 0.4;hold on
f2=scatter(x1*2,cond2,'k','filled');f1.MarkerFaceAlpha = 0.4;
ylim([0 400])
ylabel('Power (spikes^2/s^2)')
filename = strcat('C:\Users\daumail\Documents\first_year_committee_meeting\distributions_cond1_cond2_notchoff');
saveas(gcf, strcat(filename, '.png')); 
saveas(gcf, strcat(filename, '.svg')); 



  
%% Plot AUC of our two samples

    [X,Y,T,AUC]           = perfcurve([ones(length(cond1),1); repmat(2,length(cond2),1)],[cond1 cond2],1);
    NP                    = length(cond2);
    PR                    = length(cond1);
    catdat                = [cond1 cond2];
figure();
plot(X,Y)
set(gca, 'box','off')   

%% Plot  Receiver Operating Characteristics curve examples

 shufPR         = catdat(randperm(length(catdat),PR));
 shufNP         = catdat(randperm(length(catdat),NP));
[Xshuf,Yshuf,~, shufAUC]    = perfcurve([ones(PR,1); repmat(2,NP,1)],[shufPR shufNP],1);

figure();
plot(Xshuf,Yshuf)
set(gca, 'box','off')

 %% Plot randomized distribution
 

    for r       = 1:reps
       clear shufNP shufPR
       shufPR         = catdat(randperm(length(catdat),PR));
       shufNP         = catdat(randperm(length(catdat),NP));
       [~,~,~,...
       shufAUC(r)]    = perfcurve([ones(PR,1); repmat(2,NP,1)],[shufPR shufNP],1);
    end
 critT95         = quantile(shufAUC,.95);
 
    figure();
    histogram(shufAUC)
    hold on
    plot([critT95, critT95], [0, 500])
    hold on
    plot([AUC, AUC], [0, 500], 'Color', 'k')
    set(gca, 'box', 'off')
filename = strcat('C:\Users\daumail\Documents\first_year_committee_meeting\shufAUC_distribution');
saveas(gcf, strcat(filename, '.png')); 
saveas(gcf, strcat(filename, '.svg'));  %hold on
  %  histogram(shufPR)
   %  hold on
   % histogram(shufNP)
   
    critT90         = quantile(shufAUC,.90);



%% Figure 4 Plot mean power  with error bars before stimulation and during stimulation in the same analysis
%% normalize mean SUA before computing the grand cell class mean
%% plotting with power significant changes
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_origin_sup_50']; 
data_file = load(channelfilename);


channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
 peakvals = load([channeldir 'all_data_peaks']);
 sig95_idx = load( strcat(newdatadir,'roc_results95_stimonset_to1150ms.mat'));
  contrast = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

Ses = struct();
bs_data = struct();
channum = 1: length(layer);
mean_S_stim = nan(1646+128,38, length(channum));
%compute the power spectrum
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
 

clear i ;
 for i = 1:length(layer)
      if ~isempty(data_file.clean_origin_data(i).unit)
data = squeeze(data_file.clean_origin_data(i).unit(1:1901,:));
   bsl = mean(data(400:599,:));
   %stim and bl data
   norm_mean_bs = nan(length(data(:,1)),1,length(data(1,:)));
   norm_mean_bs = data(1:end,:) - bsl;
  

clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,:) ,movingwin, params); 

mean_S_stim(129:end,:,i) = nanmean(S,3);
      end

%we can also store tvec and f in a struct, but they are all identical
 end
 
time_adj = -599:-472;
x_stim = cat(2, time_adj , t*1000 -600) ;
cellclass = ['M', 'P', 'K'];

col(1,:) =[146/255 197/255 222/255] ; %--blue 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [166/255 219/255 160/255]; % -- green


for nc = 1:3
clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
%log_p_layer = zeros(length(layer),1);
%log_p_layer(layer_idx) = logical(layer_idx);
 

%here we compute the individual normalized units necessary for the variance
%and the data that we are about to plot
%for both the baseline data and the stimulus data
clear norm_chan normspec sig_su
norm_chan = nan(length(mean_S_stim(:,1,1)), length(layer_idx));
clear i;
for i = 1:length(layer_idx)
min_chan =min(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
max_chan = max(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
norm_chan(:,i) = (squeeze(mean_S_stim(:,1,layer_idx(i)))-min_chan)./(max_chan - min_chan);
end

normspec = nanmean(norm_chan,2);
ci_low = normspec(:,1) - 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
ci_high = normspec(:,1) + 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));

figure,
 h =subplot(1,1,1);
 plot(x_stim,normspec, 'LineWidth',1, 'Color',[40/255 40/255 40/255])

 hold on
 %when using ciplot, make sur there is no NAN in the vectors
 h1= ciplot(ci_low(130:end), ci_high(130:end),[-470:1174],col(nc,:),.5);
     set(h1, 'edgecolor','none') 
 hold on
 
 plot([0 0], ylim,'k')
 hold on
 plot([1150 1150], ylim,'k')
 hold on 
 
 
 cnt = 0;
  for nunit = 1:length(layer_idx)
 
   
part1 = nanmean(norm_chan(601:1175,nunit),1);
part2 = nanmean(norm_chan(1176:1750,nunit), 1);


if part1 > part2 && sig95_idx.all_sigs95(layer_idx(nunit)) ==1
    cnt = cnt +1;
     sig_su(:,cnt) = norm_chan(:,nunit); 
       % plot(x_stim,norm_chan(:, nunit)')
     %hold on
end
     
   
  end
  
  mean_sig_su = mean(sig_su,2);
  plot(x_stim, mean_sig_su,  'LineWidth',1, 'Color',[141/255 140/255 140/255] )
   xlim([-600 1250])
   ylim([-0.1 1.15])
   xticks([-500:500:1000])
 %ylim([-0.8 1.2])
     set(gca, 'linewidth',2)
      set(gca,'box','off') 
      xlabel('Time from stimulus onset(ms)')
   ylabel('Normalized Power at 4Hz(no units)')
title({sprintf('%s class cells mean power at 4Hz vs time normalized', cellclass(nc))}, 'Interpreter', 'none')

   % legend('Mean', 'Mean-1.96*sem', 'Mean+1.96*sem', 'Mean significant decrease', 'Location', 'bestoutside')
    
filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\power_plots\',contrast{2},sprintf('_%s_layer_4hz_gathered_sig_suamean_pow',cellclass(nc)));
saveas(gcf, strcat(filename, '.svg'));
saveas(gcf, strcat(filename, '.png'));
end

%% figure 6:

%% Making the data ready to plot peaks

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
    origin_dSUA = data_file.clean_origin_data(i).unit(401:1900,:); % -mean(data_file.clean_origin_data(i).unit(401:600,:),1);
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

   %% Plot three subplots: nonsig|overall mean | sig adapt

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
nlines = 5;
cmaps = struct();


cmaps(1).map =cbrewer2('Greys', nlines);
cmaps(2).map =cbrewer2('Reds', nlines);
cmaps(3).map =cbrewer2('Blues', nlines)

%{
cmaps(1).map =cbrewer2('Purples', nlines);
cmaps(2).map =cbrewer2('Oranges', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
%}
% M ylim([0.05 .75])
% P ylim([.1 .8])
% K ylim([-0.01 1])
%normylims = [[0.05 .75];[.1 .8];[-0.01 1]];
 % P ylim([-20 85])
 % M ylim([-30 135])
 % K ylim([-20 200])
ylims = [[0 160];[0 105];[0 210]];

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
 cntall = 0;
 cntnotsig = 0;
 all_mean_data = nan(4, length(layer_idx));
 sig_adapsu = nan(length(aligned_trials(:,1)),length(layer_idx));
 all_su = nan(length(aligned_trials(:,1)),length(layer_idx));
 non_sig_su =nan(length(aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
 mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cntsigadapt= cntsigadapt+1;
      sig_adapsu(:,cntsigadapt) = aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
      cntall= cntall+1;
      all_su(:,cntall) = aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
 
  if  pvalues(layer_idx(nunit),4) > .05
      
      cntnotsig= cntnotsig+1;
      non_sig_su(:,cntnotsig) = aligned_trials(:,nunit);
      
  end
      end
  end
  mean_not_sig_su = nanmean(non_sig_su,2);
  mean_all_su = nanmean(all_su,2);
  mean_sig_adapsu = nanmean(sig_adapsu,2);
 
  %counts
  %cnts = [cntnotsig; cntall; cntsigadapt];
 
  mean_origin =nanmean(aligned_trials,2);

     
     figure();
       
       cmap = flip(cmaps(nc).map) ;
         colormap(cmap); 
     
      
      h1 = subplot(1,3,1); 
      for nl = 1:4
      plot(-125:125, mean_not_sig_su(250*(nl-1)+1:250*nl+1), 'LineWidth',2,'color',cmap(nl,:));
      hold on
      end
      text(-50,50,sprintf('Number of cells: %d', cntnotsig))
      set(h1,'position',get(h1,'position').*[1 1 1.15 1])
      ylim(ylims(nc,:))
     %  ylim([-0.01 1])
     xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('non sig') 
      
      h2 = subplot(1,3,2); 
      for nl = 1:4
      plot(-125:125, mean_all_su(250*(nl-1)+1:250*nl+1), 'LineWidth',2,'color',cmap(nl,:));
      hold on
      end
      
      text(-50,50,sprintf('Number of cells: %d', cntall))
      set(h2,'position',get(h2,'position').*[1 1 1.15 1])
      ylim(ylims(nc,:))
      %ylim([-0.01 1])
      xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('overall mean')
      
      h3 = subplot(1,3,3);
      for nl = 1:4 
      plot(-125:125, mean_sig_adapsu(250*(nl-1)+1:250*nl+1),  'LineWidth',2, 'Color',cmap(nl,:) )
      hold on
      end
      text(-50,50,sprintf('Number of cells: %d', cntsigadapt))
      set(h3,'position',get(h3,'position').*[1 1 1.15 1])
     ylim(ylims(nc,:))
    % ylim([-0.01 1])
      xlim([-125 125])
      title('Adapting units mean')
      set(gca,'box','off')
      legend('peak 1', 'peak 2', 'peak 3', 'peak 4')
    
      
       % K norm ylim([-0.02 1.1])
      
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
      if pn >1 
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off';   
      end

    currfig = gcf;
    %currfig.Children(end),
    sgtitle({sprintf('%s mean responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('p_signonsigadapt_stacked_mean_origin_data_aligned_%s_cells_cbrewer2',cellclass(nc))));
   %saveas(gcf, strcat(filename, '.svg')); 
   %saveas(gcf, strcat(filename, '.png')); 
end

%%
   %% Figure 6 - Plot 1 subplot for each peak each time including nonsig - overall mean - sig adapt - sig facilitated

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
nlines = 5;
cmaps = struct();


cmaps(1).map = cbrewer2('Greys', nlines);
cmaps(2).map = cbrewer2('Reds', nlines);
cmaps(3).map = cbrewer2('Blues', nlines);


ylims = [[0 160];[0 105];[0 210]];

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
  

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))

 clear sig_adapsu sig_incsu non_sig_su mean_sig_su
 
 cntnotsig = 0;
 cntsigfac = 0;
 cntsigadapt = 0;
 cntall = 0;
 
 non_sig_su =nan(length(aligned_trials(:,1)),length(layer_idx));
 sig_adapsu = nan(length(aligned_trials(:,1)),length(layer_idx));
 sig_facsu = nan(length(aligned_trials(:,1)),length(layer_idx));
 all_su = nan(length(aligned_trials(:,1)),length(layer_idx));
 all_mean_data = nan(4, length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
          
          mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
          all_mean_data(:,nunit) = mean_data;
          
          cntall= cntall+1;
          all_su(:,cntall) = aligned_trials(:,nunit);
        
          if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
              cntsigadapt= cntsigadapt+1;
              sig_adapsu(:,cntsigadapt) = aligned_trials(:,nunit);
              % plot(x_stim,norm_chan(:, nunit)')
              %hold on
          end
           % plot(x_stim,norm_chan(:, nunit)')
          %hold on
          
          if all_mean_data(4,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
              cntsigfac = cntsigfac+1;
              sig_facsu(:,cntsigfac) = aligned_trials(:,nunit);
          end
          
          if  pvalues(layer_idx(nunit),4) > .05
              
              cntnotsig= cntnotsig+1;
              non_sig_su(:,cntnotsig) = aligned_trials(:,nunit);
              
          end
          
      end
  end
  mean_not_sigsu = nanmean(non_sig_su,2);
  mean_sig_adapsu = nanmean(sig_adapsu,2);
  mean_sig_facsu = nanmean(sig_facsu,2);
  mean_allsu = nanmean(all_su,2);
  
  gathered_means = [mean_not_sigsu, mean_sig_adapsu, mean_sig_facsu, mean_allsu];
    
  
  figure();  
         cmap = flip(cmaps(nc).map) ;
         colormap(cmap); 
    for pn = 1:4
        h1 = subplot(1,4,pn); 
      for nl = 1:length(gathered_means(1,:)) 
     
           plot(-125:125, gathered_means(250*(pn-1)+1:250*pn+1,nl), 'LineWidth',2,'color',cmap(nl,:));
           hold on
           set(gca,'box','off')
         
      end
       set(h1,'position',get(h1,'position').*[1 1 1.15 1]) 
          if pn >1 
          ax1 = gca;                   
          ax1.YAxis.Visible = 'off';   
          end
       
    end
     % text(-50,50,sprintf('Number of cells: %d', cntnotsig))
      
      ylim(ylims(nc,:))
     %  ylim([-0.01 1])
      xlim([-125 125])
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
   
     
    currfig = gcf;
    %currfig.Children(end),
    sgtitle({sprintf('%s Units subgroups behavior comparison', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
  

    legend( sprintf('Not adapting mean n= %d', cntnotsig), sprintf('Sig suppressed mean n= %d', cntsigadapt), sprintf('Sig facilitated mean n= %d', cntsigfac), sprintf('Overall Mean n= %d', cntall));
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('o_four_means_origin_data_aligned_%s_cells_cbrewer2',cellclass(nc))));
   %saveas(gcf, strcat(filename, '.svg')); 
   %saveas(gcf, strcat(filename, '.png')); 
end
%% Figure 6 - New tentative : 3D plots ; Plot 1 subplot for each peak each time including nonsig - overall mean - sig adapt - sig facilitated
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
zlims = [[0 160];[0 105];[0 210]];
%maps
nlines = 5;
cmaps = struct();


cmaps(1).map = cbrewer2('Greys', nlines);
cmaps(2).map = cbrewer2('Reds', nlines);
cmaps(3).map = cbrewer2('Blues', nlines);

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
  

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,4)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))

 clear sig_adapsu sig_incsu non_sig_su mean_sig_su
 
 cntnotsig = 0;
 cntsigfac = 0;
 cntsigadapt = 0;
 cntall = 0;
 
 non_sig_su =nan(length(aligned_trials(:,1)),length(layer_idx));
 sig_adapsu = nan(length(aligned_trials(:,1)),length(layer_idx));
 sig_facsu = nan(length(aligned_trials(:,1)),length(layer_idx));
 all_su = nan(length(aligned_trials(:,1)),length(layer_idx));
 all_mean_data = nan(4, length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
          
          mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
          all_mean_data(:,nunit) = mean_data;
          
          cntall= cntall+1;
          all_su(:,cntall) = aligned_trials(:,nunit);
        
          if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
              cntsigadapt= cntsigadapt+1;
              sig_adapsu(:,cntsigadapt) = aligned_trials(:,nunit);
              % plot(x_stim,norm_chan(:, nunit)')
              %hold on
          end
           % plot(x_stim,norm_chan(:, nunit)')
          %hold on
          
          if all_mean_data(4,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
              cntsigfac = cntsigfac+1;
              sig_facsu(:,cntsigfac) = aligned_trials(:,nunit);
          end
          
          if  pvalues(layer_idx(nunit),4) > .05
              
              cntnotsig= cntnotsig+1;
              non_sig_su(:,cntnotsig) = aligned_trials(:,nunit);
              
          end
          
      end
  end
  mean_not_sigsu = nanmean(non_sig_su,2);
  mean_sig_adapsu = nanmean(sig_adapsu,2);
  mean_sig_facsu = nanmean(sig_facsu,2);
  mean_allsu = nanmean(all_su,2);
  
 gathered_means = [ mean_allsu, mean_not_sigsu, mean_sig_adapsu, mean_sig_facsu];
 
 

  figure();  

         cmap = cmaps(nc).map ;
         cmap = cmap(5:-1:2,:);
         colormap(cmap); 
        idx = repmat(1:4,[251,1]);
         zline = nan(4,1);
        zx = nan(4,1);
        subgp = {'Overall Mean', 'Not Significant Mean', 'Sig Adapt Mean', 'Sig Facilitated Mean'};
  for ln = 1:4
           h1 = subplot(1,4,ln);
           for pn = 1:4
          % m=   mesh(-125:125,1:4,gathered_means(250*(pn-1)+1:250*pn+1,:)', 'MeshStyle','row','LineStyle', '-', 'LineWidth', 2); %'FaceAlpha', 0);
         %   m.CData = idx';
        %    m=   plot3(repmat(-125:125,[4,1])',idx, gathered_means(250*(pn-1)+1:250*pn+1,:), 'LineStyle','-','LineWidth', 2); 
            m=   plot3(-125:125, idx(:,pn), gathered_means(250*(pn-1)+1:250*pn+1,ln), 'LineStyle','-','LineWidth', 2, 'color', cmap(pn,:)); 
        hold on
           
           zline(pn) = max(gathered_means(250*(pn-1)+1:250*pn+1,ln),[],1);
           [~, xzline] = max(gathered_means(250*(pn-1)+1:250*pn+1,ln),[],1);
           zx(pn) = xzline -125; %adjusting to center the indexs on the max values
           end
           %set(gca, 'YTickLabel', sprintf('Pk%d', 1:4) )
           set(gca, 'YTick', 1:4)
           set(gca, 'YDir', 'reverse')
           ylim([1 4])
           zlim(zlims(nc,:))
          
           hold on
          
           line(zx,1:4, zline, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')
           set(gca,'box','off')
           xlabel(sprintf('%s',char(subgp(ln))));
           zlabel('Spike rate (spikes/sec)');
           view(-70,10)
          
  end
  set(gcf,'Units','inches') 
   %set(gcf,'position',[1 1 15 11])
   set(gcf,'position',[1 1 20 11])
  
filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('o_four_means_origin_data_aligned_%s_cells_3D_el10_az70',cellclass(nc))));
export_fig(gcf,strcat(filename,'.eps'));

%saveas(gcf, strcat(filename, '.eps'));    
%saveas(gcf, strcat(filename, '.svg')); 
%saveas(gcf, strcat(filename, '.png')); 
end
   

%% Other plots for Figure 6
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
nlines = 5;
cmaps = struct();

cmaps(1).map =cbrewer2('Greys', nlines);
cmaps(2).map =cbrewer2('Reds', nlines);
cmaps(3).map =cbrewer2('Blues', nlines)


%{
cmaps(1).map =cbrewer2('Blues', nlines);
cmaps(2).map =cbrewer2('Reds', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);



cmaps(1).map =cbrewer2('Purples', nlines);
cmaps(2).map =cbrewer2('Oranges', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
%}
% M ylim([0.05 .75])
% P ylim([.1 .8])
% K ylim([-0.01 1])
%normylims = [[0.05 .75];[.1 .8];[-0.01 1]];
 % P ylim([-20 85])
 % M ylim([-30 135])
 % K ylim([-20 200])
%ylims = [[-30 135];[-20 85];[-20 200]];
ylims = [[0 160];[0 105];[0 210]];

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
      sig_adapsu(:,cntsigadapt) = aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
   if all_mean_data(4,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cntsiginc= cntsiginc+1;
      sig_incsu(:,cntsiginc) = aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
 
      end
  end
 
  mean_sig_incsu = nanmean(sig_incsu,2);
  mean_sig_adapsu = nanmean(sig_adapsu,2);
 
 
  mean_origin =nanmean(aligned_trials,2);

     
     figure();
       
       cmap = flip(cmaps(nc).map) ;
      % axes('ColorOrder', cmap, 'NextPlot', 'ReplaceChildren');
       colormap(cmap); 
     
      
      h2 = subplot(1,2,1); 
      for nl = 1:4
      plot(-125:125, mean_sig_incsu(250*(nl-1)+1:250*nl+1), 'LineWidth',2,'color',cmap(nl,:));
      hold on
      end
      set(h2,'position',get(h2,'position').*[1 1 1.15 1])
      ylim(ylims(nc,:))
      %ylim([-0.01 1])
      xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('sig increase')
      
      h3 = subplot(1,2,2);
      for nl = 1:4 
      plot(-125:125, mean_sig_adapsu(250*(nl-1)+1:250*nl+1),  'LineWidth',2, 'Color',cmap(nl,:) )
      hold on
      end
      set(h3,'position',get(h3,'position').*[1 1 1.15 1])
     ylim(ylims(nc,:))
    % ylim([-0.01 1])
      xlim([-125 125])
      title('Adapting units mean')
      set(gca,'box','off')
     % legend('peak 1', 'peak 2', 'peak 3', 'peak 4')
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off';   
      
      
      
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
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\plots\',strcat(sprintf('o_sigadaptinc_stacked_origin_data_aligned_%s_cells_cbrewer2_scaled',cellclass(nc))));
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
end

%% Figure 7: Plots of the TROUGHS

%% Troughs alignment (making the data ready for plotting
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\all_units\';
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
up_dist = nan(1, length(channum),3);
max_low_dist = nan(1, length(channum));
all_locsdSUA_filtered = nan(1,length(channum),3);

for i = channum  
    if ~isempty(data_file.clean_origin_data(i).unit)
    trialidx = 1:length(data_file.clean_origin_data(i).unit(1,:));
    origin_dSUA = data_file.clean_origin_data(i).unit(401:1900,:); %- mean(data_file.clean_origin_data(i).unit(401:600,:),1);
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
    all_locsdSUA_trials = all_locsdSUA.troughs_locs(i).locs;
    
    up_dist_trials = nan(3,length(trialidx));
    clear pn
    for pn = 1:3
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
     for pn =1:3
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

%% Figure 7: Plot troughs of one single unit one after the other with 95% CI


%%store all trials of a given peak, in a matrix, across all units
%aligned_trials = nan(1204, length(layer_idx));
mean_aligned = nan(1204, 1);
stdev = nan(1204,1);
%norm_aligned_trials = nan(1204, length(layer_idx));
%norm_aligned_trials = nan(1204, length(layer_idx));

       for pn = 1:3
          % aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
           mean_aligned(250*(pn-1)+1:250*pn+1,1)= mean(suas_trials(1).aligned(max_low_dist(1)-1-124:max_low_dist(1)+125,:,pn),2);
           stdev(250*(pn-1)+1:250*pn+1,1) = squeeze(std(suas_trials(1).aligned(max_low_dist(1)-1-124:max_low_dist(1)+125,:,pn),[],2));
       end
       %normalizing with max and min of each unit
       %norm_aligned_trials(:,i) = (aligned_trials(:,i) - min(aligned_trials(:,i)))/(max(aligned_trials(:,i))-min(aligned_trials(:,i)));
     
  
  %normalizing with max and min across units
       maximum = max(mean_aligned,[],'all');
       minimum = min(mean_aligned,[],'all');
       norm_aligned_trials = (mean_aligned - minimum)/(maximum-minimum);

%aligned_trials = aligned_trials(:,~all(isnan(aligned_trials)));

%figure(); plot(1:length(aligned_trials(:,1)), aligned_trials(:,:))
%figure(); plot(1:length(aligned_trials(:,4)), norm_aligned_trials(:,:))


   h=  figure();
     
    nc =1;
      nlines = 3;
      for nl = 1:nlines  
     
      cihigh =  mean_aligned(250*(nl-1)+1:250*nl+1)+ 1.96*stdev(250*(nl-1)+1:250*nl+1)/sqrt(14);
      cilow= mean_aligned(250*(nl-1)+1:250*nl+1)- 1.96*stdev(250*(nl-1)+1:250*nl+1)/sqrt(14);
   
      h2= subplot(1,3,nl);
      plot(-125:125, mean_aligned(250*(nl-1)+1:250*nl+1), 'LineWidth',1);
      hold on
      h1= ciplot(cihigh,cilow,[-125:125],[40/255 40/255 40/255],0.1);
      ylim([0 180])
      xlim([-125 125])
      set(h1, 'edgecolor','none') 
      set(h2,'position',get(h2,'position').*[1 1 1.15 1])
      set(gca,'box','off')
      if nl == 1
            ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      end
      if nl>1
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off'; 
      end
      end

      title('Single Unit Mean')
    
      legend('Mean', '95%CI')

   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\first_year_committee_meeting\','origin_unit_example_troughs');
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
   
   %% Figure 7:  plot stacking all troughs together
   
 pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_troughs\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_troughs.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\orig_trough_values\all_units\';
troughvals = load([channeldir 'all_data_troughs']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
nlines = 5;
cmaps = struct();
%{
cmaps(1).map =cbrewer2('Blues', nlines);
cmaps(2).map =cbrewer2('Reds', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
%}
cmaps(1).map =cbrewer2('Greys', nlines);
cmaps(2).map =cbrewer2('Reds', nlines);
cmaps(3).map =cbrewer2('Blues', nlines)

ylims = [[0 140];[0 95];[0 140]];


cellclass = [ 'M', 'P', 'K'];
for nc = 1:3
layer_idx = find(strcmp(layer, cellclass(nc)));
%%store all trials of a given peak, in a matrix, across all units
%aligned_trials = nan(1204, length(layer_idx));
aligned_trials = nan(1204, length(layer_idx));

%norm_aligned_trials = nan(1204, length(layer_idx));
%norm_aligned_trials = nan(1204, length(layer_idx));

clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:3
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
 all_mean_data = nan(3, length(layer_idx));
 sig_su = nan(length(aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(troughvals.trough_vals(layer_idx(nunit)).trough)
 mean_data = nanmean(troughvals.trough_vals(layer_idx(nunit)).trough,2);
   all_mean_data(:,nunit) = mean_data;
  
      end
  end
 



   h=  figure();
     mean_origin =nanmean(aligned_trials,2);
 
        cmap = flip(cmaps(nc).map) ;
  
       colormap(cmap); 
     
      
      nlines = 3;
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
     
      set(gca,'box','off')
      legend('trough 1', 'trough 2', 'trough 3')
      % P ylim([-20 85])
      % M ylim([-35 135])
      % K ylim([-20 200])
       
      % M ylim([0.05 .75])
      % P ylim([.1 .8])
      % K ylim([-0.01 1])
      
       % K norm ylim([-0.02 1.1])
      
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
 

    currfig = gcf;
    
    title(currfig.Children(end),{sprintf('%s mean responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\plots\',strcat(sprintf('v_common_origin_data_aligned_%s_cells_cmap',cellclass(nc))));
   %saveas(gcf, strcat(filename, '.svg')); 
   %saveas(gcf, strcat(filename, '.png')); 
end

%%  Plotting the median of the troughs

  pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_troughs\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_troughs.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\orig_trough_values\all_units\';
troughvals = load([channeldir 'all_data_troughs']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
nlines = 5;
cmaps = struct();
%{
cmaps(1).map =cbrewer2('Blues', nlines);
cmaps(2).map =cbrewer2('Reds', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
%}
cmaps(1).map =cbrewer2('Greys', nlines);
cmaps(2).map =cbrewer2('Reds', nlines);
cmaps(3).map =cbrewer2('Blues', nlines)

ylims = [[0 140];[0 95];[0 140]];


cellclass = [ 'M', 'P', 'K'];
for nc = 1:3
layer_idx = find(strcmp(layer, cellclass(nc)));
%%store all trials of a given peak, in a matrix, across all units
%aligned_trials = nan(1204, length(layer_idx));
aligned_trials = nan(1204, length(layer_idx));

%norm_aligned_trials = nan(1204, length(layer_idx));
%norm_aligned_trials = nan(1204, length(layer_idx));

clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:3
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
 all_median_data = nan(3, length(layer_idx));
 sig_su = nan(length(aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(troughvals.trough_vals(layer_idx(nunit)).trough)
 median_data = nanmedian(troughvals.trough_vals(layer_idx(nunit)).trough,2);
   all_median_data(:,nunit) = median_data;
  
      end
  end
 



   h=  figure();
     median_origin =nanmedian(aligned_trials,2);
 
        cmap = flip(cmaps(nc).map) ;
  
       colormap(cmap); 
     
      
      nlines = 3;
      for nl = 1:nlines
      %cmap = jet(4); 
   
      plot(-125:125, median_origin(250*(nl-1)+1:250*nl+1), 'LineWidth',1,'color',cmap(nl,:));
      hold on
      end
      set(h,'position',get(h,'position').*[1 1 1.15 1])
      ylim(ylims(nc,:))
      xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('Overall Median')
     
      set(gca,'box','off')
      legend('trough 1', 'trough 2', 'trough 3')
      % P ylim([-20 85])
      % M ylim([-35 135])
      % K ylim([-20 200])
       
      % M ylim([0.05 .75])
      % P ylim([.1 .8])
      % K ylim([-0.01 1])
      
       % K norm ylim([-0.02 1.1])
      
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
 

    currfig = gcf;
    
    title(currfig.Children(end),{sprintf('%s median responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\plots\',strcat(sprintf('v_common_median_origin_data_aligned_%s_cells_cmap',cellclass(nc))));
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
end



%% Mean of the troughs with ci

 
      
      
 pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_troughs\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_troughs.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\orig_trough_values\all_units\';
troughvals = load([channeldir 'all_data_troughs']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
nlines = 5;
cmaps = struct();
%{
cmaps(1).map =cbrewer2('Blues', nlines);
cmaps(2).map =cbrewer2('Reds', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
%}
cmaps(1).map =cbrewer2('Greys', nlines);
cmaps(2).map =cbrewer2('Reds', nlines);
cmaps(3).map =cbrewer2('Blues', nlines)

ylims = [[0 140];[0 95];[0 140]];


cellclass = [ 'M', 'P', 'K'];
for nc = 1:3
layer_idx = find(strcmp(layer, cellclass(nc)));
%%store all trials of a given peak, in a matrix, across all units
%aligned_trials = nan(1204, length(layer_idx));
aligned_trials = nan(1204, length(layer_idx));

%norm_aligned_trials = nan(1204, length(layer_idx));
%norm_aligned_trials = nan(1204, length(layer_idx));

clear i
for i = 1:length(layer_idx)
    if ~isnan(max_low_dist(layer_idx(i)))
       for pn = 1:3
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
 all_mean_data = nan(3, length(layer_idx));
 sig_su = nan(length(aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(troughvals.trough_vals(layer_idx(nunit)).trough)
 mean_data = nanmean(troughvals.trough_vals(layer_idx(nunit)).trough,2);
   all_mean_data(:,nunit) = mean_data;
  
      end
  end
 



   h=  figure();
     mean_origin =nanmean(aligned_trials,2);
 
        cmap = flip(cmaps(nc).map) ;
  
       colormap(cmap); 
     
      
      nlines = 3;
      for nl = 1:nlines
      %cmap = jet(4); 
      cihigh =  all_mean_data(250*(nl-1)+1:250*nl+1)+ 1.96*stdev(250*(nl-1)+1:250*nl+1)/sqrt(length(all_mean_data(1,:)));
      cilow= all_mean_data(250*(nl-1)+1:250*nl+1)- 1.96*stdev(250*(nl-1)+1:250*nl+1)/sqrt(length(all_mean_data(1,:)));
   
   
      
      plot(-125:125, mean_origin(250*(nl-1)+1:250*nl+1), 'LineWidth',1,'color',cmap(nl,:));
      hold on
        
      h1= ciplot(cihigh,cilow,[-125:125],[40/255 40/255 40/255],0.1);
      
      end
      set(h,'position',get(h,'position').*[1 1 1.15 1])
      ylim(ylims(nc,:))
      xlim([-125 125])
      set(gca,'box','off')
      ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
      title('Overall Mean')
     
      set(gca,'box','off')
      legend('trough 1', 'trough 2', 'trough 3')
      % P ylim([-20 85])
      % M ylim([-35 135])
      % K ylim([-20 200])
       
      % M ylim([0.05 .75])
      % P ylim([.1 .8])
      % K ylim([-0.01 1])
      
       % K norm ylim([-0.02 1.1])
      
      %set(gca, 'linewidth',2)
      %hold on
      %plot([0 0], ylim,'k')
 

    currfig = gcf;
    
    title(currfig.Children(end),{sprintf('%s mean responses', cellclass(nc))}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\plots\',strcat(sprintf('v_common_origin_data_aligned_%s_cells_cmap',cellclass(nc))));
   %saveas(gcf, strcat(filename, '.svg')); 
   %saveas(gcf, strcat(filename, '.png')); 
end

