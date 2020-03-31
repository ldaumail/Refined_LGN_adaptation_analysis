
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

   %% plot with significant adapting units and overall mean in two separate sublots, stacking all peaks together
   
 pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_troughs\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_troughs.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\orig_trough_values\all_units\';
troughvals = load([channeldir 'all_data_troughs']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

layer_idx = find(strcmp(layer, 'K'));
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
 sig_su = nan(length(norm_aligned_trials(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
      if ~isempty(troughvals.trough_vals(layer_idx(nunit)).trough)
 mean_data = nanmean(troughvals.trough_vals(layer_idx(nunit)).trough,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(3,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),3) < .05
      cnt= cnt+1;
      sig_su(:,cnt) = norm_aligned_trials(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
      end
  end
  mean_sig_su = nanmean(sig_su,2);



     figure();
     mean_origin =nanmean(norm_aligned_trials,2);
  %red [215/255 25/255 28/255] 
 %orange [253/255 174/255 97/255] %brown [121/255 96/255 76/255]; 
 %blue [44/255 123/255 182/255]

 %brown [165/255 42/255 42/255]
     col1 = [121/255 96/255 76/255]; 
     col2 =[44/255 123/255 182/255];
     
      R = linspace(col1(1), col2(1),4);
      G = linspace(col1(2), col2(2),4);
      B = linspace(col1(3), col2(3),4);
      cmap = [R' G' B'];
      colormap(cmap); 
      
      h = subplot(1,2,1);
      
      nlines = 3;
      for nl = 1:nlines
      %cmap = jet(4); 
   
      plot(-125:125, mean_origin(250*(nl-1)+1:250*nl+1), 'LineWidth',1,'color',cmap(nl,:));
      hold on
      end
      set(h,'position',get(h,'position').*[1 1 1.15 1])
      ylim([-0.01 1])
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
      ylim([-0.01 1])
      xlim([-125 125])
      title('Adapting units mean')
      set(gca,'box','off')
      legend('peak 1', 'peak 2', 'peak 3', 'peak 4')
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
     % if pn >1 
      %ax1 = gca;                   
      %ax1.YAxis.Visible = 'off';   
     % end

    currfig = gcf;
    
    title(currfig.Children(end),{'K mean responses'}, 'Interpreter', 'none', 'FontSize', 20)
   
   % xlabel('\fontsize{14}Resolution (ms)')
    
   
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 15 11])
   
   
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\plots\',strcat('sigadapt_stacked_commonnorm_nanmean_bscorr_origin_data_aligned_Kcells_cmap_corrected'));
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 