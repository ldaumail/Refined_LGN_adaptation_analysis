newdatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020\all_units\';
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
   
   %%%%%%% align mean single units to peak of interest%%%%
      %start finding the peaks at the first non NaN location
      
        %%detect peaks locations for data aligned on peak n
        %{
      for pn = 1:4
       for len1 = 1:800
       if eval(['~isnan(mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(len1))'])
          break
       end
      end
      for len2 = len1:1550
          if eval(['~all(isnan(mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) ')) && mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(len2) < mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(len2+1)'])
          eval(['locsdSUA_filtered_pk' num2str(pn) ' = findpeaks(mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(len2:end));'])
          end
          break
      end
      end
     
      
      clear pn
      for pn = 1:4
          %detect peaks locations on data aligned on peak n
          clear len1
        for len1 =1:800
           if eval(['~isnan(mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(len1))'])
                break
           end
        end
  
        clear len2
        for len2 = len1:1550
            if eval(['~all(isnan(mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) ')) && mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(len2) < mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(len2+1)'])
             
             eval(['locsdSUA_filtered = findpeaks(mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(len2:end));'])
            
            % clear shiftidx ;
               if pn == 1
                  if length(locsdSUA_filtered.loc) == 4 || eval(['mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(locsdSUA_filtered.loc(1)+len2) >= 0.4*mean_filtered_dSUA(i).mean_peakaligned2(locsdSUA_filtered_pk2.loc(2)+len2)']) 
     
                  %store first peak location 
                  all_locsdSUA_filtered(:,i,pn) = locsdSUA_filtered.loc(1)+len2;
              
                  elseif length(locsdSUA_filtered.loc) > 4 && eval(['mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) '(locsdSUA_filtered.loc(1)+len2) < 0.4*mean_filtered_dSUA(i).mean_peakaligned2(locsdSUA_filtered_pk2.loc(2)+len2)'])

                  all_locsdSUA_filtered(:,i,pn) = locsdSUA_filtered_pk2.loc(2)+len2;    
            
                  end
               end
               if pn > 1 && (length(locsdSUA_filtered.loc) == 4 || eval('mean_filtered_dSUA(i).mean_peakaligned1(locsdSUA_filtered_pk1.loc(1)+len2) >= 0.4*mean_filtered_dSUA(i).mean_peakaligned2(locsdSUA_filtered_pk2.loc(2)+len2)'))
               
               all_locsdSUA_filtered(:,i,pn) = locsdSUA_filtered.loc(pn)+len2;
               elseif pn > 1 && (length(locsdSUA_filtered.loc) > 4 && eval('mean_filtered_dSUA(i).mean_peakaligned1(locsdSUA_filtered_pk1.loc(1)+len2) < 0.4*mean_filtered_dSUA(i).mean_peakaligned2(locsdSUA_filtered_pk2.loc(2)+len2)'))
               
               eval(['all_locsdSUA_filtered(:,i,pn) = locsdSUA_filtered_pk' num2str(pn+1) '.loc(pn+1)+len2;'])
               end
             break 
            end 
        end

    %compute the distance between the first peak and the last datapoint and store  
   %in a matrix
   eval(['up_dist(:,i,pn)= length(mean_filtered_dSUA(i).mean_peakaligned' num2str(pn) ')- all_locsdSUA_filtered(:,i,pn);'])
  
      end
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
        filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020\plots\',strcat( sprintf('trials_bscorr_origin_data_aligned_unit_%d', i)));
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
aligned_trials = nan(1000, length(layer_idx));
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
      plot(1:1000, mean_origin,'LineWidth',1, 'Color',[49/255 130/255 189/255])
      hold on
      ci_low = mean_origin - 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(aligned_trials(1,:)));
      plot(1:1000, ci_low,':', 'LineWidth',.7,'Color', [49/255 130/255 189/255])
      hold on
      ci_high = mean_origin + 1.96*std(aligned_trials,0,2, 'omitnan')./sqrt(length(aligned_trials(1,:)));
      plot(1:1000, ci_high,':', 'LineWidth',.7,'Color', [49/255 130/255 189/255])
 %blue [49/255 130/255 189/255]
 %black [222/255 45/255 38/255]
 %orange [253/255 174/255 107/255]
 %brown [165/255 42/255 42/255]
 
      set(gca, 'linewidth',2)
      set(gca,'box','off')

    title({'K mean responses'}, 'Interpreter', 'none')
    xlabel('\fontsize{14}Resolution (ms)')
    ylh = ylabel({'\fontsize{14}Spike Rate (spikes/s)'});
   
   set(gcf,'Units','inches') 
   %set(gcf,'position',[1 1 8.5 11])
   set(gcf,'position',[1 1 15 11])
   
   filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020\plots\',strcat('95ci_nanmean_bscorr_origin_data_aligned_Kcells_blue'));
   saveas(gcf, strcat(filename, '.png'));
   