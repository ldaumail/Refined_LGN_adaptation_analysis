%this script was written after get_clean_peaks_and_data.m in order to
%analyze the lmer results and plot the clean data.
%Written by Loic Daumail edited on 1/17/2022

%loading the clean data
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_SUA_sup_50']; 
data_file = load(channelfilename);


%exclude 160517, (first unit, left empty, it is a K neuron)
%Reject 180806 p1 uclust17, M cell, as doesn't seem well triggered (46)
%Reject 181207 (B) uclust22, M cell, as doesn't seem well triggered (55)
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 
 pvaluesdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
 pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
 
%% Rough plots of peaks with pvalues (not very representative, as mean unit activity)
  channum = 1: length(data_file.clean_high_SUA);

 %for n = 1:3 
for chan = 1:12:length(channum)
h = figure;
xabs = -199:1300;
idx = [1 3 5 7 9 11 2 4 6 8 10 12];
nyq = 500;
all_mean_data = nan(length(xabs), length(1:12));
clear i ;
 for i = 1:12
 
   mean_data = mean(data_file.clean_high_SUA(chan+i-1).namelist(1:1500,:),2);
   
  
   lpc       = 4.5; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   
  if ~all(isnan(mean_data))
   lpsu      = filtfilt(bwb,bwa, mean_data);
   
    sp = subplot(length(1:6), 2, idx(i));
    plot(xabs, lpsu)
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    if i == length(6)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}Spike Rate (spikes/s)'});
    end
   if i < 6 || (i >= 7)&&(i < 12)
        set(sp, 'XTick', [])
   end
      ylabelh = text(max(xabs), mean(lpsu,1), strcat(num2str(chan+i-1),' | ', layer(chan+i-1)),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 10);
   for npeak = 1:4
         for len = 231:480 %from 200 + 30 (lgn response onset) 
            if lpsu(len) < lpsu(len+1)
   locs = findpeaks(lpsu(len:1450));
        break
            end
         end
         
         if length(locs.loc) >= 4
             %adjust location to the first data point of lpsu (+len), then adjust
             %to xabs (-200)
   xlocation = locs.loc(npeak)+len-200;
         end 
            
   text(xlocation, mean(lpsu,1), strcat(num2str(sprintf('%.2f', pvalues(chan+i-1,npeak)))),'HorizontalAlignment','center','FontName', 'Arial','FontSize', 7);
   end 

 end
 end
      set(gca, 'linewidth',2)
      set(gca,'box','off')

    sgtitle({f{2}, 'all good responses, p<0.05, associated to adaptation pvalues'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
    %filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\plots\',strcat(f{2}, sprintf('x_%d_better_raw_data_peakspvalues_2dec', chan+i-1)));
    %saveas(gcf, strcat(filename, '.png'));

end

%%

 %% compute proportion of significant adaptation per peak and proportion of neurons adapting for a certain amount of 
%peak from peak 2 to 4

channeldir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
pvaluesdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected_dunnett.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

peakvals = load([channeldir 'all_raw_data_peaks']);

   cnt =0;
   for i =1:length(peakvals.peak_vals)
       if ~isempty(peakvals.peak_vals(i).peak)
           cnt = cnt+1;
       end
        
   end
 
%exclude 160517, (first unit, left empty, it is a K neuron)
%Reject 180806 p1 uclust17, M cell, as doesn't seem well triggered (46)
%Reject 181207 (B) uclust22, M cell, as doesn't seem well triggered (55)
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

 layer_idx = find(strcmp(layer, 'K'));
 
 
%{
 nan_layer_idx = layer_idx(isnan(pvalues(layer_idx,1)));
 cnt = 0;
 clear empt
  for i =1:length(layer_idx)
       if isempty(peakvals.peak_vals(layer_idx(i)).peak)
           cnt = cnt+1;
           empt(cnt) = layer_idx(i);
       end
        
   end
 %}
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};

 all_locs = nan(4,length(layer_idx));
 all_pks = nan(4,length(layer_idx));
 all_mean_data = nan(4, length(layer_idx));

  
 cntpk2 = 0;
 cntpk3 = 0;
 cntpk4 = 0;
 cntpk2pk3 = 0;
 cntpk2pk3pk4 = 0;
 cntpk3pk4 =0;
 
 cntincpk2 =0;
 cntincpk3 =0;
 cntincpk4 =0;
 
 cntnspk2 =0;
 cntnspk3 =0;
 cntnspk4=0;
 
  for nunit = 1:length(layer_idx)
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
 mean_data = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
   
   all_mean_data(:,nunit) = mean_data;

     if all_mean_data(2,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),1) < .05
         cntpk2 = cntpk2 +1;
     end
     if all_mean_data(3,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),2) < .05
         cntpk3 = cntpk3 +1;
     end
     if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cntpk4 = cntpk4 +1;
     end
     if all_mean_data(2,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),1) < .05 && ...
             all_mean_data(3,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),2) < .05
         cntpk2pk3 = cntpk2pk3 +1;
     end
     if all_mean_data(2,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),1) < .05 && ...
             all_mean_data(3,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),2) < .05 ...
             && all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cntpk2pk3pk4 = cntpk2pk3pk4 +1;
     end
       if all_mean_data(3,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),2) < .05 ...
             && all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),3) < .05 
         cntpk3pk4 = cntpk3pk4 +1;
       end
       
     if all_mean_data(2,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),1) < .05
         cntincpk2 = cntincpk2 +1;
     end
     if all_mean_data(3,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),2) < .05
         cntincpk3 = cntincpk3 +1;
     end
     
     if all_mean_data(4,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cntincpk4 = cntincpk4 +1;
     end
     if pvalues(layer_idx(nunit),1) > .05
         cntnspk2 = cntnspk2 +1;
     end
     if pvalues(layer_idx(nunit),2) > .05
         cntnspk3 = cntnspk3 +1;
     end
    if pvalues(layer_idx(nunit),3) > .05
         cntnspk4 = cntnspk4 +1;
    end
      end
  end
 all_mean_data = all_mean_data(:, ~all(isnan(all_mean_data)));
 
 percentpk2 = cntpk2*100/length(all_mean_data(1,:));
 percentpk3 = cntpk3*100/length(all_mean_data(1,:)); 
 percentpk4 = cntpk4*100/length(all_mean_data(1,:));
 
 percentpk2pk3 = cntpk2pk3*100/length(all_mean_data(1,:));
 percentpk2pk3pk4 = cntpk2pk3pk4*100/length(all_mean_data(1,:));
 percentpk3pk4 = cntpk3pk4*100/length(all_mean_data(1,:));
 
 percentincpk2 = cntincpk2*100/length(all_mean_data(1,:));
 percentincpk3 = cntincpk3*100/length(all_mean_data(1,:));
 percentincpk4 = cntincpk4*100/length(all_mean_data(1,:));
 
 percentnspk2 = cntnspk2*100/length(all_mean_data(1,:));
 percentnspk3 = cntnspk3*100/length(all_mean_data(1,:));
 percentnspk4 = cntnspk4*100/length(all_mean_data(1,:));
 
 ncells =length(all_mean_data(1,:));
 
  %% compute proportion of significant adaptation per trough and proportion of neurons adapting for a certain amount of 
%trough from trough 2  to 3

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_troughs_03032020\orig_trough_values\all_units\';
pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_troughs\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_troughs.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

troughvals = load([channeldir 'all_data_troughs']);

 
%exclude 160517, (first unit, left empty, it is a K neuron)
%Reject 180806 p1 uclust17, M cell, as doesn't seem well triggered (46)
%Reject 181207 (B) uclust22, M cell, as doesn't seem well triggered (55)
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

 layer_idx = find(strcmp(layer, 'K'));
%{
 nan_layer_idx = layer_idx(isnan(pvalues(layer_idx,1)));
 cnt = 0;
 clear empt
  for i =1:length(layer_idx)
       if isempty(peakvals.peak_vals(layer_idx(i)).peak)
           cnt = cnt+1;
           empt(cnt) = layer_idx(i);
       end
        
   end
 %}
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};

 all_locs = nan(3,length(layer_idx));
 all_trghs = nan(3,length(layer_idx));
 all_mean_data = nan(3, length(layer_idx));

  
 cntt2 = 0;
 cntt3 = 0;
 
 cntt2t3 = 0;
 
 
 cntinct2 =0;
 cntinct3 =0;

 cntnst2 =0;
 cntnst3 =0;
 
 
  for nunit = 1:length(layer_idx)
      if ~isempty(troughvals.trough_vals(layer_idx(nunit)).trough)
 mean_data = nanmean(troughvals.trough_vals(layer_idx(nunit)).trough,2);
   
   all_mean_data(:,nunit) = mean_data;

     if all_mean_data(2,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),2) < .05
         cntt2 = cntt2 +1;
     end
     if all_mean_data(3,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cntt3 = cntt3 +1;
     end
     
     if all_mean_data(2,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),2) < .05 && ...
             all_mean_data(3,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cntt2t3 = cntt2t3 +1;
     end
       
     if all_mean_data(2,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),2) < .05
         cntinct2 = cntinct2 +1;
     end
     if all_mean_data(3,nunit) > all_mean_data(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cntinct3 = cntinct3 +1;
     end
     
     if pvalues(layer_idx(nunit),2) > .05
         cntnst2 = cntnst2 +1;
     end
     if pvalues(layer_idx(nunit),3) > .05
         cntnst3 = cntnst3 +1;
     end
   
      end
  end
 all_mean_data = all_mean_data(:, ~all(isnan(all_mean_data)));
 
 percentt2 = cntt2*100/length(all_mean_data(1,:));
 percentt3 = cntt3*100/length(all_mean_data(1,:)); 
 
 percentt2t3 = cntt2t3*100/length(all_mean_data(1,:));
 
 percentinct2 = cntinct2*100/length(all_mean_data(1,:));
 percentinct3 = cntinct3*100/length(all_mean_data(1,:));

 percentnst2 = cntnst2*100/length(all_mean_data(1,:));
 percentnst3 = cntnst3*100/length(all_mean_data(1,:));

 ncells =length(all_mean_data(1,:));
 