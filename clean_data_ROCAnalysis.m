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
 for i = 1:length(channum)
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
 end

%plot the mean data only in the 5Hz range 

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
layer_idx = find(strcmp(layer, 'M'));
figure, 
normspec = (nanmean(mean_S(:,:,layer_idx),3) - min(nanmean(mean_S(:,:,layer_idx),3)))./(max(nanmean(mean_S(:,:,layer_idx),3)) - min(nanmean(mean_S(:,:,layer_idx),3)));
%time_adj_data = nan(length(tvec),1);
%time_adj_data(129:end,1) = normspec(:,1);
plot(tvec,squeeze(normspec(:,1))')
%spec = nanmean(mean_S,3);
%imagesc(tvec,sort(f),spec')
%ylim([3.9 60]); 
%set(gca,'ydir','normal')
      
title({'Mean power at 4 Hz', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus onset(ms)')
    ylabel('Power (Normalized)')
filename = strcat('C:\Users\daumail\Documents\first_year_committee_meeting\mean_power_4hz');
saveas(gcf, strcat(filename, '.png')); 
saveas(gcf, strcat(filename, '.svg')); 
figure, 
spec = nanmean(mean_S,3);
imagesc(tvec,sort(f),mean_S(:,:,1)')
ylim([2 20]); 
set(gca,'ydir','normal')
title({'Mean spectrogram', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus onset(ms)')
    ylabel('Frequency band (Hz)')
filename = strcat('C:\Users\daumail\Documents\first_year_committee_meeting\mean_powersepc_imagesc');
saveas(gcf, strcat(filename, '.png')); 
saveas(gcf, strcat(filename, '.svg')); 

%% perform the Receiver Operating Characteristics analysis

reps   = 10000;
all_sigs95 = nan(length(Ses),1);
all_sigs90 = nan(length(Ses),1);
for i = 1:length(Ses)
    if ~isempty(Ses(i).namelist2)
part1 = nanmean(squeeze(Ses(i).namelist2(1:575,1,:)), 1);
part2 = nanmean(squeeze(Ses(i).namelist2(576:1150,1,:)), 1);

    if nanmean(part1) > nanmean(part2)
    cond1               = part1;
    cond2               = part2;
    else
    cond1               = part2;
    cond2               = part1;
    end
  %{  
figure(); boxplot([cond1' cond2'],'notch','on','labels',{'cond1','cond2'}); hold on    
x=ones(length(cond1),1).*(1+(rand(length(cond1),1)-0.5)/5);
x1=ones(length(cond1),1).*(1+(rand(length(cond1),1)-0.5)/10);
f1=scatter(x,cond1,'k','filled');f1.MarkerFaceAlpha = 0.4;hold on
f2=scatter(x1*2,cond2,'k','filled');f1.MarkerFaceAlpha = 0.4;
ylim([0 400])
ylabel('Power (spikes^2/s^2)')
filename = strcat('C:\Users\daumail\Documents\first_year_committee_meeting\distributions_cond1_cond2');
saveas(gcf, strcat(filename, '.png')); 
saveas(gcf, strcat(filename, '.svg')); 
  
%}
    [X,Y,T,AUC]           = perfcurve([ones(length(cond1),1); repmat(2,length(cond2),1)],[cond1 cond2],1);
    NP                    = length(cond2);
    PR                    = length(cond1);
    catdat                = [cond1 cond2];


    for r       = 1:reps
       clear shufNP shufPR
       shufPR         = catdat(randperm(length(catdat),PR));
       shufNP         = catdat(randperm(length(catdat),NP));
       [~,~,~,...
       shufAUC(r)]    = perfcurve([ones(PR,1); repmat(2,NP,1)],[shufPR shufNP],1);
    end

    critT95         = quantile(shufAUC,.95);
    critT90         = quantile(shufAUC,.90);

    if AUC > critT95
       sig95          = 1;
    else
       sig95          = 0;
    end

    all_sigs95(i) = sig95;
    
    if AUC > critT90
       sig90        = 1;
    else
       sig90         = 0;
    end
    all_sigs90(i) = sig90;
    end
end

save( strcat(newdatadir,'roc_results95_stimonset_to1150ms.mat'), 'all_sigs95');

sig95_idx = find(all_sigs95);
sig90_idx = find(all_sigs90);


%% Analyze ROC analysis results per layer

newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_origin_sup_50']; 
data_file = load(channelfilename);

sig95_idx = load( strcat(newdatadir,'roc_results95_stimonset_to1150ms.mat'));
 
layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
layer_idx = find(strcmp(layer, 'M'));
log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
Ses = struct();
bs_data = struct();
%channum = 1: length(log_p_layer);
%add +128 ms to adjust the time to -100ms before stim onset
mean_S = nan(1174,38, length(layer_idx));
Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
xabs = -100:1301;
xpow = 0:1173;

idx = [1 3 2 4];
channum = 1: length(data_file.clean_origin_data);

clear i
clear chan
for chan =1:4:length(channum)-4
 figure(); 
    for i = 1:4
   if log_p_layer(chan+i-1) ~= 0 && ~isnan(sig95_idx.all_sigs95(chan+i-1)) 
   data = data_file.clean_origin_data(chan+i-1).unit(401:1900,:);
   bsl = mean(data(1:200,:));
   
   norm_mean_bs = data(72:end, :) - bsl;
 
   clear S namelist;
   [S,t,f]        = mtspecgramc(norm_mean_bs(:,:) ,movingwin, params); 
 

   mean_S(:,:,chan+i-1) = nanmean(S,3);
%tvec     = t*1000 + (xabs(1));
   normchan = (mean_S(:,:,chan+i-1) - min(mean_S(:,:,chan+i-1)))./(max(mean_S(:,:,chan+i-1)) - min(mean_S(:,:,chan+i-1)));
   xpow = 0:1173; 
   sp = subplot(length(1:2), 2, idx(i) );
   plot(xpow,normchan(:,1))
   hold on
   plot([0 0], ylim,'k')
   hold on
   plot([1150 1150], ylim,'k')
   xlim([-100 1174])
     if i == length(2)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}Normalized Power at 4Hz (no units)'});
     end
    

      ylabelh = text(mean(xpow)/2, max(normchan(:,1))+0.05, strcat(num2str(chan+i-1),' | ROC sign ', num2str(sig95_idx.all_sigs95(chan+i-1))),'HorizontalAlignment','left','FontName', 'Arial', 'Interpreter','none','FontSize', 10);
      
      set(gca, 'linewidth',2)
      set(gca,'box','off')
   end
    end 
    sgtitle({'DE50_NDE0_su', 'Mean M layer SUA power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
   xlabel('Time from stimulus onset(ms)')
   % ylabel('Normalized Power at 4Hz(no units)')
   
  % filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020\power_plots\',strcat('DE50_NDE0_', sprintf('mean_sua_power_4hz_M_layer_%d_stimonset_to1150ms', chan)));
   %saveas(gcf, strcat(filename, '.png')); 
   
end


%% Proportions of significant units per layer
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_origin_sup_50']; 
data_file = load(channelfilename);

sig95_idx = load( strcat(newdatadir,'roc_results95_stimonset_to1150ms.mat'));
 
layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
layer_idx = find(strcmp(layer, 'M'));
log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);




%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];

cntdecrease = 0;
cntincrease = 0;
ncells =0;
clear i ;
 for i = 1:length(layer_idx)
     if ~isempty(data_file.clean_origin_data(layer_idx(i)).unit) && ~isnan(sig95_idx.all_sigs95(layer_idx(i)))
data = squeeze(data_file.clean_origin_data(layer_idx(i)).unit(401:1900,:,:));
ncells = ncells+1;
   bsl = mean(data(1:199,:));
   
   norm_mean_bs = data(72:end, :) - bsl;
   
clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs ,movingwin, params); 
 

part1 = nanmean(nanmean(squeeze(S(1:575,:)), 1),2);

part2 = nanmean(nanmean(squeeze(S(576:1150,:)), 1),2);


if part1 > part2 && sig95_idx.all_sigs95(layer_idx(i)) ==1
    cntdecrease = cntdecrease +1;
end

if part1 < part2 && sig95_idx.all_sigs95(layer_idx(i)) ==1
    cntincrease = cntincrease +1;
end
     end
end
 
 percentdecrease = cntdecrease*100./ncells;
 percentincrease = cntincrease*100./ncells;



%% Plot mean with error bars before stimulation and during stimulation in the same analysis
%% different way to normalize the data(normalize mean SUA before computing the grand cell class mean
%% plotting with spiking activity significant changes
%{
pvaluesdir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\lmer_results\';
 pvalfilename = [pvaluesdir 'lmer_results.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\';
 peakvals = load([channeldir 'all_data_peaks']);

layer_idx = find(strcmp(layer, 'P'));
log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
 
Ses = struct();
bs_data = struct();
channum = 1: length(log_p_layer);
mean_S_stim = nan(1646+128,38, length(channum));
%compute the power spectrum
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
 

clear i ;
 for i = 1:length(channum)
data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(1:1901,:,:));
   bsl = mean(data(400:599,:));
   %stim and bl data
   norm_mean_bs = nan(length(data(:,1)),1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(1:end,:,:) - bsl;
  

clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 

mean_S_stim(129:end,:,i) = nanmean(S,3);

%we can also store tvec and f in a struct, but they are all identical
 end
 
time_adj = -99:28;
x_stim = cat(2, time_adj-500 , t*1000 -600) ;

%here we compute the individual normalized units necessary for the variance
%for both the baseline data and the stimulus data
norm_chan = nan(length(mean_S_stim(:,1,1)), length(layer_idx));
clear i;
for i = 1:length(layer_idx)
min_chan =min(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
max_chan = max(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
norm_chan(:,i) = (squeeze(mean_S_stim(:,1,layer_idx(i)))-min_chan)./(max_chan - min_chan);
end

normspec = nanmean(norm_chan,2);

figure, 
 plot(x_stim,normspec', 'LineWidth',1, 'Color',[229/255, 49/255, 90/255])
 xlim([-600 1250])
 ylim([-0.1 1])
 %green[167/255 185/255 54/255])
 %black = [24/255 23/255 23/255] )
 %pink = [229/255, 49/255, 90/255]) 
 hold on
 ci_low = normspec(:,1) - 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_low,':', 'LineWidth',1,'Color', [.40 .40 .40])
 hold on
 ci_high = normspec(:,1) + 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_high,':', 'LineWidth',1,'Color', [.40 .40 .40])
 plot([0 0], ylim,'k')
 hold on
 plot([1150 1150], ylim,'k')
 hold on 
 
 cnt = 0;
 all_mean_data = nan(4, length(layer_idx));
  for nunit = 1:length(layer_idx)
 mean_data = nanmean(peakvals.data_peaks(layer_idx(nunit)).namelist,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cnt= cnt+1;
      sig_su(:,cnt) = norm_chan(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
  
  end
  mean_sig_su = mean(sig_su,2);
  plot(x_stim, mean_sig_su,  'LineWidth',1)
  
 %ylim([-0.8 1.2])
     set(gca, 'linewidth',2)
      set(gca,'box','off') 
      xlabel('Time from stimulus onset(ms)')
   ylabel('Normalized Power at 4Hz(no units)')
title({'P class cells mean power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    legend('Mean', 'Mean-1.96*sem', 'Mean+1.96*sem', 'Mean significant decrease', 'Location', 'bestoutside')
    
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum\plots\',contrast{2},'indiv_normalized_power_freq_time_mean_95ci_P_layer_4hz_gathered_pink_sig_suamean');
saveas(gcf, strcat(filename, '.svg'));
saveas(gcf, strcat(filename, '.png'));
%export_fig(gcf, '-jpg', '-transparent');
%}

%% Plot mean with error bars before stimulation and during stimulation in the same analysis
%% different way to normalize the data(normalize mean SUA before computing the grand cell class mean
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
red = [215/255 25/255 28/255]; 
orange = [253/255 174/255 97/255];
blue = [44/255 123/255 182/255];
colors = [red; orange; blue]; 

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

figure, 
 plot(x_stim,normspec', 'LineWidth',1, 'Color',colors(nc,:))
 xlim([-600 1250])
 ylim([-0.1 1.05])
 
 hold on
 ci_low = normspec(:,1) - 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_low,':', 'LineWidth',1,'Color',colors(nc,:))
 hold on
 ci_high = normspec(:,1) + 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_high,':', 'LineWidth',1,'Color',colors(nc,:))
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
  
 %ylim([-0.8 1.2])
     set(gca, 'linewidth',2)
      set(gca,'box','off') 
      xlabel('Time from stimulus onset(ms)')
   ylabel('Normalized Power at 4Hz(no units)')
title({sprintf('%s class cells mean power at 4Hz vs time normalized', cellclass(nc))}, 'Interpreter', 'none')
    legend('Mean', 'Mean-1.96*sem', 'Mean+1.96*sem', 'Mean significant decrease', 'Location', 'bestoutside')
    
filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\power_plots\',contrast{2},sprintf('indiv_normalized_power_freq_time_mean_95ci_%s_layer_4hz_gathered_orange_sig_suamean_pow',cellclass(nc)));
%saveas(gcf, strcat(filename, '.svg'));
%saveas(gcf, strcat(filename, '.png'));
end
%export_fig(gcf, '-jpg', '-transparent');

%% save part1 and part2 for R plotting
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_origin_sup_50']; 
data_file = load(channelfilename);


layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

Ses = struct();
bs_data = struct();
channum = 1: length(layer);
mean_S_stim = nan(1646+128,38, length(channum));
S_stim = struct();
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
S_stim(i).S = nan(1646+128,1, length(S(1,1,:)));
S_stim(i).S(129:end,:) = squeeze(S(:,1,:));
      end

%we can also store tvec and f in a struct, but they are all identical
 end
 
%time_adj = -599:-472;
%x_stim = cat(2, time_adj , t*1000 -600) ;
 

%here we compute the individual normalized units necessary for the variance
%and the data that we are about to plot
%for both the baseline data and the stimulus data
norm_chan = nan(length(mean_S_stim(:,1,1)), length(channum));
clear i;
for i = channum
min_chan =min(squeeze(mean_S_stim(:,1,i)),[],1);
max_chan = max(squeeze(mean_S_stim(:,1,i)),[],1);
norm_chan(:,i) = (squeeze(mean_S_stim(:,1,i))-min_chan)./(max_chan - min_chan);
end

%{
part1 = nan(length(channum),1);
part2 = nan(length(channum),1);

  for nunit = channum

part1(nunit) = nanmean(norm_chan(601:1175,nunit),1);
part2(nunit) = nanmean(norm_chan(1176:1750,nunit), 1);

  end
 parts = [part1,part2];
filename = ['C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\part1_part2_norm_power.mat'];
save(filename, 'parts');
%}
%%save the trials distributions
parts = struct();

  for nunit = channum
      if ~isempty(S_stim(nunit).S)
parts(nunit).part1 = S_stim(nunit).S(601:1175,1,:);
parts(nunit).part2 = S_stim(nunit).S(1176:1750,1,:);
      end
  end

filename = ['C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\part1_part2_power_trials.mat'];
save(filename, 'parts');
