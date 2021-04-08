%% Figure 1: Narrow vs broad spike waveforms

gendatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [gendatadir 'refined_dataset']; 
gen_data_file = load(channelfilename);

%plot spike waveform of the retained channels for every retained signle unit
%use the peakvals to see which unit was retained in the data cleaning
%process
channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);



%spikes = nan(81,25,3);

spikes = nan(81,length(peakvals.peak_vals));

  
%for nc = 1:length(cellclass)
% clear layer_idx
%layer_idx = find(strcmp(layer, cellclass(nc)));

    for un = 1:length(peakvals.peak_vals)
         if ~isempty(peakvals.peak_vals(un).peak)
         spike_dat = gen_data_file.new_data(un).channel_data.wf.waveForms;
%mean of just one unit

        mean_data = squeeze(mean(spike_dat,2));
        chan_idx = str2double(gen_data_file.new_data(un).channel_data.chan(3:4)); 
        spikes(:,un) = mean_data(chan_idx,:);
         end
    end 
    
%end

%spikes = reshape(spikes,81,75);

%index = reshape(1:75, 3, 25).';

[~,maxLocation] = max(spikes,[],1);
[~,minLocation] = min(spikes,[],1);

wvDuration = abs(maxLocation - minLocation);



%% Waveform Duration distribution

figure();
histogram(wvDuration(find(wvDuration)).*1000/30,100)
xlabel('Waveform duration (microsecond)')


%% PCA
cenSpikes = spikes - nanmean(spikes,1); %center the data on the mean

[coeff, score, latent] = pca(cenSpikes'); 

%% plot PCA results
x = score(:,1);
y = score(:,3);
z = score(:,2);
figure();
plot3(x,y,z,'.')
hold on
yline(0)
hold on
xline(0)
xlabel('PCA1')
ylabel('PCA3')
zlabel('PCA2')
grid on

%% same but separating cell classes
layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];
cellclass = [ 'M', 'P', 'K'];
spikes = nan(81,25,3);
wvDuration = nan(25,3);
  
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
[~,maxLocation] = max(spikes(:,:,nc),[],1);
[~,minLocation] = min(spikes(:,:,nc),[],1);

wvDuration(:,nc) = maxLocation - minLocation;


end



%% Plot waveform durations per cell class
 figure();
   for nc = 1:3
wvD =wvDuration(:,nc);
     h =subplot(length(1),3,nc);
    histogram(wvD(find(wvD)).*1000/30,100)
   % hold on
   % plot(xlim, [0 0],'k')

      
       set(h,'position',get(h,'position').*[1 1 1.15 1])
       
      set(gca,'box','off')
      
      title(sprintf('Waveform duration distribution of %s cells', cellclass(nc))) 
    if nc ==1
        xlabel('Waveform Duration (microsecond)')
    end
    end
    %end
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\spike_waveform\all_units_spike_waveforms_stacked');
   %saveas(gcf, strcat(filename, '.svg')); 
   %saveas(gcf, strcat(filename, '.png')); 
