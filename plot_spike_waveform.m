
%this script was developped to plot the spike waveform


gendatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [gendatadir 'refined_dataset']; 
gen_data_file = load(channelfilename);

%% information about the data
%{
tilt = nan(71,10);
xpos = nan(71,10);
ypos = nan(71,10);
for i = 1:71
    tilt(i,1:length(unique(gen_data_file.new_data(i).channel_data.tilt))) = unique(gen_data_file.new_data(i).channel_data.tilt);
    xpos(i,1:length(unique(gen_data_file.new_data(i).channel_data.xpos))) =unique(gen_data_file.new_data(i).channel_data.xpos);
    ypos(i,1:length(unique(gen_data_file.new_data(i).channel_data.ypos))) =unique(gen_data_file.new_data(i).channel_data.ypos);

end
%}

%%

spike_dat = gen_data_file.new_data(1).channel_data.wf.waveForms;

%pltting all the channels of all retained single units
idx = [1 3 5 7 9 11 13 15 17 19 21 23 2 4 6 8 10 12 14 16 18 20 22 24];
mean_data = nan(length(-40:40), length(1:24));
clear i ;
for chan = 1:length(gen_data_file.new_data)
   h = figure;
 for i = 1:24
 
   mean_data = squeeze(mean(spike_dat,2));

    sp = subplot(length(1:12), 2, idx(i) );
    plot(-40:40, mean_data(i,:))
    hold on
    plot([0 0], ylim,'k')

    if i == length(12)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}Spike Rate (spikes/s)'});
    end
   if i < 2 || (i >= 13)&&(i < 24)
        set(sp, 'XTick', [])
   end

      set(gca, 'linewidth',2)
      set(gca,'box','off')
     % h = subplot(1,1,1); 
     %set(h,'position',get(h,'position').*[1 1 1 1.2]);
 end
end

%plot mean spike waveform of the retained channel for every retained signle unit
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

mean_spikes = squeeze(nanmean(spikes,2));
 figure();
    for cn = 1:3
      h =subplot(1,3,cn);
     
   
    plot(-40:40, mean_spikes(:,cn),  'LineWidth',2)
    hold on
    plot(xlim, [0 0],'k')

      set(gca, 'linewidth',2)
       set(h,'position',get(h,'position').*[1 1 1.15 1])
       
      set(gca,'box','off')
      %{
      if cn >1 
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off';   
      end
      %}
      title(sprintf('%s',cellclass(cn)))
    end
   filename = strcat('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\spike_waveform\all_cells_means_yaxes');
   saveas(gcf, strcat(filename, '.svg')); 
   saveas(gcf, strcat(filename, '.png')); 
      
