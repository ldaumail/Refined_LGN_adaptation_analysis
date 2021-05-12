%This script was written to analyze multiunit activity

%Developped by Loic Daumail - 05/07/2021
%datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\fft_auto_units\';
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\triggered_auto_data\';

filenames = dir(strcat(datadir, '*cinterocdrft*'));

%create new directory to save channels:
channelsdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\selected_channels\';

%1)Load MUA, see if neural response in each channel is significant
for i = 1:length(filenames)
    filename = filenames(i).name;
    STIMdMUA = load(strcat(datadir, filename));
    %trials_data = data.STIM.sdftr;
    
    %As a first step, power of the MUA response at 4Hz is calculated for all trials of each channel, to
    %determine whether the response in each channel is significant
    %preallocate power and frequency, power at 4 hz
    power = nan(length(STIMdMUA.STIM.sdftr(1,:,1)),length(STIMdMUA.STIM.sdftr(1,1,:)),1025);
    freq = nan(length(STIMdMUA.STIM.sdftr(1,:,1)),length(STIMdMUA.STIM.sdftr(1,1,:)),1025);
    fourhzpower =nan(length(STIMdMUA.STIM.sdftr(1,:,1)),length(STIMdMUA.STIM.sdftr(1,1,:)));
    
    %experiment conditions
    data_header.contrast1 = STIMdMUA.STIM.contrast == 0 & STIMdMUA.STIM.fixedc >= 0.5; %trials indices with 0 contrast in DE, and contrast >0.5 in NDE
    data_header.contrast2 = STIMdMUA.STIM.contrast >= 0.5 & STIMdMUA.STIM.fixedc == 0; %trials indices with >0.5 contrast in DE, and contrast 0 in NDE
    data_header.contrast3 = STIMdMUA.STIM.contrast >=0.5 & STIMdMUA.STIM.fixedc >= 0.5; %trials indices with >0.5 contrast in DE, and contrast >0.5 in NDE
    data_header.contrast4 = STIMdMUA.STIM.contrast ==0 & STIMdMUA.STIM.fixedc == 0; %get logical indices of trials with 0 contrast in both eyes
    f = {'_DE0_NDE50','_DE50_NDE0','_DE50_NDE50'};
    
    %preallocate stats results
    signi = nan(length(STIMdMUA.STIM.sdftr(1,:,1)), 3);
    pvalue = nan(length(STIMdMUA.STIM.sdftr(1,:,1)), 3);
    ci = nan(length(STIMdMUA.STIM.sdftr(1,:,1)),3, 2);
    stats = struct();
    
    hypo = {'_DE0_NDE50_vsB','_DE50_NDE0_vsB','_DE50_NDE50_vsB'};
   
    for chan = 1:length(STIMdMUA.STIM.sdftr(1,:,1)) % for each channel
        for tr = 1:length(STIMdMUA.STIM.sdftr(1,1,:))  %for all trials of each channel, compute power at 4Hz
            %compute fft for every trial of a channel, for every channel
            [power(chan,tr,:), freq(chan,tr,:)] = calcFFT(squeeze(STIMdMUA.STIM.sdftr(600:1700,chan,tr)));
            
            %find the index of the frequency vector closest to 4hz and point to the
            %power value of this index for every trial, and store the value in
            %fourhzpower
            [val,index] = min(abs(4-freq(chan,tr,:)));
            fourhzpower(chan,tr) = power(chan,tr,index);
            
        end
        %perform t-test accross :
        %- DE50-NDE0 vs DE0-NDE0
        %- NDE50-DE0 vs DE0-NDE0
        %- NDE50-DE50 vs DE0-NDE0
        cont_power = struct();
        cont_power.cont1_power = fourhzpower(chan,data_header.contrast1(1:length(fourhzpower(chan,:))));
        cont_power.cont2_power = fourhzpower(chan,data_header.contrast2(1:length(fourhzpower(chan,:))));
        cont_power.cont3_power = fourhzpower(chan,data_header.contrast3(1:length(fourhzpower(chan,:))));
        cont_power.cont4_power = fourhzpower(chan,data_header.contrast4(1:length(fourhzpower(chan,:))));
        
        cont_power_cell = struct2cell(cont_power);
        for testnb = 1:3
            X = cont_power_cell{testnb};
            Y = cont_power.cont4_power;
            [signi(chan, testnb), pvalue(chan,testnb)] = ttest2(X,Y);
        end
        statistics = struct();
        statistics.significance = signi;
        statistics.pvalues = pvalue;
        
        %{
    [signi(chan, testnb), pvalue(chan,testnb), ci(chan,testnb,:), stats] = ttest2(X,Y);
  
    all_stats.strcat('stats', sprintf('channel %d',chan), hypo{testnb}) = stats;
    statistics = struct();
    statistics.significance = signi;
    statistics.pvalues = pvalue;
    statistics.ConfidenceInterval = ci;
    statistics.teststats = stats;
        %}
        %if test result is significant, save the channel dMUA data and power
        %data
        clear testnb
        for testnb = 1:3
            if (signi(chan,2) ~= 0 && signi(chan,1) == 0 || (signi(chan,3) ~= 0 && signi(chan,2) ~= 0))  
                
                channel_data = struct();
                channel_data.sdftr_chan = squeeze(STIMdMUA.STIM.sdftr(:,chan,:));
                %channel_data.hypo{2}.cont_dMUA_chan = STIMdMUA.STIM.sdftr(:,chan,data_header.contrast2(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
                %channel_data.hypo{3}.cont_dMUA_chan = STIMdMUA.STIM.sdftr(:,chan,data_header.contrast3(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
                channel_data.power_chan = power(chan, :,:);
                %channel_data.hypo{2}.cont_power_chan = power(data_header.contrast2(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
                %channel_data.hypo{3}.cont_power_chan = power(data_header.contrast3(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
                
                
                channel_data.hypo{1}.cont_stats_chan.significance = statistics.significance(chan,1);
                channel_data.hypo{2}.cont_stats_chan.significance = statistics.significance(chan,2);
                channel_data.hypo{3}.cont_stats_chan.significance = statistics.significance(chan,3);
                channel_data.hypo{1}.cont_stats_chan.pvalue = statistics.pvalues(chan,1);
                channel_data.hypo{2}.cont_stats_chan.pvalue = statistics.pvalues(chan,2);
                channel_data.hypo{3}.cont_stats_chan.pvalue = statistics.pvalues(chan,3);
                channel_data.onsets = STIMdMUA.STIM.onsets;
                channel_data.offsets = STIMdMUA.STIM.offsets;
                channel_data.trstart = STIMdMUA.STIM.trstart;
                channel_data.trend = STIMdMUA.STIM.trend;
                channel_data.presnum = STIMdMUA.STIM.presnum;
                channel_data.trial = STIMdMUA.STIM.trial;
                channel_data.tilt = STIMdMUA.STIM.tilt;
                channel_data.sf = STIMdMUA.STIM.sf;
                channel_data.contrast = STIMdMUA.STIM.contrast;
                channel_data.fixedc = STIMdMUA.STIM.fixedc;
                channel_data.diameter = STIMdMUA.STIM.diameter;
                channel_data.oridist = STIMdMUA.STIM.oridist;
                channel_data.phase = STIMdMUA.STIM.phase;
                channel_data.temporal_freq = STIMdMUA.STIM.temporal_freq;
                channel_data.xpos = STIMdMUA.STIM.xpos;
                channel_data.ypos = STIMdMUA.STIM.ypos;
                if exist('STIMdMUA.STIM.fix_x')
                    channel_data.fix_x = STIMdMUA.STIM.fix_x;
                end
                if exist('STIMdMUA.STIM.fix_y')
                    channel_data.fix_y = STIMdMUA.STIM.fix_y;
                end
                channel_data.photo_on = STIMdMUA.STIM.photo_on;
                channel_data.trg_photo = STIMdMUA.STIM.trg_photo;
                channel_data.refresh = STIMdMUA.STIM.refresh;
                channel_data.measured_refresh = STIMdMUA.STIM.measured_refresh;
                channel_data.paradigm = STIMdMUA.STIM.paradigm;
                %channel_data.sdftr_chan = STIMdMUA.STIM.sdftr(:,chan,:);
                channel_data.spk_bin_chan = STIMdMUA.STIM.spk_bin(:,chan,:);
                channel_data.label_chan = STIMdMUA.STIM.label(chan);
                
                %filename(strfind(filename, '.mat')) = [];
                filename = erase(filename, '.mat');
                channelfilename = [channelsdir strcat(filename, sprintf('_channel_%d_DE',chan))];
                save(strcat(channelfilename, '.mat'), 'channel_data');
   
            else
                
                %change contrast name and data if the significance is NDE50-
                %0DE vs Blank condition as it means the DE and NDE are inverted
                if signi(chan,2) == 0 && signi(chan,1) ~= 0
                    channel_data = struct();
                    channel_data.sdftr_chan = squeeze(STIMdMUA.STIM.sdftr(:,chan,:));
                    channel_data.power_chan = power(chan, :,:);
                    channel_data.hypo{1}.cont_stats_chan.significance = statistics.significance(chan,1);
                    channel_data.hypo{2}.cont_stats_chan.significance = statistics.significance(chan,2);
                    channel_data.hypo{3}.cont_stats_chan.significance = statistics.significance(chan,3);
                    channel_data.hypo{1}.cont_stats_chan.pvalue = statistics.pvalues(chan,1);
                    channel_data.hypo{2}.cont_stats_chan.pvalue = statistics.pvalues(chan,2);
                    channel_data.hypo{3}.cont_stats_chan.pvalue = statistics.pvalues(chan,3);
                    channel_data.onsets = STIMdMUA.STIM.onsets;
                    channel_data.offsets = STIMdMUA.STIM.offsets;
                    channel_data.trstart = STIMdMUA.STIM.trstart;
                    channel_data.trend = STIMdMUA.STIM.trend;
                    channel_data.presnum = STIMdMUA.STIM.presnum;
                    channel_data.trial = STIMdMUA.STIM.trial;
                    channel_data.tilt = STIMdMUA.STIM.tilt;
                    channel_data.sf = STIMdMUA.STIM.sf;
                    channel_data.contrast = STIMdMUA.STIM.fixedc; 
                    channel_data.fixedc = STIMdMUA.STIM.contrast;
                    channel_data.diameter = STIMdMUA.STIM.diameter;
                    channel_data.oridist = STIMdMUA.STIM.oridist;
                    channel_data.phase = STIMdMUA.STIM.phase;
                    channel_data.temporal_freq = STIMdMUA.STIM.temporal_freq;
                    channel_data.xpos = STIMdMUA.STIM.xpos;
                    channel_data.ypos = STIMdMUA.STIM.ypos;
                    if exist('STIMdMUA.STIM.fix_x')
                    channel_data.fix_x = STIMdMUA.STIM.fix_x;
                    end
                    if exist('STIMdMUA.STIM.fix_y')
                        channel_data.fix_y = STIMdMUA.STIM.fix_y;
                    end
                    channel_data.photo_on = STIMdMUA.STIM.photo_on;
                    channel_data.trg_photo = STIMdMUA.STIM.trg_photo;
                    channel_data.refresh = STIMdMUA.STIM.refresh;
                    channel_data.measured_refresh = STIMdMUA.STIM.measured_refresh;
                    channel_data.paradigm = STIMdMUA.STIM.paradigm;
                    %channel_data.sdftr_chan = STIMdMUA.STIM.sdftr(:,chan,:);
                    channel_data.spk_bin_chan = STIMdMUA.STIM.spk_bin(:,chan,:);
                    channel_data.label_chan = STIMdMUA.STIM.label(chan);
  
                    %channel_data.hypo{1}.cont_dMUA_chan = channel_data.hypo{2}.cont_dMUA_chan;
                    %channel_data.hypo{2}.cont_dMUA_chan = channel_data.hypo{1}.cont_dMUA_chan;
                    %channel_data.hypo{1}.cont_power_chan = channel_data.hypo{2}.cont_power_chan;
                    %channel_data.hypo{2}.cont_power_chan = channel_data.hypo{1}.cont_power_chan;
                    channel_data.hypo{1}.cont_stats_chan.significance = statistics.significance(chan,2);
                    channel_data.hypo{2}.cont_stats_chan.significance = statistics.significance(chan,1);
                    channel_data.hypo{1}.cont_stats_chan.pvalue = channel_data.hypo{2}.cont_stats_chan.pvalue;
                    channel_data.hypo{2}.cont_stats_chan.pvalue = channel_data.hypo{1}.cont_stats_chan.pvalue;
                   
                    
                    filename = erase(filename, '.mat');
                    %filename(strfind(filename, '.mat')) = [];
                    channelfilename = [channelsdir strcat(filename, sprintf('_channel_%d_NDE',chan))];
                    save(strcat(channelfilename, '.mat'), 'channel_data');
           
                end
                 
            end   
        end
    end
end
    %save all the power values of each session
    %powerfilename = [powerdir strcat('power',STIMBRdatafile)];
    %save(powerfilename, 'power');
   
%2)Detect peaks and select the corresponding trials and the peaks
%% Smooth data, and find peak locations, to further allow us to isolate peak values of unfiltered data in order to analyze them on R and fit a LMER
 
 clear all
 %Location of multiunits data isolated in the previous section
 channelsdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\selected_channels\';
 filenames = dir(strcat(channelsdir, '*_DE*')); %only include the DE files, as contrast range differs between eayes stimulated. If we use the CRFs to label channels, we will need a wide range of contrasts ==> "DE' files

 %find the multiple contrast levels present in the data
 %{
 allContLevels =[];
 for i =1:length(filenames)
     data = load(strcat(channelsdir, filenames(i).name));
     if ~isempty(filenames(i))
    %contLevels = unique(data.channel_data.fixedc);
    contLevels = unique(data.channel_data.contrast);
    allContLevels = unique([allContLevels; contLevels]);
     end
 end
 %}
 
[peakLocs, NoFiltMultiContMUA] = peakLocsTrialSelection(channelsdir, filenames);
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\adaptation_analysis\all_channels\all_locs_data_05112021';
save(strcat(allfilename, '.mat'), 'peakLocs');
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\adaptation_analysis\all_channels\NoFiltMultiContMUA_05112021';
save(strcat(allfilename, '.mat'), 'NoFiltMultiContMUA');
 
%3) Store Peaks and peak-triggered trials
[peak_vals, peak_aligned_trials] = peaksAndPeakTrigResps(peakLocs, NoFiltMultiContMUA);

%4) Label channels based on CRFs
%5) Assess adaptation per multiunit class