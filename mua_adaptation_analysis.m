%This script was written to analyze multiunit activity

%Developped by Loic Daumail - 05/07/2021
%datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\fft_auto_units\';
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\triggered_auto_data\';

filenames = dir(strcat(datadir, '*cinterocdrft*'));

%create new directory to save channels:
channelsdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\selected_channels\';
%contLims = [0,0.1,0.3,0.5,0.7,1];
%1)Load MUA
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
                channel_data.sdftr_chan = STIMdMUA.STIM.sdftr(:,chan,:);
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
                channelfilename = [channelsdir strcat(filename, sprintf('channel_%d_DE',chan))];
                save(strcat(channelfilename, '.mat'), 'channel_data');
   
            else
                
                %change contrast name and data if the significance is NDE50-
                %0DE vs Blank condition as it means the DE and NDE are inverted
                if signi(chan,2) == 0 && signi(chan,1) ~= 0
                    channel_data = struct();
                    channel_data.sdftr_chan = STIMdMUA.STIM.sdftr(:,chan,:);
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
                    channelfilename = [channelsdir strcat(filename, sprintf('channel_%d_NDE',chan))];
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
 %find the multiple contrast levels present in the data
 
 
 allContLevels =0;
 for i =1:71
     if ~isempty(filenames(i))
    contLevels = unique(unitsData.new_data(i).channel_data.fixedc);
    allContLevels = unique([allContLevels; contLevels]);
     end
 end
 %lets create contrast limits (bins to pool different contrast levels)
contLims = [0,0.1,0.3,0.5,0.7,1];  
channum = 1: length(unitsData.new_data);
xabs = -199:1300;
nyq = 500;

%mean_filtered_dSUA = struct();


 FiltMultiContSUA =  struct();
 NoFiltMultiContSUA = struct();
 BsNoFiltMultiContSUA = struct();
%data_peaks = struct();
peakLocs = struct(); %store filtered data peak locations used to isolate peak values of unfiltered data

for n = 1:length(contLims)
    clear i
     for i = channum
        if ~isempty(filenames{i})
           filename = filenames(i);
           blankcontrast = unitsData.new_data(i).channel_data.contrast ==  0 & unitsData.new_data(i).channel_data.fixedc ==  0; %get logicacal indices of trials with 0 contrast in both eyes
           if n == 1
               contrastBin =unitsData.new_data(i).channel_data.contrast >=  0.5 & unitsData.new_data(i).channel_data.fixedc ==  0; %trials indices with 0 contrast in NDE, and contrast >0.5 in DE
               else
               if n>1
                   contrastBin = (unitsData.new_data(i).channel_data.fixedc >  contLims(n-1) & unitsData.new_data(i).channel_data.fixedc <= contLims(n))& unitsData.new_data(i).channel_data.contrast >=  0.5; 
               end
           end
        trialidx = 1:length(unitsData.new_data(i).channel_data.sdftr_chan(1,:)); %trial number of each trial for a given unit
        noFiltBs = nan(length(xabs), length(trialidx)); %to store the baseline corrected unfiltered data
        filtBs = nan(length(xabs), length(trialidx));
        origin_data = nan(length(xabs)+401, length(trialidx));
        %all_norm_lpdSUA= nan(length(xabs),length(trialidx));

        powerstim = nan(length(trialidx),1025);
        freqstim = nan(length(trialidx),1025);
        fourhzpowerstim =nan(length(trialidx),1);
       % bsl = nan(1, length(trialidx));
        mean_wnd1 = nan(1,length(trialidx));

        all_pks = nan(4,length(unitsData.new_data(i).channel_data.sdftr_chan(1,contrastBin)));

         for tridx = trialidx

                all_data = unitsData.new_data(i).channel_data.sdftr_chan(401:1900,tridx);
                origin_data(:,tridx) = unitsData.new_data(i).channel_data.sdftr_chan(:,tridx);
                noFiltBs(:,tridx) = all_data(1:end)- mean(all_data(1:200));
                
               
                lpc       = 4.5; %low pass cutoff
                lWn       = lpc/nyq;
                [bwb,bwa] = butter(4,lWn,'low');
                lpdSUA      = filtfilt(bwb,bwa, noFiltBs(:,tridx));


                filtBs(:,tridx) = lpdSUA;
                %all_norm_lpdSUA(:,tridx) = (lpdSUA - min(lpdSUA))/(max(lpdSUA)- min(lpdSUA));
                mean_wnd1(tridx) = mean(lpdSUA(201:480)); %compute mean spiking response over 280ms following stimulus onset. 250+30

                %%% power


                [powerstim(tridx,:), freqstim(tridx,:)] = calcFFT(all_data(200:1350)); %fourrier transform

                %find the index of the frequency vector closest to 4hz and point to the
                %power value of this index for every trial, and store the value in
                %fourhzpower
                [val,index] = min(abs(4-freqstim(tridx,:))); %find index closest to 4Hz
                fourhzpowerstim(tridx,1) = powerstim(tridx,index); %get power at that index, assumed to be 4Hz

         end

      %%%%%%%%%%% %reject trials below Mean + 1.96*STD in the blank condition %%%%%%
       %power related variables
       power0 = fourhzpowerstim(blankcontrast); %power of responses in blank condition
       powerDE = fourhzpowerstim(contrastBin); %power of responses with contrast stimulus >0 in DE and 0 contrast in NDE

       %spiking activity related variables
       mean_wnd1_DE =mean_wnd1(contrastBin);
       filtered_dSUA_high = filtBs(:,contrastBin);
       filtered_dSUA_blank = filtBs(:,blankcontrast);
       origin_data_high = origin_data(:,contrastBin);
       origin_data_blank = origin_data(:,blankcontrast);
       bsl_origin_data_high = noFiltBs(:,contrastBin);
       
       sua_bsl =  mean(filtered_dSUA_high(1:200,:),1);
       

       for tr = 1:length(powerDE)
          %if mean_wnd1_DE(tr) > mean(sua_bsl)+1.96*std(sua_bsl)/sqrt(length(sua_bsl))  && powerDE(tr) > mean(power0)+1.96*std(power0)/sqrt(length(power0)) %/sqrt(length(sua_bsl)) /sqrt(length(power0))
           if mean_wnd1_DE(tr) > mean(sua_bsl)+1.96*std(sua_bsl)  && powerDE(tr) > mean(power0)+1.96*std(power0) %/sqrt(length(sua_bsl)) /sqrt(length(power0))

              filtered_dSUA_high(:,tr) = filtered_dSUA_high(:,tr);
              origin_data_high(:,tr) = origin_data_high(:,tr);
              bsl_origin_data_high(:,tr) = bsl_origin_data_high(:,tr);
           else

               filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
               origin_data_high(:,tr) =  nan(length(origin_data_high(:,tr)),1);
               bsl_origin_data_high(:,tr) = nan(length(bsl_origin_data_high(:,tr)),1);
           end
       end

        %%%%%%%%%%%determine the first peak location for each trial of a given single unit %%%%%%%%
        all_locsdSUA_trials = nan(6,length(filtered_dSUA_high(1,:)));
        clear trial
        for trial = 1:length(filtered_dSUA_high(1,:))
            
            for ln = 1:550
                if filtered_dSUA_high(200+ln,trial) < filtered_dSUA_high(200+ln+1,trial) && ~all(isnan(filtered_dSUA_high(:,trial)))
                    [~,locsdSUA_trial] = findpeaks(filtered_dSUA_high(200+ln:1499,trial));
                     
                    %if peak1 is too small, peak2 becomes peak1
                    if filtered_dSUA_high(locsdSUA_trial(1)+200+ln,trial) >= 0.4*filtered_dSUA_high(locsdSUA_trial(2)+200+ln)
                        %store first peak location
                        all_locsdSUA_trials(1:length(locsdSUA_trial),trial) = locsdSUA_trial(1:end)+200+ln;
                    else
                        all_locsdSUA_trials(1:length(locsdSUA_trial(2:end)),trial) = locsdSUA_trial(2:end)+200+ln;
                        
                    end
                    
                    break
                end
            end
            
            if nnz(~isnan(all_locsdSUA_trials(:,trial))) >= 4 && ~all(isnan(all_locsdSUA_trials(:,trial)))
                %adjust location to the first data point of lpsu (+ln),
                
                all_pks(:,trial) = filtered_dSUA_high(all_locsdSUA_trials(1:4,trial), trial);
                filtered_dSUA_high(:,trial) = filtered_dSUA_high(:,trial);
                all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
                origin_data_high(:,trial) = origin_data_high(:,trial);
                bsl_origin_data_high(:,trial) = bsl_origin_data_high(:,trial);
            else
                filtered_dSUA_high(:,trial) = nan(length(filtered_dSUA_high(:,trial)),1);
                all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
                origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);
                bsl_origin_data_high(:,trial) =  nan(length(bsl_origin_data_high(:,trial)),1);
                
            end
            
            if ~all(isnan(all_locsdSUA_trials(:,trial))) && (all_locsdSUA_trials(4,trial) ~= 1500) %remove trials for which the pk4 is the last data point (= not a peak if this happens)
                %adjust location to the first data point of lpsu (+ln),
                
                all_pks(:,trial) = filtered_dSUA_high(all_locsdSUA_trials(1:4,trial), trial);
                filtered_dSUA_high(:,trial) = filtered_dSUA_high(:,trial);
                all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
                origin_data_high(:,trial) = origin_data_high(:,trial);
                bsl_origin_data_high(:,trial) = bsl_origin_data_high(:,trial);
            else
                all_pks(:,trial) = nan(length(all_pks(:,trial)),1);
                filtered_dSUA_high(:,trial) = nan(length(filtered_dSUA_high(:,trial)),1);
                all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
                origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);
                bsl_origin_data_high(:,trial) =  nan(length(bsl_origin_data_high(:,trial)),1);
                
                
            end
        end
        %{
        figure(); plot(-199:1300, filtered_dSUA_high(1:1500,:))
        hold on
        plot(all_locsdSUA_trials(1:4,1)-200, all_pks(:,1))
        set(gca,'box','off')
        %}
        %%% reject outlier peaks and the corresponding trials in
        %%% filtered_dSUA_high
        
        
        %%%%%%%%%%reject if there is a peak 1 outlier, if the max peak value in the baseline is an outlier %%%%%%%%%
        
        % First find peaks before stimulus onset

        bsl_peaks = nan(1, length(filtered_dSUA_high(1,:)));
        clear tr
        for tr = 1:length(filtered_dSUA_high(1,:))
            
            for loc = 1:200
                if filtered_dSUA_high(loc,tr) < filtered_dSUA_high(loc+1,tr) && ~all(isnan(filtered_dSUA_high(:,tr)))
                    if length(filtered_dSUA_high(loc:200,tr)) >= 3
                        if ~isempty(findpeaks(filtered_dSUA_high(loc:200,tr)))
                            [~, bsl_peak_locs] = findpeaks(filtered_dSUA_high(loc:200,tr));
                            bsl_peaks(1,tr) = max(filtered_dSUA_high(bsl_peak_locs+loc,tr));
                        else
                            bsl_peaks(1,tr) = NaN;
                        end
                    end
                    break
                end
            end
        end

        out_bsl_peaks = isoutlier(bsl_peaks);

        p1outliers = isoutlier(all_pks(1,:));
        clear tr
        for tr = 1:length(filtered_dSUA_high(1,:))
            %exclude trials
            if p1outliers(tr) == 0 && ~all(isnan(all_pks(:,tr))) && out_bsl_peaks(tr) ==0 

                filtered_dSUA_high(:,tr) = filtered_dSUA_high(:, tr);
                all_pks(:, tr) = all_pks(:,tr);
                all_locsdSUA_trials(:,tr) = all_locsdSUA_trials(:,tr);
                origin_data_high(:,tr) = origin_data_high(:, tr);
                bsl_origin_data_high(:,tr) =  bsl_origin_data_high(:,tr);


            else 
                filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
                all_pks(:,tr) = nan(length(all_pks(:,tr)),1);
                all_locsdSUA_trials(:,tr) = nan(size(all_locsdSUA_trials(:,tr)));
                origin_data_high(:,tr) = nan(length(origin_data_high(:,tr)),1);
                bsl_origin_data_high(:,tr) =  nan(length(bsl_origin_data_high(:,tr)),1);

            end
        end
       filtered_dSUA_high = filtered_dSUA_high(:,~all(isnan(filtered_dSUA_high))); % for nan - cols
       all_locsdSUA_trials =  all_locsdSUA_trials(:,~all(isnan(all_locsdSUA_trials)));
       all_pks = all_pks(:, ~all(isnan(all_pks)));
       origin_data_high = origin_data_high(:,~all(isnan(origin_data_high)));
       bsl_origin_data_high = bsl_origin_data_high(:,~all(isnan(bsl_origin_data_high)));

       filename = sprintf('x%s',char(filename));
       binNb = sprintf('bin%d', n);
       if n ==1 && length(filtered_dSUA_high(1,:)) >=10 %first bin == high contrast monocular condition will serve as an indicator of the minimum number of trials required for the analysis
           
           %eval(['peakLocs.' num2str(i) '.bin' num2str(n) ' = all_locsdSUA_trials;'])
          
           peakLocs.(filename).(binNb) = all_locsdSUA_trials; %create dynamical peak locations structures
           FiltMultiContSUA.(filename).(binNb) =  filtered_dSUA_high;
           % FiltMultiContSUA.(filename).bin0 =  filtered_dSUA_blank;
           NoFiltMultiContSUA.(filename).(binNb) = origin_data_high;
           BsNoFiltMultiContSUA.(filename).(binNb) = bsl_origin_data_high;
           %NoFiltMultiContSUA.(filename).bin0 = origin_data_blank;
           FiltMultiContSUA.(filename).cellclass = cellClass{i}; %get the cell class of selected units
           NoFiltMultiContSUA.(filename).cellclass = cellClass{i};
           BsNoFiltMultiContSUA.(filename).cellclass = cellClass{i};
       elseif n == 1 && length(filtered_dSUA_high(1,:)) <10
           
           all_pks(:,:) = [];
           FiltMultiContSUA.(filename).(binNb) =  [];
           NoFiltMultiContSUA.(filename).(binNb) = [];
           peakLocs.(filename).(binNb) = [];
        
       elseif n > 1
           peakLocs.(filename).(binNb) = all_locsdSUA_trials; %create dynamical peak locations structures
           FiltMultiContSUA.(filename).(binNb) =  filtered_dSUA_high;
           % FiltMultiContSUA.(filename).bin0 =  filtered_dSUA_blank;
           NoFiltMultiContSUA.(filename).(binNb) = origin_data_high;
           BsNoFiltMultiContSUA.(filename).(binNb) = bsl_origin_data_high;
           %NoFiltMultiContSUA.(filename).bin0 = origin_data_blank;
           FiltMultiContSUA.(filename).cellclass = cellClass{i}; %get the cell class of selected units
           NoFiltMultiContSUA.(filename).cellclass = cellClass{i};
           BsNoFiltMultiContSUA.(filename).cellclass = cellClass{i};
       end
       

     %data_peaks(i).namelist = all_pks(:,~all(isnan(all_pks)));
     %all_pks = all_pks(:,~all(isnan(all_pks)));
    %channelfilename = [unitsDir 'su_peaks_03032020_corrected\individual_units\' filename 'multiContrast'];
    %save(strcat(channelfilename, '.mat'), 'peakLocs');
        end
     end
end

%3) Store Peaks and trials
%4) Label channels based on CRFs
%5) Assess adaptation per multiunit class