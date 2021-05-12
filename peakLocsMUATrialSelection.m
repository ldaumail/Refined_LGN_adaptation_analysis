function [peakLocs, NoFiltMultiContMUA] = peakLocsMUATrialSelection(channelsdir, filenames)
%This function allows to select the peak locations of the smoothed
%(low-pass filtered) data as well as the trials they correspond to.
%Last edited by Loic Daumail on 05/10/2021
%Location of multiunits data isolated in the previous section

%lets create contrast limits (bins to pool different contrast levels in different groups)
contLims = [0,0.1,0.3,0.5,0.7,1];
channum = 1: length(filenames);
xabs = -199:1300;
nyq = 500;

%mean_filtered_dSUA = struct();

FiltMultiContMUA =  struct();
NoFiltMultiContMUA = struct();
BsNoFiltMultiContMUA = struct();
%data_peaks = struct();
peakLocs = struct(); %store filtered data peak locations used to isolate peak values of unfiltered data

    
    

for i = channum
    filename = filenames(i).name;
    MUA_datafile = load(strcat(channelsdir, filenames(i).name));
    blankcontrast = MUA_datafile.channel_data.contrast ==  0 & MUA_datafile.channel_data.fixedc ==  0; %get logical indices of trials with 0 contrast in both eyes
    filename = erase(sprintf('x%s',char(filename)),'.mat');
    for n = 1:length(contLims)
        if n == 1
            contrastBin = MUA_datafile.channel_data.contrast >=  0.5 & MUA_datafile.channel_data.fixedc ==  0; %trials indices with 0 contrast in NDE, and contrast >0.5 in DE
        else
            if n>1
                contrastBin = (MUA_datafile.channel_data.fixedc >  contLims(n-1) & MUA_datafile.channel_data.fixedc <= contLims(n))& MUA_datafile.channel_data.contrast >=  0.5;
            end
        end
        trialidx = 1:length(MUA_datafile.channel_data.sdftr_chan(1,:)); %trial number of each trial for a given unit
        origin_data = nan(length(xabs)+401, length(trialidx));
        noFiltBs = nan(length(xabs), length(trialidx)); %to store the baseline corrected unfiltered data
        filtBs = nan(length(xabs), length(trialidx)); %to store the baseline corrected filtered data
        
        powerstim = nan(length(trialidx),1025);
        freqstim = nan(length(trialidx),1025);
        fourhzpowerstim =nan(length(trialidx),1);
        % bsl = nan(1, length(trialidx));
        mean_wnd1 = nan(1,length(trialidx));
        
        all_pks = nan(4,length(MUA_datafile.channel_data.sdftr_chan(1,contrastBin)));
        
        for tridx = trialidx
            
            all_data = MUA_datafile.channel_data.sdftr_chan(401:1900,tridx);
            origin_data(:,tridx) = MUA_datafile.channel_data.sdftr_chan(:,tridx);
            noFiltBs(:,tridx) = all_data(1:end)- mean(all_data(1:200));
            
            
            lpc       = 4.5; %low pass cutoff
            lWn       = lpc/nyq;
            [bwb,bwa] = butter(4,lWn,'low');
            lpdMUA      = filtfilt(bwb,bwa, noFiltBs(:,tridx));
            
            
            filtBs(:,tridx) = lpdMUA;
            %all_norm_lpdSUA(:,tridx) = (lpdSUA - min(lpdSUA))/(max(lpdSUA)- min(lpdSUA));
            mean_wnd1(tridx) = mean(lpdMUA(201:480)); %compute mean spiking response over 280ms following stimulus onset. 250+30
            
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
        filtered_dMUA_high = filtBs(:,contrastBin);
        %filtered_dMUA_blank = filtBs(:,blankcontrast);
        origin_data_high = origin_data(:,contrastBin);
        %origin_data_blank = origin_data(:,blankcontrast);
        bsl_origin_data_high = noFiltBs(:,contrastBin);
        
        mua_bsl =  mean(filtered_dMUA_high(1:200,:),1);
        
        
        for tr = 1:length(powerDE)
            %if mean_wnd1_DE(tr) > mean(sua_bsl)+1.96*std(sua_bsl)/sqrt(length(sua_bsl))  && powerDE(tr) > mean(power0)+1.96*std(power0)/sqrt(length(power0)) %/sqrt(length(sua_bsl)) /sqrt(length(power0))
            if mean_wnd1_DE(tr) > mean(mua_bsl)+1.96*std(mua_bsl)  && powerDE(tr) > mean(power0)+1.96*std(power0) %/sqrt(length(sua_bsl)) /sqrt(length(power0))
                
                filtered_dMUA_high(:,tr) = filtered_dMUA_high(:,tr);
                origin_data_high(:,tr) = origin_data_high(:,tr);
                bsl_origin_data_high(:,tr) = bsl_origin_data_high(:,tr);
            else
                
                filtered_dMUA_high(:,tr) = nan(length(filtered_dMUA_high(:,tr)),1);
                origin_data_high(:,tr) =  nan(length(origin_data_high(:,tr)),1);
                bsl_origin_data_high(:,tr) = nan(length(bsl_origin_data_high(:,tr)),1);
            end
        end
        
        %%%%%%%%%%%determine the first peak location for each trial of a given single unit %%%%%%%%
        all_locsdMUA_trials = nan(6,length(filtered_dMUA_high(1,:)));
        clear trial
        for trial = 1:length(filtered_dMUA_high(1,:))
            
            for ln = 1:550
                if filtered_dMUA_high(200+ln,trial) < filtered_dMUA_high(200+ln+1,trial) && ~all(isnan(filtered_dMUA_high(:,trial)))
                    [~,locsdMUA_trial] = findpeaks(filtered_dMUA_high(200+ln:1499,trial));
                    
                    %if peak1 is too small, peak2 becomes peak1
                    if filtered_dMUA_high(locsdMUA_trial(1)+200+ln,trial) >= 0.4*filtered_dMUA_high(locsdMUA_trial(2)+200+ln)
                        %store first peak location
                        all_locsdMUA_trials(1:length(locsdMUA_trial),trial) = locsdMUA_trial(1:end)+200+ln;
                    else
                        all_locsdMUA_trials(1:length(locsdMUA_trial(2:end)),trial) = locsdMUA_trial(2:end)+200+ln;
                        
                    end
                    
                    break
                end
            end
            
            if nnz(~isnan(all_locsdMUA_trials(:,trial))) >= 4 && ~all(isnan(all_locsdMUA_trials(:,trial)))
                %adjust location to the first data point of lpsu (+ln),
                
                all_pks(:,trial) = filtered_dMUA_high(all_locsdMUA_trials(1:4,trial), trial);
                filtered_dMUA_high(:,trial) = filtered_dMUA_high(:,trial);
                all_locsdMUA_trials(:,trial) = all_locsdMUA_trials(:,trial);
                origin_data_high(:,trial) = origin_data_high(:,trial);
                bsl_origin_data_high(:,trial) = bsl_origin_data_high(:,trial);
            else
                filtered_dMUA_high(:,trial) = nan(length(filtered_dMUA_high(:,trial)),1);
                all_locsdMUA_trials(:,trial) = nan(size(all_locsdMUA_trials(:,trial)));
                origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);
                bsl_origin_data_high(:,trial) =  nan(length(bsl_origin_data_high(:,trial)),1);
                
            end
            
            if ~all(isnan(all_locsdMUA_trials(:,trial))) && (all_locsdMUA_trials(4,trial) ~= 1500) %remove trials for which the pk4 is the last data point (= not a peak if this happens)
                %adjust location to the first data point of lpsu (+ln),
                
                all_pks(:,trial) = filtered_dMUA_high(all_locsdMUA_trials(1:4,trial), trial);
                filtered_dMUA_high(:,trial) = filtered_dMUA_high(:,trial);
                all_locsdMUA_trials(:,trial) = all_locsdMUA_trials(:,trial);
                origin_data_high(:,trial) = origin_data_high(:,trial);
                bsl_origin_data_high(:,trial) = bsl_origin_data_high(:,trial);
            else
                all_pks(:,trial) = nan(length(all_pks(:,trial)),1);
                filtered_dMUA_high(:,trial) = nan(length(filtered_dMUA_high(:,trial)),1);
                all_locsdMUA_trials(:,trial) = nan(size(all_locsdMUA_trials(:,trial)));
                origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);
                bsl_origin_data_high(:,trial) =  nan(length(bsl_origin_data_high(:,trial)),1);
                
                
            end
        end
        %{
        figure(); plot(-199:1300, filtered_dMUA_high(1:1500,:))
        hold on
        plot(all_locsdMUA_trials(1:4,:)-200, all_pks(:,:))
        set(gca,'box','off')
        %}
        %%% reject outlier peaks and the corresponding trials in
        %%% filtered_dSUA_high
        
        
        %%%%%%%%%%reject if there is a peak 1 outlier, if the max peak value in the baseline is an outlier %%%%%%%%%
        
        % First find peaks before stimulus onset
        
        bsl_peaks = nan(1, length(filtered_dMUA_high(1,:)));
        clear tr
        for tr = 1:length(filtered_dMUA_high(1,:))
            
            for loc = 1:200
                if filtered_dMUA_high(loc,tr) < filtered_dMUA_high(loc+1,tr) && ~all(isnan(filtered_dMUA_high(:,tr)))
                    if length(filtered_dMUA_high(loc:200,tr)) >= 3
                        if ~isempty(findpeaks(filtered_dMUA_high(loc:200,tr)))
                            [~, bsl_peak_locs] = findpeaks(filtered_dMUA_high(loc:200,tr));
                            bsl_peaks(1,tr) = max(filtered_dMUA_high(bsl_peak_locs+loc,tr));
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
        for tr = 1:length(origin_data_high(1,:))
            %exclude trials
            if p1outliers(tr) == 0 && ~all(isnan(all_pks(:,tr))) && out_bsl_peaks(tr) ==0
                
                origin_data_high(:,tr) = origin_data_high(:, tr);
                all_locsdMUA_trials(:,tr) = all_locsdMUA_trials(:,tr);
                %filtered_dMUA_high(:,tr) = filtered_dMUA_high(:, tr);
                %all_pks(:, tr) = all_pks(:,tr);
                %bsl_origin_data_high(:,tr) =  bsl_origin_data_high(:,tr);
                
            else
                origin_data_high(:,tr) = nan(length(origin_data_high(:,tr)),1);
                all_locsdMUA_trials(:,tr) = nan(size(all_locsdMUA_trials(:,tr)));
                %filtered_dMUA_high(:,tr) = nan(length(filtered_dMUA_high(:,tr)),1);
                %all_pks(:,tr) = nan(length(all_pks(:,tr)),1);
                %bsl_origin_data_high(:,tr) =  nan(length(bsl_origin_data_high(:,tr)),1);
                
            end
        end
        origin_data_high = origin_data_high(:,~all(isnan(origin_data_high)));
        all_locsdMUA_trials =  all_locsdMUA_trials(:,~all(isnan(all_locsdMUA_trials)));
        
        %filtered_dMUA_high = filtered_dMUA_high(:,~all(isnan(filtered_dMUA_high))); % for nan - cols
        %all_pks = all_pks(:, ~all(isnan(all_pks)));
        %bsl_origin_data_high = bsl_origin_data_high(:,~all(isnan(bsl_origin_data_high)));
        
        
        binNb = sprintf('bin%d', n);
        if n ==1 && length(origin_data_high(1,:)) >=10 %first bin == high contrast monocular condition will serve as an indicator of the minimum number of trials required for the analysis
            
            NoFiltMultiContMUA.(filename).(binNb) = origin_data_high;
            peakLocs.(filename).(binNb) = all_locsdMUA_trials; %create dynamical peak locations structures
            %FiltMultiContMUA.(filename).(binNb) =  filtered_dMUA_high;
            %BsNoFiltMultiContMUA.(filename).(binNb) = bsl_origin_data_high;
            
        elseif n == 1 && length(origin_data_high(1,:)) <10
            NoFiltMultiContMUA.(filename).(binNb) = [];
            peakLocs.(filename).(binNb) = [];
            %all_pks(:,:) = [];
            % FiltMultiContMUA.(filename).(binNb) =  [];
            
            
        elseif n > 1 && length(origin_data_high(1,:)) >=10
            NoFiltMultiContMUA.(filename).(binNb) = origin_data_high;
            peakLocs.(filename).(binNb) = all_locsdMUA_trials; %create dynamical peak locations structures
            
            %FiltMultiContMUA.(filename).(binNb) =  filtered_dMUA_high;
            % BsNoFiltMultiContMUA.(filename).(binNb) = bsl_origin_data_high;
        elseif n > 1 && length(origin_data_high(1,:)) <10
            NoFiltMultiContMUA.(filename).(binNb) = [];
            peakLocs.(filename).(binNb) = [];
        end
        
        
        %data_peaks(i).namelist = all_pks(:,~all(isnan(all_pks)));
        %all_pks = all_pks(:,~all(isnan(all_pks)));
        %channelfilename = [unitsDir 'su_peaks_03032020_corrected\individual_units\' filename 'multiContrast'];
        %save(strcat(channelfilename, '.mat'), 'peakLocs');
        
    end
end