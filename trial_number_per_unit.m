%This script was developped to plot a distribution of the number of trials
%analyzed for each unit
%Loic Daumail 05/27/2022




%Using the updated functions to isolate peaks
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;


unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);


cellClass = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
cellClass([1,46,55]) = [];
%{
 allContLevels =0;
 alltilts =[];
 for i =1:71
     if ~isempty(filenames(i))
    contLevels = unique(unitsData.new_data(i).channel_data.fixedc);
    allContLevels = unique([allContLevels; contLevels]);
    tilts = find(unique(unitsData.new_data(i).channel_data.tilt));
    alltilts = [alltilts; tilts];
     end
 end
%}


%select trials with at least 4 peak values, of convenient quality. Keep peak locations and trial responses: 
[peakLocs, NoFiltMultiContSUA] = peakLocsTrialSelection(unitsData, filenames);

%Store Peaks and peak-triggered trials
[peak_vals, peak_aligned_trials] = peaksAndPeakTrigResps(peakLocs, NoFiltMultiContSUA);



%% Plot number of trials per unit
%get number of trials per unit in mono and binocular condition
bins = [1,6];
tn = [];
filenames = fieldnames(peak_aligned_trials);
for i = 1:length(filenames)
   if length(fieldnames(peak_aligned_trials.(filenames{i}).origin)) == 2
    for b =1:2
    bin = sprintf('bin%d',bins(b));
    tn(i,b) = nnz(~isnan(peak_aligned_trials.(filenames{i}).origin.(bin).pk1(1,:)));
    end
   end
end
tn = tn(tn(:,1) ~= 0,:);

%set up colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);
%set up data for gramm
trial_numbers = [tn(:,1);tn(:,2)];
condition = [repmat({'Monocular'},length(trial_numbers)/2,1); repmat({'Binocular'},length(trial_numbers)/2,1)];

%Plot
%1) Plot the scatter
clear g
g(1,1)=gramm('x',condition,'y',trial_numbers,'color',condition);

g(1,1).geom_jitter('width',0.4,'height',0); 
g(1,1).set_names('x','','y', '');
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]);
g(1,1).no_legend();

%Plot the mean and std for trial numbers
alltnMono = trial_numbers(strcmp(condition, 'Monocular'));
alltnBino = trial_numbers(strcmp(condition, 'Binocular'));
meanTn = nanmean([alltnMono,alltnBino],1);
shortCondition = unique(condition);
stdTns = std([alltnMono,alltnBino], 'omitnan');%/sqrt(length([alltnMono,alltnBino]));

ci_low = meanTn - stdTns;
ci_high = meanTn + stdTns;
g(1,1).update('x', shortCondition, 'y',meanTn,...
    'ymin',ci_low,'ymax',ci_high,'color',shortCondition);
%g(1,1).set_color_options('map','r');
g(1,1).geom_point('dodge',0.5);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_names('x','','y', 'Number of Trials');
g(1,1).no_legend();
g(1,1).set_point_options('base_size',7);

g(1,2) = gramm('x',ones(1,length(tn(:,1))), 'y', tn(:,1)-tn(:,2) );
g(1,2).geom_jitter('width',0.3,'height',0); 
g(1,2).axe_property('xlim', [0 2], 'ylim',[-30 30], 'XTickLabel','','XTick',''); 
g(1,2).set_names('x','','y', 'Difference Monocular-Binocular (Number of Trials)');
g.set_title('Number of trials in selected units post-processing');
g(1,2).set_color_options('map',[0.5 0.5 0.5]);


%Plot mean difference and 95%CI
meanDiff = nanmean([alltnMono-alltnBino],1);
stdDiffs = std([alltnMono-alltnBino], 'omitnan');%/sqrt(length([alltnMono-alltnBino]));
ci_low = meanDiff - stdDiffs;
ci_high = meanDiff + stdDiffs;
g(1,2).update('x', ones(1,1), 'y',meanDiff,...
    'ymin',ci_low,'ymax',ci_high);
g(1,2).geom_point('dodge',0.5);
g(1,2).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,2).axe_property('xlim', [0 2], 'ylim',[-30 30], 'XTickLabel','','XTick',''); 
%g(1,2).set_names('x','','y', '');
g(1,2).no_legend();
g(1,2).set_point_options('base_size',7);

figure('Position',[100 100 800 550]);

g.draw();

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_boxpointbar_std_numtrials_per_unit');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% Use previously isolated peaks

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']; 
peak_aligned_trials = load(channelfilename);


filenames = fieldnames(peak_aligned_trials.peak_aligned_trials);
bins = [1,6];


tn = [];
filenames = fieldnames(peak_aligned_trials.peak_aligned_trials);
for i = 1:length(filenames)
   if length(fieldnames(peak_aligned_trials.peak_aligned_trials.(filenames{i}).origin)) == 2
    for b =1:2
    bin = sprintf('bin%d',bins(b));
    tn(i,b) = nnz(~isnan(peak_aligned_trials.peak_aligned_trials.(filenames{i}).origin.(bin).pk1(1,:)));
    end
   end
end
tn = tn(tn(:,1) ~= 0,:);

%set up colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);
%set up data for gramm
trial_numbers = [tn(:,1);tn(:,2)];
condition = [repmat({'Monocular'},length(trial_numbers)/2,1); repmat({'Binocular'},length(trial_numbers)/2,1)];

%Plot
%1) Plot the scatter
clear g
g(1,1)=gramm('x',condition,'y',trial_numbers,'color',condition);

g(1,1).geom_jitter('width',0.4,'height',0); 
g(1,1).set_names('x','','y', '');
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]);
g(1,1).no_legend();

%Plot the mean and std for trial numbers
alltnMono = trial_numbers(strcmp(condition, 'Monocular'));
alltnBino = trial_numbers(strcmp(condition, 'Binocular'));
meanTn = nanmean([alltnMono,alltnBino],1);
shortCondition = unique(condition);
stdTns = std([alltnMono,alltnBino], 'omitnan');%/sqrt(length([alltnMono,alltnBino]));

ci_low = meanTn' - stdTns';
ci_high = meanTn' + stdTns';
g(1,1).update('x', shortCondition, 'y',meanTn,...
    'ymin',ci_low,'ymax',ci_high,'color',shortCondition);
%g(1,1).set_color_options('map','r');
g(1,1).geom_point('dodge',0.5);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_names('x','','y', 'Number of Trials');
g(1,1).no_legend();
g(1,1).set_point_options('base_size',7);

g(1,2) = gramm('x',ones(1,length(tn(:,1))), 'y', tn(:,1)-tn(:,2) );
g(1,2).geom_jitter('width',0.3,'height',0); 
g(1,2).axe_property('xlim', [0 2], 'ylim',[-30 30], 'XTickLabel','','XTick',''); 
g(1,2).set_names('x','','y', 'Difference Monocular-Binocular (Number of Trials)');
g.set_title('Number of trials in selected units post-processing');
g(1,2).set_color_options('map',[0.5 0.5 0.5]);


%Plot mean difference and +-std
meanDiff = nanmean([alltnMono-alltnBino],1);
stdDiffs = std([alltnMono-alltnBino], 'omitnan');%/sqrt(length([alltnMono-alltnBino]));
ci_low = meanDiff - stdDiffs;
ci_high = meanDiff + stdDiffs;
g(1,2).update('x', ones(1,1), 'y',meanDiff,...
    'ymin',ci_low,'ymax',ci_high);
g(1,2).geom_point('dodge',0.5);
g(1,2).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,2).axe_property('xlim', [0 2], 'ylim',[-20 30], 'XTickLabel','','XTick',''); 
%g(1,2).set_names('x','','y', '');
g(1,2).no_legend();
g(1,2).set_point_options('base_size',7);

figure('Position',[100 100 800 550]);

g.draw();

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_boxpointbar_std_numtrials_per_unit_old');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));


%% Assess trial number before and after preprocessing

%Load data
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;
unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);

contLims = [0,0.1,0.3,0.5,0.7,1];
channum = 1: length(filenames);
xabs = -199:1300;
nyq = 500;

%mean_filtered_dSUA = struct();

%FiltMultiContMUA =  struct();
NoFiltMultiContSUA = struct();
%BsNoFiltMultiContMUA = struct();
%data_peaks = struct();
peakLocs = struct(); %store filtered data peak locations used to isolate peak values of unfiltered data
binSpk = struct();
trialnum = struct(); %number of trials in each condition before preprocessing    
    

for i = channum
    filename = filenames{i};
    if ~isempty(filenames{i})
        SUA_data = unitsData.new_data(i);
        blankcontrast = SUA_data.channel_data.contrast ==  0 & SUA_data.channel_data.fixedc ==  0; %get logical indices of trials with 0 contrast in both eyes
        filename = erase(sprintf('x%s',char(filename)),'.mat');
        for n = 1:length(contLims)
            if n == 1
                contrastBin = SUA_data.channel_data.contrast >=  0.5 & SUA_data.channel_data.fixedc ==  0; %trials indices with 0 contrast in NDE, and contrast >0.5 in DE       
            else
                if n>1
                    contrastBin = (SUA_data.channel_data.fixedc >  contLims(n-1) & SUA_data.channel_data.fixedc <= contLims(n))& SUA_data.channel_data.contrast >=  0.5;   
                end
            end
            binNb = sprintf('bin%d', n);
            trialnum.before.(filename).(binNb) = nnz(contrastBin);  %trial num before preprocessing
            
            trialidx = 1:length(SUA_data.channel_data.sdftr_chan(1,:)); %trial number of each trial for a given unit
           origin_data = nan(length(xabs)+401, length(trialidx));
            noFiltBs = nan(length(xabs), length(trialidx)); %to store the baseline corrected unfiltered data
            filtBs = nan(length(xabs), length(trialidx)); %to store the baseline corrected filtered data

            powerstim = nan(length(trialidx),1025);
            freqstim = nan(length(trialidx),1025);
            fourhzpowerstim =nan(length(trialidx),1);
            % bsl = nan(1, length(trialidx));
            mean_wnd1 = nan(1,length(trialidx));

            all_pks = nan(4,length(SUA_data.channel_data.sdftr_chan(1,contrastBin)));

            for tridx = trialidx

                all_data = SUA_data.channel_data.sdftr_chan(401:1900,tridx);
                origin_data(:,tridx) = SUA_data.channel_data.sdftr_chan(:,tridx);
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
            trialnb = find(contrastBin);
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
                        locsdMUA_trial_struct = findpeaks_Loic(filtered_dMUA_high(200+ln:1499,trial));
                        locsdMUA_trial = locsdMUA_trial_struct.loc;

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
                            if ~isempty(findpeaks_Loic(filtered_dMUA_high(loc:200,tr)))
                                bsl_peak_locs_struct = findpeaks_Loic(filtered_dMUA_high(loc:200,tr));
                                bsl_peak_locs = bsl_peak_locs_struct.loc;
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
            trialnb = trialnb(~all(isnan(origin_data_high)));
            binary_data_high = SUA_data.channel_data.spk_bin_chan(:,trialnb); %get the binary spikes data
            origin_data_high = origin_data_high(:,~all(isnan(origin_data_high)));
            all_locsdMUA_trials =  all_locsdMUA_trials(:,~all(isnan(all_locsdMUA_trials)));

            %filtered_dMUA_high = filtered_dMUA_high(:,~all(isnan(filtered_dMUA_high))); % for nan - cols
            %all_pks = all_pks(:, ~all(isnan(all_pks)));
            %bsl_origin_data_high = bsl_origin_data_high(:,~all(isnan(bsl_origin_data_high)));


            %binNb = sprintf('bin%d', n);
            
            if n ==1 && length(origin_data_high(1,:)) >=10 %first bin == high contrast monocular condition will serve as an indicator of the minimum number of trials required for the analysis
                trialnum.after.(filename).(binNb) = nnz(trialnb);
                NoFiltMultiContSUA.(filename).(binNb) = origin_data_high;
                peakLocs.(filename).(binNb) = all_locsdMUA_trials; %create dynamical peak locations structures
                binSpk.(filename).(binNb) = binary_data_high;
                %FiltMultiContMUA.(filename).(binNb) =  filtered_dMUA_high;
                %BsNoFiltMultiContMUA.(filename).(binNb) = bsl_origin_data_high;

            elseif n == 1 && length(origin_data_high(1,:)) <10
                NoFiltMultiContSUA.(filename).(binNb) = [];
                peakLocs.(filename).(binNb) = [];
                binSpk.(filename).(binNb) = [];
                %all_pks(:,:) = [];
                % FiltMultiContMUA.(filename).(binNb) =  [];


            elseif n > 1 && length(origin_data_high(1,:)) >=10
                 trialnum.after.(filename).(binNb) = nnz(trialnb);
                NoFiltMultiContSUA.(filename).(binNb) = origin_data_high;
                peakLocs.(filename).(binNb) = all_locsdMUA_trials; %create dynamical peak locations structures
                binSpk.(filename).(binNb) = binary_data_high;
                %FiltMultiContMUA.(filename).(binNb) =  filtered_dMUA_high;
                % BsNoFiltMultiContMUA.(filename).(binNb) = bsl_origin_data_high;
            elseif n > 1 && length(origin_data_high(1,:)) <10
                NoFiltMultiContSUA.(filename).(binNb) = [];
                peakLocs.(filename).(binNb) = [];
                binSpk.(filename).(binNb) = [];
            end


            %data_peaks(i).namelist = all_pks(:,~all(isnan(all_pks)));
            %all_pks = all_pks(:,~all(isnan(all_pks)));
            %channelfilename = [unitsDir 'su_peaks_03032020_corrected\individual_units\' filename 'multiContrast'];
            %save(strcat(channelfilename, '.mat'), 'peakLocs');

        end
   end
end

savefilename = [unitsDir 'su_peaks_03032020_corrected\all_units\trial_num_before_after'];
save(strcat(savefilename, '.mat'), 'trialnum');


%% Plot the trial numbers before and after in the monocular and binocular condition
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'trial_num_before_after']; 
trialnums= load(channelfilename);


filenames = fieldnames(trialnums.trialnum.after);
times = [{'before'}, {'after'}];
bins = [1,6];
tn = [];
for i = 1:length(filenames)
    if isfield(trialnums.trialnum.after.(filenames{i}),'bin1') && isfield(trialnums.trialnum.after.(filenames{i}),'bin6')
        for t =1:2
            time = times{t};
            for b =1:2
                bin = sprintf('bin%d',bins(b));
                tn(i,t,b) = trialnums.trialnum.(time).(filenames{i}).(bin);
            end
        end
    end
end
tn = tn(tn(:,1) ~= 0,:);

%set up colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);
%set up data for gramm
trial_numbers = [tn(:,1);tn(:,2)];
condition = [repmat({'Monocular'},length(trial_numbers)/2,1); repmat({'Binocular'},length(trial_numbers)/2,1)];

%Plot
%1) Plot the scatter
clear g
g(1,1)=gramm('x',condition,'y',trial_numbers,'color',condition);

g(1,1).geom_jitter('width',0.4,'height',0); 
g(1,1).set_names('x','','y', '');
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]);
g(1,1).no_legend();

%Plot the mean and std for trial numbers
alltnMono = trial_numbers(strcmp(condition, 'Monocular'));
alltnBino = trial_numbers(strcmp(condition, 'Binocular'));
meanTn = nanmean([alltnMono,alltnBino],1);
shortCondition = unique(condition);
stdTns = std([alltnMono,alltnBino], 'omitnan');%/sqrt(length([alltnMono,alltnBino]));

ci_low = meanTn' - stdTns';
ci_high = meanTn' + stdTns';
g(1,1).update('x', shortCondition, 'y',meanTn,...
    'ymin',ci_low,'ymax',ci_high,'color',shortCondition);
%g(1,1).set_color_options('map','r');
g(1,1).geom_point('dodge',0.5);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_names('x','','y', 'Number of Trials');
g(1,1).no_legend();
g(1,1).set_point_options('base_size',7);

g(1,2) = gramm('x',ones(1,length(tn(:,1))), 'y', tn(:,1)-tn(:,2) );
g(1,2).geom_jitter('width',0.3,'height',0); 
g(1,2).axe_property('xlim', [0 2], 'ylim',[-30 30], 'XTickLabel','','XTick',''); 
g(1,2).set_names('x','','y', 'Difference Monocular-Binocular (Number of Trials)');
g.set_title('Number of trials in selected units post-processing');
g(1,2).set_color_options('map',[0.5 0.5 0.5]);


%Plot mean difference and +-std
meanDiff = nanmean([alltnMono-alltnBino],1);
stdDiffs = std([alltnMono-alltnBino], 'omitnan');%/sqrt(length([alltnMono-alltnBino]));
ci_low = meanDiff - stdDiffs;
ci_high = meanDiff + stdDiffs;
g(1,2).update('x', ones(1,1), 'y',meanDiff,...
    'ymin',ci_low,'ymax',ci_high);
g(1,2).geom_point('dodge',0.5);
g(1,2).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,2).axe_property('xlim', [0 2], 'ylim',[-20 30], 'XTickLabel','','XTick',''); 
%g(1,2).set_names('x','','y', '');
g(1,2).no_legend();
g(1,2).set_point_options('base_size',7);

figure('Position',[100 100 800 550]);

g.draw();

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_boxpointbar_std_numtrials_per_unit_old');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));



