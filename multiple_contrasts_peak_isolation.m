

%% First: use the list of single units file names that were selected in the adaptation analysis with
%high contrast
selectUnitsFilenames =load('C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;


unitsDir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);


cellClass = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
cellClass([1,46,55]) = [];


%% SECOND: isolate peak locations of smoothed data ==> use code from
%%"get_clean_peaks_and_data.m",
        % use multiple contrast levels
        % select trials the same way..?
 %% find peak locations of smoothed data, to further allow us to isolate peak values of unfiltered data in order to analyze them on R and fit a LMER
 %find the multiple contrast levels present in the data
 allContLevels =0;
 for i =1:71
     if ~isempty(filenames(i))
    contLevels = unique(unitsData.new_data(i).channel_data.contrast);
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
%data_peaks = struct();
peakLocs = struct(); %store filtered data peak locations used to isolate peak values of unfiltered data

for n = 1:length(contLims)-1
    clear i
     for i = channum
        if ~isempty(filenames{i})
       filename = filenames(i);

       blankcontrast = unitsData.new_data(i).channel_data.contrast ==  0 & unitsData.new_data(i).channel_data.fixedc ==  0;
      % highcontrast =unitsData.new_data(i).channel_data.contrast >=  0.5 & unitsData.new_data(i).channel_data.fixedc ==  0; 
       contrastBin = (unitsData.new_data(i).channel_data.contrast >  contLims(n) & unitsData.new_data(i).channel_data.contrast <= contLims(n+1))& unitsData.new_data(i).channel_data.fixedc ==  0; 

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
        mean_wnd1(tridx) = mean(lpdSUA(201:480));

        %%% power


        [powerstim(tridx,:), freqstim(tridx,:)] = calcFFT(all_data(200:1350));

        %find the index of the frequency vector closest to 4hz and point to the
        %power value of this index for every trial, and store the value in
        %fourhzpower
        [val,index] = min(abs(4-freqstim(tridx,:)));
        fourhzpowerstim(tridx,1) = powerstim(tridx,index);

         end

      %%%%%%%%%%% %reject trials below the 95%CI in the blank condition %%%%%%
       %power related variables
       power0 = fourhzpowerstim(blankcontrast); %power of responses in blank condition
       powerDE = fourhzpowerstim(contrastBin); %power of responses with contrast stimulus >0 in DE and 0 contrast in NDE

       %spiking activity related variables
       mean_wnd1_DE =mean_wnd1(contrastBin);
       filtered_dSUA_high = filtBs(:,contrastBin);
       filtered_dSUA_blank = filtBs(:,blankcontrast);
       origin_data_high = origin_data(:,contrastBin);
       origin_data_blank = origin_data(:,blankcontrast);
       %first peak location related variables 
       sua_bsl =  mean(filtered_dSUA_high(1:200,:),1);


       for tr = 1:length(powerDE)
          if mean_wnd1_DE(tr) > mean(sua_bsl)+1.96*std(sua_bsl)/sqrt(length(sua_bsl))  && powerDE(tr) > mean(power0)+1.96*std(power0)/sqrt(length(power0)) %/sqrt(length(sua_bsl)) /sqrt(length(power0))

              filtered_dSUA_high(:,tr) = filtered_dSUA_high(:,tr);
              origin_data_high(:,tr) = origin_data_high(:,tr);
           else

               filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
               origin_data_high(:,tr) =  nan(length(origin_data_high(:,tr)),1);
           end
       end

        %determine the first peak location for each trial of a given single
        %unit
        all_locsdSUA_trials = nan(6,length(filtered_dSUA_high(1,:)));
       clear trial
        for trial = 1:length(filtered_dSUA_high(1,:))

             for ln = 1:550
                 if filtered_dSUA_high(200+ln,trial) < filtered_dSUA_high(200+ln+1,trial) && ~all(isnan(filtered_dSUA_high(:,trial)))

                 locsdSUA_trial = findpeaks(filtered_dSUA_high(200+ln:1499,trial));
                 %if peak1 is too small, peak2 becomes peak1
                         if filtered_dSUA_high(locsdSUA_trial.loc(1)+200+ln,trial) >= 0.4*filtered_dSUA_high(locsdSUA_trial.loc(2)+200+ln)
                 %store first peak location 
                         all_locsdSUA_trials(1:length(locsdSUA_trial.loc),trial) = locsdSUA_trial.loc(1:end)+200+ln;
                         else
                         all_locsdSUA_trials(1:length(locsdSUA_trial.loc(2:end)),trial) = locsdSUA_trial.loc(2:end)+200+ln;

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
             else
        filtered_dSUA_high(:,trial) = nan(length(filtered_dSUA_high(:,trial)),1);
        all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
        origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);
            end 

            if ~all(isnan(all_locsdSUA_trials(:,trial))) && (all_locsdSUA_trials(4,trial) ~= 1500)
                 %adjust location to the first data point of lpsu (+ln),

        all_pks(:,trial) = filtered_dSUA_high(all_locsdSUA_trials(1:4,trial), trial);
        filtered_dSUA_high(:,trial) = filtered_dSUA_high(:,trial); 
        all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
        origin_data_high(:,trial) = origin_data_high(:,trial);
            else
        all_pks(:,trial) = nan(length(all_pks(:,trial)),1);      
        filtered_dSUA_high(:,trial) = nan(length(filtered_dSUA_high(:,trial)),1);
        all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
        origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);

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


       %reject if there is a peak 1 outlier, if the max peak value in the
       %baseline is an outlier

        % First find peaks before stimulus onset 

        bsl_peaks = nan(1, length(filtered_dSUA_high(1,:)));
        clear tr
        for tr = 1:length(filtered_dSUA_high(1,:))

             for loc = 1:200
                if filtered_dSUA_high(loc,tr) < filtered_dSUA_high(loc+1,tr) && ~all(isnan(filtered_dSUA_high(:,tr)))
               bsl_peak_locs = findpeaks(filtered_dSUA_high(loc:200,tr));
               bsl_peaks(1,tr) = max(filtered_dSUA_high(bsl_peak_locs.loc+loc,tr));
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

            else 
        filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
        all_pks(:,tr) = nan(length(all_pks(:,tr)),1);
        all_locsdSUA_trials(:,tr) = nan(size(all_locsdSUA_trials(:,tr)));
        origin_data_high(:,tr) = nan(length(origin_data_high(:,tr)),1);
            end
        end
       filtered_dSUA_high = filtered_dSUA_high(:,~all(isnan(filtered_dSUA_high))); % for nan - cols
       all_locsdSUA_trials =  all_locsdSUA_trials(:,~all(isnan(all_locsdSUA_trials)));
       all_pks = all_pks(:, ~all(isnan(all_pks)));
       origin_data_high = origin_data_high(:,~all(isnan(origin_data_high)));


      % if length(filtered_dSUA_high(1,:)) >=10

       %eval(['peakLocs.' num2str(i) '.bin' num2str(n) ' = all_locsdSUA_trials;'])
      
       filename = sprintf('x%s',char(filename));
       binNb = sprintf('bin%d', n);
       peakLocs.(filename).(binNb) = all_locsdSUA_trials; %create dynamical peak locations structures
       FiltMultiContSUA.(filename).(binNb) =  filtered_dSUA_high;
       FiltMultiContSUA.(filename).bin0 =  filtered_dSUA_blank;
       NoFiltMultiContSUA.(filename).(binNb) = origin_data_high;
       NoFiltMultiContSUA.(filename).bin0 = origin_data_blank;
       FiltMultiContSUA.(filename).cellclass = cellClass{i}; %get the cell class of selected units
       NoFiltMultiContSUA.(filename).cellclass = cellClass{i};
      % elseif length(filtered_dSUA_high(1,:)) <10  
       % all_pks(:,:) = [];
       % clean_high_SUA(i).namelist =  [];
       % clean_origin_data(i).unit = [];
       % eval(['peakLocs.unit' num2str(i) '.bin' num2str(n) '= [];']) 

       %end

     %data_peaks(i).namelist = all_pks(:,~all(isnan(all_pks)));
     %all_pks = all_pks(:,~all(isnan(all_pks)));
    %channelfilename = [unitsDir 'su_peaks_03032020_corrected\individual_units\' filename 'multiContrast'];
    %save(strcat(channelfilename, '.mat'), 'peakLocs');
        end
     end
end
allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\all_units\all_locs_data_95CI';
save(strcat(allfilename, '.mat'), 'peakLocs');


%allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\all_data_peaks'];
 %save(strcat(allfilename, '.mat'), 'data_peaks');
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\all_units\FiltMultiContSUA';
 save(strcat(allfilename, '.mat'), 'FiltMultiContSUA');
% allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\clean_SUA_locs'];
% save(strcat(allfilename, '.mat'), 'peaks_locs');
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\all_units\NoFiltMultiContSUA';
 save(strcat(allfilename, '.mat'), 'NoFiltMultiContSUA');
 

%% Third: isolate peak values of the origin data ==> use the code from
%%"get_origin_peaks"

%following the script get_clean_peaks_and_data.m (data cleaning and
%selection pipeline)
%this script was written to isolate peaks of the origin data in order to
%perform the statistical analysis of the peaks 
%some lines were commented out and replaced in order to also isolate
%normalized peak values, averaged across trials in order to plot the
%normalized average peak values of each unit on R
%Written by Loic Daumail, last edited on 6/29/2020

clear all

newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
channelfilename = [newdatadir 'FiltMultiContSUA']; 
FiltMultiContSUA = load(channelfilename);
locsfilename = [newdatadir 'all_locs_data_95CI'];
all_locsdSUA = load(locsfilename);

%gendatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
%channelfilename = [gendatadir 'refined_dataset']; 
%gen_data_file = load(channelfilename);


xabs = -199:1300;
nyq = 500;

channum = 1: length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA));
mean_origin_dSUA = struct();
mean_filtered_dSUA = struct();
suas_aligned_trials = struct();
peak_vals = struct();
%mean_peak_vals = struct();
mean_peaks =struct();
median_peaks = struct();
up_dist = nan(1, length(channum),4);
max_low_dist = struct();
all_locsdSUA_filtered = nan(1,length(channum),4);
%filenames = cell(length(channum),2);
contLims = [0,0.1,0.3,0.5,0.7,1];  
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);
for i = channum  
    for n = 1:length(contLims)-1
        
         binNb = sprintf('bin%d', n);
        filename = filenames{i};
        if ~isempty(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).(binNb))
            trialidx = 1:length(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).(binNb)(1,:));
            origin_dSUA = NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).(binNb)(401:1900,:); %- mean(data_file.clean_origin_data(i).unit(401:600,:),1);

            %create normalized origin trials data to plot average peaks for each unit with R
            %{
            norm_unit = nan(size(origin_dSUA));
            clear tr
            for tr =trialidx
                    min_unit =min(origin_dSUA(:,tr),[],1);
                    max_unit = max(origin_dSUA(:,tr),[],1);
                    norm_unit(:,tr) = (origin_dSUA(:,tr)-min_unit)./(max_unit - min_unit);
            end
            %}
            filtered_dSUA = FiltMultiContSUA.FiltMultiContSUA.(filename).(binNb);


            %determine the peak location of interest for each trial of a given single
            %unit
            all_locsdSUA_trials = all_locsdSUA.peakLocs.(filename).(binNb);

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
                   for n =trialidx
                          lower_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+1;
                          upper_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+length(xabs);

                          %origin data of the statistical analysis
                          fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = origin_dSUA(:,n);
                          %normalized data for the plotting
                         % fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = norm_unit(:,n);

                          filtered_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = filtered_dSUA(:,n);
                   end

             end
            %get the aligned data if it exists for the unit 
            suas_aligned_trials.(filename).(binNb)= fp_locked_trials;
            max_low_dist.(filename).(binNb) = max_low_dist_unit;


            clear pn
               for pn = 1:4
                   %peak data for the stats
              peak_vals.(filename).(binNb)(pn,:)= max(suas_aligned_trials.(filename).(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn), [],1);
               end
               %mean peaks for the R plots 
               mean_peaks.(filename).(binNb) = mean(peak_vals.(filename).(binNb),2);

               mean_peaks.(filename).cellclass = NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass;
               
               %median
               median_peaks.(filename).(binNb) = median(peak_vals.(filename).(binNb),2);

               median_peaks.(filename).cellclass = NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass;
         
                %peak data for the stats
            %  peak_vals(i).peak = [];
              %peak data for the R plots
             % mean_peaks(:,i) = nan(4,1);
        end
        mean_peaks.(filename).bin0 = mean(mean(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).bin0(401:1900,:),1));
        median_peaks.(filename).bin0 = median(median(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).bin0(401:1900,:),1));      
    end
 %filename = [gen_data_file.new_data(i).channel_data.filename, f{2}];
 %filename = erase(filename, '.mat');
 %filenames(i,1) = cellstr(filename);
 %filenames(i,2) = cellstr(layer(i));
 %peaks = peak_vals(i).peak;
%channelfilename = [gendatadir 'su_peaks_03032020_corrected\orig_peak_values\' filename];
%save(strcat(channelfilename, '.mat'), 'peaks');
end  
 %mean_peak_vals.peak = mean_peaks;
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\all_units\all_unfiltered_data_peaks';
 save(strcat(allfilename, '.mat'), 'mean_peaks');
 %allfilename = [gendatadir 'su_peaks_03032020_corrected\orig_peak_values\all_units\all_raw_mean_data_peaks'];
 %save(strcat(allfilename, '.mat'), 'mean_peaks');
 %savefilename = [gendatadir 'su_peaks_03032020_corrected\orig_peak_values\all_units\filenames_layers'];
 %save(strcat(savefilename, '.csv'), 'filenames');
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\all_units\all_unfiltered_median_peaks';
 save(strcat(allfilename, '.mat'), 'median_peaks');
 % post peak isolation count
 
 cellclass = {'M','P','K'};
 cnt = zeros(4,length(contLims)-1, length(cellclass));
 
for c = 1:length(cellclass)
    cellC =sprintf('%s',cellclass{c});

    for i =1:length(fieldnames(peak_vals))
        filename = filenames{i};

        if nnz(strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass,cellC))

            for bin = 1:length(contLims)-1
                binNb = sprintf('bin%d',bin);
                
                if nnz(strcmp(fieldnames(peak_vals.(filename)), binNb))
                    for pn =1:4
                        cnt(pn, bin, c) =cnt(pn, bin, c)+ numel(peak_vals.(filename).(binNb)(pn,:));


                    end
                end
            end
        end

    end
end
 
 %% Plot the contrast response curves of each cell class
 clear all
 
 estAvgCont = [0.01,0.05,0.2,0.4,0.6,0.85];
 %realAvgCont =  need to compute it from the actual contrast values
 
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);

meanPks = load([newdatadir 'all_unfiltered_data_peaks']);
 

class ={'M','P','K'};

for c = 1:length(class)
mean_pk1 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(estAvgCont)+1); 
mean_pk4 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(estAvgCont)+1); 

    for i =1:length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA))
         filename = filenames{i};
        if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
            for bin = 1:length(estAvgCont)
                binNb = sprintf('bin%d',bin);
                if nnz(strcmp(fieldnames(meanPks.mean_peaks.(filename)),binNb))
                    
                    mean_pk1(i,bin+1) = meanPks.mean_peaks.(filename).(binNb)(1);
                    mean_pk4(i,bin+1) = meanPks.mean_peaks.(filename).(binNb)(4);
                end
            end
             mean_pk1(i,1) = meanPks.mean_peaks.(filename).bin0(1);
             mean_pk4(i,1) = meanPks.mean_peaks.(filename).bin0(1);
        end
    end
    
    meanAllPk1 = nanmean(mean_pk1,1);
    meanAllPk4 = nanmean(mean_pk4,1);
    
    idxs = ~isnan(meanAllPk1);
    
    figure();
    semilogx(estAvgCont(idxs),meanAllPk1(idxs), 'Marker', 'o','col', 'r')
    %plot(estAvgCont(idxs),meanAllPk1(idxs), 'Marker', 'o','col', 'r')
    hold on
    sem1 = std(mean_pk1, 'omitnan')/sqrt(length(mean_pk1(~isnan(mean_pk1(:,1)),1)));
    h1= ciplot( meanAllPk1(idxs)+ 1.96*sem1(idxs),  meanAllPk1(idxs)- 1.96*sem1(idxs), estAvgCont(idxs),'r',0.1);
     set(h1, 'edgecolor','none') 
    hold on
    semilogx(estAvgCont(idxs),meanAllPk4(idxs), 'Marker', 'o', 'col', 'b')
    %plot(estAvgCont(idxs),meanAllPk4(idxs), 'Marker', 'o', 'col', 'b')
    hold on
    sem4 =std(mean_pk4, 'omitnan')/sqrt(length(mean_pk4(~isnan(mean_pk4(:,1)),1)));
    h2= ciplot( meanAllPk4(idxs)+ 1.96*sem4(idxs),  meanAllPk4(idxs)- 1.96*sem4(idxs), estAvgCont(idxs),'b',0.1); %[40/255 40/255 40/255]
      set(h2, 'edgecolor','none') 
  
    xlabel('Contrast')
    ylabel('Spike rate (spikes/sec)')
    title(sprintf('Average contrast response curve of %s cells',class{c}))
    legend('Pk1', '95%CI Pk1','Pk4','95%CI Pk4', 'Location', 'bestoutside')
    xlim([0 1])
    set(gca,'box','off')
    set(gca, 'linewidth',2)
    plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\plots\',strcat(sprintf('crf%s_cells_blank_log',class{c})));
   % saveas(gcf,strcat(plotdir, '.png'));
end
                
    
%% plot normalized crfs

 estAvgCont = [0.01,0.05,0.2,0.4,0.6,0.85];
 %realAvgCont =  need to compute it from the actual contrast values
 
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);

meanPks = load([newdatadir 'all_unfiltered_data_peaks']);
 

class ={'M','P','K'};

for c = 1:length(class)
mean_pk1 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(estAvgCont)+1); 
mean_pk4 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(estAvgCont)+1); 

    for i =1:length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA))
         filename = filenames{i};
        if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
             mean_pk1(i,1) = meanPks.mean_peaks.(filename).bin0(1);
             mean_pk4(i,1) = meanPks.mean_peaks.(filename).bin0(1);
            for bin = 1:length(estAvgCont)
                binNb = sprintf('bin%d',bin);
                if nnz(strcmp(fieldnames(meanPks.mean_peaks.(filename)),binNb))
                    
                    mean_pk1(i,bin+1) = meanPks.mean_peaks.(filename).(binNb)(1);
                    mean_pk4(i,bin+1) = meanPks.mean_peaks.(filename).(binNb)(4);
    
                end
            end
            
        end
    end
    
    normMeanPk1 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(estAvgCont)+1);
    normMeanPk4 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(estAvgCont)+1);
    for j =1:length(mean_pk1(:,1))
        normMeanPk1(j,:) = (mean_pk1(j,:)-min(mean_pk1(j,:),[],'omitnan'))/(max(mean_pk1(j,:),[],'omitnan')-min(mean_pk1(j,:),[],'omitnan'));
        normMeanPk4(j,:) = (mean_pk4(j,:)-min(mean_pk4(j,:)))/(max(mean_pk4(j,:))-min(mean_pk4(j,:)));
    end
        
    
    meanAllPk1 = nanmean(mean_pk1,1);
    meanAllPk4 = nanmean(mean_pk4,1);
    
    normMeanAllPk1 = (meanAllPk1 - min([meanAllPk1,meanAllPk4]))/(max([meanAllPk1,meanAllPk4])-min([meanAllPk1,meanAllPk4]));
    normMeanAllPk4 = (meanAllPk4 - min([meanAllPk1,meanAllPk4]))/(max([meanAllPk1,meanAllPk4])-min([meanAllPk1,meanAllPk4]));
   
    idxs = ~isnan(normMeanAllPk1);
    
    figure();
    %semilogx(estAvgCont(idxs),normMeanAllPk1(idxs), 'Marker', 'o','col', 'r')
    plot(estAvgCont(idxs),normMeanAllPk1(idxs), 'Marker', 'o','col', 'r')
    hold on
    sem1 = std(normMeanPk1, 'omitnan')/sqrt(length(normMeanPk1(~isnan(normMeanPk1(:,1)),1)));
    h1= ciplot( normMeanAllPk1(idxs)+ 1.96*sem1(idxs),  normMeanAllPk1(idxs)- 1.96*sem1(idxs), estAvgCont(idxs),'r',0.1);
      set(h1, 'edgecolor','none') 
    hold on
    %semilogx(estAvgCont(idxs),normMeanAllPk4(idxs), 'Marker', 'o', 'col', 'b')
    plot(estAvgCont(idxs),normMeanAllPk4(idxs), 'Marker', 'o', 'col', 'b')
    hold on
    sem4 =std(normMeanPk4, 'omitnan')/sqrt(length(normMeanPk4(~isnan(normMeanPk4(:,1)),1)));
    h2= ciplot( normMeanAllPk4(idxs)+ 1.96*sem4(idxs),  normMeanAllPk4(idxs)- 1.96*sem4(idxs), estAvgCont(idxs),'b',0.1); %[40/255 40/255 40/255]
      set(h2, 'edgecolor','none') 
  
    xlabel('Contrast')
    ylabel('Spike rate (spikes/sec)')
    title(sprintf('Average contrast response curve of %s cells',class{c}))
    legend('Pk1', '95%CI Pk1','Pk4','95%CI Pk4', 'Location', 'bestoutside')
    xlim([0 1])
    ylim([0 1])
    set(gca,'box','off')
    set(gca, 'linewidth',2)
    plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\plots\',strcat(sprintf('normalized_crf%s_cells_blank',class{c})));
   saveas(gcf,strcat(plotdir, '.png'));
end
 
 
%% plot median

estAvgCont = [0.01,0.05,0.2,0.4,0.6,0.85];
 %realAvgCont =  need to compute it from the actual contrast values
 
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);

medianPks = load([newdatadir 'all_unfiltered_median_peaks']);
 

class ={'M','P','K'};

for c = 1:length(class)
median_pk1 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(estAvgCont)+1); 
median_pk4 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(estAvgCont)+1); 

    for i =1:length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA))
         filename = filenames{i};
        if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
            for bin = 1:length(estAvgCont)
                binNb = sprintf('bin%d',bin);
                if nnz(strcmp(fieldnames(medianPks.median_peaks.(filename)),binNb))
                    
                    median_pk1(i,bin+1) = medianPks.median_peaks.(filename).(binNb)(1);
                    median_pk4(i,bin+1) = medianPks.median_peaks.(filename).(binNb)(4);
                end
            end
             median_pk1(i,1) = medianPks.median_peaks.(filename).bin0(1);
             median_pk4(i,1) = medianPks.median_peaks.(filename).bin0(1);
        end
    end
    
    medianAllPk1 = nanmedian(median_pk1,1);
    medianAllPk4 = nanmedian(median_pk4,1);
    
    idxs = ~isnan(medianAllPk1);
    
    figure();
    %semilogx(estAvgCont(idxs),medianAllPk1(idxs), 'Marker', 'o','col', 'r')
    plot(estAvgCont(idxs),medianAllPk1(idxs), 'Marker', 'o','col', 'r')
    hold on
    sem1 = std(median_pk1, 'omitnan')/sqrt(length(median_pk1(~isnan(median_pk1(:,1)),1)));
    h1= ciplot( medianAllPk1(idxs)+ 1.96*sem1(idxs),  medianAllPk1(idxs)- 1.96*sem1(idxs), estAvgCont(idxs),'r',0.1);
     set(h1, 'edgecolor','none') 
    hold on
    %semilogx(estAvgCont(idxs),medianAllPk4(idxs), 'Marker', 'o', 'col', 'b')
    plot(estAvgCont(idxs),medianAllPk4(idxs), 'Marker', 'o', 'col', 'b')
    hold on
    sem4 =std(median_pk4, 'omitnan')/sqrt(length(median_pk4(~isnan(median_pk4(:,1)),1)));
    h2= ciplot( medianAllPk4(idxs)+ 1.96*sem4(idxs),  medianAllPk4(idxs)- 1.96*sem4(idxs), estAvgCont(idxs),'b',0.1); %[40/255 40/255 40/255]
      set(h2, 'edgecolor','none') 
  
    xlabel('Contrast')
    ylabel('Spike rate (spikes/sec)')
    title(sprintf('Median contrast response curve of %s cells',class{c}))
    legend('Pk1', '95%CI Pk1','Pk4','95%CI Pk4', 'Location', 'bestoutside')
    xlim([0 1])
    set(gca,'box','off')
    set(gca, 'linewidth',2)
    plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\contrast_response_curves\plots\',strcat(sprintf('median_crf%s_cells_blank',class{c})));
    saveas(gcf,strcat(plotdir, '.png'));
end

