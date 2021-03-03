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


%% SECOND: isolate peak locations of smoothed data ==> use code from %%"get_clean_peaks_and_data.m",
        % use multiple contrast levels
        % select trials the same way..?
 %% find peak locations of smoothed data, to further allow us to isolate peak values of unfiltered data in order to analyze them on R and fit a LMER
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

           blankcontrast = unitsData.new_data(i).channel_data.contrast ==  0 & unitsData.new_data(i).channel_data.fixedc ==  0;
           if n == 1
               contrastBin =unitsData.new_data(i).channel_data.contrast >=  0.5 & unitsData.new_data(i).channel_data.fixedc ==  0; 
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
       bsl_origin_data_high = noFiltBs(:,contrastBin);
       %first peak location related variables 
       sua_bsl =  mean(filtered_dSUA_high(1:200,:),1);
       

       for tr = 1:length(powerDE)
          if mean_wnd1_DE(tr) > mean(sua_bsl)+1.96*std(sua_bsl)/sqrt(length(sua_bsl))  && powerDE(tr) > mean(power0)+1.96*std(power0)/sqrt(length(power0)) %/sqrt(length(sua_bsl)) /sqrt(length(power0))

              filtered_dSUA_high(:,tr) = filtered_dSUA_high(:,tr);
              origin_data_high(:,tr) = origin_data_high(:,tr);
              bsl_origin_data_high(:,tr) = bsl_origin_data_high(:,tr);
           else

               filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
               origin_data_high(:,tr) =  nan(length(origin_data_high(:,tr)),1);
               bsl_origin_data_high(:,tr) = nan(length(bsl_origin_data_high(:,tr)),1);
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


      % if length(filtered_dSUA_high(1,:)) >=10

       %eval(['peakLocs.' num2str(i) '.bin' num2str(n) ' = all_locsdSUA_trials;'])
      
       filename = sprintf('x%s',char(filename));
       binNb = sprintf('bin%d', n);
       peakLocs.(filename).(binNb) = all_locsdSUA_trials; %create dynamical peak locations structures
       FiltMultiContSUA.(filename).(binNb) =  filtered_dSUA_high;
      % FiltMultiContSUA.(filename).bin0 =  filtered_dSUA_blank;
       NoFiltMultiContSUA.(filename).(binNb) = origin_data_high;
       BsNoFiltMultiContSUA.(filename).(binNb) = bsl_origin_data_high;
       %NoFiltMultiContSUA.(filename).bin0 = origin_data_blank;
       FiltMultiContSUA.(filename).cellclass = cellClass{i}; %get the cell class of selected units
       NoFiltMultiContSUA.(filename).cellclass = cellClass{i};
       BsNoFiltMultiContSUA.(filename).cellclass = cellClass{i};
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
allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\all_locs_data_95CI';
save(strcat(allfilename, '.mat'), 'peakLocs');


%allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\all_data_peaks'];
 %save(strcat(allfilename, '.mat'), 'data_peaks');
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\FiltMultiContSUA';
 save(strcat(allfilename, '.mat'), 'FiltMultiContSUA');
% allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\clean_SUA_locs'];
% save(strcat(allfilename, '.mat'), 'peaks_locs');
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\NoFiltMultiContSUA';
 save(strcat(allfilename, '.mat'), 'NoFiltMultiContSUA');
 
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\BsNoFiltMultiContSUA';
 save(strcat(allfilename, '.mat'), 'BsNoFiltMultiContSUA');
 

%% Third: isolate peak values of the origin data ==> use the code from "get_origin_peaks.m"

%following the script get_clean_peaks_and_data.m (data cleaning and
%selection pipeline)
%this script was written to isolate peaks of the origin data in order to
%perform the statistical analysis of the peaks 
%some lines were commented out and replaced in order to also isolate
%normalized peak values, averaged across trials in order to plot the
%normalized average peak values of each unit on R
%Written by Loic Daumail, last edited on 6/29/2020

clear all

newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
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
%mean_filtered_dSUA = struct();
suas_aligned_trials = struct();
peak_aligned_trials = struct();
peak_vals = struct();
zscore_peak_vals = struct();
%mean_peak_vals = struct();
mean_peaks =struct();
median_peaks = struct();
zscore = struct();
up_dist = nan(1, length(channum),4);
max_low_dist = struct();
all_locsdSUA_filtered = nan(1,length(channum),4);
%filenames = cell(length(channum),2);
contLims = [0,0.1,0.3,0.5,0.7,1];  
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);
for i = channum  
    for n = 1:length(contLims)
        
         binNb = sprintf('bin%d', n);
        filename = filenames{i};
        if ~isempty(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).(binNb))
            trialidx = 1:length(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).(binNb)(1,:));
            origin_dSUA = NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).(binNb)(401:1900,:); %- mean(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).(binNb)(401:600,:),1);
            bs_origin_dSUA = NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).(binNb)(401:1900,:)- mean(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).(binNb)(401:600,:),1);
           
            %origin_dSUA = NoFiltMultiContSUA.BsNoFiltMultiContSUA.(filename).(binNb)(1:end,:);
            %create zscore origin trials data to plot average peaks for each unit with R
            
            zscore_unit = nan(size(origin_dSUA));
            clear tr
            for tr =trialidx
                    
                   zscore_unit(:,tr) = (origin_dSUA(:,tr)-mean(origin_dSUA(:,tr)))./(std(origin_dSUA(:,tr)));
            end
            
             
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
            bs_fp_locked_trials = nan(new_dist_unit,length(origin_dSUA(1,:)),4);
            zscore_fp_locked_trials = nan(new_dist_unit,length(origin_dSUA(1,:)),4);
            %filtered_fp_locked_trials = nan(new_dist_unit,length(filtered_dSUA(1,:)),4);
             clear n pn
             for pn =1:4
                   for n =trialidx
                          lower_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+1;
                          upper_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+length(xabs);

                          %origin data of the statistical analysis
                          fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = origin_dSUA(:,n);
                          %baseline corrected data
                          bs_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = bs_origin_dSUA(:,n);
                          %zscore transformed data for the plotting
                          zscore_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = zscore_unit(:,n);

                          %filtered_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = filtered_dSUA(:,n);
                   end

             end
            %get the aligned data if it exists for the unit 
            suas_aligned_trials.(filename).origin.(binNb)= fp_locked_trials;
            suas_aligned_trials.(filename).bsorigin.(binNb)= bs_fp_locked_trials;
            suas_aligned_trials.(filename).zscore.(binNb)= zscore_fp_locked_trials;
             
            max_low_dist.(filename).(binNb) = max_low_dist_unit;

            clear pn
               for pn = 1:4
                   %peak data for the stats
                  peak_vals.(filename).(binNb)(pn,:)= max(suas_aligned_trials.(filename).origin.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn), [],1);
                  zscore_peak_vals.(filename).(binNb)(pn,:)= max(suas_aligned_trials.(filename).zscore.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn), [],1);
            
                  
                  %peak aligne trials
                pkNb = sprintf('pk%d', pn);
                peak_aligned_trials.(filename).origin.(binNb).(pkNb) = suas_aligned_trials.(filename).origin.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn);
                peak_aligned_trials.(filename).bsorigin.(binNb).(pkNb) = suas_aligned_trials.(filename).bsorigin.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn);
                peak_aligned_trials.(filename).zscore.(binNb).(pkNb) = suas_aligned_trials.(filename).zscore.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn);
                peak_aligned_trials.(filename).cellclass = NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass;


               end
               
               %mean peaks for the R plots 
               mean_peaks.(filename).(binNb) = mean(peak_vals.(filename).(binNb),2);

               mean_peaks.(filename).cellclass = NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass;
               
               %median
               median_peaks.(filename).(binNb) = median(peak_vals.(filename).(binNb),2);

               median_peaks.(filename).cellclass = NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass;
               %zscore
               %zscore.(filename).(binNb) = (peak_vals.(filename).(binNb)-repmat(mean(peak_vals.(filename).(binNb),2), 1,length(peak_vals.(filename).(binNb)(1,:))))./repmat(std(peak_vals.(filename).(binNb),0,2),1,length(peak_vals.(filename).(binNb)(1,:)));
             
               zscore_peak_vals.(filename).cellclass = NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass;
   
        end
     end

end  
 %mean_peak_vals.peak = mean_peaks;
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\all_unfiltered_bsdata_peaks';
 save(strcat(allfilename, '.mat'), 'mean_peaks');
 %allfilename = [gendatadir 'su_peaks_03032020_corrected\orig_peak_values\all_units\all_raw_mean_data_peaks'];
 %save(strcat(allfilename, '.mat'), 'mean_peaks');
 %savefilename = [gendatadir 'su_peaks_03032020_corrected\orig_peak_values\all_units\filenames_layers'];
 %save(strcat(savefilename, '.csv'), 'filenames');
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\all_unfiltered_bsmedian_peaks';
 save(strcat(allfilename, '.mat'), 'median_peaks');
 
 %save z-score transformed data
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\all_unfiltered_zscore_peaks';
 save(strcat(allfilename, '.mat'), 'zscore_peak_vals');
 
 allfilename = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\all_orig_bs_zscore_trials';
 save(strcat(allfilename, '.mat'), 'peak_aligned_trials');

 % post peak isolation count
 
 cellclass = {'M','P','K'};
 scnt = zeros(4,length(contLims), length(cellclass));
 
for c = 1:length(cellclass)
    cellC =sprintf('%s',cellclass{c});

    for i =1:length(fieldnames(peak_vals))
        filename = filenames{i};

        if nnz(strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass,cellC))

            for bin = 1:length(contLims)
                binNb = sprintf('bin%d',bin);
                
                if nnz(strcmp(fieldnames(peak_vals.(filename)), binNb))
                    for pn =1:4
                        scnt(pn, bin, c) =scnt(pn, bin, c)+ numel(peak_vals.(filename).(binNb)(pn,:));


                    end
                end
            end
        end

    end
end
 
 %% Plot the comparisons of responses pkn of monocular stim vs pkn of binocular stim (varying stim contrast in NDE) of each cell class
 clear all
 
 NdeAvgCont = [0,0.85];
 %realAvgCont =  need to compute it from the actual contrast values
 
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);

meanPks = load([newdatadir 'all_unfiltered_data_peaks']);
ylims = [[0 270];[0 230];[0 250]];
 

class ={'M','P','K'};
bins =[1,6];
for c = 1:length(class)
mean_pk1 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(NdeAvgCont)); 
mean_pk4 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(NdeAvgCont)); 

    for i =1:length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA))
         filename = filenames{i};
        if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
            for bin = 1:length(NdeAvgCont)
                binNb = sprintf('bin%d',bins(bin));
                if nnz(strcmp(fieldnames(meanPks.mean_peaks.(filename)),binNb))
                    
                    mean_pk1(i,bin) = meanPks.mean_peaks.(filename).(binNb)(1);
                    mean_pk4(i,bin) = meanPks.mean_peaks.(filename).(binNb)(4);
                end
            end
        end
    end
    
    meanAllPk1 = nanmean(mean_pk1,1);
    meanAllPk4 = nanmean(mean_pk4,1);
    
    idxs = ~isnan(meanAllPk1);
    
    figure();
    h1 = boxplot([mean_pk1(:,1) mean_pk1(:,2)] ,'notch','off','labels',{'NDE = 0','NDE =0.85'}, 'Color', 'b');
    hold on
   
    h2= boxplot([mean_pk4(:,1) mean_pk4(:,2)] ,'notch','off','labels',{'NDE = 0','NDE =0.85'}, 'Color', 'r');
    
    ylim(ylims(c,:))
    legend([h1(3),h2(3)], {'Pk1','Pk4'},'Location', 'bestoutside')
    xlabel('Contrast')
    ylabel('Spike rate (spikes/sec)')
    title(sprintf('Average response of %s cells in moocular vs binocular stimulation conditions',class{c}))  
    set(gca,'box','off')
    set(gca, 'linewidth',2)
    plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat(sprintf('comparison_mono_bino_highcontrast_%s_cells',class{c})));
    saveas(gcf,strcat(plotdir, '.png'));
end
                
    
%% plot only adapting units (overal mean or all units)
 NdeAvgCont = [0,0.85];
 %realAvgCont =  need to compute it from the actual contrast values
 
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);

meanPks = load([newdatadir 'all_unfiltered_data_peaks']);
pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1)';
pvalues = pvalues(:,~all(isnan(pvalues)));

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

all_mean_data = nan(4, length(peakvals.peak_vals));
for nunit = 1:length(peakvals.peak_vals)
 if ~isempty(peakvals.peak_vals(nunit).peak)
 mean_data = nanmean(peakvals.peak_vals(nunit).peak,2);
   all_mean_data(:,nunit) = mean_data;
 end
end
all_mean_data =   all_mean_data(:,~all(isnan(all_mean_data)));

 

%% Plot difference pk1 - pk4 of adapting unitsclass ={'M','P','K'};
bins =[1,6];
for c = 1:length(class)
mean_pk1 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(NdeAvgCont)); 
mean_pk4 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(NdeAvgCont)); 
  
 
    for i =1:length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA))
        if all_mean_data(4,i) < all_mean_data(1,i) && pvalues(4,i) < .05
             filename = filenames{i};
            if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
                for bin = 1:length(NdeAvgCont)
                    binNb = sprintf('bin%d',bins(bin));
                    if nnz(strcmp(fieldnames(meanPks.mean_peaks.(filename)),binNb))

                        mean_pk1(i,bin) = meanPks.mean_peaks.(filename).(binNb)(1);
                        mean_pk4(i,bin) = meanPks.mean_peaks.(filename).(binNb)(4);
                    end
                end
            end
        end
    end
    
    meanAllPk1 = nanmean(mean_pk1,1);
    meanAllPk4 = nanmean(mean_pk4,1);
    
    idxs = find(~isnan(mean_pk1(:,1)));
    %{
    figure();
    meanAllBin1 =[meanAllPk1(1,1),meanAllPk4(1,1)]';
    mean_bin1 =[mean_pk1(:,1),mean_pk4(:,1)]';
    h1= plot(meanAllBin1,'-o', 'Color', 'b');
    %hold on
    %sem1 =std(mean_bin1,0 ,2, 'omitnan')/sqrt(length(mean_bin1(1,~isnan(mean_bin1(1,:)))));
    %ci1= ciplot( meanAllBin1(idxs)+ 1.96*sem1(idxs),  meanAllBin1(idxs)- 1.96*sem1(idxs), [1 2],'b',0.1); %[40/255 40/255 40/255]
 
    hold on
    meanAllBin2 =[meanAllPk1(1,2),meanAllPk4(1,2)]';
    mean_bin2 =[mean_pk1(:,2),mean_pk4(:,2)]';
  
   %h2= plot([meanAllPk1(1,2),meanAllPk4(1,2)]','-o','Color', 'r');
   %sem2 =std(mean_bin2,0 ,2, 'omitnan')/sqrt(length(mean_bin2(1,~isnan(mean_bin2(1,:)))));
   %ci2= ciplot( meanAllBin2(idxs)+ 1.96*sem2(idxs),  meanAllBin2(idxs)- 1.96*sem2(idxs), [1 2],'r',0.1); %[40/255 40/255 40/255]
 %}
    figure();
    
    for u =1:length(idxs)
        subplot(2, ceil(length(idxs)/2), u)
        h1 =plot([mean_pk1(idxs(u),1),mean_pk4(idxs(u),1)],'-o', 'Color', 'b');
        hold on
        h2= plot([mean_pk1(idxs(u),2),mean_pk4(idxs(u),2)],'-o', 'Color', 'r');
        xticks([0,1,2,3])
        xticklabels({'','Pk1','Pk4',''})
        xlim([0 3])
        ylabel('Spike rate (spikes/sec)')
        xlabel('Peak #')
        legend([h1(1), h2(1)],'Monocular', 'Binocular', 'location', 'bestoutside')
        set(gca,'box','off')
        set(gca, 'linewidth',2)
        title(sprintf('Monocular vs Binocular stimulation response of adapting %s cell',class{c}))  
        plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat(sprintf('comparison_mono_bino_highcontrast_adapting_mean_%s_cells',class{c})));
       % saveas(gcf,strcat(plotdir, '.png'));
    end
end 


newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);

meanPks = load([newdatadir 'all_unfiltered_bsdata_peaks']);
nBsMeanPks = load([newdatadir 'all_unfiltered_data_peaks']);
%ylims = [[0 270];[0 230];[0 250]];
 
pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1)';
pvalues = pvalues(:,~all(isnan(pvalues)));

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

all_mean_data = nan(4, length(peakvals.peak_vals));
for nunit = 1:length(peakvals.peak_vals)
 if ~isempty(peakvals.peak_vals(nunit).peak)
 mean_data = nanmean(peakvals.peak_vals(nunit).peak,2);
   all_mean_data(:,nunit) = mean_data;
 end
end
all_mean_data =   all_mean_data(:,~all(isnan(all_mean_data)));

 class ={'M','P','K'};
bins =[1,6];
for c = 1:length(class)
mean_pk1 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(bins)); 
mean_pk4 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(bins)); 
  
 
    for i =1:length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA))
        if all_mean_data(4,i) < all_mean_data(1,i) && pvalues(4,i) < .05
             filename = filenames{i};
            if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
                for bin = 1:length(bins)
                    binNb = sprintf('bin%d',bins(bin));
                    if nnz(strcmp(fieldnames(meanPks.mean_peaks.(filename)),binNb))

                        mean_pk1(i,bin) = meanPks.mean_peaks.(filename).(binNb)(1);
                        mean_pk4(i,bin) = meanPks.mean_peaks.(filename).(binNb)(4);
                    end
                end
            end
        end
    end
    
    diffBin1 = mean_pk1(:,1) -mean_pk4(:,1);
    diffBin2 = mean_pk1(:,2) -mean_pk4(:,2);
    
    meanDiffBin1 = nanmean(diffBin1);
    meanDiffBin2 = nanmean(diffBin2);
    idxs = ~isnan(diffBin1);
    
    figure();
    %x = ones(length(diffBin1(idxs)),1);
    x = ones(1,1);
    h1= plot([x,2*x]',[diffBin1(idxs), diffBin2(idxs)]', '-o');
   % h1= plot([x,2*x]',[meanDiffBin1, meanDiffBin2]', '-o');
    %hold on 
    %sem1 =std([diffBin1,diffBin2],0 ,1, 'omitnan')/sqrt(length(diffBin1(~isnan(diffBin1(:,1)))));
    %ci1= ciplot( [meanDiffBin1, meanDiffBin2]+ 1.96*sem1,  [meanDiffBin1, meanDiffBin2]- 1.96*sem1, [1 2],'b',0.1); %[40/255 40/255 40/255]
 
   % h2 = scatter(2*x,diffBin2)
    xticks([0,1,2,3])
    xticklabels({'','Monocular','Binocular',''})
    xlim([0 3])
    ylabel('Pk1-Pk4 spike rate difference (spikes/sec)')
    xlabel('Condition')
   % legend([h1(1), h1(2)],'Monocular', 'Binocular', 'location', 'bestoutside')
    set(gca,'box','off')
    set(gca, 'linewidth',2)
    title(sprintf('Monocular vs Binocular stimulation response of adapting %s cells',class{c}))  
    plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat(sprintf('comparison_mono_bino_highcontrast_adapting_meandiffBSpk1pk4_%s_cells',class{c})));
    %saveas(gcf,strcat(plotdir, '.png'));
end 


%% Plot the difference pk1 - pk4 of adapting units with z-score transformed data

newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);

zscorePks = load([newdatadir 'all_unfiltered_zscore_peaks']);

%ylims = [[0 270];[0 230];[0 250]];
 
pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1)';
pvalues = pvalues(:,~all(isnan(pvalues)));

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

all_mean_data = nan(4, length(peakvals.peak_vals));
for nunit = 1:length(peakvals.peak_vals)
 if ~isempty(peakvals.peak_vals(nunit).peak)
 mean_data = nanmean(peakvals.peak_vals(nunit).peak,2);
   all_mean_data(:,nunit) = mean_data;
 end
end
all_mean_data =   all_mean_data(:,~all(isnan(all_mean_data)));

 class ={'M','P','K'};
bins =[1,6];
for c = 1:length(class)
mean_pk1 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(bins)); 
mean_pk4 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(bins)); 
  
 
    for i =1:length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA))
        if all_mean_data(4,i) < all_mean_data(1,i) && pvalues(4,i) < .05
             filename = filenames{i};
            if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
                for bin = 1:length(bins)
                    binNb = sprintf('bin%d',bins(bin));
                    if nnz(strcmp(fieldnames(zscorePks.zscore_peak_vals.(filename)),binNb))

                        mean_pk1(i,bin) = mean(zscorePks.zscore_peak_vals.(filename).(binNb)(1,:));
                        mean_pk4(i,bin) = mean(zscorePks.zscore_peak_vals.(filename).(binNb)(4,:));
                    end
                end
            end
        end
    end
    
    diffBin1 = mean_pk1(:,1) -mean_pk4(:,1);
    diffBin2 = mean_pk1(:,2) -mean_pk4(:,2);
    
    meanDiffBin1 = nanmean(diffBin1);
    meanDiffBin2 = nanmean(diffBin2);
    idxs = ~isnan(diffBin1);
    
    figure();
    %x = ones(length(diffBin1(idxs)),1);
    x = ones(1,1);
    %h1= plot([x,2*x]',[diffBin1(idxs), diffBin2(idxs)]', '-o');
    h1= plot([x,2*x]',[meanDiffBin1, meanDiffBin2]', '-o');
    hold on 
    sem1 =std([diffBin1,diffBin2],0 ,1, 'omitnan')/sqrt(length(diffBin1(~isnan(diffBin1(:,1)))));
    ci1= ciplot( [meanDiffBin1, meanDiffBin2]+ 1.96*sem1,  [meanDiffBin1, meanDiffBin2]- 1.96*sem1, [1 2],'b',0.1); %[40/255 40/255 40/255]
 
   % h2 = scatter(2*x,diffBin2)
    xticks([0,1,2,3])
    xticklabels({'','Monocular','Binocular',''})
    xlim([0 3])
    ylabel('Pk1-Pk4 spike rate difference (z-score transformed)')
    xlabel('Condition')
   % legend([h1(1), h1(2)],'Monocular', 'Binocular', 'location', 'bestoutside')
    set(gca,'box','off')
    set(gca, 'linewidth',2)
    title(sprintf('Monocular vs Binocular stimulation response of adapting %s cells',class{c}))  
    plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat(sprintf('comparison_mono_bino_highcontrast_adapting_allmeanzscorediffpk1pk4_%s_cells',class{c})));
    %saveas(gcf,strcat(plotdir, '.png'));
end 

%% Plot the difference pk1 - pk4 of adapting units with z-score transformed data unit by unit

newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);

zscorePks = load([newdatadir 'all_unfiltered_zscore_peaks']);

%ylims = [[0 270];[0 230];[0 250]];
 
pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1)';
pvalues = pvalues(:,~all(isnan(pvalues)));

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

all_mean_data = nan(4, length(peakvals.peak_vals));
for nunit = 1:length(peakvals.peak_vals)
 if ~isempty(peakvals.peak_vals(nunit).peak)
 mean_data = nanmean(peakvals.peak_vals(nunit).peak,2);
   all_mean_data(:,nunit) = mean_data;
 end
end
all_mean_data =   all_mean_data(:,~all(isnan(all_mean_data)));

 class ={'M','P','K'};
bins =[1,6];
for c = 1:length(class)
mean_pk1 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(bins)); 
mean_pk4 = nan(length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA)),length(bins)); 
  
 
    for i =1:length(fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA))
        if all_mean_data(4,i) < all_mean_data(1,i) && pvalues(4,i) < .05
             filename = filenames{i};
            if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
                for bin = 1:length(bins)
                    binNb = sprintf('bin%d',bins(bin));
                    if nnz(strcmp(fieldnames(zscorePks.zscore_peak_vals.(filename)),binNb))

                        mean_pk1(i,bin) = mean(zscorePks.zscore_peak_vals.(filename).(binNb)(1,:));
                        mean_pk4(i,bin) = mean(zscorePks.zscore_peak_vals.(filename).(binNb)(4,:));
                    end
                end
            end
        end
    end
    
    diffBin1 = mean_pk1(:,1) -mean_pk4(:,1);
    diffBin2 = mean_pk1(:,2) -mean_pk4(:,2);

    idxs = find(~isnan(diffBin1));
    
    for u =1:length(idxs)
    figure();
    %x = ones(length(diffBin1(idxs)),1);
    x = ones(1,1);
    
    h1= plot([x,2*x]',[diffBin1(idxs(u)), diffBin2(idxs(u))]', '-o');
    %h1= plot([x,2*x]',[meanDiffBin1, meanDiffBin2]', '-o');
    %hold on 
    %sem1 =std([diffBin1,diffBin2],0 ,1, 'omitnan')/sqrt(length(diffBin1(~isnan(diffBin1(:,1)))));
    %ci1= ciplot( [meanDiffBin1, meanDiffBin2]+ 1.96*sem1,  [meanDiffBin1, meanDiffBin2]- 1.96*sem1, [1 2],'b',0.1); %[40/255 40/255 40/255]
 
   % h2 = scatter(2*x,diffBin2)
    xticks([0,1,2,3])
    xticklabels({'','Monocular','Binocular',''})
    xlim([0 3])
    ylabel('Pk1-Pk4 spike rate difference (z-score transformed)')
    xlabel('Condition')
   % legend([h1(1), h1(2)],'Monocular', 'Binocular', 'location', 'bestoutside')
    set(gca,'box','off')
    set(gca, 'linewidth',2)
    title(sprintf('Monocular vs Binocular stimulation response of adapting %s cells',class{c}))  
    plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat(sprintf('comparison_mono_bino_highcontrast_adapting_meanzscorediffpk1pk4_%s_cell_%d',class{c},u)));
    saveas(gcf,strcat(plotdir, '.png'));
    end
end 




%% Plot 3 subplots for every unit: Peak aligned origin data , bs corrected, zscore transformed 

newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';

zscorePks = load([newdatadir 'all_unfiltered_zscore_peaks']);
trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials']);

filenames = fieldnames(trialsTraces.peak_aligned_trials);

pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1)';
pvalues = pvalues(:,~all(isnan(pvalues)));

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

all_mean_data = nan(4, length(peakvals.peak_vals));
for nunit = 1:length(peakvals.peak_vals)
 if ~isempty(peakvals.peak_vals(nunit).peak)
 mean_data = nanmean(peakvals.peak_vals(nunit).peak,2);
   all_mean_data(:,nunit) = mean_data;
 end
end
all_mean_data =   all_mean_data(:,~all(isnan(all_mean_data)));

 class ={'M','P','K'};
 ylims = [[0 270];[0 230];[0 250]];
 bsylims = [[-30 250];[-20 210];[-20 230]];
 zylim = [-1 3];
bins =[1,6];

for c = 1:length(class)
    
    for i =1:length(filenames)
        if all_mean_data(4,i) < all_mean_data(1,i) && pvalues(4,i) < .05
            filename = filenames{i};
           
            if strcmp(trialsTraces.peak_aligned_trials.(filename).cellclass, class{c})
                        bin1 = sprintf('bin%d',1);
                        pk1 = sprintf('pk%d', 1);
                        
                         origTrace = repmat(nan(size(squeeze(mean(trialsTraces.peak_aligned_trials.(filename).origin.(bin1).(pk1),2)))),1,4,2);
                         bsTrace =repmat(nan(size(squeeze(mean(trialsTraces.peak_aligned_trials.(filename).bsorigin.(bin1).(pk1),2)))),1,4,2);
                         zTrace = repmat(nan(size(squeeze(mean(trialsTraces.peak_aligned_trials.(filename).zscore.(bin1).(pk1),2)))),1,4,2);
                
                for bin = 1:length(bins)
                    binNb = sprintf('bin%d',bins(bin));
                    
                        
                    for pn =1:4
                         pkNb = sprintf('pk%d', pn);
                        
                
                        if nnz(strcmp(fieldnames(zscorePks.zscore_peak_vals.(filename)),binNb))
                             origTrace(:,pn,bin) = squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).origin.(binNb).(pkNb),2)); 
                             bsTrace(:,pn,bin) = squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).bsorigin.(binNb).(pkNb),2)); 
                             zTrace(:,pn,bin) =  squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).zscore.(binNb).(pkNb),2)); 

                        end
                    end
                end
    
    
                figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);
                for p =1:4
                    %monocular
                    subplot(3,8,p)
                    plot(origTrace(:,p,1),'Color','b')
                    if p == 1
                    title(sprintf('Monocular stimulation %s',class{c}))
                    ylabel('Spike rate (spikes/s)')
                    end
                    set(gca,'box','off')
                    ylim(ylims(c,:))
                    
                    if p>1
                         ax1 = gca;                   
                         ax1.YAxis.Visible = 'off';   
                    end
                    
                    subplot(3,8,4+p)
                    plot(origTrace(:,p,2),'Color','b')
                    if p == 1
                    title(sprintf('Binocular stimulation %s',class{c}))
                    ylabel('Spike rate (spikes/s)')
                    end
                    set(gca,'box','off')
                    ylim(ylims(c,:))
                    if p>1
                         ax1 = gca;                   
                         ax1.YAxis.Visible = 'off';   
                    end
              
                    subplot(3,8,8+p)
                    plot(bsTrace(:,p,1),'Color','r')
                    if p == 1
                    title(sprintf('Baseline corrected monocular stimulation %s',class{c}))
                    
                    ylabel('Spike rate ( delta spikes/s)')
                    end
                    set(gca,'box','off')
                    ylim(bsylims(c,:))
                    
                    if p>1
                         ax1 = gca;                   
                         ax1.YAxis.Visible = 'off';   
                    end
                   
                    subplot(3,8,12+p)
                    plot(bsTrace(:,p,2),'Color','r')
                    if p == 1
                    title(sprintf('Baseline corrected binocular stimulation %s',class{c}))
                    ylabel('Spike rate ( delta spikes/s)')
                    end
                     set(gca,'box','off')
                     ylim(bsylims(c,:))
                    if p>1
                         ax1 = gca;                   
                         ax1.YAxis.Visible = 'off';   
                    end
                    
                    subplot(3,8,16+p)
                    plot(zTrace(:,p,1),'Color','k')
                    if p == 1
                    title(sprintf('Z-score transformed monocular stimulation %s',class{c}))
                    
                    ylabel('Spike rate (z-score transformed)')
                    end
                     set(gca,'box','off')
                     ylim(zylim)
                    
                    if p>1
                         ax1 = gca;                   
                         ax1.YAxis.Visible = 'off';   
                    end
                     

                    %binocular

                    subplot(3,8,20+p)
                    plot(zTrace(:,p,2),'Color','k')
                    
                    if p == 1
                    title(sprintf('Z-score transformed binocular stimulation %s',class{c}))
                    ylabel('Spike rate (z-score transformed)')
                    end 
                    
                    set(gca,'box','off')
                    ylim(zylim)
                    
                    if p>1
                         ax1 = gca;                   
                         ax1.YAxis.Visible = 'off';   
                    end
                        
                    
                end
                 
                 plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat(sprintf('mono_bino_orig_bs_zscore_%s_cell%d',class{c},i)));
                % saveas(gcf,strcat(plotdir, '.png'));
                
            end
        
        end
    end
end


%% Plot overall mean of binocular condition vs overall mean of monocular condition

newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';

zscorePks = load([newdatadir 'all_unfiltered_zscore_peaks']);
trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials']);

filenames = fieldnames(trialsTraces.peak_aligned_trials);
pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1)';
pvalues = pvalues(:,~all(isnan(pvalues)));

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

all_mean_data = nan(4, length(peakvals.peak_vals));
for nunit = 1:length(peakvals.peak_vals)
 if ~isempty(peakvals.peak_vals(nunit).peak)
 mean_data = nanmean(peakvals.peak_vals(nunit).peak,2);
   all_mean_data(:,nunit) = mean_data;
 end
end
all_mean_data =   all_mean_data(:,~all(isnan(all_mean_data)));


 class ={'M','P','K'};
%ylims = [[0 160];[0 105];[0 210]];
ylims = [[0 200];[0 140];[0 210]];
 bsylims = [[-30 250];[-20 210];[-20 230]];
 zylim = [-1 3];
bins =[1,6];
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue

for c = 1:length(class)
    
         origTrace = nan(250,4,length(filenames),2);
         bsTrace =nan(250,4,length(filenames),2);
         zTrace = nan(250,4,length(filenames),2);
    for i =1:length(filenames)
       if all_mean_data(4,i) < all_mean_data(1,i) && pvalues(4,i) < .05
            filename = filenames{i};
           
            if strcmp(trialsTraces.peak_aligned_trials.(filename).cellclass, class{c})

                for bin = 1:length(bins)
                    binNb = sprintf('bin%d',bins(bin));
                    
                        
                    for pn =1:4
                         pkNb = sprintf('pk%d', pn);
                        
                
                        if nnz(strcmp(fieldnames(zscorePks.zscore_peak_vals.(filename)),binNb))
                             origTrace(:,pn,i,bin) = squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).origin.(binNb).(pkNb),2)); 
                             bsTrace(:,pn,i,bin) = squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).bsorigin.(binNb).(pkNb),2)); 
                             zTrace(:,pn,i,bin) =  squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).zscore.(binNb).(pkNb),2)); 

                        end
                    end
                end
            end
       end
    end
    meanOrig = squeeze(nanmean(origTrace,3));
    monCi_low = meanOrig(:,:,1) - 1.96*std(origTrace(:,:,:,1),0,3, 'omitnan')./sqrt(length(find(~isnan(origTrace(1,1,:,1)))));
    monCi_high = meanOrig(:,:,1) + 1.96*std(origTrace(:,:,:,1),0,3, 'omitnan')./sqrt(length(find(~isnan(origTrace(1,1,:,1)))));
   
    binCi_low = meanOrig(:,:,2) - 1.96*std(origTrace(:,:,:,2),0,3, 'omitnan')./sqrt(length(find(~isnan(origTrace(1,1,:,1)))));
    binCi_high = meanOrig(:,:,2) + 1.96*std(origTrace(:,:,:,2),0,3, 'omitnan')./sqrt(length(find(~isnan(origTrace(1,1,:,1)))));
   
   
                figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);
                for p =1:4
                    
                    
                       h= subplot(1,4,p);
                    
                       %binocular
                       plot(-125:124,meanOrig(:,p,2),'LineWidth',2, 'Color',[40/255 40/255 40/255])
                     hold on
                      h1= ciplot(binCi_low(:,p), binCi_high(:,p),[-125:124],col(c,:),0.5);
                      set(h1, 'edgecolor','none') 
                      hold on
                         %monocular
                       plot(-125:124,meanOrig(:,p,1),'LineWidth',2, 'Color','green')
                
                    set(h,'position',get(h,'position').*[1 1 1.15 1])
                    ylim(ylims(c,:))
                    if p == 1
                        
                        ylabel('Spike rate (spikes/s)')
                    end
                    xlim([-125 125])
                    set(gca, 'linewidth',2)
                    set(gca,'box','off')
                    
                    if p>1
                         ax1 = gca;                   
                         ax1.YAxis.Visible = 'off';   
                    end               
                end
                legend('Binocular','','Monocular', 'Location', 'bestoutside')
                currfig = gcf;
                title(currfig.Children(end),sprintf('Mean response of %s cells in the binocular condition',class{c}))
                 
                 plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat(sprintf('suppressedMean_mono_bino_orig_%s_cell',class{c})));
                 saveas(gcf,strcat(plotdir, '.png'));
                 %saveas(gcf,strcat(plotdir, '.svg'));
                
end

%% Figures for the paper: scatter plot of each peak comparing binocular vs monocular


newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\';

trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials']);
filenames = fieldnames(trialsTraces.peak_aligned_trials);

%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\results_simple_Anova_diffpks.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\interaction_results_lmer_allpks.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\monobino_maineffect_results_lmer_allpks.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\mixedmodel_results_anova_allpks.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\interaction_contrast_aov.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\interaction_contrast_indepSampleTtest.csv',',',1,1)';

pvalues = dlmread('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\mixedmodel_pvals_anova_linearTrend.csv',',',1,1)';
r2 = dlmread('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\mixedmodel_r2_anova_linearTrend.csv',',',1,1)';


 class ={'M','P','K'};
%ylims = [[0 160];[0 105];[0 210]];
ylims = [[0 250];[0 200];[0 240]];
% bsylims = [[-30 250];[-20 210];[-20 230]];
 %zylim = [-1 3];
bins =[1,6];
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) =[194/255 165/255 207/255] ; %--purple
col(5,:) = [253/255 174/255 97/255]; % -- orange
col(6,:) = [166/255 219/255 160/255]; % -- green
col(7,:) = [238/255 58/255 104/255]; % -- pink




 origTrace = nan(4,length(filenames),2);
 for i =1:length(filenames)
     filename = filenames{i};
     %if strcmp(trialsTraces.peak_aligned_trials.(filename).cellclass, class{c})
     for bin = 1:length(bins)
         binNb = sprintf('bin%d',bins(bin));

         for pn =1:4
             pkNb = sprintf('pk%d', pn);

             if nnz(strcmp(fieldnames(trialsTraces.peak_aligned_trials.(filename).origin),binNb))
                 origTrace(pn,i,bin) = max(squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).origin.(binNb).(pkNb),2)));

             end
         end
     end
 end

 %% Main figure
 figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);
 
 for p =1:4
     
     
     h= subplot(1,4,p);
     clear sidx nsidx idx
     scnt =0;
     cnt =0;
     for i=1:length(pvalues)
         x =ones(1,1);
         filename = filenames{i};
         
         
         if ~isempty(trialsTraces.peak_aligned_trials.(filename).cellclass)
             cnt = cnt+1;
             idx(cnt) =i;
             % [~, maxr2] = max(r2(:,i)); %take index of best explained variance through linear trend =1, quad trend =2, cubic =3
             % if pvalues(2,i) < 0.05 && (maxr2 ==1)
             if pvalues(6,i) < 0.05
                 
                 scnt = scnt+1;
                 plot([x,2*x],[origTrace(p,i,1),origTrace(p,i,2)],'-*','LineWidth', 3, 'Color',col(scnt+1,:))
                 sidx(scnt) =i;
                 cellclass{scnt} =trialsTraces.peak_aligned_trials.(filename).cellclass;
                 
             else
                 
                 plot([x,2*x],[origTrace(p,i,1),origTrace(p,i,2)],'-*','LineWidth', 2, 'Color',[161/255,161/255,161/255])
                 
                 
             end
             
             hold on
             
             if i == length(pvalues)
                 nsidx = 1:length(origTrace(p,:,1));
                 nsidx(sidx) = 0;
                 nsidx = find(nsidx);
                 h1 =  plot([x,2*x],[mean(origTrace(p,sidx,1),2),mean(origTrace(p,sidx,2),2)],'--*','LineWidth', 6, 'Color',col(7,:));
                 hold on
                 % h2 =  plot([x,2*x],[mean(origTrace(p,nsidx,1),2),mean(origTrace(p,nsidx,2),2)],'--*','LineWidth', 6, 'Color',col(1,:));
                 h2 =  plot([x,2*x],[mean(origTrace(p,idx,1),2),mean(origTrace(p,idx,2),2)],'--*','LineWidth', 6, 'Color',col(1,:));
                 %r2sig = mean(r2,2);
                 
             end
             
         end
     end
     plot([x,2*x],[mean(origTrace(1,sidx,1),2),mean(origTrace(1,sidx,2),2)],'>','LineWidth', 6, 'Color',[61/255,61/255,61/255])
     
     set(h,'position',get(h,'position').*[1 1 1.15 1])
     ylim(ylims(1,:))
     if p == 1
         
         ylabel('Spike rate (spikes/s)')
     end
     xlim([0 3])
     xticklabels({'','', 'Monocular','','Binocular',''})
     set(gca, 'linewidth',2)
     set(gca,'box','off')
     
     if p>1
         ax1 = gca;
         ax1.YAxis.Visible = 'off';
     end
 end
 %legend([h1 h2], 'Significant mean','Population mean', 'Location', 'bestoutside')
 currfig = gcf;
 %title(currfig.Children(end),sprintf('Peak responses of %s cells in the binocular condition',class{c}))
 title(currfig.Children(end),'Peak responses of all cells in the binocular condition')
 plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat('all_pks_mono_bino_orig_identifiedcells'));
 saveas(gcf,strcat(plotdir, '.png'));
 saveas(gcf,strcat(plotdir, '.svg'));
               
%% Figure inset
%getting indices for cells that were significant within the identified
%pool of cells
cnt =0;
scnt =0;
for i=1:length(pvalues)
    filename = filenames{i};
    if ~isempty(trialsTraces.peak_aligned_trials.(filename).cellclass)
        cnt = cnt+1;
        %index of cells that were identified
        idx(cnt) =i;
        
        if pvalues(6,i) < 0.05
            
            scnt = scnt+1;
            sidx(scnt) =i;
        end
    end
end

%percent change from Pk1
%individual cells
percentTrace = nan(size(origTrace));
for n =1:length(origTrace(1,:,1))
    percentTrace(:,n,1) = 100.*(origTrace(:,n,1)-origTrace(1,n,1))./origTrace(1,n,1);
    percentTrace(:,n,2) = 100.*(origTrace(:,n,2)-origTrace(1,n,2))./origTrace(1,n,2);
end

h =figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);
x =1:4;
%h1 =  plot(x,mean(percentTrace(:,sidx,1),2),'-','LineWidth', 6, 'Color',col(7,:));
%hold on
%h2 =  plot(x,mean(percentTrace(:,sidx,2),2),'--','LineWidth', 6, 'Color',col(7,:));
%hold on
h3 =  plot(x,mean(percentTrace(:,:,1),2),'-','LineWidth', 6, 'Color',col(1,:));
hold on
h4 =  plot(x,mean(percentTrace(:,:,2),2),'-','LineWidth', 6, 'Color',col(1,:));
h3.Color(4) = 0.4;


set(h,'position',get(h,'position').*[1 1 1.15 1])
xlim([0.5 4.5])
ylim([-10 1])
hold on
hline([0 0], '-k')
ylabel('Percent Change (%)')

xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
set(gca, 'linewidth',2)
set(gca,'box','off')

legend([h3 h4],'Monocular Population mean','Binocular Population mean', 'Location', 'bestoutside')
currfig = gcf;
title(currfig.Children(end),'Monocular vs binocular condition percent change from Pk1')

plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat('all_pks_mono_bino_percentPk1_allcells_indivmeans'));
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%same type of figure but computing the percent change with the overall units mean
%percent change from Pk1
meansigTrace =squeeze(mean(origTrace(:,sidx,:),2));
meanTrace =squeeze(mean(origTrace(:,:,:),2));
percentsigTrace = 100.*(meansigTrace-meansigTrace(1,:))./meansigTrace(1,:);
percentTrace = 100.*(meanTrace-meanTrace(1,:))./meanTrace(1,:);


h =figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);
x =1:4;

%h1 =  plot(x,percentsigTrace(:,1),'-','LineWidth', 6, 'Color',col(7,:));
%hold on
%h2 =  plot(x,percentsigTrace(:,2),'-','LineWidth', 6, 'Color',col(7,:));
%hold on
h3 =  plot(x,percentTrace(:,2),'-','LineWidth', 6, 'Color',col(1,:));
hold on
h4 =  plot(x,percentTrace(:,1),'-','LineWidth', 6, 'Color',col(1,:));
%h1.Color(4) = 0.4;
h4.Color(4) = 0.4;
 

set(h,'position',get(h,'position').*[1 1 1.15 1])
xlim([0.5 4.5])
ylim([-10 1])
hold on
hline([0 0], '-k')
ylabel('Percent Change (%)')

xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
set(gca, 'linewidth',2)
set(gca,'box','off')

legend([h3 h4],'Binocular Population mean','Monocular Population mean', 'Location', 'bestoutside')
currfig = gcf;
title(currfig.Children(end),'Monocular vs binocular condition percent change from Pk1')

plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat('all_pks_mono_bino_percentPk1_allcells_population_mean'));
%saveas(gcf,strcat(plotdir, '.png'));
%saveas(gcf,strcat(plotdir, '.svg'));

%plot normalized responses relative to Pk1 (very similar)
meansigTrace =squeeze(mean(origTrace(:,sidx,:),2));
meanTrace =squeeze(mean(origTrace(:,:,:),2));
normsigTrace = 100.*(meansigTrace)./meansigTrace(1,:);
normTrace = 100.*(meanTrace)./meanTrace(1,:);


h =figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);
x =1:4;
h1 =  plot(x,normsigTrace(:,2),'-','LineWidth', 6, 'Color',col(7,:));
hold on
h2 =  plot(x,normsigTrace(:,1),'-','LineWidth', 6, 'Color',col(7,:));
hold on
h3 =  plot(x,normTrace(:,2),'-','LineWidth', 6, 'Color',col(1,:));
hold on
h4 =  plot(x,normTrace(:,1),'-','LineWidth', 6, 'Color',col(1,:));
h2.Color(4) = 0.4;
h4.Color(4) = 0.4;

set(h,'position',get(h,'position').*[1 1 1.15 1])
xlim([0.5 4.5])
ylim([80 105])
hold on
hline([1 1], '-k')
ylabel('Percent of Pk1 (%)')

xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
set(gca, 'linewidth',2)
set(gca,'box','off')

legend([h1 h2 h3 h4],'Binocular Significant mean','Monocular Significant mean','Binocular Population mean','Monocular Population mean', 'Location', 'bestoutside')
currfig = gcf;
title(currfig.Children(end),'Monocular vs binocular condition normalized response relative to Pk1')

plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat('all_pks_mono_bino_percentPk1_allcells_population_mean'));
%saveas(gcf,strcat(plotdir, '.png'));
%saveas(gcf,strcat(plotdir, '.svg'));

%% Plotting the linear regressions of monocular vs binocular of significant units

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials']);
filenames = fieldnames(trialsTraces.peak_aligned_trials);
pvalues = dlmread('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\mixedmodel_pvals_anova_linearTrend.csv',',',1,1)';
bins =[1,6];
col(1,:) =[86/255 86/255 86/255] ; %--dark grey
col(2,:) = [251/255 154/255 153/255]; % -- pink
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) =[194/255 165/255 207/255] ; %--purple
col(5,:) = [253/255 174/255 97/255]; % -- orange
col(6,:) = [166/255 219/255 160/255]; % -- green
col(7,:) = [238/255 58/255 104/255]; % -- red
origTrace = nan(4,length(filenames),2);
normTrace = nan(4,length(filenames),2);
for i =1:length(filenames)
    filename = filenames{i};
    %if strcmp(trialsTraces.peak_aligned_trials.(filename).cellclass, class{c})
    for bin = 1:length(bins)
        binNb = sprintf('bin%d',bins(bin));
        
        for pn =1:4
            pkNb = sprintf('pk%d', pn);
            
            if nnz(strcmp(fieldnames(trialsTraces.peak_aligned_trials.(filename).origin),binNb))
                origTrace(pn,i,bin) = max(squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).origin.(binNb).(pkNb),2)));
                normTrace(pn,i,bin) = origTrace(pn,i,bin)./origTrace(1,i,bin);
            end
        end
        
    end
end

%getting indices for cells that were significant within the identified
%pool of cells
ncnt =0;
scnt =0;
clear filename
for i=1:length(pvalues)
    filename = filenames{i};
    if ~isempty(trialsTraces.peak_aligned_trials.(filename).cellclass)
       
        
        if pvalues(6,i) < 0.05 
            
            scnt = scnt+1;
            sidx(scnt) =i; %index of cells with significant interaction
        else 
            if pvalues(6,i) > 0.05
         ncnt = ncnt+1;
        %index of cells that were identified
        idx(ncnt) =i;
            end
        end
    end
end

%percent change from Pk1
%individual cells
RegCoeffs = nan([2,size(origTrace,2,3)]);
t =linspace(1,4,4);
for n =1:length(origTrace(1,:,1))
    RegCoeffs(:,n,1) = polyfit(t,normTrace(:,n,1)',1);
    RegCoeffs(:,n,2) =  polyfit(t,normTrace(:,n,2)',1);
end

%Plot slopes and intercept
h =figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);
x =1:4;
h1 =subplot(1,2,1);
y = RegCoeffs(2,:,1)' +repmat((t-1),46,1).*(RegCoeffs(1,:,1)');
plot(x,y,'-','LineWidth', 4, 'Color',[180/255 180/255 180/255])
hold on

y = mean(RegCoeffs(2,sidx,1)' +repmat((t-1),length(sidx),1).*(RegCoeffs(1,sidx,1)'),1);
plot(x,y,'-','LineWidth', 8, 'Color',col(7,:))
hold on
y = mean(RegCoeffs(2,:,1)' +repmat((t-1),length(1:46),1).*(RegCoeffs(1,:,1)'),1);
plot(x,y,'-','LineWidth', 8, 'Color','k')

xlim([0.5 4.5])
ylim([0.6 1.4])
xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
title('Monocular condition')
set(gca, 'linewidth',2)
set(gca,'box','off')

h2 =subplot(1,2,2);
y = RegCoeffs(2,:,2)'+repmat((t-1),46,1).*(RegCoeffs(1,:,2)');
plot(x,y,'-','LineWidth', 4, 'Color',col(1,:))
hold on

y = mean(RegCoeffs(2,sidx,2)' +repmat((t-1),length(sidx),1).*(RegCoeffs(1,sidx,2)'),1);
plot(x,y,'-','LineWidth', 8, 'Color',col(7,:))
hold on
y = mean(RegCoeffs(2,:,2)' +repmat((t-1),length(1:46),1).*(RegCoeffs(1,:,2)'),1);
plot(x,y,'-','LineWidth', 8, 'Color','k')

ylim([0.6 1.4])
xlim([0.5 4.5])
xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
set(h,'position',get(h,'position').*[1 1 1.15 1])
ylabel('Spike Rate (normalized)')
title('Binocular condition')
set(gca, 'linewidth',2)
set(gca,'box','off')

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat('all_pks_mono_bino_linearReg_allcells'));
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

% same plot only with means
h =figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);
x =1:4;
h1 =subplot(1,2,1);

y = mean(RegCoeffs(2,sidx,1)' +repmat((t-1),length(sidx),1).*(RegCoeffs(1,sidx,1)'),1);
plot(x,y,'-','LineWidth', 8, 'Color',col(7,:))
hold on
y = mean(RegCoeffs(2,:,1)' +repmat((t-1),length(1:46),1).*(RegCoeffs(1,:,1)'),1);
plot(x,y,'-','LineWidth', 8, 'Color','k')

xlim([0.5 4.5])
ylim([0.94 1.08])
xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
title('Monocular condition')
set(gca, 'linewidth',2)
set(gca,'box','off')

h2 =subplot(1,2,2);


y = mean(RegCoeffs(2,sidx,2)' +repmat((t-1),length(sidx),1).*(RegCoeffs(1,sidx,2)'),1);
plot(x,y,'-','LineWidth', 8, 'Color',col(7,:))
hold on
y = mean(RegCoeffs(2,:,2)' +repmat((t-1),length(1:46),1).*(RegCoeffs(1,:,2)'),1);
plot(x,y,'-','LineWidth', 8, 'Color','k')

ylim([0.94 1.08])
xlim([0.5 4.5])
xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
set(h,'position',get(h,'position').*[1 1 1.15 1])
ylabel('Spike Rate (normalized)')
title('Binocular condition')
set(gca, 'linewidth',2)
set(gca,'box','off')

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat('all_pks_mono_bino_linearReg_means'));
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));
%% Plot slopes separately from intercepts
h =figure('Renderer', 'painters', 'Position', [10 10 1200 1200]);
x =1:4;
%first, plot linear change over time
h1 =subplot(3,2,1);
y = repmat((t-1),46,1).*(RegCoeffs(1,:,1)');
plot(x,y,'-','LineWidth', 2, 'Color',[180/255 180/255 180/255])
hold on

y = mean(repmat((t-1),length(sidx),1).*(RegCoeffs(1,sidx,1)'),1);
plot(x,y,'-','LineWidth', 4, 'Color',col(7,:))
hold on
y =mean(repmat((t-1),length(1:46),1).*(RegCoeffs(1,:,1)'),1);
plot(x,y,'-','LineWidth', 4, 'Color','k')

xlim([0.5 4.5])
%ylim([0.6 1.4])
xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
ylabel('Spike Rate (normalized)')
title('Linear change in the monocular condition (slope)')
set(gca, 'linewidth',2)
set(gca,'box','off')

h2 =subplot(3,2,2);
y = repmat((t-1),46,1).*(RegCoeffs(1,:,2)');
plot(x,y,'-','LineWidth', 2, 'Color',col(1,:))
hold on

y = mean(repmat((t-1),length(sidx),1).*(RegCoeffs(1,sidx,2)'),1);
plot(x,y,'-','LineWidth', 4, 'Color',col(7,:))
hold on
y = mean(repmat((t-1),length(1:46),1).*(RegCoeffs(1,:,2)'),1);
plot(x,y,'-','LineWidth', 4, 'Color','k')

%ylim([0.6 1.4])
xlim([0.5 4.5])
xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
set(h,'position',get(h,'position').*[1 1 1.15 1])
title('Linear change in the binocular condition (slope)')
set(gca, 'linewidth',2)
set(gca,'box','off')

%plot intercepts
h3 =subplot(3,2,3:4);
boxplot([RegCoeffs(2,:,1)' RegCoeffs(2,:,2)'],'notch','off','labels',{'Monocular','Binocular'});
hx = findobj(gca,'Tag','Box');
patch(get(hx(1),'XData'),get(hx(1),'YData'),[80/255 80/255 80/255],'FaceAlpha',.5);
patch(get(hx(2),'XData'),get(hx(2),'YData'),[180/255 180/255 180/255],'FaceAlpha',.5);


hold on    
x2=ones(length(RegCoeffs(2,:,1)),1).*(1+(rand(length(RegCoeffs(2,:,1)),1)-0.5)/5);
x3=ones(length(RegCoeffs(2,:,1)),1).*(1+(rand(length(RegCoeffs(2,:,1)),1)-0.5)/10);
f1=scatter(x2,RegCoeffs(2,:,1),'k','filled','jitter', 'off');f1.MarkerFaceAlpha = 0.4;
hold on
f2=scatter(x3*2,RegCoeffs(2,:,2),'k', 'filled','jitter', 'off');%f2.MarkerFaceAlpha = 0.4;
set(gca, 'linewidth',2)
set(gca,'box','off')
ylabel('Initial Spike Rate(normalized)')
title('Regression intercepts of normalized data')

%Plot everything (intercept+slope) together over time
h4 =subplot(3,2,5);
y = RegCoeffs(2,:,1)' +repmat((t-1),46,1).*(RegCoeffs(1,:,1)');
plot(x,y,'-','LineWidth', 2, 'Color',[180/255 180/255 180/255])
hold on

y = mean(RegCoeffs(2,sidx,1)' +repmat((t-1),length(sidx),1).*(RegCoeffs(1,sidx,1)'),1);
plot(x,y,'-','LineWidth', 4, 'Color',col(7,:))
hold on
y = mean(RegCoeffs(2,:,1)' +repmat((t-1),length(1:46),1).*(RegCoeffs(1,:,1)'),1);
plot(x,y,'-','LineWidth', 4, 'Color','k')

xlim([0.5 4.5])
ylim([0.6 1.4])
xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
title('Linear fit: Monocular condition')
set(gca, 'linewidth',2)
set(gca,'box','off')

h5 =subplot(3,2,6);
y = RegCoeffs(2,:,2)'+repmat((t-1),46,1).*(RegCoeffs(1,:,2)');
plot(x,y,'-','LineWidth', 2, 'Color',col(1,:))
hold on

y = mean(RegCoeffs(2,sidx,2)' +repmat((t-1),length(sidx),1).*(RegCoeffs(1,sidx,2)'),1);
plot(x,y,'-','LineWidth', 4, 'Color',col(7,:))
hold on
y = mean(RegCoeffs(2,:,2)' +repmat((t-1),length(1:46),1).*(RegCoeffs(1,:,2)'),1);
plot(x,y,'-','LineWidth', 4, 'Color','k')

ylim([0.6 1.4])
xlim([0.5 4.5])
xticklabels({'','Pk1','','Pk2','', 'Pk3','','Pk4',''})
set(h,'position',get(h,'position').*[1 1 1.15 1])
ylabel('Spike Rate (normalized)')
title('Linear fit: Binocular condition')
set(gca, 'linewidth',2)
set(gca,'box','off')

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat('all_pks_mono_bino_linearReg_allcells_slope_intercept'));
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% Plotting the slopes of monocular vs binocular



newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';

trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials']);
filenames = fieldnames(trialsTraces.peak_aligned_trials);

%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\results_simple_Anova_diffpks.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\interaction_results_lmer_allpks.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\monobino_maineffect_results_lmer_allpks.csv',',',1,1)';
pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\mixedmodel_results_anova_allpks.csv',',',1,1)';


 %class ={'M','P','K'};
%ylims = [[0 160];[0 105];[0 210]];
ylims = [[0 250];[0 200];[0 240]];
% bsylims = [[-30 250];[-20 210];[-20 230]];
 %zylim = [-1 3];
bins =[1,6];
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) = [238/255 58/255 104/255]; % -- pink

   
        % origTrace = nan(4,length(filenames),2);
         zTrace = nan(4,length(filenames),2);
         zslope = nan(length(filenames),2);
          t =linspace(1,1100,4);
   for i =1:length(filenames)
              filename = filenames{i};
       
                for bin = 1:length(bins)
                    binNb = sprintf('bin%d',bins(bin));
                    
                        
                    for pn =1:4
                         pkNb = sprintf('pk%d', pn);
                        
                
                        if nnz(strcmp(fieldnames(trialsTraces.peak_aligned_trials.(filename).origin),binNb))
                             %origTrace(pn,i,bin) = max(squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).origin.(binNb).(pkNb),2))); 
                             zTrace(pn,i,bin) =  max(squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).zscore.(binNb).(pkNb),2))); 

                            
                        end
                    end
                    slope = polyfit(t,zTrace(:,i,bin)',1);
                    zslope(i,bin) = slope(1);
                end
         
   end
  
   
 
   
                h =figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);

               
                    
                    
                     clear idx nsidx
                     scnt =0;
                     for i=1:length(pvalues)
                      x =ones(1,1);
                         if pvalues(3,i) > 0.05
                     
                       plot([x,2*x],[zslope(i,1),zslope(i,2)],'-*','LineWidth', 2, 'Color',[161/255,161/255,161/255])
                  
                         else
                       plot([x,2*x],[zslope(i,1),zslope(i,2)],'-*','LineWidth', 3, 'Color',col(4,:))
                       scnt = scnt+1;
                       sidx(scnt) =i;
                       
                         end
                       
                          hold on
                         
                         if i == length(pvalues)
                       nsidx = 1:length(origTrace(1,:,1));
                       nsidx(sidx) = 0;
                       nsidx = find(nsidx);
                        h1 =    plot([x,2*x],[mean(zslope(sidx,1)),mean(zslope(sidx,2))],'--*','LineWidth', 6, 'Color',col(4,:));
                            hold on
                        h2 =    plot([x,2*x],[mean(zslope(nsidx,1)),mean(zslope(nsidx,2))],'--*','LineWidth', 6, 'Color',col(3,:));
                   
         
                         end
                       
                             
                        
                       
                     end
                         %monocular
                      %scatter(origTrace(p,:,1),'o', 'Color',col(c,:))
                
                    set(h,'position',get(h,'position').*[1 1 1.15 1])
                    %ylim(ylims(1,:))
 
                        ylabel('Spike rate change across peaks (z-score slope)')
                   
                    xlim([0 3])
                    xticklabels({'','', 'Monocular','','Binocular',''})
                    set(gca, 'linewidth',2)
                    set(gca,'box','off')
                           
                legend([h1 h2],'Significant Interaction Mean','NS', 'Location', 'bestoutside')
                currfig = gcf;
                %title(currfig.Children(end),sprintf('Peak responses of %s cells in the binocular condition',class{c}))
                title(currfig.Children(end),'Slope across peaks of all cells in the monocular and binocular condition')
         
                 plotdir = strcat('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\plots\',strcat('all_pks_mono_bino_orig_allcells'));
                %saveas(gcf,strcat(plotdir, '.png'));
                 %saveas(gcf,strcat(plotdir, '.svg'));
                
%end

%% scatter plot of peak 1 vs Pk4 comparing binocular vs monocular


newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';

trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials']);
filenames = fieldnames(trialsTraces.peak_aligned_trials);

%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\results_simple_Anova_diffpks.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\interaction_results_lmer_allpks.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\monobino_maineffect_results_lmer_allpks.csv',',',1,1)';
%pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\mixedmodel_results_anova_allpks.csv',',',1,1)';
pvalues = dlmread('C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\interaction_contrast_aov.csv',',',1,1)';


 class ={'M','P','K'};
%ylims = [[0 160];[0 105];[0 210]];
ylims = [[0 250];[0 200];[0 240]];
% bsylims = [[-30 250];[-20 210];[-20 230]];
 %zylim = [-1 3];
bins =[1,6];
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) = [238/255 58/255 104/255]; % -- pink

%for c = 1:length(class)
    
         zTrace = nan(4,length(filenames),2);
    for i =1:length(filenames)
             filename = filenames{i};
        
                for bin = 1:length(bins)
                    binNb = sprintf('bin%d',bins(bin));
                    
                        
                    for pn =1:4
                         pkNb = sprintf('pk%d', pn);
                        
                
                        if nnz(strcmp(fieldnames(trialsTraces.peak_aligned_trials.(filename).zscore),binNb))
                               zTrace(pn,i,bin) =  max(squeeze(nanmean(trialsTraces.peak_aligned_trials.(filename).zscore.(binNb).(pkNb),2))); 

                        end
                    end
                end
    end
     
   
               h= figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);

                                         
                     clear idx nsidx
                     scnt =0;
                     for i=1:length(pvalues)
                      x =ones(1,1);
                         if pvalues(1,i) > 0.05
                     
                       %plot([x,2*x],[origTrace(1,i,1),origTrace(4,i,1)],'-*','LineWidth', 2, 'Color',[161/255,161/255,161/255])
                       %hold on
                       %plot([x,2*x],[origTrace(1,i,2),origTrace(4,i,2)],'-*','LineWidth', 2, 'Color',[161/255,161/255,161/255])
                  
                         else
                 %      plot([x,2*x],[origTrace(1,i,1),origTrace(4,i,1)],'-*','LineWidth', 3, 'Color',[161/255,161/255,161/255])
                  %     hold on
                   %     plot([x,2*x],[origTrace(1,i,2),origTrace(4,i,2)],'-*','LineWidth', 3, 'Color',[161/255,161/255,161/255])

                       scnt = scnt+1;
                       sidx(scnt) =i;
                       
                         end
                       
                          hold on
                         
                         if i == length(pvalues)
                       nsidx = 1:length(zTrace(1,:,1));
                       nsidx(sidx) = 0;
                       nsidx = find(nsidx);
                    plot([x,2*x],[mean(zTrace(1,sidx,1),2),mean(zTrace(4,sidx,1),2)],'--*','LineWidth', 6, 'Color',col(4,:))
                    hold on
                    plot([x,2*x],[mean(zTrace(1,nsidx,1),2),mean(zTrace(4,nsidx,1),2)],'--*','LineWidth', 6, 'Color',[141/255,141/255,141/255])
                    hold on
                    plot([x,2*x],[mean(zTrace(1,sidx,2),2),mean(zTrace(4,sidx,2),2)],'--*','LineWidth', 6, 'Color',col(4,:))
                    hold on
                    plot([x,2*x],[mean(zTrace(1,nsidx,2),2),mean(zTrace(4,nsidx,2),2)],'--*','LineWidth', 6, 'Color',[141/255,141/255,141/255])
                    hold on
                    plot([x,2*x],[median(zTrace(1,sidx,1),2),median(zTrace(4,sidx,1),2)],'--*','LineWidth', 6, 'Color',col(2,:))
                    hold on
                    plot([x,2*x],[median(zTrace(1,nsidx,1),2),median(zTrace(4,nsidx,1),2)],'--*','LineWidth', 6, 'Color',[161/255,161/255,161/255])
                    hold on
                    plot([x,2*x],[median(zTrace(1,sidx,2),2),median(zTrace(4,sidx,2),2)],'--*','LineWidth', 6, 'Color',col(2,:))
                    hold on
                    plot([x,2*x],[median(zTrace(1,nsidx,2),2),median(zTrace(4,nsidx,2),2)],'--*','LineWidth', 6, 'Color',[161/255,161/255,161/255])

         
                         end
                       
                             
                        
                       
                     end
              
                set(h,'position',get(h,'position').*[1 1 1.15 1])
                %ylim(ylims(1,:))
                ylabel('Spike rate (z-score)')
                xlim([0 3])
                xticklabels({'','', 'Pk1','','Pk4',''})
                set(gca, 'linewidth',2)
                set(gca,'box','off')

               
                legend('Mean Mono Significant Pk*Condition','Mean Mono NS','Mean Bino Significant Pk*Condition','Mean Bino NS',...
                    'Median Mono Significant Pk*Condition','Median Mono NS','Median Bino Significant Pk*Condition','Median Bino NS', 'Location', 'bestoutside')
                currfig = gcf;
                %title(currfig.Children(end),sprintf('Peak responses of %s cells in the binocular condition',class{c}))
                title(currfig.Children(end),'Peak 1 and Peak 4 responses of all cells in the binocular vs monocular condition')
         
   
