function [wind_peak_vals, peak_aligned_trials] = binPeakTrigResps(peakLocs,binSpkTrials, NoFiltMultiContSUA)
%This script is inspired from "get_origin_peaks.m". It isolates windows of binary spikes around peak
%locations. The script also retriggers binary trial responses to each serial cycle peak location
%Code developped and last edited on 05/11/2021
xabs = -199:1300;
nyq = 500;

channum = 1: length(fieldnames(binSpkTrials));
%mean_filtered_dSUA = struct();
bin_suas_aligned_trials = struct();
peak_aligned_trials = struct();
wind_peak_vals = struct();

%mean_peaks =struct();
%median_peaks = struct();
max_low_dist = struct();
%filenames = cell(length(channum),2);
contLims = [0,0.1,0.3,0.5,0.7,1];
filenames = fieldnames(binSpkTrials);
for i = channum
    for n = 1:length(contLims)
        
        if n ==1 || n ==6 %if we only wanna save the monocular and binocular condition
            binNb = sprintf('bin%d', n);
            filename = filenames{i};
            if isfield(binSpkTrials.(filename), binNb)
                if ~isempty(binSpkTrials.(filename).(binNb))

                    trialidx = 1:length(binSpkTrials.(filename).(binNb)(1,:));
                    bin_dSUA = binSpkTrials.(filename).(binNb)(401:1900,:); %binary data
                    origin_dSUA = NoFiltMultiContSUA.(filename).(binNb)(401:1900,:);%we also need the origin data to detect spike location on the unsmoothed data and get a window of binary data around unsmoothed spike locations
                    %filtered_dMUA = FiltMultiContMUA.(filename).(binNb);

                    %obtain the peak location of interest for each trial of a given single
                    %unit
                    all_locsdSUA_trials = peakLocs.(filename).(binNb);

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
                    bin_fp_locked_trials = nan(new_dist_unit,length(bin_dSUA(1,:)),4); 
                    fp_locked_trials = nan(new_dist_unit,length(origin_dSUA(1,:)),4);

                    clear n pn
                    for pn =1:4
                        for n =trialidx
                            lower_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+1;
                            upper_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+length(xabs);

                            %origin data of the statistical analysis
                            bin_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = bin_dSUA(:,n);
                            %binary data
                            fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = origin_dSUA(:,n);
                        end

                    end
                    %get the aligned data if it exists for the unit
                    bin_suas_aligned_trials.(filename).binarySUA.(binNb)= bin_fp_locked_trials;
                   

                    max_low_dist.(filename).(binNb) = max_low_dist_unit;

                    clear pn
                    for pn = 1:4
                        %peak data for the stats
                        [~,pnloc]= max(fp_locked_trials(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn), [],1);
                        wind_peak_vals.(filename).(binNb)(pn,:,:)= bin_fp_locked_trials(max_low_dist.(filename).(binNb)+pnloc-1-30:max_low_dist.(filename).(binNb)+pnloc-1+30,:,pn);
                        

                        %peak triggered trials
                        pkNb = sprintf('pk%d', pn);
                        peak_aligned_trials.(filename).binarySUA.(binNb).(pkNb) = bin_suas_aligned_trials.(filename).binarySUA.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn);

                    end

                    %mean peaks for the R plots
                    %mean_peaks.(filename).(binNb) = mean(wind_peak_vals.(filename).(binNb),2);

                    %median
                    %median_peaks.(filename).(binNb) = median(wind_peak_vals.(filename).(binNb),2);


                    %zscore
                    %zscore.(filename).(binNb) = (peak_vals.(filename).(binNb)-repmat(mean(peak_vals.(filename).(binNb),2), 1,length(peak_vals.(filename).(binNb)(1,:))))./repmat(std(peak_vals.(filename).(binNb),0,2),1,length(peak_vals.(filename).(binNb)(1,:)));


                end
            end
        end
    end
    
end