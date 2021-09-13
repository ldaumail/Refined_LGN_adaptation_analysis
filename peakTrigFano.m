function [slide_win_fanof, peak_aligned_trials] = peakTrigFano(peakLocs,NoFiltMultiContSUA,binSpkTrials,wsz)
%This script is inspired from "get_origin_peaks.m". It isolates windows of binary spikes around peak
%locations. The script also retriggers binary trial responses to each serial cycle peak location
%Code developped and last edited on 05/11/2021
%INPUT:
%peakLocs = peak locations of the smoothed data
%NoFiltMultiContSUA = convolved single unit activity of all selected trials
%binSpkTrials = binary data of spikes of all selected trials
%wsz = window size in ms
xabs = -199:1300;
nyq = 500;

channum = 1: length(fieldnames(binSpkTrials));
%mean_filtered_dSUA = struct();
bin_suas_aligned_trials = struct();
suas_aligned_trials = struct();
peak_aligned_trials = struct();
%wind_peak_vals = struct();
slide_win_fanof = struct();

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
                    for w = 1:length(wsz)
                        trialidx = 1:length(binSpkTrials.(filename).(binNb)(1,:));
                        bin_dSUA = binSpkTrials.(filename).(binNb)(401:1900,:); %binary data
                        origin_dSUA = NoFiltMultiContSUA.(filename).(binNb)(401:1900,:);%we also need the origin data to detect spike location on the unsmoothed data and get a window of binary data around unsmoothed spike locations
                       
                        rs_bin = binSpkTrials.(filename).(binNb)(151:400,:); %get binary spiking data at resting state
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
                        suas_aligned_trials.(filename).convolvedSUA.(binNb) = fp_locked_trials;
                        
                        max_low_dist.(filename).(binNb) = max_low_dist_unit;
                        
                        clear pn
                        for pn = 1:4
                            %peak data for the stats
                            [~,pnloc]= max(fp_locked_trials(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn), [],1);
                            %wind_peak_vals.(filename).(binNb)(pn,:,:)= bin_fp_locked_trials(max_low_dist.(filename).(binNb)+pnloc-1-30:max_low_dist.(filename).(binNb)+pnloc-1+30,:,pn);
                            %figure(); plot(fp_locked_trials(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn));
                           % figure(); scatter(pnloc,max(fp_locked_trials(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn), [],1));
                            %Compute fano factors around peaks in a range of different sliding window sizes,
                            %with 10 ms steps around each peak
                            
                            fanof = nan(length(-125:10:125-wsz(w)),1);
                            meanspkc = nan(length(-125:10:125-wsz(w)),1);
                            varspkc = nan(length(-125:10:125-wsz(w)),1);
                            for t = 1:length(-125:10:125-wsz(w))
                               % wind_vals = bin_fp_locked_trials(max_low_dist.(filename).(binNb)+pnloc-125+10*(t-1):max_low_dist.(filename).(binNb)+pnloc-125+10*(t-1)+wsz(w),:,pn);
                               wind_vals = bin_fp_locked_trials(max_low_dist.(filename).(binNb)+(-125+pnloc)-125+10*(t-1):max_low_dist.(filename).(binNb)+(-125+pnloc)-125+10*(t-1)+wsz(w),:,pn);
                                %wind_vals = bin_fp_locked_trials(max_low_dist.(filename).(binNb)-125+10*(t-1):max_low_dist.(filename).(binNb)-125+10*(t-1)+wsz(w),:,pn);
                              
                                count = sum(wind_vals,1,'omitnan');
                                meanspkc(t) = nanmean(count);
                                varspkc(t) = var(count, 'omitnan');
                                fanof(t) = varspkc(t)/meanspkc(t);
                            end
                            %figure(); plot(-125:10:125-wsz(w),fanof);
                            %hold on
                           % plot(-250:125,bin_fp_locked_trials(max_low_dist.(filename).(binNb)+(-125+pnloc)-250:max_low_dist.(filename).(binNb)+(-125+pnloc)+125,:,pn))
                          %  figure();plot(-250:125,bin_fp_locked_trials(max_low_dist.(filename).(binNb)-250:max_low_dist.(filename).(binNb)+125,:,pn))
                           
                            windSz = sprintf('wsz%d',wsz(w));
                            slide_win_fanof.(filename).fanof.(binNb).(windSz).peaks(pn,:) = fanof;
                            slide_win_fanof.(filename).meanspkc.(binNb).(windSz).peaks(pn,:) = meanspkc;
                            slide_win_fanof.(filename).varspkc.(binNb).(windSz).peaks(pn,:) = varspkc;
                            
                            
                            %peak triggered trials
                            pkNb = sprintf('pk%d', pn);
                            peak_aligned_trials.(filename).binarySUA.(binNb).(pkNb) = bin_suas_aligned_trials.(filename).binarySUA.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn);
                            peak_aligned_trials.(filename).convolvedSUA.(binNb).(pkNb) = suas_aligned_trials.(filename).convolvedSUA.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn);
                            
                        end
                            %add fano factor at baseline activity
                            fanof = nan(length(-125:10:125-wsz(w)),1);
                            meanspkc = nan(length(-125:10:125-wsz(w)),1);
                            varspkc = nan(length(-125:10:125-wsz(w)),1);
                            for t = 1:length(-125:10:125-wsz(w))
                                wind_vals = rs_bin(1+10*(t-1):10*(t-1)+wsz(w),:);
                                count = sum(wind_vals,1);
                                meanspkc(t) = nanmean(count);
                                varspkc(t) = var(count, 'omitnan');
                                fanof(t) = varspkc(t)/meanspkc(t);
                            end
                            windSz = sprintf('wsz%d',wsz(w));
                            slide_win_fanof.(filename).fanof.(binNb).(windSz).rs(:) = fanof;
                            slide_win_fanof.(filename).meanspkc.(binNb).(windSz).rs(:) = meanspkc;
                            slide_win_fanof.(filename).varspkc.(binNb).(windSz).rs(:) = varspkc;
                            
                    end
                    
                end
            end
        end
    end
    
end