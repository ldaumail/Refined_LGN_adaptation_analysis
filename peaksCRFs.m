function [peak_vals] = peaksCRFs(peakLocs,NoFiltMultiContMUA)
%This script is almost identical to "get_origin_peaks.m" as it isolates peak
%values of discretized, convolved original MUA (not smoothed). last edited
%on 05/13/2021 by Loic Daumail
%here we only get responses to monocular stimulations

xabs = -199:1300;
nyq = 500;

channum = 1: length(fieldnames(NoFiltMultiContMUA));
%mean_filtered_dSUA = struct();
muas_aligned_trials = struct();
peak_aligned_trials = struct();
peak_vals = struct();
zscore_peak_vals = struct();
%mean_peak_vals = struct();
mean_peaks =struct();
median_peaks = struct();
max_low_dist = struct();
%filenames = cell(length(channum),2);
contLims = [0,0.1,0.3,0.5,0.7,1];
filenames = fieldnames(NoFiltMultiContMUA);
for i = channum
    for n = 1:length(contLims)
        
        binNb = sprintf('bin%d', n);
        filename = filenames{i};
        if isfield(NoFiltMultiContMUA.(filename), binNb)
            if ~isempty(NoFiltMultiContMUA.(filename).(binNb))
                
                trialidx = 1:length(NoFiltMultiContMUA.(filename).(binNb)(1,:));
                origin_dMUA = NoFiltMultiContMUA.(filename).(binNb)(401:1900,:); %- mean(NoFiltMultiContMUA.(filename).(binNb)(401:600,:),1);
                bs_origin_dMUA = NoFiltMultiContMUA.(filename).(binNb)(401:1900,:)- mean(NoFiltMultiContMUA.(filename).(binNb)(401:600,:),1);
                
                %origin_dMUA = NoFiltMultiContMUA.BsNoFiltMultiContMUA.(filename).(binNb)(1:end,:);
                %create zscore origin trials data to plot average peaks for each unit with R
                
                zscore_unit = nan(size(origin_dMUA));
                clear tr
                for tr =trialidx
                    
                    zscore_unit(:,tr) = (origin_dMUA(:,tr)-mean(origin_dMUA(:,tr)))./(std(origin_dMUA(:,tr)));
                end
                
                
                %filtered_dMUA = FiltMultiContMUA.(filename).(binNb);
                
                
                %determine the peak location of interest for each trial of a given single
                %unit
                all_locsdMUA_trials = peakLocs.(filename).(binNb);
                
                up_dist_trials = nan(4,length(trialidx));
                clear pn
                for pn = 1:4
                    locs_peak = all_locsdMUA_trials(pn, :);
                    up_dist_trials(pn,:)= length(xabs)- locs_peak;
                end
                %get the max distance between the peakalign and the stimulus onset
                max_low_dist_unit = max(all_locsdMUA_trials,[],'all');
                %create new matrix with the length(max(d)+max(xabs - d))
                new_dist_unit = max_low_dist_unit + max(up_dist_trials,[],'all');
                fp_locked_trials = nan(new_dist_unit,length(origin_dMUA(1,:)),4);
                bs_fp_locked_trials = nan(new_dist_unit,length(origin_dMUA(1,:)),4);
                zscore_fp_locked_trials = nan(new_dist_unit,length(origin_dMUA(1,:)),4);
                %filtered_fp_locked_trials = nan(new_dist_unit,length(filtered_dSUA(1,:)),4);
                clear n pn
                for pn =1:4
                    for n =trialidx
                        lower_unit_bound =max_low_dist_unit-all_locsdMUA_trials(pn,n)+1;
                        upper_unit_bound =max_low_dist_unit-all_locsdMUA_trials(pn,n)+length(xabs);
                        
                        %origin data of the statistical analysis
                        fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = origin_dMUA(:,n);
                        %baseline corrected data
                        bs_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = bs_origin_dMUA(:,n);
                        %zscore transformed data for the plotting
                        zscore_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = zscore_unit(:,n);
                        
                        %filtered_fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = filtered_dSUA(:,n);
                    end
                    
                end
                %get the aligned data if it exists for the unit
                muas_aligned_trials.(filename).origin.(binNb)= fp_locked_trials;
                muas_aligned_trials.(filename).bsorigin.(binNb)= bs_fp_locked_trials;
                muas_aligned_trials.(filename).zscore.(binNb)= zscore_fp_locked_trials;
                
                max_low_dist.(filename).(binNb) = max_low_dist_unit;
                
                clear pn
                for pn = 1:4
                    %peak data for the stats
                    peak_vals.(filename).(binNb)(pn,:)= max(muas_aligned_trials.(filename).origin.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn), [],1);
                    zscore_peak_vals.(filename).(binNb)(pn,:)= max(muas_aligned_trials.(filename).zscore.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn), [],1);
                    
                    
                    %peak triggered trials
                    pkNb = sprintf('pk%d', pn);
                    peak_aligned_trials.(filename).origin.(binNb).(pkNb) = muas_aligned_trials.(filename).origin.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn);
                    peak_aligned_trials.(filename).bsorigin.(binNb).(pkNb) = muas_aligned_trials.(filename).bsorigin.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn);
                    peak_aligned_trials.(filename).zscore.(binNb).(pkNb) = muas_aligned_trials.(filename).zscore.(binNb)(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn);
                    
                    
                    
                end
                
                %mean peaks for the R plots
                %mean_peaks.(filename).(binNb) = mean(peak_vals.(filename).(binNb),2);
                
                %median
                %median_peaks.(filename).(binNb) = median(peak_vals.(filename).(binNb),2);
                
                
                %zscore
                %zscore.(filename).(binNb) = (peak_vals.(filename).(binNb)-repmat(mean(peak_vals.(filename).(binNb),2), 1,length(peak_vals.(filename).(binNb)(1,:))))./repmat(std(peak_vals.(filename).(binNb),0,2),1,length(peak_vals.(filename).(binNb)(1,:)));
                
                
            end
        end
    end
    
    
end