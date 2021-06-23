function [peak_aligned_trials] = suaPeakTrigResps(NoFiltMultiContSUA)
%This script is inspired from "get_origin_peaks.m". It retriggers convolved trial responses to each serial cycle peak location
%Code developped and last edited on 06/22/2021
%NoFiltMultiContSUA = trialsTraces.NoFiltMultiContSUA;
xabs = -199:1300;
nyq = 500;

channum = 1: length(fieldnames(NoFiltMultiContSUA));
peak_aligned_trials = struct();
max_low_dist = struct();
contLims = [0,0.1,0.3,0.5,0.7,1];
filenames = fieldnames(NoFiltMultiContSUA);
for i = channum
    for n = 1:length(contLims)
        
        if n ==1 || n ==6 %if we only wanna save the monocular and binocular condition
            binNb = sprintf('bin%d', n);
            filename = filenames{i};
            if isfield(NoFiltMultiContSUA.(filename),binNb) && isfield(NoFiltMultiContSUA.(filename).(binNb), 'neuralDat')
                if ~isempty(NoFiltMultiContSUA.(filename).(binNb).neuralDat)
                    
                    trialidx = 1:length(NoFiltMultiContSUA.(filename).(binNb).neuralDat(1,:));
                    origin_dSUA = NoFiltMultiContSUA.(filename).(binNb).neuralDat(401:1900,:);%we also need the origin data to detect spike location on the unsmoothed data and get a window of binary data around unsmoothed spike locations
                    
                    %obtain the peak location of interest for each trial of a given single
                    %unit
                    all_locsdSUA_trials = NoFiltMultiContSUA.(filename).(binNb).peaklocs;
                    
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
                    
                    clear n pn
                    for pn =1:4
                        for n =trialidx
                            lower_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+1;
                            upper_unit_bound =max_low_dist_unit-all_locsdSUA_trials(pn,n)+length(xabs);
                            
                            %origin data of the statistical analysis
                            fp_locked_trials(lower_unit_bound:upper_unit_bound,n,pn) = origin_dSUA(:,n);
                        end
                        
                    end
                    
                    max_low_dist.(filename).(binNb) = max_low_dist_unit;
                    
                    clear pn
                    for pn = 1:4
                        %peak data for the stats
                        [~,pnloc]= max(fp_locked_trials(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn), [],1);
                        %peak triggered trials
                        pkNb = sprintf('pk%d', pn);
                        %peak_aligned_trials.(filename).originSUA.(binNb).(pkNb)
                        %=
                        %fp_locked_trials(max_low_dist.(filename).(binNb)+(pnloc-125)-1-124:max_low_dist.(filename).(binNb)+(pnloc-125)-1+125,:,pn);
                        %%this one works better for binary data
                        peak_aligned_trials.(filename).originSUA.(binNb).(pkNb) = fp_locked_trials(max_low_dist.(filename).(binNb)-1-124:max_low_dist.(filename).(binNb)-1+125,:,pn); %this one works better for convolved data
                        
                    end
                    
                end
            end
        end
    end
end
