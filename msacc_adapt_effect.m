%In this script the goal is to determine whether a difference can be
%observed between trials selected without accounting for microsaccades, and
%trials selected by rejecting those that include microsaccades

indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\analysis\';
concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
%trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']); %neural data
trialsTraces =load([newdatadir 'NoFiltMultiContSUA_06212021']); %neural data +peaklocs + trial numbers obtained with  "BinocularAdaptationTrialSelection.m"
%get pvalues from lmer results with Dunnett correction
pvalues = dlmread('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\lmer_results_orig_03032020_corrected_dunnett.csv', ',', 1,1);
pk1pk4pval = pvalues(~isnan(pvalues(:,3)),3);%clear out nans
     
xfilenames = fieldnames(trialsTraces.NoFiltMultiContSUA);
cnt = 0;
eyeMovData = struct();


for i=1:length(xfilenames)
    if pk1pk4pval(i) < 0.05
        try
            xcluster = xfilenames{i};
            cluster = xcluster(2:end);
            underscore = strfind(cluster, '_');
            session =  cluster(1:underscore(2)-1);
            directory = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\',cluster,'\');
            
            xBRdatafiles = concat_filenames.(xcluster);
            eye_info =[];
            all_codes = [];
            all_times = [];
            all_analogData =[];
            for fn =1:length(xBRdatafiles)
                xBRdatafile = xBRdatafiles{fn};
                filename   = [directory xBRdatafile(2:end)];
                if exist(strcat(filename, '.bhv'),'file')
                    eye_info.(strcat(xBRdatafile,'_bhvfile')) = concatBHV(strcat(filename,'.bhv'));
                    all_codes = [all_codes, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeNumbers];
                    all_times = [all_times, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeTimes];
                    all_analogData = [all_analogData,eye_info.(strcat(xBRdatafile,'_bhvfile')).AnalogData];
                    xBaseline = eye_info.(strcat(xBRdatafile,'_bhvfile')).ScreenXresolution/4/eye_info.(strcat(xBRdatafile,'_bhvfile')).PixelsPerDegree; %since the screen monitor is split in 2 parts with the stereoscope, the center for each eye becomes the center for each side of the stereoscope (half of the half, justifying dividing by 4
                end
                
            end
            samplerate = 1000;
            %trialindex = condSelectedTrialsIdx.(xcluster);
            trialindex = trialsTraces.NoFiltMultiContSUA.(xcluster).bin1.trials; %only take trial indices of monocular stimulation
            
            ampl = [];
            veloc = [];
            ntr =0; %number of trials in which microsaccades were detected
            excludeTrials = zeros(length(trialindex),1);
            nprocTrials = zeros(length(trialindex),1);
            excluSaccs = nan(20,length(trialindex),1);
            sparedSaccs = nan(20,length(trialindex),1);
            
            for tr = 1:length(trialindex)
                if trialindex(tr) <= length(all_codes)
                    codes                 = all_codes{trialindex(tr)};
                    times                 = all_times{trialindex(tr)};
                    
                    if nnz(find( codes == 23))
                        samples = [];
                        samples(:,1) = (-1*times(codes == 23)+1) : 1 : 0 : (length(all_analogData{trialindex(tr)}.EyeSignal(:,1)) - times(codes == 23)); %trigger time points on stimulus onset time for it to be 0. Everything before that point is then negative %23 = stimulus onset time, time is measured relative to each trial onset recorded %24 = trial/stimulus offset time
                        if ~isempty(samples)
                            samples(:,2) = all_analogData{trialindex(tr)}.EyeSignal(:,1)+xBaseline; %horizontal position of the left eye in degrees baseline corrected
                            samples(:,3) = all_analogData{trialindex(tr)}.EyeSignal(:,2); %vertical position of the left eye in degrees
                            samples(:,4) = nan();
                            samples(:,5) = nan();
                            blinks = zeros(length(samples(:,1)),1);
                            recording = ClusterDetection.EyeMovRecording.Create(directory, session, samples, blinks, samplerate);
                            
                            % Runs the saccade detection
                            try
                                [saccades, stats] = recording.FindSaccades();
                                enum = ClusterDetection.SaccadeDetector.GetEnum;
                                
                                for s =1:length(saccades(:,enum.startIndex)) %loop through all microssaccades/ saccades found in one trial
                                    if  ~isempty(find(saccades(s,enum.startIndex),1)) && (saccades(s,enum.startIndex) > 0 && saccades(s,enum.startIndex) < times(codes == 24)- times(codes == 23)) %if there is at least 1 microsaccade and it occurs between stim onset and stim offset
                                        excludeTrials(tr) = excludeTrials(tr)+1;
                                        excluSaccs(s,tr) = saccades(s,enum.startIndex);
                                    else
                                        if  ~isempty(find(saccades(s,enum.startIndex),1)) && (saccades(s,enum.startIndex) < 0 || saccades(s,enum.startIndex) > times(codes == 24)- times(codes == 23))
                                            sparedSaccs(s,tr) = saccades(s,enum.startIndex);
                                        end
                                    end
                                end
                            catch
                                nprocTrials(tr) =nprocTrials(tr)+1;
                                disp(strcat({'Bad Trial',xBRdatafile}))
                            end
                            
                            % Plots a main sequence
                            ampl = [ampl; saccades(:,enum.amplitude)];
                            veloc = [veloc; saccades(:,enum.peakVelocity)];
                            eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).stats = stats;
                            eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).saccades = saccades;
                            eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).enum = enum;
                            eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).samples = samples;
                            if ~all(isnan((saccades(:,enum.startIndex))))
                                ntr =ntr+1;
                            end
                            %}
                        end
                    end
                end
            end
            %eyeMovData.(xcluster).cellclass = trialsTraces.NoFiltMultiContSUA.(xcluster).cellclass;
            %end
            %all traces without ones unprocessed by algorithm
            allTrialsTraces =trialsTraces.NoFiltMultiContSUA.(xcluster).bin1.neuralDat(:,~nprocTrials);
            %selected trials neural data excluding unprocessed ones and ones
            %with msaccs during stimulation
            selectedTrialsTraces =trialsTraces.NoFiltMultiContSUA.(xcluster).bin1.neuralDat(:,~(excludeTrials+nprocTrials));
            
            %plot mean response accross all trials above (subplot 1) and mean response with selected trials below (subplot
            %2)
            xabs = -199:1300;
            
            figure('Renderer', 'painters', 'Position', [10 10 1000 1200]);
            subplot(2,1,1)
            meanAll = nanmean(allTrialsTraces(401:1900,:),2);
            ci_high = meanAll + 1.96*std(allTrialsTraces(401:1900,:),[],2)./sqrt(length(allTrialsTraces(1,:)));
            ci_low = meanAll - 1.96*std(allTrialsTraces(401:1900,:),[],2)./sqrt(length(allTrialsTraces(1,:)));
            plot(xabs, meanAll,'linewidth',1)
            title(strcat(sprintf(' Mean response including all trials: %d trials in %s',length(allTrialsTraces(1,:)),cluster)),'Interpreter', 'none')
            hold on
            h1= ciplot( ci_high, ci_low,[-200:1299],[40/255 40/255 40/255],0.1);
            set(h1, 'edgecolor','none')
            
            %add all msaccs onset times on plot
            for t = 1:length(trialindex)
                yl = get(gca,'ylim');
                u1= zeros(size(xabs))+yl(1);
                u2= zeros(size(xabs))+yl(1);
                u1(excluSaccs(~isnan(excluSaccs(:,t)) & (excluSaccs(:,t) > -199 & excluSaccs(:,t)<1300) ,t)) = yl(2); %only keep saccades within xabs range
                u2(excluSaccs(~isnan(excluSaccs(:,t)) & (excluSaccs(:,t) > -199 & excluSaccs(:,t)<1300) ,t)+10) = yl(2);
                uOri = (cumsum(u1)-cumsum(u2));
                u =double(ischange(uOri));
                u(u ==0) = NaN;
                idc = find(~isnan(u));
                %u(idc-1) =0;
                u(idc+1) =0;
                plot(xabs, 10*u.*uOri./uOri,'k','linewidth',1)
                hold on
                uOri(uOri ==0) = NaN;
                plot(xabs,10*uOri./uOri,'k','linewidth',1)
                hold on
            end
            
            subplot(2,1,2)
            meanSel = nanmean(selectedTrialsTraces(401:1900,:),2);
            ci_high = meanSel + 1.96*std(selectedTrialsTraces(401:1900,:),[],2)./sqrt(length(selectedTrialsTraces(1,:)));
            ci_low = meanSel - 1.96*std(selectedTrialsTraces(401:1900,:),[],2)./sqrt(length(selectedTrialsTraces(1,:)));
            plot(xabs, meanSel, 'Color', 'r', 'linewidth',1)
            hold on
            h1= ciplot( ci_high, ci_low,[-200:1299],[40/255 40/255 40/255],0.1);
            set(h1, 'edgecolor','none')
            title(strcat(sprintf('Mean response excluding trials with msaccs %d trials',length(selectedTrialsTraces(1,:)))),'Interpreter', 'none')
            xlabel('Time from stimulus onset (ms)')
            ylabel('Spike rate (spikes/s)')
            
            %saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\mean_response_accounting_msaccs_',xcluster,'.png'));
            
        catch
            cnt = cnt+1;
            disp(strcat({'missing data ' xBRdatafile}))
        end
        
        if isfield(eyeMovData,xcluster)
            figure('Renderer', 'painters', 'Position', [10 10 1000 1200]);
            subplot(4,2,1)
            plot(ampl,veloc,'o')
            xlabel('Saccade amplitude (deg)');
            ylabel('Saccade peak velocity (deg/s)');
            set(gca,'box','off');
            set(gca, 'linewidth',1)
            xlim([0 2])
            title(strcat(sprintf('Microsaccades detected in %d trials of %s',ntr,cluster)),'Interpreter', 'none')
            % Plots the traces with the labeled microsaccades
            
            if ~isempty(samples)
                subplot(4,2,[3:4])
                plot(samples(:,1), samples(:,2:3),'linewidth',1);
                hold
                yl = get(gca,'ylim');
                u1= zeros(size(samples(:,1)))+yl(1);
                u2= zeros(size(samples(:,1)))+yl(1);
                u1((saccades(:,enum.startIndex))) = yl(2);
                u2(saccades(:,enum.endIndex)) = yl(2);
                uOri = (cumsum(u1)-cumsum(u2));
                u =double(ischange(uOri));
                u(u ==0) = NaN;
                idc = find(~isnan(u));
                %u(idc-1) =0;
                u(idc+1) =0;
                plot(samples(:,1), u.*uOri./uOri,'k','linewidth',1)
                hold on
                uOri(uOri ==0) = NaN;
                plot(samples(:,1),uOri./uOri,'k','linewidth',1)
                xlim([-200 1500])
                ylim([-5 5])
                set(gca,'XTickLabel',[]);
                set(gca,'box','off');
                set(gca, 'linewidth',1)
                ylabel('Eye Position (deg)');
                legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
                
                
                subplot(4,2,[5:6])
                plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1),eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,2:3),'linewidth',1)
                hold
                yl = get(gca,'ylim');
                u3= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1)))+yl(1);
                u4= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1)))+yl(1);
                u3((eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).enum.startIndex))) = yl(2);
                u4(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).enum.endIndex)) = yl(2);
                uOri = (cumsum(u3)-cumsum(u4));
                u =double(ischange(uOri));
                u(u ==0) = NaN;
                idc = find(~isnan(u));
                %u(idc-1) =0;
                u(idc+1) =0;
                plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), u.*uOri./uOri,'k','linewidth',1)
                hold on
                uOri(uOri ==0) = NaN;
                plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), uOri./uOri,'k','linewidth',1)
                xlim([-200 1500])
                ylim([-5 5])
                set(gca,'box','off');
                set(gca,'XTickLabel',[]);
                set(gca, 'linewidth',1)
                ylabel('Eye Position (deg)');
                legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
                
                subplot(4,2,[7:8])
                plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1),eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,2:3),'linewidth',1)
                hold
                yl = get(gca,'ylim');
                u5= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1)))+yl(1);
                u6= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1)))+yl(1);
                u5((eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).enum.startIndex))) = yl(2);
                u6(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).enum.endIndex)) = yl(2);
                uOri = (cumsum(u5)-cumsum(u6));
                uOri(uOri > 0)= 1;
                u =double(ischange(uOri*4));
                u(u ==0) = NaN;
                idc = find(~isnan(u));
                %u(idc-1) =0;
                u(idc+1) =0;
                plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), u.*uOri./uOri,'k','linewidth',1)
                hold on
                uOri(uOri ==0) = NaN;
                plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), uOri./uOri,'k','linewidth',1)
                xlim([-200 1500])
                ylim([-5 5])
                set(gca,'box','off');
                set(gca, 'linewidth',1)
                xlabel('Time from stimulus onset (ms)');
                ylabel('Eye Position (deg)');
                
                legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
                
            end
        end
    end
end


%% Perform same analysis with peaks triggered to cycle peaks
indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\analysis\';
concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
%trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']); %neural data
trialsTraces =load([newdatadir 'NoFiltMultiContSUA_06212021']); %neural data +peaklocs + trial numbers obtained with  "BinocularAdaptationTrialSelection.m"
% peak triggered data

peak_trig_traces = suaPeakTrigResps(trialsTraces.NoFiltMultiContSUA);

xfilenames = fieldnames(peak_trig_traces);
cnt = 0;
eyeMovData = struct();


for i=1:length(xfilenames)
    try
        xcluster = xfilenames{i};
        cluster = xcluster(2:end);
        underscore = strfind(cluster, '_');
        session =  cluster(1:underscore(2)-1);
        directory = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\',cluster,'\');
        
        xBRdatafiles = concat_filenames.(xcluster);
        eye_info =[];
        all_codes = [];
        all_times = [];
        all_analogData =[];
        for fn =1:length(xBRdatafiles)
            xBRdatafile = xBRdatafiles{fn};
            filename   = [directory xBRdatafile(2:end)];
            if exist(strcat(filename, '.bhv'),'file')
                eye_info.(strcat(xBRdatafile,'_bhvfile')) = concatBHV(strcat(filename,'.bhv'));
                all_codes = [all_codes, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeNumbers];
                all_times = [all_times, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeTimes];
                all_analogData = [all_analogData,eye_info.(strcat(xBRdatafile,'_bhvfile')).AnalogData];
                xBaseline = eye_info.(strcat(xBRdatafile,'_bhvfile')).ScreenXresolution/4/eye_info.(strcat(xBRdatafile,'_bhvfile')).PixelsPerDegree; %since the screen monitor is split in 2 parts with the stereoscope, the center for each eye becomes the center for each side of the stereoscope (half of the half, justifying dividing by 4
            end
            
        end
        samplerate = 1000;
        %trialindex = condSelectedTrialsIdx.(xcluster);
        trialindex = trialsTraces.NoFiltMultiContSUA.(xcluster).bin1.trials; %only take trial indices of monocular stimulation
        
        ampl = [];
        veloc = [];
        ntr =0; %number of trials in which microsaccades were detected
        excludeTrials = zeros(length(trialindex),1);
        nprocTrials = zeros(length(trialindex),1);
        excluSaccs = nan(20,length(trialindex),1);
        sparedSaccs = nan(20,length(trialindex),1);
        
        for tr = 1:length(trialindex)
            
            codes                 = all_codes{trialindex(tr)};
            times                 = all_times{trialindex(tr)};
            
            if nnz(find( codes == 23))
                samples = [];
                samples(:,1) = (-1*times(codes == 23)+1) : 1 : 0 : (length(all_analogData{trialindex(tr)}.EyeSignal(:,1)) - times(codes == 23)); %23 = trial onset time %24 = trial offset time
                %samples(:,1) = 1:length(all_analogData{trialindex(tr)}.EyeSignal(:,1));
                if ~isempty(samples)
                    samples(:,2) = all_analogData{trialindex(tr)}.EyeSignal(:,1)+xBaseline; %horizontal position of the left eye in degrees baseline corrected
                    samples(:,3) = all_analogData{trialindex(tr)}.EyeSignal(:,2); %vertical position of the left eye in degrees
                    samples(:,4) = nan();
                    samples(:,5) = nan();
                    blinks = zeros(length(samples(:,1)),1);
                    recording = ClusterDetection.EyeMovRecording.Create(directory, session, samples, blinks, samplerate);
                    
                    % Runs the saccade detection
                    try
                    [saccades, stats] = recording.FindSaccades();
                    enum = ClusterDetection.SaccadeDetector.GetEnum;
                    
                    for s =1:length(saccades(:,enum.startIndex)) %loop through all microssaccades/ saccades found in one trial
                        if  ~isempty(find(saccades(s,enum.startIndex),1)) && (saccades(s,enum.startIndex) > 0 && saccades(s,enum.startIndex) < times(codes == 24)- times(codes == 23)) %if there is at least 1 microsaccade and it occurs between stim onset and stim offset
                            excludeTrials(tr) = excludeTrials(tr)+1;
                            excluSaccs(s,tr) = saccades(s,enum.startIndex);
                        else
                            if  ~isempty(find(saccades(s,enum.startIndex),1)) && (saccades(s,enum.startIndex) < 0 || saccades(s,enum.startIndex) > times(codes == 24)- times(codes == 23))
                                sparedSaccs(s,tr) = saccades(s,enum.startIndex);
                            end
                        end
                    end
                    catch
                        nprocTrials(tr) =nprocTrials(tr)+1;
                        disp(strcat({'Bad Trial',xBRdatafile}))
                    end
                    %{
                 % Plots a main sequence
                ampl = [ampl; saccades(:,enum.amplitude)];
                veloc = [veloc; saccades(:,enum.peakVelocity)];
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).stats = stats;
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).saccades = saccades;
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).enum = enum;
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).samples = samples;
                if ~all(isnan((saccades(:,enum.startIndex))))
                    ntr =ntr+1;
                end
                    %}
                end
            end
        end
       
        %eyeMovData.(xcluster).cellclass = trialsTraces.NoFiltMultiContSUA.(xcluster).cellclass;
        %end
        
        %plot mean response accross all trials above (subplot 1) and mean response with selected trials below (subplot
        %2)
        xabs = -125:124;
        
        figure('Renderer', 'painters', 'Position', [10 10 1000 1200]);
        
        for p =1:4
            pkN = sprintf('pk%d',p);
            %all traces excluding unprocessed ones due to algorithm bug
            allTrialsTraces =peak_trig_traces.(xcluster).originSUA.bin1.(pkN)(:,~nprocTrials);
           
            %selected trials neural data excluding ones with msaccs during
            %stimulation
            selectedTrialsTraces =peak_trig_traces.(xcluster).originSUA.bin1.(pkN)(:,~(excludeTrials + nprocTrials));
          
            subplot(2,4,p)
            meanAll = nanmean(allTrialsTraces,2);
            ci_high = meanAll + 1.96*std(allTrialsTraces,[],2,'omitnan')./sqrt(length(allTrialsTraces(1,:)));
            ci_low = meanAll - 1.96*std(allTrialsTraces,[],2,'omitnan')./sqrt(length(allTrialsTraces(1,:)));
            plot(xabs, meanAll,'linewidth',1)
            hold on
            h1= ciplot( ci_high, ci_low,[-125:124],[40/255 40/255 40/255],0.1);
            set(h1, 'edgecolor','none')
            set(gca,'box','off')
            if p == 1
                title(strcat(sprintf(' Mean response including all trials: %d trials',length(allTrialsTraces(1,:)))),'Interpreter', 'none')
                yl = get(gca, 'ylim');
            end
            if p> 1
               ax1 = gca;
                %ax1.YAxis.Visible = 'off';
                ylim(yl);
            end
            
            
            %ylim([0 120])
            subplot(2,4,p+4)
            meanSel = nanmean(selectedTrialsTraces,2);
            ci_high = meanSel + 1.96*std(selectedTrialsTraces,[],2,'omitnan')./sqrt(length(selectedTrialsTraces(1,:)));
            ci_low = meanSel - 1.96*std(selectedTrialsTraces,[],2,'omitnan')./sqrt(length(selectedTrialsTraces(1,:)));
            plot(xabs, meanSel, 'Color', 'r', 'linewidth',1)
            hold on
            h1= ciplot( ci_high, ci_low,[-125:124],[40/255 40/255 40/255],0.1);
            set(h1, 'edgecolor','none')
            xlabel('Time from stimulus onset (ms)')
            ylabel('Spike rate (spikes/s)')
            set(gca,'box','off')
            if p == 1
                 title(strcat(sprintf('Mean response excluding trials with msaccs %d trials',length(selectedTrialsTraces(1,:)))),'Interpreter', 'none')

            %    yl2 = get(gca, 'ylim');
           end
            ylim(yl);
            if p> 1
                ax1 = gca;
                ax1.YAxis.Visible = 'off';
               
            end
            
        end
      %  saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\peaktrigg_mean_response_accounting_msaccs_',xcluster,'.png'));
        
    catch
       cnt = cnt+1;
       disp(strcat({'missing data ' xBRdatafile}))
    end
    
    
end

%% Now use knowledge from previous steps to make a good figure for the paper

indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\analysis\';
concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
%trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']); %neural data
trialsTraces =load([newdatadir 'NoFiltMultiContSUA_06212021']); %neural data +peaklocs + trial numbers obtained with  "BinocularAdaptationTrialSelection.m"
% peak triggered data

peak_trig_traces = suaPeakTrigResps(trialsTraces.NoFiltMultiContSUA);

xfilenames = fieldnames(peak_trig_traces);
cnt = 0;
eyeMovData = struct();
for i =1:length(xfilenames)
    try
        xcluster =xfilenames{i};
        cluster = xcluster(2:end);
        underscore = strfind(cluster, '_');
        session =  cluster(1:underscore(2)-1);
        directory = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\',cluster,'\');
        xBRdatafiles = concat_filenames.(xcluster);
        eye_info =[];
        all_codes = [];
        all_times = [];
        all_analogData =[];
        for fn =1:length(xBRdatafiles)
            xBRdatafile = xBRdatafiles{fn};
            filename   = [directory xBRdatafile(2:end)];
            if exist(strcat(filename, '.bhv'),'file')
                eye_info.(strcat(xBRdatafile,'_bhvfile')) = concatBHV(strcat(filename,'.bhv'));
                all_codes = [all_codes, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeNumbers];
                all_times = [all_times, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeTimes];
                all_analogData = [all_analogData,eye_info.(strcat(xBRdatafile,'_bhvfile')).AnalogData];
                xBaseline = eye_info.(strcat(xBRdatafile,'_bhvfile')).ScreenXresolution/4/eye_info.(strcat(xBRdatafile,'_bhvfile')).PixelsPerDegree; %since the screen monitor is split in 2 parts with the stereoscope, the center for each eye becomes the center for each side of the stereoscope (half of the half, justifying dividing by 4
            end
            
        end
        samplerate = 1000;
        %trialindex = condSelectedTrialsIdx.(xcluster);
        trialindex = trialsTraces.NoFiltMultiContSUA.(xcluster).bin1.trials; %only take trial indices of monocular stimulation
        
        ampl = [];
        veloc = [];
        ntr =0; %number of trials in which microsaccades were detected
        excludeTrials = zeros(length(trialindex),1);
        nprocTrials = zeros(length(trialindex),1);
        excluSaccs = nan(20,length(trialindex),1);
        sparedSaccs = nan(20,length(trialindex),1);
        
        for tr = 1:length(trialindex)
            
            codes                 = all_codes{trialindex(tr)};
            times                 = all_times{trialindex(tr)};
            
            if nnz(find( codes == 23))
                samples = [];
                samples(:,1) = (-1*times(codes == 23)+1) : 1 : 0 : (length(all_analogData{trialindex(tr)}.EyeSignal(:,1)) - times(codes == 23)); %trigger time points on stimulus onset time for it to be 0. Everything before that point is then negative %23 = stimulus onset time, time is measured relative to each trial onset recorded %24 = trial/stimulus offset time
                %samples(:,1) = 1:length(all_analogData{trialindex(tr)}.EyeSignal(:,1));
                if ~isempty(samples)
                    samples(:,2) = all_analogData{trialindex(tr)}.EyeSignal(:,1)+xBaseline; %horizontal position of the left eye in degrees baseline corrected
                    samples(:,3) = all_analogData{trialindex(tr)}.EyeSignal(:,2); %vertical position of the left eye in degrees
                    samples(:,4) = nan();
                    samples(:,5) = nan();
                    blinks = zeros(length(samples(:,1)),1);
                    recording = ClusterDetection.EyeMovRecording.Create(directory, session, samples, blinks, samplerate);
                    
                    % Runs the saccade detection
                    try
                        [saccades, stats] = recording.FindSaccades();
                        enum = ClusterDetection.SaccadeDetector.GetEnum;
                        
                        for s =1:length(saccades(:,enum.startIndex)) %loop through all microssaccades/ saccades found in one trial
                            if  ~isempty(find(saccades(s,enum.startIndex),1)) && (saccades(s,enum.startIndex) > 0 && saccades(s,enum.startIndex) < times(codes == 24)- times(codes == 23)) %if there is at least 1 microsaccade and it occurs between stim onset and stim offset
                                excludeTrials(tr) = excludeTrials(tr)+1; %trials to exclude since they have microsaccades between stim  onset and stim offset
                                excluSaccs(s,tr) = saccades(s,enum.startIndex); %save microsaccade onset time
                            else
                                if  ~isempty(find(saccades(s,enum.startIndex),1)) && (saccades(s,enum.startIndex) < 0 || saccades(s,enum.startIndex) > times(codes == 24)- times(codes == 23))
                                    sparedSaccs(s,tr) = saccades(s,enum.startIndex); %save microsaccades occurring outside of stim onset-stim offset time
                                end
                            end
                        end
                    catch
                        nprocTrials(tr) =nprocTrials(tr)+1; %trials for which the saccade detection above did not work
                        disp(strcat({'Bad Trial',xBRdatafile}))
                    end
                    
                 % Plots a main sequence
                ampl = [ampl; saccades(:,enum.amplitude)];
                veloc = [veloc; saccades(:,enum.peakVelocity)];
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).stats = stats;
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).saccades = saccades;
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).enum = enum;
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).samples = samples;
                if ~all(isnan((saccades(:,enum.startIndex))))
                    ntr =ntr+1;
                end
                   
                end
                 
            end
        end
        eyeMovData.(xcluster).nprocTrials = nprocTrials; 
        eyeMovData.(xcluster).msaccTrials = excludeTrials;
        eyeMovData.(xcluster).msaccOn = excluSaccs;
        eyeMovData.(xcluster).cellclass = trialsTraces.NoFiltMultiContSUA.(xcluster).cellclass;
        %end
        catch
        cnt = cnt+1;
        disp(strcat({'missing data ' xBRdatafile}))
    end
end

   %% plot mean response accross all trials above (subplot 1) and mean
        %response with selected trials overlayed. Subplot 2: plot
        %difference of the means
        %colors
        nlines =7;
        cmaps =cbrewer2('Oranges', nlines);

        xabs = -125:124;
   for i =1:length(xfilenames) 
       xcluster =xfilenames{i};
       %xcluster ='x160629_I_p03_uclust62_cinterocdrft_stab_fft_sig';
       %xcluster ='x180827_I_p02_uclust8_cinterocdrft_stab_fft_sig';
       
       if isfield(peak_trig_traces, xcluster) && isfield(eyeMovData,xcluster)
        figure('Renderer', 'painters', 'Position', [10 10 1000 1200]);
        for p =1:4
            pkN = sprintf('pk%d',p);
            %all traces excluding unprocessed ones due to algorithm bug
            allTrialsTraces =peak_trig_traces.(xcluster).originSUA.bin1.(pkN)(:,~eyeMovData.(xcluster).nprocTrials);
            %selected trials neural data excluding ones with msaccs during
            %stimulation
            selectedTrialsTraces =peak_trig_traces.(xcluster).originSUA.bin1.(pkN)(:,~( eyeMovData.(xcluster).msaccTrials  + eyeMovData.(xcluster).nprocTrials));
            
            subplot(2,4,p)
            meanAll = nanmean(allTrialsTraces,2);
            %ci_high = meanAll + 1.96*std(allTrialsTraces,[],2,'omitnan')./sqrt(length(allTrialsTraces(1,:)));
            %ci_low = meanAll - 1.96*std(allTrialsTraces,[],2,'omitnan')./sqrt(length(allTrialsTraces(1,:)));
            ci_high = meanAll + std(allTrialsTraces,[],2,'omitnan');
            ci_low = meanAll - std(allTrialsTraces,[],2,'omitnan');
          
            plot(xabs, meanAll,'linewidth',1,'col',[180/255 180/255 180/255])
            hold on
            h1= ciplot( ci_high, ci_low,[-125:124],[40/255 40/255 40/255],0.1);
            set(h1, 'edgecolor','none')
            hold on
            meanSel = nanmean(selectedTrialsTraces,2);
            %ci_high = meanSel + 1.96*std(selectedTrialsTraces,[],2,'omitnan')./sqrt(length(selectedTrialsTraces(1,:)));
            %ci_low = meanSel - 1.96*std(selectedTrialsTraces,[],2,'omitnan')./sqrt(length(selectedTrialsTraces(1,:)));
            ci_high = meanSel + std(selectedTrialsTraces,[],2,'omitnan');
            ci_low = meanSel - std(selectedTrialsTraces,[],2,'omitnan');
      
            plot(xabs, meanSel, 'col', cmaps(4,:), 'linewidth',1)
            hold on
            h2= ciplot( ci_high, ci_low,[-125:124],cmaps(4,:),0.1);
            set(h2, 'edgecolor','none')
            ylabel('Spike rate (spikes/s)')
            
            set(gca,'box','off')
            if p == 1
                sgtitle(strcat({sprintf(' Mean response including all trials (%d trials) vs Mean response excluding trials with msaccs (%d trials)',length(allTrialsTraces(1,:)),length(selectedTrialsTraces(1,:))), xcluster}),'Interpreter', 'none')
            end
            if p> 1
                ax1 = gca;
                ax1.YAxis.Visible = 'off';
            end
            ylim([0 270]);
            
            %ylim([0 120])
            subplot(2,4,p+4)
            diff = meanAll - meanSel;
            plot(xabs, diff, 'Color', [40/255 40/255 40/255], 'linewidth',1)
            xlabel('Time from stimulus onset (ms)')
            ylabel('Spike rate difference (spikes/s)')
            set(gca,'box','off')
            if p == 1
                title(strcat(sprintf('Mean(all trials) - Mean(msacc excluded trials) difference')),'Interpreter', 'none')
            end
            if p> 1
                ax1 = gca;
                ax1.YAxis.Visible = 'off';
            end
            ylim([-7 18]);
        end
       % saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\peaktrigg_diffmean_response_accounting_msaccs_',xcluster,'.svg'));
       end
   end
   
   
   %% plot overall mean difference of the means across single units
   xabs = -125:124;
   all_diffs = nan(250,4,length(xfilenames));
   for i =1:length(xfilenames)
       xcluster =xfilenames{i};
       if isfield(peak_trig_traces, xcluster) && isfield(eyeMovData,xcluster)
           
           for p =1:4
               pkN = sprintf('pk%d',p);
               %all traces excluding unprocessed ones due to algorithm bug
               allTrialsTraces =peak_trig_traces.(xcluster).originSUA.bin1.(pkN)(:,~eyeMovData.(xcluster).nprocTrials);
               %selected trials neural data excluding ones with msaccs during
               %stimulation
               selectedTrialsTraces =peak_trig_traces.(xcluster).originSUA.bin1.(pkN)(:,~( eyeMovData.(xcluster).msaccTrials  + eyeMovData.(xcluster).nprocTrials));
               meanAll = nanmean(allTrialsTraces,2);
               meanSel = nanmean(selectedTrialsTraces,2);
               all_diffs(:,p,i) = meanAll - meanSel;
           end
       end
   end
  
   nlines =7;
   cmaps =cbrewer2('GnBu', nlines);

   figure('Renderer', 'painters', 'Position', [10 10 1500 800]);
   meanDiffs= nanmean(all_diffs,3);
   for p =1:4
            subplot(1,4,p)
            ci_high = meanDiffs(:,p) + 1.96*std(squeeze(all_diffs(:,p,:)),[],2,'omitnan')./sqrt(size(all_diffs,3));
            ci_low = meanDiffs(:,p) - 1.96*std(squeeze(all_diffs(:,p,:)),[],2,'omitnan')./sqrt(size(all_diffs,3));
            plot(xabs, meanDiffs(:,p), 'col', cmaps(4,:), 'linewidth',1)
            hold on
            h2= ciplot( ci_high, ci_low,[-125:124],cmaps(4,:),0.3);
            set(h2, 'edgecolor','none')
             if p == 1
                sgtitle(strcat(sprintf('Mean(all trials) - Mean(msacc excluded trials) difference mean across all units')),'Interpreter', 'none')
            end
            if p> 1
                ax1 = gca;
                ax1.YAxis.Visible = 'off';
            end
            set(gca,'box','off')
            ylim([-2 2]);
            %saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\peaktrigg_diffmean_response_withandwithout_msaccs_all.svg'));
   end
   
%% Plot main sequence 

xcluster ='x180827_I_p02_uclust8_cinterocdrft_stab_fft_sig';
ampl = [];
veloc = [];
fieldn = fieldnames(eyeMovData.(xcluster));
for tr =1:length(cell2mat(strfind(fieldn,'t')))
    trIdx = fieldn{tr};
ampl = [ampl; eyeMovData.(xcluster).(trIdx).saccades(:,eyeMovData.(xcluster).(trIdx).enum.amplitude)];
veloc = [veloc; eyeMovData.(xcluster).(trIdx).saccades(:,eyeMovData.(xcluster).(trIdx).enum.peakVelocity)];
end
%{
figure('Renderer', 'painters', 'Position', [10 10 600 700]);
plot(ampl,veloc,'o')
xlabel('Saccade amplitude (deg)');
ylabel('Saccade peak velocity (deg/s)');
set(gca,'box','off');
set(gca, 'linewidth',1)
xlim([0 2])
title(strcat(sprintf('Main sequence example in %d trials',length(cell2mat(strfind(fieldn,'t'))))),'Interpreter', 'none')
%}

clear g10
figure('Position',[100 100 600 450]);
g=gramm('x',ampl,'y',veloc);
g.set_names('x','Saccade amplitude (deg)','y','Saccade peak velocity (deg/s)');
g.set_color_options('chroma',0,'lightness',60);
g.stat_glm('geom','area','disp_fit',false);
g.set_title(sprintf('Main sequence example in %d trials',length(cell2mat(strfind(fieldn,'t'))))); %Title must be provided before the first draw() call
g.draw();
snapnow;
g.update();
g.set_color_options();
g.geom_point();
g.draw();
saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\', sprintf('main_sequence_example%s',xcluster), '.svg'));


%% Frequency of msaccs accross trials broken down by animal I and B

%x = frequency/trial during stimulation
%y = number of single units
%left = I monkey
%right = B monkey

filenames = fieldnames(eyeMovData);
freq = nan(length(fieldnames(eyeMovData)), 2); %2 = 2 monkeys
icnt = 0;
bcnt = 0;
for i = 1:length(fieldnames(eyeMovData))
    xcluster = filenames{i};
    if contains(xcluster,'_I_')
        icnt = icnt+1;
        freq(i,1) = sum(eyeMovData.(xcluster).msaccTrials)/(length(find(~eyeMovData.(xcluster).nprocTrials)));
    else 
        if contains(xcluster,'_B_')
            bcnt = bcnt+1;
           freq(i,2) = sum(eyeMovData.(xcluster).msaccTrials)/(length(find(~eyeMovData.(xcluster).nprocTrials)));
 
        end
    end
    
end

%linFreq = [freq(~isnan(freq(:,1)),1);freq(~isnan(freq(:,2)),2)];
linFreq = [freq(:,1);freq(:,2)]';
%animal = [repmat('I', [length(find(~isnan(freq(:,1)))),1]);repmat('B', [length(find(~isnan(freq(:,2)))),1])];
animal = [repmat({'I'}, [length(freq(:,1)),1]);repmat({'B'}, [length(freq(:,2)),1])]';

%unit = [1:length(find(~isnan(freq(:,1)))),1:length(find(~isnan(freq(:,2))))];
unit = [1:length(freq(:,1)),1:length(freq(:,2))];
%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('GnBu', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

clear g
g(1,1)=gramm('x',linFreq,'y',unit,'color',animal);
%g(1,1).stat_bin('nbins',25,'geom','overlaid_bar');
g(1,1).stat_bin('normalization','probability','nbins',10,'geom','overlaid_bar');
%g(1,1).stat_density();
%g(1,1).set_color_options('map', [251/255 154/255 153/255;160/255 160/255 160/255]);
g(1,1).set_color_options('map',[cmaps(3).map(4,:);cmaps(1).map(4,:)]); 
g(1,1).set_names('x','Microsaccade frequency (normalized per trial)','color','Legend','row','','y','Count');
g(1,1).set_title({'Microsaccade Frequency per trial per animal'});
f = figure('Position',[100 100 800 1000]);
g.draw();
set(f,'position',get(f,'position').*[1 1 1.15 1])
saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\msacc_frequency_per_trial_I_B.svg'));

%% Replot Frequency with scatter plot
clear g
f = figure('Position',[100 100 800 1000]);
set(f,'position',get(f,'position').*[1 1 1.15 1])

  %jitter
  %gramm('x',linFreq,'y',unit,'color',animal);
g(1,1) = gramm('x',animal,'y', linFreq, 'color',animal); %[unitnb,1:length(unitnb(adapType(:,t)))]
g(1,1).geom_jitter('width',0.4,'height',0); %Scatter plot
g(1,1).set_color_options('map',[cmaps(3).map(4,:);cmaps(1).map(4,:)]); 
g(1,1).axe_property( 'xlim',[0 4] , 'ylim',[0 3.5]); 
% add confidence interval 95%

ci_low = [nanmean(linFreq(strcmp(animal, 'B'))) - std(linFreq(strcmp(animal, 'B')),0,'omitnan')/sqrt(length(linFreq(strcmp(animal, 'B')))); nanmean(linFreq(strcmp(animal, 'I'))) - std(linFreq(strcmp(animal, 'I')),0,'omitnan')/sqrt(length( linFreq(strcmp(animal, 'I')))) ];
ci_high = [nanmean(linFreq(strcmp(animal, 'B'))) + std(linFreq(strcmp(animal, 'B')),0,'omitnan')/sqrt(length(linFreq(strcmp(animal, 'B')))); nanmean(linFreq(strcmp(animal, 'I'))) + std(linFreq(strcmp(animal, 'I')),0,'omitnan')/sqrt(length( linFreq(strcmp(animal, 'I')))) ];
g(1,1).update('x',[1;2], 'y', [nanmean(linFreq(strcmp(animal, 'B'))); nanmean(linFreq(strcmp(animal, 'I')))],...
    'ymin',ci_low,'ymax',ci_high,'color',[1;2]);
g(1,1).geom_point('dodge',0.5);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).set_color_options('map',[cmaps(3).map(4,:);cmaps(1).map(4,:)]); 
g(1,1).axe_property('xlim',[0 4]); 
g(1,1).set_point_options('base_size',7);
%}
%bar
g(1,2)=gramm('x',linFreq,'y',unit, 'color',animal);
g(1,2).stat_bin('nbins',25,'geom','overlaid_bar');
g(1,2).stat_density();
g(1,2).set_color_options('map',[cmaps(3).map(4,:);cmaps(1).map(4,:)]); 
g(1,2).axe_property('xlim',[0 3.5], 'ylim', [0 6]); 
g(1,2).set_names('x','Microsaccade Frequency (normalized per trial)','color','Legend','row','','y','Count');
g(1,2).coord_flip();
g.draw();
saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\msacc_frequency_per_trial_I_B_scatter.svg'));


%% Microsaccade amplitude distributions broken down by monkey

filenames = fieldnames(eyeMovData);
ampI = []; %monkey I
ampB = []; %monkey B

icnt = 0;
bcnt = 0;
for i = 1:length(fieldnames(eyeMovData))
    xcluster = filenames{i};
    fieldn = fieldnames(eyeMovData.(xcluster));

    if contains(xcluster,'_I_')
        icnt = icnt+1;
        for tr =1:length(cell2mat(strfind(fieldn,'t')))
            trIdx = fieldn{tr};
            ampI= [ampI; eyeMovData.(xcluster).(trIdx).saccades(:,eyeMovData.(xcluster).(trIdx).enum.amplitude)]; %loop through trials
        end
            else 
        if contains(xcluster,'_B_')
            bcnt = bcnt+1;
            for tr =1:length(cell2mat(strfind(fieldn,'t')))
                trIdx = fieldn{tr};
                ampB= [ampB; eyeMovData.(xcluster).(trIdx).saccades(:,eyeMovData.(xcluster).(trIdx).enum.amplitude)]; %loop through trials
            end
        end
    end
    
end

linAmp = log([ampI(ampI<2);ampB(ampB<2)])';
animal = [repmat({'I'}, [length(ampI(ampI<2)),1]);repmat({'B'}, [length(ampB(ampB<2)),1])]';
animal = animal(isfinite(linAmp));
msacc = [1:length(ampI(ampI<2)),1:length(ampB(ampB<2))];
msacc = msacc(isfinite(linAmp));
linAmp = linAmp(isfinite(linAmp));

nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('GnBu', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

clear g
g(1,1)=gramm('x',linAmp,'y',msacc,'color',animal);
%g(1,1).stat_bin('nbins',25,'geom','overlaid_bar');
g(1,1).stat_bin('normalization','probability','nbins',50,'geom','overlaid_bar');
g(1,1).stat_density();
g(1,1).set_color_options('map',[cmaps(3).map(4,:);cmaps(1).map(4,:)]); 
g(1,1).set_names('x','Microsaccade amplitude (log)','color','Legend','row','','y','Count');
g(1,1).set_title({'Microsaccade Amplitude distribution per animal'});

f = figure('Position',[100 100 800 1000]);
g.draw();
%set(f,'position',get(f,'position').*[1 1 1.15 1])
saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\msacc_amplitudes_density_I_B.svg'));
%% Replot Amplitude with scatter plot

clear g
f = figure('Position',[100 100 800 1000]);
set(f,'position',get(f,'position').*[1 1 1.15 1])
  %jitter
g(1,1) = gramm('x',animal,'y',linAmp, 'color',animal); 
g(1,1).geom_jitter('width',0.4,'height',0); %Scatter plot
g(1,1).set_color_options('map',[cmaps(3).map(4,:);cmaps(1).map(4,:)]); 
g(1,1).axe_property( 'xlim',[0 4] , 'ylim',[-5 1]); 
% add confidence interval 95%
ci_low = [nanmean(linAmp(strcmp(animal, 'B'))) - std(linAmp(strcmp(animal, 'B')),0,'omitnan')/sqrt(length(linAmp(strcmp(animal, 'B')))); nanmean(linAmp(strcmp(animal, 'I'))) - std(linAmp(strcmp(animal, 'I')),0,'omitnan')/sqrt(length(linAmp(strcmp(animal, 'I')))) ];
ci_high = [nanmean(linAmp(strcmp(animal, 'B'))) + std(linAmp(strcmp(animal, 'B')),0,'omitnan')/sqrt(length(linAmp(strcmp(animal, 'B')))); nanmean(linAmp(strcmp(animal, 'I'))) + std(linAmp(strcmp(animal, 'I')),0,'omitnan')/sqrt(length(linAmp(strcmp(animal, 'I')))) ];
g(1,1).update('x',[1;2], 'y', [nanmean(linAmp(strcmp(animal, 'B'))); nanmean(linAmp(strcmp(animal, 'I')))],...
    'ymin',ci_low,'ymax',ci_high,'color',[1;2]);
g(1,1).geom_point('dodge',0.5);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).set_color_options('map',[cmaps(3).map(4,:);cmaps(1).map(4,:)]); 
g(1,1).axe_property('xlim',[0 4]); 
g(1,1).set_point_options('base_size',7);
%}
%bar
g(1,2)=gramm('x',linAmp,'y',msacc, 'color',animal);
%g(1,2).stat_bin('nbins',25,'geom','overlaid_bar');
g(1,2).stat_bin('normalization','probability','nbins',50,'geom','overlaid_bar');
g(1,2).stat_density();
g(1,2).set_color_options('map',[cmaps(3).map(4,:);cmaps(1).map(4,:)]); 
g(1,2).axe_property('xlim',[-5 1], 'ylim', [0 1]); 
g(1,2).set_names('x','Microsaccade Amplitude (log)','color','Legend','row','','y','Count');
g(1,2).coord_flip();
g.draw();
saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\msacc_amplitudes_density_I_B_scatter.svg'));

%% New Idea: compare frequency, amplitude, and traces of microsaccades in facilitated units vs suppressed units

%1) Frequency
filenames = fieldnames(peak_trig_traces);
freq = nan(length(fieldnames(peak_trig_traces)), 1); 

for i = 1:length(fieldnames(peak_trig_traces))
    xcluster = filenames{i};
    if isfield(peak_trig_traces, xcluster) && isfield(eyeMovData,xcluster)
       
            freq(i) = sum(eyeMovData.(xcluster).msaccTrials)/(length(find(~eyeMovData.(xcluster).nprocTrials)));

    end
end



unit = [1:length(freq(:,1)),1:length(freq(:,2))]; 

%get pvalues from lmer results with Dunnett correction
pvalues = dlmread('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\lmer_results_orig_03032020_corrected_dunnett.csv', ',', 1,1);
pk1pk4pval = pvalues(~isnan(pvalues(:,3)),3);%clear out nans
 
%get mean pk1-pk4 differences for all units
meanpk1 = nan(length(xfilenames),1);
meanpk4 = nan(length(xfilenames),1);

for i =1:length(xfilenames)
       xcluster =xfilenames{i};
       if isfield(peak_trig_traces, xcluster) && isfield(eyeMovData,xcluster)
         
               %all traces excluding unprocessed ones due to algorithm bug
               pk1Traces =peak_trig_traces.(xcluster).originSUA.bin1.pk1(:,~eyeMovData.(xcluster).nprocTrials);
               pk4Traces = peak_trig_traces.(xcluster).originSUA.bin1.pk4(:,~eyeMovData.(xcluster).nprocTrials);
               
               meanpk1(i) = max(mean(pk1Traces,2));
               meanpk4(i) = max(mean(pk4Traces,2));

       end
end

%get logicals to extract suppressed and facilitated units from the
%population
facUnits = pk1pk4pval < 0.05 & meanpk1 - meanpk4 < 0;
suppUnits = pk1pk4pval < 0.05 & meanpk1 - meanpk4 > 0;

unitBhv = cell(length(facUnits),1); %is the unit adapting or not? facilitated or suppressed?
for i =1 :length(facUnits)
    if facUnits(i) == 1
        unitBhv(i) = {'Facilitated'};
    elseif suppUnits(i) == 1
        unitBhv(i) = {'Suppressed'};
    end
end
index=cellfun(@isempty,unitBhv);
fullIndex = ~index;
linUnitBhv = [unitBhv(fullIndex); repmat({'Population'}, [length(unitBhv),1])];
linFreq = [freq(fullIndex);freq(:)]';

nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;


%plot figure
clear g
f = figure('Position',[100 100 800 1000]);
set(f,'position',get(f,'position').*[1 1 1.15 1])
  %jitter
g(1,1) = gramm('x',linUnitBhv ,'y',linFreq, 'color',linUnitBhv); 
g(1,1).geom_jitter('width',0.4,'height',0); %Scatter plot
g(1,1).set_color_options('map',[cmaps(2).map(4,:);cmaps(1).map(5,:);cmaps(3).map(5,:)]); 
%g(1,1).axe_property( 'xlim',[0 4] , 'ylim',[-5 1]); 
% add confidence interval 95%
ci_low = [nanmean(linFreq(strcmp(linUnitBhv, 'Facilitated'))) - std(linFreq(strcmp(linUnitBhv, 'Facilitated')),0,'omitnan')/sqrt(length(linFreq(strcmp(linUnitBhv, 'Facilitated'))));nanmean(linFreq(strcmp(linUnitBhv, 'Population'))) - std(linFreq(strcmp(linUnitBhv, 'Population')),0,'omitnan')/sqrt(length(linFreq(strcmp(linUnitBhv, 'Population')))); nanmean(linFreq(strcmp(linUnitBhv, 'Suppressed'))) - std(linFreq(strcmp(linUnitBhv, 'Suppressed')),0,'omitnan')/sqrt(length(linFreq(strcmp(linUnitBhv, 'Suppressed')))) ];
ci_high = [nanmean(linFreq(strcmp(linUnitBhv, 'Facilitated'))) + std(linFreq(strcmp(linUnitBhv, 'Facilitated')),0,'omitnan')/sqrt(length(linFreq(strcmp(linUnitBhv, 'Facilitated'))));nanmean(linFreq(strcmp(linUnitBhv, 'Population'))) + std(linFreq(strcmp(linUnitBhv, 'Population')),0,'omitnan')/sqrt(length(linFreq(strcmp(linUnitBhv, 'Population')))); nanmean(linFreq(strcmp(linUnitBhv, 'Suppressed'))) + std(linFreq(strcmp(linUnitBhv, 'Suppressed')),0,'omitnan')/sqrt(length(linFreq(strcmp(linUnitBhv, 'Suppressed')))) ];
g(1,1).update('x',[1;2;3], 'y', [nanmean(linFreq(strcmp(linUnitBhv, 'Facilitated')));nanmean(linFreq(strcmp(linUnitBhv, 'Population'))); nanmean(linFreq(strcmp(linUnitBhv, 'Suppressed')))],...
    'ymin',ci_low,'ymax',ci_high,'color',[1;2;3]);
g(1,1).geom_point('dodge',0.5);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).set_color_options('map',[cmaps(2).map(4,:);cmaps(1).map(5,:);cmaps(3).map(5,:)]); 
%g(1,1).axe_property('xlim',[0 4]); 
g(1,1).set_point_options('base_size',7);
g.draw();
saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\msacc_frequency_adapt_scatter.svg'));

%2) Amplitude
%3) Eye Movement Trace

%% Other idea:  Take one facilitated unit and retrigger peak responses to microsaccade onset times
%% then compare to the peak aligned average
   %get pvalues from lmer results with Dunnett correction
pvalues = dlmread('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\lmer_results_orig_03032020_corrected_dunnett.csv', ',', 1,1);
pk1pk4pval = pvalues(~isnan(pvalues(:,3)),3);%clear out nans
      

nlines =7;
        cmaps =cbrewer2('Oranges', nlines);

        xabs = -125:124;
   for i =1:length(xfilenames) 
       xcluster =xfilenames{i};
       

       if isfield(peak_trig_traces, xcluster) && isfield(eyeMovData,xcluster)
        figure('Renderer', 'painters', 'Position', [10 10 1000 1200]);
        for p =1:4
            pkN = sprintf('pk%d',p);
            %all traces excluding unprocessed ones due to algorithm bug
            allTrialsTraces =peak_trig_traces.(xcluster).originSUA.bin1.(pkN)(:,~eyeMovData.(xcluster).nprocTrials);
            %selected trials neural data retriggered to msaccs onset times
            %1)retrigger msac onset times within each peak response window
            %2)retrigger peak responses to msac onset times
            selectedTrialsTraces =peak_trig_traces.(xcluster).originSUA.bin1.(pkN)(:,~eyeMovData.(xcluster).nprocTrials);
            
            subplot(2,4,p)
            meanAll = nanmean(allTrialsTraces,2);
            %ci_high = meanAll + 1.96*std(allTrialsTraces,[],2,'omitnan')./sqrt(length(allTrialsTraces(1,:)));
            %ci_low = meanAll - 1.96*std(allTrialsTraces,[],2,'omitnan')./sqrt(length(allTrialsTraces(1,:)));
            ci_high = meanAll + std(allTrialsTraces,[],2,'omitnan');
            ci_low = meanAll - std(allTrialsTraces,[],2,'omitnan');
          
            plot(xabs, meanAll,'linewidth',1,'col',[180/255 180/255 180/255])
            hold on
            h1= ciplot( ci_high, ci_low,[-125:124],[40/255 40/255 40/255],0.1);
            set(h1, 'edgecolor','none')
            hold on
            meanSel = nanmean(selectedTrialsTraces,2);
            %ci_high = meanSel + 1.96*std(selectedTrialsTraces,[],2,'omitnan')./sqrt(length(selectedTrialsTraces(1,:)));
            %ci_low = meanSel - 1.96*std(selectedTrialsTraces,[],2,'omitnan')./sqrt(length(selectedTrialsTraces(1,:)));
            ci_high = meanSel + std(selectedTrialsTraces,[],2,'omitnan');
            ci_low = meanSel - std(selectedTrialsTraces,[],2,'omitnan');
      
            plot(xabs, meanSel, 'col', cmaps(4,:), 'linewidth',1)
            hold on
            h2= ciplot( ci_high, ci_low,[-125:124],cmaps(4,:),0.1);
            set(h2, 'edgecolor','none')
            ylabel('Spike rate (spikes/s)')
            
            set(gca,'box','off')
            if p == 1
                sgtitle(strcat({sprintf(' Mean response including all trials (%d trials) vs Mean response excluding trials with msaccs (%d trials)',length(allTrialsTraces(1,:)),length(selectedTrialsTraces(1,:))), xcluster}),'Interpreter', 'none')
            end
            if p> 1
                ax1 = gca;
                ax1.YAxis.Visible = 'off';
            end
            ylim([0 270]);
            
            %ylim([0 120])
            subplot(2,4,p+4)
            diff = meanAll - meanSel;
            plot(xabs, diff, 'Color', [40/255 40/255 40/255], 'linewidth',1)
            xlabel('Time from stimulus onset (ms)')
            ylabel('Spike rate difference (spikes/s)')
            set(gca,'box','off')
            if p == 1
                title(strcat(sprintf('Mean(all trials) - Mean(msacc excluded trials) difference')),'Interpreter', 'none')
            end
            if p> 1
                ax1 = gca;
                ax1.YAxis.Visible = 'off';
            end
            ylim([-7 18]);
        end
       % saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\peaktrigg_diffmean_response_accounting_msaccs_',xcluster,'.svg'));
       end
   end
