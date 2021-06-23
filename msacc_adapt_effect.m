%In this script the goal is to determine whether a difference can be
%observed between trials selected without accounting for microsaccades, and
%trials selected by rejecting those that include microsaccades

indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\analysis\';
concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
%trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']); %neural data
trialsTraces =load([newdatadir 'NoFiltMultiContSUA_06212021']); %neural data +peaklocs + trial numbers obtained with  "BinocularAdaptationTrialSelection.m"

xfilenames = fieldnames(trialsTraces.NoFiltMultiContSUA);
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
        title(strcat(sprintf(' Mean response including all trials: %d trials',length(allTrialsTraces(1,:)))),'Interpreter', 'none')
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
        
        saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\mean_response_accounting_msaccs_',xcluster,'.png'));
  
    catch
       cnt = cnt+1;
       disp(strcat({'missing data ' xBRdatafile}))
    end
    
    %{
       figure('Renderer', 'painters', 'Position', [10 10 1000 1200]);
        subplot(4,2,1)
        plot(ampl,veloc,'o')
        xlabel('Saccade amplitude (deg)');
        ylabel('Saccade peak velocity (deg/s)');
        set(gca,'box','off');
        set(gca, 'linewidth',1)
        xlim([0 2])
        title(strcat(sprintf('Microsaccades detected in %d trials',ntr)),'Interpreter', 'none')
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
    %}
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
        saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\plots\peaktrigg_mean_response_accounting_msaccs_',xcluster,'.png'));
        
    catch
       cnt = cnt+1;
       disp(strcat({'missing data ' xBRdatafile}))
    end
    
    
end
