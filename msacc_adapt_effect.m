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


for i=31:length(xfilenames)
    %try
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
    
    for tr = 1:length(trialindex)
        
        %if trialindex(tr) <= length(eye_info.CodeNumbers)
        codes                 = all_codes{trialindex(tr)};
        times                 = all_times{trialindex(tr)};
        
        
        %timestamps = zeros(length( eye_info.AnalogData{1,1}.EyeSignal(:,1)),1);
        %timestamps(trialindex) = trialindex;
        if nnz(find( codes == 23))
            samples = [];
            samples(:,1) = (-1*times(codes == 23)+1) : 1 : 0 : (length(all_analogData{trialindex(tr)}.EyeSignal(:,1)) - times(codes == 23));
            %(-1*times(codes == 23)+1) : 1 : 0 : (length(eye_info.AnalogData{1,1}.EyeSignal(:,1)) - times(codes == 23)); %timestamps of the recording in miliseconds
            %1:length(eye_info.AnalogData{1,1}.EyeSignal(:,1)); %timestamps of the recording in miliseconds
            if ~isempty(samples)
                samples(:,2) = all_analogData{trialindex(tr)}.EyeSignal(:,1)+xBaseline; %horizontal position of the left eye in degrees baseline corrected
                samples(:,3) = all_analogData{trialindex(tr)}.EyeSignal(:,2); %vertical position of the left eye in degrees
                samples(:,4) = nan();
                samples(:,5) = nan();
                blinks = zeros(length(samples(:,1)),1);
                recording = ClusterDetection.EyeMovRecording.Create(directory, session, samples, blinks, samplerate);
                
                % Runs the saccade detection
                [saccades, stats] = recording.FindSaccades();
                
                % Plots a main sequence
                enum = ClusterDetection.SaccadeDetector.GetEnum;
                ampl = [ampl; saccades(:,enum.amplitude)];
                veloc = [veloc; saccades(:,enum.peakVelocity)];
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).saccades = saccades;
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).enum = enum;
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).stats = stats;
                eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).samples = samples;
                if ~all(isnan((saccades(:,enum.startIndex))))
                    ntr =ntr+1;
                end
                
            end
        end
    end
    eyeMovData.(xcluster).cellclass = trialsTraces.NoFiltMultiContSUA.(xcluster).cellclass;
    %end
    
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
end