%This script was developped to assess the number of trials per session
%the change in recording depth within a single unit file
%the duration of a session
%Loic -10/12/2022


%% This section is only necessary to transfer the nev and bhv files we need
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_pooled_09292021';
selected_data = load(strcat(allfilename, '.mat'));
selected_data =selected_data.selected_data_pooled;
xFilenames = fieldnames(selected_data);
PATHS = {'E:\LGN_data_TEBA\rig021\', 'E:\LGN_data_TEBA\rig021_2\', 'E:\LGN_data_TEBA\rig022\', 'E:\LGN_data_TEBA\rig022_2\'};
for p =1:length(PATHS)
    for i = 1:length(xFilenames)
        selectxFile = xFilenames{i};
        selectDate = selectxFile(2:9);
        bhvDir = strcat(PATHS{p}, selectDate);
        if exist(bhvDir, 'dir')
            sourcedir = bhvDir;
            %fbhv = dir(strcat(sourcedir,'\*cinterocdrft*.bhv'));
            fbhv = dir(strcat(sourcedir,'\*cinterocdrft*.nev'));
            fbhvNames = {fbhv.name};
            %destdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\session_trial_count\bhv_all\',selectDate);
            destdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\session_trial_count\nev_all\',selectDate);

            mkdir(destdir)
            if  ~isempty(fbhvNames)
                for n =1:length(fbhvNames)
                    [status, msg] = copyfile( strcat(sourcedir,'\',fbhv(n).name),  destdir);
                    %[status, msg] = copyfile( strcat(sourcedir,'\*cinterocdrft*.nev'),  destdir);
                    %[status, msg] = copyfile( strcat(sourcedir,'\*.gCINTEROCDRFTGrating_di'),  destdir);
                    %[status, msg] = copyfile( strcat(sourcedir,'\*cinterocdrft*.bhv'),  destdir);
                end
            end
        end
    end
end


%% total number of trials per session (including the non successful ones)


indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\analysis\';
concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
% newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
% trialsTraces =load([newdatadir 'NoFiltMultiContSUA_06212021']); %neural data +peaklocs + trial numbers obtained with  "BinocularAdaptationTrialSelection.m"
% xfilenames = fieldnames(trialsTraces.NoFiltMultiContSUA);
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_pooled_09292021';
selected_data = load(strcat(allfilename, '.mat'));
selected_data =selected_data.selected_data_pooled;
xFilenames = fieldnames(selected_data);

cnt =0;
trialCnt = [];
for i=1:length(xFilenames)
    try
        xcluster = xFilenames{i};
        cluster = xcluster(2:end);
        underscore = strfind(cluster, '_');
        session =  cluster(1:underscore(2)-1);
        directory = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\session_trial_count\bhv_all\',session, '\');
        
        %xBRdatafiles = concat_filenames.(xcluster);
        BRdatafiles = dir(strcat(directory,'\*cinterocdrft*.bhv'));
        BRdataFnames = {BRdatafiles.name};
        eye_info =struct();
        all_codes = [];
        %   all_times = [];
        %      all_analogData =[];
        for n =1:length(BRdataFnames)
            BRdatafile = BRdatafiles(n).name;
            filename   = [directory BRdatafile];
            if exist(filename,'file')
                dot = strfind(BRdatafile, '.');
                fieldn =  BRdatafile(1:dot-1);
                eye_info.(strcat('x',fieldn)) = concatBHV(filename);
                if  ~isempty(eye_info.(strcat('x',fieldn)))
                    all_codes = [all_codes, eye_info.(strcat('x',fieldn)).CodeNumbers];
                    %all_times = [all_times, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeTimes];
                    %all_analogData = [all_analogData,eye_info.(strcat(xBRdatafile,'_bhvfile')).AnalogData];
                    % xBaseline = eye_info.(strcat(xBRdatafile,'_bhvfile')).ScreenXresolution/4/eye_info.(strcat(xBRdatafile,'_bhvfile')).PixelsPerDegree; %since the screen monitor is split in 2 parts with the stereoscope, the center for each eye becomes the center for each side of the stereoscope (half of the half, justifying dividing by 4
                end
            end
        end
        trialCnt = [trialCnt, length(all_codes)];
    catch
        cnt = cnt+1;
        disp(strcat({'missing data ' BRdatafile}))
    end
end

trialCnt = trialCnt(find(trialCnt));

avgTrC = mean(trialCnt);
CITrC = 1.96*std(trialCnt)/sqrt(length(trialCnt));

%% Do the same with nev files
% NEV             = openNEV([filename '.nev'],'noread','overwrite');
% EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;        % we don't know why we have to subtract 128 but we do
% EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                 % in samples 
% EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000); % convert to ms 
% [pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);          % sorts codes, samps or times into trials
% 
% 
% allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_pooled_09292021';
% selected_data = load(strcat(allfilename, '.mat'));
% selected_data =selected_data.selected_data_pooled;
% xFilenames = fieldnames(selected_data);
% 
% cnt =0;
% trialCnt = [];
% for i=1:length(xFilenames)
%     try
%         xcluster = xFilenames{i};
%         cluster = xcluster(2:end);
%         underscore = strfind(cluster, '_');
%         session =  cluster(1:underscore(2)-1);
%         directory = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\session_trial_count\bhv_all\',session, '\');
%         
%         %xBRdatafiles = concat_filenames.(xcluster);
%         BRdatafiles = dir(strcat(directory,'\*cinterocdrft*.bhv'));
%         BRdataFnames = {BRdatafiles.name};
%         eye_info =struct();
%         all_codes = [];
%         %   all_times = [];
%         %      all_analogData =[];
%         for n =1:length(BRdataFnames)
%             BRdatafile = BRdatafiles(n).name;
%             filename   = [directory BRdatafile];
%             if exist(filename,'file')
%                 dot = strfind(BRdatafile, '.');
%                 fieldn =  BRdatafile(1:dot-1);
%                 eye_info.(strcat('x',fieldn)) = concatBHV(filename);
%                 if  ~isempty(eye_info.(strcat('x',fieldn)))
%                     all_codes = [all_codes, eye_info.(strcat('x',fieldn)).CodeNumbers];
%                     %all_times = [all_times, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeTimes];
%                     %all_analogData = [all_analogData,eye_info.(strcat(xBRdatafile,'_bhvfile')).AnalogData];
%                     % xBaseline = eye_info.(strcat(xBRdatafile,'_bhvfile')).ScreenXresolution/4/eye_info.(strcat(xBRdatafile,'_bhvfile')).PixelsPerDegree; %since the screen monitor is split in 2 parts with the stereoscope, the center for each eye becomes the center for each side of the stereoscope (half of the half, justifying dividing by 4
%                 end
%             end
%         end
%         trialCnt = [trialCnt, length(all_codes)];
%     catch
%         cnt = cnt+1;
%         disp(strcat({'missing data ' BRdatafile}))
%     end
% end

%% Count total number of trials per single unit recorded

% %contLims = [0,0.1,0.3,0.5,0.7,1];

sourcedir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\';
funits = dir(strcat(sourcedir,'\*cinterocdrft*.mat'));
funitsNames = {funits.name};
trials = struct();
 tn = [];
for i =1:length(funitsNames)
    suaDat = load(strcat(sourcedir, funitsNames{i}));
    filename = funitsNames{i};
    dot = strfind(funitsNames{i}, '.');
    fname =  filename(1:dot-1);
    trials.(strcat('x',fname)) = suaDat.STIM.trial;
     tn = [tn, length(suaDat.STIM.trial)];
end
tn = tn(tn(:,1) ~= 0,:);
mean(tn)
std(tn)
%count the number of trials in each condition (high monocular, high
%binocular contrast)
sourcedir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\';
funits = dir(strcat(sourcedir,'\*cinterocdrft*.mat'));
funitsNames = {funits.name};

 tnmono = [];
 tnbino = [];
for i =1:length(funitsNames)
    suaDat = load(strcat(sourcedir, funitsNames{i}));
    filename = funitsNames{i};
    monohighcontrast = suaDat.STIM.contrast >=  0.5 & suaDat.STIM.fixedc ==  0; %contrast >=0.5 in DE, contrast = 0 in NDE
    monotn = nnz(find(monohighcontrast));
    tnmono = [tnmono, monotn];
    binohighcontrast = suaDat.STIM.contrast >=  0.5 & suaDat.STIM.fixedc >=  0.5;
    binotn = nnz(find(binohighcontrast));
    tnbino = [tnbino, binotn];
end
tnmono = tnmono(tnmono(:,1) ~= 0,:);
tnbino = tnbino(tnbino(:,1) ~= 0,:);
mean(tnmono)
std(tnmono)
mean(tnbino)
std(tnbino)


%%  Assess: (with the 71 units post power check)
%1) number of days that the 71 single units represent

sourcedir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\';
funits = dir(strcat(sourcedir,'\*cinterocdrft*.mat'));
funitsNames = {funits.name};
dates = [];
I_dates = [];
B_dates = [];
for i =1:length(funitsNames)

    filename = funitsNames{i};
    fname =  filename(1:6);
    dates = [dates , str2double(fname)];
    mfname = filename(1:8);
    if contains(mfname, 'I')
        I_dates = [I_dates, {mfname}];
    elseif contains(mfname, 'B')
        B_dates = [B_dates, {mfname}];
    end
end
days = unique(dates);
Idays = unique(I_dates);
Bdays = unique(B_dates);
%2) number of files per day that they represent

%below code snippet adapted from "C:\Users\daumail\OneDrive -
%Vanderbilt\Documents\loic_code_042021\single_units_analysis\refined_analysis\s_potential_analysis\get_selected_trials_idx.m"
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [newdatadir 'refined_dataset']; 
data_file = load(channelfilename);
f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
channum = 1: length(data_file.new_data);
unitNames =[]; %first step = get units name
 clear i
 for i = channum
     
     filename = [data_file.new_data(i).channel_data.filename, f{2}];
     filename = erase(filename, '.mat');
     unitNames = [unitNames, {filename}];
 end
savename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\analysis\stim_penetrations_10162022';
save(strcat(savename, '.mat'), 'unitNames');
 
% use code adapted from C:\Users\daumail\OneDrive - Vanderbilt\Documents\loic_code_042021\single_units_analysis\refined_analysis\microsaccades\get_concat_filenames.m
%to obtain the filenames that constitue the single units of interest
ns6dir = 'E:\PC_all_data_2020_to_April2021\LGN_data\single_units\s_potentials\data\ns6_selected_units\';
%ssdir
ssdir = 'E:\PC_all_data_2020_to_April2021\LGN_data\single_units\s_potentials\data\kilosorted_files\';
selected_name_list = dir(ssdir);  %ssdir content

%load file with selected trial indexes for a given penetration with
%penetration name
indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\analysis\';
stim_penetrations = load( [indexdir, 'stim_penetrations_10162022']);

SingUnitConcatFiles = struct();
for pnt = 1:length(stim_penetrations.unitNames) 
        penetration_name = char(erase(stim_penetrations.unitNames(pnt), 'matDE50_NDE0'));
        underscore = strfind(penetration_name, '_');
            for ss =3:length(selected_name_list)
                if contains(erase(selected_name_list(ss).name,'_ss.mat'),penetration_name(1:underscore(3)-1))
                    ss_file = load(strcat(ssdir,selected_name_list(ss).name));
                    ss_file_fieldnames = fieldnames(ss_file.ss);
                    cinterocdrft_names = ss_file_fieldnames(contains(ss_file_fieldnames,'cinterocdrft'));
                    %SingUnitConcatFiles.(strcat('x',penetration_name(1:underscore(2)-1))) = cinterocdrft_names;
                    SingUnitConcatFiles.(strcat('x',penetration_name)) = cinterocdrft_names;
                end
            end
end
save(strcat(indexdir,'concat_origfilenames_completenames_10162022.mat'),'-struct', 'SingUnitConcatFiles')

%3) number of units per file


%4) Duration
