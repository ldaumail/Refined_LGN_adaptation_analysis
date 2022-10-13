%This script was developped to assess the number of trials per session
%Loic -10/12/2022


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
            %mkdir(strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul\', selectDate))
            sourcedir = bhvDir;
            fbhv = dir(strcat(sourcedir,'\*cinterocdrft*.bhv'));
            fbhvNames = {fbhv.name};
            destdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\session_trial_count\bhv_all\',selectDate);
            %mkdir(destdir)
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


%% total number of trials per session


% indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\microsaccades_adaptation_analysis\analysis\';
% concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
% newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
% trialsTraces =load([newdatadir 'NoFiltMultiContSUA_06212021']); %neural data +peaklocs + trial numbers obtained with  "BinocularAdaptationTrialSelection.m"
% xfilenames = fieldnames(trialsTraces.NoFiltMultiContSUA);
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
        for fn =1:length(BRdataFnames)
            BRdatafile = BRdatafiles(n).name;
            filename   = [directory BRdatafile];
            if exist(filename,'file')
                dot = strfind(BRdatafile, '.');
                fieldn =  BRdatafile(1:dot-1);
                eye_info.(strcat('x',fieldn)) = concatBHV(filename);
                all_codes = [all_codes, eye_info.(strcat('x',fieldn)).CodeNumbers];
                %all_times = [all_times, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeTimes];
                %all_analogData = [all_analogData,eye_info.(strcat(xBRdatafile,'_bhvfile')).AnalogData];
               % xBaseline = eye_info.(strcat(xBRdatafile,'_bhvfile')).ScreenXresolution/4/eye_info.(strcat(xBRdatafile,'_bhvfile')).PixelsPerDegree; %since the screen monitor is split in 2 parts with the stereoscope, the center for each eye becomes the center for each side of the stereoscope (half of the half, justifying dividing by 4
            end
        end
        trialCnt = [trialCnt, length(all_codes)];
    catch
        cnt = cnt+1;
        disp(strcat({'missing data ' BRdatafile}))
    end
end