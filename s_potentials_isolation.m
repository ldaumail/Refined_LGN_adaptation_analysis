%this script was developped in order to investigate the s-potentials
%of retinal ganglion cells potentially present in the data around the 666Hz
%frequency.
%script developped by Loic Daumail -08/10/2020









%% other strategy to load data

%load packages used to load .ns6 file
npmkdir    = 'C:\Users\daumail\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\daumail\Documents\bootcamp-selected\nbanalysis\'; 
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))
addpath 'C:\Users\daumail\Documents\loic_code'

%ns6 directory 
ns6dir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\ns6_selected_units\';

%ssdir
ssdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\kilosorted_files\';
%ssdir content
selected_name_list = dir(ssdir);  

%load file with selected trial indexes for a given penetration with
%penetration name
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'selected_trials_idx']);

%load selected penetrations files list
%selected_penetrations = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\selected_orig_units_filenames';
%penetrations_names = textscan( fopen(strcat(selected_penetrations, '.txt')), '%s', 'HeaderLines', 1);
%penetrations_names = cell2struct(penetrations_names, 'penetrations_names');

for pnt = 1:length(selected_trials_idx.logicals)
 
    if ~isempty(selected_trials_idx.logicals(pnt).idx)
       
        penetration_name = erase(selected_trials_idx.logicals(pnt).penetration, 'matDE50_NDE0');
        STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);
    
     %1)get the stim onset times for all trials of the given penetration
        STIM_onsets = STIM_file.STIM.photo_on;
     %only keep the selected trials onsets
        selected_STIM_onsets = cell2mat(STIM_onsets(selected_trials_idx.logicals(pnt).idx));
        selected_STIM_onsets = selected_STIM_onsets(1:4:end);
     %contact index/letter
        contact = STIM_file.STIM.chan;
        
     %2) get appropriate ns6 file of penetration 
     %isolate session date/animal for ns6 folder name
        underscore = strfind(penetration_name, '_');
        session =  penetration_name(1:underscore(2)-1);
    
        
      %load ns6 file of penetration of interest  
        %loop through ss files until we find an ss file with same
        %session and penetration as the penetration of interest
        for ss =1:length(selected_name_list)
            if contains(erase(selected_name_list(ss).name,'_ss.mat'),penetration_name(1:underscore(3)-1))
                ss_file = load(strcat(ssdir,selected_name_list(ss).name));
                ss_file_fieldnames = fieldnames(ss_file.ss);

            
                 cinterocdrft_names = ss_file_fieldnames(contains(ss_file_fieldnames,'cinterocdrft'));
                 for cint =1:length(cinterocdrft_names)
                     ns6_filename = char(cinterocdrft_names(cint));            
                     clear ext NS_header banks neural 
                     % Read in NS Header
                     ext          = 'ns6'; 
                     NS_Header    = openNSx(strcat(ns6dir,session,'\',ns6_filename(2:end),'.',ext),'noread');
                    el  = 'eD';
                    % get basic info about recorded data
                    neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},el(2)); % logicals where contact bank name matches electrode of interest
                    N.neural     = sum(neural); % number of neural channels 
                    NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label}; %get labels
                    Fs           = NS_Header.MetaTags.SamplingFreq; % get sampling frequency
                    nyq          = Fs/2; 
                    r            = Fs/1000; 
                    r2           = Fs/15000;

                    % counters
                    clear nct
                    nct = 0;

                    tic 
                    % process data electrode by electrode
                    for e = 1:length(neural)

                        if neural(e) == 1    % why? because neural is a vector of logicals, 1 = contacts we want

                            nct = nct+1;

                            % open data for this channel. 
                            clear NS DAT
                            electrode = sprintf('c:%u',e);
                            NS        = openNSx(strcat(ns6dir,session,'\',ns6_filename(2:end),'.',ext),electrode,'read','uV');
                            DAT       = NS.Data; NS.Data = [];  % this is the whole signal on one channel, 30 kHz!


                            % preallocate data matrices 
                            if nct == 1
                                N.samples = length(DAT); 
                                MUA       = zeros(ceil(N.samples/r),N.neural); % preallocating for downsampled data
                                sLFP       = zeros(ceil(N.samples/r2),N.neural);
                            end

                            % extract the aMUA. 
                              clear hpc hWn bwb bwa hpMUA
                            hpc       = 750; %high pass cutoff
                            hWn       = hpc/nyq;
                            [bwb,bwa] = butter(4,hWn,'high');
                            hpMUA      = filtfilt(bwb,bwa,DAT);  %low pass filter 

                            clear lpc lWn bwb bwa
                            lpc       = 5000; %low pass cutoff
                            lWn       = lpc/nyq;
                            [bwb,bwa] = butter(4,lWn,'low');
                            hpMUA      = abs(filtfilt(bwb,bwa,hpMUA));  %low pass filter 
                            
                            
                             clear lpc lWn bwb bwa lpMUA
                            lpc       = 200; %low pass cutoff
                            lWn       = lpc/nyq;
                            [bwb,bwa] = butter(4,lWn,'low');
                            lpMUA      = filtfilt(bwb,bwa,hpMUA);  %low pass filter 
                           
                            % extract the sLFP 
                            clear hpc hWn bwb bwa hpsLFP
                            hpc       = 250;  %high pass cutoff
                            hWn       = hpc/nyq;
                            [bwb,bwa] = butter(4,hWn,'high');
                            hpsLFP     = filtfilt(bwb,bwa,DAT); %high pass filter

                            % low pass at 750 Hz and rectify 
                            clear lpc lWn bwb bwa 
                            lpc       = 750;  % cutoff
                            lWn       = lpc/nyq;
                            [bwb,bwa] = butter(4,lWn,'low');
                            lpsLFP     = abs(filtfilt(bwb,bwa,hpsLFP)); %low pass filter &rectify

                           

                            % decimate analog MUA (aMUA) to get 1kHz samp freq
                            MUA(:,nct) = decimate(lpMUA,r);
                            
                            % decimate sLFP to get 20kHz samp freq
                            sLFP(:,nct) = decimate(lpsLFP,r2); 
                            
                            clear DAT 

                        end

                    end
                    toc
                    % Warning! THESE DATA ARE IN THE SAME ORDER AS THE BR PINS, NOT THE ORDER OF THE PROBE
                    %As the .ns6 data was retrieved using openNSx() and not getLFP(), we need
                    %to sort the channels ourselves. With getLFP(), sorting is done
                    %automatically
                    % sort data from top of electrode to bottom.

                       % get indices
                    idx = zeros(1,length(NeuralLabels));
                    for i = 1:length(NeuralLabels)

                       Str  = cell2mat(NeuralLabels(i));
                       Key   = 'eD';
                       Str(strfind(Str, '%02d')) = [];

                       Index = strfind(Str, Key);
                       idx(1, i) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);

                    end
                    Srt_sLFP = nan(length(sLFP(:,1)), length(NeuralLabels));
                    Srt_aMUA = nan(length(MUA(:,1)), length(NeuralLabels));


                    Srt_sLFP(:,idx) = sLFP(:,:); 
                    Srt_aMUA(:,idx) = MUA(:,:);
                    sortedLabels = NeuralLabels(idx);
                    % calculate CSD before triggering to trials OR on the trial data BUT not on
                    % the mean LFP. 

                    CSD = mod_iCSD(Srt_sLFP')';  % this function takes LFP in channels x samples so let's transpose LFP and then flip it right back 
                                            % feed in units of microV and get back units of
                                            % nA/mm^3
                    % pad array if you want to keep the matrix the same size on the channel
                    % dimension as the other matrices

                    CSD = padarray(CSD,[0 1],NaN,'replicate');
                     % trigger the neural data to the event codes of interest
                    pre   = -500;
                    post  = 1500; 
                    
                    %need to scale up the triggering onset and offset times
                    %for data sampled af FS = 15kHz, ==> 15 X more than FS
                    %at 1000 Hz
                    pre2  = -7500;
                    post2 = 22500;

                    STIM_file.STIM.aMUA  = trigData(Srt_aMUA,floor(selected_STIM_onsets./30),-pre,post); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units
                    STIM_file.STIM.sCSD  = trigData(CSD,floor(selected_STIM_onsets./2),-pre2,post2); 
                    STIM_file.STIM.sLFP = trigData(Srt_sLFP,floor(selected_STIM_onsets./2),-pre2,post2); 
                    STIM_aMUA = STIM_file.STIM.aMUA;
                    STIM_sLFP = STIM_file.STIM.sLFP;
                    STIM_sCSD = STIM_file.STIM.sCSD;


                    trigdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\trig_data';

                    save(strcat(trigdir,'\',ns6_filename, '_rectified_1000hz_aMUA.mat'), 'STIM_aMUA');
                    save(strcat(trigdir,'\',ns6_filename, '_rectified_15khz_sCSD.mat'), 'STIM_sCSD');
                    save(strcat(trigdir,'\',ns6_filename, '_rectified_15khz_sLFP.mat'), 'STIM_sLFP');
                    contrast = '_domsup50_nondom0';
                    LimePlot(STIM_aMUA, ns6_filename, ns6dir, contrast)
                  %  LimePlot(STIM_sCSD, ns6_filename, ns6dir, contrast)
                   % LimePlot(STIM_sLFP, ns6_filename, ns6dir, contrast)
                   % data = cat(4,STIM_sLFP,  STIM_sCSD, STIM_aMUA);
                    plotdir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\plots\';
                  %  newLimePlot(data, ns6_filename, plotdir, contrast)
                    for n = 1:length(STIM_aMUA(1,1,:))
                   figure();
                   subplot(2,1,1)
                   x = -500:1500;
                   plot(x,STIM_aMUA(:,14,n), 'Color', 'b')
                   legend('aMUA')
                   xlim([-500 1500])
                   
                   subplot(2,1,2)
                   x= -7500:22500;
                   plot(x,STIM_sLFP(:,14,n), 'Color', 'r')
                   legend('s potential')
                   xlim([-7500 22500])
                    end
                    
                 end
            end
        end
    end
end
 

ns6_filename = 'x160602_I_cinterocdrft012';
trigdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\trig_data\';

STIM_aMUA =load(strcat(trigdir, ns6_filename, '_rectified_1000hz_aMUA.mat'));
STIM_sLFP =load(strcat(trigdir, ns6_filename, '_15khz_sLFP.mat'));
STIM_sCSD =load(strcat(trigdir, ns6_filename, '_15khz_sCSD.mat'));

for n = 1:length(STIM_aMUA.STIM_aMUA(1,1,:))
    figure();
    subplot(2,1,1)
    x = -500:1500;
    plot(x,STIM_aMUA.STIM_aMUA(:,14,n), 'Color', 'b')
    legend('aMUA')
    xlim([-500 1500])

    subplot(2,1,2)
    x= -7500:22500;
    plot(x,STIM_sLFP.STIM_sLFP(:,14,n), 'Color', 'r')
    legend('s potential')
    xlim([-7500 22500])
end


%% Compute spectrogram

%Ses = struct();
bs_data = struct();
%channum = 1: length(data_file.clean_origin_data);
%mean_S = nan(1174,38, length(channum));
%mean_S = nan(1174,38);
xabs = -7500:22500;

%filtered_dMUA = nan(length(xabs), length(channum));
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [250 750];

clear i ;
  %for i = 1:length(channum)
     if ~isempty(STIM_sCSD.STIM_sCSD)
data = squeeze(STIM_sLFP.STIM_sLFP(1:30001,14,:));
   bsl = mean(data(1:3000,:));
   
   norm_mean_bs = data(72:end, :) - bsl;
  % namelist1(1,1:length(sprintf('chan_%d',i))) = sprintf('chan_%d',i);
   %bs_data(i).namelist1 = norm_mean_bs;
clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,:) ,movingwin, params); 
 
%namelist2(1,1:length(sprintf('S_%d',i))) = sprintf('S_%d',i);
%Ses(i).namelist2 = S;
mean_S = nanmean(S,3);

tvec = t*1000 -129;

%time_adj = 1:128;
%tvec = cat(2, time_adj , t*1000) ;
%we can also store tvec and f in a struct, but they are all identical
     end
  %end

  figure, 

imagesc(tvec,sort(f),mean_S(:,:)')
%ylim([2 20]); 
set(gca,'ydir','normal')
title({'Mean spectrogram', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus onset(ms)')
    ylabel('Frequency band (Hz)')
    
    
    for tr =1:length(S(1,1,:))
    figure, 

imagesc(tvec,sort(f),S(:,:, tr)')
%ylim([2 20]); 
set(gca,'ydir','normal')
title({'Spectrogram', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus onset(ms)')
    ylabel('Frequency band (Hz)')
    end
    
  %%
%{
%%
%load selected trials idx
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'selected_trials_idx']);

%load packages used to load .ns6 file
npmkdir    = 'C:\Users\daumail\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\daumail\Documents\bootcamp-selected\nbanalysis\'; 
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))


% get the selected dates for ns6 files
ns6dir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\ns6_selected_units\';
selected_name_list = dir(ns6dir);      

%STIMULUS ONSETS
%LOAD STIM FILE OF INTEREST TO GET THE STIMULUS ONSET TIMES OF SELECTED
%TRIALS
% get the session names of the selected STIM files corresponding to the
% selected dates of ns6 files to load the STIM files of interest
selected_penetrations = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\selected_orig_units_filenames';
penetrations_names = textscan( fopen(strcat(selected_penetrations, '.txt')), '%s', 'HeaderLines', 1);
penetrations_names = cell2struct(penetrations_names, 'penetrations_names');
  
  %load the STIM file of interest
    if strfind(penetrations_names.penetrations_names{1,1}, selected_name_list(3).name)  
  stim_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetrations_names.penetrations_names{1,1}]);
    end
     %get the stim onset times for all trials
stim_onsets = stim_file.STIM.photo_on;
     %only keep the selected trials onsets
selected_stim_onsets = cell2mat(stim_onsets(selected_trials_idx.logicals(1).idx));
selected_stim_onsets = selected_stim_onsets(1:4:end);
 
%PENETRATIONS NAMES of NS6 FILES OF INTEREST
%LOAD SS (KILOSORTED) FILE TO GET THE PENETRATIONS OF INTEREST
 %load ss file (we need those files to determine which penetration had
 %a single unit (there could be multiple penetrations per session)
ssdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\kilosorted_files\';
   %load the ss file with penetrations of interest    
ss_file_names = dir(ssdir);
   if strfind(ss_file_names(8).name, selected_name_list(3).name)
 ss_file = load( [ssdir, ss_file_names(8).name]);
   end
%LOAD NS6 FILE OF PENETRATION OF INTETEREST
 %get ns6 penetration file of interest:
ss_file_fieldnames = fieldnames(ss_file.ss);
ns6_filename = ss_file_fieldnames{10};





%% LOAD LFP and analog MUA with NS6 file ( you can follow these steps to load with ns2 as well but beware of sampling freq.)

 %load ns6 file
clear ext NS_header banks neural 
  % Read in NS Header
  ext          = 'ns6'; 
  NS_Header    = openNSx(strcat(ns6dir,selected_name_list(3).name,'\',ns6_filename(2:end),'.',ext),'noread');

%  let's break the whole thing down. 
el  = 'eD';
% get basic info about recorded data
neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},el(2)); % logicals where contact bank name matches electrode of interest
N.neural     = sum(neural); % number of neural channels 
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label}; %get labels
Fs           = NS_Header.MetaTags.SamplingFreq; % get sampling frequency
nyq          = Fs/2; 
r            = Fs/1000; 

% counters
clear nct
nct = 0;

tic 
% process data electrode by electrode
for e = 1:length(neural)
    
    if neural(e) == 1    % why? because neural is a vector of logicals, 1 = contacts we want

        nct = nct+1;
        
        % open data for this channel. 
        clear NS DAT
        electrode = sprintf('c:%u',e);
        NS        = openNSx(strcat(ns6dir,selected_name_list(3).name,'\',ns6_filename(2:end),'.',ext),electrode,'read','uV');
        DAT       = NS.Data; NS.Data = [];  % this is the whole signal on one channel, 30 kHz!
        
        
        % preallocate data matrices 
        if nct == 1
            N.samples = length(DAT); 
            LFP       = zeros(ceil(N.samples/r),N.neural); % preallocating for downsampled data
            MUA       = zeros(ceil(N.samples/r),N.neural);
        end
        
        % extract the LFP. 
        clear lpc lWn bwb bwa lpLFP
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        lpLFP      = filtfilt(bwb,bwa,DAT);  %low pass filter 
        
        % extract the MUA:
        clear hpc hWn bwb bwa hpMUA
        hpc       = 750;  %high pass cutoff
        hWn       = hpc/nyq;
        [bwb,bwa] = butter(4,hWn,'high');
        hpMUA     = filtfilt(bwb,bwa,DAT); %high pass filter
        
        % low pass at 5000 Hz and rectify 
        clear lpc lWn bwb bwa 
        lpc       = 5000;  % cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        hpMUA     = abs(filtfilt(bwb,bwa,hpMUA)); %low pass filter &rectify
        
        % low pass filter at x Hz. 
        clear lpc lWn bwb bwa lpMUA
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low'); 
        lpMUA     = filtfilt(bwb,bwa,hpMUA);  %low pass filter to smooth
        
        
                           
      
        % decimate both LFP and analog MUA (aMUA) to get 1kHz samp freq
        MUA(:,nct) = decimate(lpMUA,r); 
        LFP(:,nct) = decimate(lpLFP,r); 
        
        clear DAT 
        
    end
    
end
toc

%%
% Warning! THESE DATA ARE IN THE SAME ORDER AS THE BR PINS, NOT THE ORDER OF THE PROBE
%As the .ns6 data was retrieved using openNSx() and not getLFP(), we need
%to sort the channels ourselves. With getLFP(), sorting is done
%automatically
% sort data from top of electrode to bottom.

   % get indices
idx = zeros(1,length(NeuralLabels));
for i = 1:length(NeuralLabels)
    
   Str  = cell2mat(NeuralLabels(i));
   Key   = 'eD';
   Str(strfind(Str, '%02d')) = [];
   
   Index = strfind(Str, Key);
   idx(1, i) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
  
end
Srt_MUA = nan(length(MUA(:,1)), length(NeuralLabels));
Srt_LFP = nan(length(MUA(:,1)), length(NeuralLabels));


Srt_MUA(:,idx) = MUA(:,:); 
Srt_LFP(:,idx) = LFP(:,:);
sortedLabels = NeuralLabels(idx);

%% calculate CSD 
% calculate CSD before triggering to trials OR on the trial data BUT not on
% the mean LFP. 

CSD = mod_iCSD(Srt_LFP')';  % this function takes LFP in channels x samples so let's transpose LFP and then flip it right back 
                        % feed in units of microV and get back units of
                        % nA/mm^3
% pad array if you want to keep the matrix the same size on the channel
% dimension as the other matrices

CSD = padarray(CSD,[0 1],NaN,'replicate');


%%
% trigger the neural data to the event codes of interest
pre   = -500;
post  = 1500; 

stim_file.STIM.LFP  = trigData(Srt_LFP,floor(selected_stim_onsets./30),-pre,post); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units
stim_file.STIM.CSD  = trigData(CSD,floor(stim_file.STIM.onsets./30),-pre,post); 
stim_file.STIM.aMUA = trigData(Srt_MUA,floor(selected_stim_onsets./30),-pre,post); 
STIM_LFP = stim_file.STIM.LFP;
STIM_aMUA = stim_file.STIM.aMUA;
STIM_CSD = stim_file.STIM.CSD;

%save(strcat(filename, 'LFP.mat'), 'STIM_LFP');
%save(strcat(filename, 'CSD.mat'), 'STIM_CSD');
%save(strcat(filename, 'aMUA.mat'), 'STIM_aMUA');


%% Plot data
contrast = '_domsup50_nondom0';
LimePlot(STIM_LFP, ns6_filename, ns6dir, contrast)
LimePlot(STIM_CSD, ns6_filename, ns6dir, contrast)
LimePlot(STIM_aMUA, ns6_filename, ns6dir, contrast)

data = cat(4,STIM_LFP, STIM_aMUA, STIM_CSD);
newLimePlot(data, BRdatafile, filename, contrast)
%}