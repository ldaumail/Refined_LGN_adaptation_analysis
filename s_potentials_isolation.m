
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
channelfilename = [newdatadir 'clean_origin_sup_50']; 
data_file = load(channelfilename);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

%%
for i =1:length(layer)
    if ~isempty(data_file.clean_origin_data(i).unit) && ~isempty(layer(i))
   
   raw_unit = data_file.new_data(i).channel_data.sdftr_chan
    end
end

%%

if strcmp(getenv('USER'),'maierav')
    npmkdir    = '/Users/alex 1/Desktop/LAB/Loic/NPMK-master/'; 
    nbanadir   = '/Users/alex 1/Desktop/LAB/Loic/nbanalysis/'; 
 
    directory  = '/Users/alex 1/Desktop/LAB/LoicLGNinfo_4LD/';
    BRdatafile = '190119_B_cinterocdrft002';
else
    npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
    nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 
 
    directory  = 'C:\Users\maier\Documents\LGNinfo_4LD-20190826T172747Z-001\LGNinfo_4LD\';
    BRdatafile = '190119_B_cinterocdrft002';
end

npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\LGN_190124_B_cinterocdrft_data\';
BRdatafile = '190326_B_cinterocdrft001';
filename   = [directory BRdatafile]; 

addpath(genpath(directory))
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

%% STEP ONE: LOAD STIMULUS CONDITIONS with text file 

patterns   = {'rforidrft','rfsfdrft','posdisparitydrft','disparitydrft','cinterocdrft','coneinterocdrft','conedrft', ...
                'colorflicker','bwflicker','rfori','rfsize','cinteroc','color','rfsf','mcosinteroc','dotmapping'}; 

for p = 1:length(patterns)

   pattern      = patterns{p}; 
  
   if any(strfind(BRdatafile,pattern))
       startlog = strfind(BRdatafile,pattern); 
       if ~isequal(BRdatafile(startlog:end-3),pattern),continue
       else
       match    = patterns{p}; 
       end
   end
   
end

if isequal(match,'dotmapping')
ext  = '.gDotsXY_di';
else
ext  = ['.g' upper(match) 'Grating_di']; 
end

if contains(ext,'DRFT')
      grating     = readgDRFTGrating([filename ext]); % from nbanalysis (or even MLAnalysisOnline--might be out of date)
elseif contains(ext,'Dots')
      grating     = readgDotsXY([filename ext]);
else
      grating     = readgGrating([filename ext]);
end

%% STEP TWO: LOAD EVENT TIMES/CODES

NEV             = openNEV([filename '.nev'],'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;        % we don't know why we have to subtract 128 but we do
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                 % in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000); % convert to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);          % sorts codes, samps or times into trials


% So far all of these data are from EVERY trial, including trials where
% animal breaks fixation. Lets get rid of those and make a new structure
% with the grating info and the stimulus onsets 

STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); % this is in nbanalysis. definitely double check it before you use it. 




%% STEP THREE: LOAD NEURAL DATA

%% LOAD LFP with NS2 file

clear ext 

ext = 'ns2';
el  = 'eC';
lfp = getLFP(filename,ext,el,'ascending');        % if you use this function you need to know the bank + sort direction 
                                                  % this function can be found in MLAnalysisOnline or nbanalysis
                                                  % if your sort direction input is correct, the data come out in order
                                                  % from top--> bottom
%% here's a quick way to get  bank info + sort direction if you don't have it in an Excel sheet or notes: 

% Read in NS Header

NS_Header    = openNSx(strcat(filename,'.',ext),'noread');
banks        = unique({NS_Header.ElectrodesInfo.ConnectorBank}); banks(ismember(banks,'E')) = []; % bank E is BNC cable inputs = eye tracking

for b = 1:length(banks)
    clear neural label 
    neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},banks{b}); 
    firstlabel   = cell2mat({NS_Header.ElectrodesInfo(find(neural,1,'first')).Label}); 
    if str2double(firstlabel(3:4)) < 2
        sortdirection = 'ascending'; 
    else
        sortdirection = 'descending'; 
    end
end

%}                                                                                   
%% LOAD LFP and analog MUA with NS6 file ( you can follow these steps to load with ns2 as well but beware of sampling freq.)
%  let's break the whole thing down. 

clear ext NS_header banks neural 

% Read in NS Header
ext          = 'ns6'; 
NS_Header    = openNSx(strcat(filename,'.',ext),'noread');

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
        NS        = openNSx(strcat(filename,'.',ext),electrode,'read','uV');
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


% Warning! THESE DATA ARE IN THE SAME ORDER AS THE BR PINS, NOT THE ORDER OF THE PROBE
%As the .ns6 data was retrieved using openNSx() and not getLFP(), we need
%to sort the channels ourselves. With getLFP(), sorting is done
%automatically
% sort data from top of electrode to bottom.

   % get indices
idx = zeros(1,length(NeuralLabels));
for i = 1:length(NeuralLabels)
    
   Str  = cell2mat(NeuralLabels(i));
   Key   = 'eC';
   Str(strfind(Str, '%02d')) = [];
   
   Index = strfind(Str, Key);
   idx(1, i) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
  
end
% Sort data thank to indices created above
% better  code here, but need to rename variables
%{
switch sortdirection
    case 'ascending'
        MUA = MUA(:,idx);
        LFP = LFP(:,idx);
        sortedLabels = NeuralLabels(idx); 
    case 'descending'
        MUA = MUA(:,flipud(idx));
        LFP = LFP(:,flipud(idx));
        sortedLabels = NeuralLabels(fliplr(idx)); 
end
%}


Srt_MUA = nan(length(MUA(:,1)), length(NeuralLabels));
Srt_LFP = nan(length(MUA(:,1)), length(NeuralLabels));


Srt_MUA(:,idx) = MUA(:,:); 
Srt_LFP(:,idx) = LFP(:,:);
sortedLabels = NeuralLabels(idx);