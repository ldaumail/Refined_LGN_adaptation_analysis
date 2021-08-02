
%code developped to obtain stimulus size statistics
%of selected single units
%Loic Daumail 8/2/2021

%get filenames where the data is located
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\s_potentials_analysis\analysis\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName; %get filenames of selected single units

unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir); %get single units data

contLims = [0,0.1,0.3,0.5,0.7,1];
channum = 1: length(filenames);
stimDiameter = struct();

for i = channum
    filename = filenames{i};
    if ~isempty(filenames{i})
        SUA_data = unitsData.new_data(i);
        blankcontrast = SUA_data.channel_data.contrast ==  0 & SUA_data.channel_data.fixedc ==  0; %get logical indices of trials with 0 contrast in both eyes
        filename = erase(sprintf('x%s',char(filename)),'.mat');
        for n = 1:length(contLims)
            if n == 1
                contrastBin = SUA_data.channel_data.contrast >=  0.5 & SUA_data.channel_data.fixedc ==  0; %trials indices with 0 contrast in NDE, and contrast >0.5 in DE
            else
                if n>1
                    contrastBin = (SUA_data.channel_data.fixedc >  contLims(n-1) & SUA_data.channel_data.fixedc <= contLims(n))& SUA_data.channel_data.contrast >=  0.5;
                end
            end
            stimDiameter.(filename).(sprintf('bin%d', n)) = SUA_data.channel_data.diameter(contrastBin);
        end
    end
end


%stim size stats by monkey

selectFilenames = fieldnames(stimDiameter);
%diamIBin1 = [];
%diamIBin2 = [];
allDiamI = [];
%diamBBin1 = [];
%diamBBin2 = [];
allDiamB = [];
cntI = 0;
cntB = 0;

for i =1:length(fieldnames(stimDiameter))
    filename = selectFilenames{i};
    if contains(filename, '_I_')
        cntI =cntI+1;
       % diamIBin1 = [diamIBin1; stimDiameter.(filename).bin1];
        %diamIBin2 = [diamIBin2; stimDiameter.(filename).bin6];
        allDiamI = [allDiamI; stimDiameter.(filename).bin1;stimDiameter.(filename).bin6];
    else
        if contains(filename, '_B_')
             cntB =cntB+1;
           % diamBBin1 = [diamBBin1; stimDiameter.(filename).bin1];
           % diamBBin2 = [diamBBin2; stimDiameter.(filename).bin6];
            allDiamB = [allDiamB; stimDiameter.(filename).bin1;stimDiameter.(filename).bin6];

        end
    end
end
%{
 meanDiamI_bin1 = mean(diamIBin1);
 meanDiamI_bin2 = mean(diamIBin2);  
  
 stdDiamI_bin1 = std(diamIBin1);
 stdDiamI_bin2 = std(diamIBin2);
 
 minDiamI_bin1 = min(diamIBin1);
 maxDiamI_bin1 = max(diamIBin1);
 
 meanDiamB_bin1 = mean(diamBBin1);  
 meanDiamB_bin2 = mean(diamBBin2);  
  
 stdDiamB_bin1 = std(diamBBin1);
 stdDiamB_bin2 = std(diamBBin2);
 
 minDiamB_bin1 = min(diamBBin1);
 maxDiamB_bin1 = max(diamBBin1);
 %}
 %%all stims (bin1+bin2) stats
 meanDiamI = mean(allDiamI);  
 medianDiamI = median(allDiamI); 
 stdDiamI = std(allDiamI);
 minDiamI = min(allDiamI);
 maxDiamI = max(allDiamI);
 
 meanDiamB = mean(allDiamB); 
 medianDiamB = median(allDiamB); 
 stdDiamB = std(allDiamB);
 minDiamB = min(allDiamB);
 maxDiamB = max(allDiamB);
 
 
 