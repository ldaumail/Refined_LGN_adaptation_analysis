% This script was written in order to compute the effect sizes of
% significantly adapting units. 
%This was applied for both suppressed and facilitated units
% This was applied for both spiking activity and the power data
%Written by Loic Daumail edited on 6/16/2020


%% first of all = cohens d for the spiking activity analysis
pvaluesdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected_dunnett.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_raw_data_peaks']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];
sig_inccd = nan(25,3,3);
 sig_adapcd = nan(25,3,3);
 %non_sig_su =nan(25,3);
  
for nc = 1:length(cellclass)
 clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
   for pn =2:4
 all_mean_peaks = nan(4, length(layer_idx));
         cntsigadapt = 0;
         cntsiginc = 0;
         cntnotsig = 0;
  for nunit = 1:length(layer_idx)
       
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
        

            mean_peaks = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
            all_mean_peaks(:,nunit) = mean_peaks;
            if all_mean_peaks(4,nunit) > all_mean_peaks(1,nunit) && pvalues(layer_idx(nunit),3) < .05
                
                cntsiginc= cntsiginc+1;
                %should compute the cohen's d with var or std of the
                %difference between peak1 and peak4
                var1 =var(peakvals.peak_vals(layer_idx(nunit)).peak(1,:)-peakvals.peak_vals(layer_idx(nunit)).peak(pn,:));
                
               % ntrials =length(peakvals.peak_vals(layer_idx(nunit)).peak(1,:));
            
                sig_inccd(cntsiginc,pn-1,nc) = (all_mean_peaks(1,nunit)-all_mean_peaks(pn,nunit))/(sqrt(var1)); 
     
            
            elseif all_mean_peaks(4,nunit) < all_mean_peaks(1,nunit) && pvalues(layer_idx(nunit),3) < .05
                cntsigadapt= cntsigadapt+1;
                var1 =var(peakvals.peak_vals(layer_idx(nunit)).peak(1,:)-peakvals.peak_vals(layer_idx(nunit)).peak(pn,:));
            
                %ntrials =length(peakvals.peak_vals(layer_idx(nunit)).peak(1,:));
                sig_adapcd(cntsigadapt,pn-1,nc) = (all_mean_peaks(1,nunit)-all_mean_peaks(pn,nunit))/(sqrt(var1));
               
      
            end
      end
  end
   end
end

mean_adapcd =nanmean(sig_adapcd,1);
mean_incd =nanmean(sig_inccd,1);


%% plot cohen's D




%% compute the percent of amplitude decrease
 pvaluesdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected_dunnett.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];
sig_ainc = nan(25,3,3);
 sig_aadap = nan(25,3,3);
 %non_sig_su =nan(25,3);
  
for nc = 1:length(cellclass)
 clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
   for pn =2:4
 all_mean_peaks = nan(4, length(layer_idx));
         cntsigadapt = 0;
         cntsiginc = 0;
       
  for nunit = 1:length(layer_idx)
       
      if ~isempty(peakvals.peak_vals(layer_idx(nunit)).peak)
        

            mean_peaks = nanmean(peakvals.peak_vals(layer_idx(nunit)).peak,2);
            all_mean_peaks(:,nunit) = mean_peaks;
            if all_mean_peaks(4,nunit) > all_mean_peaks(1,nunit) && pvalues(layer_idx(nunit),3) < .05
                
                cntsiginc= cntsiginc+1;
               
                sig_ainc(cntsiginc,pn-1,nc) = 100*(all_mean_peaks(1,nunit)-all_mean_peaks(pn,nunit))/all_mean_peaks(1,nunit); 
     
            
            elseif all_mean_peaks(4,nunit) < all_mean_peaks(1,nunit) && pvalues(layer_idx(nunit),3) < .05
                cntsigadapt= cntsigadapt+1;
                sig_aadap(cntsigadapt,pn-1,nc) = 100*(all_mean_peaks(1,nunit)-all_mean_peaks(pn,nunit))/all_mean_peaks(1,nunit);
               
      
            end
      end
  end
   end
end

mean_aadap =nanmean(sig_aadap,1);
mean_ain =nanmean(sig_ainc,1);

%% cohens d for the analysis of the power
 pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
pvalfilename = [pvaluesdir 'roc_results95_stimonset_to1150ms'];
pvalues = load(pvalfilename);

channeldir = '\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
partsvals = load([channeldir '\part1_part2_power_trials']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];
sig_inccd = nan(25,3);
 sig_adapcd = nan(25,3);
 %non_sig_su =nan(25,3);
  
for nc = 1:length(cellclass)
 clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
  
 all_mean_parts = nan(2, length(layer_idx));
         cntsigadapt = 0;
         cntsiginc = 0;
         cntnotsig = 0;
  for nunit = 1:length(layer_idx)
       
      if ~isempty(partsvals.parts(layer_idx(nunit)).part1)
        

             mean_part1 = nanmean(nanmean(partsvals.parts(layer_idx(nunit)).part1,1),3);
             mean_part2 = nanmean(nanmean(partsvals.parts(layer_idx(nunit)).part2,1),3);
            all_mean_parts(:,nunit) = [mean_part1, mean_part2];
            if all_mean_parts(2,nunit) > all_mean_parts(1,nunit) && pvalues.all_sigs95(layer_idx(nunit)) == 1
                
                cntsiginc= cntsiginc+1;
                var1 =var(nanmean(partsvals.parts(layer_idx(nunit)).part1,1)-nanmean(partsvals.parts(layer_idx(nunit)).part2,1),[],3);
                
                %ntrials =length(squeeze(nanmean(partsvals.parts(layer_idx(nunit)).part1,1)));
                sig_inccd(cntsiginc,nc) = (all_mean_parts(1,nunit)-all_mean_parts(2,nunit))/(sqrt(var1)); 
     
            
            elseif all_mean_parts(2,nunit) < all_mean_parts(1,nunit) && pvalues.all_sigs95(layer_idx(nunit)) ==1
                cntsigadapt= cntsigadapt+1;
                var1 =var(nanmean(partsvals.parts(layer_idx(nunit)).part1,1)-nanmean(partsvals.parts(layer_idx(nunit)).part2,1),[],3);
               % ntrials =length(squeeze(nanmean(partsvals.parts(layer_idx(nunit)).part1,1)));
                sig_adapcd(cntsigadapt,nc) = (all_mean_parts(1,nunit)-all_mean_parts(2,nunit))/(sqrt(var1)); 
               
      
            end
      end
  end
   
end

mean_adapcd =nanmean(sig_adapcd,1);
mean_incd =nanmean(sig_inccd,1);

%% percent change fore the analysis on the power

%% cohens d for the analysis of the power
 pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\all_units\';
pvalfilename = [pvaluesdir 'roc_results95_stimonset_to1150ms'];
pvalues = load(pvalfilename);

channeldir = '\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
partsvals = load([channeldir '\part1_part2_power_trials']);

layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];

cellclass = [ 'M', 'P', 'K'];
sig_inccd = nan(25,3);
 sig_adapcd = nan(25,3);
 %non_sig_su =nan(25,3);
  
for nc = 1:length(cellclass)
 clear layer_idx
layer_idx = find(strcmp(layer, cellclass(nc)));
  
 all_mean_parts = nan(2, length(layer_idx));
         cntsigadapt = 0;
         cntsiginc = 0;
         cntnotsig = 0;
  for nunit = 1:length(layer_idx)
       
      if ~isempty(partsvals.parts(layer_idx(nunit)).part1)
        

             mean_part1 = nanmean(nanmean(partsvals.parts(layer_idx(nunit)).part1,1),3);
             mean_part2 = nanmean(nanmean(partsvals.parts(layer_idx(nunit)).part2,1),3);
            all_mean_parts(:,nunit) = [mean_part1, mean_part2];
            if all_mean_parts(2,nunit) > all_mean_parts(1,nunit) && pvalues.all_sigs95(layer_idx(nunit)) == 1
                
                cntsiginc= cntsiginc+1;
                
                sig_inccd(cntsiginc,nc) = 100*(all_mean_parts(1,nunit)-all_mean_parts(2,nunit))/(all_mean_parts(1,nunit)); 
     
            
            elseif all_mean_parts(2,nunit) < all_mean_parts(1,nunit) && pvalues.all_sigs95(layer_idx(nunit)) ==1
                cntsigadapt= cntsigadapt+1;
               
                sig_adapcd(cntsigadapt,nc) = 100*(all_mean_parts(1,nunit)-all_mean_parts(2,nunit))/(all_mean_parts(1,nunit)); 
               
      
            end
      end
  end
   
end

mean_adapcd =nanmean(sig_adapcd,1);
mean_incd =nanmean(sig_inccd,1);
