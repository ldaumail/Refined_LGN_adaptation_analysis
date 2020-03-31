
 pvaluesdir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
pvalfilename = [pvaluesdir 'lmer_results_orig_03032020_corrected.csv'];
pvalues = dlmread(pvalfilename, ',', 1,1);

channeldir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\';
peakvals = load([channeldir 'all_data_peaks']);

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
            if all_mean_peaks(pn,nunit) > all_mean_peaks(1,nunit) && pvalues(layer_idx(nunit),pn) < .05
                
                cntsiginc= cntsiginc+1;
                var1 =var(peakvals.peak_vals(layer_idx(nunit)).peak(1,:));
                var2 =var(peakvals.peak_vals(layer_idx(nunit)).peak(pn,:));
                ntrials =length(peakvals.peak_vals(layer_idx(nunit)).peak(1,:));
                sig_inccd(cntsiginc,pn-1,nc) = (all_mean_peaks(1,nunit)-all_mean_peaks(pn,nunit))/(sqrt((var1+var2)/ntrials)); 
     
            
            elseif all_mean_peaks(pn,nunit) < all_mean_peaks(1,nunit) && pvalues(layer_idx(nunit),pn) < .05
                cntsigadapt= cntsigadapt+1;
                var1 =var(peakvals.peak_vals(layer_idx(nunit)).peak(1,:));
                var2 =var(peakvals.peak_vals(layer_idx(nunit)).peak(pn,:));
                ntrials =length(peakvals.peak_vals(layer_idx(nunit)).peak(1,:));
                sig_adapcd(cntsigadapt,pn-1,nc) = (all_mean_peaks(1,nunit)-all_mean_peaks(pn,nunit))/(sqrt((var1+var2)/ntrials));
               
      
            end
      end
  end
   end
end

mean_adapcd =nanmean(sig_adapcd,1);
mean_incd =nanmean(sig_inccd,1);

