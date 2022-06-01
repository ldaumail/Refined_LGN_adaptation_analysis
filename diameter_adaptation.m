%This script was developped to assess the influence of grating diameter on
%adaptation
%Loic Daumail 05/06/2022

%With updated functions
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;


unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);


cellClass = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
cellClass([1,46,55]) = [];
 allContLevels =0;
 alltilts =[];
 for i =1:71
     if ~isempty(filenames(i))
    contLevels = unique(unitsData.new_data(i).channel_data.fixedc);
    allContLevels = unique([allContLevels; contLevels]);
    tilts = find(unique(unitsData.new_data(i).channel_data.tilt));
    alltilts = [alltilts; tilts];
     end
 end
%select trials with at least 4 peak values, of convenient quality. Keep peak locations and trial responses: 
[peakLocs, NoFiltMultiContSUA] = peakLocsTrialSelection(unitsData, filenames);

%Store Peaks and peak-triggered trials
[peak_vals, peak_aligned_trials] = peaksAndPeakTrigResps(peakLocs, NoFiltMultiContSUA);

%Get grating diameter
grating_diameters =[];
 for i =1:length(filenames)
     if ~isempty(filenames{i})
    diameter = unique(unitsData.new_data(i).channel_data.diameter);
    grating_diameters = [grating_diameters; diameter];
     end
 end

 %assess adaptation effect in monocular condition
 dp1p4 = struct('bin',[]);
 normdp1p4 = struct('bin',[]);
 midx = struct('bin',[]);
 meanpk = struct('bin',[]);
 filenames = fieldnames(peak_vals);
 contLims = [0,0.1,0.3,0.5,0.7,1];
 clear i
 for i = 1:length(fieldnames(peak_vals))
     for  n =1:length(contLims)
         if n ==1 || n ==6
         binNb = sprintf('bin%d', n);
          if isfield(peak_vals.(filenames{i}), binNb)
         dp1p4.bin.(binNb)(i) = -mean(peak_vals.(filenames{i}).(binNb)(1,:) - peak_vals.(filenames{i}).(binNb)(4,:));
         normdp1p4.bin.(binNb)(i) = -mean(peak_vals.(filenames{i}).(binNb)(1,:) - peak_vals.(filenames{i}).(binNb)(4,:))/mean(peak_vals.(filenames{i}).(binNb)(1,:));
         midx.bin.(binNb)(i) = -2*(mean(peak_vals.(filenames{i}).(binNb)(1,:) - peak_vals.(filenames{i}).(binNb)(4,:)))/(mean(peak_vals.(filenames{i}).(binNb)(1,:) + peak_vals.(filenames{i}).(binNb)(4,:)));
         meanpk.(binNb) = mean(peak_vals.(filenames{i}).(binNb),2);
          end
         end
     end
 end
 
 
 figure()
 for  n =1:length(contLims)
     if n ==1 || n ==6
         binNb = sprintf('bin%d', n);
         if isfield(midx.bin, binNb)
             plot(grating_diameters, midx.bin.(binNb), 'o')
             hold on
         end
     end
 end
  

%% Plot adaptation as a function of diameter with gramm
x=squeeze(grating_diameters');
y=squeeze(normdp1p4.bin.(binNb));
figure('Position',[100 100 800 500]);
clear g
g=gramm('x',x,'y',y);
g.stat_glm();
g.axe_property('xlim',[0 9]);
g.set_names('x','Diameter','y','(P4-P1)/P4','column','');
g(1,1).geom_point();

g.set_title('Adaptation index as a function of grating diameter');
g.draw();

%% Plot adaptation as a function of diameter with polyfit

nlines = 7;
cmap =flip(cbrewer2('Blues', nlines));
figure();
x=squeeze(grating_diameters');
binNb = sprintf('bin%d', 1); %plot monocular distribution
y=squeeze(midx.bin.(binNb));
coeffs = polyfit(x,y,1);
f = polyval(coeffs,x);
plot(x, y,'o',x, f,'-','Color',cmap(4,:),'MarkerSize',2, 'MarkerFaceColor',cmap(4,:),'linewidth',2)
xlim([0 9])
text(max(x)/1.3,max(y)/20, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)))
hold on
binNb = sprintf('bin%d', 6); %plot binocular condition
y1=squeeze(midx.bin.(binNb));
coeffs1 = polyfit(x,y1,1);
f1 = polyval(coeffs1,x);
plot(x, y1,'o',x, f1,'-','Color',cmap(1,:),'MarkerSize',2, 'MarkerFaceColor',cmap(1,:),'linewidth',2)
text(max(x)/1.3,max(y1)/2, sprintf('y = %.2f + %.2f*x', round(coeffs1(2),2), round(coeffs1(1),2)))

colormap(cmap);
hold on
set(gca,'box','off')
set(gca, 'linewidth', 2)
title('Adaptation index as a function of grating diameter');
xlabel('Diameter (deg)')
ylabel('Adaptation strength')
legend('','Monocular','','Binocular')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\adaptationidx_vs_diameter');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

 %% test significance of slope versus 0

    xtest = x;
    ytest = y1;
    linreg = fitlm(xtest,ytest);
    [linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linreg);
    Rsquared(1,1) = linreg.Rsquared.Ordinary;

    
   %% Add trends of adapting single units
 
   %get pvalues from lmer results with Dunnett correction
   pvalues = dlmread('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\lmer_results_orig_03032020_corrected_dunnett.csv', ',', 1,1);
   
   %clear out nans
   pk1pk4pval = pvalues(~isnan(pvalues(:,3)),3);
   
   %plot results
   
nlines = 7;
cmap =flip(cbrewer2('Greys', nlines));
cmap2 = flip(cbrewer2('Blues', nlines));
cmap3 = flip(cbrewer2('Reds', nlines));
figure('Renderer', 'painters', 'Position', [300 300 1200 1000]);
x=squeeze(grating_diameters');

%plot monocular distribution
y=squeeze(midx.bin.bin1);
coeffs = polyfit(x,y,1);
f = polyval(coeffs,x);
plot(x, y,'o',x, f,'-','Color',cmap(4,:),'MarkerSize',5, 'MarkerFaceColor',cmap(4,:),'linewidth',2)
xlim([0 9])
text(max(x)/1.3,max(y)/20, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)))
hold on

%plot binocular condition
y1=squeeze(midx.bin.bin6);
coeffs1 = polyfit(x,y1,1);
f1 = polyval(coeffs1,x);
plot(x, y1,'o',x, f1,'-','Color',cmap(1,:),'MarkerSize',5, 'MarkerFaceColor',cmap(1,:),'linewidth',2)
text(max(x)/1.3,max(y1)/2, sprintf('y = %.2f + %.2f*x', round(coeffs1(2),2), round(coeffs1(1),2)))
hold on

%Monocular
%plot linreg of y vals with pvaln<0.05 that are suppressed (monocular)
ysig = y(pk1pk4pval<0.05 & midx.bin.bin1'<0);
xsig = x(pk1pk4pval<0.05 & midx.bin.bin1'<0);
coeffsig = polyfit(xsig,ysig,1);
f = polyval(coeffsig,xsig);
plot(xsig, ysig,'o',xsig, f,'-','Color',cmap2(4,:),'MarkerSize',2, 'MarkerFaceColor',cmap2(4,:),'linewidth',2)
xlim([0 9])
text(max(xsig)/1.3,max(ysig)/20, sprintf('ysig = %.2f + %.2f*x', round(coeffsig(2),2), round(coeffsig(1),2)))
hold on  

%plot linreg of y vals with pvaln<0.05 that are facilitated (monocular)
yfac = y(pk1pk4pval<0.05 & midx.bin.bin1'>0);
xfac = x(pk1pk4pval<0.05 & midx.bin.bin1'>0);
coeffac = polyfit(xfac,yfac,1);
f = polyval(coeffac,xfac);
plot(xfac, yfac,'o',xfac, f,'-','Color',cmap3(4,:),'MarkerSize',2, 'MarkerFaceColor',cmap3(4,:),'linewidth',2)
xlim([0 9])
text(max(xfac)/1.3,max(yfac)/20, sprintf('yfac = %.2f + %.2f*x', round(coeffac(2),2), round(coeffac(1),2)))
hold on  

%Binocular
%keep y vals with pvaln<0.05 that are suppressed (binocular)
y1sig = y1(pk1pk4pval<0.05 & midx.bin.bin1'<0); %keep using bin1 since the test results were performed on monocular stimulation condition
coeff1sig = polyfit(xsig,y1sig,1);
f = polyval(coeff1sig,xsig);
plot(xsig, y1sig,'o',xsig, f,'-','Color',cmap2(1,:),'MarkerSize',2, 'MarkerFaceColor',cmap2(1,:),'linewidth',2)
xlim([0 9])
text(max(xsig)/1.3,max(y1sig)/20, sprintf('ysigbino = %.2f + %.2f*x', round(coeff1sig(2),2), round(coeff1sig(1),2)))
hold on  

%keep y vals with pvaln<0.05 that are facilitated (binocular)
y1fac = y1(pk1pk4pval<0.05 & midx.bin.bin1'>0);
coeff1fac = polyfit(xfac,y1fac,1);
f = polyval(coeff1fac,xfac);
plot(xfac, y1fac,'o',xfac, f,'-','Color',cmap3(1,:),'MarkerSize',2, 'MarkerFaceColor',cmap3(1,:),'linewidth',2)
xlim([0 9])
text(max(xfac)/1.3,max(y1fac)/20, sprintf('yfacbino = %.2f + %.2f*x', round(coeff1fac(2),2), round(coeff1fac(1),2)))
hold on  


colormap(cmap);
hold on
set(gca,'box','off')
set(gca, 'linewidth', 2)
title('Adaptation index as a function of grating diameter');
xlabel('Diameter (deg)')
ylabel('Adaptation strength')
legend('','Monocular','','Binocular','','Monocular Suppressed','','Monocular Facilitated','','Binocular Suppressed','','Binocular Facilitated')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\adaptationidx_vs_diameter_signif');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));
   
   %% test significance of slope versus 0 for suppressed units

    xtest = xfac;
    ytest = yfac;
    linreg = fitlm(xtest,ytest);
    [linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linreg) %r =numerator degrees of freedom, linreg.DFE = denominator degrees of freedom.
    linreg.DFE
    
    
    
    Rsquared(1,1) = linreg.Rsquared.Ordinary

    
    %% Redo the whole analysis with previously isolated peak values
    
    
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']; 
peak_aligned_trials = load(channelfilename);

%Get grating diameter
unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);

filenames = fieldnames(peak_aligned_trials.peak_aligned_trials);
grating_diameters =[];
 for i =1:length(filenames)
     if ~isempty(filenames{i})
    diameter = unique(unitsData.new_data(i).channel_data.diameter);
    grating_diameters = [grating_diameters; diameter];
     end
 end

%assess adaptation effect in monocular condition
 dp1p4 = struct('bin',[]);
 normdp1p4 = struct('bin',[]);
 midx = struct('bin',[]);
 meanpk = struct('bin',[]);
 %filenames = fieldnames(peak_vals);

 bins = [1,6];
 contLims = [0,0.1,0.3,0.5,0.7,1];
 clear i
 for i = 1:length(filenames)
     for  n =1:length(contLims)
          if length(fieldnames(peak_aligned_trials.peak_aligned_trials.(filenames{i}).origin)) == 2
              
              for b =bins
                  binNb = sprintf('bin%d', b);
                  peak_vals = [];
                  %  isfield(peak_aligned_trials.peak_aligned_trials.(filenames{i}).origin, binNb)
                  for p =1:4
                      pk = sprintf('pk%d', p);
                      peak_vals = [peak_vals; max(peak_aligned_trials.peak_aligned_trials.(filenames{i}).origin.(binNb).(pk),[],1)];
                  end
                  
                  dp1p4.bin.(binNb)(i) = -mean(peak_vals(1,:) - peak_vals(4,:)); %Pk1 - Pk4
                  normdp1p4.bin.(binNb)(i) = -mean(peak_vals(1,:) - peak_vals(4,:))/mean(peak_vals(1,:)); %(Pk1 - Pk4)/Pk1
                  midx.bin.(binNb)(i) = mean(-2*(peak_vals(1,:) - peak_vals(4,:))/(peak_vals(1,:) + peak_vals(4,:))); %mean(2*(Pk1 - Pk4)/(Pk1 + Pk4))
                  meanpk.(binNb) = mean(peak_vals,2); %mean across all peaks
              end
          end
     end
 end
 
 
 figure()
 for  n =1:length(contLims)
     if n ==1 || n ==6
         binNb = sprintf('bin%d', n);
         if isfield(midx.bin, binNb)
             plot(grating_diameters, midx.bin.(binNb), 'o')
             hold on
         end
     end
 end
  
 %% Plot trends 
 
   %get pvalues from lmer results with Dunnett correction
   pvalues = dlmread('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\lmer_results_orig_03032020_corrected_dunnett.csv', ',', 1,1);
   
   %clear out nans
   pk1pk4pval = pvalues(~isnan(pvalues(:,3)),3);
   
   pk1pk4pval = pk1pk4pval(midx.bin.(binNb) ~= 0); % if midx.bin.(binNb) = 0, this means there was no difference between Pk1 and Pk4, this is impossible, this means binN was empty 
   %plot results
   
nlines = 7;
cmap =flip(cbrewer2('Greys', nlines));
cmap2 = flip(cbrewer2('Blues', nlines));
cmap3 = flip(cbrewer2('Reds', nlines));
figure('Renderer', 'painters', 'Position', [300 300 1200 1000]);
x=squeeze(grating_diameters(midx.bin.bin1 ~= 0)');

%plot monocular distribution
y=squeeze(midx.bin.bin1(midx.bin.bin1 ~= 0));
coeffs = polyfit(x,y,1);
f = polyval(coeffs,x);
plot(x, y,'o',x, f,'-','Color',cmap(4,:),'MarkerSize',5, 'MarkerFaceColor',cmap(4,:),'linewidth',2)
xlim([0 9])
text(max(x)/1.3,max(y)/20, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)))
hold on

%plot binocular condition
y1=squeeze(midx.bin.bin6(midx.bin.bin1 ~= 0));
coeffs1 = polyfit(x,y1,1);
f1 = polyval(coeffs1,x);
plot(x, y1,'o',x, f1,'-','Color',cmap(1,:),'MarkerSize',5, 'MarkerFaceColor',cmap(1,:),'linewidth',2)
text(max(x)/1.3,max(y1)/2, sprintf('y = %.2f + %.2f*x', round(coeffs1(2),2), round(coeffs1(1),2)))
hold on

%Monocular
%plot linreg of y vals with pvaln<0.05 that are suppressed (monocular)
ysig = y(pk1pk4pval<0.05 & midx.bin.bin1(midx.bin.bin1 ~= 0)'<0);
xsig = x(pk1pk4pval<0.05 & midx.bin.bin1(midx.bin.bin1 ~= 0)'<0);
coeffsig = polyfit(xsig,ysig,1);
f = polyval(coeffsig,xsig);
plot(xsig, ysig,'o',xsig, f,'-','Color',cmap2(4,:),'MarkerSize',2, 'MarkerFaceColor',cmap2(4,:),'linewidth',2)
xlim([0 9])
text(max(xsig)/1.3,max(ysig)/20, sprintf('ysig = %.2f + %.2f*x', round(coeffsig(2),2), round(coeffsig(1),2)))
hold on  

%plot linreg of y vals with pvaln<0.05 that are facilitated (monocular)
yfac = y(pk1pk4pval<0.05 & midx.bin.bin1(midx.bin.bin1 ~= 0)'>0);
xfac = x(pk1pk4pval<0.05 & midx.bin.bin1(midx.bin.bin1 ~= 0)'>0);
coeffac = polyfit(xfac,yfac,1);
f = polyval(coeffac,xfac);
plot(xfac, yfac,'o',xfac, f,'-','Color',cmap3(4,:),'MarkerSize',2, 'MarkerFaceColor',cmap3(4,:),'linewidth',2)
xlim([0 9])
text(max(xfac)/1.3,max(yfac)/20, sprintf('yfac = %.2f + %.2f*x', round(coeffac(2),2), round(coeffac(1),2)))
hold on  

%Binocular
%keep y vals with pvaln<0.05 that are suppressed (binocular)
y1sig = y1(pk1pk4pval<0.05 & midx.bin.bin1(midx.bin.bin1 ~= 0)'<0); %keep using bin1 since the test results were performed on monocular stimulation condition
coeff1sig = polyfit(xsig,y1sig,1);
f = polyval(coeff1sig,xsig);
plot(xsig, y1sig,'o',xsig, f,'-','Color',cmap2(1,:),'MarkerSize',2, 'MarkerFaceColor',cmap2(1,:),'linewidth',2)
xlim([0 9])
text(max(xsig)/1.3,max(y1sig)/20, sprintf('ysigbino = %.2f + %.2f*x', round(coeff1sig(2),2), round(coeff1sig(1),2)))
hold on  

%keep y vals with pvaln<0.05 that are facilitated (binocular)
y1fac = y1(pk1pk4pval<0.05 & midx.bin.bin1(midx.bin.bin1 ~= 0)'>0);
coeff1fac = polyfit(xfac,y1fac,1);
f = polyval(coeff1fac,xfac);
plot(xfac, y1fac,'o',xfac, f,'-','Color',cmap3(1,:),'MarkerSize',2, 'MarkerFaceColor',cmap3(1,:),'linewidth',2)
xlim([0 9])
text(max(xfac)/1.3,max(y1fac)/20, sprintf('yfacbino = %.2f + %.2f*x', round(coeff1fac(2),2), round(coeff1fac(1),2)))
hold on  


colormap(cmap);
hold on
set(gca,'box','off')
set(gca, 'linewidth', 2)
title('Adaptation index as a function of grating diameter');
xlabel('Diameter (deg)')
ylabel('Adaptation strength')
legend('','Monocular','','Binocular','','Monocular Suppressed','','Monocular Facilitated','','Binocular Suppressed','','Binocular Facilitated')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\adaptationidx_vs_diameter_signif_old');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));
   
   %% test significance of slope versus 0 for suppressed units

    xtest = xfac;
    ytest = y1fac;
    linreg = fitlm(xtest,ytest);
    [linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linreg) %r =numerator degrees of freedom, linreg.DFE = denominator degrees of freedom.
    linreg.DFE
    length(xtest)
    
    
    Rsquared(1,1) = linreg.Rsquared.Ordinary

 