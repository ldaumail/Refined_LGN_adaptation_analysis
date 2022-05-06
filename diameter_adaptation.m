%This script was developped to assess the influence of grating diameter on
%adaptation
%Loic Daumail 05/06/2022
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
y=squeeze(normdp1p4.bin.(binNb));
coeffs = polyfit(x,y,1);
f = polyval(coeffs,x);
plot(x, y,'o',x, f,'-','Color',cmap(4,:),'MarkerSize',2, 'MarkerFaceColor',cmap(4,:),'linewidth',2)
xlim([0 9])
text(max(x)/1.3,max(y)/20, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)))
hold on
binNb = sprintf('bin%d', 6); %plot binocular condition
y1=squeeze(normdp1p4.bin.(binNb));
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
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\adaptation_vs_diameter');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

 %% test significance of slope versus 0

    xtest = x;
    ytest = y1;
    linreg = fitlm(xtest,ytest);
    [linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linreg);
    Rsquared(1,1) = linreg.Rsquared.Ordinary;
