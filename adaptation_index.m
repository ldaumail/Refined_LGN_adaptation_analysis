% This script intends to compute an adaptation index for every single unit
% both in the monocular condition and the binocular condition and to plot both
% monocular and binocular conditions
% This script follows the processing steps of "BinocularAdaptationTrialSelection.m"

% Developped by Loic Daumail, started 04-08-2021, last edited 06-02-2022


%1) Get mean peak response values for all units in both mono vs bino
%conditions ===> cf binocular adaptation analysis
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA_05022021']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);
meanPks = load([newdatadir 'all_unfiltered_data_peaks_05022021']); %peak values obtained with BinocularAdaptationTrialSelection.m


NdeAvgCont = [0,0.85];
class ={'M','P','K'};
bins =[1,6];

mean_pk1 = nan(length(fieldnames( meanPks.mean_peaks)),length(NdeAvgCont),length(class));
mean_pk4 = nan(length(fieldnames( meanPks.mean_peaks)),length(NdeAvgCont),length(class));

 for c = 1:length(class)   
    for i =1:length(fieldnames( meanPks.mean_peaks))
        filename = filenames{i};
        if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
            if nnz(strcmp(fieldnames(meanPks.mean_peaks.(filename)),'bin1')) && nnz(strcmp(fieldnames(meanPks.mean_peaks.(filename)),'bin6'))  
                for bin = 1:length(NdeAvgCont)
                    binNb = sprintf('bin%d',bins(bin));    
                    mean_pk1(i,bin,c) = meanPks.mean_peaks.(filename).(binNb)(1);
                    mean_pk4(i,bin,c) = meanPks.mean_peaks.(filename).(binNb)(4);
                end
            end
        end
    end
    
    %meanAllPk1 = nanmean(mean_pk1,1);
    %meanAllPk4 = nanmean(mean_pk4,1);
    %idxs(c) = ~isnan(meanAllPk1(:,c));
 end
 
 %2) Compute index
 adapt_idx = nan(length(fieldnames( meanPks.mean_peaks)),length(NdeAvgCont),length(class));
 for i = 1:length(mean_pk1(1,1,:))
    for j = 1:length(mean_pk1(1,:,1))
        for k = 1:length(mean_pk1(:,1,1))
            adapt_idx(k,j,i) = 2*(mean_pk1(k,j,i) - mean_pk4(k,j,i))/(mean_pk1(k,j,i) + mean_pk4(k,j,i));
        end
    end
 end
 
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue

 
%% overlay mono/bino histograms in same axes

%first need to convert data in the right structure for gramm toolbox
%standards
unit =  nan(2*length(adapt_idx(:,1,1)),1);
condition = cell(2*length(adapt_idx(:,1,1)),1);
cellclass = cell(2*length(adapt_idx(:,1,1)),1);
index =  nan(2*length(adapt_idx(:,1,1)),1);

for k = 1:length(adapt_idx(:,1,1))
    for j = 1:length(adapt_idx(1,:,1))
        for i = 1:length(adapt_idx(1,1,:))
            
            if ~isnan(adapt_idx(k,j,i))
                
                if j == 1
                    unit(k) = k;
                    condition{k} = 'Monocular';
                    index(k) = adapt_idx(k,j,i);
                    if i == 1
                        cellclass{k} = class{1};
                    else
                        if i == 2
                            cellclass{k} = class{2};
                        else
                            if i == 3
                                cellclass{k} = class{3};
                            end
                        end
                    end
                    %{
                    if isempty(condition{k})
                        condition{k} = char();
                        cellclass{k} = char();
                    end
                    %}
                else
                    if j == 2
                        unit(length(adapt_idx(:,1,1))+k) = k;
                        condition{length(adapt_idx(:,1,1))+k} = 'Binocular';
                        index(length(adapt_idx(:,1,1))+k) = adapt_idx(k,j,i);
                        if i == 1
                            cellclass{length(adapt_idx(:,1,1))+ k} = class{1};
                        else
                            if i == 2
                                cellclass{length(adapt_idx(:,1,1))+ k} = class{2};
                            else
                                if i == 3
                                    cellclass{length(adapt_idx(:,1,1))+ k} = class{3};
                                end
                            end
                        end
                        %{
                        if isempty(condition{length(adapt_idx(:,1,1))+ k})
                            condition{length(adapt_idx(:,1,1))+ k} = char();
                            cellclass{length(adapt_idx(:,1,1))+ k} = char();
                        end
                        %}
                    end
                    
                end
            end
            
        end
    end
end

for n = 1:length(condition)
    if isempty(condition{n})
        condition{n} = char();
        cellclass{n} = char();
        
    end
end

%% Create a table to store data
FileName = [filenames; filenames];

T = table(FileName, unit, cellclass, condition, index, 'VariableNames',{'File Name','Unit Number','Cell Class','Condition', 'Adaptation Index'}); 
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\data\AdaptationIndexData_05022021';
writetable(T,strcat(allfilename, '.csv'));


%% Raster + Histogram + density function
%(1) subplots for individual cell classes
clear g
g(1,1)=gramm('x',index,'y',unit,'color',condition);
g(1,2)=copy(g(1));
g(1,3)=copy(g(1));
%g(2,2)=copy(g(1));

%Raw data as raster plot
g(1,1).facet_grid(cellclass,[]);
g(1,1).geom_raster();
g(1,1).no_legend();
g(1,1).set_names('x','Adaptation Index','color','Legend','row','','y','');
%Histogram
g(1,2).facet_grid(cellclass,[]);
g(1,2).stat_bin('nbins',25);
g(1,2).set_title({'Adaptation index distributions of each cell class', 'in the monocular and binocular conditions'});
g(1,2).set_names('x','Adaptation Index','color','Legend','row','','y','Count');


%Kernel smoothing density estimate
g(1,3).facet_grid(cellclass,[]);
g(1,3).stat_density();
g(1,3).no_legend();
g(1,3).set_names('x','Adaptation Index','color','Legend','row','','y','Count');

f = figure('Position',[100 100 1400 550]);
g.draw();
set(f,'position',get(f,'position').*[1 1 1.15 1])
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\plots\hist_rast_sdf_mono_bino_cells_05022021');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));



%(2) Plot all cells together


clear g
g(1,1)=gramm('x',index,'y',unit,'color',condition);
g(1,2)=copy(g(1));
g(1,3)=copy(g(1));


%Raw data as raster plot
%g(1,1).facet_grid(cellclass,[]);
g(1,1).geom_raster();
g(1,1).no_legend();
g(1,1).set_names('x','Adaptation Index','color','Legend','row','','y','');
%Histogram
%g(1,2).facet_grid(cellclass,[]);
g(1,2).stat_bin('nbins',8);
g(1,2).set_title({'Adaptation index distribution across all cells', 'in the monocular and binocular conditions'});
g(1,2).set_names('x','Adaptation Index','color','Legend','row','','y','Count');


%Kernel smoothing density estimate
%g(1,3).facet_grid(cellclass,[]);
g(1,3).stat_density();
g(1,3).no_legend();
g(1,3).set_names('x','Adaptation Index','color','Legend','row','','y','Count');

f = figure('Position',[100 100 1400 550]);
g.draw();
set(f,'position',get(f,'position').*[1 1 1.15 1])
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\plots\hist_rast_sdf_mono_bino_allcells_combined');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));


%% same plot as the previous one with histogram and density function in the same panel
clear g
g(2,1)=gramm('x',index,'y',unit,'color',condition, 'subset', strcmp(cellclass,'M'));
g(1,1)=gramm('x',index,'y',unit,'color',condition);

g(1,1).stat_bin('nbins',25);
g(1,1).stat_density();
g(1,1).set_names('x','Adaptation Index','color','Legend','row','','y','Count');
g(1,1).set_title({'Adaptation index distribution across all cells in the monocular and binocular conditions'});


g(2,1).stat_bin('nbins',15);
g(2,1).stat_density();
g(2,1).set_names('x','Adaptation Index','color','Legend','row','','y','Count');
%g(2,1).set_color_options('hue_range',[-100 100],'chroma',60);%,'legend','separate');
%g(2,1).set_color_options('hue_range',[-40 40],'chroma',30,'lightness',90);
g(2,1).set_title({'Adaptation index distribution of M cells, in the monocular and binocular conditions'});


f = figure('Position',[100 100 800 1400]);
%g.set_title({'Adaptation index distribution across all cells in the monocular and binocular conditions'});
g.draw();
set(f,'position',get(f,'position').*[1 1 1.15 1])
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\plots\hist_sdf_mono_bino_allcells_combined_and_M_taller_05022021');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));


%% SAME  plot, only including the whole population
%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

clear g
g(1,1)=gramm('x',index,'y',unit,'color',condition);

g(1,1).stat_bin('nbins',25,'geom','overlaid_bar');
g(1,1).stat_density();
%g(1,1).set_color_options('map', [251/255 154/255 153/255;160/255 160/255 160/255]);
g(1,1).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]); 
g(1,1).set_names('x','Adaptation Index','color','Legend','row','','y','Count');
g(1,1).set_title({'Adaptation index distribution across all cells in the monocular and binocular conditions'});


f = figure('Position',[100 100 800 1000]);
%g.set_title({'Adaptation index distribution across all cells in the monocular and binocular conditions'});
g.draw();
set(f,'position',get(f,'position').*[1 1 1.15 1])
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\plots\hist_sdf_mono_bino_allcells_adaptindex_purpOr');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));
%% Same plot as above, with scatter plots and vertical

clear g
f = figure('Position',[100 100 800 1000]);
set(f,'position',get(f,'position').*[1 1 1.15 1])
  %jitter
 %g(1,1)=gramm('x',index,'y',unit,'color',condition);
g(1,1) = gramm('x',condition,'y',index, 'color',condition); 
g(1,1).geom_jitter('width',0.4,'height',0); %Scatter plot
g(1,1).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]); 
g(1,1).axe_property( 'xlim',[0 4] , 'ylim',[-0.3 0.5]); 
% add confidence interval 95%
ci_low = [nanmean(index(strcmp(condition, 'Monocular'))) - std(index(strcmp(condition, 'Monocular')),0,'omitnan')/sqrt(length(index(strcmp(condition, 'Monocular')))); nanmean(index(strcmp(condition, 'Binocular'))) - std(index(strcmp(condition, 'Binocular')),0,'omitnan')/sqrt(length(index(strcmp(condition, 'Binocular'))))];
ci_high = [nanmean(index(strcmp(condition, 'Monocular'))) + std(index(strcmp(condition, 'Monocular')),0,'omitnan')/sqrt(length(index(strcmp(condition, 'Monocular')))); nanmean(index(strcmp(condition, 'Binocular'))) + std(index(strcmp(condition, 'Binocular')),0,'omitnan')/sqrt(length(index(strcmp(condition, 'Binocular'))))];
g(1,1).update('x',[1;2], 'y', [nanmean(index(strcmp(condition, 'Monocular'))); nanmean(index(strcmp(condition, 'Binocular')))],...
    'ymin',ci_low,'ymax',ci_high,'color',[1;2]);
g(1,1).geom_point('dodge',0.5);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]); 
g(1,1).axe_property('xlim',[0 4]); 
g(1,1).set_point_options('base_size',7);
%}
%bar
g(1,2)=gramm('x',index,'y',unit, 'color',condition);
g(1,2).stat_bin('nbins',25,'geom','overlaid_bar');
%g(1,2).stat_bin('normalization','probability','nbins',20,'geom','overlaid_bar');
g(1,2).stat_density();
g(1,2).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]); 
g(1,2).axe_property('xlim',[-0.3 0.5], 'ylim', [0 5]); 
g(1,2).set_names('x','Adaptation index','color','Legend','row','','y','Count');
g(1,2).coord_flip();
g.draw();
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\plots\hist_scat_mono_bino_allcells_adaptindex_purpOr');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% Plot the difference of adaptation indices across adaptation conditions
%Use trial selected data from other trial selection ('same as adaptation_and_binocular_interaction_plots.m', that has more single unit data


newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']; 
peak_aligned_trials = load(channelfilename);

%Normalize data differently so Pk1 mono also has a normal distribution
filenames = fieldnames(peak_aligned_trials.peak_aligned_trials);
bins = [1,6];

mean_peaks = nan(250, 4, 2,length(filenames) );
for i = 1: length(filenames)
    filename = filenames{i};
    if length(fieldnames(peak_aligned_trials.peak_aligned_trials.(filename).origin)) == 2
        for b = 1:2
            binN = sprintf('bin%d',bins(b));
             %compute mean peak responses
             for p = 1:4
                 pkN = sprintf('pk%d',p);
                 mean_peaks(:,p,b,i) = mean(peak_aligned_trials.peak_aligned_trials.(filename).origin.(binN).(pkN),2);
             end
        end
    end
end
%store norm peak values
 peaks = nan(4,2,length(filenames));
for i = 1:length(filenames)
    for b = 1:2
        for p = 1:4
            peaks(p,b,i) = max(mean_peaks(:,p,b,i));
        end
    end
end

%compute adaptation index
adapt_idx = nan(2,length(filenames));
for i =1:length(filenames)
    for b =1:2
        adapt_idx(b,i) = 2*(peaks(1,b,i)-peaks(4,b,i))/(peaks(1,b,i)+peaks(4,b,i));
    end
end

index = adapt_idx(:,~isnan(adapt_idx(1,:)));

idx_diff = index(1,:) - index(2,:);

unitnb = 1:length(idx_diff);
%get pvalues from lmer results with Dunnett correction
pvalues = dlmread('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\lmer_results_orig_03032020_corrected_dunnett.csv', ',', 1,1);

%clear out nans
pk1pk4pval = pvalues(~isnan(pvalues(:,3)),3);

pk1pk4pval = pk1pk4pval(~isnan(adapt_idx(1,:)));

suppsig = pk1pk4pval(pk1pk4pval < 0.05 & index(1,:)' > 0);
facsig = pk1pk4pval(pk1pk4pval < 0.05 & index(1,:)' < 0);


nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('Blues', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);

colors = [cmaps(2).map(7,:); cmaps(1).map(7,:)];
adapType = [pk1pk4pval < 0.05 & index(1,:)' > 0,pk1pk4pval < 0.05 & index(1,:)' < 0 ];
%colormap(cmap);
clear g
f = figure('Position',[100 100 800 1000]);
set(f,'position',get(f,'position').*[1 1 1.15 1])
for t = 1:size(adapType,2)
 condition = [repmat({'all'},length(idx_diff),1); repmat({'sig'},length(idx_diff(unitnb(adapType(:,t)))),1)];
  %jitter
g(1,2*(t-1)+1) = gramm('x',condition,'y', [idx_diff,idx_diff(adapType(:,t))]', 'color',condition); %[unitnb,1:length(unitnb(adapType(:,t)))]
g(1,2*(t-1)+1).geom_jitter('width',0.4,'height',0); %Scatter plot
g(1,2*(t-1)+1).set_color_options('map',[cmaps(3).map(3,:);colors(t,:)]); 
g(1,2*(t-1)+1).axe_property( 'xlim',[0 4] , 'ylim',[-1 1]); 
% add confidence interval 95%
ci_low = [mean(idx_diff) - std(idx_diff)/sqrt(length(idx_diff)); mean(idx_diff(adapType(:,t))) - std(idx_diff(adapType(:,t)))/sqrt(length(idx_diff(adapType(:,t)))) ];
ci_high = [mean(idx_diff) + std(idx_diff)/sqrt(length(idx_diff)); mean(idx_diff(adapType(:,t))) + std(idx_diff(adapType(:,t)))/sqrt(length(idx_diff(adapType(:,t)))) ];
g(1,2*(t-1)+1).update('x',[1;2], 'y', [mean(idx_diff); mean(idx_diff(adapType(:,t)))],...
    'ymin',ci_low,'ymax',ci_high,'color',[1;2]);
g(1,2*(t-1)+1).geom_point('dodge',0.5);
g(1,2*(t-1)+1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,2*(t-1)+1).set_color_options('map',[cmaps(3).map(3,:);colors(t,:)]); 
g(1,2*(t-1)+1).axe_property('xlim',[0 4]); 
g(1,2*(t-1)+1).set_point_options('base_size',7);

%bar
g(1,2*(t-1)+2)=gramm('x',idx_diff,'y',unitnb);
g(1,2*(t-1)+2).stat_bin('nbins',25,'geom','overlaid_bar');
g(1,2*(t-1)+2).stat_density();
%g(1,1).set_color_options('map', [251/255 154/255 153/255;160/255 160/255 160/255]);
g(1,2*(t-1)+2).set_color_options('map',cmaps(3).map(3,:)); 
g(1,2*(t-1)+2).axe_property('xlim',[-1 1], 'ylim', [0 6]); 
g(1,2*(t-1)+2).set_names('x','Adaptation Index Difference','color','Legend','row','','y','Count');

g(1,2*(t-1)+2).update('x',idx_diff(adapType(:,t)), 'y',unitnb(adapType(:,t)))
g(1,2*(t-1)+2).stat_bin('nbins',25,'geom','overlaid_bar');
g(1,2*(t-1)+2).stat_density();
g(1,2*(t-1)+2).set_color_options('map',colors(t,:)); 
g(1,2*(t-1)+2).axe_property('xlim',[-1 1], 'ylim', [0 6]); 

end
g(1,2).coord_flip();
g(1,4).coord_flip();
g.set_title({'Adaptation index difference distribution across all cells'});
g.draw();
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\plots\hist_sdf_mono_bino_allcells_adaptindex_difference');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));


%% Test for significant difference between stimulation conditions

%1) Population level
[hpop, ppop, cipop, statspop] = ttest(idx_diff);
%2) Suppressed units
[hsup, psup, cisup, statsup] = ttest(idx_diff(adapType(:,1)));
%% Plots of mon/bin distributions of indices
nbins = [15, 20; 15, 15; 15, 3];
xlims = [-0.2 0.4;-0.25 0.3;0 0.4];
for c =1:3
 f = figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);;
 subplot(1,2,1)
 histogram(adapt_idx(:,1,c),nbins(c,1), 'EdgeColor', 'none', 'FaceColor',col(c,:));
 xlim(xlims(c,:))
 ylim([0 6])
 set(gca, 'Box', 'off')
 set(gca, 'linewidth',2)
 set(gca, 'fontsize', 24)
 title({'Adaptation index distribution of',sprintf('the %s cell class in the monocular condition',char(class(c)))},'FontSize', 24)

 hold on
 subplot(1,2,2)
  h = histogram(adapt_idx(:,2,c),nbins(c,2), 'EdgeColor', 'none', 'FaceColor', [180/255 180/255 180/255]);
 %alpha(h,1)

 xlim(xlims(c,:))
 ylim([0 6])
 set(gca, 'Box', 'off')
 set(gca, 'linewidth',2)
 set(gca, 'fontsize', 24)
title({'Adaptation index distribution of',sprintf('the %s cell class in the binocular condition',char(class(c)))},'FontSize', 24)
set(f,'position',get(f,'position').*[1 1 1.15 1])
%plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\plots\',strcat(sprintf('mono_bino_%s_cell',class{c})));
%saveas(gcf,strcat(plotdir, '.png'));
%saveas(gcf,strcat(plotdir, '.svg'));

                
end


%%
%% Example data from gramm tutorial
%data = load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\MATLAB\gramm-master\example_data.mat');
load 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\MATLAB\gramm-master\example_data.mat';
% Create a gramm object, provide x (year of production) and y (fuel economy) data,
% color grouping data (number of cylinders) and select a subset of the data
g=gramm('x',cars.Model_Year,'y',cars.MPG,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5);
%%% 

%% Tests difference between two distributions (monocular/binocular)

%%All units of all cell classes combined
cond1 = index(strcmp(condition, 'Monocular'));
cond2 = index(strcmp(condition, 'Binocular'));
%two sample Kolmogorov Smirnov test
[h,p,ks2stat] = kstest2(cond1,cond2);

%Qi square goodness of fit test
[h1,p1,stats1] = chi2gof(cond1);
[h2,p2,stats2] = chi2gof(cond2);

%testing normality: 
%Shapiro-Wilk
[H1, pValue1, W1] = swtest(cond1, 0.05);
[H2, pValue2, W2]= swtest(cond2, 0.05);


%Look at M cells only 
cond1 = index(strcmp(condition, 'Monocular') & strcmp(cellclass, 'M'));
cond2 = index(strcmp(condition, 'Binocular') & strcmp(cellclass, 'M'));

%Qi square goodness of fit test
[h1,p1,stats1] = chi2gof(cond1);
[h2,p2,stats2] = chi2gof(cond2);
%%Shapiro-Wilk
[H1, pValue1, W1] = swtest(cond1, 0.05);
[H2, pValue2, W2]= swtest(cond2, 0.05);


%% Plot scatter plot of adapt index in mono vs bino condition to create clusters

%all cells
cond1 = index(strcmp(condition, 'Monocular'));
cond2 = index(strcmp(condition, 'Binocular'));


figure();
for i = 1:length(cond1)
    if cond2(i) < cond1(i)
        scatter(cond1(i), cond2(i), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
        hold on
    else
        scatter(cond1(i), cond2(i), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
        hold on
    end
    
end
hold on
xf = [min(cond1), max(cond1)];
plot(xf, polyval(polyfit(cond1,cond2,1), xf));

xlim([-0.3 0.4])
ylim([-0.3 0.4])
xlabel('Monocular adaptation index')
ylabel('Binocular adaptation index')
title('All cells adaptation indices')

%M cells
cond1 = index(strcmp(condition, 'Monocular') & strcmp(cellclass, 'M'));
cond2 = index(strcmp(condition, 'Binocular') & strcmp(cellclass, 'M'));

figure();
for i = 1:length(cond1)
    if cond2(i) < cond1(i)
        scatter(cond1(i), cond2(i), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
        hold on
    else
        scatter(cond1(i), cond2(i), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
        hold on
    end
    
end
hold on
xf = [min(cond1), max(cond1)];
plot(xf, polyval(polyfit(cond1,cond2,1), xf));

xlabel('Monocular adaptation index')
ylabel('Binocular adaptation index')
title('M cells adaptation indices')

%P cells
cond1 = index(strcmp(condition, 'Monocular') & strcmp(cellclass, 'P'));
cond2 = index(strcmp(condition, 'Binocular') & strcmp(cellclass, 'P'));

figure();
for i = 1:length(cond1)
    if cond2(i) < cond1(i)
        scatter(cond1(i), cond2(i), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
        hold on
    else
        scatter(cond1(i), cond2(i), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
    end
    
end
hold on
xf = [min(cond1), max(cond1)];
plot(xf, polyval(polyfit(cond1,cond2,1), xf));

xlabel('Monocular adaptation index')
ylabel('Binocular adaptation index')
title('P cells adaptation indices')


%K cells
cond1 = index(strcmp(condition, 'Monocular') & strcmp(cellclass, 'K'));
cond2 = index(strcmp(condition, 'Binocular') & strcmp(cellclass, 'K'));

figure();
for i = 1:length(cond1)
    if cond2(i) < cond1(i)
        scatter(cond1(i), cond2(i), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
        hold on
    else
        scatter(cond1(i), cond2(i), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
    end
    
end
xlabel('Monocular adaptation index')
ylabel('Binocular adaptation index')
title('K cells adaptation indices')

%% Compute adaptation index difference in the binocular versus monocular condition







%{
%% ROC analyis of distributions in the monocular versus binocular condition

reps   = 10000;
all_sigs95 = nan(length(Ses),1);
all_sigs90 = nan(length(Ses),1);

part1 = index(strcmp(condition, 'Monocular'));
part2 = index(strcmp(condition, 'Binocular'));

if nanmean(part1) > nanmean(part2)
    cond1               = part1;
    cond2               = part2;
else
    cond1               = part2;
    cond2               = part1;
end
%% Plot cond1 and cond2 distributions  in box/scatter plots
figure(); boxplot([cond1 cond2],'notch','off','labels',{'Binocular','Monocular'}); hold on    %here we now know that cond1 = part2 = binocular/ cond2 = part1 = Monocular
x=ones(length(cond1),1).*(1+(rand(length(cond1),1)-0.5)/5);
x1=ones(length(cond1),1).*(1+(rand(length(cond1),1)-0.5)/10);
f1=scatter(x,cond1,'k','filled');f1.MarkerFaceAlpha = 0.4;hold on
f2=scatter(x1*2,cond2,'k','filled');f1.MarkerFaceAlpha = 0.4;
%ylim([0 400])
ylabel('Adaptation Index')

%% Plot AUC of our two samples

[X,Y,T,AUC]           = perfcurve([ones(length(cond1),1); repmat(2,length(cond2),1)],[cond1' cond2'],1);
NP                    = length(cond2);
PR                    = length(cond1);
catdat                = [cond1' cond2'];
figure();
plot(X,Y)
set(gca, 'box','off')

%% Plot Receiver Operating Characteristics curve examples

 shufPR         = catdat(randperm(length(catdat),PR));
 shufNP         = catdat(randperm(length(catdat),NP));
[Xshuf,Yshuf,~, shufAUC]    = perfcurve([ones(PR,1); repmat(2,NP,1)],[shufPR shufNP],1);

figure();
plot(Xshuf,Yshuf)
set(gca, 'box','off')

 %% Plot randomized distribution
 
 NP                    = length(cond2);
 PR                    = length(cond1);
 catdat                = [cond1' cond2'];
 for r       = 1:reps
     clear shufNP shufPR
     shufPR         = catdat(randperm(length(catdat),PR));
     shufNP         = catdat(randperm(length(catdat),NP));
     [~,~,~,...
         shufAUC(r)]    = perfcurve([ones(PR,1); repmat(2,NP,1)],[shufPR shufNP],1);
 end
 critT95         = quantile(shufAUC,.95);
 
 figure();
 histogram(shufAUC)
 hold on
 plot([critT95, critT95], [0, 900], 'linewidth', 2)
 hold on
 plot([AUC, AUC], [0, 900], 'Color', 'k','linewidth', 2)
 set(gca, 'box', 'off')

 critT90         = quantile(shufAUC,.90);
%}