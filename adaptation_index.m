% This script intends to compute an adaptation index for every single unit
% both in the monocular condition and the binocular condition and to plot both
% monocular and binocular conditions
% This script follows the processing steps of BinocularAdaptationTrialSelection.m

% Developped by Loic Daumail, started 04-08-2021


%1) Get mean peak response values for all units in both mono vs bino
%conditions ===> cf binocular adaptation analysis
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'NoFiltMultiContSUA']; 
NoFiltMultiContSUA = load(channelfilename);
filenames = fieldnames(NoFiltMultiContSUA.NoFiltMultiContSUA);
meanPks = load([newdatadir 'all_unfiltered_data_peaks']); %peak values obtained with BinocularAdaptationTrialSelection.m


NdeAvgCont = [0,0.85];
class ={'M','P','K'};
bins =[1,6];

mean_pk1 = nan(length(fieldnames( meanPks.mean_peaks)),length(NdeAvgCont),length(class));
mean_pk4 = nan(length(fieldnames( meanPks.mean_peaks)),length(NdeAvgCont),length(class));

 for c = 1:length(class)   
    for i =1:length(fieldnames( meanPks.mean_peaks))
        filename = filenames{i};
        if strcmp(NoFiltMultiContSUA.NoFiltMultiContSUA.(filename).cellclass, class{c})
            for bin = 1:length(NdeAvgCont)
                binNb = sprintf('bin%d',bins(bin));
                if nnz(strcmp(fieldnames(meanPks.mean_peaks.(filename)),binNb))
                    
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
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\plots\hist_rast_sdf_mono_bino_cells');
%saveas(gcf,strcat(plotdir, '.png'));
%saveas(gcf,strcat(plotdir, '.svg'));



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
g(2,1).set_color_options('hue_range',[-100 100],'chroma',60);%,'legend','separate');
%g(2,1).set_color_options('hue_range',[-40 40],'chroma',30,'lightness',90);
g(2,1).set_title({'Adaptation index distribution of M cells, in the monocular and binocular conditions'});


f = figure('Position',[100 100 800 1400]);
%g.set_title({'Adaptation index distribution across all cells in the monocular and binocular conditions'});
g.draw();
set(f,'position',get(f,'position').*[1 1 1.15 1])
%plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\adaptation_index\plots\hist_sdf_mono_bino_allcells_combined_and_M_taller');
%saveas(gcf,strcat(plotdir, '.png'));
%saveas(gcf,strcat(plotdir, '.svg'));


%%
%% Example data from gramm tutorial
%data = load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\MATLAB\gramm-master\example_data.mat');
load 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\MATLAB\gramm-master\example_data.mat';
% Create a gramm object, provide x (year of production) and y (fuel economy) data,
% color grouping data (number of cylinders) and select a subset of the data
g=gramm('x',cars.Model_Year,'y',cars.MPG,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5);
%%% 

%% Tests of goodness of fit


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