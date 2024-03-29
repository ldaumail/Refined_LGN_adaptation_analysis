%This script is developped following trial selection with
%"BinocularAdaptationTrialSelection.m" and adaptation and binocular interaction analysis 
% with "lmer_peaks_binocular_adaptation.R"
%Last edited on 05-25-2022 by Loic Daumail

%% Part 1 : Single Units Analysis (Part 2 for MUA in other script in MUA section)

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']; 
peak_aligned_trials = load(channelfilename);

%1) Plot all normalized data in mono vs bino condition

%We first need to normalize all mean responses. Since trials are stored 4 times as each time is triggered to one peak,  
filenames = fieldnames(peak_aligned_trials.peak_aligned_trials);
bins = [1,6];
norm_aligned_resps = nan(250, 4, 2,length(filenames) );
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
     %normalize in relation to all mean peaks in the monocular condition
         %norm_aligned_resps(:,:,b,i) = (mean_peaks(:,:,b,i) - min(mean_peaks(:,1,1,i),[], 'all'))./(max(mean_peaks(:,1,1,i), [], 'all') - min(mean_peaks(:,1,1,i),[], 'all'));
         norm_aligned_resps(:,:,b,i) = (mean_peaks(:,:,b,i))./(max(mean_peaks(:,1,1,i), [], 'all'));
        end
    end
end

%plot overall mean in binocular and monocular condition. Add confidence
%interval
figure();
     mean_unit = squeeze(nanmean(norm_aligned_resps,4));
     stdev = squeeze(std(norm_aligned_resps,[],4, 'omitnan'));
     for pn =1:4
         h = subplot(1,4,pn);
         %first plot monocular condition
         plot(-125:124, mean_unit(:,pn,1));
         hold on
         h1= ciplot( mean_unit(:,pn,1)+ 1.96*stdev(:,pn,1)/sqrt(length(norm_aligned_resps(1,1,1,:))), mean_unit(:,pn,1)-1.96*stdev(:,pn,1)/sqrt(length(norm_aligned_resps(1,1,1,:))),[-125:124],[40/255 40/255 40/255],0.1);
         set(h1, 'edgecolor','none')
         hold on 
         plot(-125:124, mean_unit(:,pn,2), 'r');
         hold on
         h2= ciplot( mean_unit(:,pn,2)+ 1.96*stdev(:,pn,2)/sqrt(length(norm_aligned_resps(1,1,1,:))), mean_unit(:,pn,2)-1.96*stdev(:,pn,2)/sqrt(length(norm_aligned_resps(1,1,1,:))),[-125:124],'r',0.1);
         set(h2, 'edgecolor','none')
       
         set(h,'position',get(h,'position').*[1 1 1.15 1])
         ylim([0 1])
         xlim([-125 125])
         set(gca,'box','off')
         set(gca, 'linewidth',2)
         ylabel({'\fontsize{14}Spike Rate (normalized)'});
         if pn > 1
             ax1 = gca;
             ax1.YAxis.Visible = 'off';
         end
     end
     
     sgtitle({'Mean responses in the monocular condition','vs binocular condition'}, 'Interpreter', 'none')
     xlabel('Time from cycle peak (ms)')
     set(gcf,'Units','inches')
     
%brainstorm with other ways like Gramm
%% plot horizontal histograms for each peak

%store norm peak values
 peaks = nan(4,2,length(filenames));
for i = 1:length(filenames)
    for b = 1:2
        for p = 1:4
            peaks(p,b,i) = max(norm_aligned_resps(:,p,b,i));
        end
    end
end
%%find a way to melt matrix in
%%p1p1p1p1....p2p2p2p2....p3p3p3p3....p4p4p4p4p4 fashion
peakvals = [squeeze(peaks(1,:,:))';squeeze(peaks(2,:,:))';squeeze(peaks(3,:,:))';squeeze(peaks(4,:,:))'];
%peakvals = reshape(peaks, [length(peaks(:,1,1))*length(peaks(1,1,:)), length(peaks(1,:,1))]);
linPeakVals = reshape(peakvals, [2*length(peakvals(:,1)),1]); 
condition = [repmat({'Monocular'},length(peakvals(:,1)),1); repmat({'Binocular'},length(peakvals(:,1)),1)];
unit = repmat(1:length(peaks(1,1,:)),1,8)';
peakLabel = repmat([repmat({'Pk1'}, length(peaks(1,1,:)),1);repmat({'Pk2'}, length(peaks(1,1,:)),1);repmat({'Pk3'}, length(peaks(1,1,:)),1);repmat({'Pk4'}, length(peaks(1,1,:)),1)],2,1);
%col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
%col(2,:) = [251/255 154/255 153/255]; % -- red
%col(3,:) = [146/255 197/255 222/255]; % -- blue


clear g
figure('Position',[100 100 1400 600]);
g(1,1)=gramm('x',linPeakVals, 'color',condition);
g(1,1).set_names('x','Spiking activity (Normalized)', 'column', '');
g(1,1).facet_grid([],peakLabel); %Provide facets
g(1,1).geom_jitter('width',0,'height',0.2);
%g(1,1).stat_bin('fill','transparent'); 
g(1,1).stat_bin('geom','overlaid_bar');%Histogram

g(1,1).set_color_options('map', [251/255 154/255 153/255;160/255 160/255 160/255]);
g(1,1).set_title('Population peak responses in the binocular and monocular conditions');
g(1,1).axe_property('xlim',[0.4 1.8],'ylim',[-2 15]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten

%Set global axe properties
g.axe_property('TickDir','out');
g.coord_flip();
%g.set_title({'Adaptation index distribution across all cells in the monocular and binocular conditions'});
g.draw();
%plot(g.results.stat_bin().child_axe_handle,[-2 -2],[0 50],'k:','LineWidth',2)

%set(f,'position',get(f,'position').*[1 1 1.15 1])
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\hist_mono_bino_allcells_peaks');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% Plot jitter scatter plot + horizontal histogram on the side for each peak
%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);


clear g
figure('Position',[100 100 1400 600]);
for p =1:4
%Create a scatter plot
meanpMono = nanmean(linPeakVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p))));
meanpBino = nanmean(linPeakVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p))));
g(1,2*(p-1)+1)=gramm('x',linPeakVals,'color',condition,'subset', strcmp(peakLabel,sprintf('Pk%d',p)));
g(1,2*(p-1)+1).set_names('x','Spiking activity (Normalized)','column','');
g(1,2*(p-1)+1).geom_jitter('width',0,'height',0.2); %Scatter plot
g(1,2*(p-1)+1).geom_vline('xintercept', meanpMono, 'style', '-or');
g(1,2*(p-1)+1).geom_vline('xintercept', meanpBino, 'style', '-b');
%g(1,2*(p-1)+1).axe_property('Ygrid','on', 'ylim',[0.3 1.7],'YTickLabel','','YTick',''); 
g(1,2*(p-1)+1).axe_property('xlim', [0.6 1.3], 'ylim',[0.3 1.7]); 
g(1,2*(p-1)+1).no_legend();

%Create y data histogram on the right
g(1,2*(p-1)+2)=gramm('x',linPeakVals,'color',condition, 'subset', strcmp(peakLabel,sprintf('Pk%d',p)));
%g(2,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
%    'legend',false,...
%    'margin_height',[0.1 0.02],...
%    'margin_width',[0.02 0.05],...
%    'redraw',false);
g(1,2*(p-1)+2).geom_vline('xintercept', meanpMono, 'style', '-or');
g(1,2*(p-1)+2).geom_vline('xintercept', meanpBino, 'style', '-b');
g(1,2*(p-1)+2).set_names('x','');
g(1,2*(p-1)+2).stat_bin('geom','overlaid_bar'); %histogram
g(1,2*(p-1)+2).axe_property('xlim',[0.3 1.7],'ylim',[-2 15],'XTickLabel','','XTick',''); 
g(1,2*(p-1)+2).no_legend();

end
%Set global axe properties
g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.axe_property('xlim',[0.4 1.8]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g.coord_flip();
g.set_title('Population peak responses in the binocular and monocular conditions');
g.set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]);
%g.set_color_options('map', [251/255 154/255 153/255;160/255 160/255 160/255]);
g.draw();
%set(findobj(gcf, 'type','axes'), 'Visible','off')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_hist_mono_bino_allcells_peaks_purpOr');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));


%% %% Plot jitter scatter plot + point bar +-1.96sem plot on the side for each peak

%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

%plot
clear g
figure('Position',[100 100 1400 600]);
for p =1:4
%Create a scatter plot
meanpMono = nanmean(linPeakVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p))));
meanpBino = nanmean(linPeakVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p))));
g(1,2*(p-1)+1)=gramm('x',linPeakVals,'color',condition,'subset', strcmp(peakLabel,sprintf('Pk%d',p)));
g(1,2*(p-1)+1).set_names('x','Spiking activity (Normalized)','column','');
g(1,2*(p-1)+1).geom_jitter('width',0,'height',0.2); %Scatter plot
%g(1,2*(p-1)+1).geom_vline('xintercept', meanpMono, 'style', '-k');
%g(1,2*(p-1)+1).geom_vline('xintercept', meanpBino, 'style', '-p');
%g(1,2*(p-1)+1).axe_property('Ygrid','on', 'ylim',[0.3 1.7],'YTickLabel','','YTick',''); 
%g(1,2*(p-1)+1).axe_property('xlim', [0.6 1.3], 'ylim',[0.3 1.7]); 
g(1,2*(p-1)+1).no_legend();
g(1,2*(p-1)+1).axe_property('xlim',[0.4 1.8],'ylim',[0.6 1.3],'YTickLabel','','YTick',''); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g(1,2*(p-1)+1).coord_flip();
%Create point box plot on the right

shortCondition = {'Monocular', 'Binocular'};
meanPks = [meanpMono, meanpBino];
semMonoPks = std(linPeakVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p))), 'omitnan')/sqrt(length(linPeakVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p)))));
semBinoPks = std(linPeakVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p))), 'omitnan')/sqrt(length(linPeakVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p)))));
semPks = [semMonoPks, semBinoPks];

ci_low = meanPks - 1.96*semPks;
ci_high = meanPks + 1.96*semPks;

g(1,2*(p-1)+2)=gramm('x', shortCondition, 'y',meanPks,...
    'ymin',ci_low,'ymax',ci_high,'color',shortCondition);
g(1,2*(p-1)+2).set_names('color','Condition','y','Mean Peak spike rate (Normalized)');
%g(1,2*(p-1)+2).set_color_options('map',cmap);
g(1,2*(p-1)+2).geom_point('dodge',0.2);
g(1,2*(p-1)+2).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,2*(p-1)+2).axe_property('ylim',[0.4 1.8],'XTickLabel','','XTick',''); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g(1,2*(p-1)+2).no_legend();

end
%Set global axe properties
g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);

g.set_title('Population peak responses in the binocular and monocular conditions');
g.set_color_options('map',[cmap(3,:);160/255 160/255 160/255]);
g.draw();
%set(findobj(gcf, 'type','axes'), 'Visible','off')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_boxpointbar_mono_bino_allcells_peaks_std');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% %% Plot jitter scatter plot + point bar +-1.96sem plot on above for each peak

%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

%plot
clear g
figure('Position',[100 100 800 800]);
for p =1:4
%Create a scatter plot
meanpMono = nanmean(linPeakVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p))));
meanpBino = nanmean(linPeakVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p))));
g(1,p)=gramm('x',linPeakVals,'color',condition,'subset', strcmp(peakLabel,sprintf('Pk%d',p)));
g(1,p).set_names('x','Spiking activity (Normalized)','column','');
g(1,p).geom_jitter('width',0,'height',0.1); %Scatter plot
%g(1,2*(p-1)+1).geom_vline('xintercept', meanpMono, 'style', '-k');
%g(1,2*(p-1)+1).geom_vline('xintercept', meanpBino, 'style', '-p');
%g(1,2*(p-1)+1).axe_property('Ygrid','on', 'ylim',[0.3 1.7],'YTickLabel','','YTick',''); 
%g(1,2*(p-1)+1).axe_property('xlim', [0.6 1.3], 'ylim',[0.3 1.7]); 
g(1,p).no_legend();
g(1,p).axe_property('xlim',[0.4 1.8],'ylim',[0.6 1.3],'YTickLabel','','YTick',''); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g(1,p).coord_flip();

%Create point box plot below

shortCondition = {'Monocular', 'Binocular'};
meanPks = [meanpMono, meanpBino];
semMonoPks = std(linPeakVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p))), 'omitnan')/sqrt(length(linPeakVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p)))));
semBinoPks = std(linPeakVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p))), 'omitnan')/sqrt(length(linPeakVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p)))));
semPks = [semMonoPks, semBinoPks];

ci_low = meanPks - 1.96*semPks;
ci_high = meanPks + 1.96*semPks;

g(2,p)=gramm('x', shortCondition, 'y',meanPks,...
    'ymin',ci_low,'ymax',ci_high,'color',shortCondition);
g(2,p).set_names('color','Condition','y','Mean Peak spike rate (Normalized)');
%g(1,2*(p-1)+2).set_color_options('map',cmap);
g(2,p).geom_point('dodge',0.2);
g(2,p).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(2,p).axe_property('ylim',[0.85 1.1],'XTickLabel','','XTick',''); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g(2,p).no_legend();

end
%Set global axe properties
g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);

g.set_title('Population peak responses in the binocular and monocular conditions');
g.set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]);
g.draw();
%set(findobj(gcf, 'type','axes'), 'Visible','off')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_boxpointbar_mono_bino_allcells_peaks_95ci');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%%
%% %% Plot jitter scatter plot + point bar +-1.96sem plot overlay for each peak combining both conditions

%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('Blues', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

%plot
clear g
figure('Position',[100 100 800 800]);

%Create a scatter plot
%combine mono and bino data
pksMono = linPeakVals( strcmp(condition, 'Monocular'));
pksBino = linPeakVals( strcmp(condition, 'Binocular'));
meanpks = nanmean([pksMono,pksBino],2);
g(1,1)=gramm('x',categorical(peakLabel(1:length(peakLabel)/2)),'y',meanpks);
g(1,1).geom_jitter('width',0.3,'height',0.2);
g(1,1).set_names('color','Peak','y','Spiking activity (Normalized)','x','Peak #');
g(1,1).no_legend();
g(1,1).axe_property('ylim',[0.4 1.8]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g(1,1).set_color_options('map',cmap(4,:));

%Create point box plot overlay
allpksMono = linPeakVals(strcmp(condition, 'Monocular'));
allpksBino = linPeakVals(strcmp(condition, 'Binocular'));
meanUnitPks = nanmean([allpksMono,allpksBino],2);
meanPks = nan(4,1);
semPks = nan(4,1);
for p =1:4
meanPks(p) = nanmean(meanUnitPks( strcmp(peakLabel(1:length(peakLabel)/2),sprintf('Pk%d',p))));
semPks(p) = std(meanUnitPks(strcmp(peakLabel(1:length(peakLabel)/2),sprintf('Pk%d',p))), 'omitnan')/sqrt(length(meanUnitPks(strcmp(peakLabel(1:length(peakLabel)/2),sprintf('Pk%d',p)))));
end
ci_low = meanPks - 1.96*semPks;
ci_high = meanPks + 1.96*semPks;
shortPks = {'Pk1','Pk2','Pk3','Pk4'};
g(1,1).update('x',categorical(shortPks), 'y',meanPks,...
    'ymin',ci_low,'ymax',ci_high);
g(1,1).set_color_options('map',cmaps(1).map(4,:));
g(1,1).geom_point('dodge',0.2);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
%Set global axe properties
%g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.set_title('Population peak responses');
g.draw();
%set(findobj(gcf, 'type','axes'), 'Visible','off')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_boxpointbar_conditions_fused_allcells_peaks_95ci');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% %% Plot violin plot + point bar +-1.96sem plot on above for each peak

%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

%plot
clear g
figure('Position',[100 100 800 800]);
for p =1:4

%Create point box plot +violin plot
pkLinVals = linPeakVals( strcmp(peakLabel,sprintf('Pk%d',p)));
g(1,p)=gramm('x', condition( strcmp(peakLabel,sprintf('Pk%d',p))), 'y',pkLinVals,'color',condition( strcmp(peakLabel,sprintf('Pk%d',p))));
g(1,p).stat_violin('fill','transparent' );
g(1,p).stat_boxplot('width',0.15 );
g(1,p).axe_property('ylim',[0.4 1.75],'XTickLabel','','XTick',''); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g(1,p).no_legend();

end
%Set global axe properties
g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.set_title('Population peak responses in the binocular and monocular conditions');
g.set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]);
g.draw();
%set(findobj(gcf, 'type','axes'), 'Visible','off')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\violin_mono_bino_allcells_peaks');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% Stats on variance
peaks = nan(2,4,length(filenames));
for i = 1:length(filenames)
    for b = 1:2
        for p = 1:4
            peaks(b,p,i) = max(mean_peaks(:,p,b,i));
        end
    end
end

norm_peaks = nan(2,4,length(filenames));
for i = 1:length(filenames)
    for b = 1:2
        for p = 1:4
            %norm_peaks(b,p,i) = (peaks(b,p,i) - min(peaks(b,:,i)))/(max(peaks(b,:,i))-min(peaks(b,:,i)));
            norm_peaks(b,p,i) = peaks(b,p,i)/peaks(1,1,i);

        end
    end
end
%first try to test variance on peaks in the monocular condition
%set each peak as a columns
monoPeaks = squeeze(norm_peaks(1,:,:))';
binoPeaks = squeeze(norm_peaks(2,:,:))';
[pval,stats] = vartestn(binoPeaks);

%% Plot jitter scatter plot + horizontal histogram on the side for each mono peak - bino peak difference
diffLinPeakVals = linPeakVals(strcmp(condition, 'Monocular')) - linPeakVals(strcmp(condition, 'Binocular'));
diffPeakLabel = peakLabel(1:length(peakLabel)/2);

clear g
figure('Position',[100 100 1400 600]);
for p =1:4
%Create a scatter plot
meanpDiff = nanmean(diffLinPeakVals(strcmp(diffPeakLabel, sprintf('Pk%d',p)))); %compute difference between monocular and binocular condition

g(1,2*(p-1)+1)=gramm('x',diffLinPeakVals,'subset', strcmp(diffPeakLabel,sprintf('Pk%d',p)));
g(1,2*(p-1)+1).set_names('x','Spiking activity difference (Normalized)','column','');
g(1,2*(p-1)+1).geom_jitter('width',0,'height',0.2); %Scatter plot
g(1,2*(p-1)+1).geom_vline('xintercept', meanpDiff, 'style', '-k');

%g(1,2*(p-1)+1).axe_property('Ygrid','on', 'ylim',[0.3 1.7],'YTickLabel','','YTick',''); 
g(1,2*(p-1)+1).axe_property('ylim',[0.3 1.7]); 
g(1,2*(p-1)+1).no_legend();

%Create y data histogram on the right
g(1,2*(p-1)+2)=gramm('x',diffLinPeakVals,'subset', strcmp(diffPeakLabel,sprintf('Pk%d',p)));
%g(2,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
%    'legend',false,...
%    'margin_height',[0.1 0.02],...
%    'margin_width',[0.02 0.05],...
%    'redraw',false);
g(1,2*(p-1)+2).geom_vline('xintercept', meanpDiff, 'style', '-k');

g(1,2*(p-1)+2).set_names('x','');
g(1,2*(p-1)+2).stat_bin('geom','overlaid_bar', 'nbins', round(std(diffLinPeakVals(strcmp(diffPeakLabel, sprintf('Pk%d',p))),'omitnan')*100)); %histogram
%g(1,2*(p-1)+2).axe_property('xlim',[0 0.5],'ylim',[-2 15],'XTickLabel','','XTick',''); 
g(1,2*(p-1)+2).no_legend();

end
%Set global axe properties
g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.axe_property('xlim',[-0.75 0.45]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g.coord_flip();
g.set_title('Population peak response of monocular - binocular conditions differences');
g.set_color_options('map','d3_10');
%g.set_color_options('map', [251/255 154/255 153/255;160/255 160/255 160/255]);
g.draw();
%set(findobj(gcf, 'type','axes'), 'Visible','off')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_hist_mono_bino_allcells_peakdiffs');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%}
%%
%2) Only plot significantly modulated units that show adaptation in the
%monocular condition
%% Box plot with jitter of significant binocularly modulated units that show suppressive adaptation
%select data of those specific units  using the table created in "lmer_peaks_binocular_adaptation.R" script
%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
filename = 'summary_table_pvalues_normmono_meanpks_mono_bino_dunnett';
T = readtable(strcat(newdatadir, filename, '.csv'));

T(ismember(T.Wpvalue,'NA'),:)=[]; %remove rows with NAs
%T(T.c4>T.c1,:)=[];
T.Wpvalue = str2double(T.Wpvalue);
goodrows = T.Wpvalue < 0.05; %keep units with significant binocular modulation (turns out they also all show adaptation in either of the mono or binocular condition)
sigBinoT = T(goodrows,:);

%modify Var1
sigBinoT.Var1 = repmat(1:4, 1,2)';
stackSigBinoT = stack(sigBinoT,{'Pk1','Pk2','Pk3','Pk4'},'NewDataVariableName','PeakResp');

%plot
clear g
figure('Position',[100 100 1000 600]);
g(1,1)=gramm('x',stackSigBinoT.PeakResp_Indicator ,'y',str2double(stackSigBinoT.PeakResp),'color',stackSigBinoT.Condition);
%g(1,1)=gramm('x',str2double(stackSigBinoT.PeakResp),'color',stackSigBinoT.Condition);

g(1,1).set_names('x','Spiking activity (Normalized)', 'column', '');
%g(1,1).facet_grid([],grp2idx(stackSigBinoT.PeakResp_Indicator)); %Provide facets
%g(1,1).set_color_options('map', [251/255 154/255 153/255;160/255 160/255 160/255]);
g(1,1).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]); %We make it light grey

g(1,1).geom_jitter('width',0,'height',0); %Scatter plot
g(1,1).stat_boxplot();
%g(1,1).stat_density();
%g(1,1).stat_violin('fill','transparent');
%g(1,1).axe_property('xlim',[0.4 1.8],'ylim',[-2 15]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g.set_title({'Peak responses in the binocular and monocular conditions of','binocularly modulated adapting units'});
%Set global axe properties
g.axe_property('TickDir','out');
%g.coord_flip();
g.draw();
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\box_mono_bino_modadaptcells_peaks_purpOr_dunnett');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% 
%% Box plot with jitter of significant binocularly modulated units that show suppressive adaptation
%select data of those specific units  using the table created in "lmer_peaks_binocular_adaptation.R"

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
filename = 'summary_table_pvalues_normmono_meanpks_mono_bino';
T = readtable(strcat(newdatadir, filename, '.csv'));

T(ismember(T.Wpvalue,'NA'),:)=[]; %remove rows with NAs
%T(T.c4>T.c1,:)=[];
T.Wpvalue = str2double(T.Wpvalue);


stackT = stack(T,{'Pk1','Pk2','Pk3','Pk4'},'NewDataVariableName','PeakResp');

%plot
clear g
figure('Position',[100 100 1000 600]);
g(1,1)=gramm('x',stackT.PeakResp_Indicator ,'y',str2double(stackT.PeakResp),'color',stackT.Condition);
%g(1,1)=gramm('x',str2double(stackSigBinoT.PeakResp),'color',stackSigBinoT.Condition);
g(1,1).set_names('y','Spiking activity (Normalized)', 'column', '');
%g(1,1).facet_grid([],grp2idx(stackSigBinoT.PeakResp_Indicator)); %Provide facets
g(1,1).set_color_options('map', [251/255 154/255 153/255;160/255 160/255 160/255]);
%g(1,1).set_color_options('chroma',0,'lightness',75); %We make it light grey
g(1,1).geom_jitter('width',0,'height',0); %Scatter plot
g(1,1).stat_boxplot();
%g(1,1).stat_density();
%g(1,1).stat_violin('fill','transparent');
%g(1,1).axe_property('xlim',[0.4 1.8],'ylim',[-2 15]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g.set_title({'Peak responses in the binocular and monocular conditions of','single units population'});
%Set global axe properties
g.axe_property('TickDir','out');
%g.coord_flip();
g.draw();
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\box_mono_bino_allcells_peaks');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));



%% Quick stats to compare the mean distributions of each peak between conditions

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']; 
peak_aligned_trials = load(channelfilename);

%Normalize data differently so Pk1 mono also has a normal distribution
filenames = fieldnames(peak_aligned_trials.peak_aligned_trials);
bins = [1,6];
norm_aligned_resps = nan(250, 4, 2,length(filenames) );
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
     %normalize in relation to all mean peaks in the monocular condition
         norm_aligned_resps(:,:,b,i) = (mean_peaks(:,:,b,i) - min(mean_peaks(:,:,1,i),[], 'all'))./(max(mean_peaks(:,:,1,i), [], 'all') - min(mean_peaks(:,:,1,i),[], 'all'));
         %norm_aligned_resps(:,:,b,i) = (mean_peaks(:,:,b,i))./(max(mean_peaks(:,1,1,i), [], 'all'));
        end
    end
end
%store norm peak values
 peaks = nan(4,2,length(filenames));
for i = 1:length(filenames)
    for b = 1:2
        for p = 1:4
            peaks(p,b,i) = max(norm_aligned_resps(:,p,b,i));
        end
    end
end
%%find a way to melt matrix in
%%p1p1p1p1....p2p2p2p2....p3p3p3p3....p4p4p4p4p4 fashion
peakvals = [squeeze(peaks(1,:,:))';squeeze(peaks(2,:,:))';squeeze(peaks(3,:,:))';squeeze(peaks(4,:,:))'];
%peakvals = reshape(peaks, [length(peaks(:,1,1))*length(peaks(1,1,:)), length(peaks(1,:,1))]);
linPeakVals = reshape(peakvals, [2*length(peakvals(:,1)),1]); 
condition = [repmat({'Monocular'},length(peakvals(:,1)),1); repmat({'Binocular'},length(peakvals(:,1)),1)];
unit = repmat(1:length(peaks(1,1,:)),1,8)';
peakLabel = repmat([repmat({'Pk1'}, length(peaks(1,1,:)),1);repmat({'Pk2'}, length(peaks(1,1,:)),1);repmat({'Pk3'}, length(peaks(1,1,:)),1);repmat({'Pk4'}, length(peaks(1,1,:)),1)],2,1);

for p =1:4
    meansMono =linPeakVals(strcmp(condition, 'Monocular')& strcmp(peakLabel,sprintf('Pk%d',p)));
    meansBino =linPeakVals(strcmp(condition, 'Binocular')& strcmp(peakLabel,sprintf('Pk%d',p))); 
    [ttestRes(p), pval(p),~,stats(p).stats] = ttest(meansMono,meansBino);
    WtestRes(p) = signrank(meansMono,meansBino);
end

%% Correlational stats between peaks

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

%compute correlation coeff Pearson's r between differences btween peaks

%1) Compute difference between peaks
peaks_diff = nan(3,2,length(filenames));
for i = 1:length(filenames)
    for b = 1:2
        for p = 2:4
            peaks_diff(p-1,b,i) = (peaks(1,b,i)- peaks(p,b,i)); %/peaks(1,b,i);
        end
    end
end
peaks_diff = peaks_diff(:,:,~isnan(peaks_diff(1,1,:)));


%% Plot the correlation trends between subsequent peaks (Pk1-Pk3 vs Pk1-Pk4)(Pk1-Pk2 vs Pk1-Pk3)
%1) Plot Pk1-Pk2 vs Pk1-Pk3 + overlay Pk1-Pk3 vs Pk1-Pk4
%nlines = 7;
%cmap =flip(cbrewer2('Blues', nlines));
%colormap(cmap);
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =flip(cbrewer2('BuPu', nlines));

%cmap = flip(cmaps(2).map) ;
colormap(cmaps(2).map);

colmaps = [cmaps(1).map(4,:);cmaps(2).map(3,:)]; %invert order of colors since we plot monocular condition first
titles = {'Monocular'; 'Binocular'};
figure('Position',[100 100 800 400],'Color',[1 1 1]);   
for b =1:2
ax = subplot(1,2,b);
    hold on
    % Monocular condition
    x1 = squeeze(peaks_diff(1,b,:)); %Pk1-Pk2
    y1 =  squeeze(peaks_diff(2,b,:)); %Pk1-Pk3

    %linear regression
    coeffs1 = polyfit(x1(isfinite(x1) & isfinite(y1)),y1(isfinite(x1) & isfinite(y1)),1);
    f1 = polyval(coeffs1,x1);
    plot(x1, y1,'o',x1, f1,'-','Color',[160/255 160/255 160/255],'MarkerSize',2, 'MarkerFaceColor',[160/255 160/255 160/255],'linewidth',2)
    %xlim([0 10])
    %ylim([0 4.5])
    text(max(x1)/1.3,max(y1)/20, sprintf('y1 = %.2f + %.2f*x', round(coeffs1(2),2), round(coeffs1(1),2)))

    x2 = squeeze(peaks_diff(2,b,:));
    y2 = squeeze(peaks_diff(3,b,:));
    
    % Keep the same color for the statistics
    coeffs2 = polyfit(x2(isfinite(x2) & isfinite(y2)),y2(isfinite(x2) & isfinite(y2)),1);
    f2 = polyval(coeffs2,x2);
    plot(x2, y2,'o',x2, f2,'-','Color',colmaps(b,:),'MarkerSize',2, 'MarkerFaceColor',colmaps(b,:),'linewidth',2)
    %xlim([0 10])
    %ylim([-0.7 0.8])
    text(max(x2)/1.3,max(y2)/20, sprintf('y2 = %.2f + %.2f*x', round(coeffs2(2),2), round(coeffs2(1),2)))

    hold on
    set(gca, 'linewidth',2)
    set(gca,'box','off')
    legend('','Pk1-Pk3 = f(Pk1-Pk2)','','Pk1-Pk4=f(Pk1-Pk3)')
    title(titles(b))
end
 
 plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\trends_across_peaks_normalized');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));
 
%% Need to make statistical tests on slopes vs 0 for above linear fits
xtest = x2;
ytest = y2;
linreg = fitlm(xtest,ytest);
[linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linreg); %r =numerator degrees of freedom
length(xtest)
linPvalue(1,1)
F(1,1)

%% Assess difference betwen slopes of linear regressions performed above (previous section) 
%% between monocular and binocular stimulation

xm = squeeze(peaks_diff(1,1,:)); %Pk1-Pk2 mono
ym =  squeeze(peaks_diff(2,1,:)); %Pk1-Pk3 mono

xb = squeeze(peaks_diff(1,2,:)); %Pk1-Pk2 bino
yb =  squeeze(peaks_diff(2,2,:)); %Pk1-Pk3 bino

linm = fitlm(xm,ym);
linb = fitlm(xb,yb);
[linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linm,[0 1],linb.Coefficients{2,1}); %[0 1] = origin*0+1*slopemono = slopebino?
linPvalue(1,1)
F(1,1)

xm = squeeze(peaks_diff(2,1,:)); %Pk1-Pk3 mono
ym =  squeeze(peaks_diff(3,1,:)); %Pk1-Pk4 mono

xb = squeeze(peaks_diff(2,2,:)); %Pk1-Pk3 bino
yb =  squeeze(peaks_diff(3,2,:)); %Pk1-Pk4 bino

linm = fitlm(xm,ym);
linb = fitlm(xb,yb);
[linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linm,[0 1],linb.Coefficients{2,1}); %[0 1] = origin*0+1*slopemono = slopebino?
linPvalue(1,1)
F(1,1)


%% 2)Compute Pearson's r correlations

corrs = nan(2,2);

    for b = 1:2
        corrs(1,b) = corr(squeeze(peaks_diff(1,b,:)), squeeze(peaks_diff(2,b,:))); %compare pk1-pk2 vs pk1-pk3
        corrs(2,b) = corr(squeeze(peaks_diff(2,b,:)), squeeze(peaks_diff(3,b,:))); %compare pk1-pk3 vs pk1-pk4
    end



