%This script is developped following trial selection with
%"BinocularAdaptationTrialSelection.m" and adaptation and binocular interaction analysis 
% with "mua_lmer_peaks_binocular_adaptation.R"

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
         norm_aligned_resps(:,:,b,i) = (mean_peaks(:,:,b,i) - min(mean_peaks(:,1,1,i),[], 'all'))./(max(mean_peaks(:,1,1,i), [], 'all') - min(mean_peaks(:,1,1,i),[], 'all'));
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
%g(1,1).set_color_options('chroma',0,'lightness',75); %We make it light grey
g(1,1).set_title('Population peak responses in the binocular and monocular conditions');
g(1,1).axe_property('xlim',[0.4 1.8],'ylim',[-2 15]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten

%Set global axe properties
g.axe_property('TickDir','out');
g.coord_flip();
%g.set_title({'Adaptation index distribution across all cells in the monocular and binocular conditions'});
g.draw();

%set(f,'position',get(f,'position').*[1 1 1.15 1])
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\hist_mono_bino_allcells_peaks');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% Plot jitter scatter plot + horizontal histogram on the side for each peak
%{

clear g
figure('Position',[100 100 1400 600]);
for p =1:4
%Create a scatter plot
g(1,2*(n-1)+1)=gramm('x',linPeakVals,'color',condition);
g(1,2*(n-1)+1).set_names('x','Spiking activity (Normalized)','column','');
g(1,2*(n-1)+1).geom_jitter('width',0.4,'height',0.3); %Scatter plot
g(1,2*(n-1)+1).axe_property('Ygrid','on', 'ylim',[0.3 1.7],'YTickLabel','','YTick',''); 
g(1,2*(n-1)+1).coord_flip();
g(1,2*(n-1)+1).draw();

%Create y data histogram on the right
g(1,2*(n-1)+2)=gramm('x',linPeakVals,'color',condition);
%g(2,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
%    'legend',false,...
%    'margin_height',[0.1 0.02],...
%    'margin_width',[0.02 0.05],...
%    'redraw',false);
g(1,2*(n-1)+2).set_names('x','');
g(1,2*(n-1)+2).stat_bin('geom','stacked_bar','fill','transparent','nbins',15); %histogram
g(1,2*(n-1)+2).coord_flip();
g(1,2*(n-1)+2).draw();

end

%Set global axe properties
g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.set_title('Population peak responses in the binocular and monocular conditions');
g.set_color_options('map','d3_10');
g.draw();


%% violin plot with jitter

clear g
figure('Position',[100 100 1400 600]);
g(1,1)=gramm('x',unit ,'y',linPeakVals,'color',condition);
%g(1,1).facet_grid([],peakLabel); %Provide facets
%g(1,1).geom_jitter('width',0.4,'height',0.3); %Scatter plot
g(1,1).stat_violin('fill','transparent');
g(1,1).set_names('y','Spiking activity (Normalized)','column','');
%g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.set_title('Population peak responses in the binocular and monocular conditions');
%g.set_color_options('map','d3_10');
g.draw();



%%
%2) Only plot significantly modulated units that show adaptation in the
%monocular condition

%}