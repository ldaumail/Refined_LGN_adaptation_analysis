%This script was developped to analyze the noise present in single units data

%get filenames where the data is located
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\s_potentials_analysis\analysis\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;

unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);

%get binary data
[binSpkTrials,NoFiltMultiContSUA, peakLocs] = binTrialSelection(unitsData, filenames);

allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\binary_trials_data_06022021';
save(strcat(allfilename, '.mat'), 'binSpkTrials');
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\origin_trials_data_06022021';
save(strcat(allfilename, '.mat'), 'NoFiltMultiContSUA');
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\peak_locs_data_06022021';
save(strcat(allfilename, '.mat'), 'peakLocs');

%{
%to load older data
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
peakLocs = load(strcat(datadir, 'all_locs_data_95CI_05022021'));
peakLocs = peakLocs.peakLocs;
NoFiltMultiContSUA = load(strcat(datadir,'NoFiltMultiContSUA_05022021'));
NoFiltMultiContSUA = NoFiltMultiContSUA.NoFiltMultiContSUA;
%}


%Load peakLocs, NoFiltMultiContSUA, binSpkTrials
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
peakLocs = load(strcat(datadir, 'peak_locs_data_06022021'));
peakLocs = peakLocs.peakLocs;
NoFiltMultiContSUA = load(strcat(datadir,'origin_trials_data_06022021'));
NoFiltMultiContSUA = NoFiltMultiContSUA.NoFiltMultiContSUA;
binSpkTrials = load(strcat(datadir,'binary_trials_data_06022021'));
binSpkTrials = binSpkTrials.binSpkTrials;

%align binary trials to peak locations for each peak
[wind_peak_vals, peak_aligned_trials] = binPeakTrigResps(peakLocs, NoFiltMultiContSUA, binSpkTrials);
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\peak_aligned_binary_trials_data_06022021';
save(strcat(allfilename, '.mat'), 'peak_aligned_trials');
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\binary_peak_winds30_06022021';
save(strcat(allfilename, '.mat'), 'wind_peak_vals');

%load the binary peak values data
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
wind_peak_vals = load(strcat(datadir, 'binary_peak_winds30_06022021'));
wind_peak_vals = wind_peak_vals.wind_peak_vals;


filenames = fieldnames(wind_peak_vals);
bins = [1,6];

mean_winds = nan(4, 2,length(filenames) ); %Mean spike # per window for each peak in each unit
var_winds = nan(4, 2,length(filenames) ); %Variance of spike counts for each peak in each unit
fanoF = nan(4, 2,length(filenames) );


for i = 1: length(filenames)
    filename = filenames{i};
    if length(fieldnames(wind_peak_vals.(filename))) == 2
        for b = 1:2
            binN = sprintf('bin%d',bins(b));
             %compute mean peak responses
             for p = 1:4
                 
                 count = sum(squeeze(wind_peak_vals.(filename).(binN)(p,:,:)),1);
                 mean_winds(p,b,i) = sum(count)/size(squeeze(wind_peak_vals.(filename).(binN)(p,:,:)),2);
                 fanoF(p,b,i) = var(count)/mean_winds(p,b,i);
             end

        end
    end
end

%%find a way to melt matrix in
%%p1p1p1p1....p2p2p2p2....p3p3p3p3....p4p4p4p4p4 fashion
peakvals = [squeeze(mean_winds(1,:,:))';squeeze(mean_winds(2,:,:))';squeeze(mean_winds(3,:,:))';squeeze(mean_winds(4,:,:))'];
%peakvals = reshape(peaks, [length(peaks(:,1,1))*length(peaks(1,1,:)), length(peaks(1,:,1))]);
linPeakVals = reshape(peakvals, [2*length(peakvals(:,1)),1]); 
condition = [repmat({'Monocular'},length(peakvals(:,1)),1); repmat({'Binocular'},length(peakvals(:,1)),1)];
unit = repmat(1:length(mean_winds(1,1,:)),1,8)';
peakLabel = repmat([repmat({'Pk1'}, length(mean_winds(1,1,:)),1);repmat({'Pk2'}, length(mean_winds(1,1,:)),1);repmat({'Pk3'}, length(mean_winds(1,1,:)),1);repmat({'Pk4'}, length(mean_winds(1,1,:)),1)],2,1);

%fano vals
fanovals = [squeeze(fanoF(1,:,:))';squeeze(fanoF(2,:,:))';squeeze(fanoF(3,:,:))';squeeze(fanoF(4,:,:))'];
linFanoVals = reshape(fanovals, [2*length(fanovals(:,1)),1]); 


%compute fano factor for each peak
%Fano factor = variance/ mean
clear g
figure('Position',[100 100 1400 600]);
for p =1:4
%Create a scatter plot
%meanpMono(p) = nanmean(linPeakVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p))));
%meanpBino(p) = nanmean(linPeakVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p))));

%varpMono(p) = var(linPeakVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p))),'omitnan');
%varpBino(p) = var(linPeakVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p))),'omitnan');

fanofMono(p) = nanmean(linFanoVals( strcmp(condition, 'Monocular') & strcmp(peakLabel,sprintf('Pk%d',p))));
fanofBino(p) = nanmean(linFanoVals( strcmp(condition, 'Binocular') & strcmp(peakLabel,sprintf('Pk%d',p))));

g(1,2*(p-1)+1)=gramm('x',linFanoVals,'color',condition,'subset', strcmp(peakLabel,sprintf('Pk%d',p)));
g(1,2*(p-1)+1).set_names('x','Fano Factor','column','');
g(1,2*(p-1)+1).geom_jitter('width',0,'height',0.2); %Scatter plot
g(1,2*(p-1)+1).geom_vline('xintercept', fanofMono(p), 'style', '-k');
g(1,2*(p-1)+1).geom_vline('xintercept', fanofBino(p), 'style', '-p');

%g(1,2*(p-1)+1).axe_property('Ygrid','on', 'ylim',[0.3 1.7],'YTickLabel','','YTick',''); 
g(1,2*(p-1)+1).axe_property('ylim', [0.6 1.3]); 
g(1,2*(p-1)+1).no_legend();

%Create y data histogram on the right
g(1,2*(p-1)+2)=gramm('x',linFanoVals,'color',condition, 'subset', strcmp(peakLabel,sprintf('Pk%d',p)));
%g(2,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
%    'legend',false,...
%    'margin_height',[0.1 0.02],...
%    'margin_width',[0.02 0.05],...
%    'redraw',false);
g(1,2*(p-1)+2).geom_vline('xintercept', fanofMono(p), 'style', '-k');
g(1,2*(p-1)+2).geom_vline('xintercept', fanofBino(p), 'style', '-p');

g(1,2*(p-1)+2).set_names('x','');
g(1,2*(p-1)+2).stat_bin('geom','overlaid_bar'); %histogram
%g(1,2*(p-1)+2).axe_property('xlim',[0.3 1.7],'ylim',[-2 15],'XTickLabel','','XTick',''); 
g(1,2*(p-1)+2).no_legend();

end
%Set global axe properties
g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.axe_property('xlim',[0 9.2]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g.coord_flip();
g.set_title('Population peak fano factors in the binocular and monocular conditions');
%g.set_color_options('map','d3_10');
g.set_color_options('map', [251/255 154/255 153/255;160/255 160/255 160/255]);
g.draw();
%set(findobj(gcf, 'type','axes'), 'Visible','off')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_hist_mono_bino_allcells_fanof');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% Fano Factor with a sliding window 
%Load peakLocs, NoFiltMultiContSUA, binSpkTrials
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
peakLocs = load(strcat(datadir, 'peak_locs_data_06022021'));
peakLocs = peakLocs.peakLocs;
NoFiltMultiContSUA = load(strcat(datadir,'origin_trials_data_06022021'));
NoFiltMultiContSUA = NoFiltMultiContSUA.NoFiltMultiContSUA;
binSpkTrials = load(strcat(datadir,'binary_trials_data_06022021'));
binSpkTrials = binSpkTrials.binSpkTrials;
wsz = [20,30,50,70,100]; %window size
%align binary trials to peak locations for each peak
[slide_win_fanof, peak_aligned_trials] = peakTrigFano(peakLocs, NoFiltMultiContSUA, binSpkTrials,wsz);

%set fanof in matrix format
filenames = fieldnames(slide_win_fanof);
bins = [1,6];
all_fanofs = nan(length(-125:10:125-wsz(1)),5,2,length(wsz),length(filenames)); %5 = 4 pks +1 resting state
len = nan(length(wsz),1);
for i =1:length(filenames)
    filename = filenames{i};
    if length(fieldnames(slide_win_fanof.(filename).fanof)) == 2 
        for b =1:2
            binN = sprintf('bin%d',bins(b));
            %fill up matrix with FF of peaks
            for p =1:4
                for w =1:length(wsz)
                    windSz = sprintf('wsz%d',wsz(w));
                    len(w) = length(slide_win_fanof.(filename).fanof.(binN).(windSz).peaks(p,:));
                    all_fanofs(1:len(w),p+1,b,w,i) =  slide_win_fanof.(filename).fanof.(binN).(windSz).peaks(p,:);
                end
            end
            %fill up first  matrix column with FF of resting state 
              for w =1:length(wsz)
                    windSz = sprintf('wsz%d',wsz(w));
                    len(w) = length(slide_win_fanof.(filename).fanof.(binN).(windSz).rs(:));
                    all_fanofs(1:len(w),1,b,w,i) =  slide_win_fanof.(filename).fanof.(binN).(windSz).rs(:);
              end
        end
    end
    
    
end

all_fanofs = squeeze(nanmean(all_fanofs(:,:,:,:,:),3)); %computing the mean of between mono and binocular to reduce error size
col(1,:) =[194/255 165/255 207/255] ; %--purple
col(2,:) = [253/255 174/255 97/255]; % -- orange
%col(3,:) = [166/255 219/255 160/255]; % -- green


figure();

mean_fanof =nanmean(all_fanofs,4);
for p = 1:5
    h =subplot(1,5,p);
    for w =1:length(wsz)
        %ci_low = mean_fanof(:,p,w) - 1.96*std(squeeze(all_fanofs(:,p,w,:)),0,2, 'omitnan')./sqrt(length(find(~isnan(all_fanofs(1,p,w,:)))));
        %ci_high = mean_fanof(:,p,w) + 1.96*std(squeeze(all_fanofs(:,p,w,:)),0,2, 'omitnan')./sqrt(length(find(~isnan(all_fanofs(1,p,w,:)))));
        
        plot(-125+wsz(w)/2:10:125-wsz(w)/2,mean_fanof(1:len(w),p,w),'LineWidth',2)%'Color',[40/255 40/255 40/255] )
        hold on
        %h1= ciplot(ci_low, ci_high,[-125+wsz/2:10:125-wsz/2],col(1,:),0.5);
        %set(h1, 'edgecolor','none')
    end
    ylim([0 3])
    
    
    % K norm ylim([-0.02 1.1])
    xlim([-125 125])
    set(gca, 'linewidth',2)
    set(gca,'box','off')
    %set(gca, 'linewidth',2)
    %hold on
    %plot([0 0], ylim,'k')
    if p ==1
        ylabel({'\fontsize{14}Fano Factor'})
    end
    if p >1
        ax1 = gca;
        ax1.YAxis.Visible = 'off';
    end
    if p ==5
        
        legend('wsz = 20ms', 'wsz = 30ms','wsz = 50ms', 'wsz = 70ms', 'wsz = 100ms')
    end
end


%brown [165/255 42/255 42/255]
%currfig = gcf;
%title(currfig.Children(end),'Mean Sliding Window (50 ms/ 10 ms step) Fano Factor', 'FontSize', 20)
sgtitle('Mean Sliding Window (10 ms steps) Fano Factor', 'FontSize', 20)

% xlabel('\fontsize{14}Resolution (ms)')


set(gcf,'Units','inches')
set(gcf,'position',[1 1 15 11])


plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\',sprintf('all_sliding_window_mean_fanof'));
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));


%% %%%%%%%% raster plot of example single unit %%%%
 %raster plot of trials binary spikes
 %for i =1:length(filenames)
 
 
 for i =40:length(filenames)
     filename = filenames{i};
     
     figure();
     
     for p = 1:4
         h =subplot(1,4,p);
         
         pkN = sprintf('pk%d',p);
         if isfield(peak_aligned_trials.(filename).binarySUA,'bin1')
             for tr =1:length(peak_aligned_trials.(filename).binarySUA.bin1.(pkN)(1,:))
                 spikeTimes =find(peak_aligned_trials.(filename).binarySUA.bin1.(pkN)(:,tr)==1)';
                 x = repmat(spikeTimes,3,1);
                 y = nan(size(x));
                 
                 if ~isempty(y)
                     y(1,:) = tr-1;
                     y(2,:) = tr;
                 end
                 plot(x-500,y,'Color','k')
                 hold on
                 
             end
             set(gca, 'linewidth',2)
             set(gca,'box','off')
             xlabel('Time (ms)')
             ylabel('Trial number')
             if p >1
                 ax1 = gca;
                 ax1.YAxis.Visible = 'off';
             end
         end
     end
     
     sgtitle({'Binary spike data of an example single unit','aligned to each peak'})
     set(gcf,'position',get(gcf,'position').*[1 1 1.15 1])
     plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\', sprintf('raster_plots_single_unit_example%d',i));
     saveas(gcf,strcat(plotdir, '.png'));
     saveas(gcf,strcat(plotdir, '.svg'));
     
     
     
 end


    %%%%%%   %%%%%%%%%
%% Fano Factor as a function of the mean spike count

%Load peakLocs, NoFiltMultiContSUA, binSpkTrials
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
peakLocs = load(strcat(datadir, 'peak_locs_data_06022021'));
peakLocs = peakLocs.peakLocs;
NoFiltMultiContSUA = load(strcat(datadir,'origin_trials_data_06022021'));
NoFiltMultiContSUA = NoFiltMultiContSUA.NoFiltMultiContSUA;
binSpkTrials = load(strcat(datadir,'binary_trials_data_06022021'));
binSpkTrials = binSpkTrials.binSpkTrials;

wsz = [20,30,50,70,100]; %window size

%align binary trials to peak locations for each peak
[slide_win_fanof, peak_aligned_trials] = peakTrigFano(peakLocs, NoFiltMultiContSUA, binSpkTrials, wsz);

%set fanof in matrix format and get peak values
filenames = fieldnames(slide_win_fanof); 
bins = [1,6];
all_fanofs = nan(length(-125:10:125-wsz(1)), 4,2,length(wsz),length(filenames));
peak_vals = nan(4,2,length(wsz),length(filenames));
var_vals = nan(4,2,length(wsz),length(filenames));
peak_fanofs = nan(4,2,length(wsz),length(filenames));

for i =1:length(filenames)
    filename = filenames{i};
    if length(fieldnames(slide_win_fanof.(filename).fanof)) == 2
        for b =1:2
            binN = sprintf('bin%d',bins(b));
            for p =1:4
                for w =1:length(wsz)
                    windSz = sprintf('wsz%d',wsz(w));
                    len(w) = length(slide_win_fanof.(filename).fanof.(binN).(windSz)(p,:));
                    all_fanofs(1:len(w),p,b,w,i) =  slide_win_fanof.(filename).fanof.(binN).(windSz)(p,:);
                    peak_fanofs(p,b,w,i) =  slide_win_fanof.(filename).fanof.(binN).(windSz)(p,round(len(w)/2));
                    peak_vals(p,b,w,i) = slide_win_fanof.(filename).meanspkc.(binN).(windSz)(p,round(len(w)/2));
                    var_vals(p,b,w,i) = slide_win_fanof.(filename).varspkc.(binN).(windSz)(p,round(len(w)/2));
                    
                    % pkN = sprintf('pk%d',p);
                    % for tr = 1:length(peak_aligned_trials.(filename).convolvedSUA.(binN)(1,:))
                    %    peak_vals(p,b,i) = nanmean(max(peak_aligned_trials.(filename).convolvedSUA.(binN).(pkN),[],1));
                    % end
                end
            end
        end
    end
end

%merge mono and bino since they are not different
mpeak_fanofs = squeeze(nanmean(peak_fanofs,2));
mpeak_vals = squeeze(nanmean(peak_vals,2));    
mvar_vals = squeeze(nanmean(var_vals,2));              

%

col(1,:) =[194/255 165/255 207/255] ; %--purple
col(2,:) = [253/255 174/255 97/255]; % -- orange


figure('Position',[100 100 1400 600]);

for p =1:4
    h =subplot(1,4,p);
    
        x =squeeze(mpeak_vals(p,3,:));
        y =squeeze(mpeak_fanofs(p,3,:));
        %scatter(x, y, 15, col(b,:),'filled')
        %hold on
        %hline = refline(0, nanmean(squeeze(peak_fanofs(p,b,:))));
        %hline.Color = col(b,:);
        %hline.LineWidth = 2;
        coeffs = polyfit(x(isfinite(x) & isfinite(y)),y(isfinite(x) & isfinite(y)),1); 
        f = polyval(coeffs,x); 
        plot(x, y,'o',x, f,'-','Color',col(1,:))
        hold on
        %format short
    text(max(x)/12,max(y)/2, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)));
    set(gca, 'linewidth',2)
    set(gca,'box','off')
    xlim([0 4])
    ylim([0 5])
    if p ==1 
    xlabel('Mean Spike Count in a 50ms Window')
    ylabel('Fano Factor')
    end
    title(sprintf('Peak %d', p))
    if p >1
        ax1 = gca;
       ax1.YAxis.Visible = 'off';
       
    end
end
sgtitle('Fano Factor as a function of mean spike count across trials of single units')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\scatter_plot_fanof_vs_spkcmean_regression');
     saveas(gcf,strcat(plotdir, '.png'));
     saveas(gcf,strcat(plotdir, '.svg'));

     

%% fano factor vs variance
figure('Position',[100 100 1400 600]);

for p =1:4
    h =subplot(1,4,p);
    
        x =squeeze(mvar_vals(p,3,:));
        y =squeeze(mpeak_fanofs(p,3,:));
        %scatter(x, y, 15, col(b,:),'filled')
        %hold on
        %hline = refline(0, nanmean(squeeze(peak_fanofs(p,b,:))));
        %hline.Color = col(b,:);
        %hline.LineWidth = 2;
        coeffs = polyfit(x(isfinite(x) & isfinite(y)),y(isfinite(x) & isfinite(y)),1); 
        f = polyval(coeffs,x); 
        plot(x, y,'o',x, f,'-','Color',col(2,:))
        hold on
        %format short
    text(max(x)/12,max(y)/2, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)));
    set(gca, 'linewidth',2)
    set(gca,'box','off')
    xlim([0 4])
    ylim([0 5])
    if p ==1 
    xlabel('Spike Count Variance in a 50ms Window')
    ylabel('Fano Factor')
    end
    title(sprintf('Peak %d', p))
    if p >1
        ax1 = gca;
       ax1.YAxis.Visible = 'off';
       
    end
end
sgtitle('Fano Factor as a function of spike count variance across trials of single units')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\scatter_plot_fanof_vs_spkcvar_regression');
     saveas(gcf,strcat(plotdir, '.png'));
     saveas(gcf,strcat(plotdir, '.svg'));

%% Mean vs variance

figure('Position',[100 100 1400 600]);

for p =1:4
    h =subplot(1,4,p);
    
    x =squeeze(mpeak_vals(p,3,:));
    y =squeeze(mvar_vals(p,3,:));
    %scatter(x, y, 15, col(b,:),'filled')
    %hold on
    %hline = refline(0, nanmean(squeeze(peak_fanofs(p,b,:))));
    %hline.Color = col(b,:);
    %hline.LineWidth = 2;
    coeffs = polyfit(x(isfinite(x) & isfinite(y)),y(isfinite(x) & isfinite(y)),1);
    f = polyval(coeffs,x);
    plot(x, y,'o',x, f,'-','Color',col(1,:))
    hold on
    %format short
    text(max(x)/12,max(y)/2, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)));
    set(gca, 'linewidth',2)
    set(gca,'box','off')
    xlim([0 4])
    ylim([0 4])
    if p ==1
        xlabel('Mean Spike Count in a 50ms Window')
        ylabel('Spike Count Variance')
    end
    title(sprintf('Peak %d', p))
    if p >1
        ax1 = gca;
        ax1.YAxis.Visible = 'off';
        
    end
end
sgtitle('Fano Factor as a function of spike count variance across trials of single units')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\scatter_plot_spkcvar_vs_spkcmean_regression');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));


%% Variance only


figure('Position',[100 100 1400 800]);

for p =1:4
    h =subplot(1,4,p);
    y =squeeze(mvar_vals(p,3,:));
    boxplot(y,'Color',col(1,:))
    hold on
    set(gca,'box','off')
    ylim([0 4])
    set(gca, 'linewidth',2)
    set(gca,'box','off')
    
        xlabel(sprintf('Peak %d', p))
    if p ==1
        ylabel('Spike Count Variance')
    end
    title(sprintf('Peak %d', p))
    if p >1
        ax1 = gca;
        ax1.YAxis.Visible = 'off';
        
    end
    plot(nanmean(y), 'dg')
    hold off
    text(1,nanmean(y), sprintf('<y> = %.2f', nanmean(y)));
  
end
sgtitle('Spike count variance of single units across trials')
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\boxplot_spkcvar');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

