%This script was developped to analyze the variability quenching present in single units data
%last edited by Loic Daumail on 08-18-2021

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

%%Prepare binary data (no need to repeat that if already done once)
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

%% compute Fano Factor and start plotting data
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
peak_vals = nan(4,2,length(wsz),length(filenames)); %get mean peak spike counts
var_vals = nan(5,2,length(wsz),length(filenames)); %get trial-to-trial variance

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
                    var_vals(p+1,b,w,i) = slide_win_fanof.(filename).varspkc.(binN).(windSz).peaks(p,round(len(w)/2));
                    peak_vals(p+1,b,w,i) = slide_win_fanof.(filename).meanspkc.(binN).(windSz).peaks(p,round(len(w)/2));
                   
                end
            end
            %fill up first  matrix column with FF of resting state 
              for w =1:length(wsz)
                    windSz = sprintf('wsz%d',wsz(w));
                    len(w) = length(slide_win_fanof.(filename).fanof.(binN).(windSz).rs(:));
                    all_fanofs(1:len(w),1,b,w,i) =  slide_win_fanof.(filename).fanof.(binN).(windSz).rs(:);
                    var_vals(1,b,w,i) = slide_win_fanof.(filename).varspkc.(binN).(windSz).rs(round(len(w)/2));
                    peak_vals(1,b,w,i) = slide_win_fanof.(filename).meanspkc.(binN).(windSz).rs(round(len(w)/2));
                  
              end
        end
    end
    
    
end

%merge mono and bino since they are not different
all_fanofs = squeeze(nanmean(all_fanofs(:,:,:,:,:),3)); %computing the mean of between mono and binocular to reduce error size
mvar_vals = squeeze(nanmean(var_vals,2));
mpeak_vals = squeeze(nanmean(peak_vals,2));    


outliers = nan(size(mvar_vals));
for p = 1:size(mvar_vals,1)
    for w =1:size(mvar_vals,2)
        outliers(p,w,:) = ~isoutlier(mvar_vals(p,w,:));
    end
end

select_mvar = mvar_vals.*outliers; %zero out outliers
select_mvar(select_mvar == 0) = NaN; %replace zeros by nans
 
% use variables below in plot code above to replot them
mvar_vals = select_mvar;

mfanofs = mvar_vals./mpeak_vals;

col(1,:) =[146/255 197/255 222/255] ; %--blue 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [166/255 219/255 160/255]; % -- green

col(4,:) =[194/255 165/255 207/255] ; %--purple
col(5,:) = [253/255 174/255 97/255]; % -- orange
%col(3,:) = [166/255 219/255 160/255]; % -- green
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('Greys', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Blues', nlines);


         
figure();
mean_fanof =nanmean(all_fanofs,4);
for p = 1:5
    h =subplot(1,5,p);
    %cmap = flip(cmaps(2).map) ;
    cmap = cmaps(2).map ;
    colormap(cmap);
    if p ==1
        for w =1:length(wsz)
            
            plot(-450+wsz(w)/2:10:-200-wsz(w)/2,mean_fanof(1:len(w),p,w),'LineWidth',2, 'col',cmap(w,:))%'Color',[40/255 40/255 40/255] )
            hold on
        end
        xlim([-375 -275])
    end
    if p ~= 1
        for w =1:length(wsz)
            %ci_low = mean_fanof(:,p,w) - 1.96*std(squeeze(all_fanofs(:,p,w,:)),0,2, 'omitnan')./sqrt(length(find(~isnan(all_fanofs(1,p,w,:)))));
            %ci_high = mean_fanof(:,p,w) + 1.96*std(squeeze(all_fanofs(:,p,w,:)),0,2, 'omitnan')./sqrt(length(find(~isnan(all_fanofs(1,p,w,:)))));
            
            plot(-125+wsz(w)/2:10:125-wsz(w)/2,mean_fanof(1:len(w),p,w),'LineWidth',2, 'col',cmap(w,:))%'Color',[40/255 40/255 40/255] )
            hold on
            %h1= ciplot(ci_low, ci_high,[-125+wsz/2:10:125-wsz/2],col(1,:),0.5);
            %set(h1, 'edgecolor','none')
            
        end
        xlim([-50 50])
    end
    ylim([0.2 1.8])
    
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

 sgtitle('Mean Sliding Window (10 ms steps) Fano Factor', 'FontSize', 20)
 set(gcf,'Units','inches')
 set(gcf,'position',[1 1 15 11])
 plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\',sprintf('all_sliding_window_mean_fanof_rsstimonset_right_plusminus50ms'));
 saveas(gcf,strcat(plotdir, '.png'));
 saveas(gcf,strcat(plotdir, '.svg'));

 %% Fano Factor as a function of peak in a point/box plot with mean and 95% CI
 
%data for plot
halflen = round(len./2);
plot_dat = nan(size(all_fanofs,4),size(all_fanofs,3),length(halflen));
for w = 1:length(halflen)
  % plot_dat(:,:,w) = squeeze(all_fanofs(halflen(w),:,w,:))'; 
   plot_dat(:,:,w) = squeeze(mfanofs(:,w,:))'; %for plotting data without the variance outliers
end
%make data ready for gramm


meanFanofs = reshape(nanmean(plot_dat,1),[size(plot_dat,2)*size(plot_dat,3),1]);
stdFanofs = reshape(std(plot_dat,[],1,'omitnan'),[size(plot_dat,2)*size(plot_dat,3),1]);
ci_high = meanFanofs + 1.96*stdFanofs/sqrt(size(plot_dat,1));
ci_low = meanFanofs - 1.96*stdFanofs/sqrt(size(plot_dat,1));
peakLabel = repmat({'Baseline State'; 'Pk1';'Pk2';'Pk3';'Pk4'},size(plot_dat,3),1);
windowSz = [repmat({'20ms'}, size(plot_dat,2),1);repmat({'30ms'}, size(plot_dat,2),1);repmat({'50ms'}, size(plot_dat,2),1);repmat({'70ms'}, size(plot_dat,2),1);repmat({'x100ms'}, size(plot_dat,2),1)];

%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(1).map) ;
colormap(cmap);

%points and error bars plot
clear g

g(1,1)=gramm('x',categorical(peakLabel),'y',meanFanofs,...
    'ymin',ci_low,'ymax',ci_high,'color',windowSz);
g(1,1).set_names('color','Window Size','y','Fano Factor','x','Peak #');
g(1,1).set_color_options('map',cmap);
g(1,1).geom_point('dodge',0.2);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).set_title('Mean Sliding Window (10 ms steps) Fano Factor');
%g(1,1).axe_property('ylim',[0.2 1.8]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten

figure('Position',[100 100 800 450]);
g.draw();

 plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\',sprintf('all_sliding_window_mean_fanof_rsstimonset_point_ci_excludeVarOutliers'));
 saveas(gcf,strcat(plotdir, '.png'));
 saveas(gcf,strcat(plotdir, '.svg'));
 
%% Perform stats on fano factor distributions: multiple pair-wise comparisons

w =3; %for window size = 50 ms
stat_dat = squeeze(mfanofs(:,w,:))'; %data without the variance outliers %data dimensions dim 1 = units, dim 2 = peaks, dim 3 = window  

%test for normality

%(1) Just take a look at the data if we can se any non-normality in the
%distributions
%prepare data for jitter
peakLabel = [repmat({'Baseline State'}, size(stat_dat,1),1);repmat({'Pk1'}, size(stat_dat,1),1);repmat({'Pk2'}, size(stat_dat,1),1);repmat({'Pk3'}, size(stat_dat,1),1);repmat({'Pk4'}, size(stat_dat,1),1)];
linFfs = reshape(stat_dat, size(stat_dat,1)*size(stat_dat,2), 1); 

%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('Blues', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

%points and error bars plot
clear g

%jitter

g(1,1) = gramm('x',categorical(peakLabel),'y',linFfs);
g(1,1).geom_jitter('width',0,'height',0.2);
g(1,1).set_names('y','Fano Factor','x','Peak #');
g(1,1).axe_property('ylim',[0 2]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g(1,1).set_color_options('map',cmap(4,:));
g(1,1).set_title('Fano Factor distributions across peaks');

figure('Position',[100 100 800 450]);
g.draw();

%(2) Use Shapiro-Wilk test to test for normality on each peak
labs = unique(peakLabel);
alpha = 0.05;
for i = 1:length(unique(peakLabel))
    pkN = labs{i};
    x = linFfs(strcmp(peakLabel, pkN));
[H(i), pValue(i), W(i)] = swtest(x, alpha);
end

%(3) Perform pairwise comparisons corrected for FWER = 0.05
%since distributions aren't all normally distributed, lets use Wilcoxon
%signed rank tests
%here we want to perform 7 two-sided tests, with Holm's method for
%correction for multiple comparisons

%(1) Prestim vs Pk1
x = stat_dat(:,1);
y = stat_dat(:,2);
[p1,~,stats1] = signrank(x,y);
%(2) Prestim vs Pk4
x = stat_dat(:,1);
y = stat_dat(:,5);
[p2,~,stats2] = signrank(x,y);
%(3) Pk1 vs Pk2
x = stat_dat(:,2);
y = stat_dat(:,3);
p3 = signrank(x,y);
%(4) Pk1 vs Pk3
x = stat_dat(:,2);
y = stat_dat(:,4);
p4 = signrank(x,y);
%(5) Pk1 vs Pk4
x = stat_dat(:,2);
y = stat_dat(:,5);
p5 = signrank(x,y);
%(6) Pk2 vs Pk3
x = stat_dat(:,3);
y = stat_dat(:,4);
p6 = signrank(x,y);
%(7) Pk3 vs Pk4
x = stat_dat(:,4);
y = stat_dat(:,5);
[p7,~,stats7] = signrank(x,y);

%% Perform similar analysis on significantly adapting units
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
filename = 'summary_table_pvalues_normmono_meanpks_mono_bino';
T = readtable(strcat(newdatadir, filename, '.csv'));
%T(ismember(T.Wpvalue,'NA'),:)=[]; %remove rows with NAs
Tmono = T(ismember(T.Condition, 'Monocular'),:);
suppressed =  ismember(char(Tmono.Pk1Pk4supress), 'True');
goodrows = (Tmono.Pk1Pk4pvalue <0.05) & suppressed(:,1); %Gert rows of suppressed (significant) single units
Tsuppress = Tmono(goodrows,:);

idx = str2num(char(Tsuppress.Var1)); 

%plot fano factor distribution of suppressed units
w =3; %for window size = 50 ms
halflen = round(len./2); %middle sliding value for each window size
adapt_dat = squeeze(all_fanofs(halflen(w),:,w,idx))'; %fano factor values (we didn't exclude outliers here, we want all adapting neurons FF values%data dimensions dim 1 = units, dim 2 = peaks, dim 3 = window  
%95%CI
meanFfs = nanmean(adapt_dat);
stdFfs = std(adapt_dat,[],'omitnan');
ci_high = meanFfs + 1.96*stdFfs/sqrt(length(adapt_dat));
ci_low = meanFfs - 1.96*stdFfs/sqrt(length(adapt_dat));;
%prepare data for jitter
peakLabel = [repmat({'Baseline State'}, size(adapt_dat,1),1);repmat({'Pk1'}, size(adapt_dat,1),1);repmat({'Pk2'}, size(adapt_dat,1),1);repmat({'Pk3'}, size(adapt_dat,1),1);repmat({'Pk4'}, size(adapt_dat,1),1)];
linFfs = reshape(adapt_dat, size(adapt_dat,1)*size(adapt_dat,2), 1); 

%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('Blues', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

%points and error bars plot
clear g

%jitter

g(1,1) = gramm('x',categorical(peakLabel),'y',linFfs);
g(1,1).geom_jitter('width',0.2,'height',0.2);
g(1,1).set_names('y','Mean Fano Factor','x','Peak #');
g(1,1).axe_property('ylim',[0 1.5]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g(1,1).set_color_options('map',cmap(4,:));
g(1,1).set_title('Mean 50 ms Sliding Window (10 ms steps) Fano Factor');

%point bar plot

g(1,1).update('x',categorical(unique(peakLabel)),'y',meanFfs,...
    'ymin',ci_low,'ymax',ci_high);
g(1,1).set_color_options('map',cmaps(1).map(4,:));
g(1,1).geom_point('dodge',0.2);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).axe_property('ylim',[0 1.5]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten

figure('Position',[100 100 800 450]);
g.draw();
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\',sprintf('adapt_sliding_window50ms_mean_fanofactor_rsstimonset_jitter_point_ci'));
 saveas(gcf,strcat(plotdir, '.png'));
 saveas(gcf,strcat(plotdir, '.svg'));
 
 %perform stats on adapting single units
 %(2) Use Shapiro-Wilk test to test for normality on each peak
labs = unique(peakLabel);
alpha = 0.05;
for i = 1:length(unique(peakLabel))
    pkN = labs{i};
    x = linFfs(strcmp(peakLabel, pkN));
[H(i), pValue(i), W(i)] = swtest(x, alpha);
end

%(3) Perform pairwise comparisons corrected for FWER = 0.05
%since distributions aren't all normally distributed, lets use Wilcoxon
%signed rank tests
%here we want to perform 7 two-sided tests, with Holm's method for
%correction for multiple comparisons

%(1) Prestim vs Pk1
x = adapt_dat(:,1);
y = adapt_dat(:,2);
[p1, ~,stats1] = signrank(x,y);
%(2) Prestim vs Pk4
x = adapt_dat(:,1);
y = adapt_dat(:,5);
p2 = signrank(x,y);
%(3) Pk1 vs Pk2
x = adapt_dat(:,2);
y = adapt_dat(:,3);
p3 = signrank(x,y);
%(4) Pk1 vs Pk3
x = adapt_dat(:,2);
y = adapt_dat(:,4);
p4 = signrank(x,y);
%(5) Pk1 vs Pk4
x = adapt_dat(:,2);
y = adapt_dat(:,5);
p5 = signrank(x,y);
%(6) Pk2 vs Pk3
x = adapt_dat(:,3);
y = adapt_dat(:,4);
p6 = signrank(x,y);
%(7) Pk3 vs Pk4
x = adapt_dat(:,4);
y = adapt_dat(:,5);
[p7,~,stats7] = signrank(x,y);

%% Test for a non linear trend
%modelfun = @(b,x)(b(1)+b(2)*exp(b(3)*x));

%modelfun = @(b,x)(b(1)+b(2)*x.^b(3));
b = [0.1,0.3,0.2, 0.2];
%modelfun = @(b,x)(b(1)+b(2)*x+b(3)*x.^b(4)); %polynomial
%modelfun = @(b,x)(b(1)+b(2)*x+b(3)*x.^2); %polynomial with only 3 coefficients to estimate
%modelfun = @(b,x)(b(1)+b(2)* (x-(b(3))).^2);
%b = [.2,.2];
%modelfun = @(b,x)(b(1)+b(2)*log(x));

adapt_dat2 = adapt_dat(:,2:end);
adapt_dat3 = adapt_dat2(isfinite(adapt_dat2)); 
%fit model
rng('default') % for reproducibility

x = [1*ones(9,1); 2*ones(9,1); 3*ones(9,1); 4*ones(9,1)];
y = adapt_dat3;

%use nlinfit
%{
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
beta0 = [2;2;2];
[beta,R,J,CovB,MSE] = nlinfit(x,y,modelfun,beta0,opts);

%get 95%CI
xrange = min(x):.01:max(x);
[ypred,delta] = nlpredci(modelfun,xrange,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower = ypred - delta;
upper = ypred + delta;
%}
%use more recent fitnlm function
tabxy = table(x,y);

%modelfun = @(b,x)(b(1)+b(2)* (x-(b(3))).^2); %parabolic
%modelfun = @(b,x)(b(1)+b(2)* (x)); %linear
modelfun = @(b,x)(b(1)+b(2)* (x).^(b(3)));%logarithm

beta0 = [2;2;2];
mdl = fitnlm(tabxy,modelfun,beta0);
AIC = mdl.ModelCriterion.AIC;

%get 95%CI
xrange = min(x):.01:max(x);
[ypred,delta] = nlpredci(modelfun,xrange,mdl.Coefficients.Estimate,mdl.Residuals.Raw,'Covar',mdl.CoefficientCovariance,...
                         'MSE',mdl.MSE,'SimOpt','on');
lower = ypred - delta;
upper = ypred + delta;

figure()
plot(x,y,'ko') % observed data
hold on
plot(xrange,ypred,'k','LineWidth',2)
plot(xrange,[lower;upper],'r--','LineWidth',1.5)
hold on
plot([1 2 3 4], meanFfs(2:5),'b-o')

set(gca,'box','off')
set(gca, 'linewidth', 2)
title(sprintf('x = %.2d + %.2d*x^(%.2d)', mdl.Coefficients.Estimate),'Interpreter', 'none');

xlim([0 5])
ylim([0 1])
xlabel('Peak #')
ylabel('Fano Factor')
text(2.5,.5,sprintf('AIC = %d', AIC));
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\',sprintf('nlfit_logarithm'));
 saveas(gcf,strcat(plotdir, '.png'));
 saveas(gcf,strcat(plotdir, '.svg'));

%% Plot population linear regressions of mean versus trial-to-trial variance comparing peaks with resting state
%window size = 50ms
peakLabel = repmat({'Baseline State';'Pk1';'Pk2';'Pk3';'Pk4'}, size(mvar_vals,3),1);
linVars = reshape(squeeze(mvar_vals(:,3,:)), size(mvar_vals,1)*size(mvar_vals,3), 1); 
linMeans = reshape(mpeak_vals(:,3,:), size(mpeak_vals,1)*size(mpeak_vals,3), 1);

nlines = 7;
cmap =flip(cbrewer2('Blues', nlines));
colormap(cmap);


figure('Position',[100 100 800 400],'Color',[1 1 1]);
% Define groups
peakLab = unique(peakLabel); % Based on data

% Loop over groups
for p = 2:length(peakLab) % External loop on the axes
    
    % Axes creation
    ax = subplot(1,length(peakLab)-1,p-1);
    hold on
    % Pre-stimulation Data selection
    sel = strcmp(peakLabel,peakLab{1}) &...
        ~isnan(linVars);
    x1 = linMeans(sel);
    y1 = linVars(sel);
    % Plotting of raw data
    %linear regression
    coeffs1 = polyfit(x1(isfinite(x1) & isfinite(y1)),y1(isfinite(x1) & isfinite(y1)),1);
    f1 = polyval(coeffs1,x1);
    plot(x1, y1,'o',x1, f1,'-','Color',[160/255 160/255 160/255],'MarkerSize',2, 'MarkerFaceColor',[160/255 160/255 160/255],'linewidth',2)
    xlim([0 10])
    ylim([0 4.5])
    % Peak Data selection
    sel = strcmp(peakLabel,peakLab{p}) &...
        ~isnan(linVars);
    x = linMeans(sel);
    y = linVars(sel);
    % Plotting of raw data
    % Keep the same color for the statistics
    coeffs = polyfit(x(isfinite(x) & isfinite(y)),y(isfinite(x) & isfinite(y)),1);
    f = polyval(coeffs,x);
    plot(x, y,'o',x, f,'-','Color',cmap(4,:),'MarkerSize',2, 'MarkerFaceColor',cmap(4,:),'linewidth',2)
    xlim([0 10])
    ylim([0 4.5])
    hold on
    %format short
    %text(max(x)/12,max(y)/2, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)));
    %ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    set(gca, 'linewidth',2)
    set(gca,'box','off')
    % Statistics (linear fit and plotting)
    %{
        b = [ones(sum(sel),1) linMeans(sel)] \ ...
			linVars(sel);
        x_fit = [min(linMeans(sel)) ...
			max(linVars(sel))];
        plot(x_fit, x_fit * b(2) + b(1),'LineWidth',1.5);
    %}
    % Axes legends
    title(['Label: ' peakLab{p}]);
    xlabel('Mean spike count');
    ylabel('Spike count variance');
end

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\',sprintf('scatter_linreg_50ms_window_mean_vs_variance_rsstimonset_excludeVarOutliers'));
 saveas(gcf,strcat(plotdir, '.png'));
 saveas(gcf,strcat(plotdir, '.svg'));

 %% Test significance of slope relative to 0
 peakLabel = repmat({'Baseline State';'Pk1';'Pk2';'Pk3';'Pk4'}, size(mvar_vals,3),1);
linVars = reshape(squeeze(mvar_vals(:,3,:)), size(mvar_vals,1)*size(mvar_vals,3), 1); 
linMeans = reshape(mpeak_vals(:,3,:), size(mpeak_vals,1)*size(mpeak_vals,3), 1);
peakLab = unique(peakLabel); % Based on data

linPvalue = zeros(5,1);
Rsquared = zeros(5,1);
slope = zeros(5,1);
% Loop over groups
for p = 2:length(peakLab) % External loop on the axes

    % Pre-stimulation Data selection
    sel = strcmp(peakLabel,peakLab{1}) &...
        ~isnan(linVars);
    x1 = linMeans(sel);
    y1 = linVars(sel);
    linreg = fitlm(x1,y1);
    [linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linreg);
    Rsquared(1,1) = linreg.Rsquared.Ordinary;
    coeffs1 = polyfit(x1(isfinite(x1) & isfinite(y1)),y1(isfinite(x1) & isfinite(y1)),1);
    slope(1,1) = coeffs1(1);
    % Plotting of raw data
    %linear regression
    
    %f1 = polyval(coeffs1,x1);
    %plot(x1, y1,'o',x1, f1,'-','Color',[160/255 160/255 160/255],'MarkerSize',2, 'MarkerFaceColor',[160/255 160/255 160/255],'linewidth',2)

    % Peak Data selection
    sel = strcmp(peakLabel,peakLab{p}) &...
        ~isnan(linVars);
    x = linMeans(sel);
    y = linVars(sel);
    % Plotting of raw data
    linreg2 = fitlm(x,y);
    [linPvalue(p,1),F(p,1), r(p,1)] = coefTest(linreg2);
    Rsquared(p,1) = linreg2.Rsquared.Ordinary;
    coeffs = polyfit(x(isfinite(x) & isfinite(y)),y(isfinite(x) & isfinite(y)),1);
    slope(p,1) = coeffs(1);
   % f = polyval(coeffs,x);
   % plot(x, y,'o',x, f,'-','Color',cmap(4,:),'MarkerSize',2, 'MarkerFaceColor',cmap(4,:),'linewidth',2)
   
end

 
%% Plot linear regressions of mean versus trial-to-trial variance comparing peaks with resting state 
%for ADAPTING neurons
%window size = 50ms
adapt_vars = squeeze(mvar_vals(:,w,idx))'; 
adapt_mpeaks = squeeze(mpeak_vals(:,w,idx))'; 
peakLabel = repmat({'Baseline State';'Pk1';'Pk2';'Pk3';'Pk4'}, length(idx),1);
linVars = reshape(adapt_vars, size(adapt_vars,1)*size(adapt_vars,2), 1); 
linMeans = reshape(adapt_mpeaks, size(adapt_mpeaks,1)*size(adapt_mpeaks,2), 1);

nlines = 7;
cmap =flip(cbrewer2('Blues', nlines));
colormap(cmap);


figure('Position',[100 100 800 400],'Color',[1 1 1]);
% Define groups
peakLab = unique(peakLabel); % Based on data

% Loop over groups
for p = 2:length(peakLab) % External loop on the axes
    
    % Axes creation
    ax = subplot(1,length(peakLab)-1,p-1);
    hold on
    % Pre-stimulation Data selection
    sel = strcmp(peakLabel,peakLab{1}) &...
        ~isnan(linVars);
    x1 = linMeans(sel);
    y1 = linVars(sel);
    % Plotting of raw data
    %linear regression
    coeffs1 = polyfit(x1(isfinite(x1) & isfinite(y1)),y1(isfinite(x1) & isfinite(y1)),1);
    f1 = polyval(coeffs1,x1);
    plot(x1, y1,'o',x1, f1,'-','Color',[160/255 160/255 160/255],'MarkerSize',2, 'MarkerFaceColor',[160/255 160/255 160/255],'linewidth',2)
    xlim([0 10])
    ylim([0 6])
    % Peak Data selection
    sel = strcmp(peakLabel,peakLab{p}) &...
        ~isnan(linVars);
    x = linMeans(sel);
    y = linVars(sel);
    % Plotting of raw data
    % Keep the same color for the statistics
    coeffs = polyfit(x(isfinite(x) & isfinite(y)),y(isfinite(x) & isfinite(y)),1);
    f = polyval(coeffs,x);
    plot(x, y,'o',x, f,'-','Color',cmap(4,:),'MarkerSize',2, 'MarkerFaceColor',cmap(4,:),'linewidth',2)
    xlim([0 10])
    ylim([0 6])
    hold on
    %format short
    %text(max(x)/12,max(y)/2, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)));
    %ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    set(gca, 'linewidth',2)
    set(gca,'box','off')
    % Statistics (linear fit and plotting)
    %{
        b = [ones(sum(sel),1) linMeans(sel)] \ ...
			linVars(sel);
        x_fit = [min(linMeans(sel)) ...
			max(linVars(sel))];
        plot(x_fit, x_fit * b(2) + b(1),'LineWidth',1.5);
    %}
    % Axes legends
    title(['Label: ' peakLab{p}]);
    xlabel('Mean spike count');
    ylabel('Spike count variance');
end

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\',sprintf('scatter_linreg_50ms_window_mean_vs_variance_rsstimonset_excludeVarOutliers'));
 saveas(gcf,strcat(plotdir, '.png'));
 saveas(gcf,strcat(plotdir, '.svg'));


%% Variance (population) point plot 
%{
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
%all_fanofs = nan(length(-125:10:125-wsz(1)), 4,2,length(wsz),length(filenames));
peak_vals = nan(4,2,length(wsz),length(filenames));
var_vals = nan(5,2,length(wsz),length(filenames));
%peak_fanofs = nan(4,2,length(wsz),length(filenames));

for i =1:length(filenames)
    filename = filenames{i};
    if length(fieldnames(slide_win_fanof.(filename).fanof)) == 2
        for b =1:2
            binN = sprintf('bin%d',bins(b));
            for p =1:4
                for w =1:length(wsz)
                    windSz = sprintf('wsz%d',wsz(w));
                    len(w) = length(slide_win_fanof.(filename).fanof.(binN).(windSz).peaks(p,:));
                    var_vals(p+1,b,w,i) = slide_win_fanof.(filename).varspkc.(binN).(windSz).peaks(p,round(len(w)/2));
                    peak_vals(p+1,b,w,i) = slide_win_fanof.(filename).meanspkc.(binN).(windSz).peaks(p,round(len(w)/2));
                  
                end
            end
            %fill up first  matrix column with FF of resting state 
              for w =1:length(wsz)
                    windSz = sprintf('wsz%d',wsz(w));
                    len(w) = length(slide_win_fanof.(filename).varspkc.(binN).(windSz).rs(:));
                    var_vals(1,b,w,i) = slide_win_fanof.(filename).varspkc.(binN).(windSz).rs(round(len(w)/2));
                    peak_vals(1,b,w,i) = slide_win_fanof.(filename).meanspkc.(binN).(windSz).rs(round(len(w)/2));
                  
              end
        end
    end
end

%merge mono and bino since they are not different
mvar_vals = squeeze(nanmean(var_vals,2));
mpeak_vals = squeeze(nanmean(peak_vals,2));    
%}
%prepare data for gramm
%prepare data for bar/point
meanVars = reshape(nanmean(mvar_vals,3),[size(mvar_vals,1)*size(mvar_vals,2),1]);
stdVars = reshape(std(mvar_vals,[],3,'omitnan'),[size(mvar_vals,1)*size(mvar_vals,2),1]);
ci_high = meanVars + 1.96*stdVars/sqrt(size(mvar_vals,3));
ci_low = meanVars - 1.96*stdVars/sqrt(size(mvar_vals,3));
peakLabel = repmat({'Baseline State'; 'Pk1';'Pk2';'Pk3';'Pk4'},size(mvar_vals,2),1);
windowSz = [repmat({'20ms'}, size(mvar_vals,1),1);repmat({'30ms'}, size(mvar_vals,1),1);repmat({'50ms'}, size(mvar_vals,1),1);repmat({'70ms'}, size(mvar_vals,1),1);repmat({'x100ms'}, size(mvar_vals,1),1)];

%prepare data for jitter
longPeakLabel = repmat(repmat({'Baseline State';'Pk1';'Pk2';'Pk3';'Pk4'}, size(mvar_vals,2),1),size(mvar_vals,3),1);
longWindowSz = repmat([repmat({'20ms'}, size(mvar_vals,1),1);repmat({'30ms'}, size(mvar_vals,1),1);repmat({'50ms'}, size(mvar_vals,1),1);repmat({'70ms'}, size(mvar_vals,1),1);repmat({'x100ms'}, size(mvar_vals,1),1)], size(mvar_vals,3),1);
linVars = reshape(reshape(mvar_vals, [size(mvar_vals,1)*size(mvar_vals,2), size(mvar_vals,3)]), [size(mvar_vals,1)*size(mvar_vals,2)*size(mvar_vals,3),1]); 


%colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('Blues', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

%points and error bars plot
clear g

%jitter

g(1,1) = gramm('x',categorical(longPeakLabel),'y',linVars, 'subset',strcmp(longWindowSz, '50ms') );
g(1,1).geom_jitter('width',0,'height',0.2);
g(1,1).set_names('color','Window Size','y','Mean Variance','x','Peak #');
g(1,1).axe_property('ylim',[0 7]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
g(1,1).set_color_options('map',cmap(4,:));
g(1,1).set_title('Mean 50 ms Sliding Window (10 ms steps) Variance');

%point bar plot
g(1,1).update('x',categorical(peakLabel),'y',meanVars,...
    'ymin',ci_low,'ymax',ci_high,'subset',strcmp(windowSz, '50ms'));
g(1,1).set_color_options('map',cmaps(1).map(4,:));
g(1,1).geom_point('dodge',0.2);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).axe_property('ylim',[0 7]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten


figure('Position',[100 100 800 450]);
g.draw();
plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\',sprintf('all_sliding_window_mean_variance_rsstimonset_jitter_point_ci_no_outlier'));
 saveas(gcf,strcat(plotdir, '.png'));
 saveas(gcf,strcat(plotdir, '.svg'));


%% Monocular vs Binocular condition comparison
%% Plot linear regressions of mean versus trial-to-trial variance comparing peaks in the monocular vs binocular condition
%merge mono and bino since they are not different

%all_fanofs = squeeze(nanmean(all_fanofs(:,:,:,:,:),3)); %computing the mean of between mono and binocular to reduce error size
%mvar_vals = squeeze(nanmean(var_vals,2));
%mpeak_vals = squeeze(nanmean(peak_vals,2));

outliers = nan(size(var_vals));
for b =1:2
    for p = 1:size(var_vals,1)
        for w =1:size(var_vals,3)
            outliers(p,b,w,:) = ~isoutlier(var_vals(p,b,w,:));
        end
    end
end
select_var = var_vals.*outliers; %zero out outliers
select_var(select_var == 0) = NaN; %replace zeros by nans
% use variables below in plot code above to replot them
svar_vals = select_var;
sfanofs = svar_vals./peak_vals;
%window size = 50ms

peakLabel = repmat({'Baseline State';'Pk1';'Pk2';'Pk3';'Pk4'}, size(svar_vals,4),2);
linVars = reshape(squeeze(svar_vals(:,:,3,:)), size(svar_vals,1)*size(svar_vals,4), 2); 
linMeans = reshape(peak_vals(:,:,3,:), size(peak_vals,1)*size(peak_vals,4), 2);

nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);

figure('Position',[100 100 800 400],'Color',[1 1 1]);
% Define groups
peakLab = unique(peakLabel); % Based on data


    % Loop over groups
    for p = 2:length(peakLab) % External loop on the axes
        % Axes creation
        ax = subplot(1,length(peakLab)-1,p-1);
        hold on
        % Binocular Data selection
        sel = strcmp(peakLabel(:,2),peakLab{p}) &...
            ~isnan(linVars(:,2));
        x1 = linMeans(sel,2);
        y1 = linVars(sel,2);
        % Plotting of Bino data
        %linear regression
        coeffs1 = polyfit(x1(isfinite(x1) & isfinite(y1)),y1(isfinite(x1) & isfinite(y1)),1);
        f1 = polyval(coeffs1,x1);
        plot(x1, y1,'o',x1, f1,'-','Color',cmap(3,:),'MarkerSize',2, 'MarkerFaceColor',cmap(3,:),'linewidth',2)
        xlim([0 10])
        ylim([0 4.5])
        % Monocular Peak Data selection
        sel = strcmp(peakLabel(:,1),peakLab{p}) &...
            ~isnan(linVars(:,1));
        x = linMeans(sel,1);
        y = linVars(sel,1);
        % Plotting of raw data
        % Keep the same color for the statistics
        coeffs = polyfit(x(isfinite(x) & isfinite(y)),y(isfinite(x) & isfinite(y)),1);
        f = polyval(coeffs,x);
        plot(x, y,'o',x, f,'-','Color',cmaps(1).map(4,:),'MarkerSize',2, 'MarkerFaceColor',cmaps(1).map(4,:),'linewidth',2)
        xlim([0 10])
        ylim([0 4.5])
        hold on
        %format short
        %text(max(x)/12,max(y)/2, sprintf('y = %.2f + %.2f*x', round(coeffs(2),2), round(coeffs(1),2)));
        %ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
        set(gca, 'linewidth',2)
        set(gca,'box','off')
        % Statistics (linear fit and plotting)
        %{
        b = [ones(sum(sel),1) linMeans(sel)] \ ...
			linVars(sel);
        x_fit = [min(linMeans(sel)) ...
			max(linVars(sel))];
        plot(x_fit, x_fit * b(2) + b(1),'LineWidth',1.5);
        %}
        % Axes legends
        title(['Label: ' peakLab{p}]);
        xlabel('Mean spike count');
        ylabel('Spike count variance');
    end

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\',sprintf('scatter_linreg_50ms_window_mean_vs_variance_mono_bino_excludeVarOutliers'));
 saveas(gcf,strcat(plotdir, '.png'));
 saveas(gcf,strcat(plotdir, '.svg'));



%% %%%%%%%% raster plot of example single unit %%%%
 %raster plot of trials binary spikesaligned to each peak
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
                 plot(x-125,y,'Color','k')
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
    % saveas(gcf,strcat(plotdir, '.png'));
    % saveas(gcf,strcat(plotdir, '.svg'));
     
     
     
 end
 
%% raster plots with trials triggered to stimulus onset time
datadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
binSpkTrials = load(strcat(datadir,'binary_trials_data_06022021'));
binSpkTrials = binSpkTrials.binSpkTrials;
filenames = fieldnames(binSpkTrials);
 
for i =40:length(filenames)
    filename = filenames{i};
    figure();
    
    if isfield(binSpkTrials.(filename),'bin1')
        for tr =1:length(binSpkTrials.(filename).bin1(1,:))
            spikeTimes =find(binSpkTrials.(filename).bin1(:,tr)==1)';
            x = repmat(spikeTimes,3,1);
            y = nan(size(x));
            
            if ~isempty(y)
                y(1,:) = tr-1;
                y(2,:) = tr;
            end
            plot(x-600,y,'Color','k')
            hold on
            
        end
        set(gca, 'linewidth',2)
        set(gca,'box','off')
        xlabel('Time (ms)')
        ylabel('Trial number')
        
    end

     
     sgtitle({'Binary spike data of an example single unit','aligned to each peak'})
     set(gcf,'position',get(gcf,'position').*[1 1 1.15 1])
     plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\noise_suppression\plots\', sprintf('stimtrigg_raster_plots_single_unit_example%d',i));
     saveas(gcf,strcat(plotdir, '.png'));
     saveas(gcf,strcat(plotdir, '.svg'));
     
     
     
 end

%%

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
                    len(w) = length(slide_win_fanof.(filename).fanof.(binN).(windSz).peaks(p,:));
                    all_fanofs(1:len(w),p,b,w,i) =  slide_win_fanof.(filename).fanof.(binN).(windSz).peaks(p,:);
                    peak_fanofs(p,b,w,i) =  slide_win_fanof.(filename).fanof.(binN).(windSz).peaks(p,round(len(w)/2));
                    peak_vals(p,b,w,i) = slide_win_fanof.(filename).meanspkc.(binN).(windSz).peaks(p,round(len(w)/2));
                    var_vals(p,b,w,i) = slide_win_fanof.(filename).varspkc.(binN).(windSz).peaks(p,round(len(w)/2));
                    
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
    %xlim([0 4])
    %ylim([0 5])
    if p ==1 
    xlabel('Mean Spike Count in a 50ms Window')
    ylabel('Fano Factor')
    end
    title(sprintf('Peak %d', p))
    if p >1
        ax1 = gca;
     %  ax1.YAxis.Visible = 'off';
       
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
    %xlim([0 4])
    %ylim([0 5])
    if p ==1 
    xlabel('Spike Count Variance in a 50ms Window')
    ylabel('Fano Factor')
    end
    title(sprintf('Peak %d', p))
    if p >1
        ax1 = gca;
     %  ax1.YAxis.Visible = 'off';
       
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
    xlim([0 10])
    ylim([0 10])
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

