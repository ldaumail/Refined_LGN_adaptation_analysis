%This script was developped to plot a distribution of the number of trials
%analyzed for each unit
%Loic Daumail 05/27/2022


%Using the updated functions to isolate peaks
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;


unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);


cellClass = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
cellClass([1,46,55]) = [];
%{
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
%}
%select trials with at least 4 peak values, of convenient quality. Keep peak locations and trial responses: 
[peakLocs, NoFiltMultiContSUA] = peakLocsTrialSelection(unitsData, filenames);

%Store Peaks and peak-triggered trials
[peak_vals, peak_aligned_trials] = peaksAndPeakTrigResps(peakLocs, NoFiltMultiContSUA);



%% Plot number of trials per unit
%get number of trials per unit in mono and binocular condition
bins = [1,6];
tn = [];
filenames = fieldnames(peak_aligned_trials);
for i = 1:length(filenames)
   if length(fieldnames(peak_aligned_trials.(filenames{i}).origin)) == 2
    for b =1:2
    bin = sprintf('bin%d',bins(b));
    tn(i,b) = nnz(~isnan(peak_aligned_trials.(filenames{i}).origin.(bin).pk1(1,:)));
    end
   end
end
tn = tn(tn(:,1) ~= 0,:);

%set up colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);
%set up data for gramm
trial_numbers = [tn(:,1);tn(:,2)];
condition = [repmat({'Monocular'},length(trial_numbers)/2,1); repmat({'Binocular'},length(trial_numbers)/2,1)];

%Plot
%1) Plot the scatter
clear g
g(1,1)=gramm('x',condition,'y',trial_numbers,'color',condition);

g(1,1).geom_jitter('width',0.4,'height',0); 
g(1,1).set_names('x','','y', '');
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]);
g(1,1).no_legend();

%Plot the mean and std for trial numbers
alltnMono = trial_numbers(strcmp(condition, 'Monocular'));
alltnBino = trial_numbers(strcmp(condition, 'Binocular'));
meanTn = nanmean([alltnMono,alltnBino],1);
shortCondition = unique(condition);
stdTns = std([alltnMono,alltnBino], 'omitnan');%/sqrt(length([alltnMono,alltnBino]));

ci_low = meanTn - stdTns;
ci_high = meanTn + stdTns;
g(1,1).update('x', shortCondition, 'y',meanTn,...
    'ymin',ci_low,'ymax',ci_high,'color',shortCondition);
%g(1,1).set_color_options('map','r');
g(1,1).geom_point('dodge',0.5);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_names('x','','y', 'Number of Trials');
g(1,1).no_legend();
g(1,1).set_point_options('base_size',7);

g(1,2) = gramm('x',ones(1,length(tn(:,1))), 'y', tn(:,1)-tn(:,2) );
g(1,2).geom_jitter('width',0.3,'height',0); 
g(1,2).axe_property('xlim', [0 2], 'ylim',[-30 30], 'XTickLabel','','XTick',''); 
g(1,2).set_names('x','','y', 'Difference Monocular-Binocular (Number of Trials)');
g.set_title('Number of trials in selected units post-processing');
g(1,2).set_color_options('map',[0.5 0.5 0.5]);


%Plot mean difference and 95%CI
meanDiff = nanmean([alltnMono-alltnBino],1);
stdDiffs = std([alltnMono-alltnBino], 'omitnan');%/sqrt(length([alltnMono-alltnBino]));
ci_low = meanDiff - stdDiffs;
ci_high = meanDiff + stdDiffs;
g(1,2).update('x', ones(1,1), 'y',meanDiff,...
    'ymin',ci_low,'ymax',ci_high);
g(1,2).geom_point('dodge',0.5);
g(1,2).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,2).axe_property('xlim', [0 2], 'ylim',[-30 30], 'XTickLabel','','XTick',''); 
%g(1,2).set_names('x','','y', '');
g(1,2).no_legend();
g(1,2).set_point_options('base_size',7);

figure('Position',[100 100 800 550]);

g.draw();

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_boxpointbar_std_numtrials_per_unit');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

%% Use previously isolated peaks

newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\';
channelfilename = [newdatadir 'all_orig_bs_zscore_trials_05022021_mono_bino']; 
peak_aligned_trials = load(channelfilename);


filenames = fieldnames(peak_aligned_trials.peak_aligned_trials);
bins = [1,6];


tn = [];
filenames = fieldnames(peak_aligned_trials.peak_aligned_trials);
for i = 1:length(filenames)
   if length(fieldnames(peak_aligned_trials.peak_aligned_trials.(filenames{i}).origin)) == 2
    for b =1:2
    bin = sprintf('bin%d',bins(b));
    tn(i,b) = nnz(~isnan(peak_aligned_trials.peak_aligned_trials.(filenames{i}).origin.(bin).pk1(1,:)));
    end
   end
end
tn = tn(tn(:,1) ~= 0,:);

%set up colors
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('BuPu', nlines);
cmaps(3).map =cbrewer2('Greens', nlines);
cmap = flip(cmaps(2).map) ;
colormap(cmap);
%set up data for gramm
trial_numbers = [tn(:,1);tn(:,2)];
condition = [repmat({'Monocular'},length(trial_numbers)/2,1); repmat({'Binocular'},length(trial_numbers)/2,1)];

%Plot
%1) Plot the scatter
clear g
g(1,1)=gramm('x',condition,'y',trial_numbers,'color',condition);

g(1,1).geom_jitter('width',0.4,'height',0); 
g(1,1).set_names('x','','y', '');
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_color_options('map',[cmap(3,:);cmaps(1).map(4,:)]);
g(1,1).no_legend();

%Plot the mean and std for trial numbers
alltnMono = trial_numbers(strcmp(condition, 'Monocular'));
alltnBino = trial_numbers(strcmp(condition, 'Binocular'));
meanTn = nanmean([alltnMono,alltnBino],1);
shortCondition = unique(condition);
stdTns = std([alltnMono,alltnBino], 'omitnan');%/sqrt(length([alltnMono,alltnBino]));

ci_low = meanTn - stdTns;
ci_high = meanTn + stdTns;
g(1,1).update('x', shortCondition, 'y',meanTn,...
    'ymin',ci_low,'ymax',ci_high,'color',shortCondition);
%g(1,1).set_color_options('map','r');
g(1,1).geom_point('dodge',0.5);
g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,1).axe_property('xlim', [0 4], 'ylim',[0 50],'XTickLabel','','XTick',''); 
g(1,1).set_names('x','','y', 'Number of Trials');
g(1,1).no_legend();
g(1,1).set_point_options('base_size',7);

g(1,2) = gramm('x',ones(1,length(tn(:,1))), 'y', tn(:,1)-tn(:,2) );
g(1,2).geom_jitter('width',0.3,'height',0); 
g(1,2).axe_property('xlim', [0 2], 'ylim',[-30 30], 'XTickLabel','','XTick',''); 
g(1,2).set_names('x','','y', 'Difference Monocular-Binocular (Number of Trials)');
g.set_title('Number of trials in selected units post-processing');
g(1,2).set_color_options('map',[0.5 0.5 0.5]);


%Plot mean difference and 95%CI
meanDiff = nanmean([alltnMono-alltnBino],1);
stdDiffs = std([alltnMono-alltnBino], 'omitnan');%/sqrt(length([alltnMono-alltnBino]));
ci_low = meanDiff - stdDiffs;
ci_high = meanDiff + stdDiffs;
g(1,2).update('x', ones(1,1), 'y',meanDiff,...
    'ymin',ci_low,'ymax',ci_high);
g(1,2).geom_point('dodge',0.5);
g(1,2).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
g(1,2).axe_property('xlim', [0 2], 'ylim',[-20 30], 'XTickLabel','','XTick',''); 
%g(1,2).set_names('x','','y', '');
g(1,2).no_legend();
g(1,2).set_point_options('base_size',7);

figure('Position',[100 100 800 550]);

g.draw();

plotdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\plots\jitter_boxpointbar_std_numtrials_per_unit_old');
saveas(gcf,strcat(plotdir, '.png'));
saveas(gcf,strcat(plotdir, '.svg'));

