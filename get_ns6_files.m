
        
   %%
   
%metadata = load('C:\Users\daumail\Documents\LGN_data\single_units\kilosorted_files\160602_I_p01_ss.mat');
sessionslist = load('C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\selected_units_sessions.mat');
charsessions =cell(length(sessionslist.selectedfilenames),1);
for i = 1:length(sessionslist.selectedfilenames)
    if ~isempty(sessionslist.selectedfilenames{i})
    charsessions(i) =sessionslist.selectedfilenames(i);
    else
     charsessions(i) = cellstr(string);       
    end
end

Unique_sessions =uniqueStrCell(charsessions);
 for i = 1:length(Unique_sessions)-1 
     %Big Drobo 2 - Backup\Drobo2\data\rig021
     %r1a\maierlab\DATA\NEUROPHYS\rig021\
     myFolder=strcat('Y:\r1a\maierlab\DATA\NEUROPHYS\rig021\', Unique_sessions{i+1});
     if exist(myFolder, 'dir')
mkdir(strcat('C:\Users\daumail\Documents\LGN_data\single_units\ns2_selected_units\', Unique_sessions{i+1}))
sourcedir=strcat('Y:\r1a\maierlab\DATA\NEUROPHYS\rig021\', Unique_sessions{i+1});
destdir=strcat('C:\Users\daumail\Documents\LGN_data\single_units\ns2_selected_units\',Unique_sessions{i+1});
     
status = copyfile( strcat(sourcedir,'\*.ns2'),  destdir);

     end
 end
%'/Volumes/Drobo2/data/neurophys/KiloSort-ed/LGNgrouped/' dates{d} 

%% Kacie'S Code
%{
bigdrobopath           = '/volumes/BigDrobo2/Drobo2/data/';
dates                  = unique(DataList.Datestr(:));
missing                = []; 
specific_date          = find(ismember(dates,{'190326'}));

for d =  specific_date' 


    fprintf('d : %u\n',d); 

    clear these_idfilelist
    these_id               =  find(ismember(DataList.Datestr,dates{d}));
    for p =[these_id]'

        clearvars -except DataList these_id p d targetpath flag_checkforexisting dates monocdates  bigdrobopath missing extHDpath

        %try
         depth = sprintf('_p%02u',sum(ismember(DataList.Datestr(1:p-1),DataList.Datestr(p))) + 1);
        if DataList.Drobo(p) == 1
            drname = sprintf('/volumes/drobo/data/neurophys/rig%03u/%s_%s/',DataList.Rig(p),DataList.Datestr{p},DataList.Monkey{p});
        elseif DataList.Drobo(p) == 2 & (datenum(dates{d},'yymmdd') < datenum('05/05/2018','mm/dd/yyyy'))
            drname = sprintf('/volumes/drobo2/data/neurophys/rig%03u/%s_%s/',DataList.Rig(p),DataList.Datestr{p},DataList.Monkey{p});

        else
            drname = sprintf('%s/rig%03u/%s_%s/',bigdrobopath,DataList.Rig(p),DataList.Datestr{p},DataList.Monkey{p});
        end

%}