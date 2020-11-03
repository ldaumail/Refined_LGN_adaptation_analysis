
%ns6 directory 
ns6dir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\ns6_selected_units\';

%ssdir
ssdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\kilosorted_files\';
%ssdir content
selected_name_list = dir(ssdir);  

%load file with selected trial indexes for a given penetration with
%penetration name
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'selected_trials_idx']);

%load selected penetrations files list
%selected_penetrations = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\selected_orig_units_filenames';
%penetrations_names = textscan( fopen(strcat(selected_penetrations, '.txt')), '%s', 'HeaderLines', 1);
%penetrations_names = cell2struct(penetrations_names, 'penetrations_names');
SingUnitns6Files = cell(71,3);
for pnt = 1:length(selected_trials_idx.logicals)
    if isnumeric(SingUnitns6Files{pnt,1})
    SingUnitns6Files{pnt,1} = num2str(SingUnitns6Files{pnt,1});
    SingUnitns6Files{pnt,2} = num2str(SingUnitns6Files{pnt,2});
        if ~isempty(selected_trials_idx.logicals(pnt).idx)
        penetration_name = erase(selected_trials_idx.logicals(pnt).penetration, 'matDE50_NDE0');
        STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);

        % get appropriate ns6 file of penetration 
        %isolate session date/animal for ns6 folder name
        underscore = strfind(penetration_name, '_');
        session =  penetration_name(1:underscore(2)-1);
        %load ns6 file of penetration of interest  
        %loop through ss files until we find an ss file with same
        %session and penetration as the penetration of interest to get the
        %right file name for the .NS6 files of interest
            for ss =1:length(selected_name_list)
                if contains(erase(selected_name_list(ss).name,'_ss.mat'),penetration_name(1:underscore(3)-1))
                    ss_file = load(strcat(ssdir,selected_name_list(ss).name));
                    ss_file_fieldnames = fieldnames(ss_file.ss);

                    cinterocdrft_names = ss_file_fieldnames(contains(ss_file_fieldnames,'cinterocdrft'));

                    for cint =1:length(cinterocdrft_names)
                        ns6_filename = char(cinterocdrft_names(cint));
                        if  ~contains(strjoin(SingUnitns6Files(1:pnt,1)), ns6_filename)
                            SingUnitns6Files{pnt,1} = ns6_filename(2:end);
                            SingUnitns6Files{pnt,2} = penetration_name;

                        end
                    end

                end
            end
                SingUnitns6Files{pnt,3} = STIM_file.STIM.chan;

                else
                    if isempty(selected_trials_idx.logicals(pnt).idx)
                    SingUnitns6Files{pnt,1} = char();
                    end
        end
    end
end
metadata = struct();
metadata.ns6FileName = SingUnitns6Files(:,1);
metadata.STIMFileName = SingUnitns6Files(:,2);
metadata.channel = SingUnitns6Files(:,3);




save(strcat(indexdir, 'single_units_ns6_metadata.mat'),'-struct', 'metadata'),  



