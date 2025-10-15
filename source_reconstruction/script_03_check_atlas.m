%% source reconstruction fieldtrip
% https://www.youtube.com/watch?v=Ez72OFjSABs&ab_channel=fieldtriptoolbox

clear variables; close all;

path_utilities = '/path_to_local/scripts/source_reconstruction/utilities';
addpath(path_utilities);
addpath('/path_to_software/fieldtrip-lite-2020083');
path_mri = '/path_to_local/data/MRI';
path_head = '/path_to_local/data/head_shape';
path_meg = '/path_to_local/data/MEG';
path_preproc_data = '/path_to_local/results/preprocessing/M3';
path_results = '/path_to_local/results/source_reconstruction';
path_atlas = '/path_to_software/atlas/Schaefer_2018/MNI';

atlas_parcels = [200];

% settings
ft_defaults;

% n = 15 (zthr 10) - exclude 2 subjects without MRI (14 and 41)
subjects = [5 7 18 23 25 31 34 37 40 45 47 53 61 201 202];
sessions = 1:2;

%% check linear registrations

for sub_i = 1:length(subjects)

    subject = subjects(sub_i);

    for ses_i = 1:length(sessions)

        session = sessions(ses_i);

        % set paths
        if subject < 10
            subject_ID = ['sub_0' num2str(subject)];
            subjectID = ['sub0' num2str(subject)];
        else
            subject_ID = ['sub_' num2str(subject)];
            subjectID = ['sub' num2str(subject)];
        end
        
        session_ID = ['sess_0' num2str(session)];

        disp([subject_ID ' ' session_ID]);

        for number_parcels_i = 1:length(atlas_parcels)

            number_parcels = atlas_parcels(number_parcels_i);

            for run_i = 0:3
        
                if run_i == 0
        
                    filename = 'localizer';
        
                else
        
                    filename = ['task' num2str(run_i)];
        
                end
        
                t1 = [path_results '/' subject_ID '/' session_ID '/t1_aln_' filename '.nii.gz'];
                atlas_individual = [path_results '/' subject_ID '/' session_ID '/atlas_' filename '_Schaefer2018_' num2str(number_parcels) 'Parcels_7Networks_order_FSLMNI152_1mm.nii.gz'];
    
                system(['fsleyes ' t1 ' ' atlas_individual ' --cmap cortical -a 70']);

            end

        end

    end

end


%% check nonlinear registrations

for sub_i = 1:length(subjects)

    subject = subjects(sub_i);

    for ses_i = 1:length(sessions)

        session = sessions(ses_i);

        % set paths
        if subject < 10
            subject_ID = ['sub_0' num2str(subject)];
            subjectID = ['sub0' num2str(subject)];
        else
            subject_ID = ['sub_' num2str(subject)];
            subjectID = ['sub' num2str(subject)];
        end
        
        session_ID = ['sess_0' num2str(session)];

        disp([subject_ID ' ' session_ID]);

        for number_parcels_i = 1:length(atlas_parcels)

            number_parcels = atlas_parcels(number_parcels_i);

            for run_i = 0:3
        
                if run_i == 0
        
                    filename = 'localizer';
        
                else
        
                    filename = ['task' num2str(run_i)];
        
                end
        
                t1 = [path_results '/' subject_ID '/' session_ID '/t1_aln_' filename '.nii.gz'];
                atlas_individual = [path_results '/' subject_ID '/' session_ID '/atlas_nonlinear_' filename '_Schaefer2018_' num2str(number_parcels) 'Parcels_7Networks_order_FSLMNI152_1mm.nii.gz'];
    
                system(['fsleyes ' t1 ' ' atlas_individual ' --cmap cortical -a 70']);

            end

        end

    end

end

