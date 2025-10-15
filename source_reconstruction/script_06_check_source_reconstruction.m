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

% n = 15 - exclude 2 subjects without MRI (14 and 41)
subjects = [5 7 18 23 25 31 34 37 40 45 47 53 61 201 202];
sessions = 1:2;


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

            for run_i = 0:1
        
                if run_i == 0
        
                    filename = 'localizer';
        
                else
        
                    filename = 'task';
        
                end
        
                load([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(number_parcels) '_' filename '.mat']); % data_atlas

                tmp = data_atlas.sources_per_parcel{randperm(size(data_atlas.trial,2),1)};
                tmp(tmp > 0) = [];

                disp(['    ' num2str(length(tmp)) ' parcels with no sources - ' filename ' ' num2str(number_parcels)]);

            end

        end

    end

end

%% check dimensionality

atlas_parcels = [200];

source_dimensionality_localizer = nan(length(subjects)*2, length(atlas_parcels), 2);   % third dimension: trials with deficient rank, median rank of deficient
source_dimensionality_task = nan(length(subjects)*2, length(atlas_parcels), 2);

for sub_i = 1:length(subjects)

    subject = subjects(sub_i);

    for ses_i = 1:length(sessions)

        session = sessions(ses_i);

        %% load data

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

        %% symmetric orthogonalization
        % https://github.com/OHBA-analysis/MEG-ROI-nets/blob/master/%2BROInets/run_individual_network_analysis.m
        % line 405

        for atlas_i = 1:length(atlas_parcels)

            %% localizer

            % load data
            load([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(atlas_i)) '_localizer.mat']); % data_atlas
        
            deficient_rank = [];

            for trial_i = 1:size(data_atlas.trial,2)

                % force full-rank
                if rank(data_atlas.trial{trial_i}) < size(data_atlas.trial{trial_i},1)
                    
                    deficient_rank  = [deficient_rank rank(data_atlas.trial{trial_i})];

                end
                    
            end

            disp(['atlas ' num2str(atlas_parcels(atlas_i)) ' localizer: ' num2str(length(deficient_rank)) ' trials / rank: ' num2str(median(deficient_rank))])

            source_dimensionality_localizer(sub_i+(length(subjects)*(ses_i-1)), atlas_i, 1) = length(deficient_rank);
            source_dimensionality_localizer(sub_i+(length(subjects)*(ses_i-1)), atlas_i, 2) = median(deficient_rank);

        end

        for atlas_i = 1:length(atlas_parcels)

            %% task

            % load data
            load([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(atlas_i)) '_task.mat']); % data_atlas
        
            deficient_rank = [];
            
            for trial_i = 1:size(data_atlas.trial,2)

                % force full-rank
                if rank(data_atlas.trial{trial_i}) < size(data_atlas.trial{trial_i},1)
                    
                    deficient_rank  = [deficient_rank rank(data_atlas.trial{trial_i})];

                end
                    
            end

            disp(['atlas ' num2str(atlas_parcels(atlas_i)) ' task: ' num2str(length(deficient_rank)) ' trials / rank: ' num2str(median(deficient_rank))])

            source_dimensionality_task(sub_i+(length(subjects)*(ses_i-1)), atlas_i, 1) = length(deficient_rank);
            source_dimensionality_task(sub_i+(length(subjects)*(ses_i-1)), atlas_i, 2) = median(deficient_rank);

        end

    end

end

source_dimensionality = [];
source_dimensionality.localizer = source_dimensionality_localizer;
source_dimensionality.task = source_dimensionality_task;

save([path_results '/source_dimensionality.mat'],'source_dimensionality');


disp('Analysis done');
