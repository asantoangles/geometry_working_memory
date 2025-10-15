%% source reconstruction fieldtrip
% https://www.youtube.com/watch?v=Ez72OFjSABs&ab_channel=fieldtriptoolbox

% extract timecourse of virtual channels (sources)
% https://www.fieldtriptoolbox.org/workshop/paris2019/handson_sourceanalysis/

clear variables; close all;

path_utilities = '/path_to_local/scripts/source_reconstruction/utilities';
addpath(path_utilities);
addpath('/path_to_software/fieldtrip-lite-2020083');
addpath /path_to_software/MEG-ROI-nets-master


path_mri = '/path_to_local/data/MRI';
path_head = '/path_to_local/data/head_shape';
path_meg = '/path_to_local/data/MEG';
path_preproc_data = '/path_to_local/results/preprocessing/M3';
path_results = '/path_to_local/results/source_reconstruction';
path_scripts = '/path_to_local/scripts/source_reconstruction';
path_atlas = '/path_to_software/atlas/Schaefer_2018/MNI';

atlas_parcels = [200];

% settings
ft_defaults;

% n = 15 - exclude 2 subjects without MRI (14 and 41)
subjects = [5 7 18 23 25 31 34 37 40 45 47 53 61 201 202];
sessions = 1:2;

%% symmetric orthogonalization

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

            if ~isfile([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(atlas_i)) '_localizer_zthr25_ortho.mat'])

                % load data
                load([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(atlas_i)) '_localizer_zthr25.mat']); % data_atlas
            
                for trial_i = 1:size(data_atlas.trial,2)
    
                    % force full-rank
                    if rank(data_atlas.trial{trial_i}) < size(data_atlas.trial{trial_i},1)
        
                        A = data_atlas.trial{trial_i};
        
                        % Perform Singular Value Decomposition (SVD)
                        [U, S, V] = svd(A, 'econ');
                        
                        % Get the singular values
                        singular_values = diag(S);
                        
                        % Tolerance for determining zero or near-zero singular values
                        tolerance = 1e-10;
                        
                        % Find indices of near-zero singular values
                        near_zero_indices = find(abs(singular_values) < tolerance);
                        
                        % Modify near-zero singular values to some small value
                        singular_values(near_zero_indices) = 1e-10;
                        
                        % Reconstruct the modified matrix
                        A_modified = U * diag(singular_values) * V';
    
                        data_atlas.trial{trial_i} = A_modified;
    
                    end
    
                    data_atlas.trial{trial_i} = ROInets.remove_source_leakage(ROInets.demean(data_atlas.trial{trial_i},2), 'symmetric');
                        
                end
    
                save([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(atlas_i)) '_localizer_zthr25_ortho.mat'], 'data_atlas'); % data_atlas
    
                clear data_atlas

            end

            %% task

            if ~isfile([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(atlas_i)) '_task_zthr25_ortho.mat'])

                % load data
                load([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(atlas_i)) '_task_zthr25.mat']); % data_atlas
            
                for trial_i = 1:size(data_atlas.trial,2)
    
                    % force full-rank
                    if rank(data_atlas.trial{trial_i}) < size(data_atlas.trial{trial_i},1)
        
                        A = data_atlas.trial{trial_i};
        
                        % Perform Singular Value Decomposition (SVD)
                        [U, S, V] = svd(A, 'econ');
                        
                        % Get the singular values
                        singular_values = diag(S);
                        
                        % Tolerance for determining zero or near-zero singular values
                        tolerance = 1e-10;
                        
                        % Find indices of near-zero singular values
                        near_zero_indices = find(abs(singular_values) < tolerance);
                        
                        % Modify near-zero singular values to some small value
                        singular_values(near_zero_indices) = 1e-10;
                        
                        % Reconstruct the modified matrix
                        A_modified = U * diag(singular_values) * V';
    
                        data_atlas.trial{trial_i} = A_modified;
    
                    end
    
                    data_atlas.trial{trial_i} = ROInets.remove_source_leakage(ROInets.demean(data_atlas.trial{trial_i},2), 'symmetric');
                        
                end
    
                save([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(atlas_i)) '_task_zthr25_ortho.mat'], 'data_atlas'); % data_atlas
    
                clear data_atlas

            end

        end

    end

end

disp('Analysis done');
