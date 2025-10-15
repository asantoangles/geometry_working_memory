%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = '/path_to_local/results/source_reconstruction';
    path_results = ['/path_to_local/results/source_geometry_lm/' folder_geometry];
    addpath('/path_to_local/scripts/source_geometry_lm/utilities')
    functions_data ={@mean};
    hpc = 0;
    folder_original = 'LB20';
    path_results_original = ['/path_to_local/results/source_geometry_lm/' folder_original];
end

% time windows
segments_start = 1:300:3701;
segments_end = 300:300:4000;
delay_segments = [segments_start; segments_end]';

iterations_bootstrap = 1000;
number_parcels = 200;
sessions = 1:2;
components = 1:3;

number_locations = 8;

performance = {'correct_trials', 'incorrect_trials'};
events = {'delay'};

%% geometry of memory representations 

for event_i = 1:length(events)

    if strcmp(events{event_i}, 'delay')

        time_segments = delay_segments;

    elseif strcmp(events{event_i}, 'stim_resolved')
        
        time_segments = stim_resolved_segments;

    end

    for sub_i = 1:length(subjects)
    
        subject = subjects(sub_i);

        % set paths
        if subject < 10
            subject_ID = ['sub_0' num2str(subject)];
            subjectID = ['sub0' num2str(subject)];
        else
            subject_ID = ['sub_' num2str(subject)];
            subjectID = ['sub' num2str(subject)];
        end
        
        disp(' '); disp(subject_ID);
    
        % loop over functions
        for fun_i = functions_data

            fun_i = fun_i{1};
            
            % subset by time windows within delay period
            for delay_i = 1:size(time_segments, 1)

                %% geometric analysis
        
                % geometric analysis for each performance
                for perf_i = 1:length(performance)

                    % loop over sequence length used
                    for sequence_length = [3 4 34]
    
                        if sequence_length == 3
                            sequence_length_filename = 'length3';
                            ranks = 3;
                        elseif sequence_length == 4
                            sequence_length_filename = 'length4';
                            ranks = 4;
                        elseif sequence_length == 34
                            sequence_length_filename = 'lengthall';
                            ranks = 4;
                        end

                        % load X matrix
                        load([path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                            events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) ...
                            '/X_matrix_' sequence_length_filename '.mat']); % X
                        
                        % load X matrix for correct trials where to compute eigenvectors
                        X_correct_trials = load([path_results_original '/' subject_ID '/' func2str(fun_i) '/correct_trials/' ...
                            events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) ...
                            '/X_matrix_' sequence_length_filename '.mat']);
                        X_correct_trials = X_correct_trials.X;

                        %% randomization

                        if ~isfolder([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                                events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))])
                            mkdir([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                                events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);
                        end

                        rng('shuffle');

                        X_surrogate = X;
                        for j = 1:size(X_surrogate, 2)
                            X_surrogate(:, j) = X_surrogate(randperm(size(X_surrogate, 1)), j);
                        end
                        X = X_surrogate;

                        X_surrogate = X_correct_trials;
                        for j = 1:size(X_surrogate, 2)
                            X_surrogate(:, j) = X_surrogate(randperm(size(X_surrogate, 1)), j);
                        end
                        X_correct_trials = X_surrogate;

                        %% compute PCA, best-fitting planes, principal angles and variance metrics

                        if ~isfile([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/vaf_' sequence_length_filename '.mat'])
                        
                            % pca, best-fitting planes and angles - pca
                            % on all trials data, and PC scores (data
                            % projection on PCs) done with perforance
                            % specific X matrix
                            [z_k, z_k_exp, planes] = compute_plane_components_alltrials(X_correct_trials, X, components);
    
                            % compute angle and vaf
                            [angle, vaf] = compute_angle_vaf(planes, ranks);
    
                            low_dim_space = [];
                            low_dim_space.z_k = z_k;
                            low_dim_space.plane_r1 = planes.plane_r1;
                            low_dim_space.plane_r2 = planes.plane_r2;
                            low_dim_space.plane_r3 = planes.plane_r3;
                            if sequence_length ~= 3
                                low_dim_space.plane_r4 = planes.plane_r4;
                            end
    
                            save([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/low_dim_space_' sequence_length_filename '.mat'], 'low_dim_space');
    
                            save([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/angle_' sequence_length_filename '.mat'], 'angle');
    
                            save([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/vaf_' sequence_length_filename '.mat'], 'vaf');
    
                            variance_explained = [];
                            variance_explained.z_k = z_k_exp;
                            variance_explained.explained_r1 = planes.explained_r1;
                            variance_explained.explained_r2 = planes.explained_r2;
                            variance_explained.explained_r3 = planes.explained_r3;
                            if sequence_length ~= 3
                                variance_explained.explained_r4 = planes.explained_r4;
                            end
    
                            save([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/variance_explained_' sequence_length_filename '.mat'], 'variance_explained');
        
                        end
    
    
                        %% separability
        
                        if ~isfile([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                            '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                            '/separability_' sequence_length_filename '.mat'])
    
                            load([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/low_dim_space_' sequence_length_filename '.mat']);
        
                            separability = [];
                            separability.volume = cell(4,1);
                            separability.distance_by_separation = cell(4,1);
            
                            for rank_i = 1:ranks
                                    
                                rank = (number_locations*(rank_i-1))+(1:number_locations);
                                points = zscore(low_dim_space.z_k(rank, 1:3));
    
                                [separability.volume{rank_i}, separability.distance_by_separation{rank_i}] = compute_distance(points);
            
                            end
            
                            save([path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/separability_' sequence_length_filename '.mat'], 'separability');
        
                        end
                
                    end

                end

            end

        end

    end

end

disp('Analysis done');
