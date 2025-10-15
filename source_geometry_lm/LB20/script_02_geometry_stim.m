%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = '/path_to_local/results/source_reconstruction';
    path_results = ['/path_to_local/results/source_geometry_lm/' folder_geometry];
    addpath('/path_to_local/scripts/source_geometry_lm/utilities')
    functions_data ={@mean};
end

% time windows
segments_start = 1:300:3701;
segments_end = 300:300:4000;
delay_segments = [segments_start; segments_end]';

stim_resolved_segments = [];
stim_resolved_segments(end+1,:) = [1 300];
stim_resolved_segments(end+1,:) = [400 700];
stim_resolved_segments(end+1,:) = [800 1100];
stim_resolved_segments(end+1,:) = [1200 1500];

iterations_bootstrap = 1000;
number_parcels = 200;
sessions = 1:2;
components = 1:3;

number_locations = 8;

performance = {'correct_trials', 'incorrect_trials'};
events = {'stim_resolved_refined'};

%% geometry of memory representations 

for event_i = 1:length(events)

    if strcmp(events{event_i}, 'stim_resolved_refined')
        
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
    
        %% pool data across sessions, separately for performance
        for perf_i = 1:length(performance)

            for ses_i = 1:length(sessions)
        
                session = sessions(ses_i);
                        
                session_ID = ['sess_0' num2str(session)];
                                    
                % load data
                load([path_inputs '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(number_parcels) '_task_zthr25_ortho.mat']); % data_atlas
                data_atlas.label = data_atlas.label(1:number_parcels);
                data_task = data_atlas;
                clear data_atlas
    
                if strcmp(performance{perf_i}, 'correct_trials')
    
                    % subset correct trials
                    subset_trials = data_task.trialinfo(:,3);
                    subset_trials = find(subset_trials == 1); % correct
                    data_task.trial = data_task.trial(subset_trials);
                    data_task.time = data_task.time(subset_trials);
                    data_task.trialinfo = data_task.trialinfo(subset_trials, :);                  
    
                elseif strcmp(performance{perf_i}, 'incorrect_trials')
    
                    % subset correct trials
                    subset_trials = data_task.trialinfo(:,3);
                    subset_trials = find(subset_trials == 0); % incorrect
                    data_task.trial = data_task.trial(subset_trials);
                    data_task.time = data_task.time(subset_trials);
                    data_task.trialinfo = data_task.trialinfo(subset_trials, :);
    
                end
    
                % merge data_task
    
                if ses_i == 1
    
                    data_task_allsessions = data_task;
    
                else
    
                    for trial_i = 1:size(data_task.trial,2)
    
                        data_task_allsessions.trial{1,end+1} = data_task.trial{trial_i};
                        data_task_allsessions.time{1,end+1} = data_task.time{trial_i};
    
                    end
    
                    data_task_allsessions.trialinfo = [data_task_allsessions.trialinfo; data_task.trialinfo];
    
                    data_task = data_task_allsessions; clear data_task_allsessions
    
                end
    
            end
    
            % data_baseline
    
            data_baseline = data_task;
            for trial_i = 1:size(data_baseline.trial, 2)
                data_baseline.time{trial_i} = data_baseline.time{trial_i}(1:1000);
                data_baseline.trial{trial_i} = data_baseline.trial{trial_i}(:,1:1000);
            end
    
            if strcmp(events{event_i}, 'stim_resolved_refined') % same than stim
        
                time_segments = stim_resolved_segments;
    
                % subset stim period
                for trial_i = 1:size(data_task.trial, 2)
                    data_task.time{trial_i} = data_task.time{trial_i}(1:1500);
                    data_task.trial{trial_i} = data_task.trial{trial_i}(:,1001:2500);
                end
            
            end

            if strcmp(performance{perf_i}, 'correct_trials')

                data_task_whole_delay_correct_trials = data_task;
                data_baseline_correct_trials = data_baseline;

            elseif strcmp(performance{perf_i}, 'incorrect_trials')

                data_task_whole_delay_incorrect_trials = data_task;
                data_baseline_incorrect_trials = data_baseline;

            end

        end

        % loop over functions
        for fun_i = functions_data

            fun_i = fun_i{1};
            
            % subset by time windows within delay period
            for delay_i = 1:size(time_segments, 1)

                disp([events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);

                % correct_trials
                data_task_correct_trials = data_task_whole_delay_correct_trials;
                for trial_i = 1:size(data_task_correct_trials.trial, 2)
                    data_task_correct_trials.time{trial_i} = data_task_correct_trials.time{trial_i}(1:time_segments(delay_i,2));
                    data_task_correct_trials.trial{trial_i} = data_task_correct_trials.trial{trial_i}(:,time_segments(delay_i,1):time_segments(delay_i,2));
                end
                
                % incorrect_trials
                data_task_incorrect_trials = data_task_whole_delay_incorrect_trials;
                for trial_i = 1:size(data_task_incorrect_trials.trial, 2)
                    data_task_incorrect_trials.time{trial_i} = data_task_incorrect_trials.time{trial_i}(1:time_segments(delay_i,2));
                    data_task_incorrect_trials.trial{trial_i} = data_task_incorrect_trials.trial{trial_i}(:,time_segments(delay_i,1):time_segments(delay_i,2));
                end
    
                data_task_correct_trials_lengthall = data_task_correct_trials;
                data_task_incorrect_trials_lengthall = data_task_incorrect_trials;

                %% geometric analysis
        
                % geometric analysis for each performance
                for perf_i = 1:length(performance)

                    if strcmp('correct_trials', performance{perf_i})

                        data_task = data_task_correct_trials_lengthall;
                        data_baseline = data_baseline_correct_trials;

                    elseif strcmp('incorrect_trials', performance{perf_i})

                        data_task = data_task_incorrect_trials_lengthall;
                        data_baseline = data_baseline_incorrect_trials;

                    end

                    data_task_before_subset_length = data_task;
                    data_baseline_before_subset_length = data_baseline;

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

                        data_task = data_task_before_subset_length;
                        data_baseline = data_baseline_before_subset_length;

                        % subset data by sequence length
                        if sequence_length == 3
                            idx = find(data_task.trialinfo(:, 7) == 0);
                            data_task.trial = data_task.trial(idx);
                            data_task.time = data_task.time(idx);
                            data_task.trialinfo = data_task.trialinfo(idx,:);
                        elseif sequence_length == 4
                            idx = find(data_task.trialinfo(:, 7) ~= 0);
                            data_task.trial = data_task.trial(idx);
                            data_task.time = data_task.time(idx);
                            data_task.trialinfo = data_task.trialinfo(idx,:);
                        end

                        % subset data_baseline by sequence length
                        if sequence_length == 3
                            idx = find(data_baseline.trialinfo(:, 7) == 0);
                            data_baseline.trial = data_baseline.trial(idx);
                            data_baseline.time = data_baseline.time(idx);
                            data_baseline.trialinfo = data_baseline.trialinfo(idx,:);
                        elseif sequence_length == 4
                            idx = find(data_baseline.trialinfo(:, 7) ~= 0);
                            data_baseline.trial = data_baseline.trial(idx);
                            data_baseline.time = data_baseline.time(idx);
                            data_baseline.trialinfo = data_baseline.trialinfo(idx,:);
                        end

                        % compute X matrix to compute PC scores 
                        if number_locations == 8
    
                            X = compute_X_lm_baseline_corrected(data_task, data_baseline, sequence_length, number_locations);                            
                            X = X - mean(X);
    
                        elseif number_locations == 4
    
                            X1 = compute_X_lm_baseline_corrected(data_task, data_baseline, sequence_length, number_locations); 
                            X2 = compute_X_lm_baseline_corrected(data_task, data_baseline, sequence_length, number_locations, 1); 
                            X = (X1 + X2) ./ 2;
                            X = X - mean(X);
    
                        end

                        % NOTE change here relative to main script
                        % subset X matrix depending on stimuli presented
                        if sequence_length == 3
                            if time_segments(delay_i,1) == 1 && time_segments(delay_i,2) == 300
                                stim_order = 1;
                                X = X(1:(number_locations*1),:);
                            elseif time_segments(delay_i,1) == 400 && time_segments(delay_i,2) == 700
                                stim_order = 1;
                                X = X(1:(number_locations*1),:);
                            elseif time_segments(delay_i,1) == 800 && time_segments(delay_i,2) == 1100
                                stim_order = 2;
                                X = X(1:(number_locations*2),:);
                            elseif time_segments(delay_i,1) == 1200 && time_segments(delay_i,2) == 1500
                                stim_order = 3;
                                X = X(1:(number_locations*3),:);
                            end

                        else
                            if time_segments(delay_i,1) == 1 && time_segments(delay_i,2) == 300
                                X = X(1:(number_locations*1),:);
                                stim_order = 1;
                            elseif time_segments(delay_i,1) == 400 && time_segments(delay_i,2) == 700
                                X = X(1:(number_locations*2),:);
                                stim_order = 2;
                            elseif time_segments(delay_i,1) == 800 && time_segments(delay_i,2) == 1100
                                X = X(1:(number_locations*3),:);
                                stim_order = 3;
                            elseif time_segments(delay_i,1) == 1200 && time_segments(delay_i,2) == 1500
                                X = X(1:(number_locations*4),:);
                                stim_order = 4;
                            end
                        end


                        if ~isfolder([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                                events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))])
                            mkdir([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                                events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);
                        end
    
                        save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                                events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) ...
                            '/X_matrix_' sequence_length_filename '.mat'], ...
                            'X', '-v7.3');

                        % load X matrix for correct trials where to compute
                        % eigenvectors (delay)
                        X_correct_trials = load([path_results '/' subject_ID '/' func2str(fun_i) ...
                            '/correct_trials/delay_1to4000/X_matrix_' sequence_length_filename '.mat']);
                        X_correct_trials = X_correct_trials.X;
                        
                        %% compute PCA, best-fitting planes, principal angles and variance metrics
        
                        if ~isfile([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/vaf_' sequence_length_filename '.mat'])
                        
                            % pca, best-fitting planes and angles - pca
                            % on all trials data, and PC scores (data
                            % projection on PCs) done with perforance
                            % specific X matrix
                            [z_k, z_k_exp, planes] = compute_plane_components_alltrials_stim_refined(X_correct_trials, X, components, stim_order);
    
                            % compute angle and vaf
                            [angle, vaf] = compute_angle_vaf(planes, stim_order);
    
                            low_dim_space = [];
                            low_dim_space.z_k = z_k;
                            low_dim_space.plane_r1 = planes.plane_r1;

                            if stim_order > 1
                                low_dim_space.plane_r2 = planes.plane_r2;
                            end

                            if stim_order > 2
                                low_dim_space.plane_r3 = planes.plane_r3;
                            end

                            if stim_order > 3
                                low_dim_space.plane_r4 = planes.plane_r4;
                            end
    
                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/low_dim_space_' sequence_length_filename '.mat'], 'low_dim_space');
    
                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/angle_' sequence_length_filename '.mat'], 'angle');
    
                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/vaf_' sequence_length_filename '.mat'], 'vaf');
    
                            variance_explained = [];
                            variance_explained.z_k = z_k_exp;
                            variance_explained.explained_r1 = planes.explained_r1;

                            if stim_order > 1
                                variance_explained.explained_r2 = planes.explained_r2;
                            end

                            if stim_order > 2
                                variance_explained.explained_r3 = planes.explained_r3;
                            end

                            if stim_order > 3
                                variance_explained.explained_r4 = planes.explained_r4;
                            end
    
                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/variance_explained_' sequence_length_filename '.mat'], 'variance_explained');
        
                        end
    
    
                        %% separability
        
                        if ~isfile([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                            '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                            '/separability_' sequence_length_filename '.mat'])
    
                            load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/low_dim_space_' sequence_length_filename '.mat']);
        
                            separability = [];
                            separability.volume = cell(4,1);
                            separability.distance_by_separation = cell(4,1);
            
                            for rank_i = 1:stim_order
                                    
                                rank = (number_locations*(rank_i-1))+(1:number_locations);
                                points = zscore(low_dim_space.z_k(rank, 1:3));
    
                                [separability.volume{rank_i}, separability.distance_by_separation{rank_i}] = compute_distance(points);
            
                            end
            
                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/separability_' sequence_length_filename '.mat'], 'separability');
        
                        end
    
                        %% bootstrapp procedure to compute within-subspace metrics
    
                        % split trials within-session, compute planes for
                        % each split, and calculate angle/vaf between
                        % splits
    
                        if ~isfile([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/vaf_bootstrap_' sequence_length_filename '.mat'])
    
                            angle_bootstrap_alliter = nan(stim_order, stim_order, iterations_bootstrap);
                            vaf_bootstrap_alliter = nan(stim_order, stim_order, iterations_bootstrap);
    
                            for iter_i = 1:iterations_bootstrap
    
                                % split data
                                trials = 1:size(data_task.trial,2);
                                idx_split1 = randperm(length(trials), round(length(trials)/2));
                                idx_split2 = setdiff(trials, idx_split1);
    
                                data_task_split1 = data_task;
                                data_task_split1.trial = data_task_split1.trial(idx_split1);
                                data_baseline_split1 = data_baseline;
                                data_baseline_split1.trial = data_baseline_split1.trial(idx_split1);
    
                                data_task_split2 = data_task;
                                data_task_split2.trial = data_task_split2.trial(idx_split2);
                                data_baseline_split2 = data_baseline;
                                data_baseline_split2.trial = data_baseline_split2.trial(idx_split2);
    
                                if number_locations == 4
    
                                    X_split1_1 = compute_X_lm_baseline_corrected(data_task_split1, data_baseline_split1, sequence_length, number_locations, 1);
                                    X_split1_2 = compute_X_lm_baseline_corrected(data_task_split1, data_baseline_split1, sequence_length, number_locations);
                                    X_split1 = (X_split1_1 + X_split1_2) / 2;
                                    X_split1 = X_split1 - mean(X_split1);
    
                                    X_split2_1 = compute_X_lm_baseline_corrected(data_task_split2, data_baseline_split2, sequence_length, number_locations, 1);
                                    X_split2_2 = compute_X_lm_baseline_corrected(data_task_split2, data_baseline_split2, sequence_length, number_locations);
                                    X_split2 = (X_split2_1 + X_split2_2) / 2;
                                    X_split2 = X_split2 - mean(X_split2);
    
                                elseif number_locations == 8
    
                                    X_split1 = compute_X_lm_baseline_corrected(data_task_split1, data_baseline_split1, sequence_length, number_locations);
                                    X_split1 = X_split1 - mean(X_split1);
    
                                    X_split2 = compute_X_lm_baseline_corrected(data_task_split2, data_baseline_split2, sequence_length, number_locations);
                                    X_split2 = X_split2 - mean(X_split2);
    
                                end

                                % NOTE change here relative to main script
                                % subset X matrix depending on stimuli presented
                                if sequence_length == 3
                                    if stim_order == 1
                                        X_split1 = X_split1(1:(number_locations*1),:);
                                        X_split2 = X_split2(1:(number_locations*1),:);
                                    elseif stim_order == 2
                                        X_split1 = X_split1(1:(number_locations*2),:);
                                        X_split2 = X_split2(1:(number_locations*2),:);
                                    elseif stim_order == 3
                                        X_split1 = X_split1(1:(number_locations*3),:);
                                        X_split2 = X_split2(1:(number_locations*3),:);
                                    end
        
                                else
                                    if stim_order == 1
                                        X_split1 = X_split1(1:(number_locations*1),:);
                                        X_split2 = X_split2(1:(number_locations*1),:);
                                    elseif stim_order == 2
                                        X_split1 = X_split1(1:(number_locations*2),:);
                                        X_split2 = X_split2(1:(number_locations*2),:);
                                    elseif stim_order == 3
                                        X_split1 = X_split1(1:(number_locations*3),:);
                                        X_split2 = X_split2(1:(number_locations*3),:);
                                    elseif stim_order == 4
                                        X_split1 = X_split1(1:(number_locations*4),:);
                                        X_split2 = X_split2(1:(number_locations*4),:);
                                    end
                                end

    
                                % compute plane within
                                [planes_split1, planes_split2, Y_r, score] = compute_plane_components_within_alltrials_stim_refined(X_correct_trials, X_split1, X_split2, components, stim_order);
    
                                % compute angle and vaf within
                                [angle_bootstrap_alliter(:,:,iter_i), vaf_bootstrap_alliter(:,:,iter_i)] = compute_angle_vaf_within(planes_split1, planes_split2, Y_r, score, stim_order);
        
                            end
    
                            % save outputs    
                            angle_bootstrap = [];
                            angle_bootstrap.median = median(angle_bootstrap_alliter,3);
                            angle_bootstrap.std = std(angle_bootstrap_alliter,1,3);
    
                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/angle_bootstrap_' sequence_length_filename '.mat'], 'angle_bootstrap', '-v7.3');
    
                            vaf_bootstrap = [];
                            vaf_bootstrap.median = median(vaf_bootstrap_alliter,3);
                            vaf_bootstrap.std = std(vaf_bootstrap_alliter,1,3);
    
                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/vaf_bootstrap_' sequence_length_filename '.mat'], 'vaf_bootstrap', '-v7.3');
                            
                        end
            
                    end

                end

            end

        end

    end

end

disp('Analysis done');
