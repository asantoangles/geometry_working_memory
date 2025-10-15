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

performance = {'correct_trials'};
iterations_bootstrap = 1000;
permutations = 1000;
sessions = 1:2;
number_parcels = 200;
components = 1:3;

number_locations = 8;

events = {'delay'};
delay_segments = delay_segments(subset_i,:);

%% balanced number trials - correct trials only 

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

    for event_i = 1:length(events)
    
        if strcmp(events{event_i}, 'delay')
    
            time_segments = delay_segments;
        
        elseif strcmp(events{event_i}, 'stim_resolved')
            
            time_segments = stim_resolved_segments;
    
        end

        for perf_i = 1:length(performance)

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
    
                number_incorrect_trials = 0;
                number_correct_trials = 0;
                
                for ses_i = 1:length(sessions)
            
                    session = sessions(ses_i);
                                    
                    session_ID = ['sess_0' num2str(session)];
                                        
                    % load data
                    load([path_inputs '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(number_parcels) '_task_zthr25_ortho.mat']); % data_atlas
                    data_atlas.label = data_atlas.label(1:number_parcels);
                    data_task = data_atlas;
                    clear data_atlas

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
    
                    % find number of incorrect trials
                    subset_incorrect_trials = data_task.trialinfo(:,3);
                    subset_incorrect_trials = find(subset_incorrect_trials == 0); % incorrect
                    number_incorrect_trials = number_incorrect_trials + length(subset_incorrect_trials);
    
                    % find number of correct trials
                    subset_correct_trials = data_task.trialinfo(:,3);
                    subset_correct_trials = find(subset_correct_trials == 1); % correct
                    number_correct_trials = number_correct_trials + length(subset_correct_trials);

                    if strcmp(performance{perf_i}, 'correct_trials')
        
                        % subset correct trials
                        subset_trials = data_task.trialinfo(:,3);
                        subset_trials = find(subset_trials == 1); % correct
                        data_task.trial = data_task.trial(subset_trials);
                        data_task.time = data_task.time(subset_trials);
                        data_task.trialinfo = data_task.trialinfo(subset_trials, :);
        
                    elseif strcmp(performance{perf_i}, 'incorrect_trials')
        
                        % subset incorrect trials
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
                
                if strcmp(events{event_i}, 'delay')
            
                    time_segments = delay_segments;
    
                    % subset delay period
                    for trial_i = 1:size(data_task.trial, 2)
                        data_task.time{trial_i} = data_task.time{trial_i}(1:4000);
                        data_task.trial{trial_i} = data_task.trial{trial_i}(:,(end-3999):end);
                    end
            
                elseif strcmp(events{event_i}, 'stim_resolved')
            
                    time_segments = stim_resolved_segments;
    
                    % subset stim period
                    for trial_i = 1:size(data_task.trial, 2)
                        data_task.time{trial_i} = data_task.time{trial_i}(1:1500);
                        data_task.trial{trial_i} = data_task.trial{trial_i}(:,1001:2500);
                    end
                    
                end
        
                data_task_whole_delay = data_task;
                data_baseline_before_perm = data_baseline;
        
                % loop over functions
                for fun_i = functions_data
    
                    fun_i = fun_i{1};
                    
                    % subset by time windows within delay period
                    for delay_i = 1:size(time_segments, 1)
    
                        disp([events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);
        
                        data_task_subset_time = data_task_whole_delay;

                        for trial_i = 1:size(data_task_subset_time.trial, 2)
                            data_task_subset_time.time{trial_i} = data_task_subset_time.time{trial_i}(1:time_segments(delay_i,2));
                            data_task_subset_time.trial{trial_i} = data_task_subset_time.trial{trial_i}(:,time_segments(delay_i,1):time_segments(delay_i,2));
                        end

                        if ~isfolder([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))])
                            mkdir([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);
                        end

                    
                        %% compute PCA, best-fitting planes, principal angles and variance metrics
        
                        if ~isfile([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/separability_perm_' sequence_length_filename '.mat'])
        
                            angle_perm = nan(ranks,ranks,permutations);
                            vaf_perm = nan(ranks,ranks,permutations);
                            euclidean_distances_perm = zeros(number_locations * ranks, permutations);
                            separability_perm = cell(permutations,1);
                            separability_allpairs_perm = cell(permutations,1);

                            % compute matrix X for correct trials

                            if number_locations == 8
    
                                X_before_perm = compute_X_lm_baseline_corrected(data_task_subset_time, data_baseline_before_perm, sequence_length, number_locations);                            
                                X_before_perm = X_before_perm - mean(X_before_perm);
    
                            elseif number_locations == 4
    
                                X1 = compute_X_lm_baseline_corrected(data_task_subset_time, data_baseline_before_perm, sequence_length, number_locations); 
                                X2 = compute_X_lm_baseline_corrected(data_task_subset_time, data_baseline_before_perm, sequence_length, number_locations, 1); 
                                X_before_perm = (X1 + X2) ./ 2;
                                X_before_perm = X_before_perm - mean(X_before_perm);
    
                            end

                            for perm_i = 1:permutations

                                data_task = data_task_subset_time;

                                % subset trials
                                number_trials_bootstrap = number_incorrect_trials;
                                rng('shuffle');
                                idx = randperm(size(data_task.trial,2), number_trials_bootstrap);
                                data_task.trial = data_task.trial(idx);
                                data_task.time = data_task.time(idx);
                                data_task.trialinfo = data_task.trialinfo(idx, :);

                                data_baseline = data_baseline_before_perm;
                                data_baseline.trial = data_baseline.trial(idx);
                                data_baseline.time = data_baseline.time(idx);
                                data_baseline.trialinfo = data_baseline.trialinfo(idx, :);

                                % compute matrix X
                                if number_locations == 8
        
                                    X = compute_X_lm_baseline_corrected(data_task, data_baseline, sequence_length, number_locations);                            
                                    X = X - mean(X);
        
                                elseif number_locations == 4
        
                                    X1 = compute_X_lm_baseline_corrected(data_task, data_baseline, sequence_length, number_locations); 
                                    X2 = compute_X_lm_baseline_corrected(data_task, data_baseline, sequence_length, number_locations, 1); 
                                    X = (X1 + X2) ./ 2;
                                    X = X - mean(X);
        
                                end
                    
                                % pca, best-fitting planes and angles
                                [z_k, z_k_exp, planes] = compute_plane_components_alltrials(X_before_perm, X, components);
    
                                % compute angle and vaf
                                [angle_perm(:,:,perm_i), vaf_perm(:,:,perm_i)] = compute_angle_vaf(planes, ranks);

                                low_dim_space = [];
                                low_dim_space.z_k = z_k;
                                low_dim_space.plane_r1 = planes.plane_r1;
                                low_dim_space.plane_r2 = planes.plane_r2;
                                low_dim_space.plane_r3 = planes.plane_r3;
                                if sequence_length ~= 3
                                    low_dim_space.plane_r4 = planes.plane_r4;
                                end

                                %% separability

                                separability = [];
                                separability.volume = cell(4,1);
                                separability.distance_by_separation = cell(4,1);
                
                                for rank_i = 1:ranks
                                        
                                    rank = (number_locations*(rank_i-1))+(1:number_locations);
                                    points = zscore(low_dim_space.z_k(rank, 1:3));
    
                                    [separability.volume{rank_i}, separability.distance_by_separation{rank_i}] = compute_distance(points);
                
                                end

                                separability_perm{perm_i} = separability;

                            end

                            % save
                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/angle_perm_' sequence_length_filename '.mat'], 'angle_perm');

                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/vaf_perm_' sequence_length_filename '.mat'], 'vaf_perm');

                            save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/separability_perm_' sequence_length_filename '.mat'], 'separability_perm');
                                    
                        end
    
                    end
    
                end
    
            end
        
        end
    
    end

end

disp('Analysis done');
