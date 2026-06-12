%% GERE project

folder_outputs = 'LB23_decode_rank_highdim';

% path
if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = [1:15];
else
    path_root = '/path_to_hpc';
    subjects = [1:15];
end

path_inputs = [path_root '/results/source_reconstruction'];
path_results = [path_root '/github/results/source_geometry_lm/' folder_outputs];
path_results_template = [path_root '/github/results/source_geometry_lm/X_matrix'];
addpath([path_root '/github/scripts/source_geometry_lm/utilities'])

% settings
segments_start = 1:300:3701;
segments_end = 300:300:4000;
delay_segments = [segments_start; segments_end]';

stim_segments = [];
stim_segments(end+1,:) = [1 300];
stim_segments(end+1,:) = [400 700];
stim_segments(end+1,:) = [800 1100];
stim_segments(end+1,:) = [1200 1500];

number_parcels = 200;
sessions = 1:2;
number_locations = 8;
performance = {'correct_trials'};
events = {'stim', 'delay'};
functions_data ={@mean};
fun_i = @mean;
iterations_bootstrap = 500;
number_perm_decoding = 500;
components = 1:3;

% variables to define in sbatch scripts
% iterations_vector = 1
% length_vector = 3

% loop over sequence lengths
for sequence_length = length_vector

    if sequence_length == 3
        seq_name = 'length3';
    else
        seq_name = 'length4';
    end

    for global_iter_i = iterations_vector
    
        if ~isfile([path_results '/decoding_results_' seq_name '_globaliter' num2str(global_iter_i) '.mat'])
    
            % outputs subject-level
            z_all_subjects = nan(length(subjects), size(delay_segments, 1) + 1);
            accuracy_all_subjects = nan(length(subjects), size(delay_segments, 1) + 1);
            accuracy_perm_all_subjects = nan(length(subjects), size(delay_segments, 1) + 1);
        
            % loop over subjects
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
            
                segment_idx = 0;
                
                % loop over events
                for event_i = 1:length(events)
                            
                    % loop over performance
                    for perf_i = 1:length(performance)
                    
                        is_correct  = strcmp(performance{perf_i}, 'correct_trials');
                        target_value = is_correct;                 % correct=1, incorrect=0
                    
                        %%% merge sessions
            
                        data_task_allsessions = [];
                    
                        for ses_i = 1:length(sessions)
                    
                            session    = sessions(ses_i);
                            session_ID = ['sess_0' num2str(session)];
                    
                            % Load data
                            load([path_inputs '/' subject_ID '/' session_ID ...
                                '/data_nonlinear_atlas' num2str(number_parcels) '_task_zthr25_ortho.mat']);
                    
                            % Trim labels
                            data_atlas.label = data_atlas.label(1:number_parcels);
                            data_task = data_atlas; 
                            clear data_atlas;
                    
                            % Subset correct/incorrect
                            idx = find(data_task.trialinfo(:,3) == target_value);
                            data_task.trial     = data_task.trial(idx);
                            data_task.time      = data_task.time(idx);
                            data_task.trialinfo = data_task.trialinfo(idx,:);
                    
                            % Merge sessions
                            if ses_i == 1
                                data_task_allsessions = data_task;
                            else
                                data_task_allsessions.trial     = [data_task_allsessions.trial,     data_task.trial];
                                data_task_allsessions.time      = [data_task_allsessions.time,      data_task.time];
                                data_task_allsessions.trialinfo = [data_task_allsessions.trialinfo; data_task.trialinfo];
                            end
                        end
                    
                        %%% baseline
                        data_baseline_perf = data_task_allsessions;
                    
                        for trial_i = 1:length(data_baseline_perf.trial)
                            data_baseline_perf.time{trial_i}  = data_baseline_perf.time{trial_i}(1:1000);
                            data_baseline_perf.trial{trial_i} = data_baseline_perf.trial{trial_i}(:,1:1000);
                        end
                    
                        %%% delay
                        data_task_event = data_task_allsessions;
                    
                        if strcmp(events{event_i}, 'stim')
                    
                            time_segments = stim_segments(sequence_length,:);
                    
                            for trial_i = 1:length(data_task_event.trial)
                                data_task_event.time{trial_i}  = data_task_event.time{trial_i}(1:1500);
                                data_task_event.trial{trial_i} = data_task_event.trial{trial_i}(:,1001:2500);
                            end
                    
                        else  % delay
                            
                            time_segments = delay_segments;
                    
                            for trial_i = 1:length(data_task_event.trial)
                                L = length(data_task_event.time{trial_i});
                                data_task_event.time{trial_i}  = data_task_event.time{trial_i}(1:4000);
                                data_task_event.trial{trial_i} = data_task_event.trial{trial_i}(:,L-3999:L);
                            end
                        end
                                            
                        % loop over time windows
                        for seg_i = 1:size(time_segments, 1)
                
                            seg_start = time_segments(seg_i, 1);
                            seg_end   = time_segments(seg_i, 2);

                            segment_idx = segment_idx + 1;

                            % set stim_order variable
                            if strcmp(events{event_i}, 'stim')
            
                                if sequence_length == 3
                                    if time_segments(seg_i,1) == 1 && time_segments(seg_i,2) == 300
                                        stim_order = 1;
                                    elseif time_segments(seg_i,1) == 400 && time_segments(seg_i,2) == 700
                                        stim_order = 2;
                                    elseif time_segments(seg_i,1) == 800 && time_segments(seg_i,2) == 1100
                                        stim_order = 3;
                                    elseif time_segments(seg_i,1) == 1200 && time_segments(seg_i,2) == 1500
                                        stim_order = 3;
                                    end
                                elseif sequence_length == 4
                                    if time_segments(seg_i,1) == 1 && time_segments(seg_i,2) == 300
                                        stim_order = 1;
                                    elseif time_segments(seg_i,1) == 400 && time_segments(seg_i,2) == 700
                                        stim_order = 2;
                                    elseif time_segments(seg_i,1) == 800 && time_segments(seg_i,2) == 1100
                                        stim_order = 3;
                                    elseif time_segments(seg_i,1) == 1200 && time_segments(seg_i,2) == 1500
                                        stim_order = 4;
                                    end
                                end                            
            
                            else
                                stim_order = sequence_length;
                            end
                        
                            % Crop segment
                            data_segment = data_task_event;
                            for trial_i = 1:length(data_segment.trial)
                                data_segment.time{trial_i}  = data_segment.time{trial_i}(1:seg_end);
                                data_segment.trial{trial_i} = data_segment.trial{trial_i}(:,seg_start:seg_end);
                            end
                
                            % Create output directory
                            outdir = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                                events{event_i} '_' num2str(seg_start) 'to' num2str(seg_end)];
            
                            if ~isfolder(outdir)
                                mkdir(outdir);
                            end
            
                            % Subset the segment trials
                            if sequence_length == 3
                                idx = find(data_segment.trialinfo(:,7) == 0);
                            else
                                idx = find(data_segment.trialinfo(:,7) ~= 0);
                            end

                            % Subset data_task
                            data_task = data_segment;
                            data_task.trial     = data_task.trial(idx);
                            data_task.time      = data_task.time(idx);
                            data_task.trialinfo = data_task.trialinfo(idx,:);
            
                            % Subset baseline
                            data_baseline = data_baseline_perf;
                            data_baseline.trial     = data_baseline.trial(idx);
                            data_baseline.time      = data_baseline.time(idx);
                            data_baseline.trialinfo = data_baseline.trialinfo(idx,:);
        
                            % load X matrix for dimensions
        
                            filename_X_matrix = [path_results_template '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                                events{event_i} '_' num2str(seg_start) 'to' num2str(seg_end) '/X_matrix_' seq_name '.mat'];
                            load(filename_X_matrix); % X
                
                            %% X matrices half split train/test, and resampling within each split
        
                            % split data
                            trials = 1:size(data_task.trial,2);
    
                            rng('shuffle');
        
                            idx_split_train = randperm(length(trials), round(length(trials)/2));
                            data_task_train = data_task;
                            data_task_train.trial = data_task_train.trial(idx_split_train);
                            data_task_train.time = data_task_train.time(idx_split_train);
                            data_baseline_train = data_baseline;
                            data_baseline_train.trial = data_baseline_train.trial(idx_split_train);
                            data_baseline_train.time = data_baseline_train.time(idx_split_train);
        
                            idx_split_test = setdiff(trials, idx_split_train);
                            data_task_test = data_task;
                            data_task_test.trial = data_task_test.trial(idx_split_test);
                            data_task_test.time = data_task_test.time(idx_split_test);
                            data_baseline_test = data_baseline;
                            data_baseline_test.trial = data_baseline_test.trial(idx_split_test);
                            data_baseline_test.time = data_baseline_test.time(idx_split_test);
        
                            X_decode_rank = [];
                            X_decode_rank.X_train = nan(size(X,1), size(X,2), iterations_bootstrap);
                            X_decode_rank.X_test = nan(size(X,1), size(X,2), iterations_bootstrap);
        
                            for iter_i = 1:iterations_bootstrap
        
                                % resample with replacement - train
                                idx_iter = randi(size(data_task_train.trial,2),size(data_task_train.trial,2),1);
                                data_task_train_iter = data_task_train;
                                data_task_train_iter.trial = data_task_train_iter.trial(idx_iter);
                                data_task_train_iter.time = data_task_train_iter.time(idx_iter);
                                data_baseline_train_iter = data_baseline_train;
                                data_baseline_train_iter.trial = data_baseline_train_iter.trial(idx_iter);
                                data_baseline_train_iter.time = data_baseline_train_iter.time(idx_iter);
        
                                % resample with replacement - test
                                idx_iter = randi(size(data_task_test.trial,2),size(data_task_test.trial,2),1);
                                data_task_test_iter = data_task_test;
                                data_task_test_iter.trial = data_task_test_iter.trial(idx_iter);
                                data_task_test_iter.time = data_task_test_iter.time(idx_iter);
                                data_baseline_test_iter = data_baseline_test;
                                data_baseline_test_iter.trial = data_baseline_test_iter.trial(idx_iter);
                                data_baseline_test_iter.time = data_baseline_test_iter.time(idx_iter);
        
        
                                if number_locations == 4
        
                                    X_split1_1 = compute_X_matrix_lm(data_task_train_iter, data_baseline_train, sequence_length, number_locations, 1);
                                    X_split1_2 = compute_X_matrix_lm(data_task_train_iter, data_baseline_train, sequence_length, number_locations);
                                    X_split1 = (X_split1_1 + X_split1_2) / 2;
                                    X_split1 = X_split1 - mean(X_split1);
        
                                    X_split2_1 = compute_X_matrix_lm(data_task_test_iter, data_baseline_test, sequence_length, number_locations, 1);
                                    X_split2_2 = compute_X_matrix_lm(data_task_test_iter, data_baseline_test, sequence_length, number_locations);
                                    X_split2 = (X_split2_1 + X_split2_2) / 2;
                                    X_split2 = X_split2 - mean(X_split2);
        
                                elseif number_locations == 8
        
                                    X_split1 = compute_X_matrix_lm(data_task_train_iter, data_baseline_train, sequence_length, number_locations);
                                    X_split1 = X_split1 - mean(X_split1);
        
                                    X_split2 = compute_X_matrix_lm(data_task_test_iter, data_baseline_test, sequence_length, number_locations);
                                    X_split2 = X_split2 - mean(X_split2);
        
                                end
        
                                % subset X matrix of stim period
                        
                                if strcmp(events{event_i}, 'stim')
        
                                    if sequence_length == 3
                                        if time_segments(seg_i,1) == 1 && time_segments(seg_i,2) == 300
                                            X_split1 = X_split1(1:(number_locations*1),:);
                                            X_split2 = X_split2(1:(number_locations*1),:);
                                        elseif time_segments(seg_i,1) == 400 && time_segments(seg_i,2) == 700
                                            X_split1 = X_split1(1:(number_locations*2),:);
                                            X_split2 = X_split2(1:(number_locations*2),:);
                                        elseif time_segments(seg_i,1) == 800 && time_segments(seg_i,2) == 1100
                                            X_split1 = X_split1(1:(number_locations*3),:);
                                            X_split2 = X_split2(1:(number_locations*3),:);
                                        elseif time_segments(seg_i,1) == 1200 && time_segments(seg_i,2) == 1500
                                            X_split1 = [];
                                            X_split2 = [];
                                        end
            
                                    else
                                        if time_segments(seg_i,1) == 1 && time_segments(seg_i,2) == 300
                                            X_split1 = X_split1(1:(number_locations*1),:);
                                            X_split2 = X_split2(1:(number_locations*1),:);
                                        elseif time_segments(seg_i,1) == 400 && time_segments(seg_i,2) == 700
                                            X_split1 = X_split1(1:(number_locations*2),:);
                                            X_split2 = X_split2(1:(number_locations*2),:);
                                        elseif time_segments(seg_i,1) == 800 && time_segments(seg_i,2) == 1100
                                            X_split1 = X_split1(1:(number_locations*3),:);
                                            X_split2 = X_split2(1:(number_locations*3),:);
                                        elseif time_segments(seg_i,1) == 1200 && time_segments(seg_i,2) == 1500
                                            X_split1 = X_split1(1:(number_locations*4),:);
                                            X_split2 = X_split2(1:(number_locations*4),:);
                                        end
                                    end
        
                                end
        
                                X_decode_rank.X_train(:,:,iter_i) = X_split1;
                                X_decode_rank.X_test(:,:,iter_i) = X_split2;
        
                            end
                                        
                            %% compute subspaces (pc_scores)
                
                            for train_test_i = 1:2
        
                                if train_test_i == 1
        
                                    X_allperm = X_decode_rank.X_train;
        
                                elseif train_test_i == 2
        
                                    X_allperm = X_decode_rank.X_test;
        
                                end
        
                                X_all = [];
                                Y_all = [];
                            
                                for perm_i = 1:size(X_allperm, 3)
                                    
                                    % =========================================================
                                    % compute subspaces (z_k)
                                    % =========================================================
                                    X = squeeze(X_allperm(:,:,perm_i));

                                    z_k = X; % assign full data matrix to subspaces (z_k) variable
                                
                                    % =========================================================
                                    % vectorize subspaces (rank-specific)
                                    % =========================================================
                                    nLoc = number_locations;
                                
                                    rank1 = z_k(1:nLoc, :);             rank1 = rank1(:)';
                                    rank2 = z_k(nLoc+1:2*nLoc, :);      rank2 = rank2(:)';
                                    rank3 = z_k(2*nLoc+1:3*nLoc, :);    rank3 = rank3(:)';
                                
                                    X_all = [X_all; rank1; rank2; rank3];
                                    Y_all = [Y_all; 1; 2; 3];
                                
                                    % optional 4th rank
                                    if sequence_length == 4
                                        rank4 = z_k(3*nLoc+1:4*nLoc, :);
                                        rank4 = rank4(:)';
                                
                                        X_all = [X_all; rank4];
                                        Y_all = [Y_all; 4];
                                    end
        
                                end
        
                                if train_test_i == 1
        
                                    X_train_all = X_all;
                                    Y_train_all = Y_all;
        
                                elseif train_test_i == 2
        
                                    X_test_all = X_all;
                                    Y_test_all = Y_all;
        
                                end
                                                        
                            end
        
                            %% decoding
        
                            % =========================================================
                            % MULTICLASS DECODING (LOGISTIC REGRESSION) + PERMUTATION TEST
                            % =========================================================
        
                            % -------------------------
                            % REAL ACCURACY (LOGISTIC REGRESSION)
                            % -------------------------
                            t = templateLinear('Learner','logistic');
                            
                            % Train on full training set
                            mdl = fitcecoc(X_train_all, Y_train_all, 'Learners', t);
                            
                            % Predict on test set
                            pred = predict(mdl, X_test_all);
                            
                            % Compute accuracy
                            acc_emp = mean(pred == Y_test_all);
        
                            % -------------------------
                            % PERMUTATION TEST
                            % -------------------------
                            t = templateLinear('Learner','logistic');
                            
                            acc_perm = zeros(number_perm_decoding,1);
                            
                            for p = 1:number_perm_decoding
                                
                                % -------------------------
                                % shuffle TRAIN labels only
                                % -------------------------
                                Y_shuff = Y_train_all(randperm(length(Y_train_all)));
                                
                                % -------------------------
                                % train on shuffled labels
                                % -------------------------
                                mdl = fitcecoc(X_train_all, Y_shuff, 'Learners', t);
                                
                                % -------------------------
                                % test on TRUE test labels
                                % -------------------------
                                pred = predict(mdl, X_test_all);
                                
                                acc_perm(p) = mean(pred == Y_test_all);
                                
                            end
                            
                            accuracy_perm_all_subjects(sub_i, segment_idx) = mean(acc_perm);
                                
                            % -------------------------
                            % P-VALUE
                            % -------------------------
                            p_value = mean(acc_perm >= acc_emp);
                            
                            % fprintf('Permutation p-value: %.4f\n', p_value);
                                            
                            % -------------------------
                            % STATS
                            % -------------------------
                            mu_perm = mean(acc_perm);
                            sd_perm = std(acc_perm);
                            z_subj = (acc_emp - mu_perm) / sd_perm;
            
                            z_all_subjects(sub_i, segment_idx) = z_subj;
                            accuracy_all_subjects(sub_i, segment_idx) = acc_emp;
        
        
                        end
        
                    end
                
                end
        
            end
    
            decoding_results = [];
            decoding_results.z_all_subjects = z_all_subjects;
            decoding_results.accuracy_all_subjects = accuracy_all_subjects;
            decoding_results.accuracy_perm_all_subjects = accuracy_perm_all_subjects;
        
            save([path_results '/decoding_results_' seq_name '_globaliter' num2str(global_iter_i) '.mat'], 'decoding_results');
    
        end

    end

end

disp('Analysis done');
