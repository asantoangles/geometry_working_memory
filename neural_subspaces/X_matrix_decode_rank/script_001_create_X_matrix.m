%% GERE project

folder_outputs = 'X_matrix_decode_rank';

% path
if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = [1:15];
else
    path_root = '/path_to_hpc';
end

path_inputs = [path_root '/results/source_reconstruction'];
path_results = [path_root '/github/results/source_geometry_lm/' folder_outputs];
path_results_template = [path_root '/github/results/source_geometry_lm/X_matrix'];
addpath([path_root '/github/scripts/source_geometry_lm/utilities'])

% settings
segments_start = 1:300:3701;
segments_end = 300:300:4000;
delay_segments = [segments_start; segments_end]';
delay_segments(end+1,:) = [1 4000];
delay_segments(end+1,:) = [1 2000];
delay_segments(end+1,:) = [2001 4000];

stim_segments = [];
stim_segments(end+1,:) = [1 300];
stim_segments(end+1,:) = [400 700];
stim_segments(end+1,:) = [800 1100];
stim_segments(end+1,:) = [1200 1500];

number_parcels = 200;
sessions = 1:2;
number_locations = 8;
performance = {'correct_trials'};
events = {'delay', 'stim'};
functions_data ={@mean};
iterations_bootstrap = 1000;
fun_i = @mean;

%% create X matrices 

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

    % loop over events
    for event_i = 1:length(events)
    
        if strcmp(events{event_i}, 'delay')
    
            time_segments = delay_segments;
    
        elseif strcmp(events{event_i}, 'stim')
            
            time_segments = stim_segments;
    
        end

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
        
                time_segments = stim_segments;
        
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
        
            %%% create X matrices
                
            % loop over time windows
            for seg_i = 1:size(time_segments, 1)
    
                seg_start = time_segments(seg_i, 1);
                seg_end   = time_segments(seg_i, 2);
            
                % Crop segment
                data_segment = data_task_event;
                for trial_i = 1:length(data_segment.trial)
                    data_segment.time{trial_i}  = data_segment.time{trial_i}(1:seg_end);
                    data_segment.trial{trial_i} = data_segment.trial{trial_i}(:,seg_start:seg_end);
                end
    
                % loop over sequence lengths
                for sequence_length = [3 4]
    
                    if sequence_length == 3
                        seq_name = 'length3';
                        idx = find(data_segment.trialinfo(:,7) == 0);
                    else
                        seq_name = 'length4';
                        idx = find(data_segment.trialinfo(:,7) ~= 0);
                    end

                    % Create output directory
                    outdir = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                        events{event_i} '_' num2str(seg_start) 'to' num2str(seg_end)];
    
                    if ~isfolder(outdir)
                        mkdir(outdir);
                    end
    
                    % Subset the segment trials
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
    
                    filename_outputs = [outdir '/X_matrix_' seq_name '_decode_rank.mat'];
                    
                    if ~isfile(filename_outputs)

                        % split data
                        trials = 1:size(data_task.trial,2);

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
    
                        save(filename_outputs, 'X_decode_rank', '-v7.3');
                            
                    end

                end

            end
        
        end

    end

end

disp('Analysis done');
