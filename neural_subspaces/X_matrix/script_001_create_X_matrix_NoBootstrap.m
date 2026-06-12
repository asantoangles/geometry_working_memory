%% GERE project

folder_outputs = 'X_matrix';

% path
if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = [1:15];
else
    path_root = '/path_to_hpc';
end

path_inputs = [path_root '/results/source_reconstruction'];
path_results = [path_root '/github/results/source_geometry_lm/' folder_outputs];
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
performance = {'correct_trials', 'incorrect_trials'};
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
        
                    %% X matrix

                    filename_outputs = [outdir '/X_matrix_' seq_name '.mat'];
                       
                    if ~isfile(filename_outputs)
                    
                        % Compute X
                        if number_locations == 8
        
                            X = compute_X_matrix_lm( ...
                                data_task, data_baseline, sequence_length, number_locations);
                            X = X - mean(X);

                        else  % number_locations == 4
        
                            X1 = compute_X_matrix_lm( ...
                                data_task, data_baseline, sequence_length, number_locations);
                            X2 = compute_X_matrix_lm( ...
                                data_task, data_baseline, sequence_length, number_locations, 1);
        
                            X = (X1 + X2) ./ 2;
                            X = X - mean(X);

                        end
        
                        % subset X matrix of stim period

                        if strcmp(events{event_i}, 'stim')

                            if sequence_length == 3
                                if time_segments(seg_i,1) == 1 && time_segments(seg_i,2) == 300
                                    X = X(1:(number_locations*1),:);
                                elseif time_segments(seg_i,1) == 400 && time_segments(seg_i,2) == 700
                                    X = X(1:(number_locations*2),:);
                                elseif time_segments(seg_i,1) == 800 && time_segments(seg_i,2) == 1100
                                    X = X(1:(number_locations*3),:);
                                elseif time_segments(seg_i,1) == 1200 && time_segments(seg_i,2) == 1500
                                    X = [];
                                end
    
                            else
                                if time_segments(seg_i,1) == 1 && time_segments(seg_i,2) == 300
                                    X = X(1:(number_locations*1),:);
                                elseif time_segments(seg_i,1) == 400 && time_segments(seg_i,2) == 700
                                    X = X(1:(number_locations*2),:);
                                elseif time_segments(seg_i,1) == 800 && time_segments(seg_i,2) == 1100
                                    X = X(1:(number_locations*3),:);
                                elseif time_segments(seg_i,1) == 1200 && time_segments(seg_i,2) == 1500
                                    X = X(1:(number_locations*4),:);
                                end
                            end

                        end
                                
                        % Save X matrix
                        save(filename_outputs, 'X', '-v7.3');

                    else

                        load(filename_outputs); % X

                    end

                end

            end
        
        end

    end

end

disp('Analysis done');
