%% GERE project
% replay

% path to data
if isfolder('/path_to_local')
    path_data = '/path_to_local/results/preprocessing/M3';
    path_inputs = ['/path_to_local/results/mvpa/' folder_decoding];
    path_results = ['/path_to_local/results/replay/' folder_replay];
    path_sequences = '/path_to_local/results/replay/sequences';
    classifiers = {'l2_e-2'};
    inputs = {'bins'};
    state_spaces = {'prob'};
    addpath('/path_to_local/scripts/replay/utilities')
end

number_stimuli = 8;
folders = {'task_stim2delay_balanced_trial_predictions_averaged'...
           'localizer2delay_balanced_trial_predictions_averaged'};
folder_output = 'TDLM_sequence';
time_windows = {'whole_delay'};
sessions = 1:2;

%% replay

errors = cell(0);

for folder_i = 1:length(folders)

    predictions_folder = folders{folder_i};

    try

        for time_window_i = 1:length(time_windows)
            
            for state_i = 1:length(state_spaces)
                                    
                for input_i = 1:length(inputs)
        
                    if strcmp(inputs{input_i}, 'timepoints')
                        nlags = 1000; % number of different lags (500 medians that shifts include [1...500] timepoints of the period of interest (delay period)
                    elseif strcmp(inputs{input_i}, 'bins')
                        nlags = 200;
                    end
                    
                    if strcmp(inputs{input_i}, 'timepoints')
                        if strcmp(time_windows{time_window_i}, 'whole_delay')
                            time_window = 1:4000;
                        elseif strcmp(time_windows{time_window_i}, 'early_delay')
                            time_window = 1:2000;
                        elseif strcmp(time_windows{time_window_i}, 'late_delay')
                            time_window = 2001:4000;
                        end
                    elseif strcmp(inputs{input_i}, 'bins')
                        if strcmp(time_windows{time_window_i}, 'whole_delay')
                            time_window = 1:791;
                        elseif strcmp(time_windows{time_window_i}, 'early_delay')
                            time_window = 1:400;
                        elseif strcmp(time_windows{time_window_i}, 'late_delay')
                            time_window = 401:791;
                        end
                    end
        
                    for class_i = 1:length(classifiers)
            
                        classifier = classifiers{class_i};

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

                            if ~isfolder([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}])
                                mkdir([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}]);
                            end

                            X_l3 = nan(0,3);
                            X_l4 = nan(0,4);
                            X_lall = nan(0,4);

                            X_l3_trials = zeros(1,2); % trials per session
                            X_l4_trials = zeros(1,2);
                            X_lall_trials = zeros(1,2);
                        
                            for ses_i = 1:length(sessions)
                        
                                session = sessions(ses_i);
                                                                                
                                session_ID = ['sess_0' num2str(session)];
        
                                disp([subject_ID ' ' session_ID]);
                            
                                %% load data_task
                            
                                % load
                                load([path_data '/' subject_ID '/sess_0' num2str(session)...
                                        '/GERE_task_zthr25_wholetrial_baseline1sec_postICA.mat']); % data
                        
                                data_task = data;
                        
                                clear data
                    
                                % subset correct trials
                                subset_trials = data_task.trialinfo(:,3);
                                subset_trials = find(subset_trials == 1);
                                data_task.trial = data_task.trial(subset_trials);
                                data_task.time = data_task.time(subset_trials);
                                data_task.sampleinfo = data_task.sampleinfo(subset_trials, :);
                                data_task.trialinfo = data_task.trialinfo(subset_trials, :);
                                                        
                                % load predictions
                                load([path_inputs '/' subject_ID '/' session_ID '/' predictions_folder '/' inputs{input_i} '/' classifier ... 
                                    '/predictions_singletrial.mat']); % predictions
                
                                % subset time window
                                predictions.dval = predictions.dval(:,time_window,:);
                
                                % sigmoid function
                                predictions.prob = (1 ./ (1 + exp(-predictions.dval)));
        
                                if ~isfolder([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                              folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures'])
                                    mkdir([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                              folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures']);
                                end
                                        
                                %% concatenate trials

                                for trial_i = 1:size(predictions.prob, 1)
                        
                                    % subset predictions for trial in the
                                    % sequence
                                    sequence = data_task.trialinfo(trial_i, 4:7);
                                    sequence = sequence(sequence > 0);
                
                                    %% TDLM - sequence
        
                                    % settings TDLM
                                    nstates = length(sequence);
                                    nbins=nlags+1;
        
                                    if strcmp(state_spaces{state_i}, 'prob')
                            
                                        % input - decoded state space
                                        X = squeeze(predictions.prob(trial_i,:,sequence)); % time by states
        
                                    elseif strcmp(state_spaces{state_i}, 'dval')
                            
                                        % input - decoded state space
                                        X = squeeze(predictions.dval(trial_i,:,sequence)); % time by states
        
                                    end

                                    % append
                                    
                                    if size(X,2) == 3
                                        X_l3 = [X; X_l3];
                                        X = [X zeros(size(X,1),1)];
                                        X_lall = [X; X_lall];

                                        X_l3_trials(1,ses_i) = X_l3_trials(1, ses_i) + 1;
                                        X_lall_trials(1, ses_i) = X_lall_trials(1, ses_i) + 1;

                                    else
                                        X_l4 = [X; X_l4];
                                        X_lall = [X; X_lall];

                                        X_l4_trials(1, ses_i) = X_l4_trials(1, ses_i) + 1;
                                        X_lall_trials(1, ses_i) = X_lall_trials(1, ses_i) + 1;
                                    end

                                end
    
                            end

                            % pool trials of lengths 3 and 4, but deleting the
                            % forth column
                            X_l34 = X_lall(:,1:3);
                            X_l34_trials = X_lall_trials;

                            %% TDLM
                                        
                            % X = state space

                            % transition matrix from sequence (deleting
                            % empty rows/columns, and ordering by
                            % appearence in sequence, e.g. [2 5 4] is
                            % encoded as [1 2 3], where [1==2, 2==5,
                            % 3==4]
                            nstates = 3;
                            T_l3 = zeros(nstates,nstates);
                            for item_i = 2:nstates
                                T_l3(item_i-1, item_i) = 1;
                            end

                            nstates = 4;
                            T_l4 = zeros(nstates,nstates);
                            for item_i = 2:nstates
                                T_l4(item_i-1, item_i) = 1;
                            end

                            % number shuffle for null distribution
                            nshuf_l3 = 6; % maximum permutations of 4 states
                            nshuf_l4 = 24; % maximum permutations of 3 states
        
                            % TDLM
                            if ~isfile([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msb_l3.mat'])

                                [Msf_l3, Msb_l3] = TDLMlme2_2ndGLM_concatenate_nocontamination(X_l3,T_l3,3,nlags,nshuf_l3, X_l3_trials);
                                save([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msf_l3.mat'], 'Msf_l3', '-v7.3');
                                save([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msb_l3.mat'], 'Msb_l3', '-v7.3');
    
                            end

                            if ~isfile([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msb_l4.mat'])

                                [Msf_l4, Msb_l4] = TDLMlme2_2ndGLM_concatenate_nocontamination(X_l4,T_l4,4,nlags,nshuf_l4, X_l4_trials);
                                save([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msf_l4.mat'], 'Msf_l4', '-v7.3');
                                save([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msb_l4.mat'], 'Msb_l4', '-v7.3');
    
                            end

                            if ~isfile([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msb_l34.mat'])
                                
                                [Msf_l34, Msb_l34] = TDLMlme2_2ndGLM_concatenate_nocontamination(X_l34,T_l3,3,nlags,nshuf_l3, X_l34_trials);
                                save([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msf_l34.mat'], 'Msf_l34', '-v7.3');
                                save([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msb_l34.mat'], 'Msb_l34', '-v7.3');
            
                            end

                        end
            
                    end
                
                end
            
            end
        
        end

    catch

        errors{end+1} = predictions_folder;

    end
    
end

disp('Analysis done');

