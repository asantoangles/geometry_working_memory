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
folders = {'task_stim_balanced2task_stim_trial_predictions_averaged' 'localizer2task_stim_trial_predictions_averaged'};
folder_output = 'TDLM_sequence';

%% replay

errors = cell(0);

for folder_i = 1:length(folders)

    predictions_folder = folders{folder_i};

    try

        for state_i = 1:length(state_spaces)
        
            % loop over trials
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
        
                    if ~size(data_task.trial, 2) == 0
        
                        %% set input   
                        for input_i = 1:length(inputs)
                
                            if strcmp(inputs{input_i}, 'timepoints')
                                nlags = 500; % number of different lags (500 medians that shifts include [1...500] timepoints of the period of interest (delay period)
                            elseif strcmp(inputs{input_i}, 'bins')
                                nlags = 100;
                            end
                
                            %% loop across classifiers
                            for class_i = 1:length(classifiers)
                    
                                classifier = classifiers{class_i};
        
                                disp(classifier);
                
                                if ~isfolder([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i}])
                                    mkdir([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i}]);
                                end
                        
                                % load predictions
                                load([path_inputs '/' subject_ID '/' session_ID '/' predictions_folder '/' inputs{input_i} '/' classifier ... 
                                    '/predictions_singletrial.mat']); % predictions

%                                 % subset time in predictions to include
%                                 % only stimulus presentation
%                                 interval_stim = [101:400 501:800 901:1200 1301:1600];
%                                 predictions.dval = predictions.dval(:,interval_stim,:);

                                % sigmoid function
                                predictions.prob = (1 ./ (1 + exp(-predictions.dval)));
        
                                Msf_alltrials = cell(0);
                                Msb_alltrials = cell(0);
                                Msz_alltrials = cell(0);
        
                                if ~isfile([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_length4.mat'])
                                    
                                    %% sequenceness and null distributions
                                    
                                    for trial_i = 1:size(predictions.prob, 1)
                    
                                        disp([subject_ID ' - ' session_ID ' - trial ' num2str(trial_i)])
        
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
                            
                                        %% TDLM
                                        
                                        % X = state space
        
                                        % transition matrix from sequence (deleting
                                        % empty rows/columns, and ordering by
                                        % appearence in sequence, e.g. [2 5 4] is
                                        % encoded as [1 2 3], where [1==2, 2==5,
                                        % 3==4]
                                        T = zeros(nstates,nstates);
                                        for item_i = 2:length(sequence)
                                            T(item_i-1, item_i) = 1;
                                        end
        
                                        % number shuffle for null distribution
                                        nshuf = 24; % maximum permutations of 4 states
        
                                        [Msf, Msb] = TDLM_2ndGLM(X,T,nstates,nlags,nshuf);
                                
                                        Msf_alltrials{end+1} = Msf;
                                        Msb_alltrials{end+1} = Msb;
                                        Msz_alltrials{end+1} = Msf - Msb;
        
                                    end
                    
                                    %% average across permutations for each time lag
        
                                    % empirical sequenceness
                                    Msf_empirical = Msf_alltrials;
                                    Msb_empirical = Msb_alltrials;
                                    Msz_empirical = Msz_alltrials;
        
                                    for trial_i = 1:size(Msf_alltrials, 2)
                                        Msf_empirical{trial_i} = squeeze(Msf_empirical{trial_i}(1,1,:));
                                        Msb_empirical{trial_i} = squeeze(Msb_empirical{trial_i}(1,1,:));
                                        Msz_empirical{trial_i} = squeeze(Msz_empirical{trial_i}(1,1,:));
                                    end
        
                                    % empty matrix for permutations
                                    Msf = nan(size(Msf_alltrials{1}, 2), size(Msf_alltrials{1}, 3), size(Msf_alltrials, 2)); % permutation, time lag, trial
                                    Msb = nan(size(Msb_alltrials{1}, 2), size(Msb_alltrials{1}, 3), size(Msb_alltrials, 2)); % permutation, time lag, trial
                                    Msz = nan(size(Msz_alltrials{1}, 2), size(Msz_alltrials{1}, 3), size(Msz_alltrials, 2)); % permutation, time lag, trial
                                    
                                    % pool trials
                                    for trial_i = 1:size(Msf_alltrials, 2)
                                        
                                        Msf(:,:,trial_i) = squeeze(Msf_alltrials{trial_i});
                                        Msb(:,:,trial_i) = squeeze(Msb_alltrials{trial_i});
                                        Msz(:,:,trial_i) = squeeze(Msz_alltrials{trial_i});
                                    
                                    end
        
                                    % subset by sequence length (and remove first
                                    % permutation, that is, empirical)
                                    Msf_length3 = nan(5, size(Msf, 2), size(Msf, 3)); % permutations by lags by trials
                                    Msb_length3 = nan(5, size(Msf, 2), size(Msf, 3));
                                    Msz_length3 = nan(5, size(Msf, 2), size(Msf, 3));
        
                                    Msf_length4 = nan(23, size(Msf, 2), size(Msf, 3));
                                    Msb_length4 = nan(23, size(Msf, 2), size(Msf, 3));
                                    Msz_length4 = nan(23, size(Msf, 2), size(Msf, 3));
        
                                    Msf_all = nan(23, size(Msf, 2), size(Msf, 3));
                                    Msb_all = nan(23, size(Msf, 2), size(Msf, 3));
                                    Msz_all = nan(23, size(Msf, 2), size(Msf, 3));
        
                                    trials_length3 = [];
                                    trials_length4 = [];
                                    for trial_i = 1:size(Msf_alltrials, 2)
                                        if isnan(Msf(10,1,trial_i)) % length 3
                                            trials_length3 = [trials_length3 trial_i];
                                            Msf_length3(:,:,trial_i) = Msf(2:6,:,trial_i);
                                            Msb_length3(:,:,trial_i) = Msb(2:6,:,trial_i);
                                            Msz_length3(:,:,trial_i) = Msz(2:6,:,trial_i);
                                        else
                                            trials_length4 = [trials_length4 trial_i];
                                            Msf_length4(:,:,trial_i) = Msf(2:24,:,trial_i);
                                            Msb_length4(:,:,trial_i) = Msb(2:24,:,trial_i);
                                            Msz_length4(:,:,trial_i) = Msz(2:24,:,trial_i);
                                        end
                                        Msf_all(:,:,trial_i) = Msf(2:24,:,trial_i);
                                        Msb_all(:,:,trial_i) = Msb(2:24,:,trial_i);
                                        Msz_all(:,:,trial_i) = Msz(2:24,:,trial_i);
                                    end
                            
                                    Msf_length3 = Msf_length3(:,:,trials_length3);
                                    Msb_length3 = Msb_length3(:,:,trials_length3);
                                    Msz_length3 = Msz_length3(:,:,trials_length3);
                                    
                                    Msf_length4 = Msf_length4(:,:,trials_length4);
                                    Msb_length4 = Msb_length4(:,:,trials_length4);
                                    Msz_length4 = Msz_length4(:,:,trials_length4);
                                    
                                    % average permutations across trials
                                    Msf_all = median(Msf_all, 3, 'omitna');
                                    Msb_all = median(Msb_all, 3, 'omitna');
                                    Msz_all = median(Msz_all, 3, 'omitna');
        
                                    Msf_length3 = median(Msf_length3, 3);
                                    Msb_length3 = median(Msb_length3, 3);
                                    Msz_length3 = median(Msz_length3, 3);
                                    
                                    Msf_length4 = median(Msf_length4, 3);
                                    Msb_length4 = median(Msb_length4, 3);
                                    Msz_length4 = median(Msz_length4, 3);
        
                                    %% save outputs
        
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_empirical.mat'], 'Msf_empirical', '-v7.3');
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_empirical.mat'], 'Msb_empirical', '-v7.3');
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_empirical.mat'], 'Msz_empirical', '-v7.3');
                            
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_all.mat'], 'Msf_all', '-v7.3');
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_all.mat'], 'Msb_all', '-v7.3');
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_all.mat'], 'Msz_all', '-v7.3');
                
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_length3.mat'], 'Msf_length3', '-v7.3');
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_length3.mat'], 'Msb_length3', '-v7.3');
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_length3.mat'], 'Msz_length3', '-v7.3');
        
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_length4.mat'], 'Msf_length4', '-v7.3');
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_length4.mat'], 'Msb_length4', '-v7.3');
                                    save([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_length4.mat'], 'Msz_length4', '-v7.3');
        
                                end
                
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

%% pool across subjects

for folder_i = 1:length(folders)

    predictions_folder = folders{folder_i};

    try
        
        for state_i = 1:length(state_spaces)
        
            %% set input   
            for input_i = 1:length(inputs)
        
                if strcmp(inputs{input_i}, 'timepoints')
                    nlags = 500; % this variable can be shorter than 500, but no longer
                elseif strcmp(inputs{input_i}, 'bins')
                    nlags = 100;
                end
        
                %% loop across classifiers
                for class_i = 1:length(classifiers)
        
                    classifier = classifiers{class_i};
        
                    disp([state_spaces{state_i} ' ' inputs{input_i} ' ' classifier]);
        
                    if ~isfile([path_results '/group_results/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/data_individual_plots.mat'])
            
                        % empty structures for data
                        Msz_allsub = nan(nlags, 0);
                        Msf_allsub = nan(nlags, 0);
                        Msb_allsub = nan(nlags, 0);
                        Msz_l3_allsub = nan(nlags, 0);
                        Msf_l3_allsub = nan(nlags, 0);
                        Msb_l3_allsub = nan(nlags, 0);
                        Msz_l4_allsub = nan(nlags, 0);
                        Msf_l4_allsub = nan(nlags, 0);
                        Msb_l4_allsub = nan(nlags, 0);
            
                        null_Msz_thr_pos_allsub = [];
                        null_Msz_thr_neg_allsub = [];
                        null_Msf_thr_pos_allsub = [];
                        null_Msf_thr_neg_allsub = [];
                        null_Msb_thr_pos_allsub = [];
                        null_Msb_thr_neg_allsub = [];
            
                        null_Msz_thr_pos_l3_allsub = [];
                        null_Msz_thr_neg_l3_allsub = [];
                        null_Msf_thr_pos_l3_allsub = [];
                        null_Msf_thr_neg_l3_allsub = [];
                        null_Msb_thr_pos_l3_allsub = [];
                        null_Msb_thr_neg_l3_allsub = [];
            
                        null_Msz_thr_pos_l4_allsub = [];
                        null_Msz_thr_neg_l4_allsub = [];
                        null_Msf_thr_pos_l4_allsub = [];
                        null_Msf_thr_neg_l4_allsub = [];
                        null_Msb_thr_pos_l4_allsub = [];
                        null_Msb_thr_neg_l4_allsub = [];
            
                        % loop over trials
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
                        
                                disp([subject_ID ' - ' session_ID]); 
                        
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
            
                                %% null distribution computed at trial level, and averaged across trials for each time lag,
                                % then, for each permutation, take the
                                % minimum/maximum across time lags, to get two null
                                % distributions (minimum and maximum), and
                                % threshold at 95% to get the statistical bound
                                % see Schwartenbeck 2023
            
                                %% length all
                
                                % load sequenceness - empirical
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_empirical.mat']); % Msf_empirical
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_empirical.mat']); % Msb_empirical
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_empirical.mat']); % Msz_empirical
                
                                % average sequences across trials (for each time lag)
                                Msf_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
                                Msb_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
                                Msz_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
                
                                for trial_i = 1:size(Msf_tmp,2)
                            
                                    Msz_tmp(:,trial_i) = Msz_empirical{trial_i};
                                    Msf_tmp(:,trial_i) = Msf_empirical{trial_i};
                                    Msb_tmp(:,trial_i) = Msb_empirical{trial_i};
                
                                end
                
                                Msz = median(Msz_tmp, 2);
                                Msf = median(Msf_tmp, 2);
                                Msb = median(Msb_tmp, 2);
            
                                % pool across subjects
                                Msz_allsub = [Msz_allsub Msz];
                                Msf_allsub = [Msf_allsub Msf];
                                Msb_allsub = [Msb_allsub Msb];
                
                                %%% load null distribution - length all
            
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_all.mat']); % Msf_all
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_all.mat']); % Msb_all
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_all.mat']); % Msz_all
                
                                % omit na
                                Msf_all(isnan(Msf_all(:,1)),:) = [];
                                Msb_all(isnan(Msb_all(:,1)),:) = [];
                                Msz_all(isnan(Msz_all(:,1)),:) = [];
                                
                                % maximum (first column) and minimum (second column) across time lags for each permutation
                                null_Msf = [max(Msf_all, [],2) min(Msf_all, [],2)];
                                null_Msb = [max(Msb_all, [],2) min(Msb_all, [],2)];
                                null_Msz = [max(Msz_all, [],2) min(Msz_all, [],2)];
                
                                % threshold (P_fwe = 0.05)
                                null_Msz_thr_pos = percentile(null_Msz(:,1), 95);
                                null_Msz_thr_neg = percentile(null_Msz(:,2), 5);
                                null_Msf_thr_pos = percentile(null_Msf(:,1), 95);
                                null_Msf_thr_neg = percentile(null_Msf(:,2), 5);
                                null_Msb_thr_pos = percentile(null_Msb(:,1), 95);
                                null_Msb_thr_neg = percentile(null_Msb(:,2), 5);
            
                                null_Msz_thr_pos_allsub = [null_Msz_thr_pos_allsub null_Msz_thr_pos];
                                null_Msz_thr_neg_allsub = [null_Msz_thr_neg_allsub null_Msz_thr_neg];
                                null_Msf_thr_pos_allsub = [null_Msf_thr_pos_allsub null_Msf_thr_pos];
                                null_Msf_thr_neg_allsub = [null_Msf_thr_neg_allsub null_Msf_thr_neg];
                                null_Msb_thr_pos_allsub = [null_Msb_thr_pos_allsub null_Msb_thr_pos];
                                null_Msb_thr_neg_allsub = [null_Msb_thr_neg_allsub null_Msb_thr_neg];
            
                                %% length 3
            
                                % load sequenceness - empirical
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_empirical.mat']); % Msf_empirical
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_empirical.mat']); % Msb_empirical
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_empirical.mat']); % Msz_empirical
            
                                % find idx of trials of length 3 and 4
                                trials_length3 = [];
                                trials_length4 = [];
                                for trial_i = 1:size(data_task.trial, 2)
            
                                    if data_task.trialinfo(trial_i, 7) == 0
            
                                        trials_length3 = [trials_length3 trial_i];
            
                                    else
            
                                        trials_length4 = [trials_length4 trial_i];
            
                                    end
                                end
            
                                % subset trials - length 3
                                Msf_empirical = Msf_empirical(1, trials_length3);
                                Msb_empirical = Msb_empirical(1, trials_length3);
                                Msz_empirical = Msz_empirical(1, trials_length3);
                
                                % average sequences across trials (for each time lag)
                                Msz_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
                                Msf_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
                                Msb_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
                
                                for trial_i = 1:size(Msf_tmp,2)
                            
                                    Msz_tmp(:,trial_i) = Msz_empirical{trial_i};
                                    Msf_tmp(:,trial_i) = Msf_empirical{trial_i};
                                    Msb_tmp(:,trial_i) = Msb_empirical{trial_i};
                
                                end
        
                                Msz = median(Msz_tmp, 2);
                                Msf = median(Msf_tmp, 2);
                                Msb = median(Msb_tmp, 2);
            
                                % pool across subjects
                                Msz_l3_allsub = [Msz_l3_allsub Msz];
                                Msf_l3_allsub = [Msf_l3_allsub Msf];
                                Msb_l3_allsub = [Msb_l3_allsub Msb];
                
                                %%% load null distribution 
            
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_length3.mat']); % Msf_length3
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_length3.mat']); % Msb_length3
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_length3.mat']); % Msz_length3
                
                                % omit na
                                Msf_length3(isnan(Msf_length3(:,1)),:) = [];
                                Msb_length3(isnan(Msb_length3(:,1)),:) = [];
                                Msz_length3(isnan(Msz_length3(:,1)),:) = [];
                                
                                % maximum (first column) and minimum (second column) across time lags for each permutation
                                null_Msf = [max(Msf_length3, [],2) min(Msf_length3, [],2)];
                                null_Msb = [max(Msb_length3, [],2) min(Msb_length3, [],2)];
                                null_Msz = [max(Msz_length3, [],2) min(Msz_length3, [],2)];
                
                                % threshold (P_fwe = 0.05)
                                null_Msz_thr_pos = percentile(null_Msz(:,1), 95);
                                null_Msz_thr_neg = percentile(null_Msz(:,2), 5);
                                null_Msf_thr_pos = percentile(null_Msf(:,1), 95);
                                null_Msf_thr_neg = percentile(null_Msf(:,2), 5);
                                null_Msb_thr_pos = percentile(null_Msb(:,1), 95);
                                null_Msb_thr_neg = percentile(null_Msb(:,2), 5);
            
                                null_Msz_thr_pos_l3_allsub = [null_Msz_thr_pos_l3_allsub null_Msz_thr_pos];
                                null_Msz_thr_neg_l3_allsub = [null_Msz_thr_neg_l3_allsub null_Msz_thr_neg];
                                null_Msf_thr_pos_l3_allsub = [null_Msf_thr_pos_l3_allsub null_Msf_thr_pos];
                                null_Msf_thr_neg_l3_allsub = [null_Msf_thr_neg_l3_allsub null_Msf_thr_neg];
                                null_Msb_thr_pos_l3_allsub = [null_Msb_thr_pos_l3_allsub null_Msb_thr_pos];
                                null_Msb_thr_neg_l3_allsub = [null_Msb_thr_neg_l3_allsub null_Msb_thr_neg];
            
                                %% length 4
                
                                % load sequenceness - empirical
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_empirical.mat']); % Msf_empirical
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_empirical.mat']); % Msb_empirical
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_empirical.mat']); % Msz_empirical
            
                                % find idx of trials of length 3 and 4
                                trials_length3 = [];
                                trials_length4 = [];
                                for trial_i = 1:size(data_task.trial, 2)
            
                                    if data_task.trialinfo(trial_i, 7) == 0
            
                                        trials_length3 = [trials_length3 trial_i];
            
                                    else
            
                                        trials_length4 = [trials_length4 trial_i];
            
                                    end
                                end
            
                                % subset trials - length 3
                                Msf_empirical = Msf_empirical(1, trials_length4);
                                Msb_empirical = Msb_empirical(1, trials_length4);
                                Msz_empirical = Msz_empirical(1, trials_length4);
                
                                % average sequences across trials (for each time lag)
                                Msz_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
                                Msf_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
                                Msb_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
                
                                for trial_i = 1:size(Msf_tmp,2)
                            
                                    Msz_tmp(:,trial_i) = Msz_empirical{trial_i};
                                    Msf_tmp(:,trial_i) = Msf_empirical{trial_i};
                                    Msb_tmp(:,trial_i) = Msb_empirical{trial_i};
                
                                end
                
                                Msz = median(Msz_tmp, 2);
                                Msf = median(Msf_tmp, 2);
                                Msb = median(Msb_tmp, 2);
            
                                % pool across subjects
                                Msz_l4_allsub = [Msz_l4_allsub Msz];
                                Msf_l4_allsub = [Msf_l4_allsub Msf];
                                Msb_l4_allsub = [Msb_l4_allsub Msb];
                
                                %%% load null distribution 
            
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msf_length4.mat']); % Msf_length4
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msb_length4.mat']); % Msb_length4
                                load([path_results '/' subject_ID '/' session_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/Msz_length4.mat']); % Msz_length4
                
                                % omit na
                                Msf_length4(isnan(Msf_length4(:,1)),:) = [];
                                Msb_length4(isnan(Msb_length4(:,1)),:) = [];
                                Msz_length4(isnan(Msz_length4(:,1)),:) = [];
                                
                                % maximum (first column) and minimum (second column) across time lags for each permutation
                                null_Msf = [max(Msf_length4, [],2) min(Msf_length4, [],2)];
                                null_Msb = [max(Msb_length4, [],2) min(Msb_length4, [],2)];
                                null_Msz = [max(Msz_length4, [],2) min(Msz_length4, [],2)];
                
                                % threshold (P_fwe = 0.05)
                                null_Msz_thr_pos = percentile(null_Msz(:,1), 95);
                                null_Msz_thr_neg = percentile(null_Msz(:,2), 5);
                                null_Msf_thr_pos = percentile(null_Msf(:,1), 95);
                                null_Msf_thr_neg = percentile(null_Msf(:,2), 5);
                                null_Msb_thr_pos = percentile(null_Msb(:,1), 95);
                                null_Msb_thr_neg = percentile(null_Msb(:,2), 5);
            
                                null_Msz_thr_pos_l4_allsub = [null_Msz_thr_pos_l4_allsub null_Msz_thr_pos];
                                null_Msz_thr_neg_l4_allsub = [null_Msz_thr_neg_l4_allsub null_Msz_thr_neg];
                                null_Msf_thr_pos_l4_allsub = [null_Msf_thr_pos_l4_allsub null_Msf_thr_pos];
                                null_Msf_thr_neg_l4_allsub = [null_Msf_thr_neg_l4_allsub null_Msf_thr_neg];
                                null_Msb_thr_pos_l4_allsub = [null_Msb_thr_pos_l4_allsub null_Msb_thr_pos];
                                null_Msb_thr_neg_l4_allsub = [null_Msb_thr_neg_l4_allsub null_Msb_thr_neg];
            
                
                
                            end
                
                        end
            
                        if ~isfolder([path_results '/group_results/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures'])
                            mkdir([path_results '/group_results/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures']);
                        end
            
                        % save
                        save([path_results '/group_results/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/data_individual_plots.mat'],...
                            'Msz_allsub', 'Msf_allsub', 'Msb_allsub',...
                            'null_Msz_thr_pos_allsub', 'null_Msf_thr_pos_allsub', 'null_Msb_thr_pos_allsub',...
                            'null_Msz_thr_neg_allsub', 'null_Msf_thr_neg_allsub', 'null_Msb_thr_neg_allsub',...
                            'Msz_l3_allsub', 'Msf_l3_allsub', 'Msb_l3_allsub',...
                            'null_Msz_thr_pos_l3_allsub', 'null_Msf_thr_pos_l3_allsub', 'null_Msb_thr_pos_l3_allsub',...
                            'null_Msz_thr_neg_l3_allsub', 'null_Msf_thr_neg_l3_allsub', 'null_Msb_thr_neg_l3_allsub',...
                            'Msz_l4_allsub', 'Msf_l4_allsub', 'Msb_l4_allsub',...
                            'null_Msz_thr_pos_l4_allsub', 'null_Msf_thr_pos_l4_allsub', 'null_Msb_thr_pos_l4_allsub',...
                            'null_Msz_thr_neg_l4_allsub', 'null_Msf_thr_neg_l4_allsub', 'null_Msb_thr_neg_l4_allsub',...
                            '-v7.3');
        
                    end
        
                end
        
            end
        
        end
        
    catch
    end
    
end

disp('Analysis done');