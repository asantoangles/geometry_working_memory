%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = ['/path_to_local/results/preprocessing/' folder_preproc];
    path_results = ['/path_to_local/results/mvpa/' folder_decoding];
end

number_stimuli = 8;

%% compute decoded state space (subject-specific) for the delay period
% with the classifier trained at the time with the highest performance averaged across locations

for sub_i = 1:length(subjects)

    subject = subjects(sub_i);

    % loop over sessions
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

        disp([subject_ID ' - ' session_ID]); 
    
        % loop across classifiers
        for class_i = 1:length(classifiers)

            classifier = classifiers{class_i};

            for input_i = 1:length(inputs)

                if ~isfile([path_results '/' subject_ID '/' session_ID '/task_stim2delay_balanced_trial_predictions_averaged/' inputs{input_i} '/' classifier...
                                '/predictions_singletrial.mat'])

                    if strcmp(inputs{input_i}, 'timepoints')
                        time_delay = 4000;
                    elseif strcmp(inputs{input_i}, 'bins')
                        time_delay = 791;
                    end
        
                    if ~isfolder([path_results '/' subject_ID '/' session_ID '/task_stim2delay_balanced_trial_predictions_averaged/' inputs{input_i} '/' classifier])
                        mkdir([path_results '/' subject_ID '/' session_ID '/task_stim2delay_balanced_trial_predictions_averaged/' inputs{input_i} '/' classifier]);
                    end
    
                    % average best time
                    best_times = [];
                    for stim_i = 1:number_stimuli
    
                        load([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier...
                                '/result_location' num2str(stim_i) '.mat']); % result
        
                        result = mv_select_result(result, 'auc');
            
                        [val, best_time] = max(result.perf, [], 1);
    
                        best_times = [best_times best_time];
    
                    end
    
                    best_time = round(median(best_times)); clear best_times
    
                    for stim_i = 1:number_stimuli
    
                        disp([classifier ' - ' inputs{input_i} ' - stimulus ' num2str(stim_i)])
    
                        % load
                        load([path_inputs '/' subject_ID '/sess_0' num2str(session)...
                                '/GERE_task_zthr' num2str(zthr) '_wholetrial_baseline1sec_postICA.mat']); % data_task
    
                        data_task = data;
    
                        clear data
    
                        % subset correct trials
                        subset_trials = data_task.trialinfo(:,3);
                        subset_trials = find(subset_trials == 1);
                        data_task.trial = data_task.trial(subset_trials);
                        data_task.time = data_task.time(subset_trials);
                        data_task.sampleinfo = data_task.sampleinfo(subset_trials, :);
                        data_task.trialinfo = data_task.trialinfo(subset_trials, :);
            
                        %% prepare task data - delay
                        
                        % X_test is variable with trial x channels x timepoints in task
                        % time_test is vector with the timepoints in task
                    
                        if strcmp(inputs{input_i}, 'bins')
        
                            % create bins
                            timepoints = 4000; % delay 4 seconds
                            bin_length = 50;
                            bin_step = 5;
                            bin_number = round((timepoints) / bin_step) - (bin_step*2)+1;
                        
                            for bin_i = 1:bin_number
                                if bin_i == 1
                                    bins = 1:bin_length;
                                else
                                    bins = [bins; bins(end,:)+5];
                                end
                            end
        
                            % subset data_task
                            X = nan(size(data_task.trial, 2), size(data_task.label, 1), bin_number);
                            for trial_i = 1:size(data_task.trialinfo, 1)
                                all_timepoints = (size(data_task.trial{trial_i}, 2) - 3999):size(data_task.trial{trial_i}, 2);
                                for bin_i = 1:bin_number
                                    X(trial_i, :, bin_i) = mean(data_task.trial{trial_i}(:, all_timepoints(bins(bin_i,:))), 2)';
                                end
                            end
                    
                            X_test = X;
                        
                            % time
                            time_test = squeeze(mean(bins'));
        
                            clear X
        
                        elseif strcmp(inputs{input_i}, 'timepoints')
                        
                            timepoints = 4000; % delay 4 seconds
        
                            % subset data_task with only delay period
                            for trial_i = 1:size(data_task.trial, 2)
                                data_task.time{trial_i} = data_task.time{trial_i}(1, (end-timepoints+1):end);
                                data_task.trial{trial_i} = data_task.trial{trial_i}(:, (end-timepoints+1):end);
                            end
        
                            data_task.sampleinfo(:,1) = data_task.sampleinfo(:,2) - 4000;
        
                            % reset data_task
                            X = nan(size(data_task.trial, 2), size(data_task.label, 1), timepoints);
                            for trial_i = 1:size(data_task.trialinfo, 1)
                                X(trial_i, :, :) = data_task.trial{trial_i};
                            end
        
                            X_test = X;
        
                            time_test = 1:timepoints;
        
                            clear X
        
                        end                        
                        
                        %% prepare data - task stimulation
        
                        load([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier...
                                '/result_location' num2str(stim_i) '.mat']); % result
        
                        result = mv_select_result(result, 'auc');
    
                        % load
                        load([path_inputs '/' subject_ID '/sess_0' num2str(session)...
                                '/GERE_task_zthr' num2str(zthr) '_wholetrial_baseline1sec_postICA.mat']); % data_task
    
                        data_task = data;
    
                        clear data
    
                        % subset correct trials
                        subset_trials = data_task.trialinfo(:,3);
                        subset_trials = find(subset_trials == 1);
                        data_task.trial = data_task.trial(subset_trials);
                        data_task.time = data_task.time(subset_trials);
                        data_task.sampleinfo = data_task.sampleinfo(subset_trials, :);
                        data_task.trialinfo = data_task.trialinfo(subset_trials, :);
                        
                        % X is variable with trial x channels x timepoints in task
                        % time_test is vector with the timepoints in task
                    
                        if strcmp(inputs{input_i}, 'bins')
                
                            % create bins
                            timepoints = 1700; % 100 ms baseline + (300 ms * 4 stimuli = 1200) + (100 ms * 3 isi = 300) + 100 ms delay
                            bin_length = 50;
                            bin_step = 5;
                            bin_number = round((timepoints) / bin_step) - (bin_step*2)+1;
                        
                            for bin_i = 1:bin_number
                                if bin_i == 1
                                    bins = 1:bin_length;
                                else
                                    bins = [bins; bins(end,:)+5];
                                end
                            end
                
                            % subset data_task
                            X = nan(size(data_task.trial, 2), size(data_task.label, 1), bin_number);
                            for trial_i = 1:size(data_task.trialinfo, 1)
                                all_timepoints = 901:2600;
                                for bin_i = 1:bin_number
                                    X(trial_i, :, bin_i) = mean(data_task.trial{trial_i}(:, all_timepoints(bins(bin_i,:))), 2)';
                                end
                            end
                            
                            % time
                            time = squeeze(mean(bins'));
        
                            % format X to accomodate order sequence
                            X_oldformat = X;
                            X = nan(size(X,1), size(X, 2), 91, 4); % third dimension: time comprising 100 ms isi + (300 ms stim * 4) + (100 ms isi * 3) + 100 ms delay
                                                                    % forth dimension: order sequence
            
                            for trial_i = 1:size(X, 1)
            
                                X(trial_i, :,:,1) = X_oldformat(trial_i, :, 1:91);          % 1:500 in timepoints
                                X(trial_i, :,:,2) = X_oldformat(trial_i, :, 81:171);        % 401:900 in timepoints
                                X(trial_i, :,:,3) = X_oldformat(trial_i, :, 161:251);       % 801:1300 in timepoints
            
                                if data_task.trialinfo(trial_i, 7) > 0
                                    X(trial_i, :,:,4) = X_oldformat(trial_i, :, 241:331); % 1201:1700 in timepoints
                                end
            
                            end
                        
                        end           
        
                        % subset X for trials of length 4
        
                        X_order4 = nan(0, size(X, 2), size(X, 3), 4);
        
                        for trial_i = 1:size(X, 1)
        
                            if data_task.trialinfo(trial_i, 7) > 0
        
                                X_order4(end+1,:,:,:) = X(trial_i,:,:,:);
        
                            end
        
                        end
        
                        % reformat X and clabel
        
                        X = [squeeze(X(:,:,:,1)); squeeze(X(:,:,:,2)); squeeze(X(:,:,:,3)); squeeze(X_order4(:,:,:,4))];
        
                        % set labels
                        clabel_train = data_task.trialinfo(:,4:7);
                        clabel_train = [clabel_train(:,1); clabel_train(:,2); clabel_train(:,3); clabel_train(:,4)];
                        clabel_train = clabel_train(clabel_train>0);
    
                        clabel_train(clabel_train == stim_i) = 100;
                        clabel_train(clabel_train < 100) = 2;
                        clabel_train(clabel_train == 100) = 1;
        
                        % subset best time
                        X_besttime_1tp = X(:,:,best_time);
                    
                        % set train data (repeating best time classifier)
                        X_train = nan(size(X_besttime_1tp, 1), size(X_besttime_1tp, 2), length(time_test));
                        for time_i = 1:size(X_train, 3)
                            X_train(:,:,time_i) = X_besttime_1tp;
                        end
    
                        clabel_test = nan(size(X_test, 1), 1);
        
                        for trial_i = 1:length(clabel_test)
        
                            if data_task.trialinfo(trial_i, 4) == stim_i || ...
                               data_task.trialinfo(trial_i, 5) == stim_i || ...
                               data_task.trialinfo(trial_i, 6) == stim_i || ...
                               data_task.trialinfo(trial_i, 7) == stim_i
        
                                clabel_test(trial_i) = 1;
        
                            else
        
                                clabel_test(trial_i) = 2;
        
                            end
        
                        end
                    
                                    
                        %% compute weights
                        cfg                         = [];
                        if strcmp(classifier, 'l1_e-1')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l1';
                            cfg.hyperparameter.lambda      = 0.1;
                        elseif strcmp(classifier, 'l2_e-3')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l2';
                            cfg.hyperparameter.lambda      = 0.001;
                        elseif strcmp(classifier, 'l2_e-4')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l2';
                            cfg.hyperparameter.lambda      = 0.0001;
                        elseif strcmp(classifier, 'l1_e-3')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l1';
                            cfg.hyperparameter.lambda      = 0.001;
                        elseif strcmp(classifier, 'l1_e-4')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l1';
                            cfg.hyperparameter.lambda      = 0.0001;
                        elseif strcmp(classifier, 'l1_e-2')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l1';
                            cfg.hyperparameter.lambda      = 0.01;
                        elseif strcmp(classifier, 'l2')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l2';
                        elseif strcmp(classifier, 'l1')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l1';
                        elseif strcmp(classifier, 'l2_e-1')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l2';
                            cfg.hyperparameter.lambda      = 0.1;
                        elseif strcmp(classifier, 'l2_e-2')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l2';
                            cfg.hyperparameter.lambda      = 0.01;
                        elseif strcmp(classifier, 'l1_e-5')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l1';
                            cfg.hyperparameter.lambda      = 0.00001;
                        elseif strcmp(classifier, 'l1_e-6')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l1';
                            cfg.hyperparameter.lambda      = 0.000001;
                        elseif strcmp(classifier, 'l1_e-7')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l1';
                            cfg.hyperparameter.lambda      = 0.0000001;
                        elseif strcmp(classifier, 'l1_adaptive')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l1_adaptive';
                            cfg.hyperparameter.lambda      = 0.000001; % lambda for prob < 0.5, lambda/100 for prob > 0.5
                        elseif strcmp(classifier, 'l2_adaptive')
                            cfg.classifier              = 'logreg';
                            cfg.hyperparameter.reg      = 'l2_adaptive';
                            cfg.hyperparameter.lambda      = 0.01; % lambda for prob < 0.5, lambda/100 for prob > 0.5
                        else
                            cfg.classifier              = classifier;
                        end
                        cfg.metric                  = {'dval'};
                        cfg.preprocess              = {'zscore'}; % nested preprocessing
                        cfg.save                    = {'model_param'};
                        [~, result, ~]   = mv_classify(cfg, X_train(:,:,1), clabel_train, X_test(1,:,1), clabel_test(1));
    
                        cf = result.model_param{1};

                        save([path_results '/' subject_ID '/' session_ID '/task_stim2delay_balanced_trial_predictions_averaged/' inputs{input_i} '/' classifier...
                                '/weights_location' num2str(stim_i) '.mat'], 'cf', '-v7.3');
    
                        %% preprocessing X_test
    
                        pparam = mv_get_preprocess_param('zscore');
                        [~, X_preproc] = mv_preprocess_zscore(pparam, X_test);
    
                        %% compute predictions
    
                        dval_all = nan(size(X_preproc, 1), size(X_preproc, 3));
    
                        for trial_i = 1:size(dval_all,1)
                            
                            %% single trial predictions
                        
                            for time_i = 1:size(X_preproc, 3)
            
                                X_trial = squeeze(X_preproc(trial_i, :, time_i));
    
                                if strcmp(classifier, 'svm')
    
                                    [~,dval] = test_svm(cf, X_trial);
    
                                else
    
                                    [~,dval] = test_logreg(cf, X_trial);
    
                                end
    
                                dval_all(trial_i, time_i) = dval;
    
                            end
                                            
                        end
                                    
                        predictions_stim = dval_all;
    
                        filename_output = [path_results '/' subject_ID '/' session_ID '/task_stim2delay_balanced_trial_predictions_averaged/' inputs{input_i} '/' classifier...
                                '/predictions_trial_predictions_dval_stim' num2str(stim_i) '.mat'];
        
                        save(filename_output, 'predictions_stim', '-v7.3');
            
                    end
    
                    % merge predictions
                    predictions = [];
                    predictions.info = {'dimensions: trial / time (delay) / stimuli'};
                    predictions.dval = nan(size(X_test, 1), size(X_test, 3), number_stimuli);
    
                    for stim_i = 1:number_stimuli
    
                        filename_output = [path_results '/' subject_ID '/' session_ID '/task_stim2delay_balanced_trial_predictions_averaged/' inputs{input_i} '/' classifier...
                                '/predictions_trial_predictions_dval_stim' num2str(stim_i) '.mat'];
    
                        load(filename_output)
    
                        predictions.dval(:,:,stim_i) = predictions_stim;
    
                    end
    
                    save([path_results '/' subject_ID '/' session_ID '/task_stim2delay_balanced_trial_predictions_averaged/' inputs{input_i} '/' classifier...
                                '/predictions_singletrial.mat'], 'predictions', '-v7.3');
    
                    system(['rm ' path_results '/' subject_ID '/' session_ID '/task_stim2delay_balanced_trial_predictions_averaged/' inputs{input_i} '/' classifier...
                                '/predictions_trial_predictions_dval_stim*']);

                end

            end
            
        end
        
    end

end


disp('analysis done');
