%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = ['/path_to_local/results/preprocessing/' folder_preproc];
    path_results = ['/path_to_local/results/mvpa/' folder_decoding];
    classifiers = {'logreg', 'l2_e-1', 'l2_e-2', 'l2_e-3', 'l2_e-4', 'l1_e-1', 'l1_e-2', 'l1_e-3', 'l1_e-4'};
    inputs = {'bins'};
end

number_stimuli = 8;
nk              = 5;
nrepeat         = 10;

%% train classifiers task

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

        if ~isfolder([path_results '/' subject_ID '/' session_ID '/task_stim_balanced'])
            mkdir([path_results '/' subject_ID '/' session_ID '/task_stim_balanced']);
        end

        disp([subject_ID ' - ' session_ID]); 

        %% prepare task data

        % load
        load([path_inputs '/' subject_ID '/sess_0' num2str(session)...
                '/GERE_task_zthr' num2str(zthr) '_wholetrial_baseline1sec_postICA.mat']);

        data_task = data;

        clear data

        % subset correct trials
        subset_trials = data_task.trialinfo(:,3);
        subset_trials = find(subset_trials == 1);
        data_task.trial = data_task.trial(subset_trials);
        data_task.time = data_task.time(subset_trials);
        data_task.sampleinfo = data_task.sampleinfo(subset_trials, :);
        data_task.trialinfo = data_task.trialinfo(subset_trials, :);

        for input_i = 1:length(inputs)
    
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
            clabel = data_task.trialinfo(:,4:7);
            clabel = [clabel(:,1); clabel(:,2); clabel(:,3); clabel(:,4)];
            clabel = clabel(clabel>0);

            % loop across classifiers
            for class_i = 1:length(classifiers)
    
                classifier = classifiers{class_i};
    
                if ~isfolder([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier])
                    mkdir([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier]);
                end
    
                %% classify stim

                for stim_i = 1:number_stimuli

                    try

                        % decoding

                        if ~isfile([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier '/result_location' num2str(stim_i) '.mat'])
    
                            clabel_stim = clabel;

                            clabel_stim(clabel_stim == stim_i) = 100;
                            idx_location1 = find(clabel_stim == 1); % location 1
                            idx_location2 = find(clabel_stim == 2); % location 2...
                            idx_location3 = find(clabel_stim == 3);
                            idx_location4 = find(clabel_stim == 4);
                            idx_location5 = find(clabel_stim == 5);
                            idx_location6 = find(clabel_stim == 6);
                            idx_location7 = find(clabel_stim == 7);
                            idx_location8 = find(clabel_stim == 8);
                            clabel_stim(clabel_stim == 100) = stim_i;
    
                            % clabel_binary contains 1 for positive
                            % example, 2 for all negative examples (other
                            % locations)
                            clabel_binary = clabel_stim;
                            clabel_binary(clabel_binary == stim_i) = 100;
                            clabel_binary(clabel_binary < 100) = 2;
                            clabel_binary(clabel_binary == 100) = 1;
    
                            % indices of positive examples
                            idx_label1 = find(clabel_stim == stim_i); 
    
                            % balance number classes for negative examples
                            idx_label2_balanced = nan(size(idx_label1));
                            stimxlocation = round(length(idx_label2_balanced)/7);
    
                            % assign equal number of random indices per
                            % location
                            idx = [];
                            for loc_i = 0:7
    
                                try
                                    tmp = randperm(length(eval(['idx_location' num2str(loc_i+1)])), stimxlocation);
                                    idx = [idx eval(['idx_location' num2str(loc_i+1) '([' num2str(tmp) '],1)'])'];
    
                                catch
                                end
    
                            end
    
                            % complete empty indices
                            for idx_i = 1:length(idx_label2_balanced)
    
                                if idx_i > length(idx)
                                    no_repetition = 0;
                                    while ~no_repetition
                                        idx_new = randi(length(clabel));
                                        if ismember(idx_new, idx)
                                            idx = [idx idx_new];
                                            no_repetition = 1;
                                        end
                                    end
                                end
    
                            end
    
                            idx_label2_balanced = idx';
                                         
                            % pool indices
                            idx_all = [idx_label1; idx_label2_balanced];
    
                            % create data and clabel
                            X_undersample = X(idx_all,:,:);
                            clabel_undersample = clabel_binary(idx_all,1);
                            
                            % classification across time
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
                            elseif strcmp(classifier, 'l2_vector')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l2';
                                cfg.hyperparameter.lambda   = [10 1 0.1 0.01 0.001];
                            elseif strcmp(classifier, 'l1_vector')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l1';
                                cfg.hyperparameter.lambda   = [10 1 0.1 0.01 0.001];
                            elseif strcmp(classifier, 'l2_e-1')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l2';
                                cfg.hyperparameter.lambda      = 0.1;
                            elseif strcmp(classifier, 'l2_e-2')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l2';
                                cfg.hyperparameter.lambda      = 0.01;
                            elseif strcmp(classifier, 'l1_e-3')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l1';
                                cfg.hyperparameter.lambda      = 0.001;
                            elseif strcmp(classifier, 'l1_e-4')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l1';
                                cfg.hyperparameter.lambda      = 0.0001;
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
                            cfg.metric                  = {'auc' 'accuracy'};
                            cfg.k                       = nk;
                            cfg.repeat                  = nrepeat;
                            cfg.preprocess              = {'zscore'}; % nested preprocessing
                            [~, result, ~]   = mv_classify_across_time(cfg, X_undersample, clabel_undersample);

                            result.idx_undersample = idx_all;
        
                            % save
                            save([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier...
                                '/result_location' num2str(stim_i) '.mat'], 'result');

                        end
        
                        % select result by auc
                        result = mv_select_result(result, 'auc');
    
                        % plot auc
                        mv_plot_result(result, 1:length(result.perf));
    
                        % save
                        saveas(gcf, [path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier '/figure_location' num2str(stim_i) '.fig'])
                        fig = openfig([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier '/figure_location' num2str(stim_i) '.fig']);
                        saveas(fig, [path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier '/figure_location' num2str(stim_i) '.jpeg']);
                        close all

                        % time-generalization

                        if ~isfile([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier...
                                '/timextime_location' num2str(stim_i) '.mat'])

                            % load result.idx_undersample
                            load([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier...
                                '/result_location' num2str(stim_i) '.mat']); % result

                            clabel_stim = clabel;

                            clabel_stim(clabel_stim == stim_i) = 100;
                            idx_location1 = find(clabel_stim == 1); % location 1
                            idx_location2 = find(clabel_stim == 2); % location 2...
                            idx_location3 = find(clabel_stim == 3);
                            idx_location4 = find(clabel_stim == 4);
                            idx_location5 = find(clabel_stim == 5);
                            idx_location6 = find(clabel_stim == 6);
                            idx_location7 = find(clabel_stim == 7);
                            idx_location8 = find(clabel_stim == 8);
                            clabel_stim(clabel_stim == 100) = stim_i;
    
                            % clabel_binary contains 1 for positive
                            % example, 2 for all negative examples (other
                            % locations)
                            clabel_binary = clabel_stim;
                            clabel_binary(clabel_binary == stim_i) = 100;
                            clabel_binary(clabel_binary < 100) = 2;
                            clabel_binary(clabel_binary == 100) = 1;
    
                            % indices
                            idx_all = result.idx_undersample;
    
                            % create data and clabel
                            X_undersample = X(idx_all,:,:);
                            clabel_undersample = clabel_binary(idx_all,1);
                            
                            % classification across time
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
                            elseif strcmp(classifier, 'l2_vector')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l2';
                                cfg.hyperparameter.lambda   = [10 1 0.1 0.01 0.001];
                            elseif strcmp(classifier, 'l1_vector')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l1';
                                cfg.hyperparameter.lambda   = [10 1 0.1 0.01 0.001];
                            elseif strcmp(classifier, 'l2_e-1')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l2';
                                cfg.hyperparameter.lambda      = 0.1;
                            elseif strcmp(classifier, 'l2_e-2')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l2';
                                cfg.hyperparameter.lambda      = 0.01;
                            elseif strcmp(classifier, 'l1_e-3')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l1';
                                cfg.hyperparameter.lambda      = 0.001;
                            elseif strcmp(classifier, 'l1_e-4')
                                cfg.classifier              = 'logreg';
                                cfg.hyperparameter.reg      = 'l1';
                                cfg.hyperparameter.lambda      = 0.0001;
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
                            cfg.metric                  = {'auc' 'accuracy'};
                            cfg.k                       = nk;
                            cfg.repeat                  = nrepeat;
                            cfg.preprocess              = {'zscore'}; % nested preprocessing
                            [~, timextime, ~]           = mv_classify_timextime(cfg, X_undersample, clabel_undersample);
        
                            % save
                            save([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier...
                                '/timextime_location' num2str(stim_i) '.mat'], 'timextime');
                            
                            % plot time generalization - auc
                            close all
                            imagesc(1:length(result.perf), 1:length(result.perf), timextime.perf{1});
                            set(gca, 'YDir', 'normal');
                            colorbar; grid on;
                            ylabel('Training time'), xlabel('Test time');
                            saveas(gcf, [path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier '/figure_timextime_location' num2str(stim_i) '.fig'])
                            fig = openfig([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier '/figure_timextime_location' num2str(stim_i) '.fig']);
                            saveas(fig, [path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier '/figure_timextime_location' num2str(stim_i) '.jpeg']);
                            close all
    
                        end

                    catch
                    end

                end
                
            end

        end

    end

end

disp('Analysis done');