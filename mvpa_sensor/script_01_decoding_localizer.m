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

%% train classifiers localizer

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

        if ~isfolder([path_results '/' subject_ID '/' session_ID])
            mkdir([path_results '/' subject_ID '/' session_ID]);
        end

        disp([subject_ID ' - ' session_ID]); 

        if isfile([path_inputs '/' subject_ID '/sess_0' num2str(session)...
                '/GERE_localizer_validtrialszthr' num2str(zthr) '_postICA.mat'])

            % load data localizer
            load([path_inputs '/' subject_ID '/sess_0' num2str(session)...
                    '/GERE_localizer_validtrialszthr' num2str(zthr) '_postICA.mat']); % data
    
            data_localizer = data;
    
            clear data
    
            %% set input   
            for input_i = 1:length(inputs)
                
                if strcmp(inputs{input_i}, 'bins')
                
                    % X is variable with trial x channels x timepoints 
                    % (average in windows of 50 ms, in steps of 5 ms)
                
                    % create bins
                    timepoints = size(data_localizer.time{1}, 2);
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
                
                    % subset data
                    X = nan(size(data_localizer.trial, 2), size(data_localizer.label, 1), bin_number);
                    for trial_i = 1:size(data_localizer.trialinfo, 1)
                        for bin_i = 1:bin_number
                            for sens_i = 1:size(data_localizer.label, 1)
                                X(trial_i, sens_i, bin_i) = mean(data_localizer.trial{trial_i}(sens_i, bins(bin_i,:)));
                            end
                        end
                    end
                
                    % time
                    time = mean(bins');
    
                    X_all = X;
                                
                end
    
                %% loop across classifiers
                for class_i = 1:length(classifiers)
        
                    classifier = classifiers{class_i};
    
                    if ~isfolder([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier])
                        mkdir([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier]);
                    end
                            
                    % loop over stimuli
                    for stim_i = 1:number_stimuli
    
                        X = X_all;
                        clabel_binary = data_localizer.trialinfo(:,3);
                        clabel_binary(clabel_binary == stim_i) = 100;
                        clabel_binary(clabel_binary < 100) = 2;
                        clabel_binary(clabel_binary == 100) = 1;
                        
                        %% balanced number of samples
    
                        disp([inputs{input_i} ' ' classifier ' balanced - stimulus ' num2str(stim_i)]);
    
                        if ~isfile([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier...
                                '/result_balanced_location' num2str(stim_i) '.mat'])
    
                            % clabel contains the spatial location [1-8] in each trial
                            clabel = data_localizer.trialinfo(:,3);
                            clabel(clabel == stim_i) = 100;
                            idx_location1 = find(clabel == 1); % location 1
                            idx_location2 = find(clabel == 2); % location 2...
                            idx_location3 = find(clabel == 3);
                            idx_location4 = find(clabel == 4);
                            idx_location5 = find(clabel == 5);
                            idx_location6 = find(clabel == 6);
                            idx_location7 = find(clabel == 7);
                            idx_location8 = find(clabel == 8);
                            clabel(clabel == 100) = stim_i;
    
                            % clabel_binary contains 1 for positive
                            % example, 2 for all negative examples (other
                            % locations)
                            clabel_binary = data_localizer.trialinfo(:,3);
                            clabel_binary(clabel_binary == stim_i) = 100;
                            clabel_binary(clabel_binary < 100) = 2;
                            clabel_binary(clabel_binary == 100) = 1;
    
                            % indices of positive examples
                            idx_label1 = find(clabel == stim_i); 
    
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
                            X_undersample = X_all(idx_all,:,:);
                            clabel_undersample = clabel_binary(idx_all,1);
    
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
                            cfg.metric                  = {'auc' 'accuracy' 'dval' 'tval'};
                            cfg.k                       = nk;
                            cfg.repeat                  = nrepeat;
                            cfg.preprocess              = {'zscore'}; % nested preprocessing
                            [~, result, ~]   = mv_classify_across_time(cfg, X_undersample, clabel_undersample);
    
                            result.idx_undersample = idx_all;
    
                            % save
                            save([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier...
                                '/result_balanced_location' num2str(stim_i) '.mat'], 'result');
    
                            % select result by auc
                            result = mv_select_result(result, 'auc');
    
                            % plot auc
                            mv_plot_result(result, time);
    
                            % save
                            saveas(gcf, [path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier '/figure_balanced_location' num2str(stim_i) '.fig'])
                            fig = openfig([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier '/figure_balanced_location' num2str(stim_i) '.fig']);
                            saveas(fig, [path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier '/figure_balanced_location' num2str(stim_i) '.jpeg']);
                            close all

                        end

                        if ~isfile([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier...
                                '/figure_timextime_balanced_location' num2str(stim_i) '.mat'])

                            % load result.idx_undersample
                            load([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier...
                                '/result_balanced_location' num2str(stim_i) '.mat']); % result

                            % clabel contains the spatial location [1-8] in each trial
                            clabel = data_localizer.trialinfo(:,3);
                            clabel(clabel == stim_i) = 100;
                            idx_location1 = find(clabel == 1); % location 1
                            idx_location2 = find(clabel == 2); % location 2...
                            idx_location3 = find(clabel == 3);
                            idx_location4 = find(clabel == 4);
                            idx_location5 = find(clabel == 5);
                            idx_location6 = find(clabel == 6);
                            idx_location7 = find(clabel == 7);
                            idx_location8 = find(clabel == 8);
                            clabel(clabel == 100) = stim_i;
    
                            % clabel_binary contains 1 for positive
                            % example, 2 for all negative examples (other
                            % locations)
                            clabel_binary = data_localizer.trialinfo(:,3);
                            clabel_binary(clabel_binary == stim_i) = 100;
                            clabel_binary(clabel_binary < 100) = 2;
                            clabel_binary(clabel_binary == 100) = 1;
                                                 
                            % indices
                            idx_all = result.idx_undersample;
    
                            % create data and clabel
                            X_undersample = X_all(idx_all,:,:);
                            clabel_undersample = clabel_binary(idx_all,1);

                            % time generalization
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
                            cfg.metric                  = {'auc' 'accuracy' 'dval' 'tval'};
                            cfg.k                       = nk;
                            cfg.repeat                  = nrepeat;
                            cfg.preprocess              = {'zscore'}; % nested preprocessing
                            [~, timextime, ~]           = mv_classify_timextime(cfg, X_undersample, clabel_undersample);
    
                            % save
                            save([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier...
                                '/timextime_balanced_location' num2str(stim_i) '.mat'], 'timextime');
                            
                            % plot time generalization - auc
                            close all
                            imagesc(time, time, timextime.perf{1});
                            set(gca, 'YDir', 'normal');
                            colorbar; grid on;
                            ylabel('Training time'), xlabel('Test time');
                            saveas(gcf, [path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier '/figure_timextime_balanced_location' num2str(stim_i) '.fig'])
                            fig = openfig([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier '/figure_timextime_balanced_location' num2str(stim_i) '.fig']);
                            saveas(fig, [path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier '/figure_timextime_balanced_location' num2str(stim_i) '.jpeg']);
                            close all
    
                        end
    
                    end
                    
                end
            
            end

        end
        
    end

end

disp('Analysis done')