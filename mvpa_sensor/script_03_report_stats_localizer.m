%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = ['/path_to_local/results/preprocessing/' folder_preproc];
    path_results = ['/path_to_local/results/mvpa/' folder_decoding];
    classifiers = {'logreg', 'l2_e-1', 'l2_e-2', 'l2_e-3', 'l2_e-4', 'l1_e-1', 'l1_e-2', 'l1_e-3', 'l1_e-4'};
    inputs = {'bins'};
end

number_stimuli = 8;

%% report stats localizer

for input_i = 1:length(inputs)

    if strcmp(inputs{input_i}, 'bins')
        time = 25:5:24+(71)*5;
    end
    
    % loop across classifiers
    for class_i = 1:length(classifiers)

        classifier = classifiers{class_i};
        
        %% pool results across subjects
        results_across_subjects = cell(0);

        for sub_i = 1:length(subjects)
        
            subject = subjects(sub_i);
        
            for ses_i = 1:length(sessions)
        
                session = sessions(ses_i);
        
                %% load results
        
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
        
                results_across_stimuli = cell(1,number_stimuli);
                        
                % loop over stimuli
                for stim_i = 1:number_stimuli
                    
                    %% balanced number of samples

                    % result 
                    load([path_results '/' subject_ID '/' session_ID '/localizer_twoclass/' inputs{input_i} '/' classifier...
                        '/result_balanced_location' num2str(stim_i) '.mat']); % result

                    if isnan(result.perf{1}(1))
                        disp(['nan stim' num2str(stim_i) ' ' subject_ID ' ' session_ID]);
                    end

                    results_across_stimuli{stim_i} = result;

                end

                % plot average across locations
                result_average = mv_combine_results(results_across_stimuli, 'average');
                result_average = mv_select_result(result_average, 'auc');

                results_across_subjects{end+1} = result_average;
                    
            end
        
        end

        result_average = mv_combine_results(results_across_subjects, 'average');
    
        %% reports

        %%% peak

        [~,idx_group] = max(result_average.perf{1});

        peak_vector = [];

        for sub_i = 1:length(subjects)

            [~,idx] = max(results_across_subjects{sub_i}.perf);

            peak_vector = [peak_vector idx];

        end
        
        disp(['peak: ' num2str(time(idx_group)) ' ± ' num2str(time(round(std(peak_vector))))]);

        %%% auc

        auc_vector = [];

        for sub_i = 1:length(subjects)

            [value,~] = max(results_across_subjects{sub_i}.perf);

            auc_vector = [auc_vector value];

        end

        disp(['auc: ' num2str(mean(auc_vector)) ' ± ' num2str(std(auc_vector))]);
        
    end
    
end



