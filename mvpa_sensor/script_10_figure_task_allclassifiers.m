%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = ['/path_to_local/results/preprocessing/' folder_preproc];
    path_results = ['/path_to_local/results/mvpa/' folder_decoding];
end

number_stimuli = 8;
classifiers = {'logreg', 'l2_e-1', 'l2_e-2', 'l2_e-3', 'l2_e-4', 'l1_e-1', 'l1_e-2', 'l1_e-3', 'l1_e-4'};
inputs = {'bins'};

%% figure decoding task all classifiers

for input_i = 1:length(inputs)

    if strcmp(inputs{input_i}, 'bins')
        time = 25:5:24+(91)*5;
    end

    performance_allclassifiers = nan(length(classifiers), length(time));
        
    % pool performance across classifiers
    for class_i = 1:length(classifiers)

        classifier = classifiers{class_i};
        
        %% classification across time

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
                    load([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier...
                        '/result_location' num2str(stim_i) '.mat']); % result

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

        performance_allclassifiers(class_i,:) = result_average.perf{1}';

    end

    % plot
    figure
    for class_i = 1:length(classifiers)

        if strcmp(classifiers{class_i},'logreg')

            color = [0 0 0]; % black
            trace = '-';

        elseif strcmp(classifiers{class_i},'l2_e-1')

            color = [0.7 0 0]; % red
            trace = '-';

        elseif strcmp(classifiers{class_i},'l2_e-2')

            color = [0 0.7 0]; % green
            trace = '-';

        elseif strcmp(classifiers{class_i},'l2_e-3')

            color = [0 0 0.7]; % blue
            trace = '-';

        elseif strcmp(classifiers{class_i},'l2_e-4')

            color = [1 0.5 0]; % orange
            trace = '-';

        elseif strcmp(classifiers{class_i},'l1_e-1')

            color = [0.7 0 0]; % red
            trace = '--';

        elseif strcmp(classifiers{class_i},'l1_e-2')

            color = [0 0.7 0]; % green
            trace = '--';

        elseif strcmp(classifiers{class_i},'l1_e-3')

            color = [0 0 0.7]; % blue
            trace = '--';

        elseif strcmp(classifiers{class_i},'l1_e-4')

            color = [1 0.5 0]; % orange
            trace = '--';

        end

        % subset 1 to 300 ms since stimulus onset
        subset_auc = 21:71;
        plot(25:5:275, performance_allclassifiers(class_i,subset_auc), trace, 'Color', color, 'LineWidth', 2); hold on

    end

    legend('logreg', 'l2 e-1', 'l2 e-2', 'l2 e-3', 'l2 e-4', 'l1 e-1', 'l1 e-2', 'l1 e-3', 'l1 e-4');
    ylabel('decoding accuracy (auc)');
    xlabel('time (ms)');
    title('Decoding of location during task stimulation');
    ylim([0.48 0.66])

    saveas(gcf, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/figure_allclassifiers.fig']);
    fig = openfig([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/figure_allclassifiers.fig']);
    saveas(fig, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/figure_allclassifiers.jpeg']);
    close all

end
    
