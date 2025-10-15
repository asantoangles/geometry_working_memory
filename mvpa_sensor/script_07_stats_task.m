
%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = ['/path_to_local/results/preprocessing/' folder_preproc];
    path_results = ['/path_to_local/results/mvpa/' folder_decoding];
end

number_stimuli = 8;
classifiers = {'l2_e-1'};
inputs = {'bins'};

%% stats task

for input_i = 1:length(inputs)

    if strcmp(inputs{input_i}, 'bins')
        time = 25:5:24+(91)*5;
    else
        time = 1:500;
    end
        
    % loop across classifiers
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


        %% group statistics within-group (different than auc=0.5)

        % For each subject and every time point, we have calculated AUC values. We
        % will now perform a cluster permutation test. See Maris & Oostenveld's
        % paper or the FieldTrip tutorials (https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/) 
        % for explanations on cluster permutation tests
        cfg_stat = [];
        cfg_stat.metric          = 'auc';
        cfg_stat.test            = 'permutation';
        cfg_stat.correctm        = 'cluster';  % correction method is cluster
        cfg_stat.n_permutations  = 5000;
        
        % Clusterstatistic is the actual statistic used for the clustertest.
        % Normally the default value 'maxum' is used, we are setting it here
        % explicitly for clarity. Maxsum adds up the statistics calculated at each
        % time point (the latter are set below using cfg_stat.statistic)
        cfg_stat.clusterstatistic = 'maxsum';
        cfg_stat.alpha           = 0.005; % use standard significance threshold of 5%
        
        % Level 2 stats design: we have to choose between within-subject and
        % between-subjects. Between-subjects is relevant when there is two
        % different experimental groups (eg patients vs controls) and we want to
        % investigate whether their MVPA results are significantly different. Here,
        % we have only one group and we want to see whether the AUC is
        % significantly different from a null value, hence the statistical design
        % is within-subject
        cfg_stat.design          = 'within';

        % cfg_stat.statistic defines how the difference between the AUC values and
        % the null is calculated at each time point (across subjects). 
        % We can choose t-test or its nonparametric counterpart Wilcoxon test. We
        % choose Wilcoxon here.
        cfg_stat.statistic       = 'wilcoxon';

        % The null value for AUC (corresponding to a random classifier) is 0.5
        cfg_stat.null            = 0.5;

        % clustercritval is a parameter that needs to be selected by the user. In a
        % Level 2 (group level) test, it represents the critical cutoff value for
        % the statistic. Here, we selected Wilcoxon, so clustercritval corresponds
        % to the cutoff value for the z-statistic which is obtained by a normal
        % approximation
        cfg_stat.clustercritval  = 2.56;
        % z-val = 1.65 corresponds to uncorrected p-value = 0.1
        % z-val = 1.96 corresponds to uncorrected p-value = 0.05
        % z-val = 2.58 corresponds to uncorrected p-value = 0.01

        if ~isfile([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                    '/figure_stats_withingroup_perm' num2str(cfg_stat.n_permutations) '.jpeg'])

            try
                stat_level2 = mv_statistics(cfg_stat, results_across_subjects);
            catch
                stat_level2 = [];
            end
    
            if ~isfolder([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier])
                mkdir([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier]);
            end
    
            save([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                '/stats_withingroup_perm' num2str(cfg_stat.n_permutations) '.mat'], 'stat_level2');
            
            try
    
                % plot the grand average result again and indicate the cluster in bold
                mv_plot_result(result_average, time, 'mask', stat_level2.mask);
                saveas(gcf, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                    '/figure_stats_withingroup_perm' num2str(cfg_stat.n_permutations) '.fig']);
                close all
    
            catch
    
                % plot the grand average result again and indicate the cluster in bold
                mv_plot_result(result_average, time);
                saveas(gcf, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                    '/figure_stats_withingroup_perm' num2str(cfg_stat.n_permutations) '.fig']);
                close all
    
            end
    
            fig = openfig([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                    '/figure_stats_withingroup_perm' num2str(cfg_stat.n_permutations) '.fig']);

            % Add a horizontal line at y = some_value
            some_value = 0.5;  % Change this value accordingly
            plot(xlim, [some_value, some_value], '--', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
    
            xline(100, '-', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
            xline(400, '-', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');

            title(['task stim ' classifiers{class_i} ' ' inputs{input_i}]);
    
            hold off;


            saveas(fig, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                    '/figure_stats_withingroup_perm' num2str(cfg_stat.n_permutations) '.jpeg']);
            close all

        end

        %% time generalization

        try

            if ~isfile([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                '/figure_timextime_average_locations.jpeg'])
    
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
                                     
                            load([path_results '/' subject_ID '/' session_ID '/task_stim_balanced/' inputs{input_i} '/' classifier...
                                '/timextime_location' num2str(stim_i) '.mat']); % timextime
        
                            results_across_stimuli{stim_i} = timextime;
        
                        end
        
                        % plot average across locations
                        result_average = mv_combine_results(results_across_stimuli, 'average');
                        result_average = mv_select_result(result_average, 'auc');
        
                        results_across_subjects{end+1} = result_average;
                            
                    end
                
                end
        
                result_average = mv_combine_results(results_across_subjects, 'average');
                        
                % plot average across locations
                close all
                imagesc(1:size(result_average.perf{1},1), 1:size(result_average.perf{1},1), result_average.perf{1});
                set(gca, 'YDir', 'normal');
                colorbar; grid on;
                ylabel('Training time'), xlabel('Test time');
                title(['task stim ' classifiers{class_i} ' averaged locations / subjects']);
                saveas(gcf, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                    '/figure_timextime_average_locations.fig'])
                fig = openfig([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                    '/figure_timextime_average_locations.fig']);                
                saveas(fig, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                    '/figure_timextime_average_locations.jpeg']);
                close all
    
            end

        catch
        end

    end

end
    
disp('Analysis done');

