%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = ['/path_to_local/results/preprocessing/' folder_preproc];
    path_results = ['/path_to_local/results/mvpa/' folder_decoding];
end

number_stimuli = 8;
classifiers = {'l2_e-1'};
inputs = {'bins'};

%% figure decoding task best classifier

for input_i = 1:length(inputs)

    if strcmp(inputs{input_i}, 'bins')
        time = 25:5:24+(91)*5;
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


        cfg_stat = [];
        cfg_stat.n_permutations  = 5000;
    
        load([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
            '/stats_withingroup_perm' num2str(cfg_stat.n_permutations) '.mat']);

        % subset 300 ms of stimulus presentation
        subset_idx = 21:71;

        time = time(subset_idx);
        time = time - 100;
        result_average.perf{1} = result_average.perf{1}(subset_idx);
        result_average.perf_std{1} = result_average.perf_std{1}(subset_idx);
        
        
        try

            % plot the grand average result again and indicate the cluster in bold
            mv_plot_result(result_average, time, 'mask', stat_level2.mask(subset_idx));
            grid off
            yticks([0.5 0.55 0.6 0.65 0.7 0.75]);
            xticks([100 200]);
            title(' ');
            xlabel('time (ms)');
            ylabel('decoding performance (auc)');
%             xlim([25 375]);
            
            saveas(gcf, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                '/figure_stats_withingroup_perm' num2str(cfg_stat.n_permutations) '_v2.fig']);
            close all

        catch

            % plot the grand average result again and indicate the cluster in bold
            mv_plot_result(result_average, time);
            grid off
            yticks([0.5 0.55 0.6 0.65 0.7 0.75]);
            xticks([100 200]);
            title(' ');
            xlabel('time (ms)');
            ylabel('decoding performance (auc)');
%             xlim([25 375]);

            saveas(gcf, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                '/figure_stats_withingroup_perm' num2str(cfg_stat.n_permutations) '_v2.fig']);
            close all

        end

        fig = openfig([path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                '/figure_stats_withingroup_perm' num2str(cfg_stat.n_permutations) '_v2.fig']);
        saveas(fig, [path_results '/group_results/task_stim_balanced/' inputs{input_i} '/' classifier...
                '/figure_stats_withingroup_perm' num2str(cfg_stat.n_permutations) '_v2.jpeg']);
        close all

    end

end
    
disp('Analysis done');
