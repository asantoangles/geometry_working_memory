%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = '/path_to_local/results/source_reconstruction';
    path_results = ['/path_to_local/results/source_geometry_lm/' folder_geometry];
    addpath('/path_to_local/scripts/source_geometry_lm/utilities')
    functions_data ={@mean};
end

% settings loops
segments_start = 1:300:3701;
segments_end = 300:300:4000;
delay_segments = [segments_start; segments_end]';

stim_resolved_segments = [];
stim_resolved_segments(end+1,:) = [1 300];
stim_resolved_segments(end+1,:) = [400 700];
stim_resolved_segments(end+1,:) = [800 1100];
stim_resolved_segments(end+1,:) = [1200 1500];

performance = {'correct_trials'};
iterations_bootstrap = 10;

number_locations = 8;

components = 1:3;

sessions = 1:2;

number_sensors = 200;

number_randomizatons = 1000;
subset_eigenvalues = 20;

%% empirical - surrogate / length 3 and 4 in one figure - delay

events = {'delay'};

for event_i = 1:length(events)

    if strcmp(events{event_i}, 'delay')

        time_segments = delay_segments;

    elseif strcmp(events{event_i}, 'baseline')

        time_segments = baseline_segments;

    elseif strcmp(events{event_i}, 'stim')
        
        time_segments = stim_segments;

    end

    for perf_i = 1:length(performance)
            
        % loop over functions
        for fun_i = functions_data

            fun_i = fun_i{1};
                    

            path_figure = [path_results '/group_results/variance_explained/' func2str(fun_i) '/' performance{perf_i} '/delay_avg'];

            if ~isfolder(path_figure)
                mkdir(path_figure);
            end
                
            filename_figure = [path_figure '/figure_empirical_surrogate.png'];

            if ~isfile(filename_figure)

                % Create the plot
                figure('Position', [100, 100, 250, 400]);  % [left, bottom, width, height]
                hold on;
                
                x = 1:number_sensors;

                %% empirical

                % loop over sequence length used
                for sequence_length = [3 4]
    
                    if sequence_length == 3
                        sequence_length_filename = 'length3';
                        ranks = 3;
                    elseif sequence_length == 4
                        sequence_length_filename = 'length4';
                        ranks = 4;
                    elseif sequence_length == 34
                        sequence_length_filename = 'lengthall';
                        ranks = 4;
                    end

                    % variance explained
                    variance_explained_allsubjects = nan(length(subjects), number_sensors);

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

                        % subset by time windows within delay period
                        for delay_i = 1:size(time_segments, 1)
                                                                                
                            % variance explained
                            load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/variance_explained_' sequence_length_filename '.mat']); % variance_explained

                            if delay_i == 1
    
                                var_explained = 100 * (variance_explained.z_k / sum(variance_explained.z_k));

                            else

                                var_explained = [var_explained, 100 * (variance_explained.z_k / sum(variance_explained.z_k))];

                            end
                                
                        end

                        var_explained = mean(var_explained,2);

                        variance_explained_allsubjects(sub_i, :) = mean(var_explained,2);

                    end

                    if sequence_length == 3
                        variance_explained_allsubjects_l3 = variance_explained_allsubjects;
                    else
                        variance_explained_allsubjects_l4 = variance_explained_allsubjects;
                    end

                end

                % accumulated
                for col_i = 2:size(variance_explained_allsubjects_l3,2)
                    variance_explained_allsubjects_l3(:,col_i) = variance_explained_allsubjects_l3(:,col_i) + variance_explained_allsubjects_l3(:,col_i-1);
                end
                for col_i = 2:size(variance_explained_allsubjects_l4,2)
                    variance_explained_allsubjects_l4(:,col_i) = variance_explained_allsubjects_l4(:,col_i) + variance_explained_allsubjects_l4(:,col_i-1);
                end

                % variance explained
                disp('empirical');
                k = 3;
                disp(['variance explained by ' num2str(k) ' PCs in length 3: ' num2str(mean(variance_explained_allsubjects_l3(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l3(:,k))) ])
                disp(['variance explained by ' num2str(k) ' PCs in length 4: ' num2str(mean(variance_explained_allsubjects_l4(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l4(:,k))) ])
                k = 5;
                disp(['variance explained by ' num2str(k) ' PCs in length 3: ' num2str(mean(variance_explained_allsubjects_l3(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l3(:,k))) ])
                disp(['variance explained by ' num2str(k) ' PCs in length 4: ' num2str(mean(variance_explained_allsubjects_l4(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l4(:,k))) ])

                % plot
                var_mean_l3 = mean(variance_explained_allsubjects_l3);
                var_std_l3 = std(variance_explained_allsubjects_l3);
                var_mean_l4 = mean(variance_explained_allsubjects_l4);
                var_std_l4 = std(variance_explained_allsubjects_l4);

                % Calculate the upper and lower bounds of the shaded area
                upper_bound_l3 = var_mean_l3 + var_std_l3;
                lower_bound_l3 = var_mean_l3 - var_std_l3;
                upper_bound_l4 = var_mean_l4 + var_std_l4;
                lower_bound_l4 = var_mean_l4 - var_std_l4;
                    
                % subset plots
                x = x(1:subset_eigenvalues);
                var_mean_l3 = var_mean_l3(1:subset_eigenvalues);
                var_std_l3 = var_std_l3(1:subset_eigenvalues);
                upper_bound_l3 = upper_bound_l3(1:subset_eigenvalues);
                lower_bound_l3 = lower_bound_l3(1:subset_eigenvalues);
                var_mean_l4 = var_mean_l4(1:subset_eigenvalues);
                var_std_l4 = var_std_l4(1:subset_eigenvalues);
                upper_bound_l4 = upper_bound_l4(1:subset_eigenvalues);
                lower_bound_l4 = lower_bound_l4(1:subset_eigenvalues);

                % Plot length 3
                fill([x, fliplr(x)], [upper_bound_l3, fliplr(lower_bound_l3)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                plot(x, var_mean_l3, 'r', 'LineWidth', 1.2);
                hold on

                % Plot length 4
                fill([x, fliplr(x)], [upper_bound_l4, fliplr(lower_bound_l4)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                plot(x, var_mean_l4, 'b', 'LineWidth', 1.2);
                hold on

                %% surrogate

                % loop over sequence length used
                for sequence_length = [3 4]
    
                    if sequence_length == 3
                        sequence_length_filename = 'length3';
                        ranks = 3;
                    elseif sequence_length == 4
                        sequence_length_filename = 'length4';
                        ranks = 4;
                    elseif sequence_length == 34
                        sequence_length_filename = 'lengthall';
                        ranks = 4;
                    end

                    % variance explained
                    variance_explained_allsubjects = nan(length(subjects), number_sensors);

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
                                    
                        for rand_i = 1:number_randomizatons
                                    
                            % subset by time windows within delay period
                            for delay_i = 1:size(time_segments, 1)
                                                                                    
                                % variance explained
                                load([path_results '_surrogate/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/variance_explained_' sequence_length_filename '.mat']); % variance_explained
    
                                if delay_i == 1
        
                                    var_explained_tmp = 100 * (variance_explained.z_k / sum(variance_explained.z_k));
    
                                else
    
                                    var_explained_tmp = [var_explained_tmp, 100 * (variance_explained.z_k / sum(variance_explained.z_k))];
    
                                end
                                    
                            end

                            % compute variance explained
                            if rand_i == 1
                                var_explained = mean(var_explained_tmp,2);
                            else
                                var_explained = [var_explained mean(var_explained_tmp,2)];
                            end

                        end

                        variance_explained_allsubjects(sub_i, :) = mean(var_explained,2);
                        
                    end

                    if sequence_length == 3
                        variance_explained_allsubjects_l3 = variance_explained_allsubjects;
                    else
                        variance_explained_allsubjects_l4 = variance_explained_allsubjects;
                    end

                end

                % accumulated
                for col_i = 2:size(variance_explained_allsubjects_l3,2)
                    variance_explained_allsubjects_l3(:,col_i) = variance_explained_allsubjects_l3(:,col_i) + variance_explained_allsubjects_l3(:,col_i-1);
                end
                for col_i = 2:size(variance_explained_allsubjects_l4,2)
                    variance_explained_allsubjects_l4(:,col_i) = variance_explained_allsubjects_l4(:,col_i) + variance_explained_allsubjects_l4(:,col_i-1);
                end

                % variance explained
                disp('surrogate');
                k = 3;
                disp(['variance explained by ' num2str(k) ' PCs in length 3: ' num2str(mean(variance_explained_allsubjects_l3(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l3(:,k))) ])
                disp(['variance explained by ' num2str(k) ' PCs in length 4: ' num2str(mean(variance_explained_allsubjects_l4(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l4(:,k))) ])
                k = 5;
                disp(['variance explained by ' num2str(k) ' PCs in length 3: ' num2str(mean(variance_explained_allsubjects_l3(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l3(:,k))) ])
                disp(['variance explained by ' num2str(k) ' PCs in length 4: ' num2str(mean(variance_explained_allsubjects_l4(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l4(:,k))) ])
                
                % plot
                var_mean_l3 = mean(variance_explained_allsubjects_l3);
                var_std_l3 = std(variance_explained_allsubjects_l3);
                var_mean_l4 = mean(variance_explained_allsubjects_l4);
                var_std_l4 = std(variance_explained_allsubjects_l4);

                % Calculate the upper and lower bounds of the shaded area
                upper_bound_l3 = var_mean_l3 + var_std_l3;
                lower_bound_l3 = var_mean_l3 - var_std_l3;
                upper_bound_l4 = var_mean_l4 + var_std_l4;
                lower_bound_l4 = var_mean_l4 - var_std_l4;
                    
                % subset plots
                x = x(1:subset_eigenvalues);
                var_mean_l3 = var_mean_l3(1:subset_eigenvalues);
                var_std_l3 = var_std_l3(1:subset_eigenvalues);
                upper_bound_l3 = upper_bound_l3(1:subset_eigenvalues);
                lower_bound_l3 = lower_bound_l3(1:subset_eigenvalues);
                var_mean_l4 = var_mean_l4(1:subset_eigenvalues);
                var_std_l4 = var_std_l4(1:subset_eigenvalues);
                upper_bound_l4 = upper_bound_l4(1:subset_eigenvalues);
                lower_bound_l4 = lower_bound_l4(1:subset_eigenvalues);

                % Plot length 3
                fill([x, fliplr(x)], [upper_bound_l3, fliplr(lower_bound_l3)], [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                plot(x, var_mean_l3, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
                hold on

                % Plot length 4
                fill([x, fliplr(x)], [upper_bound_l4, fliplr(lower_bound_l4)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                plot(x, var_mean_l4, 'k', 'LineWidth', 1.2);
                hold on

                ylim([0 100]);
                yticks(0:25:100);

                ylabel('variance explained (%)', 'FontSize', 12);
                xlabel('eigenvectors', 'FontSize', 12);
                
                %% legend and save
                legend('', 'Empirical length 3', '', 'Empirical length 4', '', 'Surrogate length 3', '', 'Surrogate length 4', 'Location', 'southeast');

                ax = gca;  % Get current axis
                ax.FontSize = 12;  % Set font size of the axis ticks

                if ~isfolder(path_figure)
                    mkdir(path_figure);
                end
                
                saveas(gcf, filename_figure);

            end

            close all

        end
    
    end

end

%% empirical - surrogate / length 3 and 4 in one figure - stim

events = {'stim_resolved_refined'};

for event_i = 1:length(events)

    if strcmp(events{event_i}, 'delay')

        time_segments = delay_segments;

    elseif strcmp(events{event_i}, 'baseline')

        time_segments = baseline_segments;

    elseif strcmp(events{event_i}, 'stim_resolved_refined')
        
        time_segments = stim_resolved_segments;

    end

    for perf_i = 1:length(performance)
            
        % loop over functions
        for fun_i = functions_data

            fun_i = fun_i{1};
                    

            path_figure = [path_results '/group_results/variance_explained/' func2str(fun_i) '/' performance{perf_i} '/stim_avg'];

            if ~isfolder(path_figure)
                mkdir(path_figure);
            end
                
            filename_figure = [path_figure '/figure_empirical_surrogate.png'];

            if ~isfile(filename_figure)

                % Create the plot
                figure('Position', [100, 100, 250, 400]);  % [left, bottom, width, height]
                hold on;
                
                x = 1:number_sensors;

                %% empirical

                % loop over sequence length used
                for sequence_length = [3 4]
    
                    if sequence_length == 3
                        sequence_length_filename = 'length3';
                        ranks = 3;
                    elseif sequence_length == 4
                        sequence_length_filename = 'length4';
                        ranks = 4;
                    elseif sequence_length == 34
                        sequence_length_filename = 'lengthall';
                        ranks = 4;
                    end

                    % variance explained
                    variance_explained_allsubjects = nan(length(subjects), number_sensors);

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

                        % subset by time windows within delay period
                        for delay_i = 1:size(time_segments, 1)
                                                                                
                            % variance explained
                            load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                '/variance_explained_' sequence_length_filename '.mat']); % variance_explained

                            if delay_i == 1
    
                                var_explained = 100 * (variance_explained.z_k / sum(variance_explained.z_k));

                            else

                                var_explained = [var_explained, 100 * (variance_explained.z_k / sum(variance_explained.z_k))];

                            end
                                
                        end

                        var_explained = mean(var_explained,2);

                        variance_explained_allsubjects(sub_i, :) = mean(var_explained,2);

                    end

                    if sequence_length == 3
                        variance_explained_allsubjects_l3 = variance_explained_allsubjects;
                    else
                        variance_explained_allsubjects_l4 = variance_explained_allsubjects;
                    end

                end

                % accumulated
                for col_i = 2:size(variance_explained_allsubjects_l3,2)
                    variance_explained_allsubjects_l3(:,col_i) = variance_explained_allsubjects_l3(:,col_i) + variance_explained_allsubjects_l3(:,col_i-1);
                end
                for col_i = 2:size(variance_explained_allsubjects_l4,2)
                    variance_explained_allsubjects_l4(:,col_i) = variance_explained_allsubjects_l4(:,col_i) + variance_explained_allsubjects_l4(:,col_i-1);
                end

                % variance explained
                disp('empirical');
                k = 3;
                disp(['variance explained by ' num2str(k) ' PCs in length 3: ' num2str(mean(variance_explained_allsubjects_l3(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l3(:,k))) ])
                disp(['variance explained by ' num2str(k) ' PCs in length 4: ' num2str(mean(variance_explained_allsubjects_l4(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l4(:,k))) ])
                k = 5;
                disp(['variance explained by ' num2str(k) ' PCs in length 3: ' num2str(mean(variance_explained_allsubjects_l3(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l3(:,k))) ])
                disp(['variance explained by ' num2str(k) ' PCs in length 4: ' num2str(mean(variance_explained_allsubjects_l4(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l4(:,k))) ])

                % plot
                var_mean_l3 = mean(variance_explained_allsubjects_l3);
                var_std_l3 = std(variance_explained_allsubjects_l3);
                var_mean_l4 = mean(variance_explained_allsubjects_l4);
                var_std_l4 = std(variance_explained_allsubjects_l4);

                % Calculate the upper and lower bounds of the shaded area
                upper_bound_l3 = var_mean_l3 + var_std_l3;
                lower_bound_l3 = var_mean_l3 - var_std_l3;
                upper_bound_l4 = var_mean_l4 + var_std_l4;
                lower_bound_l4 = var_mean_l4 - var_std_l4;
                    
                % subset plots
                x = x(1:subset_eigenvalues);
                var_mean_l3 = var_mean_l3(1:subset_eigenvalues);
                var_std_l3 = var_std_l3(1:subset_eigenvalues);
                upper_bound_l3 = upper_bound_l3(1:subset_eigenvalues);
                lower_bound_l3 = lower_bound_l3(1:subset_eigenvalues);
                var_mean_l4 = var_mean_l4(1:subset_eigenvalues);
                var_std_l4 = var_std_l4(1:subset_eigenvalues);
                upper_bound_l4 = upper_bound_l4(1:subset_eigenvalues);
                lower_bound_l4 = lower_bound_l4(1:subset_eigenvalues);

                % Plot length 3
                fill([x, fliplr(x)], [upper_bound_l3, fliplr(lower_bound_l3)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                plot(x, var_mean_l3, 'r', 'LineWidth', 1.2);
                hold on

                % Plot length 4
                fill([x, fliplr(x)], [upper_bound_l4, fliplr(lower_bound_l4)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                plot(x, var_mean_l4, 'b', 'LineWidth', 1.2);
                hold on

                %% surrogate

                % loop over sequence length used
                for sequence_length = [3 4]
    
                    if sequence_length == 3
                        sequence_length_filename = 'length3';
                        ranks = 3;
                    elseif sequence_length == 4
                        sequence_length_filename = 'length4';
                        ranks = 4;
                    elseif sequence_length == 34
                        sequence_length_filename = 'lengthall';
                        ranks = 4;
                    end

                    % variance explained
                    variance_explained_allsubjects = nan(length(subjects), number_sensors);

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
                                    
                        for rand_i = 1:number_randomizatons
                                    
                            % subset by time windows within delay period
                            for delay_i = 1:size(time_segments, 1)
                                                                                    
                                % variance explained
                                load([path_results '_surrogate/stim_rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/variance_explained_' sequence_length_filename '.mat']); % variance_explained
    
                                if delay_i == 1
        
                                    var_explained_tmp = 100 * (variance_explained.z_k / sum(variance_explained.z_k));
    
                                else
    
                                    var_explained_tmp = [var_explained_tmp, 100 * (variance_explained.z_k / sum(variance_explained.z_k))];
    
                                end
                                    
                            end

                            % compute variance explained
                            if rand_i == 1
                                var_explained = mean(var_explained_tmp,2);
                            else
                                var_explained = [var_explained mean(var_explained_tmp,2)];
                            end

                        end

                        variance_explained_allsubjects(sub_i, :) = mean(var_explained,2);
                        
                    end

                    if sequence_length == 3
                        variance_explained_allsubjects_l3 = variance_explained_allsubjects;
                    else
                        variance_explained_allsubjects_l4 = variance_explained_allsubjects;
                    end

                end

                % accumulated
                for col_i = 2:size(variance_explained_allsubjects_l3,2)
                    variance_explained_allsubjects_l3(:,col_i) = variance_explained_allsubjects_l3(:,col_i) + variance_explained_allsubjects_l3(:,col_i-1);
                end
                for col_i = 2:size(variance_explained_allsubjects_l4,2)
                    variance_explained_allsubjects_l4(:,col_i) = variance_explained_allsubjects_l4(:,col_i) + variance_explained_allsubjects_l4(:,col_i-1);
                end

                % variance explained
                disp('surrogate');
                k = 3;
                disp(['variance explained by ' num2str(k) ' PCs in length 3: ' num2str(mean(variance_explained_allsubjects_l3(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l3(:,k))) ])
                disp(['variance explained by ' num2str(k) ' PCs in length 4: ' num2str(mean(variance_explained_allsubjects_l4(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l4(:,k))) ])
                k = 5;
                disp(['variance explained by ' num2str(k) ' PCs in length 3: ' num2str(mean(variance_explained_allsubjects_l3(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l3(:,k))) ])
                disp(['variance explained by ' num2str(k) ' PCs in length 4: ' num2str(mean(variance_explained_allsubjects_l4(:,k))) ' ± ' num2str(std(variance_explained_allsubjects_l4(:,k))) ])
                
                % plot
                var_mean_l3 = mean(variance_explained_allsubjects_l3);
                var_std_l3 = std(variance_explained_allsubjects_l3);
                var_mean_l4 = mean(variance_explained_allsubjects_l4);
                var_std_l4 = std(variance_explained_allsubjects_l4);

                % Calculate the upper and lower bounds of the shaded area
                upper_bound_l3 = var_mean_l3 + var_std_l3;
                lower_bound_l3 = var_mean_l3 - var_std_l3;
                upper_bound_l4 = var_mean_l4 + var_std_l4;
                lower_bound_l4 = var_mean_l4 - var_std_l4;
                    
                % subset plots
                x = x(1:subset_eigenvalues);
                var_mean_l3 = var_mean_l3(1:subset_eigenvalues);
                var_std_l3 = var_std_l3(1:subset_eigenvalues);
                upper_bound_l3 = upper_bound_l3(1:subset_eigenvalues);
                lower_bound_l3 = lower_bound_l3(1:subset_eigenvalues);
                var_mean_l4 = var_mean_l4(1:subset_eigenvalues);
                var_std_l4 = var_std_l4(1:subset_eigenvalues);
                upper_bound_l4 = upper_bound_l4(1:subset_eigenvalues);
                lower_bound_l4 = lower_bound_l4(1:subset_eigenvalues);

                % Plot length 3
                fill([x, fliplr(x)], [upper_bound_l3, fliplr(lower_bound_l3)], [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                plot(x, var_mean_l3, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
                hold on

                % Plot length 4
                fill([x, fliplr(x)], [upper_bound_l4, fliplr(lower_bound_l4)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                plot(x, var_mean_l4, 'k', 'LineWidth', 1.2);
                hold on

                ylim([0 100]);
                yticks(0:25:100);

                ylabel('variance explained (%)', 'FontSize', 12);
                xlabel('eigenvectors', 'FontSize', 12);
                
                %% legend and save
                legend('', 'Empirical length 3', '', 'Empirical length 4', '', 'Surrogate length 3', '', 'Surrogate length 4', 'Location', 'southeast');

                ax = gca;  % Get current axis
                ax.FontSize = 12;  % Set font size of the axis ticks

                if ~isfolder(path_figure)
                    mkdir(path_figure);
                end
                
                saveas(gcf, filename_figure);

            end

            close all

        end
    
    end

end


disp('Analysis done');
