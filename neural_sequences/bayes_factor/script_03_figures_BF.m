%% GERE project
% replay

if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = 1:17;
end

path_data = [path_root, '/github/results/replay/data'];
addpath([path_root, '/github/scripts/replay/utilities'])

classifier = 'l2_e-1';
nlags = 100;
number_stimuli = 8;
sessions = 1:2;

seq_lengths = [3 4];

replay_folders = {'sequence_concat_delay', 'sequence_concat_stim', 'sequence_trial_avg_delay', 'sequence_trial_avg_stim'};

% ATTENTION
replay_folders = {'sequence_concat_delay' 'sequence_trial_avg_delay'};
% ATTENTION


%% bayes factor

% loop over replay approaches
for folder_replay_i = 1:length(replay_folders)

    folder_replay = replay_folders{folder_replay_i};
    path_results = [path_root, '/github/results/replay/' folder_replay];

    if strcmp(folder_replay((end-4):end), 'delay')
        folders = {'task_stim2delay' 'localizer2delay'};
    else
        folders = {'task_stim2task_stim' 'localizer2task_stim'};
    end

    % loop over folders
    for folder_i = 1:length(folders)
    
        path_outputs = [path_root, '/github/results/replay/' folder_replay '/group_results/' classifier '/' folders{folder_i} '/bayes_factor'];
    
        if ~isfolder(path_outputs)
            mkdir(path_outputs);
        end
            
        %% bayes factor
    
        % loop over sequence lengths
        for length_i = seq_lengths
    
            %% ---------------------------------------------------------
            % Plot BF10 and BF01 across lags
            %
            % BF10  -> solid lines
            % BF01  -> dashed lines
            % forward  -> blue
            % backward -> red
            %% ---------------------------------------------------------
            
            figure; hold on;
            
            for fwd_bwd_i = 1:2 % 1:2 % ATTENTION
            
                if fwd_bwd_i == 1
                    seq_direction = 'f';
                    seq_name = 'forward';
            
                    line_color = [0 0 1]; % blue
            
                else
                    seq_direction = 'b';
                    seq_name = 'backward';
            
                    line_color = [1 0 0]; % red

                end

                disp(' ')
                disp('######################################################')
                disp([folder_replay ' - length ' num2str(length_i) ' - ' seq_name]);
                disp('######################################################')

                % -----------------------------------------------------
                % load BF results
                % -----------------------------------------------------
                load([path_outputs '/bf10_lags_' seq_name '_length' num2str(length_i) '.mat']); % bf10

                bf10_clean = bf10;
                
                % replace zeros (and any negatives just in case) with small floor value
                bf10_clean(bf10_clean <= 0) = 1e-4;
                
                % log transform
                log_bf = 2 * log(bf10_clean);

                % -----------------------------------------------------
                % plot
                % -----------------------------------------------------
                plot(log_bf, ...
                     'LineWidth', 3, ...
                     'Color', line_color, ...
                     'LineStyle', '-');



                %% =========================================================
                % Summary statistics for Bayes Factors across lags
                %% =========================================================
                
                lags = 1:nlags;
                lags_real = linspace(1,500,length(lags));
                
                mean_bf = mean(log_bf);
                std_bf  = std(log_bf);
                
                median_bf = median(log_bf);
                
                iqr_bf = iqr(log_bf);
                
                q25_bf = prctile(log_bf,25);
                q75_bf = prctile(log_bf,75);
                
                [min_bf, idx_min] = min(log_bf);
                [max_bf, idx_max] = max(log_bf);
                
                % Times of min/max peaks in rescaled lag axis
                time_min_bf = lags_real(idx_min);
                time_max_bf = lags_real(idx_max);
                
                %% =========================================================
                % Find intervals where logBF > 10
                %% =========================================================
                
                threshold = 10;
                
                pos_mask = log_bf > threshold;
                
                d_pos = diff([0 pos_mask 0]);
                
                start_idx_pos = find(d_pos == 1);
                end_idx_pos   = find(d_pos == -1) - 1;
                
                pos_intervals = [lags_real(start_idx_pos)' ...
                                 lags_real(end_idx_pos)'];
                
                %% =========================================================
                % Find intervals where logBF < -10
                %% =========================================================
                
                neg_mask = log_bf < -threshold;
                
                d_neg = diff([0 neg_mask 0]);
                
                start_idx_neg = find(d_neg == 1);
                end_idx_neg   = find(d_neg == -1) - 1;
                
                neg_intervals = [lags_real(start_idx_neg)' ...
                                 lags_real(end_idx_neg)'];
                
                %% ---------------------------------------------------------
                % Print summary
                %% ---------------------------------------------------------
                
                fprintf('Mean logBF = %.1f ± %.1f [%.1f %.1f]\n', ...
                    mean_bf, std_bf, min_bf, max_bf);
                
                fprintf('Min peak: logBF = %.1f at %.0f ms\n', ...
                        min_bf, time_min_bf);
                
                fprintf('Max peak: logBF = %.1f at %.0f ms\n', ...
                         max_bf, time_max_bf);
                
                
                %% Positive intervals
                if isempty(pos_intervals)
                    fprintf('No intervals with logBF > %.1f\n', threshold);
                else
                    fprintf('Intervals with logBF > %.1f:\n', threshold);
                    for i = 1:size(pos_intervals,1)
                        fprintf('  %.1f to %.1f\n', ...
                            pos_intervals(i,1), pos_intervals(i,2));
                    end
                end
                
                %% Negative intervals
                if isempty(neg_intervals)
                    fprintf('No intervals with logBF < -%.1f\n', threshold);
                else
                    fprintf('Intervals with logBF < -%.1f:\n', threshold);
                    for i = 1:size(neg_intervals,1)
                        fprintf('  %.1f to %.1f\n', ...
                            neg_intervals(i,1), neg_intervals(i,2));
                    end
                end

            end

            % Very strong evidence (|log BF| = 6)
            yline( 6, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
            yline(-6, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
            
            % Moderate evidence (|log BF| = 2)
            yline( 2, ':',  'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);
            yline(-2, ':',  'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);
            
            % Neutral reference (log BF = 0)
            yline(0, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);

            %% ---------------------------------------------------------
            % figure formatting
            %% ---------------------------------------------------------

            % ylabel('2 log_e(BF)');
            
            yl = ylim;
            ylim([min(yl(1), -10), max(yl(2), 10)]);

            xticks([]);
            % xticks([1 50 100]);
            % xticklabels({'0','250','500'});
            % xlabel('time lag (ms)');

            % legend({'BF Forward', 'BF Backward'}, 'Location', 'best');
            
            box off;
            set(gca, 'FontSize', 15);

            set(gcf, 'Units', 'centimeters');
            set(gcf, 'Position', [0 0 8 6]); % [x y width height]
            
            % save
            exportgraphics(gcf, [path_outputs '/figure_length' num2str(length_i) '.png'], 'Resolution', 300);


        
        end
    
    end

end

close all

disp('Analysis done');

