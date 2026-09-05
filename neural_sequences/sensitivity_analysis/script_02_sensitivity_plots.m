%% GERE project

% path to data
if isfolder('/path_to_local')
    path_root = '/path_to_local';
end

path_trialinfo = [path_root '/results/preprocessing/trialinfo'];
path_predictions = [path_root '/github/results/neural_sequences/data/predictions'];
path_sensitivity = [path_root '/github/results/neural_sequences/sensitivity_analysis'];
addpath([path_root '/github/scripts/neural_sequences/utilities']);

number_stimuli = 8;
classifiers = {'l2_e-1'}; classifier = classifiers{1};
time_delay = 791; % bins
subjects = 1:17;
sessions = 1:2;
number_iterations = 50;
nlags = 200;

folder_replay = 'sequence_trial_avg_delay';
folders = {'task_stim2delay', 'localizer2delay'};
amplitudes = [1.5 1.4 1.3 1.2 1.1 1 0.9 0.8 0.7 0.6 0.5];
lags = [100, 200, 300]; % milliseconds
percentiles = [85, 90, 95];


%% Load sensitivity results

sensitivity_data = nan( ...
    length(percentiles), ...
    length(folders), ...
    length(seq_lengths), ...
    length(amplitudes), ...
    length(lags), ...
    number_iterations);

for percentile_i = 1:length(percentiles)

    percentile = percentiles(percentile_i);

    for folder_i = 1:length(folders)

        folder_space = folders{folder_i};

        for sequence_length_i = 1:length(seq_lengths)

            % Actual sequence length (3 or 4)
            sequence_length = seq_lengths(sequence_length_i);

            for amplitude_i = 1:length(amplitudes)

                injected_percentage = amplitudes(amplitude_i);

                for lag_i = 1:length(lags)

                    lag_ms = lags(lag_i);

                    for iter_i = 1:number_iterations

                        outdir = [ ...
                            path_sensitivity '/' ...
                            folder_space '/percentile_' num2str(percentile) ...
                            '/injected_' num2str(injected_percentage * 100) ...
                            '/lag_' num2str(lag_ms) ...
                            '/iter_' num2str(iter_i)];

                        %% Load significance result

                        filename = [outdir ...
                            '/injected_significance_length' ...
                            num2str(sequence_length) '.mat'];

                        load(filename);

                        %% Extract appropriate sequence-length result

                        if sequence_length == 3

                            inj_sign = injected_significance_length3;

                        elseif sequence_length == 4

                            inj_sign = injected_significance_length4;

                        end

                        %% Store result

                        sensitivity_data( ...
                            percentile_i, ...
                            folder_i, ...
                            sequence_length_i, ...
                            amplitude_i, ...
                            lag_i, ...
                            iter_i) = inj_sign;

                    end

                end

            end

        end

    end

end

%% Plot sensitivity separately for L3 and L4
% Probability of significant TDLM effect
% Averaged across lags and iterations

% Sort amplitudes from smallest to largest
[amplitudes_sorted, sort_idx] = sort(amplitudes);

% Evenly spaced x-axis positions
x_positions = 1:length(amplitudes_sorted);

for folder_i = 1:length(folders)

    folder_space = folders{folder_i};

    if strcmp('task_stim2delay', folder_space)
        folder_space = 'task2delay';
    end

    for sequence_length_i = 1:length(seq_lengths)

        sequence_length = seq_lengths(sequence_length_i);

        figure;
        hold on;

        for percentile_i = 1:length(percentiles)

            % amplitude x lag x iteration
            data = squeeze(sensitivity_data( ...
                percentile_i, ...
                folder_i, ...
                sequence_length_i, ...
                :, :, :));

            % Sort amplitudes and corresponding sensitivity values
            data = data(sort_idx,:,:);

            % Average across lags and iterations
            detection_probability = squeeze(mean(data, [2 3], 'omitnan'));

            plot(x_positions, detection_probability, ...
                '-o', ...
                'LineWidth', 4, ...
                'MarkerSize', 6, ...
                'DisplayName', ...
                [num2str(percentiles(percentile_i)) 'th percentile']);

            xline(6, '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

        end

        xlabel('Injected amplitude (× empirical in-vs-out effect)');
        ylabel('Probability of significant TDLM effect');
        
        % title([folder_space ' - sequence length ' num2str(sequence_length)]);
        
        xticks(x_positions);
        xticklabels(string(amplitudes_sorted));
        
        ylim([0 1]);
        yticks(0:0.1:1);
        
        legend('Location', 'southeast');
        box off;
        set(gca, 'FontSize', 20);


        % -------------------------
        % Save figure
        % -------------------------
        outdir_fig = [path_sensitivity '/figures'];

        if ~isfolder(outdir_fig)
            mkdir(outdir_fig);
        end

        filename_fig = [outdir_fig ...
            '/sensitivity_' folder_space '_length' num2str(sequence_length) '.png'];

        exportgraphics(gcf, filename_fig, 'Resolution', 300);
        close(gcf);

    end

end

disp('  Analysis done');
