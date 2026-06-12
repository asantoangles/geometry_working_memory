%% GERE project
% replay

folder_replay = 'sequence_trial_avg_delay';

if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = 1:17;
end

path_data = [path_root, '/github/results/replay/data'];
path_results = [path_root, '/github/results/replay/supplementary_exclude1/' folder_replay];
addpath([path_root, '/github/scripts/replay/utilities'])

classifier = 'l2_e-1';
nlags = 100; % for plotting
number_stimuli = 8;
folders = {'task_stim2delay' 'localizer2delay'};
sessions = 1:2;

%% replay

for folder_i = 1:length(folders)

    disp(folders{folder_i})

    outdir_fig = [path_results '/group_results/' classifier '/' folders{folder_i} '/figures'];
    
    if ~isfolder(outdir_fig)
        mkdir(outdir_fig);
    end
    
    %% pool data

    % preallocate final outputs
    % dimensions:
    % (subject-session) x nlags
    
    n_subjectsessions = length(subjects) * length(sessions);
    
    Msf_l3_allsubjects = nan(n_subjectsessions, 6, nlags);
    Msf_l4_allsubjects = nan(n_subjectsessions, 24, nlags);
    
    Msb_l3_allsubjects = nan(n_subjectsessions, 6, nlags);
    Msb_l4_allsubjects = nan(n_subjectsessions, 24, nlags);
    
    row_i = 1;
    
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
    
        for ses_i = 1:length(sessions)
    
            session = sessions(ses_i);
    
            session_ID = ['sess_0' num2str(session)];
    
            indir = [path_data '/sequences/' subject_ID '/' session_ID];
    
            load([indir '/sequences.mat']); % sequences
    
            for trial_i = size(sequences.trialinfo,1):-1:1

                % subset predictions for trial in the
                % sequence
                sequence = sequences.trialinfo(trial_i, 4:7);
                sequence = sequence(sequence > 0);

                % exclude trials with consecutive
                sorted_sequence = unique(sort(sequence));
                exclude_trial = 0;

                for item_i = 2:length(sorted_sequence)

                    if abs(sorted_sequence(item_i) - sorted_sequence(item_i - 1)) == 1
                        exclude_trial = 1;
                    end
                    if sorted_sequence(1) == 1 && sorted_sequence(end) == 8
                        exclude_trial = 1;
                    end
                    if sorted_sequence(end) == 1 && sorted_sequence(1) == 8
                        exclude_trial = 1;
                    end

                end

                if exclude_trial == 1

                    sequences.trialinfo(trial_i, :) = [];

                end

            end

            outdir = [path_results '/' subject_ID '/' session_ID '/bins/' classifier '/TDLM_sequence/' folders{folder_i} '/prob/whole_delay'];
    
            % load sequenceness - empirical
            load([outdir '/Msf_empirical.mat']); % Msf_empirical
            load([outdir '/Msb_empirical.mat']); % Msb_empirical
    
            % ---------------------------------------------------------
            % find idx of trials of length 3 and 4
            % ---------------------------------------------------------
            trials_length3 = [];
            trials_length4 = [];
    
            for trial_i = 1:size(sequences.trialinfo,1)
    
                if sequences.trialinfo(trial_i,7) == 0
                    trials_length3 = [trials_length3 trial_i];
                else
                    trials_length4 = [trials_length4 trial_i];
                end
    
            end
    
            % ---------------------------------------------------------
            % convert selected trials to matrices
            % output: nTrials x nlags
            % ---------------------------------------------------------
            fwd_trial_lag = nan(length(trials_length3), nlags);
            bwd_trial_lag = nan(length(trials_length3), nlags);
            for trial_i = 1:length(trials_length3)
                fwd_trial_lag(trial_i,:) = Msf_empirical{trials_length3(trial_i)}(1:nlags);
                bwd_trial_lag(trial_i,:) = Msb_empirical{trials_length3(trial_i)}(1:nlags);
            end
            Msf_l3 = mean(fwd_trial_lag, 1);
            Msb_l3 = mean(bwd_trial_lag, 1);
            
            fwd_trial_lag = nan(length(trials_length4), nlags);
            bwd_trial_lag = nan(length(trials_length4), nlags);
            for trial_i = 1:length(trials_length4)
                fwd_trial_lag(trial_i,:) = Msf_empirical{trials_length4(trial_i)}(1:nlags);
                bwd_trial_lag(trial_i,:) = Msb_empirical{trials_length4(trial_i)}(1:nlags);
            end
            Msf_l4 = mean(fwd_trial_lag, 1);
            Msb_l4 = mean(bwd_trial_lag, 1);

            % ---------------------------------------------------------
            % average across trials
            % result: 1 x nlags
            % ---------------------------------------------------------
            Msf_l3_allsubjects(row_i, 1, :) = Msf_l3;
            Msb_l3_allsubjects(row_i, 1, :) = Msb_l3;
    
            Msf_l4_allsubjects(row_i, 1, :) = Msf_l4;
            Msb_l4_allsubjects(row_i, 1, :) = Msb_l4;


            % ---------------------------------------------------------
            % null distribution
            % ---------------------------------------------------------

            % load sequenceness - empirical
            load([outdir '/Msf_length3.mat']); % Msf_length3
            load([outdir '/Msb_length3.mat']); % Msb_length3

            load([outdir '/Msf_length4.mat']); % Msf_length4
            load([outdir '/Msb_length4.mat']); % Msb_length4

            Msf_l3_allsubjects(row_i, 2:end, :) = reshape(Msf_length3(:,1:nlags), [1 size(Msf_length3,1) nlags]);
            Msb_l3_allsubjects(row_i, 2:end, :) = reshape(Msb_length3(:,1:nlags), [1 size(Msb_length3,1) nlags]);
    
            Msf_l4_allsubjects(row_i, 2:end, :) = reshape(Msf_length4(:,1:nlags), [1 size(Msf_length4,1) nlags]);
            Msb_l4_allsubjects(row_i, 2:end, :) = reshape(Msb_length4(:,1:nlags), [1 size(Msb_length4,1) nlags]);
    
            row_i = row_i + 1;
    

        end
    
    end

    % peak
    tmp = reshape(Msf_l3_allsubjects, 2, [], size(Msf_l3_allsubjects,2), size(Msf_l3_allsubjects,3));
    Msf_l3_allsubjects_peak = squeeze(mean(tmp,1));
    Msf_l3_allsubjects_peak = squeeze(Msf_l3_allsubjects_peak(:,1,:));

    tmp = reshape(Msf_l4_allsubjects, 2, [], size(Msf_l4_allsubjects,2), size(Msf_l4_allsubjects,3));
    Msf_l4_allsubjects_peak = squeeze(mean(tmp,1));
    Msf_l4_allsubjects_peak = squeeze(Msf_l4_allsubjects_peak(:,1,:));
    
    tmp = reshape(Msb_l3_allsubjects, 2, [], size(Msb_l3_allsubjects,2), size(Msb_l3_allsubjects,3));
    Msb_l3_allsubjects_peak = squeeze(mean(tmp,1));
    Msb_l3_allsubjects_peak = squeeze(Msb_l3_allsubjects_peak(:,1,:));

    tmp = reshape(Msb_l4_allsubjects, 2, [], size(Msb_l4_allsubjects,2), size(Msb_l4_allsubjects,3));
    Msb_l4_allsubjects_peak = squeeze(mean(tmp,1));
    Msb_l4_allsubjects_peak = squeeze(Msb_l4_allsubjects_peak(:,1,:));

    %% plot l3

    null_Msf_thr_pos = mean(max(Msf_l3_allsubjects(:,2:end,:),[],2:3));
    null_Msf_thr_neg = mean(min(Msf_l3_allsubjects(:,2:end,:),[],2:3));
    null_Msb_thr_pos = mean(max(Msb_l3_allsubjects(:,2:end,:),[],2:3));
    null_Msb_thr_neg = mean(min(Msb_l3_allsubjects(:,2:end,:),[],2:3));

    % avoid full overlap of threshold lines
    if round(null_Msf_thr_pos,3) == round(null_Msb_thr_pos,3)
        null_Msf_thr_pos = null_Msf_thr_pos + 0.001;
    end

    if round(null_Msf_thr_neg,3) == round(null_Msb_thr_neg,3)
        null_Msf_thr_neg = null_Msf_thr_neg + 0.001;
    end

    % mean/median across subjects
    Msf_mean = squeeze(mean(Msf_l3_allsubjects(:,1,:), 1));
    Msb_mean = squeeze(mean(Msb_l3_allsubjects(:,1,:), 1));
    
    % SEM across subjects
    Msf_sem = squeeze(std(Msf_l3_allsubjects(:,1,:), [], 1)) ./ ...
              sqrt(size(Msf_l3_allsubjects,1));
    
    Msb_sem = squeeze(std(Msb_l3_allsubjects(:,1,:), [], 1)) ./ ...
              sqrt(size(Msb_l3_allsubjects,1));
    
    x = (1:nlags) * 5;
    
    figure
    
    % ---------------------------------------------------------
    % forward sequenceness
    % ---------------------------------------------------------
    fill([x fliplr(x)], ...
         [Msf_mean + Msf_sem; flipud(Msf_mean - Msf_sem)]', ...
         [0 0 0.8], ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none');
    hold on
    
    plot(x, Msf_mean, ...
        'Color', [0 0 0.8], ...
        'LineWidth', 4);
    
    % ---------------------------------------------------------
    % backward sequenceness
    % ---------------------------------------------------------
    fill([x fliplr(x)], ...
         [Msb_mean + Msb_sem; flipud(Msb_mean - Msb_sem)]', ...
         [0.8 0 0], ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none');
    
    plot(x, Msb_mean, ...
        'Color', [0.8 0 0], ...
        'LineWidth', 4);
           
    % ---------------------------------------------------------
    % significance threshold
    % ---------------------------------------------------------
    
    yline(null_Msf_thr_pos,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
    yline(null_Msf_thr_neg,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
    yline(null_Msb_thr_pos,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
    yline(null_Msb_thr_neg,'--', 'Color', [0.8 0 0], 'LineWidth', 3);

    % ---------------------------------------------------------
    % subject peak points - forward
    % ---------------------------------------------------------
    
    nSubjects = size(Msf_l3_allsubjects_peak,1);
    
    for sub_i = 1:nSubjects
    
        % peak value and index
        [peak_val, peak_idx] = max(Msf_l3_allsubjects_peak(sub_i,:));
    
        % convert lag index to time
        peak_time = x(peak_idx);
    
        % plot point
        scatter(peak_time, peak_val, ...
            80, ...
            [0 0 0.8], ...
            'filled', ...
            'MarkerFaceAlpha', 0.7);
    
    end
    
    % ---------------------------------------------------------
    % subject peak points - backward
    % ---------------------------------------------------------
    
    for sub_i = 1:nSubjects
    
        % peak value and index
        [peak_val, peak_idx] = max(Msb_l3_allsubjects_peak(sub_i,:));
    
        % convert lag index to time
        peak_time = x(peak_idx);
    
        % plot point
        scatter(peak_time, peak_val, ...
            80, ...
            [0.8 0 0], ...
            'filled', ...
            'MarkerFaceAlpha', 0.7);
    
    end
    
    % ---------------------------------------------------------
    % settings plot
    % ---------------------------------------------------------
    yticks(-0.20:0.10:0.20);
    xticks([]);

    % Increase the font size of the tick labels
    ax = gca; % Get current axes                            
    ax.FontSize = 20;
    
    saveas(gcf, [outdir_fig '/sequenceness_length3.fig']);
    fig = openfig([outdir_fig '/sequenceness_length3.fig']);
    saveas(fig, [outdir_fig '/sequenceness_length3.jpeg']);
    system(['rm ' outdir_fig '/sequenceness_length3.fig']);


    %% plot l4

    null_Msf_thr_pos = mean(max(Msf_l4_allsubjects(:,2:end,:),[],2:3));
    null_Msf_thr_neg = mean(min(Msf_l4_allsubjects(:,2:end,:),[],2:3));
    null_Msb_thr_pos = mean(max(Msb_l4_allsubjects(:,2:end,:),[],2:3));
    null_Msb_thr_neg = mean(min(Msb_l4_allsubjects(:,2:end,:),[],2:3));

    % avoid full overlap of threshold lines
    if round(null_Msf_thr_pos,3) == round(null_Msb_thr_pos,3)
        null_Msf_thr_pos = null_Msf_thr_pos + 0.001;
    end

    if round(null_Msf_thr_neg,3) == round(null_Msb_thr_neg,3)
        null_Msf_thr_neg = null_Msf_thr_neg + 0.001;
    end

    % mean/median across subjects
    Msf_mean = squeeze(mean(Msf_l4_allsubjects(:,1,:), 1));
    Msb_mean = squeeze(mean(Msb_l4_allsubjects(:,1,:), 1));
    
    % SEM across subjects
    Msf_sem = squeeze(std(Msf_l4_allsubjects(:,1,:), [], 1)) ./ ...
              sqrt(size(Msf_l4_allsubjects,1));
    
    Msb_sem = squeeze(std(Msb_l4_allsubjects(:,1,:), [], 1)) ./ ...
              sqrt(size(Msb_l4_allsubjects,1));
    
    x = (1:nlags) * 5;
    
    figure
    
    % ---------------------------------------------------------
    % forward sequenceness
    % ---------------------------------------------------------
    fill([x fliplr(x)], ...
         [Msf_mean + Msf_sem; flipud(Msf_mean - Msf_sem)]', ...
         [0 0 0.8], ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none');
    hold on
    
    plot(x, Msf_mean, ...
        'Color', [0 0 0.8], ...
        'LineWidth', 4);
    
    % ---------------------------------------------------------
    % backward sequenceness
    % ---------------------------------------------------------
    fill([x fliplr(x)], ...
         [Msb_mean + Msb_sem; flipud(Msb_mean - Msb_sem)]', ...
         [0.8 0 0], ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none');
    
    plot(x, Msb_mean, ...
        'Color', [0.8 0 0], ...
        'LineWidth', 4);
           
    % ---------------------------------------------------------
    % significance threshold
    % ---------------------------------------------------------
    
    yline(null_Msf_thr_pos,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
    yline(null_Msf_thr_neg,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
    yline(null_Msb_thr_pos,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
    yline(null_Msb_thr_neg,'--', 'Color', [0.8 0 0], 'LineWidth', 3);

    % ---------------------------------------------------------
    % subject peak points - forward
    % ---------------------------------------------------------
    
    nSubjects = size(Msf_l4_allsubjects_peak,1);
    
    for sub_i = 1:nSubjects
    
        % peak value and index
        [peak_val, peak_idx] = max(Msf_l4_allsubjects_peak(sub_i,:));
    
        % convert lag index to time
        peak_time = x(peak_idx);
    
        % plot point
        scatter(peak_time, peak_val, ...
            80, ...
            [0 0 0.8], ...
            'filled', ...
            'MarkerFaceAlpha', 0.7);
    
    end
    
    % ---------------------------------------------------------
    % subject peak points - backward
    % ---------------------------------------------------------
    
    for sub_i = 1:nSubjects
    
        % peak value and index
        [peak_val, peak_idx] = max(Msb_l4_allsubjects_peak(sub_i,:));
    
        % convert lag index to time
        peak_time = x(peak_idx);
    
        % plot point
        scatter(peak_time, peak_val, ...
            80, ...
            [0.8 0 0], ...
            'filled', ...
            'MarkerFaceAlpha', 0.7);
    
    end
    
    % ---------------------------------------------------------
    % settings plot
    % ---------------------------------------------------------
    yticks(-0.20:0.10:0.20);
    xticks([]);

    % Increase the font size of the tick labels
    ax = gca; % Get current axes                            
    ax.FontSize = 20;
    
    saveas(gcf, [outdir_fig '/sequenceness_length4.fig']);
    fig = openfig([outdir_fig '/sequenceness_length4.fig']);
    saveas(fig, [outdir_fig '/sequenceness_length4.jpeg']);
    system(['rm ' outdir_fig '/sequenceness_length4.fig']);

end

disp('Analysis done');