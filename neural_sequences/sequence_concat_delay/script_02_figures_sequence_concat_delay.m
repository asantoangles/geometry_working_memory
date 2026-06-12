%% GERE project
% replay

folder_replay = 'sequence_concat_delay';

if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = 1:17;
end

path_data = [path_root, '/github/results/replay/data'];
path_results = [path_root, '/github/results/replay/' folder_replay];
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

        outdir = [path_results '/' subject_ID '/' classifier '/' folders{folder_i}];
            
        load([outdir '/Msf_l3.mat']);
        load([outdir '/Msf_l4.mat']);
        load([outdir '/Msb_l3.mat']);
        load([outdir '/Msb_l4.mat']);                            

        Msf_l3 = Msf_l3(:,:,1:(end/2));
        Msb_l3 = Msb_l3(:,:,1:(end/2));
        Msf_l4 = Msf_l4(:,:,1:(end/2));
        Msb_l4 = Msb_l4(:,:,1:(end/2));

        if sub_i == 1

            Msf_l3_allsubjects = Msf_l3;
            Msb_l3_allsubjects = Msb_l3;
            Msf_l4_allsubjects = Msf_l4;
            Msb_l4_allsubjects = Msb_l4;

        else

            Msf_l3_allsubjects(end+1,:,:) = Msf_l3;
            Msb_l3_allsubjects(end+1,:,:) = Msb_l3;
            Msf_l4_allsubjects(end+1,:,:) = Msf_l4;
            Msb_l4_allsubjects(end+1,:,:) = Msb_l4;

        end

    end

    % peak
    Msf_l3_allsubjects_peak = squeeze(Msf_l3_allsubjects(:,1,:));
    Msf_l4_allsubjects_peak = squeeze(Msf_l4_allsubjects(:,1,:));
    Msb_l3_allsubjects_peak = squeeze(Msb_l3_allsubjects(:,1,:));
    Msb_l4_allsubjects_peak = squeeze(Msb_l4_allsubjects(:,1,:));

    %% plot l3

    null_Msf_thr_pos = median(max(Msf_l3_allsubjects(:,2:end,:),[],2:3));
    null_Msf_thr_neg = median(min(Msf_l3_allsubjects(:,2:end,:),[],2:3));
    null_Msb_thr_pos = median(max(Msb_l3_allsubjects(:,2:end,:),[],2:3));
    null_Msb_thr_neg = median(min(Msb_l3_allsubjects(:,2:end,:),[],2:3));

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
    yticks(-0.03:0.03:0.03);
    xticks([]);

    % Increase the font size of the tick labels
    ax = gca; % Get current axes                            
    ax.FontSize = 20;
    
    saveas(gcf, [outdir_fig '/sequenceness_length3.fig']);
    fig = openfig([outdir_fig '/sequenceness_length3.fig']);
    saveas(fig, [outdir_fig '/sequenceness_length3.jpeg']);
    system(['rm ' outdir_fig '/sequenceness_length3.fig']);


    %% plot l4

    null_Msf_thr_pos = median(max(Msf_l4_allsubjects(:,2:end,:),[],2:3));
    null_Msf_thr_neg = median(min(Msf_l4_allsubjects(:,2:end,:),[],2:3));
    null_Msb_thr_pos = median(max(Msb_l4_allsubjects(:,2:end,:),[],2:3));
    null_Msb_thr_neg = median(min(Msb_l4_allsubjects(:,2:end,:),[],2:3));

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
    yticks(-0.03:0.03:0.03);
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