%% GERE project

folder_outputs = 'LB23_decode_rank_highdim';

% path
if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = [1:15];
else
    path_root = '/path_to_hpc';
end

path_inputs = [path_root '/results/source_reconstruction'];
path_results = [path_root '/github/results/source_geometry_lm/' folder_outputs];
addpath([path_root '/github/scripts/source_geometry_lm/utilities'])

path_figures = [path_results '/figures'];
if ~isfolder(path_figures)
    mkdir(path_figures);
end

% settings
segments_start = 1:300:3701;
segments_end = 300:300:4000;
delay_segments = [segments_start; segments_end]';

stim_segments = [];
stim_segments(end+1,:) = [1 300];
stim_segments(end+1,:) = [400 700];
stim_segments(end+1,:) = [800 1100];
stim_segments(end+1,:) = [1200 1500];

number_parcels = 200;
sessions = 1:2;
number_locations = 8;
events = {'stim', 'delay'};
fun_i = @mean;
components = 1:3;
nFolds = 5;
nPerm  = 500; 

% loop over sequence length used
for sequence_length = [3 4]

    if sequence_length == 3
        seq_name = 'length3';
    elseif sequence_length == 4
        seq_name = 'length4';
    end

    disp(' ');
    disp('########')
    disp(seq_name);
    disp('########')

    %% =========================================================
    % LOAD DATA
    %% =========================================================

    number_iterations = 100;
    
    accuracy_all_subjects = cell(number_iterations, 1);
    accuracy_perm_all_subjects = cell(number_iterations, 1);
    
    % --------------------------------------------------
    % Load all iterations
    % --------------------------------------------------
    for iter_i = 1:number_iterations
    
        load([path_results '/decoding_results_' seq_name '_globaliter' num2str(iter_i) '.mat']);
    
        accuracy_all_subjects{iter_i} = decoding_results.accuracy_all_subjects;          % (sub × time)
        accuracy_perm_all_subjects{iter_i} = decoding_results.accuracy_perm_all_subjects; % (sub × time × perm)
    
    end
    
    % --------------------------------------------------
    % Concatenate
    % --------------------------------------------------
    acc_emp_all  = cat(3, accuracy_all_subjects{:});      % (sub × time × iter)
    acc_perm_all = cat(3, accuracy_perm_all_subjects{:}); % (sub × time × iter)
    
    [nSub, nTime, ~] = size(acc_emp_all);
    nPerm = size(acc_perm_all, 3);
    
    % --------------------------------------------------
    % Preallocate
    % --------------------------------------------------
    z = nan(1, nTime);
    p = nan(1, nTime);
    emp_mean = nan(1, nTime);
    mu_null = nan(1, nTime);
    sigma_null = nan(1, nTime);
    disp(' ');

    % --------------------------------------------------
    % Loop over time windows
    % --------------------------------------------------
    for t = 1:nTime
    
        % empirical: (sub × iter)
        emp_tmp = squeeze(acc_emp_all(:, t, :)); 
        
        % null: (sub × perm × iter)
        perm_tmp = squeeze(acc_perm_all(:, t, :, :)); 
    
        % vectorize
        emp_vec  = emp_tmp(:);
        perm_vec = perm_tmp(:);
    
        % remove NaNs
        emp_vec  = emp_vec(~isnan(emp_vec));
        perm_vec = perm_vec(~isnan(perm_vec));
    
        % compute stats
        emp_mean(t)   = mean(emp_vec);
        mu_null(t)    = mean(perm_vec);
        sigma_null(t) = std(perm_vec);
    
        z(t) = (emp_mean(t) - mu_null(t)) / sigma_null(t);
        p(t) = mean(perm_vec >= emp_mean(t));
    
    end

    disp(['mean(z) = ' num2str(mean(z))]);
    disp(['std(z) = ' num2str(std(z))]);
    disp(['min(z) = ' num2str(min(z))]);
    disp(['max(z) = ' num2str(max(z))]);


    %% =========================================================
    % PREPARE EMPIRICAL DATA
    %% =========================================================
    
    % acc_emp_all: (subjects × time × iterations)
    nIter = size(acc_emp_all, 3);
    
    % --------------------------------------------------
    % Mean across iterations (subject-level preserved)
    % --------------------------------------------------
    emp_mean_sub = squeeze(mean(acc_emp_all, 3, 'omitnan')); % (sub × time)
    
    mean_emp = mean(emp_mean_sub, 1, 'omitnan');
    
    % --------------------------------------------------
    % SEM: average of subject SEM across iterations
    % --------------------------------------------------
    sem_emp_iter = nan(nIter, length(time_axis));
    
    for it = 1:nIter
        
        tmp = squeeze(acc_emp_all(:,:,it)); % (sub × time)
        
        sem_emp_iter(it,:) = std(tmp, 0, 1, 'omitnan') ./ sqrt(size(tmp,1));
    end
    
    sem_emp = mean(sem_emp_iter, 1, 'omitnan');
    
    %% =========================================================
    % PERMUTATION (NULL) — time-resolved
    %% =========================================================
    
    % acc_perm_all: (subjects × time × perm × iter)
    
    nTime = size(acc_emp_all, 2);
    nIter = size(acc_perm_all, 4);
    
    mean_perm = nan(1, nTime);
    sem_perm  = nan(1, nTime);
    
    for t = 1:nTime
    
        % --------------------------------------------------
        % extract null for this time point
        % (subjects × perm × iter)
        % --------------------------------------------------
        tmp = squeeze(acc_perm_all(:, t, :, :));
    
        % vectorize across subjects, perm, iterations
        perm_vec = tmp(:);
        perm_vec = perm_vec(~isnan(perm_vec));
    
        % --------------------------------------------------
        % null statistics per time point
        % --------------------------------------------------
        mu_null    = mean(perm_vec);
        sigma_null = std(perm_vec);
    
        mean_perm(t) = mu_null;
        sem_perm(t)  = sigma_null / sqrt(length(perm_vec));
    
    end    

    %% =========================================================
    % MULTIPLE COMPARISONS (FDR CORRECTION)
    %% =========================================================
    
    adj_p = mafdr(p, 'BHFDR', true);
    sig_timepoints = adj_p < 0.05;
    
    %% =========================================================
    % FIGURE
    %% =========================================================
    
    figure; hold on;
    
    % --- empirical (red shaded) ---
    fill([time_axis fliplr(time_axis)], ...
         [mean_emp + sem_emp fliplr(mean_emp - sem_emp)], ...
         [1 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    plot(time_axis, mean_emp, 'r', 'LineWidth', 2);
    
    % --- permutation (blue shaded) ---
    fill([time_axis fliplr(time_axis)], ...
         [mean_perm + sem_perm fliplr(mean_perm - sem_perm)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    plot(time_axis, mean_perm, 'b', 'LineWidth', 2);
    
    % --- significant timepoints (black bar on top) ---
    ymax = max(mean_emp + sem_emp);
    y_sig = ymax + 0.01;
    
    sig = sig_timepoints(:)';
    
    % find start/end of contiguous significant segments
    d = diff([0 sig 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    
    y_sig = max(mean_emp + sem_emp) + 0.01;
    
    for i = 1:length(starts)
        
        x_start = time_axis(starts(i));
        x_end   = time_axis(ends(i));
        
        plot([x_start x_end], [y_sig y_sig], ...
            'k', 'LineWidth', 3);
    end
    
    %% =========================================================
    % FORMATTING
    %% =========================================================
    
    y_min = min(mean_perm - sem_perm);
    y_max = max(mean_emp + sem_emp);
    ylim([0.995 * y_min y_max * 1.005]);
    
    xlabel('Time (ms)');
    ylabel('Decoding accuracy');
    title(['Rank decoding - ' seq_name]);
    
    xticks(1:length(time_axis));
    xticklabels({'stim','150','450','750','1050','1350','1650','1950','2250','2550','2850','3150','3450','3750','4050'});
    xtickangle(45);
    
    xlim([1 length(time_axis)]);
    
    box on;

    %% =========================================================
    % SAVE
    %% =========================================================
    
    % File name
    outname_png = [path_figures '/figure_decoding_rank_' seq_name '.png'];
    
    % Make renderer clean for vector export
    set(gcf, 'Renderer', 'painters');
    
    % Save high-resolution PNG
    exportgraphics(gcf, outname_png, 'Resolution', 300);
    
end

disp('Analysis done');
