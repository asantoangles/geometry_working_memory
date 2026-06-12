%% GERE project

folder_outputs = 'LB23';

% path
if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = [1:15];
else
    path_root = '/path_to_hpc';
end

path_results = [path_root '/github/results/source_geometry_lm/' folder_outputs];
addpath([path_root '/github/scripts/source_geometry_lm/utilities'])

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
performance = {'correct_trials'};
fun_i = @mean;
components = 1:3;

%%% merged time segments with event type

% event type: 1 = stim, 2 = delay
stim_labels  = ones(size(stim_segments,1),1);     % 1 = stim
delay_labels = 2 * ones(size(delay_segments,1),1); % 2 = delay

% combine segments
time_segments_all = [stim_segments, stim_labels; 
                     delay_segments, delay_labels];

time_segments_all_length3 = time_segments_all;
time_segments_all_length3(4,:) = [];
time_segments_all_length4 = time_segments_all;


% time_segments_all columns:
% column 1 = start, column 2 = end, column 3 = event type (1=stim, 2=delay)

%% subspace similarity metrics 

% loop over performance
for perf_i = 1:length(performance)

    disp(performance{perf_i});

    % loop over sequence length used
    for sequence_length = [3 4]

        if sequence_length == 3
            seq_name = 'length3';
            time_segments_all = time_segments_all_length3;
        elseif sequence_length == 4
            seq_name = 'length4';
            time_segments_all = time_segments_all_length4;
        end
        
        disp(' ');

        disp(seq_name)


        % outputs group
        frobenius_norm_ranks_ij_allsubjects = nan(size(time_segments_all,1), size(time_segments_all,1), sequence_length, sequence_length, length(subjects));
        procrustes_distance_ranks_ij_allsubjects = nan(size(time_segments_all,1), size(time_segments_all,1), sequence_length, sequence_length, length(subjects));
        % dimensions: time, time, rank i in distance_ij, rank j in distance_ij, subject

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
                
            outdir = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/subspaces_similarity'];
            
            if ~isfolder(outdir)
                mkdir(outdir);
            end
    
            % loop over time windows
            for seg_i = 1:size(time_segments_all, 1)
        
                seg_start = time_segments_all(seg_i, 1);
                seg_end   = time_segments_all(seg_i, 2);

                if time_segments_all(seg_i, 3) == 1 % stim
                    event = 'stim';
                else
                    event = 'delay';
                end

                event_i = event;

                indir = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                    event '_' num2str(seg_start) 'to' num2str(seg_end)];

                % set stim_order variable
                if time_segments_all(seg_i, 3) == 1 % stim

                    if sequence_length == 3
                        if seg_start == 1 && seg_end == 300
                            stim_order = 1;
                        elseif seg_start == 400 && seg_end == 700
                            stim_order = 2;
                        elseif seg_start == 800 && seg_end == 1100
                            stim_order = 3;
                        elseif seg_start == 1200 && seg_end == 1500
                            stim_order = 3;
                        end
                    elseif sequence_length == 4
                        if seg_start == 1 && seg_end == 300
                            stim_order = 1;
                        elseif seg_start == 400 && seg_end == 700
                            stim_order = 2;
                        elseif seg_start == 800 && seg_end == 1100
                            stim_order = 3;
                        elseif seg_start == 1200 && seg_end == 1500
                            stim_order = 4;
                        end
                    end                            

                else
                    stim_order = sequence_length;
                end

                stim_order_i = stim_order;

                % load pc_scores
                load([indir '/low_dim_space_' seq_name '.mat']); % low_dim_space

                pc_scores_i = low_dim_space.z_k; % (conditions by PCs)

                % loop over time windows
                for seg_ii = 1:size(time_segments_all, 1)

                    if seg_ii >= seg_i
            
                        seg_start = time_segments_all(seg_ii, 1);
                        seg_end   = time_segments_all(seg_ii, 2);
    
                        if time_segments_all(seg_ii, 3) == 1 % stim
                            event = 'stim';
                        else
                            event = 'delay';
                        end
                        
                        indir = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                            event '_' num2str(seg_start) 'to' num2str(seg_end)];
    
                        % set stim_order variable
                        if time_segments_all(seg_ii, 3) == 1 % stim
    
                            if sequence_length == 3
                                if seg_start == 1 && seg_end == 300
                                    stim_order = 1;
                                elseif seg_start == 400 && seg_end == 700
                                    stim_order = 2;
                                elseif seg_start == 800 && seg_end == 1100
                                    stim_order = 3;
                                elseif seg_start == 1200 && seg_end == 1500
                                    stim_order = 3;
                                end
                            elseif sequence_length == 4
                                if seg_start == 1 && seg_end == 300
                                    stim_order = 1;
                                elseif seg_start == 400 && seg_end == 700
                                    stim_order = 2;
                                elseif seg_start == 800 && seg_end == 1100
                                    stim_order = 3;
                                elseif seg_start == 1200 && seg_end == 1500
                                    stim_order = 4;
                                end
                            end                            
    
                        else
                            stim_order = sequence_length;
                        end
    
                        stim_order_ii = stim_order;
    
                        % load pc_scores
                        load([indir '/low_dim_space_' seq_name '.mat']); % low_dim_space
    
                        pc_scores_ii = low_dim_space.z_k; % (conditions by PCs)
    
                        %% compute subspace metrics for each rank
                        
                        for rank_i = 1:sequence_length

                            for rank_j = 1:sequence_length
    
                                % default to zero
                                dFro  = NaN;
                                dProc  = NaN;
                                                            
                                try
        
                                    start_row = (rank_i - 1)*number_locations + 1;
                                    end_row   = rank_i*number_locations;
                                    X1 = pc_scores_i(start_row:end_row, :);

                                    start_row = (rank_j - 1)*number_locations + 1;
                                    end_row   = rank_j*number_locations;
                                    X2 = pc_scores_ii(start_row:end_row, :);
                    
                                    [dFro, dProc] = subspace_metrics(X1, X2);
        
                                end
                    
                                %%% store group
                                frobenius_norm_ranks_ij_allsubjects(seg_i, seg_ii, rank_i, rank_j, sub_i) = dFro;
                                procrustes_distance_ranks_ij_allsubjects(seg_i, seg_ii, rank_i, rank_j, sub_i) = dProc;
                                frobenius_norm_ranks_ij_allsubjects(seg_ii, seg_i, rank_i, rank_j, sub_i) = dFro;
                                procrustes_distance_ranks_ij_allsubjects(seg_ii, seg_i, rank_i, rank_j, sub_i) = dProc;

                            end
    
                        end

                    end

                end

            end

        end

        %% raw distances - all subjects - all ranks
        
        outdir_group = [path_results '/subspaces_similarity/group_level'];
        if ~isfolder(outdir_group)
            mkdir(outdir_group);
        end

        save([outdir_group '/frobenius_norm_ranks_ij_allsubjects_' seq_name '.mat'], 'frobenius_norm_ranks_ij_allsubjects');
        save([outdir_group '/procrustes_distance_ranks_ij_allsubjects_' seq_name '.mat'], 'procrustes_distance_ranks_ij_allsubjects');
        


        %% average across ranks

        %%% distance between same rank
                
        frobenius_norm_ranks_ij_allsubjects_within = zeros(size(frobenius_norm_ranks_ij_allsubjects, 1), size(frobenius_norm_ranks_ij_allsubjects, 2), sequence_length, length(subjects));
        for r = 1:sequence_length
            frobenius_norm_ranks_ij_allsubjects_within(:,:,r,:) = ...
                squeeze(frobenius_norm_ranks_ij_allsubjects(:,:,r,r,:));
        end
        frobenius_norm_ranks_ij_allsubjects_within = squeeze(mean(frobenius_norm_ranks_ij_allsubjects_within, 3, 'omitnan'));
        save([outdir_group '/frobenius_norm_ranks_ij_allsubjects_within_' seq_name '.mat'], 'frobenius_norm_ranks_ij_allsubjects_within');



        procrustes_distance_ranks_ij_allsubjects_within = zeros(size(frobenius_norm_ranks_ij_allsubjects, 1), size(frobenius_norm_ranks_ij_allsubjects, 2), sequence_length, length(subjects));
        for r = 1:sequence_length
            procrustes_distance_ranks_ij_allsubjects_within(:,:,r,:) = ...
                squeeze(procrustes_distance_ranks_ij_allsubjects(:,:,r,r,:));
        end
        procrustes_distance_ranks_ij_allsubjects_within = squeeze(mean(procrustes_distance_ranks_ij_allsubjects_within, 3, 'omitnan'));
        save([outdir_group '/procrustes_distance_ranks_ij_allsubjects_within_' seq_name '.mat'], 'procrustes_distance_ranks_ij_allsubjects_within');

        %%% distance between different ranks
        
        num_cross = sequence_length^2 - sequence_length;   % 3^2 - 3 = 6
        frobenius_norm_ranks_ij_allsubjects_between = zeros(size(frobenius_norm_ranks_ij_allsubjects, 1), size(frobenius_norm_ranks_ij_allsubjects, 2), num_cross, length(subjects));
        pair_idx = 1;
        for i = 1:sequence_length
            for j = 1:sequence_length
                if i ~= j
                    frobenius_norm_ranks_ij_allsubjects_between(:,:,pair_idx,:) = ...
                        squeeze(frobenius_norm_ranks_ij_allsubjects(:,:,i,j,:));
                    pair_idx = pair_idx + 1;
                end
            end
        end
        frobenius_norm_ranks_ij_allsubjects_between = squeeze(mean(frobenius_norm_ranks_ij_allsubjects_between, 3, 'omitnan'));
        A_mean = mean(frobenius_norm_ranks_ij_allsubjects_between(~isnan(frobenius_norm_ranks_ij_allsubjects_between)));   % mean of non-NaN values
        frobenius_norm_ranks_ij_allsubjects_between(isnan(frobenius_norm_ranks_ij_allsubjects_between)) = A_mean;          % replace NaNs
        save([outdir_group '/frobenius_norm_ranks_ij_allsubjects_between_' seq_name '.mat'], 'frobenius_norm_ranks_ij_allsubjects_between');

        num_cross = sequence_length^2 - sequence_length;   % 3^2 - 3 = 6
        procrustes_distance_ranks_ij_allsubjects_between = zeros(size(frobenius_norm_ranks_ij_allsubjects, 1), size(frobenius_norm_ranks_ij_allsubjects, 2), num_cross, length(subjects));        pair_idx = 1;
        for i = 1:sequence_length
            for j = 1:sequence_length
                if i ~= j
                    procrustes_distance_ranks_ij_allsubjects_between(:,:,pair_idx,:) = ...
                        squeeze(procrustes_distance_ranks_ij_allsubjects(:,:,i,j,:));
                    pair_idx = pair_idx + 1;
                end
            end
        end
        procrustes_distance_ranks_ij_allsubjects_between = squeeze(mean(procrustes_distance_ranks_ij_allsubjects_between, 3, 'omitnan'));
        A_mean = mean(procrustes_distance_ranks_ij_allsubjects_between(~isnan(procrustes_distance_ranks_ij_allsubjects_between)));   % mean of non-NaN values
        procrustes_distance_ranks_ij_allsubjects_between(isnan(procrustes_distance_ranks_ij_allsubjects_between)) = A_mean;          % replace NaNs
        save([outdir_group '/procrustes_distance_ranks_ij_allsubjects_between_' seq_name '.mat'], 'procrustes_distance_ranks_ij_allsubjects_between');



        
        %% hierarchical models

        % first level: between > within (within-rank more stable than between-rank are more stable over time)
        %       here, the higher the t-values, the more stable within-rank

        % second level: delay > stim
        %       here, the higher the t-values, the more stable during the
        %       delay period

        % frobenius_norm_ranks_ij_allsubjects_within
        % procrustes_distance_ranks_ij_allsubjects_within
        % frobenius_norm_ranks_ij_allsubjects_between
        % procrustes_distance_ranks_ij_allsubjects_between

        for distance_i = 1:2

            disp(' ');

            if distance_i == 1

                X_within = frobenius_norm_ranks_ij_allsubjects_within;
                X_between = frobenius_norm_ranks_ij_allsubjects_between;
                disp('frobenius_norm')
                
            else


                X_within = procrustes_distance_ranks_ij_allsubjects_within;
                X_between = procrustes_distance_ranks_ij_allsubjects_between;
                disp('procrustes_distance')

            end

            [I, J, K] = size(X_within);
            T = nan(I, J);
            
            for i = 1:I
                for j = 1:J
                    
                    % vectorize
                    x = squeeze(X_between(i, j, :));   % data across k
                    y = squeeze(X_within(i, j, :));   % data across k
    
                    % omit na
                    x(isnan(x)) = mean(x, 'omitna');
                    y(isnan(y)) = mean(y, 'omitna');
    
                    % t-test
                    [h, p, ~, stats] = ttest(x, y);
                    
                    T(i,j) = stats.tstat;
                end
            end
    
            % set to zero
            T(1,1) = 0;
    
            %%% one-sample t-test (between > within) - all matrix
            disp(' ');
            disp('one-sample t-test (between > within) - all matrix');
            [h, p, ci, stats] = ttest(vec(T));
            fprintf('t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

%             %%% one-sample t-test (between > within) - on-diagonal
%             disp(' ');
%             disp('one-sample t-test (between > within) - on-diagonal');
%             [h, p, ci, stats] = ttest(diag(T));
%             fprintf('t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);
% 
%             %%% one-sample t-test (between > within) - upper triangle
%             disp(' ');
%             disp('one-sample t-test (between > within) - upper triangle');
%             A = vec(triu(T, 1));
%             A(A == 0) = [];
%             [h, p, ci, stats] = ttest(A);   % upper triangle, exclude diagonal);
%             fprintf('t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);
    
            %%% one-sample t-test (between > within) - stim
            disp(' ');
            disp('one-sample t-test (between > within) - stim');
            n = size(T,1);
            idx1 = 1:sequence_length;
            G_stim = T(idx1, idx1);
            G_stim = triu(G_stim, 1);   % upper triangle, exclude diagonal
            G_stim = G_stim(:);
            G_stim = G_stim(G_stim ~= 0);          % remove zeros
            [h, p, ci, stats] = ttest(G_stim);   % upper triangle, exclude diagonal);
            fprintf('t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            %%% one-sample t-test (between > within) - delay
            disp(' ');
            disp('one-sample t-test (between > within) - delay');
            n = size(T,1);
            idx2 = (sequence_length+1):n;
            G_delay = T(idx2, idx2);
            G_delay = triu(G_delay, 1);   % upper triangle
            G_delay = G_delay(:);
            G_delay = G_delay(G_delay ~= 0);          % remove zeros
            [h, p, ci, stats] = ttest(G_delay);   % upper triangle, exclude diagonal);
            fprintf('t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            %%% t-test delay > stim (between > within) - upper triangle
            disp(' ');
            disp('t-test delay > stim (between > within) - upper triangle')
            n = size(T,1);
            idx1 = 1:sequence_length;
            G_stim = T(idx1, idx1);
            G_stim = triu(G_stim, 1);   % upper triangle, exclude diagonal
            G_stim = G_stim(:);
            G_stim = G_stim(G_stim ~= 0);          % remove zeros
            idx2 = (sequence_length+1):n;
            G_delay = T(idx2, idx2);
            G_delay = triu(G_delay, 1);   % upper triangle
            G_delay = G_delay(:);
            G_delay = G_delay(G_delay ~= 0);          % remove zeros
            [h, p, ci, stats] = ttest2(G_delay, G_stim);
            fprintf('t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            %%% time consistency of subspaces during delay
            
            disp(' ');
            disp('time consistency of subspaces during delay');
            
            % indices for delay
            n = size(T,1);
            delay_idx = (sequence_length+1):n;
            num_delay = length(delay_idx);
            
            % extract delay–delay submatrix
            T_delay = T(delay_idx, delay_idx);
            
            % build vectors of:
            %   D = distances
            %   Δt = temporal separation
            D = [];
            Delta_t = [];
            
            for a = 1:num_delay
                for b = a+1:num_delay    % upper triangle only
                    D(end+1) = T_delay(a,b);
                    Delta_t(end+1) = abs(delay_idx(a) - delay_idx(b));
                end
            end
                                    
            % Regression
            X = [ones(length(Delta_t),1), Delta_t(:)];
            [b, ~, ~, ~, stats_reg] = regress(D(:), X);
            
            % Extract t-value for the slope term
            % stats_reg contains: [R2, F, p, error variance]
            % To get t, compute from F-statistic for the slope (1 df in numerator):
            t_slope = sqrt(stats_reg(2));   % since F(1,df) = t^2
            
            fprintf('t = %.4f, p = %.6f, \n', t_slope, stats_reg(3));

        end

    end

end

disp('Analysis done');
