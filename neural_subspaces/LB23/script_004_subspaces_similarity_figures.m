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

        %% load distances

        outdir_group = [path_results '/subspaces_similarity/group_level'];

        %%% distance between same rank
        load([outdir_group '/frobenius_norm_ranks_ij_allsubjects_within_' seq_name '.mat']);        % frobenius_norm_ranks_ij_allsubjects_within
        load([outdir_group '/procrustes_distance_ranks_ij_allsubjects_within_' seq_name '.mat']);   % procrustes_distance_ranks_ij_allsubjects_within

        %%% distance between different ranks
        load([outdir_group '/frobenius_norm_ranks_ij_allsubjects_between_' seq_name '.mat']);       % frobenius_norm_ranks_ij_allsubjects_between
        load([outdir_group '/procrustes_distance_ranks_ij_allsubjects_between_' seq_name '.mat']);  % procrustes_distance_ranks_ij_allsubjects_between

        %% average across subjects
        frobenius_norm_ranks_ij_allsubjects_within = mean(frobenius_norm_ranks_ij_allsubjects_within, 3);
        procrustes_distance_ranks_ij_allsubjects_within = mean(procrustes_distance_ranks_ij_allsubjects_within, 3);
        frobenius_norm_ranks_ij_allsubjects_between = mean(frobenius_norm_ranks_ij_allsubjects_between, 3);
        procrustes_distance_ranks_ij_allsubjects_between = mean(procrustes_distance_ranks_ij_allsubjects_between, 3);

        %% stack between and within in one matrix
    
        frob_matrix = nan(size(frobenius_norm_ranks_ij_allsubjects_within));
        proc_matrix = nan(size(frobenius_norm_ranks_ij_allsubjects_within));

        for i = 1:size(frob_matrix, 1)
            for j = 1:size(frob_matrix, 2)

                if i < j
                    frob_matrix(i, j) = frobenius_norm_ranks_ij_allsubjects_within(i, j);
                    proc_matrix(i, j) = procrustes_distance_ranks_ij_allsubjects_within(i, j);
                else
                    frob_matrix(i, j) = frobenius_norm_ranks_ij_allsubjects_between(i, j);
                    proc_matrix(i, j) = procrustes_distance_ranks_ij_allsubjects_between(i, j);
                end


            end
        end

        %% plot

        outdir_figures = [path_results '/subspaces_similarity/figures'];
        if ~isfolder(outdir_figures)
            mkdir(outdir_figures);
        end

        for fig_i = 1:2

            if fig_i == 1

                % fig_title = 'Frobenius norm - upper (within-rank), lower (between-ranks)';
                M = frob_matrix;
                filename = [outdir_figures '/figure_frobenius_' seq_name '.png'];

            else

                % fig_title = 'Procrustes distance - upper (within-rank), lower (between-ranks)';
                M = proc_matrix;
                filename = [outdir_figures '/figure_procrustes_' seq_name '.png'];

            end
            
            % Plot heatmap
            figure;
            imagesc(M);
            colorbar;
            
            % Ensure correct axis orientation (important!)
            axis xy;   % makes y-axis increase upward (so bottom-left is (1,1))
            
            % color map
            colormap(turbo);
            
            % Color limits
            caxis([min(M(:)) max(M(:))]);
            
            hold on;
            
            % Draw separating lines (adjusted for flipped matrix)
            if sequence_length == 3
                xline(3.5, 'k', 'LineWidth', 2);
                yline(3.5, 'k', 'LineWidth', 2);
            else
                xline(4.5, 'k', 'LineWidth', 2);
                yline(4.5, 'k', 'LineWidth', 2);
            end
            
            hold off;
            
            % Axis labels
            xlabel('Time (ms)');
            ylabel('Time (ms)');
            
            % -------------------------
            % Custom tick labels
            % -------------------------
            
            n = size(M,1);
                    
            tick_vals = strings(1, n);
            
            tick_vals(1:sequence_length) = "s" + string(1:sequence_length);
            
            for i = (sequence_length+1):n
                val = 0.150 + (i - (sequence_length+1)) * 0.300;
                tick_vals(i) = sprintf('%.3f', val);
            end
    
            % From sequence_length+1 onward: 0.150, 0.450, 0.750, ...
            for i = (sequence_length+1):n
                tick_vals(i) = 0.150 + (i - (sequence_length+1)) * 0.300;
            end       
            
            % Show only indices 5 and 10
            if sequence_length == 3
                idx = [1 2 3 5 7 9 11 13 15];
            else
                idx = [1 2 3 4 6 8 10 12 14 16];
            end
            
            xticks(idx);
            yticks(idx);
            
            xticklabels(string(tick_vals(idx)));
            yticklabels(string(tick_vals(idx)));

            % save
            exportgraphics(gcf, filename, 'ContentType', 'vector');

        end

        
    end

end

disp('Analysis done');
