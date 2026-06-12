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

path_figures = [path_results '/rank_specificity'];
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
performance = {'correct_trials'};
fun_i = @mean;
components = 1:3;
events = {'stim', 'delay'};

%% procrustes/frobenius distance - within time windows

% loop over sequence length used
for sequence_length = [3 4]

    if sequence_length == 3
        seq_name = 'length3';
        seq_name_fig = 'Length 3';
    elseif sequence_length == 4
        seq_name = 'length4';
        seq_name_fig = 'Length 4';
    end
    
    %% procrustes_distance_ranks_ij_allsubjects
    
    outdir_group = [path_results '/subspaces_similarity/group_level'];
    if ~isfolder(outdir_group)
        mkdir(outdir_group);
    end

    load([outdir_group '/frobenius_norm_ranks_ij_allsubjects_' seq_name '.mat']);       % frobenius_norm_ranks_ij_allsubjects
    load([outdir_group '/procrustes_distance_ranks_ij_allsubjects_' seq_name '.mat']);  % procrustes_distance_ranks_ij_allsubjects
    % dimensions: time, time, rank i in distance_ij, rank j in distance_ij, subject

    for distance_metric_i = 1:2

        if distance_metric_i == 1
            distance_tensor = frobenius_norm_ranks_ij_allsubjects;
            name_yaxis = 'Frobenius norm';
            filename_figure = 'frobenius';
        else
            distance_tensor = procrustes_distance_ranks_ij_allsubjects;
            name_yaxis = 'Procrustes distance';
            filename_figure = 'procrustes';
        end

        % Dimensions
        [nT1, nT2, nRi, nRj, nSub] = size(distance_tensor);
        % time, time, rank, rank, subject
        
        % --------------------------------------------------
        % 1) Get upper-triangular indices (time1 < time2)
        % --------------------------------------------------
        mask_diag = logical(eye(nT1, nT2)); % only diagonal (t1 = t2)
        nPairs = sum(mask_diag(:));     

        % --------------------------------------------------
        % 2) Extract upper triangle and vectorize
        % Result: samples × rank_i × rank_j
        % samples = (time pairs) × subject
        % --------------------------------------------------
        dist_vec = nan(nPairs*nSub, nRi, nRj);
        
        cnt = 0;
        
        for sub = 1:nSub
            
            % Extract this subject: (time × time × rank_i × rank_j)
            tmp = distance_tensor(:,:,:,:,sub);
            
            % Loop over rank dims and extract upper triangle
            tmp_vec = nan(nPairs, nRi, nRj);
            
            for i = 1:nRi
                for j = 1:nRj
                    
                    mat = tmp(:,:,i,j); % time × time
                    tmp_vec(:,i,j) = mat(mask_diag);
                    
                end
            end
            
            % store
            dist_vec(cnt + (1:nPairs),:,:) = tmp_vec;
            cnt = cnt + nPairs;
        end
        
        nSamples = size(dist_vec,1);
        
        % --------------------------------------------------
        % 3) Subject index
        % --------------------------------------------------
        sub_ids = repelem((1:nSub)', nPairs);
    
        % --------------------------------------------------
        % 4) Extract distances for all rank separations
        % --------------------------------------------------
        rank_dists = 1:(nRi-1);
        nD = length(rank_dists);
        
        vals = nan(nSamples, nD);
        
        for di = 1:nD
            d = rank_dists(di);
            
            for s = 1:nSamples
                
                tmp = [];
                
                for i = 1:nRi
                    for j = 1:nRj
                        
                        if i < j && abs(i - j) == d
                            tmp = [tmp; dist_vec(s,i,j)];
                        end
                        
                    end
                end
                
                vals(s,di) = mean(tmp, 'omitnan');
                
            end
        end
    
        % --------------------------------------------------
        % 5) Long format table
        % --------------------------------------------------
        Y = vals(:);
        
        RankDist = repelem(rank_dists', nSamples);
        Subject  = repmat(sub_ids, nD, 1);
        
        tbl = table(Y, RankDist, Subject);
        
        tbl.Subject = categorical(tbl.Subject);
        tbl.RankDist = categorical(tbl.RankDist);
    
        % --------------------------------------------------
        % 6) Linear model (subject as covariate)
        % --------------------------------------------------
        mdl = fitlm(tbl, 'Y ~ RankDist + Subject');
            
        disp(mdl);
        
        % Extract p-value for rank separation effect
        p = mdl.Coefficients.pValue(2); % RankDist effect
        tval = mdl.Coefficients.tStat(2);
        
        fprintf('Effect of rank separation: t = %.3f, p = %.4e\n', tval, p);
        
        % --------------------------------------------------
        % 7) Figure
        % --------------------------------------------------
        figure('Color','w','Position',[100 100 600 450]);
        hold on;
        
        fs = 18;
        ms = 40;
        lw_box = 1.8;
        lw_line = 1.5;
        
        % -------------------------
        % Boxplot
        % -------------------------
        boxplot(vals, ...
            'Labels', arrayfun(@(x) sprintf('d = %d',x), rank_dists, 'UniformOutput', false), ...
            'Widths', 0.5, ...
            'Symbol','');
        
        set(gca, 'LineWidth', 1.5, ...
                 'FontSize', fs, ...
                 'FontName', 'Arial');
        
        h = findobj(gca,'Tag','Box');
        for j = 1:length(h)
            patch(get(h(j),'XData'), get(h(j),'YData'), ...
                [0.8 0.8 0.8], ...
                'FaceAlpha', 0.4, ...
                'EdgeColor', 'k', ...
                'LineWidth', lw_box);
        end
        
        % -------------------------
        % Light samples
        % -------------------------
        nPlot = min(4000, nSamples);
        idx = randperm(nSamples, nPlot);
        
        for di = 1:nD
            scatter(ones(nPlot,1)*di, vals(idx,di), 12, ...
                'k', 'filled', 'MarkerFaceAlpha', 0.08);
        end
        
        % -------------------------
        % Subject means
        % -------------------------
        vals_sub = nan(nSub, nD);
        
        for s = 1:nSub
            idx_s = sub_ids == s;
            
            for di = 1:nD
                vals_sub(s,di) = mean(vals(idx_s,di), 'omitnan');
            end
        end
        
        % plot subject means + lines
        for s = 1:nSub
            scatter(1:nD, vals_sub(s,:), ms, ...
                'filled', 'MarkerFaceColor',[0 0.45 0.74]);
            
            plot(1:nD, vals_sub(s,:), '-', ...
                'Color', [0.6 0.6 0.6], ...
                'LineWidth', lw_line);
        end
        
        % axis
        if distance_metric_i == 1
            ylim([0.8 1.8])
        else
            ylim([0.4 1])
        end
        xlim([0.5 nD + 0.5]);
        xticks(1:nD);
        
        xticklabels(arrayfun(@(x) sprintf('Rank separation = %d', x), ...
            rank_dists, 'UniformOutput', false));
        
        set(gca, 'TickLength', [0 0]); % cleaner look   
            
        % -------------------------
        % Labels
        % -------------------------
        ylabel(name_yaxis, 'FontSize', fs+2);
        
        xticklabels(arrayfun(@(x) sprintf('Δ rank = %d', x), rank_dists, 'UniformOutput', false));
        
        if p > 0.05
            ttl = sprintf('t = %.1f, p > 0.05', round(tval, 2));
        elseif p < 0.0001
            ttl = sprintf('t = %.1f, p < 0.0001', round(tval, 2));
        else
            ttl = sprintf('t = %.1f, p = %.4f', round(tval, 2), p);
        end
        
        title(ttl, 'FontSize', fs+2, 'FontWeight', 'bold');

        grid on;
        box on;

        % -------------------------
        % Save
        % -------------------------

        % Set high-resolution export settings
        set(gcf, 'Renderer', 'painters');
                
        % Save high-resolution PNG
        exportgraphics(gcf, [path_figures '/figure_withintime_' filename_figure '_length' num2str(sequence_length) '.png'], 'Resolution', 300);

    end

end



%% procrustes/frobenius distance - between time windows

if 1 == 2

    % loop over sequence length used
    for sequence_length = [3 4]
    
        if sequence_length == 3
            seq_name = 'length3';
            seq_name_fig = 'Length 3';
        elseif sequence_length == 4
            seq_name = 'length4';
            seq_name_fig = 'Length 4';
        end
        
        %% procrustes_distance_ranks_ij_allsubjects
        
        outdir_group = [path_results '/subspaces_similarity/group_level'];
        if ~isfolder(outdir_group)
            mkdir(outdir_group);
        end
    
        load([outdir_group '/frobenius_norm_ranks_ij_allsubjects_' seq_name '.mat']);       % frobenius_norm_ranks_ij_allsubjects
        load([outdir_group '/procrustes_distance_ranks_ij_allsubjects_' seq_name '.mat']);  % procrustes_distance_ranks_ij_allsubjects
        % dimensions: time, time, rank i in distance_ij, rank j in distance_ij, subject
    
        for distance_metric_i = 1:2
    
            if distance_metric_i == 1
                distance_tensor = frobenius_norm_ranks_ij_allsubjects;
                name_yaxis = 'Frobenius norm';
                filename_figure = 'frobenius';
            else
                distance_tensor = procrustes_distance_ranks_ij_allsubjects;
                name_yaxis = 'Procrustes distance';
                filename_figure = 'procrustes';
            end
    
            % Dimensions
            [nT1, nT2, nRi, nRj, nSub] = size(distance_tensor);
            
            % --------------------------------------------------
            % 1) Get upper-triangular indices (time1 < time2)
            % --------------------------------------------------
            mask_ut = triu(true(nT1, nT2), 1); % exclude diagonal
            nPairs = sum(mask_ut(:));
            
            % --------------------------------------------------
            % 2) Extract upper triangle and vectorize
            % Result: samples × rank_i × rank_j
            % samples = (time pairs) × subject
            % --------------------------------------------------
            dist_vec = nan(nPairs*nSub, nRi, nRj);
            
            cnt = 0;
            
            for sub = 1:nSub
                
                % Extract this subject: (time × time × rank_i × rank_j)
                tmp = distance_tensor(:,:,:,:,sub);
                
                % Loop over rank dims and extract upper triangle
                tmp_vec = nan(nPairs, nRi, nRj);
                
                for i = 1:nRi
                    for j = 1:nRj
                        
                        mat = tmp(:,:,i,j); % time × time
                        tmp_vec(:,i,j) = mat(mask_ut);
                        
                    end
                end
                
                % store
                dist_vec(cnt + (1:nPairs),:,:) = tmp_vec;
                cnt = cnt + nPairs;
            end
            
            nSamples = size(dist_vec,1);
            
            % --------------------------------------------------
            % 3) Subject index
            % --------------------------------------------------
            sub_ids = repelem((1:nSub)', nPairs);
        
            % --------------------------------------------------
            % 4) Extract distances for all rank separations
            % --------------------------------------------------
            rank_dists = 1:(nRi-1);
            nD = length(rank_dists);
            
            vals = nan(nSamples, nD);
            
            for di = 1:nD
                d = rank_dists(di);
                
                for s = 1:nSamples
                    
                    tmp = [];
                    
                    for i = 1:nRi
                        for j = 1:nRj
                            
                            if i < j && abs(i - j) == d
                                tmp = [tmp; dist_vec(s,i,j)];
                            end
                            
                        end
                    end
                    
                    vals(s,di) = mean(tmp, 'omitnan');
                    
                end
            end
        
            % --------------------------------------------------
            % 5) Long format table
            % --------------------------------------------------
            Y = vals(:);
            
            RankDist = repelem(rank_dists', nSamples);
            Subject  = repmat(sub_ids, nD, 1);
            
            tbl = table(Y, RankDist, Subject);
            
            tbl.Subject = categorical(tbl.Subject);
            tbl.RankDist = categorical(tbl.RankDist);
        
            % --------------------------------------------------
            % 6) Linear model (subject as covariate)
            % --------------------------------------------------
            mdl = fitlm(tbl, 'Y ~ RankDist + Subject');
                
            disp(mdl);
            
            % Extract p-value for rank separation effect
            p = mdl.Coefficients.pValue(2); % RankDist effect
            tval = mdl.Coefficients.tStat(2);
            
            fprintf('Effect of rank separation: t = %.3f, p = %.4e\n', tval, p);
            
            % --------------------------------------------------
            % 7) Figure
            % --------------------------------------------------
            figure('Color','w','Position',[100 100 600 450]);
            hold on;
            
            fs = 18;
            ms = 40;
            lw_box = 1.8;
            lw_line = 1.5;
            
            % -------------------------
            % Boxplot
            % -------------------------
            boxplot(vals, ...
                'Labels', arrayfun(@(x) sprintf('d = %d',x), rank_dists, 'UniformOutput', false), ...
                'Widths', 0.5, ...
                'Symbol','');
            
            set(gca, 'LineWidth', 1.5, ...
                     'FontSize', fs, ...
                     'FontName', 'Arial');
            
            h = findobj(gca,'Tag','Box');
            for j = 1:length(h)
                patch(get(h(j),'XData'), get(h(j),'YData'), ...
                    [0.8 0.8 0.8], ...
                    'FaceAlpha', 0.4, ...
                    'EdgeColor', 'k', ...
                    'LineWidth', lw_box);
            end
            
            % -------------------------
            % Light samples
            % -------------------------
            nPlot = min(4000, nSamples);
            idx = randperm(nSamples, nPlot);
            
            for di = 1:nD
                scatter(ones(nPlot,1)*di, vals(idx,di), 12, ...
                    'k', 'filled', 'MarkerFaceAlpha', 0.08);
            end
            
            % -------------------------
            % Subject means
            % -------------------------
            vals_sub = nan(nSub, nD);
            
            for s = 1:nSub
                idx_s = sub_ids == s;
                
                for di = 1:nD
                    vals_sub(s,di) = mean(vals(idx_s,di), 'omitnan');
                end
            end
            
            % plot subject means + lines
            for s = 1:nSub
                scatter(1:nD, vals_sub(s,:), ms, ...
                    'filled', 'MarkerFaceColor',[0 0.45 0.74]);
                
                plot(1:nD, vals_sub(s,:), '-', ...
                    'Color', [0.6 0.6 0.6], ...
                    'LineWidth', lw_line);
            end
            
            % axis
            ylim([0 Inf])
            xlim([0.5 nD + 0.5]);
            xticks(1:nD);
            
            xticklabels(arrayfun(@(x) sprintf('Rank separation = %d', x), ...
                rank_dists, 'UniformOutput', false));
            
            set(gca, 'TickLength', [0 0]); % cleaner look   
                
            % -------------------------
            % Labels
            % -------------------------
            ylabel(name_yaxis, 'FontSize', fs+2);
            
            xticklabels(arrayfun(@(x) sprintf('Δ rank = %d', x), rank_dists, 'UniformOutput', false));
            
            if p > 0.05
                ttl = sprintf('t = %.1f, p > 0.05', round(tval, 2));
            elseif p < 0.0001
                ttl = sprintf('t = %.1f, p < 0.0001', round(tval, 2));
            else
                ttl = sprintf('t = %.1f, p = %.4f', round(tval, 2), p);
            end
            
            title(ttl, 'FontSize', fs+2, 'FontWeight', 'bold');
    
            grid on;
            box on;
    
            % -------------------------
            % Save
            % -------------------------
    
            % Set high-resolution export settings
            set(gcf, 'Renderer', 'painters');
                    
            % Save high-resolution PNG
            exportgraphics(gcf, [path_figures '/figure_betweentime_' filename_figure '_length' num2str(sequence_length) '.png'], 'Resolution', 300);
    
        end
    
    end

end

%% PA/VAF

% --------------------------------------------------
% LOOP OVER ALIGNMENT METRICS
% --------------------------------------------------
for alignment_metric = 1:2

    if alignment_metric == 1
        name_yaxis = 'PA';
        filename_figure = 'angle';
    else
        name_yaxis = 'VAF';
        filename_figure = 'vaf';
    end

    % --------------------------------------------------
    % LOOP OVER SEQUENCE LENGTH
    % --------------------------------------------------
    for sequence_length = [3 4]
    
        if sequence_length == 3
            seq_name = 'length3';
        elseif sequence_length == 4
            seq_name = 'length4';
        end
    
        % --------------------------------------------------
        % INITIALIZE CONCATENATION
        % (subjects × rank_dist × segments)
        % --------------------------------------------------
        vals_all = [];
        
        % --------------------------------------------------
        % LOOP OVER EVENTS
        % --------------------------------------------------
        for event_i = 1:length(events)
        
            if strcmp(events{event_i}, 'delay')
                time_segments = delay_segments;
            elseif strcmp(events{event_i}, 'stim')
                time_segments = stim_segments(sequence_length,:);
            end
    
            % --------------------------------------------------
            % LOOP OVER TIME WINDOWS (NOW ONLY FOR COLLECTION)
            % --------------------------------------------------
            for seg_i = 1:size(time_segments, 1)
            
                seg_start = time_segments(seg_i, 1);
                seg_end   = time_segments(seg_i, 2);
    
                vals_sub = [];
    
                % --------------------------------------------------
                % LOOP OVER SUBJECTS
                % --------------------------------------------------
                for sub_i = 1:length(subjects)
                
                    subject = subjects(sub_i);
                
                    if subject < 10
                        subject_ID = ['sub_0' num2str(subject)];
                    else
                        subject_ID = ['sub_' num2str(subject)];
                    end
                
                    indir = [path_results '/' subject_ID '/' func2str(fun_i) '/correct_trials/' ...
                        events{event_i} '_' num2str(seg_start) 'to' num2str(seg_end)];
    
                    load([indir '/' filename_figure '_' seq_name '.mat']); % or vaf

                    if strcmp(filename_figure, 'angle')
                        alignment_matrix = angle;
                    else
                        alignment_matrix = vaf;
                    end
        
                    [nRi, nRj] = size(alignment_matrix);
    
                    rank_dists = 1:(nRi-1);
                    nD = length(rank_dists);
    
                    vals_tmp = nan(1, nD);
    
                    for di = 1:nD
                        d = rank_dists(di);
    
                        vals = [];
    
                        for i = 1:nRi
                            for j = 1:nRj
                                if i < j && abs(i - j) == d
                                    vals = [vals; alignment_matrix(i,j)];
                                end
                            end
                        end
    
                        vals_tmp(di) = mean(vals, 'omitnan');
                    end
    
                    vals_sub(sub_i, :) = vals_tmp;
    
                end % subjects
    
                % store: (subjects × rank_dist × segment)
                vals_all(:,:,seg_i) = vals_sub;
    
            end % segments
        end % events
    
        % --------------------------------------------------
        % CONCATENATE ACROSS SEGMENTS
        % → treat segments as repeated measures
        % --------------------------------------------------
        [nSub, nD, nSeg] = size(vals_all);
    
        % reshape: (subject × segment) → samples
        vals_concat = reshape(vals_all, nSub*nSeg, nD);
    
        % subject index (repeat per segment)
        sub_ids = repelem((1:nSub)', nSeg);
    
        % --------------------------------------------------
        % STATS
        % --------------------------------------------------
        Y = vals_concat(:);
    
        RankDist = repelem((1:nD)', nSub*nSeg);
        Subject  = repmat(sub_ids, nD, 1);
    
        tbl = table(Y, RankDist, Subject);
    
        tbl.Subject = categorical(tbl.Subject);
        tbl.RankDist = categorical(tbl.RankDist);
    
        mdl = fitlm(tbl, 'Y ~ RankDist + Subject');
    
        p = mdl.Coefficients.pValue(2);
        tval = mdl.Coefficients.tStat(2);
    
        % --------------------------------------------------
        % PLOT
        % --------------------------------------------------
        figure('Color','w','Position',[100 100 600 450]);
        hold on;
        
        fs = 18;
        ms = 40;
        lw_box = 1.8;
        lw_line = 1.5;
        
        % -------------------------
        % Boxplot
        % -------------------------
        boxplot(vals_concat, ...
            'Labels', arrayfun(@(x) sprintf('d = %d',x), rank_dists, 'UniformOutput', false), ...
            'Widths', 0.5, ...
            'Symbol','');
        
        set(gca, 'LineWidth', 1.5, ...
                 'FontSize', fs, ...
                 'FontName', 'Arial');
        
        % nicer boxes
        h = findobj(gca,'Tag','Box');
        for j = 1:length(h)
            patch(get(h(j),'XData'), get(h(j),'YData'), ...
                [0.8 0.8 0.8], ...
                'FaceAlpha', 0.4, ...
                'EdgeColor', 'k', ...
                'LineWidth', lw_box);
        end
        
        % -------------------------
        % Light samples (all segments)
        % -------------------------
        nSamples = size(vals_concat,1);
        nPlot = min(4000, nSamples);
        idx = randperm(nSamples, nPlot);
        
        for di = 1:nD
            scatter(ones(nPlot,1)*di, vals_concat(idx,di), 12, ...
                'k', 'filled', 'MarkerFaceAlpha', 0.08);
        end
        
        % -------------------------
        % Subject means (across segments)
        % -------------------------
        vals_sub_mean = squeeze(mean(vals_all, 3, 'omitnan')); % (subjects × rank_dist)
        
        for s = 1:nSub
            scatter(1:nD, vals_sub_mean(s,:), ms, ...
                'filled', 'MarkerFaceColor',[0 0.45 0.74]);
        
            plot(1:nD, vals_sub_mean(s,:), '-', ...
                'Color', [0.6 0.6 0.6], ...
                'LineWidth', lw_line);
        end
        
        % -------------------------
        % Axis formatting (FIXED)
        % -------------------------
        if alignment_metric == 1
            ylim([10 100])
        else
            ylim([0.2 1])
        end
        xlim([0.5 nD + 0.5]);
        xticks(1:nD);
        
        xticklabels(arrayfun(@(x) sprintf('Δ rank = %d', x), ...
            rank_dists, 'UniformOutput', false));
        
        set(gca, 'TickLength', [0 0]);
        
        % -------------------------
        % Labels
        % -------------------------
        ylabel(name_yaxis, 'FontSize', fs+2);
        
        if p > 0.05
            ttl = sprintf('t = %.1f, p > 0.05', round(tval, 2));
        elseif p < 0.0001
            ttl = sprintf('t = %.1f, p < 0.0001', round(tval, 2));
        else
            ttl = sprintf('t = %.1f, p = %.4f', round(tval, 2), p);
        end
        
        title(ttl, 'FontSize', fs+2, 'FontWeight', 'bold');

        grid on;
        box on;    
        
        % -------------------------
        % Save
        % -------------------------
    
        % Set high-resolution export settings
        set(gcf, 'Renderer', 'painters');
                
        % Save high-resolution PNG
        exportgraphics(gcf, [path_figures '/figure_' filename_figure '_length' num2str(sequence_length) '.png'], 'Resolution', 300);
    
    end % sequence length


end

close all

disp('Analysis done');
