%% GERE project

folder_inputs = 'X_matrix';
folder_outputs = 'LB25';

% path
if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = [1:15];
else
    path_root = '/path_to_hpc';
end

path_inputs = [path_root '/github/results/source_geometry_lm/' folder_inputs];
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
performance = {'correct_trials', 'incorrect_trials'};
events = {'delay', 'stim'};
fun_i = @mean;
components = 1:3;
iterations_bootstrap = 1000;

%% geometric analysis of neural subspaces 

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
    
    disp(subject_ID);

    % loop over evens
    for event_i = 1:length(events)
    
        if strcmp(events{event_i}, 'delay')
            time_segments = delay_segments;
        elseif strcmp(events{event_i}, 'stim')
            time_segments = stim_segments;
        end
            
        % loop over time windows
        for seg_i = 1:size(time_segments, 1)
    
            seg_start = time_segments(seg_i, 1);
            seg_end   = time_segments(seg_i, 2);

            % loop over performance
            for perf_i = 1:length(performance)

                % loop over sequence length used
                for sequence_length = [3 4]

                    if sequence_length == 3
                        seq_name = 'length3';
                    elseif sequence_length == 4
                        seq_name = 'length4';
                    end

                    indir = [path_inputs '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                        events{event_i} '_' num2str(seg_start) 'to' num2str(seg_end)];

                    indir_correct = [path_inputs '/' subject_ID '/' func2str(fun_i) '/correct_trials/delay_1801to2100'];
                    
                    outdir = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                        events{event_i} '_' num2str(seg_start) 'to' num2str(seg_end)];
                    
                    if ~isfolder(outdir)
                        mkdir(outdir);
                    end

                    % load X matrix for correct trials (irrespective of perf_i value)
                    % used to compute eigvecs_correct
                    X_correct_trials = load([indir_correct '/X_matrix_' seq_name '.mat']);
                    X_correct_trials = X_correct_trials.X;

                    % load X matrix for correct or incorrect trials
                    % depending on the value of perf_i
                    % used to compute pc_scores by projecting X to the eigvecs_correct
                    load([indir '/X_matrix_' seq_name '.mat']); % X

                    % set stim_order variable
                    if strcmp(events{event_i}, 'stim')

                        if sequence_length == 3
                            if time_segments(seg_i,1) == 1 && time_segments(seg_i,2) == 300
                                stim_order = 1;
                            elseif time_segments(seg_i,1) == 400 && time_segments(seg_i,2) == 700
                                stim_order = 2;
                            elseif time_segments(seg_i,1) == 800 && time_segments(seg_i,2) == 1100
                                stim_order = 3;
                            elseif time_segments(seg_i,1) == 1200 && time_segments(seg_i,2) == 1500
                                stim_order = 3;
                            end
                        elseif sequence_length == 4
                            if time_segments(seg_i,1) == 1 && time_segments(seg_i,2) == 300
                                stim_order = 1;
                            elseif time_segments(seg_i,1) == 400 && time_segments(seg_i,2) == 700
                                stim_order = 2;
                            elseif time_segments(seg_i,1) == 800 && time_segments(seg_i,2) == 1100
                                stim_order = 3;
                            elseif time_segments(seg_i,1) == 1200 && time_segments(seg_i,2) == 1500
                                stim_order = 4;
                            end
                        end                            

                    else
                        stim_order = sequence_length;
                    end

                    if ~isempty(X)

                        %% compute eigvecs
    
                        if ~isfile([outdir '/eigvecs_' seq_name '.mat'])
    
                            eigvecs = compute_eigvecs(X);
    
                            save([outdir '/eigvecs_' seq_name '.mat'], 'eigvecs');
    
                        end
    
                        %% compute PCA, best-fitting planes, principal angles and variance metrics
        
                        if ~isfile([outdir '/vaf_' seq_name '.mat'])
                                                    
                            % pca, best-fitting planes and angles - pca
                            % on all trials data, and PC scores (data
                            % projection on PCs) done with perforance
                            % specific X matrix
    
                            if strcmp(events{event_i}, 'delay')
                                
                                [z_k, z_k_exp, planes] = compute_subspaces(X_correct_trials, X, components, number_locations);
    
                            elseif strcmp(events{event_i}, 'stim')
    
                                [z_k, z_k_exp, planes] = compute_subspaces_stim(X_correct_trials, X, components, stim_order, number_locations);
    
                            end
    
                            % compute angle and vaf
                            [angle, vaf] = compute_alignment_metrics(planes, stim_order);
    
                            low_dim_space = [];
                            low_dim_space.z_k = z_k;
                            low_dim_space.plane_r1 = planes.plane_r1;
    
                            if stim_order > 1
                                low_dim_space.plane_r2 = planes.plane_r2;
                            end
    
                            if stim_order > 2
                                low_dim_space.plane_r3 = planes.plane_r3;
                            end
    
                            if stim_order > 3
                                low_dim_space.plane_r4 = planes.plane_r4;
                            end
        
                            save([outdir '/low_dim_space_' seq_name '.mat'], 'low_dim_space');
                            save([outdir '/angle_' seq_name '.mat'], 'angle');
                            save([outdir '/vaf_' seq_name '.mat'], 'vaf');
    
                            variance_explained = [];
                            variance_explained.z_k = z_k_exp;
                            variance_explained.explained_r1 = planes.explained_r1;
    
                            if stim_order > 1
                                variance_explained.explained_r2 = planes.explained_r2;
                            end
    
                            if stim_order > 2
                                variance_explained.explained_r3 = planes.explained_r3;
                            end
    
                            if stim_order > 3
                                variance_explained.explained_r4 = planes.explained_r4;
                            end
    
                            save([outdir '/variance_explained_' seq_name '.mat'], 'variance_explained');
        
                        end
    
    
                        %% separability of locations
        
                        if ~isfile([outdir '/separability_' seq_name '.mat'])
    
                            load([outdir '/low_dim_space_' seq_name '.mat']);
        
                            separability = [];
                            separability.volume = cell(4,1);
                            separability.distance_by_separation = cell(4,1);
                    
                            for rank_i = 1:stim_order
                                    
                                rank = (number_locations*(rank_i-1))+(1:number_locations);
                                points = zscore(low_dim_space.z_k(rank, components));
    
                                [separability.volume{rank_i}, separability.distance_by_separation{rank_i}] = compute_separability_metrics(points);
            
                            end
            
                            save([outdir '/separability_' seq_name '.mat'], 'separability');
        
                        end
    
                        %% bootstrapp procedure to compute within-subspace metrics
    
                        % split trials within-session, compute planes for
                        % each split, and calculate angle/vaf between
                        % splits
    
                        if strcmp(performance{perf_i}, 'correct_trials')
                        
                            if ~isfile([outdir '/vaf_bootstrap_' seq_name '.mat'])
        
                                % load X_boostrap
                                filename = [indir '/X_matrix_' seq_name '_bootstrap.mat']; % X_boostrap
                                load(filename); 
                                
                                angle_bootstrap_alliter = nan(stim_order, stim_order, iterations_bootstrap);
                                vaf_bootstrap_alliter = nan(stim_order, stim_order, iterations_bootstrap);
        
                                for iter_i = 1:iterations_bootstrap
        
                                    if strcmp(events{event_i}, 'delay')
                                        
                                        [planes_split1, planes_split2, Y_r, score] = compute_subspaces_within(X_correct_trials, X_boostrap.X_split1(:,:,iter_i), X_boostrap.X_split2(:,:,iter_i), components, number_locations);
            
                                    elseif strcmp(events{event_i}, 'stim')
            
                                        [planes_split1, planes_split2, Y_r, score] = compute_subspaces_within_stim(X_correct_trials, X_boostrap.X_split1(:,:,iter_i), X_boostrap.X_split2(:,:,iter_i), components, stim_order);
            
                                    end
        
                                    % compute angle and vaf within
                                    [angle_bootstrap_alliter(:,:,iter_i), vaf_bootstrap_alliter(:,:,iter_i)] = compute_alignment_metrics_within(X_boostrap.X_split1(:,:,iter_i), planes_split1, planes_split2, Y_r, score, stim_order);
            
                                end
        
                                % save outputs
        
                                angle_bootstrap = [];
                                angle_bootstrap.median = median(angle_bootstrap_alliter,3);
                                save([outdir '/angle_bootstrap_' seq_name '.mat'], 'angle_bootstrap', '-v7.3');
        
                                vaf_bootstrap = [];
                                vaf_bootstrap.median = median(vaf_bootstrap_alliter,3);
                                save([outdir '/vaf_bootstrap_' seq_name '.mat'], 'vaf_bootstrap', '-v7.3');
                                
                            end
    
                        end

                    end
    
                end

            end

        end

    end

end

disp('Analysis done');
