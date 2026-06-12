%% GERE project

folder_inputs = 'X_matrix';
folder_outputs = 'LB25_controlled_resampling';

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

folder_X_matrix = 'X_matrix';
path_X_matrix = [path_root '/github/results/source_geometry_lm/' folder_X_matrix];

folder_original = 'LB25';
path_inputs_original = [path_root '/github/results/source_geometry_lm/' folder_original];

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
events = {'delay', 'stim'};
fun_i = @mean;
components = 1:3;
iterations_bootstrap = 1000;

%% geometric variables across iterations

disp('geometric variables across iterations')

performance = {'correct_trials'};

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

                    % subset X matrix depending on stimuli presented
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
                
                        else
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
                
                    %% compute PCA, best-fitting planes, principal angles and variance metrics
    
                    if ~isfile([outdir '/separability_perm_' seq_name '.mat'])
                            
                        angle_perm = nan(stim_order,stim_order,iterations_bootstrap);
                        vaf_perm = nan(stim_order,stim_order,iterations_bootstrap);
                        euclidean_distances_perm = zeros(number_locations * stim_order, iterations_bootstrap);
                        separability_perm = cell(iterations_bootstrap,1);
                        
                        % load X matrix
                        load([indir_correct '/X_matrix_' seq_name '.mat']); % X
                        X_before_perm = X;
                        clear X

                        if ~isempty(X_before_perm)

                            if ~isfile([outdir '/separability_perm_' seq_name '.mat'])

                                % load X matrix - controlled resampling
                                filename = [indir '/X_matrix_' seq_name '_resampling.mat'];
                                load(filename); % X_resampling

                                if ~isempty(X_resampling)
        
                                    for iter_i = 1:iterations_bootstrap
            
                                        % subset X_bootstrap iteration
                                        X = X_resampling(:,:,iter_i);
                            
                                        % pca, best-fitting planes and angles
                                        
                                        if strcmp(events{event_i}, 'delay')
                                            
                                            [z_k, z_k_exp, planes] = compute_subspaces(X_before_perm, X, components, number_locations);
                
                                        elseif strcmp(events{event_i}, 'stim')
                
                                            [z_k, z_k_exp, planes] = compute_subspaces_stim(X_before_perm, X, components, stim_order, number_locations);
                
                                        end
                                            
                                        % compute angle and vaf
                                        [angle_perm(:,:,iter_i), vaf_perm(:,:,iter_i)] = compute_alignment_metrics(planes, stim_order);
            
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
            
                                        %% separability
            
                                        separability = [];
                                        separability.volume = cell(4,1);
                                        separability.distance_by_separation = cell(4,1);
                        
                                        for rank_i = 1:stim_order
                                                
                                            rank = (number_locations*(rank_i-1))+(1:number_locations);
                                            points = zscore(low_dim_space.z_k(rank, 1:3));
                        
                                            [separability.volume{rank_i}, separability.distance_by_separation{rank_i}] = compute_separability_metrics(points);
                                            
                                        end
            
                                        separability_perm{iter_i} = separability;
            
            
                                    end
            
                                    % save
                                    save([outdir '/angle_perm_' seq_name '.mat'], 'angle_perm');
                                    save([outdir '/vaf_perm_' seq_name '.mat'], 'vaf_perm');
                                    save([outdir '/separability_perm_' seq_name '.mat'], 'separability_perm');

                                end
                                
                            end

                        end
                                                        
                    end

                end
        
            end
        
        end
    
    end

end

%% pool iterations

disp('pool iterations')

performance = {'correct_trials', 'incorrect_trials'};

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

                    outdir = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                        events{event_i} '_' num2str(seg_start) 'to' num2str(seg_end)];

                    if ~isfolder(outdir)
                        mkdir(outdir)
                    end

                    indir_original = [path_inputs_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                        events{event_i} '_' num2str(seg_start) 'to' num2str(seg_end)];

                    indir_X_matrix = [path_X_matrix '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                        events{event_i} '_' num2str(seg_start) 'to' num2str(seg_end)];
                        
                    % subset X matrix depending on stimuli presented
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
                
                        else
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


                    % load X matrix
                    load([indir_X_matrix '/X_matrix_' seq_name '.mat']); % X

                    if ~isempty(X)
                    
                        if strcmp(performance{perf_i}, 'correct_trials')
    
                            %% angle
                            load([outdir '/angle_perm_' seq_name '.mat']); % angle_perm
    
                            angle = mean(angle_perm,3);
    
                            save([outdir '/angle_' seq_name '.mat'], 'angle');
    
                            %% vaf
                            load([outdir '/vaf_perm_' seq_name '.mat']); % vaf_perm
    
                            vaf = mean(vaf_perm,3);
    
                            save([outdir '/vaf_' seq_name '.mat'], 'vaf');
        
                            %% separability
                            load([outdir '/separability_perm_' seq_name '.mat']); % separability_perm
    
                            separability = separability_perm{1};
    
                            for iter_i = 1:iterations_bootstrap
    
                                for rank_i = 1:stim_order
    
                                    if iter_i == 1
    
                                        separability.volume{rank_i} = separability_perm{iter_i}.volume{rank_i};
                                        separability.distance_by_separation{rank_i} = separability_perm{iter_i}.distance_by_separation{rank_i};
    
                                    else
    
                                        separability.volume{rank_i} = [separability.volume{rank_i} separability_perm{iter_i}.volume{rank_i}];
                                        separability.distance_by_separation{rank_i} = [separability.distance_by_separation{rank_i}; separability_perm{iter_i}.distance_by_separation{rank_i}];
    
                                    end
                            
                                end
    
                            end
    
                            for rank_i = 1:stim_order
    
                                separability.volume{rank_i} = mean(separability.volume{rank_i});
                                separability.distance_by_separation{rank_i} = mean(separability.distance_by_separation{rank_i},1);
    
                            end
    
                            save([outdir '/separability_' seq_name '.mat'], 'separability');
    
    
                            %% transfer remaining files
                            source = [indir_original '/angle_bootstrap_' seq_name '.mat'];
                            target = [outdir '/angle_bootstrap_' seq_name '.mat'];
                            system(['cp ' source ' ' target]);
    
                            source = [indir_original '/vaf_bootstrap_' seq_name '.mat'];
                            target = [outdir '/vaf_bootstrap_' seq_name '.mat'];
                            system(['cp ' source ' ' target]);
    
                        elseif strcmp(performance{perf_i}, 'incorrect_trials')
    
                            %% transfer files
                            source = [indir_original '/angle_' seq_name '.mat'];
                            target = [outdir '/angle_' seq_name '.mat'];
                            system(['cp ' source ' ' target]);
    
                            source = [indir_original '/vaf_' seq_name '.mat'];
                            target = [outdir '/vaf_' seq_name '.mat'];
                            system(['cp ' source ' ' target]);
    
                            source = [indir_original '/separability_' seq_name '.mat'];
                            target = [outdir '/separability_' seq_name '.mat'];
                            system(['cp ' source ' ' target]);
    
                        end

                    end
        
                end

            end
            
        end
    
    end

end

disp('Analysis done');
