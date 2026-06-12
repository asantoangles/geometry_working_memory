%% GERE project

folder_inputs = 'X_matrix';
folder_outputs = 'LB21_random';

% path
if isfolder('/path_to_local')
    path_root = '/path_to_local';
else
    path_root = '/path_to_hpc';
end

path_inputs = [path_root '/github/results/source_geometry_lm/X_matrix'];
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
events = {'delay', 'stim'};
fun_i = @mean;
components = 1:3;
iterations_bootstrap = 1000;

subjects = [1:15];

%% geometry of memory representations 

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
    
    disp(' '); disp(subject_ID);

    for event_i = 1:length(events)
    
        if strcmp(events{event_i}, 'delay')
    
            time_segments = delay_segments;
    
        elseif strcmp(events{event_i}, 'stim')
            
            time_segments = stim_segments;
    
        end
                
        % loop over time windows
        for seg_i = 1:size(time_segments, 1)

            %% geometric analysis
    
            % loop over performance
            for perf_i = 1:length(performance)

                % loop over sequence length used
                for sequence_length = [3 4]

                    if sequence_length == 3
                        seq_name = 'length3';
                        ranks = 3;
                    elseif sequence_length == 4
                        seq_name = 'length4';
                        ranks = 4;
                    end

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
                    
                    indir = [path_inputs '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                        events{event_i} '_' num2str(time_segments(seg_i,1)) 'to' num2str(time_segments(seg_i,2))];
                    
                    indir_correct = [path_inputs '/' subject_ID '/' func2str(fun_i) '/correct_trials/delay_1to4000'];
                    
                    outdir = [path_results '/rand' num2str(random_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} '/' ...
                            events{event_i} '_' num2str(time_segments(seg_i,1)) 'to' num2str(time_segments(seg_i,2))];

                    if ~isfolder(outdir)
                        mkdir(outdir);
                    end

                    % load X matrix
                    load([indir '/X_matrix_' seq_name '.mat']); % X
                    
                    % load X matrix for correct trials where to compute eigenvectors
                    X_correct_trials = load([indir_correct '/X_matrix_' seq_name '.mat']);
                    X_correct_trials = X_correct_trials.X;

                    %% randomization

                    if ~isempty(X)

                        rng('shuffle');
    
                        X = X(randperm(size(X, 1)), randperm(size(X, 2)));
                        X_correct_trials = X_correct_trials(randperm(size(X_correct_trials, 1)), randperm(size(X_correct_trials, 2)));
    
                        %% compute PCA, best-fitting planes, principal angles and variance metrics
    
                        if ~isfile([outdir '/vaf_' seq_name '.mat'])
                        
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

                    end
            
                end

            end

        end

    end

end

disp('Analysis done');
