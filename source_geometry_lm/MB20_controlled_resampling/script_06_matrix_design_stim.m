%% GERE project
% geometry projection

% path to data
if isfolder('/path_to_local')
    path_inputs = '/path_to_local/results/source_reconstruction';
    path_results = ['/path_to_local/results/source_geometry_lm/' folder_geometry];
    addpath('/path_to_local/scripts/source_geometry_lm/utilities')
    folder_original = 'MB20';
    path_results_original = ['/path_to_local/results/source_geometry_lm/' folder_original];
    subjects = [5 7 18 23 25 31 34 37 40 45 47 53 61 201 202];
    sessions = 1:2;
end

% settings loops
functions_data ={@mean};
performance = {'correct_trials' 'incorrect_trials'};

% time windows
stim_resolved_segments = [];
stim_resolved_segments(end+1,:) = [1 300];
stim_resolved_segments(end+1,:) = [400 700];
stim_resolved_segments(end+1,:) = [800 1100];
stim_resolved_segments(end+1,:) = [1200 1500];

%%%%%%%%%%%%%%%%%%
%% separability %%
%%%%%%%%%%%%%%%%%%

%%% stim

time_segments = stim_resolved_segments;

% loop over functions
for fun_i = functions_data

    fun_i = fun_i{1};

    % subset by time windows within delay period
    for delay_i = 1:size(time_segments, 1)

        disp(['stim_resolved_refined_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);

        % loop over sequence length used
        for sequence_length = [3 4 34]

            if sequence_length == 3
                sequence_length_filename = 'length3';
                ranks = 3;
            elseif sequence_length == 4
                sequence_length_filename = 'length4';
                ranks = 4;
            elseif sequence_length == 34
                sequence_length_filename = 'lengthall';
                ranks = 4;
            end

            % NOTE change here relative to main script
            % subset X matrix depending on stimuli presented
            if sequence_length == 3
                if time_segments(delay_i,1) == 1 && time_segments(delay_i,2) == 300
                    stim_order = 1;
                elseif time_segments(delay_i,1) == 400 && time_segments(delay_i,2) == 700
                    stim_order = 1;
                elseif time_segments(delay_i,1) == 800 && time_segments(delay_i,2) == 1100
                    stim_order = 2;
                elseif time_segments(delay_i,1) == 1200 && time_segments(delay_i,2) == 1500
                    stim_order = 3;
                end

            else
                if time_segments(delay_i,1) == 1 && time_segments(delay_i,2) == 300
                    stim_order = 1;
                elseif time_segments(delay_i,1) == 400 && time_segments(delay_i,2) == 700
                    stim_order = 2;
                elseif time_segments(delay_i,1) == 800 && time_segments(delay_i,2) == 1100
                    stim_order = 3;
                elseif time_segments(delay_i,1) == 1200 && time_segments(delay_i,2) == 1500
                    stim_order = 4;
                end
            end

            for perf_i = 1:length(performance)

                if ~isfile([path_results '/group_results/volume/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'...
                    '/design_matrix_stim_resolved_refined_' sequence_length_filename '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) '.txt'])

                    if ~isfolder([path_results '/group_results/volume/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'])
                        mkdir([path_results '/group_results/volume/' func2str(fun_i) '/' performance{perf_i} '/design_matrix']);
                    end
                    if ~isfolder([path_results '/group_results/distance/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'])
                        mkdir([path_results '/group_results/distance/' func2str(fun_i) '/' performance{perf_i} '/design_matrix']);
                    end

                    matrix_design_volume = nan(0,5); % columns: (1) intercept, (2) variable of interest (e.g. PA), (3) condition, (4) subject, (5) session
                    matrix_design_dist = nan(0,5); % columns: (1) intercept, (2) variable of interest (e.g. PA), (3) condition, (4) subject, (5) session
    
                    for sub_i = 1:length(subjects)
                    
                        subject = subjects(sub_i);
                                    
                        %% load data
                
                        % set paths
                        if subject < 10
                            subject_ID = ['sub_0' num2str(subject)];
                            subjectID = ['sub0' num2str(subject)];
                        else
                            subject_ID = ['sub_' num2str(subject)];
                            subjectID = ['sub' num2str(subject)];
                        end
                            
                        % delay
                        load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/stim_resolved_refined_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                            '/separability_' sequence_length_filename '.mat']); % separability

                        if stim_order == 1

                            volume_delay = [separability.volume{1}];
                            distance_delay = [separability.distance_by_separation{1}];

                        end

                        if stim_order == 2
                            
                            volume_delay = [separability.volume{1} separability.volume{2}];
                            distance_delay = [separability.distance_by_separation{1}; separability.distance_by_separation{2}];

                        end


                        if stim_order == 3

                            volume_delay = [separability.volume{1} separability.volume{2} separability.volume{3}];
                            distance_delay = [separability.distance_by_separation{1}; separability.distance_by_separation{2}; separability.distance_by_separation{3}];

                        end

                        if stim_order == 4

                            volume_delay = [separability.volume{1} separability.volume{2} separability.volume{3}, separability.volume{4}];
                            distance_delay = [separability.distance_by_separation{1}; separability.distance_by_separation{2}; separability.distance_by_separation{3}; separability.distance_by_separation{4}];

                        end

                        % volume
                        for row_i = 1:length(volume_delay)
                            matrix_design_volume(size(matrix_design_volume,1)+1,2) = volume_delay(row_i);
                            matrix_design_volume(size(matrix_design_volume,1),1) = 1;
                            matrix_design_volume(size(matrix_design_volume,1),3) = 1;
                            matrix_design_volume(size(matrix_design_volume,1),4) = sub_i;
                        end

                        % dist
                        distance_delay = distance_delay(:);

                        for row_i = 1:length(distance_delay)
                            matrix_design_dist(size(matrix_design_dist,1)+1,2) = distance_delay(row_i);
                            matrix_design_dist(size(matrix_design_dist,1),1) = 1;
                            matrix_design_dist(size(matrix_design_dist,1),3) = 1;
                            matrix_design_dist(size(matrix_design_dist,1),4) = sub_i;
                        end                        
                        
                    end

                    %% volume

                    matrix_design = matrix_design_volume;
    
                    % response variable
                    response = matrix_design(:,2);
    
                    % random intercepts
                    condition = categorical(matrix_design(:,3));
                    subject = categorical(matrix_design(:,4));
                    session = categorical(matrix_design(:,5));
    
                    % Create a table from the data
                    data = table(response, condition, subject, session);
    
                    % save data matrix
                    writetable(data, [path_results '/group_results/volume/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'...
                        '/design_matrix_stim_resolved_refined_' sequence_length_filename '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) '.txt']);
    
                    %% distance

                    matrix_design = matrix_design_dist;
    
                    % response variable
                    response = matrix_design(:,2);
    
                    % random intercepts
                    condition = categorical(matrix_design(:,3));
                    subject = categorical(matrix_design(:,4));
                    session = categorical(matrix_design(:,5));
    
                    % Create a table from the data
                    data = table(response, condition, subject, session);
    
                    % save data matrix
                    writetable(data, [path_results '/group_results/distance/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'...
                        '/design_matrix_stim_resolved_refined_' sequence_length_filename '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) '.txt']);
    
                    
                end
                
            end

        end

    end

end




%%%%%%%%%%%%%%%%%%%%%
%% principal angle %%
%%%%%%%%%%%%%%%%%%%%%


%%% stim 

time_segments = stim_resolved_segments;

% loop over functions
for fun_i = functions_data

    fun_i = fun_i{1};

    % subset by time windows within delay period
    for delay_i = 1:size(time_segments, 1)

        disp(['stim_resolved_refined_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);

        % loop over sequence length used
        for sequence_length = [3 4]

            if sequence_length == 3
                sequence_length_filename = 'length3';
                ranks = 3;
            elseif sequence_length == 4
                sequence_length_filename = 'length4';
                ranks = 4;
            elseif sequence_length == 34
                sequence_length_filename = 'lengthall';
                ranks = 4;
            end

            % NOTE change here relative to main script
            % subset X matrix depending on stimuli presented
            if sequence_length == 3
                if time_segments(delay_i,1) == 1 && time_segments(delay_i,2) == 300
                    stim_order = 1;
                elseif time_segments(delay_i,1) == 400 && time_segments(delay_i,2) == 700
                    stim_order = 1;
                elseif time_segments(delay_i,1) == 800 && time_segments(delay_i,2) == 1100
                    stim_order = 2;
                elseif time_segments(delay_i,1) == 1200 && time_segments(delay_i,2) == 1500
                    stim_order = 3;
                end

            else
                if time_segments(delay_i,1) == 1 && time_segments(delay_i,2) == 300
                    stim_order = 1;
                elseif time_segments(delay_i,1) == 400 && time_segments(delay_i,2) == 700
                    stim_order = 2;
                elseif time_segments(delay_i,1) == 800 && time_segments(delay_i,2) == 1100
                    stim_order = 3;
                elseif time_segments(delay_i,1) == 1200 && time_segments(delay_i,2) == 1500
                    stim_order = 4;
                end
            end

            for perf_i = 1:length(performance)

                if ~isfile([path_results '/group_results/principal_angle/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'...
                    '/between_within_stim_resolved_refined_' sequence_length_filename '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) '.txt'])

                    matrix_design = nan(0,4); % columns: (1) intercept, (2) variable of interest (e.g. PA), (3) subject, (4) session
                    between_within = nan(0,4); % columns: (1) variable of interest (PAs), (2) between [1] or within [2], (3) subject, (4) session
    
                    % pool data into the matrix design
                    for sub_i = 1:length(subjects)
                    
                        subject = subjects(sub_i);
                                    
                        %% load data
                
                        % set paths
                        if subject < 10
                            subject_ID = ['sub_0' num2str(subject)];
                            subjectID = ['sub0' num2str(subject)];
                        else
                            subject_ID = ['sub_' num2str(subject)];
                            subjectID = ['sub' num2str(subject)];
                        end
                            
                        % load data
                        load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/stim_resolved_refined_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                            '/angle_' sequence_length_filename '.mat']); % angle

                        load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/stim_resolved_refined_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                            '/angle_bootstrap_' sequence_length_filename '.mat']); % angle_bootstrap

                        % between minus within PAs
                        data = angle - median([diag(angle_bootstrap.median)]);

                        % DISCARDED CODE
                        % data = zeros(size(angle));
                        % for row_i = 1:size(angle,1)
                        %     for col_i = 1:size(angle,1)
                        %         if row_i ~= col_i
                        %             data(row_i, col_i) = angle(row_i, col_i) - median([angle_bootstrap.median(row_i, row_i) angle_bootstrap.median(col_i, col_i)]);
                        %         end
                        %     end
                        % end

                        % set design matrix
                        if size(data,1) > 1
                            data(logical(eye(size(data)))) = 0;
                            data = triu(data); data = data(:);
                            data(data == 0) = [];
                        end

                        first_row = size(matrix_design,1)+1;

                        for item_i = 1:length(data)
                            matrix_design(end+1,2) = data(item_i);
                        end

                        last_row = size(matrix_design,1);

                        matrix_design(first_row:last_row,1) = 1;
                        matrix_design(first_row:last_row,3) = sub_i;
                        

                        % store between and within angles
                        between = angle;
                        if size(between,1) > 1
                            between(logical(eye(size(between)))) = 0;
                            between = triu(between); between = between(:);
                            between(between == 0) = [];
                        end

                        within = diag(angle_bootstrap.median);

                        first_row = size(between_within,1)+1;
                        tmp = [between; within];
                        last_row = first_row + length(tmp) - 1;

                        between_within(first_row:last_row,1) = tmp;
                        between_within(first_row:last_row,2) = [ones(length(between),1); 2*ones(length(within),1)];
                        between_within(first_row:last_row,3) = sub_i;
                        
    
                    end
    
                    % response variable
                    response = matrix_design(:,2);
    
                    % random intercepts
                    subject = categorical(matrix_design(:,3));
                    session = categorical(matrix_design(:,4));
    
                    % Create a table from the data
                    data = table(response, subject, session);
    
                    % save data matrix
                    writetable(data, [path_results '/group_results/principal_angle/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'...
                        '/design_matrix_stim_resolved_refined_' sequence_length_filename '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) '.txt']);
    
                    % save between_within
                    writetable(table(between_within), [path_results '/group_results/principal_angle/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'...
                        '/between_within_stim_resolved_refined_' sequence_length_filename '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) '.txt']);
                    

                end

            end

        end

    end

end



%%%%%%%%%
%% vaf %%
%%%%%%%%%

% It equals 0 if the two subspaces are orthogonal and 
% equals 1 if they completely overlap with each other.

% lme intercept tests whether vaf(between - within) is
% different from zero

% if t-stat < 0, vaf_between < vaf_within, so 
% between-subsapces is more orthogonal than within-subspaces


%%% stim 

time_segments = stim_resolved_segments;

% loop over functions
for fun_i = functions_data

    fun_i = fun_i{1};

    % subset by time windows within delay period
    for delay_i = 1:size(time_segments, 1)

        disp(['stim_resolved_refined_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);

        % loop over sequence length used
        for sequence_length = [3 4]

            if sequence_length == 3
                sequence_length_filename = 'length3';
                ranks = 3;
            elseif sequence_length == 4
                sequence_length_filename = 'length4';
                ranks = 4;
            elseif sequence_length == 34
                sequence_length_filename = 'lengthall';
                ranks = 4;
            end

            % NOTE change here relative to main script
            % subset X matrix depending on stimuli presented
            if sequence_length == 3
                if time_segments(delay_i,1) == 1 && time_segments(delay_i,2) == 300
                    stim_order = 1;
                elseif time_segments(delay_i,1) == 400 && time_segments(delay_i,2) == 700
                    stim_order = 1;
                elseif time_segments(delay_i,1) == 800 && time_segments(delay_i,2) == 1100
                    stim_order = 2;
                elseif time_segments(delay_i,1) == 1200 && time_segments(delay_i,2) == 1500
                    stim_order = 3;
                end

            else
                if time_segments(delay_i,1) == 1 && time_segments(delay_i,2) == 300
                    stim_order = 1;
                elseif time_segments(delay_i,1) == 400 && time_segments(delay_i,2) == 700
                    stim_order = 2;
                elseif time_segments(delay_i,1) == 800 && time_segments(delay_i,2) == 1100
                    stim_order = 3;
                elseif time_segments(delay_i,1) == 1200 && time_segments(delay_i,2) == 1500
                    stim_order = 4;
                end
            end

            for perf_i = 1:length(performance)
                
                if ~isfile([path_results '/group_results/vaf/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'...
                    '/between_within_stim_resolved_refined_' sequence_length_filename '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) '.txt'])

                    matrix_design = nan(0,4); % columns: (1) intercept, (2) variable of interest (e.g. PA), (3) subject, (4) session
                    between_within = nan(0,4); % columns: (1) variable of interest (PAs), (2) between [1] or within [2], (3) subject, (4) session
    
                    % pool data into the matrix design
                    for sub_i = 1:length(subjects)
                    
                        subject = subjects(sub_i);
                    
                        %% load data
                
                        % set paths
                        if subject < 10
                            subject_ID = ['sub_0' num2str(subject)];
                            subjectID = ['sub0' num2str(subject)];
                        else
                            subject_ID = ['sub_' num2str(subject)];
                            subjectID = ['sub' num2str(subject)];
                        end
                            
                        % load data
                        load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/stim_resolved_refined_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                            '/vaf_' sequence_length_filename '.mat']); % vaf

                        load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/stim_resolved_refined_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                            '/vaf_bootstrap_' sequence_length_filename '.mat']); % vaf_bootstrap

                        % between minus within PAs
                        data = vaf - median([diag(vaf_bootstrap.median)]);

                        % DISCARDED CODE
                        % data = zeros(size(vaf));
                        % for row_i = 1:size(vaf,1)
                        %     for col_i = 1:size(vaf,1)
                        %         if row_i ~= col_i
                        %             data(row_i, col_i) = vaf(row_i, col_i) - median([vaf_bootstrap.median(row_i, row_i) vaf_bootstrap.median(col_i, col_i)]);
                        %         end
                        %     end
                        % end

                        % set design matrix
                        if size(data,1) > 1
                            data(logical(eye(size(data)))) = 0;
                            data = triu(data); data = data(:);
                            data(data == 0) = [];
                        end
                        

                        first_row = size(matrix_design,1)+1;

                        for item_i = 1:length(data)
                            matrix_design(end+1,2) = data(item_i);
                        end

                        last_row = size(matrix_design,1);

                        matrix_design(first_row:last_row,1) = 1;
                        matrix_design(first_row:last_row,3) = sub_i;
                        

                        % store between and within angles
                        between = vaf;
                        if size(between,1) > 1
                            between(logical(eye(size(between)))) = 0;
                            between = triu(between); between = between(:);
                            between(between == 0) = [];
                        end

                        within = diag(vaf_bootstrap.median);

                        first_row = size(between_within,1)+1;
                        tmp = [between; within];
                        last_row = first_row + length(tmp) - 1;

                        between_within(first_row:last_row,1) = tmp;
                        between_within(first_row:last_row,2) = [ones(length(between),1); 2*ones(length(within),1)];
                        between_within(first_row:last_row,3) = sub_i;
                        
    
                    end
    
                    % response variable
                    response = matrix_design(:,2);
    
                    % random intercepts
                    subject = categorical(matrix_design(:,3));
                    session = categorical(matrix_design(:,4));
    
                    % Create a table from the data
                    data = table(response, subject, session);
    
                    % save data matrix
                    writetable(data, [path_results '/group_results/vaf/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'...
                        '/design_matrix_stim_resolved_refined_' sequence_length_filename '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) '.txt']);
    
                    % save between_within
                    writetable(table(between_within), [path_results '/group_results/vaf/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'...
                        '/between_within_stim_resolved_refined_' sequence_length_filename '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2)) '.txt']);
                    
                end

            end

        end

    end

end


disp('Analysis done');
