%% GERE project

folder_outputs = 'NB23_PC1andPC2_random';

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
events = {'delay', 'stim'};
fun_i = @mean;
components = 1:3;
iterations_bootstrap = 1000;




%%%%%%%%%%%%%%%%%%
%% separability %%
%%%%%%%%%%%%%%%%%%

disp('separability');

%%% delay

time_segments = delay_segments;

% loop over time windows
for time_i = 1:size(time_segments, 1)

    % loop over sequence length used
    for sequence_length = [3 4]

        if sequence_length == 3
            seq_name = 'length3';
            ranks = 3;
        elseif sequence_length == 4
            seq_name = 'length4';
            ranks = 4;
        end

        for perf_i = 1:length(performance)

            outdir_volume = [path_results '/group_results/volume/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'];
            outdir_distance = [path_results '/group_results/distance/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'];

            if ~isfolder(outdir_volume)
                mkdir(outdir_volume);
            end

            if ~isfolder(outdir_distance)
                mkdir(outdir_distance);
            end

            if ~isfile([outdir_volume '/delay_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt'])

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

                    % pool iterations_bootstrap
                    for rand_i = 1:iterations_bootstrap
                        
                        % delay
                        load([path_results '/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/delay_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2))...
                            '/separability_' seq_name '.mat']); % separability

                        if ranks == 3

                            volume_delay = [separability.volume{1} separability.volume{2} separability.volume{3}];
                            distance_delay = [separability.distance_by_separation{1}; separability.distance_by_separation{2}; separability.distance_by_separation{3}];

                        else

                            volume_delay = [separability.volume{1} separability.volume{2} separability.volume{3}, separability.volume{4}];
                            distance_delay = [separability.distance_by_separation{1}; separability.distance_by_separation{2}; separability.distance_by_separation{3}; separability.distance_by_separation{4}];

                        end

                        if rand_i == 1

                            distance_delay_all = distance_delay(:);
                            volume_delay_all = volume_delay(:);

                        else

                            distance_delay_all = [distance_delay_all distance_delay(:)];

                            if size(volume_delay, 2) == size(volume_delay_all, 1)
                                volume_delay_all = [volume_delay_all volume_delay(:)];
                            end

                        end

                    end

                    distance_delay = mean(distance_delay_all,2);
                    volume_delay = mean(volume_delay_all, 2);

                    % volume
                    for row_i = 1:length(volume_delay)
                        matrix_design_volume(size(matrix_design_volume,1)+1,2) = volume_delay(row_i);
                        matrix_design_volume(size(matrix_design_volume,1),1) = 1;
                        matrix_design_volume(size(matrix_design_volume,1),3) = 1;
                        matrix_design_volume(size(matrix_design_volume,1),4) = sub_i;
                    end

                    % dist
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
                writetable(data, [outdir_volume '/delay_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt']);

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
                writetable(data, [outdir_distance '/delay_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt']);
                    
            end
            
        end

    end

end

%%% stim

time_segments = stim_segments;

% loop over time windows
for time_i = 1:size(time_segments, 1)

    % loop over sequence length used
    for sequence_length = [3 4]

        if sequence_length == 3
            seq_name = 'length3';
            ranks = 3;
        elseif sequence_length == 4
            seq_name = 'length4';
            ranks = 4;
        end

        % subset X matrix depending on stimuli presented
        if sequence_length == 3
            if time_segments(time_i,1) == 1 && time_segments(time_i,2) == 300
                stim_order = 1;
            elseif time_segments(time_i,1) == 400 && time_segments(time_i,2) == 700
                stim_order = 2;
            elseif time_segments(time_i,1) == 800 && time_segments(time_i,2) == 1100
                stim_order = 3;
            elseif time_segments(time_i,1) == 1200 && time_segments(time_i,2) == 1500
                stim_order = 3;
            end

        else
            if time_segments(time_i,1) == 1 && time_segments(time_i,2) == 300
                stim_order = 1;
            elseif time_segments(time_i,1) == 400 && time_segments(time_i,2) == 700
                stim_order = 2;
            elseif time_segments(time_i,1) == 800 && time_segments(time_i,2) == 1100
                stim_order = 3;
            elseif time_segments(time_i,1) == 1200 && time_segments(time_i,2) == 1500
                stim_order = 4;
            end
        end
        
        for perf_i = 1:length(performance)

            outdir_volume = [path_results '/group_results/volume/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'];
            outdir_distance = [path_results '/group_results/distance/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'];

            if ~isfolder(outdir_volume)
                mkdir(outdir_volume);
            end

            if ~isfolder(outdir_distance)
                mkdir(outdir_distance);
            end
            
            if ~isfile([outdir_volume '/stim_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt'])

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

                    distance_stim_all = [];
                    volume_stim_all = [];
                    
                    % pool iterations_bootstrap
                    for rand_i = 1:iterations_bootstrap

                        if isfile([path_results '/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/stim_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2))...
                            '/separability_' seq_name '.mat'])
                        
                            % delay
                            load([path_results '/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                                '/stim_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2))...
                                '/separability_' seq_name '.mat']); % separability
    
                            if stim_order == 1
                                volume_stim = [separability.volume{1}];
                                distance_stim = [separability.distance_by_separation{1}];
                            end
        
                            if stim_order == 2
                                volume_stim = [separability.volume{1} separability.volume{2}];
                                distance_stim = [separability.distance_by_separation{1}; separability.distance_by_separation{2}];
                            end
        
                            if stim_order == 3
                                volume_stim = [separability.volume{1} separability.volume{2} separability.volume{3}];
                                distance_stim = [separability.distance_by_separation{1}; separability.distance_by_separation{2}; separability.distance_by_separation{3}];
                            end
        
                            if stim_order == 4
                                volume_stim = [separability.volume{1} separability.volume{2} separability.volume{3}, separability.volume{4}];
                                distance_stim = [separability.distance_by_separation{1}; separability.distance_by_separation{2}; separability.distance_by_separation{3}; separability.distance_by_separation{4}];
                            end
    
                            if rand_i == 1
    
                                distance_stim_all = distance_stim(:);
                                volume_stim_all = volume_stim(:);
    
                            else
    
                                distance_stim_all = [distance_stim_all distance_stim(:)];
                                if size(volume_stim, 2) == size(volume_stim_all, 1)
                                    volume_stim_all = [volume_stim_all volume_stim(:)];
                                end
    
                            end

                        end

                    end

                    distance_stim = mean(distance_stim_all,2);
                    volume_stim = mean(volume_stim_all, 2);

                    % volume
                    for row_i = 1:length(volume_stim)
                        matrix_design_volume(size(matrix_design_volume,1)+1,2) = volume_stim(row_i);
                        matrix_design_volume(size(matrix_design_volume,1),1) = 1;
                        matrix_design_volume(size(matrix_design_volume,1),3) = 1;
                        matrix_design_volume(size(matrix_design_volume,1),4) = sub_i;
                    end

                    % dist
                    for row_i = 1:length(distance_stim)
                        matrix_design_dist(size(matrix_design_dist,1)+1,2) = distance_stim(row_i);
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
                writetable(data, [outdir_volume '/stim_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt']);

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
                writetable(data, [outdir_distance '/stim_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt']);
                    
            end
            
        end

    end

end




%%%%%%%%%%%%%%%%%%%%%
%% principal angle %%
%%%%%%%%%%%%%%%%%%%%%

disp('principal angle');

%%% delay 

time_segments = delay_segments;

% loop over time windows
for time_i = 1:size(time_segments, 1)

    % loop over sequence length used
    for sequence_length = [3 4]

        if sequence_length == 3
            seq_name = 'length3';
            ranks = 3;
        elseif sequence_length == 4
            seq_name = 'length4';
            ranks = 4;
        end

        for perf_i = 1:length(performance)

            outdir = [path_results '/group_results/principal_angle/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'];

            if ~isfolder(outdir)
                mkdir(outdir);
            end

            if ~isfile([outdir '/delay_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt'])

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

                    % pool iterations_bootstrap
                    for rand_i = 1:iterations_bootstrap
                        
                        % load data
                        load([path_results '/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/delay_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2))...
                            '/angle_' seq_name '.mat']); % angle
                                
                        if rand_i == 1

                            angle_all = angle;

                        else

                            angle_all = angle_all + angle;

                        end

                    end

                    angle = angle_all / iterations_bootstrap;
                        
                    % store between angles
                    between = angle;
                    between(logical(eye(size(between)))) = 0;
                    between = triu(between); between = between(:);
                    between(between == 0) = [];

                    first_row = size(between_within,1)+1;
                    tmp = [between];
                    last_row = first_row + length(tmp) - 1;

                    between_within(first_row:last_row,1) = tmp;
                    between_within(first_row:last_row,2) = [ones(length(between),1)];
                    between_within(first_row:last_row,3) = sub_i;

                end
    
                % save between_within
                writetable(table(between_within), [outdir '/delay_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt']);
                
            end

        end

    end

end



%%% stim

time_segments = stim_segments;

% loop over time windows
for time_i = 1:size(time_segments, 1)

    % loop over sequence length used
    for sequence_length = [3 4]

        if sequence_length == 3
            seq_name = 'length3';
            ranks = 3;
        elseif sequence_length == 4
            seq_name = 'length4';
            ranks = 4;
        end
        
        for perf_i = 1:length(performance)

            outdir = [path_results '/group_results/principal_angle/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'];

            if ~isfolder(outdir)
                mkdir(outdir);
            end

            if ~isfile([outdir '/stim_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt'])

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

                    angle_all = [];

                    % pool iterations_bootstrap
                    for rand_i = 1:iterations_bootstrap

                        if isfile([path_results '/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/stim_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2))...
                            '/angle_' seq_name '.mat'])
                        
                            % load data
                            load([path_results '/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                                '/stim_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2))...
                                '/angle_' seq_name '.mat']); % angle
                                    
                            if rand_i == 1
    
                                angle_all = angle;
    
                            else
    
                                angle_all = angle_all + angle;
    
                            end

                        end

                    end

                    angle = angle_all / iterations_bootstrap;
                        

                    % store between angles
                    between = angle;
                    between(logical(eye(size(between)))) = 0;
                    between = triu(between); between = between(:);
                    between(between == 0) = [];

                    first_row = size(between_within,1)+1;
                    tmp = [between];
                    last_row = first_row + length(tmp) - 1;

                    between_within(first_row:last_row,1) = tmp;
                    between_within(first_row:last_row,2) = [ones(length(between),1)];
                    between_within(first_row:last_row,3) = sub_i;

                end
    
                % save between_within
                writetable(table(between_within), [outdir '/stim_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt']);
                
            end

        end

    end

end




%%%%%%%%%
%% vaf %%
%%%%%%%%%

disp('vaf');

% It equals 0 if the two subspaces are orthogonal and 
% equals 1 if they completely overlap with each other.

% lme intercept tests whether vaf(between - within) is
% different from zero

% if t-stat < 0, vaf_between < vaf_within, so 
% between-subsapces is more orthogonal than within-subspaces


%%% delay 

time_segments = delay_segments;

% loop over time windows
for time_i = 1:size(time_segments, 1)

    % loop over sequence length used
    for sequence_length = [3 4]

        if sequence_length == 3
            seq_name = 'length3';
            ranks = 3;
        elseif sequence_length == 4
            seq_name = 'length4';
            ranks = 4;
        end

        for perf_i = 1:length(performance)

            outdir = [path_results '/group_results/vaf/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'];

            if ~isfolder(outdir)
                mkdir(outdir);
            end

            if ~isfile([outdir '/delay_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt'])

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

                    % pool iterations_bootstrap
                    for rand_i = 1:iterations_bootstrap
                        
                        % load data
                        load([path_results '/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/delay_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2))...
                            '/vaf_' seq_name '.mat']); % vaf
                                
                        if rand_i == 1

                            vaf_all = vaf;

                        else

                            vaf_all = vaf_all + vaf;

                        end

                    end

                    vaf = vaf_all / iterations_bootstrap;
                        
                    % store between angles
                    between = vaf;
                    between(logical(eye(size(between)))) = 0;
                    between = triu(between); between = between(:);
                    between(between == 0) = [];

                    first_row = size(between_within,1)+1;
                    tmp = [between];
                    last_row = first_row + length(tmp) - 1;

                    between_within(first_row:last_row,1) = tmp;
                    between_within(first_row:last_row,2) = [ones(length(between),1)];
                    between_within(first_row:last_row,3) = sub_i;

                end
    
                % save between_within
                writetable(table(between_within), [outdir '/delay_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt']);
                
            end

        end

    end

end


%%% stim 

time_segments = stim_segments;

% loop over time windows
for time_i = 1:size(time_segments, 1)

    % loop over sequence length used
    for sequence_length = [3 4]

        if sequence_length == 3
            seq_name = 'length3';
            ranks = 3;
        elseif sequence_length == 4
            seq_name = 'length4';
            ranks = 4;
        end

        for perf_i = 1:length(performance)

            outdir = [path_results '/group_results/vaf/' func2str(fun_i) '/' performance{perf_i} '/design_matrix'];

            if ~isfolder(outdir)
                mkdir(outdir);
            end

            if ~isfile([outdir '/stim_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt'])

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

                    vaf_all = [];

                    % pool iterations_bootstrap
                    for rand_i = 1:iterations_bootstrap

                        if isfile([path_results '/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                            '/stim_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2))...
                            '/vaf_' seq_name '.mat'])
                        
                            % load data
                            load([path_results '/rand' num2str(rand_i) '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i} ...
                                '/stim_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2))...
                                '/vaf_' seq_name '.mat']); % vaf
                                    
                            if rand_i == 1
    
                                vaf_all = vaf;
    
                            else
    
                                vaf_all = vaf_all + vaf;
    
                            end

                        end

                    end

                    vaf = vaf_all / iterations_bootstrap;
                        

                    % store between angles
                    between = vaf;
                    between(logical(eye(size(between)))) = 0;
                    between = triu(between); between = between(:);
                    between(between == 0) = [];

                    first_row = size(between_within,1)+1;
                    tmp = [between];
                    last_row = first_row + length(tmp) - 1;

                    between_within(first_row:last_row,1) = tmp;
                    between_within(first_row:last_row,2) = [ones(length(between),1)];
                    between_within(first_row:last_row,3) = sub_i;

                end
    
                % save between_within
                writetable(table(between_within), [outdir '/stim_' seq_name '_' num2str(time_segments(time_i,1)) 'to' num2str(time_segments(time_i,2)) '.txt']);
                
            end

        end

    end

end



disp('Analysis done');
