%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = '/path_to_local/results/source_reconstruction';
    path_results = ['/path_to_local/results/source_geometry_lm/' folder_geometry '/comp1to' num2str(last_comp_i)];
    addpath('/path_to_local/scripts/source_geometry_lm/utilities')
    functions_data ={@mean};
end

% settings loops
stim_resolved_segments = [];
stim_resolved_segments(end+1,:) = [1 300];
stim_resolved_segments(end+1,:) = [400 700];
stim_resolved_segments(end+1,:) = [800 1100];
stim_resolved_segments(end+1,:) = [1200 1500];

events = {'stim_resolved_refined'};
performance = {'correct_trials', 'incorrect_trials'};
iterations_bootstrap = 1000;
permutations = 1000;
components = 1:3;
sessions = 1:2;

number_locations = 8;

%% pool permutations stim 

for last_comp_i = [6 8]

    path_results = ['/path_to_local/results/source_geometry_lm/' folder_geometry '/comp1to' num2str(last_comp_i)];

    folder_original = 'LB20_allPCs_accumulated';
    path_results_original = ['/path_to_local/results/source_geometry_lm/' folder_original '/comp1to' num2str(last_comp_i)];

    for event_i = 1:length(events)
    
        if strcmp(events{event_i}, 'stim_resolved_refined')
            
            time_segments = stim_resolved_segments;
    
        end
        
        for perf_i = 1:length(performance)
    
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
                    
                % loop over functions
                for fun_i = functions_data
        
                    fun_i = fun_i{1};
                    
                    % subset by time windows within delay period
                    for delay_i = 1:size(time_segments, 1)
    
                        disp([events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);
                
                        if ~isfolder([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))])
                            mkdir([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);
                        end
    
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
                            
                            %% pool permutations
    
                            if strcmp(performance{perf_i}, 'correct_trials')
    
                                %%% angle
                                load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_perm_' sequence_length_filename '.mat']); % angle_perm
        
                                angle = mean(angle_perm,3);
        
                                save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_' sequence_length_filename '.mat'], 'angle');

                                %%% angle min
                                load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_min_perm_' sequence_length_filename '.mat']); % angle_min_perm
        
                                angle_min = mean(angle_min_perm,3);
        
                                save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_min_' sequence_length_filename '.mat'], 'angle_min');
        
        
                                %%% vaf
                                load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/vaf_perm_' sequence_length_filename '.mat']); % vaf_perm
        
                                vaf = mean(vaf_perm,3);
        
                                save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/vaf_' sequence_length_filename '.mat'], 'vaf');
            
                                %%% separability
                                load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_perm_' sequence_length_filename '.mat']); % separability_perm
        
                                separability = separability_perm{1};
        
                                for perm_i = 1:permutations
        
                                    for rank_i = 1:stim_order
        
                                        if perm_i == 1
        
                                            separability.volume{rank_i} = separability_perm{perm_i}.volume{rank_i};
                                            separability.distance_by_separation{rank_i} = separability_perm{perm_i}.distance_by_separation{rank_i};
        
                                        else

                                            try
        
                                                separability.volume{rank_i} = [separability.volume{rank_i} separability_perm{perm_i}.volume{rank_i}];
                                                separability.distance_by_separation{rank_i} = [separability.distance_by_separation{rank_i}; separability_perm{perm_i}.distance_by_separation{rank_i}];
        
                                            catch
                                            end

                                        end
                                
                                    end
        
                                end
        
                                for rank_i = 1:stim_order
        
                                    separability.volume{rank_i} = mean(separability.volume{rank_i});
                                    separability.distance_by_separation{rank_i} = mean(separability.distance_by_separation{rank_i},1);
        
                                end
        
                                save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_' sequence_length_filename '.mat'], 'separability');
        
        
                                % transfer files
                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_bootstrap_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_bootstrap_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);
        
                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_min_bootstrap_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_min_bootstrap_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);
                                
                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/vaf_bootstrap_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/vaf_bootstrap_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);
    
                            elseif strcmp(performance{perf_i}, 'incorrect_trials')
    
                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);

                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_min_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_min_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);
    
                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/vaf_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/vaf_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);
    
                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);
    
                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_allpairs_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_allpairs_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);
    
                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_bootstrap_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_bootstrap_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);

                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_min_bootstrap_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/angle_min_bootstrap_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);
        
                                source = [path_results_original '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/vaf_bootstrap_' sequence_length_filename '.mat'];
                                target = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/vaf_bootstrap_' sequence_length_filename '.mat'];
                                system(['cp ' source ' ' target]);
    
                            end
                
                        end
        
                    end
        
                end
            
            end
        
        end
    
    end

end

disp('Analysis done');
