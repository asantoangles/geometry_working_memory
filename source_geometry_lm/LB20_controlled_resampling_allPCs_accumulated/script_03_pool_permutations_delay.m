%% GERE project

% path to data
if isfolder('/path_to_local')
    path_inputs = '/path_to_local/results/source_reconstruction';
    path_results = ['/path_to_local/results/source_geometry_lm/' folder_geometry '/comp1to' num2str(last_comp_i)];
    addpath('/path_to_local/scripts/source_geometry_lm/utilities')
    functions_data ={@mean};
end

% time windows
segments_start = 1:300:3701;
segments_end = 300:300:4000;
delay_segments = [segments_start; segments_end]';

events = {'delay'};
performance = {'correct_trials', 'incorrect_trials'};
iterations_bootstrap = 1000;
permutations = 1000;
sessions = 1:2;
number_locations = 8;


%% pool permutations delay 

for last_comp_i = [6 8]

    path_results = ['/path_to_local/results/source_geometry_lm/' folder_geometry '/comp1to' num2str(last_comp_i)];

    folder_original = 'LB20_allPCs_accumulated';
    path_results_original = ['/path_to_local/results/source_geometry_lm/' folder_original '/comp1to' num2str(last_comp_i)];

    for event_i = 1:length(events)
    
        if strcmp(events{event_i}, 'delay')
    
            time_segments = delay_segments;
    
        elseif strcmp(events{event_i}, 'baseline')
    
            time_segments = baseline_segments;
    
        elseif strcmp(events{event_i}, 'stim_resolved')
            
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
    
                            if ~isfolder([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))])
                                mkdir([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))]);
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
        
                                    for rank_i = 1:ranks
        
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
        
                                for rank_i = 1:ranks
        
                                    separability.volume{rank_i} = mean(separability.volume{rank_i});
                                    separability.distance_by_separation{rank_i} = mean(separability.distance_by_separation{rank_i},1);
        
                                end
        
                                save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_' sequence_length_filename '.mat'], 'separability');
        
                                %%% separability all pairs
        
                                old_file = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_alpairs_perm_' sequence_length_filename '.mat'];
                                new_file = [path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_allpairs_perm_' sequence_length_filename '.mat'];
        
                                if isfile(old_file)
                                    system(['mv ' old_file ' ' new_file]);
                                end
        
                                load([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_allpairs_perm_' sequence_length_filename '.mat']); % separability_allpairs_perm
        
        
                                separability_allpairs = [];
                                separability_allpairs.distance_by_separation = separability_allpairs_perm{1};
                                for rank_i = 1:ranks
                                    separability_allpairs.distance_by_separation{rank_i} = nan(number_locations, number_locations,permutations);
                                end
        
                                for perm_i = 1:permutations
        
                                    for rank_i = 1:ranks
        
                                        separability_allpairs.distance_by_separation{rank_i}(:,:,perm_i) = separability_allpairs_perm{perm_i}{rank_i};
                                
                                    end
        
                                end
        
                                for rank_i = 1:ranks
                                    separability_allpairs.distance_by_separation{rank_i} = squeeze(mean(separability_allpairs.distance_by_separation{rank_i},3));
                                end
        
                                save([path_results '/' subject_ID '/' func2str(fun_i) '/' performance{perf_i}...
                                    '/' events{event_i} '_' num2str(time_segments(delay_i,1)) 'to' num2str(time_segments(delay_i,2))...
                                    '/separability_allpairs_' sequence_length_filename '.mat'], 'separability_allpairs');
        
        
                                % transfer remaining files from C1 to C1_balance
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
