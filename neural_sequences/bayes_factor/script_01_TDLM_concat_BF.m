%% GERE project
% replay

if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = 1:17;
end

path_data = [path_root, '/github/results/replay/data'];
addpath([path_root, '/github/scripts/replay/utilities'])

classifier = 'l2_e-1';
nlags = 100;
number_stimuli = 8;
sessions = 1:2;

seq_lengths = [3 4];

replay_folders = {'sequence_concat_delay', 'sequence_concat_stim'};

%% bayes factor

% loop over replay approaches
for folder_replay_i = 1:length(replay_folders)

    folder_replay = replay_folders{folder_replay_i};
    path_results = [path_root, '/github/results/replay/' folder_replay];

    if strcmp(folder_replay, 'sequence_concat_delay')
        folders = {'task_stim2delay' 'localizer2delay'};
    else
        folders = {'task_stim2task_stim' 'localizer2task_stim'};
    end

    % loop over folders
    for folder_i = 1:length(folders)
    
        path_outputs = [path_root, '/github/results/replay/' folder_replay '/group_results/' classifier '/' folders{folder_i} '/bayes_factor'];
    
        if ~isfolder(path_outputs)
            mkdir(path_outputs);
        end
    
        %% compute empirical prior from surrogate data across lengths and lags
        
        % ---------------------------------------------------------
        % collect pooled null values
        % ---------------------------------------------------------
        null_all = [];
        
        for length_i = seq_lengths
        
            % -----------------------------------------------------
            % forward / backward
            % -----------------------------------------------------
            for fwd_bwd_i = 1:2
        
                if fwd_bwd_i == 1
                    seq_direction = 'f';
                else
                    seq_direction = 'b';
                end
        
                % -------------------------------------------------
                % subjects
                % -------------------------------------------------
                for sub_i = 1:length(subjects)
        
                    subject = subjects(sub_i);
        
                    % -------------------------------------------------
                    % subject ID
                    % -------------------------------------------------
                    if subject < 10
                        subject_ID = ['sub_0' num2str(subject)];
                    else
                        subject_ID = ['sub_' num2str(subject)];
                    end
        
                    % -------------------------------------------------
                    % paths
                    % -------------------------------------------------
                    outdir = [path_results '/' subject_ID '/' ...
                              classifier '/' folders{folder_i}];
        
                    % -------------------------------------------------
                    % load sequenceness
                    % -------------------------------------------------
                    load([outdir '/Ms' seq_direction '_l' ...
                          num2str(length_i) '.mat']);
        
                    % -------------------------------------------------
                    % select variable
                    % -------------------------------------------------
                    if seq_direction == 'f'
        
                        if length_i == 3
                            sequenceness = Msf_l3;
                        else
                            sequenceness = Msf_l4;
                        end
        
                    else
        
                        if length_i == 3
                            sequenceness = Msb_l3;
                        else
                            sequenceness = Msb_l4;
                        end
        
                    end
        
                    % -------------------------------------------------
                    % surrogate/null values
                    %
                    % size:
                    %   [nPerms-1 x nLags]
                    % -------------------------------------------------
                    null_vals = squeeze( ...
                        sequenceness(1,2:end,1:nlags));
        
                    % -------------------------------------------------
                    % append to pooled null
                    %
                    % reshape into vector
                    % -------------------------------------------------
                    null_all = [null_all; null_vals(:)];
        
                end
            end
        end
    
        % ---------------------------------------------------------
        % GLOBAL null SD
        % ---------------------------------------------------------
        global_null_sd = std(null_all);
        
        %% bayes factor
    
        % loop over sequence lengths
        for length_i = seq_lengths
    
            % loop over fwd or bwd sequenceness
            for fwd_bwd_i = 1:2
    
                if fwd_bwd_i == 1
                    seq_direction = 'f';
                    seq_name = 'forward';
                else
                    seq_direction = 'b';
                    seq_name = 'backward';
                end
                
                % empirical sequenceness
                sequenceness_allsubjects = nan(length(subjects), nlags);
                                        
                for sub_i = 1:length(subjects)
                
                    subject = subjects(sub_i);
                
                    % -----------------------------------------------------
                    % subject ID
                    % -----------------------------------------------------
                    if subject < 10
                        subject_ID = ['sub_0' num2str(subject)];
                        subjectID = ['sub0' num2str(subject)];
                    else
                        subject_ID = ['sub_' num2str(subject)];
                        subjectID = ['sub' num2str(subject)];
                    end
                            
                    % -----------------------------------------------------
                    % paths
                    % -----------------------------------------------------
                    outdir = [path_results '/' subject_ID '/' classifier '/' folders{folder_i}];
                
                    % -----------------------------------------------------
                    % load TDLM results
                    % -----------------------------------------------------
                    load([outdir '/Ms' seq_direction '_l' num2str(length_i) '.mat']);   % variable
    
                    if seq_direction == 'f'
                        if length_i == 3 
                            sequenceness = Msf_l3;
                        else
                            sequenceness = Msf_l4;
                        end
                    else
                        if length_i == 3 
                            sequenceness = Msb_l3;
                        else
                            sequenceness = Msb_l4;
                        end
                    end
    
                    % -----------------------------------------------------
                    % empirical sequenceness
                    %
                    % Msf_l3(1,1,:) = empirical
                    % -----------------------------------------------------
                    empirical = squeeze(sequenceness(1,1,1:nlags));
                
                    sequenceness_allsubjects(sub_i,:) = empirical;
                            
                end
                
                        
                %% BF across lags
    
                bf_lags = nan(1, nlags);
    
                for lag_i = 1:nlags
                        
                    % BF with empirical scale (std of null distribution)
                    [bf_lags(lag_i), ~] = bf.ttest(sequenceness_allsubjects(:,lag_i), 'scale', global_null_sd, 'tail', 'right');
    
                end

                bf10 = bf_lags;
                bf01 = 1./bf_lags;
    
                % save
                save([path_outputs '/bf10_lags_' seq_name '_length' num2str(length_i) '.mat'], 'bf10');
                save([path_outputs '/bf01_lags_' seq_name '_length' num2str(length_i) '.mat'], 'bf01');
    
            end
    
        end
    
    end

end

disp('Analysis done');

