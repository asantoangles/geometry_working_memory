%% automatic artifact detection (trial-wise)
% https://www.fieldtriptoolbox.org/tutorial/automatic_artifact_rejection/

% loop over z-score threshold of artifact detection to select the appropiate
% trade-off between clean data and statistical power

% path to outputs
path_results = ['/path_to_local/results/preprocessing/' folder];

%% loop over trials

threshold_muscle_blinks = [25];

for sub_i = 1:length(subjects)

    subject = subjects(sub_i);

    for ses_i = 1:length(sessions)

        session = sessions(ses_i);

        % set paths
        if subject < 10
            subject_ID = ['sub_0' num2str(subject)];
            subjectID = ['sub0' num2str(subject)];
        else
            subject_ID = ['sub_' num2str(subject)];
            subjectID = ['sub' num2str(subject)];
        end
        
        session_ID = ['sess_0' num2str(session)];

        disp([subject_ID ' - ' session_ID]);

        clear survivor_trials_table

        % load output if exists (for the if statement below)
        if isfile([path_results '/' subject_ID '/' session_ID '/GERE_task_artifact_thr.mat'])
            load([path_results '/' subject_ID '/' session_ID '/GERE_task_artifact_thr.mat'])
        end

        % run if input exists
        if isfile([path_results '/' subject_ID '/' session_ID '/GERE_task_block1_wholetrial.mat'])

            % run if output does not exist
            if ~isfile([path_results '/' subject_ID '/' session_ID '/GERE_task_artifact_thr_20to30.mat'])
                
                %% settings loops
                
                trial_events = {'wholetrial' 'wholetrial_baseline1sec' 'delay' 'baseline_1sec'};
                
                %% loops
                
                % emtpy matrix for results
                rownames        = {'wholetrial - block 1',...
                                   'wholetrial - block 2',...    
                                   'wholetrial - block 3',...    
                                   'wholetrial_baseline1sec - block 1',...
                                   'wholetrial_baseline1sec - block 2',...    
                                   'wholetrial_baseline1sec - block 3',... 
                                   'delay - block 1',...    
                                   'delay - block 2',...    
                                   'delay - block 3',...    
                                   'baseline_1sec - block 1',...    
                                   'baseline_1sec - block 2',...    
                                   'baseline_1sec - block 3'};
                colnames = {'20' '25' '30'};
                
                survivor_trials = nan(length(rownames), length(colnames));
                
                % loop over trial events
                for trialevent_i = 1:length(trial_events)
                
                    % loop over blocks
                    number_blocks = 3;
                    for block_i = 1:number_blocks

                        data_locked_file = [path_results '/' subject_ID '/' session_ID '/GERE_task_'...
                         'block' num2str(block_i) '_' trial_events{trialevent_i} '.mat'];

                        if isfile(data_locked_file)
                    
                            % set location on output
                            if trialevent_i == 1
                                location = block_i;
                            elseif trialevent_i == 2
                                location = block_i + 3;
                            elseif trialevent_i == 3
                                location = block_i + 6;
                            elseif trialevent_i == 4
                                location = block_i + 9;
                            end
        
                            % loop over thresholds
                            for thr_i = 1:length(threshold_muscle_blinks)

                                try

                                    load([path_results '/' subject_ID '/' session_ID '/GERE_task_block' num2str(block_i) ... 
                                        '_no_artifacts_zthr' num2str(threshold_muscle_blinks(thr_i)) '_' trial_events{trialevent_i} '.mat']);

                                    if ~isempty(data_no_artifacts)
                                        data_no_artifacts.trialinfo(1,1) + 10;
                                    end

                                catch

                                    % load data
                                    load(data_locked_file); % data_locked
                        
                                    %% Jump artifact detection
                                    
                                    % channel selection, cutoff and padding
                                    cfg                                 = [];
                                    cfg.artfctdef.zvalue.channel        = 'AG*';
                                    cfg.artfctdef.zvalue.cutoff         = threshold_muscle_blinks(thr_i);
                                    cfg.artfctdef.zvalue.trlpadding     = 0;
                                    cfg.artfctdef.zvalue.artpadding     = 0;
                                    cfg.artfctdef.zvalue.fltpadding     = 0;
                                    
                                    % algorithmic parameters
                                    cfg.artfctdef.zvalue.cumulative     = 'yes';
                                    cfg.artfctdef.zvalue.medianfilter   = 'yes';
                                    cfg.artfctdef.zvalue.medianfiltord  = 9;
                                    cfg.artfctdef.zvalue.absdiff        = 'yes';
                                    
                                    [~, artifact_jump] = ft_artifact_zvalue(cfg, data_locked);
                                    
                                    %% Detection of muscle artifacts
                                    
                                    % channel selection, cutoff and padding
                                    cfg                               = [];
                                    cfg.artfctdef.zvalue.channel      = 'AG*';
                                    cfg.artfctdef.zvalue.cutoff       = threshold_muscle_blinks(thr_i);
                                    cfg.artfctdef.zvalue.trlpadding   = 0;
                                    cfg.artfctdef.zvalue.fltpadding   = 0;
                                    cfg.artfctdef.zvalue.artpadding   = 0.1;
                                    
                                    % algorithmic parameters
                                    cfg.artfctdef.zvalue.bpfilter     = 'yes';
                                    cfg.artfctdef.zvalue.bpfreq       = [110 140];
                                    cfg.artfctdef.zvalue.bpfiltord    = 9;
                                    cfg.artfctdef.zvalue.bpfilttype   = 'but';
                                    cfg.artfctdef.zvalue.hilbert      = 'yes';
                                    cfg.artfctdef.zvalue.boxcar       = 0.2;
                                    
                                    [~, artifact_muscle] = ft_artifact_zvalue(cfg, data_locked);
                                    
                                    %% Detection of EOG artifacts
                                    
                                    % channel selection, cutoff and padding
                                    cfg                              = [];
                                    cfg.artfctdef.zvalue.channel     = {'AG087' 'AG130'};
                                    cfg.artfctdef.zvalue.cutoff      = threshold_muscle_blinks(thr_i);
                                    cfg.artfctdef.zvalue.trlpadding  = 0;
                                    cfg.artfctdef.zvalue.artpadding  = 0.1;
                                    cfg.artfctdef.zvalue.fltpadding  = 0;
                                    
                                    % algorithmic parameters
                                    cfg.artfctdef.zvalue.bpfilter   = 'yes';
                                    cfg.artfctdef.zvalue.bpfilttype = 'but';
                                    cfg.artfctdef.zvalue.bpfreq     = [2 15];
                                    cfg.artfctdef.zvalue.bpfiltord  = 4;
                                    cfg.artfctdef.zvalue.hilbert    = 'yes';
                                    
                                    [~, artifact_eog] = ft_artifact_zvalue(cfg, data_locked);
                        
                                    %% artifact removal (just to count removed trials)
            
                                    try
            
                                        cfg                           = [];
                                        cfg.artfctdef.reject          = 'complete';
                                        cfg.artfctdef.eog.artifact    = artifact_eog;
                                        cfg.artfctdef.jump.artifact   = artifact_jump;
                                        cfg.artfctdef.muscle.artifact = artifact_muscle;
                                        data_no_artifacts = ft_rejectartifact(cfg, data_locked);

                                    catch

                                        data_no_artifacts = [];
                                        
                                    end
            
                                    % save
                                    save([path_results '/' subject_ID '/' session_ID '/GERE_task_block' num2str(block_i) ... 
                                        '_no_artifacts_zthr' num2str(threshold_muscle_blinks(thr_i)) '_' trial_events{trialevent_i} '.mat'], 'data_no_artifacts');
                    
                                end

                                if ~isempty(data_no_artifacts)

                                    % survival trials
                                    survivor_trials(location, thr_i) = size(data_no_artifacts.trial, 2);

                                else

                                    survivor_trials(location, thr_i) = 0;

                                end
                                
                            
                            end
                    
                        end

                    end
                
                end
                
                % add row/column names
                survivor_trials_table = array2table(survivor_trials, 'RowNames', rownames, 'VariableNames', colnames);
        
                % save
                save([path_results '/' subject_ID '/' session_ID '/GERE_task_artifact_thr_20to30.mat'], 'survivor_trials_table');
    
            end

        end

    end

end



