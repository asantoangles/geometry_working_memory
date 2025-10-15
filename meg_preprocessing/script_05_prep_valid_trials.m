%% count outputs

% path to outputs
path_results = ['/path_to_local/results/preprocessing/' folder];

%% subset data_locked with valid trials

zthreshold = [25];

for zthr_i = 1:length(zthreshold)

    zthr = zthreshold(zthr_i);

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
    
            if isfile([path_results '/' subject_ID '/' session_ID ...
                        '/GERE_localizer_rejected_artifacts_zthr' num2str(zthr) '.mat'])
    
                %% localizer
                % create data_loscked with valid trials
                % excluding trials with dim/color stimulation and response
        
                if ~isfile([path_results '/' subject_ID '/' session_ID ...
                    '/GERE_localizer_validtrialszthr' num2str(zthr) '.mat'])
        
                    % load data_no_artifacts
                    load([path_results '/' subject_ID '/' session_ID ...
                        '/GERE_localizer_rejected_artifacts_zthr' num2str(zthr) '.mat']); % data_no_artifacts
                    
                    % load data_locked
                    load([path_results '/' subject_ID '/' session_ID ...
                        '/GERE_localizer.mat']); % data_locked
            
                    % find valid trials
                    subset_trials = [];
                    for trial_i = 1:size(data_locked.trialinfo, 1)
                        if data_locked.trialinfo(trial_i, 2) == 1
                            if ismember(trial_i, data_no_artifacts.trialinfo(:,1))
                                subset_trials = [subset_trials trial_i];
                            end
                        end
                    end
                   
                    % subset data_locked with valid trials from data_no_artifacts
                    data_locked.trial = data_locked.trial(1,subset_trials);
                    data_locked.time = data_locked.time(1,subset_trials);
                    data_locked.sampleinfo = data_locked.sampleinfo(subset_trials,:);
                    data_locked.trialinfo = data_locked.trialinfo(subset_trials,:);
        
                    data_validtrials = data_locked;
            
                    % save
                    save([path_results '/' subject_ID '/' session_ID ...
                        '/GERE_localizer_validtrialszthr' num2str(zthr) '.mat'], 'data_validtrials');
        
                end
                
                %% task - wholetrial_baseline1sec
                % create data_locked with valid trials
                % (including incorrect responses)
                % from 'zthr 10 during delay', but
                % including wholetrial_baseline1sec, as the input for ICA
        
                for block_i = 1:3
        
                    if ~isfile([path_results '/' subject_ID '/' session_ID ...
                                    '/GERE_task_block' num2str(block_i) '_wholetrial_baseline1sec_validtrialszthr'...
                                    num2str(zthr) '.mat'])
        
                        if isfile([path_results '/' subject_ID '/' session_ID ...
                            '/GERE_task_block' num2str(block_i) '_wholetrial_baseline1sec.mat'])
                
                            % load data_locked
                            load([path_results '/' subject_ID '/' session_ID ...
                                '/GERE_task_block' num2str(block_i) '_wholetrial_baseline1sec.mat']); % data_locked
        
                            % load data_no_artifacts
                            load([path_results '/' subject_ID '/' session_ID ...
                                '/GERE_task_block' num2str(block_i) '_no_artifacts_zthr' num2str(zthr) '_delay.mat']); % data_no_artifacts
            
    %                         if isfield(data_no_artifacts, 'trialinfo')
                
                                % find valid trials
                                subset_trials = [];
                                for trial_i = 1:size(data_locked.trialinfo, 1)
                                    if data_locked.trialinfo(trial_i, 2) == 1
                                        if ismember(trial_i, data_no_artifacts.trialinfo(:,1))
                                            subset_trials = [subset_trials trial_i];
                                        end
                                    end
                                end
                               
                                % subset data_locked with valid trials from data_no_artifacts
                                data_locked.trial = data_locked.trial(1,subset_trials);
                                data_locked.time = data_locked.time(1,subset_trials);
                                data_locked.sampleinfo = data_locked.sampleinfo(subset_trials,:);
                                data_locked.trialinfo = data_locked.trialinfo(subset_trials,:);
        
                                data_validtrials = data_locked;
    
                                %%% add info about baseline without artifacts
                                
                                % load data from baseline (1 second)
                                load([path_results '/' subject_ID '/' session_ID ...
                                    '/GERE_task_block' num2str(block_i) '_no_artifacts_zthr'...
                                    num2str(zthr) '_baseline_1sec.mat']); % data_no_artifacts
                
                                % legend
                                data_validtrials.trialinfolegend{12, 1} = 'valid baseline 1 second';
                                data_validtrials.trialinfo(:,12) = zeros(size(data_validtrials.trialinfo, 1), 1);
                    
                                for trial_i = 1:size(data_validtrials.trial, 2)
                        
                                    if ismember(data_validtrials.trialinfo(trial_i, 1), data_no_artifacts.trialinfo(:,1))
                    
                                        data_validtrials.trialinfo(trial_i,12) = 1;
                    
                                    end
                    
                                end
                                
                                % save
                                save([path_results '/' subject_ID '/' session_ID ...
                                    '/GERE_task_block' num2str(block_i) '_wholetrial_baseline1sec_validtrialszthr'...
                                    num2str(zthr) '.mat'], 'data_validtrials');
            
    %                         end
        
                        end
        
                    end
        
                end
        
            end
    
        end
    
    end

end