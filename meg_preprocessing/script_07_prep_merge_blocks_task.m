%% preprocessing MEG data for WM task [GERE project]

% concatenate meg sensor data across runs within session

% path to data
path_results = ['/path_to_local/results/preprocessing/' folder];
path_preprocessing = '/path_to_local/results/preprocessing/triggers';

%% loop over trials

trialevent = 'wholetrial_baseline1sec';
zthr = 25;

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

        if ~isfile([path_results '/' subject_ID '/sess_0' num2str(session)...
                    '/GERE_task_allblocks' ...
                    '_' trialevent '_validtrialszthr'...
                    num2str(zthr) '.mat'])

            data_validtrials_all = [];
                                    
            for block_i = 1:3
    
                data_locked_file = [path_results '/' subject_ID '/' session_ID ...
                    '/GERE_task_block' num2str(block_i) ...
                    '_' trialevent '_validtrialszthr'...
                    num2str(zthr) '.mat'];

                if isfile(data_locked_file)
    
                    % load data
                    load(data_locked_file); % data_validtrials
        
                    if ~isempty(data_validtrials)
                                    
                        if isempty(data_validtrials_all)
        
                            data_validtrials_all = data_validtrials;
        
                        else
        
                            last_trial = size(data_validtrials_all.trial, 2);
        
                            for trial_i = 1:size(data_validtrials.trial, 2)
        
                                data_validtrials_all.time{last_trial + trial_i} = data_validtrials.time{trial_i};
                                data_validtrials_all.trial{last_trial + trial_i} = data_validtrials.trial{trial_i};
                                data_validtrials_all.sampleinfo(last_trial + trial_i,:) = data_validtrials.sampleinfo(trial_i,:);
                                data_validtrials_all.trialinfo(last_trial + trial_i,:) = data_validtrials.trialinfo(trial_i,:);
        
                            end
        
                        end
        
                    end

                end
    
            end
    
            data_task = data_validtrials_all;

            if ~isempty(data_task)
    
                save([path_results '/' subject_ID '/sess_0' num2str(session)...
                        '/GERE_task_allblocks' ...
                        '_' trialevent '_validtrialszthr'...
                        num2str(zthr) '.mat'], 'data_task', '-v7.3');

            end

        end
            
    end

end
