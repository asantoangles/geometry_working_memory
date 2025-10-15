%% Extracting condition-specific trials

% this is taking 2 seconds as baseline, instead of 2.2 seconds as coded in ptb, 
% discarding the first 200 ms in the beginning, because of overlap with
% motor response to initiate trial
% then, trial finishes at response onset
% see triggers_trial_events.info about trial events

% from ft_definetrial.m
% The trial definition "trl" is an Nx3 matrix, N is the number of trials.
% The first column contains the sample-indices of the begin of each trial
% relative to the begin of the raw data, the second column contains the
% sample-indices of the end of each trial, and the third column contains
% the offset of the trigger with respect to the trial. An offset of 0
% means that the first sample of the trial corresponds to the trigger. A
% positive offset indicates that the first sample is later than the trigger,
% a negative offset indicates that the trial begins before the trigger.

% path to outputs
path_outputs = ['/path_to_local/results/preprocessing/' folder];
path_preprocessing = '/path_to_local/results/preprocessing/triggers';

if ~isfolder(path_outputs)
    mkdir(path_outputs);
end

if ~isfolder(path_preprocessing)
    mkdir(path_preprocessing);
end


% path to data
path_inputs = '/path_to_local/data/MEG';

%% loop over subjects

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

        % outputs folder
        if ~isfolder([path_outputs '/' subject_ID '/' session_ID])
            mkdir([path_outputs '/' subject_ID '/' session_ID]);
        end

        if ~isfolder([path_preprocessing '/' subject_ID '/' session_ID])
            mkdir([path_preprocessing '/' subject_ID '/' session_ID]);
        end

        % copy renamed file
        meg_con_original = [path_inputs '/' subject_ID '/sess_0' num2str(session) '/GERE_' subjectID '_sess' num2str(session) '*_01_analysis_01.con'];
        meg_data = [path_preprocessing '/' subject_ID '/sess_0' num2str(session) '/GERE_' subjectID '_sess' num2str(session) '_localizer.con'];
        
        if ~isfile(meg_data)
            system(['cp ' meg_con_original ' ' meg_data]);
        end

        if isfile(meg_data)
                            
            %% segment continuous data of localizer into trials 
            % each trial is defined as the stimulus presentation
                    
            % load stimulus timings
            load([path_preprocessing '/' subject_ID '/' session_ID '/triggers_localizer.mat']); % triggers_time
            
            % order triggers
            triggers_localizer = nan(2, size(triggers_time.trigger_location_1, 2)*8);
            
            triggers_localizer(1,:) = [triggers_time.trigger_location_1 triggers_time.trigger_location_2...
                                       triggers_time.trigger_location_3 triggers_time.trigger_location_4...
                                       triggers_time.trigger_location_5 triggers_time.trigger_location_6...
                                       triggers_time.trigger_location_7 triggers_time.trigger_location_8];
            
            triggers_localizer(2,:) = [ones(1, size(triggers_time.trigger_location_1, 2)) * 1 ...
                                       ones(1, size(triggers_time.trigger_location_1, 2)) * 2 ...
                                       ones(1, size(triggers_time.trigger_location_1, 2)) * 3 ...
                                       ones(1, size(triggers_time.trigger_location_1, 2)) * 4 ...
                                       ones(1, size(triggers_time.trigger_location_1, 2)) * 5 ...
                                       ones(1, size(triggers_time.trigger_location_1, 2)) * 6 ...
                                       ones(1, size(triggers_time.trigger_location_1, 2)) * 7 ...
                                       ones(1, size(triggers_time.trigger_location_1, 2)) * 8];
            
            triggers_localizer = sortrows(triggers_localizer')';

            data_locked_file = [path_outputs '/' subject_ID '/' session_ID '/GERE_localizer.mat'];
            
            %%% each stimulus is one trial - stimulus present (400 ms)
            if ~isfile(data_locked_file)
                
                % create trl
                number_trials = size(triggers_localizer, 2);
                trl_events = zeros(number_trials, 3);
                trl_events(:, 1) = triggers_localizer(1, :); % stimulus onset
                trl_events(:, 2) = triggers_localizer(1, :) + 399; % stimulus offset
                
                % define trials based on trl, and demean 
                cfg             = [];
                cfg.dataset     = meg_data;
                cfg.trl         = trl_events;
                cfg.demean      = 'yes';    
                cfg.channel     = {'AG*'};  
                data_locked     = ft_preprocessing(cfg);
                
                % save
                save(data_locked_file, 'data_locked', '-v7.3');

            end            

            data_locked_file = [path_outputs '/' subject_ID '/' session_ID '/GERE_localizer_isi.mat'];
            
            %%% each isi is one trial - stimulus present (400 ms)
            if ~isfile(data_locked_file)
                
                % create trl
                number_trials = size(triggers_localizer, 2);
                trl_events = zeros(number_trials, 3);
                trl_events(:, 1) = triggers_localizer(1, :) + 399; % isi onset
                trl_events(:, 2) = triggers_localizer(1, :) + 598; % isi offset
                
                % define trials based on trl, and demean 
                cfg             = [];
                cfg.dataset     = meg_data;
                cfg.trl         = trl_events;
                cfg.demean      = 'yes';    
                cfg.channel     = {'AG*'};  
                data_locked     = ft_preprocessing(cfg);
                
                % save
                save(data_locked_file, 'data_locked', '-v7.3');

            end            

            data_locked_file = [path_outputs '/' subject_ID '/' session_ID '/GERE_localizer_stimisi.mat'];
            
            %%% each isi is one trial - stimulus present (400 ms)
            if ~isfile(data_locked_file)
                
                % create trl
                number_trials = size(triggers_localizer, 2);
                trl_events = zeros(number_trials, 3);
                trl_events(:, 1) = triggers_localizer(1, :); % stimulus onset
                trl_events(:, 2) = triggers_localizer(1, :) + 599; % isi offset
                
                % define trials based on trl, and demean 
                cfg             = [];
                cfg.dataset     = meg_data;
                cfg.trl         = trl_events;
                cfg.demean      = 'yes';    
                cfg.channel     = {'AG*'};  
                data_locked     = ft_preprocessing(cfg);
                
                % save
                save(data_locked_file, 'data_locked', '-v7.3');

            end            

            
            %% segment continuous data of WM task into trials
            
            number_blocks = 3;

            % exceptions
            if subject == 20
                if session == 1
                    number_blocks = 2;
                end
            end

            % exception
            if subject == 51
                if session == 2
                    number_blocks = 2;
                end
            end
            
            for block_i = 1:number_blocks

                meg_con_original = [path_inputs '/' subject_ID '/sess_0' num2str(session) '/GERE_' subjectID '_sess' num2str(session) '*_0' num2str(block_i+1) '_analysis_01.con'];
                
                % copy renamed file
                meg_data = [path_preprocessing '/' subject_ID '/sess_0' num2str(session) '/GERE_' subjectID '_sess' num2str(session) '_block' num2str(block_i) '.con'];
                if ~isfile(meg_data)
                    system(['cp ' meg_con_original ' ' meg_data]);
                end
                            
                % load trial timings
                load([path_preprocessing '/' subject_ID '/' session_ID '/triggers_task_block' num2str(block_i) '.mat']); % triggers_trial_events

                %% wholetrial - 1 second baseline

                filename = 'wholetrial_baseline1sec';

                data_locked_file = [path_outputs '/' subject_ID '/' session_ID ...
                    '/GERE_task_block' num2str(block_i) '_' filename '.mat'];

                if ~isfile(data_locked_file)
                
                    % create trl
                    number_trials = size(triggers_trial_events.data, 1);
                    trl_events = zeros(number_trials, 3);
                    trl_events(:, 1) = triggers_trial_events.data(:, 1) + 1200; % baseline (1 second)
                    trl_events(:, 2) = triggers_trial_events.data(:, 4); % response onset
                    
                    % define trials based on trl, and demean 
                    cfg             = [];
                    cfg.dataset     = meg_data;
                    cfg.trl         = trl_events;
                    cfg.demean      = 'yes';    
                    cfg.channel     = {'AG*'};  
                    data_locked     = ft_preprocessing(cfg);
                
                    % save
                    save(data_locked_file, 'data_locked', '-v7.3');

                end
            
                %% baseline - 1 second
            
                filename = 'baseline_1sec';

                data_locked_file = [path_outputs '/' subject_ID '/' session_ID ...
                    '/GERE_task_block' num2str(block_i) '_' filename '.mat'];

                if ~isfile(data_locked_file)

                    % create trl
                    number_trials = size(triggers_trial_events.data, 1);
                    trl_events = zeros(number_trials, 3);
                    trl_events(:, 1) = triggers_trial_events.data(:, 1) + 1200; % baseline (2 seconds)
                    trl_events(:, 2) = triggers_trial_events.data(:, 2); % stimulus onset
                    
                    % define trials based on trl, and demean 
                    cfg             = [];
                    cfg.dataset     = meg_data;
                    cfg.trl         = trl_events;
                    cfg.demean      = 'yes';    
                    cfg.channel     = {'AG*'};  
                    data_locked     = ft_preprocessing(cfg);
                
                    % save
                    save(data_locked_file, 'data_locked', '-v7.3');

                end
            
                %% delay
            
                filename = 'delay';

                data_locked_file = [path_outputs '/' subject_ID '/' session_ID ...
                    '/GERE_task_block' num2str(block_i) '_' filename '.mat'];
            
                if ~isfile(data_locked_file)
                
                    % create trl
                    number_trials = size(triggers_trial_events.data, 1);
                    trl_events = zeros(number_trials, 3);
                    trl_events(:, 1) = triggers_trial_events.data(:, 3); % delay onset
                    trl_events(:, 2) = triggers_trial_events.data(:, 4); % response onset
                    
                    % define trials based on trl, and demean 
                    cfg             = [];
                    cfg.dataset     = meg_data;
                    cfg.trl         = trl_events;
                    cfg.demean      = 'yes';    
                    cfg.channel     = {'AG*'};  
                    data_locked     = ft_preprocessing(cfg);
                
                    % save
                    save(data_locked_file, 'data_locked', '-v7.3');

                end
            
            end

        end

    end

end
