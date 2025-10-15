%% additional information in data_locked
% add valid trials and correct responses in sampleinfo field

% path to outputs
path_outputs = ['/path_to_local/results/preprocessing/' folder];
path_preprocessing = '/path_to_local/results/preprocessing/triggers';
path_behavior = '/path_to_local/data/behavior';

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
        sessionID = ['sess0' num2str(session)];

        disp([subject_ID ' - ' session_ID]); 

        data_locked_file = [path_outputs '/' subject_ID '/' session_ID '/GERE_localizer.mat'];
        
        if isfile(data_locked_file)

            load(data_locked_file);

            if ~isfield(data_locked, 'trialinfo')

                % reset sampleinfo
                data_locked.sampleinfo = data_locked.sampleinfo(:,1:2);
                if isfield(data_locked, 'sampleinfolegend')
                    data_locked = rmfield(data_locked, 'sampleinfolegend');
                end
    
                % new field with trial information
                data_locked.trialinfo = nan(size(data_locked.sampleinfo, 1), 2);
                data_locked.trialinfolegend = {'trial number'; 'valid trials (no dim/color task or response)'; 'stimulus'};
    
                % trial number
                data_locked.trialinfo(:, 1) = 1:(size(data_locked.trialinfo, 1));
    
                % valid trials (exclude responses and changing color)
                load([path_behavior '/' subject_ID '/' session_ID '/functional_localizer.mat']); % functional_localizer
                load([path_behavior '/' subject_ID '/' session_ID '/responses_localizer.mat']); % responses
                tmp = ones(size(data_locked.trialinfo, 1), 1);
                tmp = tmp - responses';
                tmp = tmp - functional_localizer.data(2,:)';
                tmp = tmp - functional_localizer.data(3,:)';
                tmp(tmp < 0) = 0;
                data_locked.trialinfo(:, 2) = tmp;
                clear tmp;
    
                % stimulus location
                data_locked.trialinfo(:, 3) = functional_localizer.data(1, :);
    
                % save
                save(data_locked_file, 'data_locked', '-v7.3');

            end

        end            


        %%% localizer isi
        data_locked_file = [path_outputs '/' subject_ID '/' session_ID '/GERE_localizer_isi.mat'];
        
        if isfile(data_locked_file)

            load(data_locked_file);

            if ~isfield(data_locked, 'trialinfo')

                % reset sampleinfo
                data_locked.sampleinfo = data_locked.sampleinfo(:,1:2);
                if isfield(data_locked, 'sampleinfolegend')
                    data_locked = rmfield(data_locked, 'sampleinfolegend');
                end
    
                % new field with trial information
                data_locked.trialinfo = nan(size(data_locked.sampleinfo, 1), 2);
                data_locked.trialinfolegend = {'trial number'; 'valid trials (no dim/color task or response)'; 'stimulus'};
    
                % trial number
                data_locked.trialinfo(:, 1) = 1:(size(data_locked.trialinfo, 1));
    
                % valid trials (exclude responses and changing color)
                load([path_behavior '/' subject_ID '/' session_ID '/functional_localizer.mat']); % functional_localizer
                load([path_behavior '/' subject_ID '/' session_ID '/responses_localizer.mat']); % responses
                tmp = ones(size(data_locked.trialinfo, 1), 1);
                tmp = tmp - responses';
                tmp = tmp - functional_localizer.data(2,:)';
                tmp = tmp - functional_localizer.data(3,:)';
                tmp(tmp < 0) = 0;
                data_locked.trialinfo(:, 2) = tmp;
                clear tmp;
    
                % stimulus location
                data_locked.trialinfo(:, 3) = functional_localizer.data(1, :);
    
                % save
                save(data_locked_file, 'data_locked', '-v7.3');

            end

        end   


        %%% localizer stimisi
        data_locked_file = [path_outputs '/' subject_ID '/' session_ID '/GERE_localizer_stimisi.mat'];
        
        if isfile(data_locked_file)

            load(data_locked_file);

            if ~isfield(data_locked, 'trialinfo')

                % reset sampleinfo
                data_locked.sampleinfo = data_locked.sampleinfo(:,1:2);
                if isfield(data_locked, 'sampleinfolegend')
                    data_locked = rmfield(data_locked, 'sampleinfolegend');
                end
    
                % new field with trial information
                data_locked.trialinfo = nan(size(data_locked.sampleinfo, 1), 2);
                data_locked.trialinfolegend = {'trial number'; 'valid trials (no dim/color task or response)'; 'stimulus'};
    
                % trial number
                data_locked.trialinfo(:, 1) = 1:(size(data_locked.trialinfo, 1));
    
                % valid trials (exclude responses and changing color)
                load([path_behavior '/' subject_ID '/' session_ID '/functional_localizer.mat']); % functional_localizer
                load([path_behavior '/' subject_ID '/' session_ID '/responses_localizer.mat']); % responses
                tmp = ones(size(data_locked.trialinfo, 1), 1);
                tmp = tmp - responses';
                tmp = tmp - functional_localizer.data(2,:)';
                tmp = tmp - functional_localizer.data(3,:)';
                tmp(tmp < 0) = 0;
                data_locked.trialinfo(:, 2) = tmp;
                clear tmp;
    
                % stimulus location
                data_locked.trialinfo(:, 3) = functional_localizer.data(1, :);
    
                % save
                save(data_locked_file, 'data_locked', '-v7.3');

            end

        end   

        
        %% segment continuous data of WM task into trials
        
        number_blocks = 3;
        
        for block_i = 1:number_blocks

            intervals = {'wholetrial' 'wholetrial_baseline1sec' 'delay' 'baseline_1sec'};
            for inter_i = 1:length(intervals)

                filename = intervals{inter_i};

                data_locked_file = [path_outputs '/' subject_ID '/' session_ID ...
                    '/GERE_task_block' num2str(block_i) '_' filename '.mat'];

                if isfile(data_locked_file)
                
                    load(data_locked_file)

                    if ~isfield(data_locked, 'trialinfo')

                        % reset sampleinfo
                        data_locked.sampleinfo = data_locked.sampleinfo(:,1:2);
                        try
                            data_locked = rmfield(data_locked, 'sampleinfolegend');
                        catch
                        end
    
                        % new field with trial information
                        data_locked.trialinfo = nan(size(data_locked.sampleinfo, 1), 11);
                        data_locked.trialinfolegend = {'trial number'; 'valid trials (no dim task / response)'; 'correct responses';...
                            'stim1'; 'stim2'; 'stim3'; 'stim4';...
                            'resp1'; 'resp2'; 'resp3'; 'resp4'};
            
                        % trial number
                        data_locked.trialinfo(:, 1) = 1:(size(data_locked.trialinfo, 1));
        
                        % valid trials
                        load([path_behavior '/' subject_ID '/' session_ID '/block_0' num2str(block_i) '/task_trials_dimfixation.mat']); % dim
                        
                        % adjust length (for exceptions describe in
                        % script_02_trigger_info_sequential_WM.m
                        tmp = (~dim.dim_task_trial)';
                        tmp = tmp(1:(size(data_locked.trial, 2)));
                        
                        data_locked.trialinfo(:, 2) = tmp';
    
                        clear tmp
        
                        % correct responses
                        load([path_behavior '/' subject_ID '/' session_ID '/block_0' num2str(block_i) '/task_trial_presentation.mat']); % sequences
                        load([path_behavior '/' subject_ID '/' session_ID '/block_0' num2str(block_i) '/responses_task.mat']); % responses
                        correct_length3 = 0; correct_length4 = 0;
                        number_trials_length3 = 0; number_trials_length4 = 0;
                        responses = responses(find(~isnan(responses(:,1))),:);
                        sequences = sequences(find(~isnan(responses(:,1))),:);
        
                        tmp = zeros(size(responses, 1), 1);
                        
                        for trial_i = 1:size(sequences, 1)
                        
                            if sequences(trial_i, end) ~= 0 % length 4
                                        
                                if sequences(trial_i, :) == responses(trial_i, :)
                        
                                    tmp(trial_i) = 1;
                        
                                end
                        
                            else % length 3
                                        
                                if sequences(trial_i, 1:3) == responses(trial_i, 1:3)
                        
                                    tmp(trial_i) = 1;
                        
                                end
                        
                            end
                        
                        end
        
                        % adjust length (for exceptions describe in
                        % script_02_trigger_info_sequential_WM.m
                        tmp = tmp(1:(size(data_locked.trial, 2)),:);
    
                        data_locked.trialinfo(:, 3) = tmp;
    
                        clear tmp
    
                        % stimulus presentation
                        tmp = sequences(1:(size(data_locked.trial, 2)),:);
                        data_locked.trialinfo(:,4:7) = tmp;
                        clear tmp
    
                        % responses
                        tmp = responses(1:(size(data_locked.trial, 2)),:);
                        data_locked.trialinfo(:,8:11) = tmp;
                                        
                        % save
                        save(data_locked_file, 'data_locked', '-v7.3');

                    end
    
                end

            end
        
        end

    end

end
