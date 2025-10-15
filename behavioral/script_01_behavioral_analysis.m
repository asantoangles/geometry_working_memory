%% behavioral analysis of sequential working memory task

% extract behavioral responses from psychotoolbox outputs

path = '/path_to_local/data/behavior';

for sub_i = 1:length(subjects)

    subject = subjects(sub_i);

    for ses_i = 1:length(sessions)

        session = sessions(ses_i);

        if subject < 10
            subjectID = ['0' num2str(subject)];
        else
            subjectID = ['' num2str(subject)];
        end
        
        sessionID = ['0' num2str(session)];
                
        %% functional localizer - counting for the delayed responses
        
        % load data
        load([path '/sub_' subjectID '/sess_' sessionID '/functional_localizer.mat']); % functional_localizer
        load([path '/sub_' subjectID '/sess_' sessionID '/responses_localizer.mat']); % responses
        
        unique_responses = responses;
        
        for stim_i = 1:length(unique_responses)
            if unique_responses(stim_i) == 1
                unique_responses(stim_i + 1) = 0;
            end
        end
        
        total_responses = sum(unique_responses);
        
        % count as response the previous stimulus at which the response was
        % recorded (because of delayed response)
        for stim_i = 2:length(responses)
            if responses(stim_i) == 1
                responses(stim_i - 1) = 1;
            end
        end
        
        disp(' ');
        disp('Localizer - dim flickering fixation point');
        
        % color task
        total_hits = functional_localizer.data(2,:) + responses;
        hit = length(total_hits(total_hits == 2));
        total_dim = sum(functional_localizer.data(2,:));
        disp(['Color: hits = ' num2str(hit) '/' num2str(total_dim)]);
        
        % dim task
        total_hits = functional_localizer.data(3,:) + responses;
        hit = length(total_hits(total_hits == 2));
        total_dim = sum(functional_localizer.data(2,:));
        disp(['Flickering: hits = ' num2str(hit) '/' num2str(total_dim)]);
        
        % false positives
        total_hits = functional_localizer.data(2,:) + functional_localizer.data(3,:) + responses;
        false_positives = sum(total_hits(total_hits == 2))/2;
        disp(['Total responses = ' num2str(total_responses)]);
        
        disp(' ');
        
        
        %% main task
        
        blocks = 3;
        
        %% subject-specific corrections
        
        if subject == 20
            if session == 1
                blocks = 2;
            end
        end

        if subject == 51
            if session == 2
                blocks = 2;
            end
        end

        if subject == 49
            if session == 1
                blocks = 2;
            end
        end
        
        for block_i = 1:blocks
        
            %% behavioral responses to main task
        
            % load data
            load([path '/sub_' subjectID '/sess_' sessionID '/block_0' num2str(block_i) '/task_trial_presentation.mat']); % sequences
            load([path '/sub_' subjectID '/sess_' sessionID '/block_0' num2str(block_i) '/responses_task.mat']); % responses
        
            % number of correct sequences, depending on sequence length
            correct_length3 = 0;
            correct_length4 = 0;
            
            number_trials_length3 = 0;
            number_trials_length4 = 0;
            
            % remove empty responses
            responses = responses(find(~isnan(responses(:,1))),:);
            sequences = sequences(find(~isnan(responses(:,1))),:);
            
            for trial_i = 1:size(sequences, 1)
            
                if sequences(trial_i, end) ~= 0 % length 4
            
                    number_trials_length4 = number_trials_length4 + 1;
            
                    if sequences(trial_i, :) == responses(trial_i, :)
            
                        correct_length4 = correct_length4 + 1;
            
                    end
            
                else % length 3
            
                    number_trials_length3 = number_trials_length3 + 1;
            
                    if sequences(trial_i, 1:3) == responses(trial_i, 1:3)
            
                        correct_length3 = correct_length3 + 1;
            
                    end
            
                end
            
            end
            
            disp(' ');
            disp(['WM task - block ' num2str(block_i)]);
            disp(['Correct responses (length 3): ' num2str(correct_length3) ' out of ' num2str(number_trials_length3) ' [' num2str(round(correct_length3/number_trials_length3, 2)*100) '%]']);
            disp(['Correct responses (length 4): ' num2str(correct_length4) ' out of ' num2str(number_trials_length4) ' [' num2str(round(correct_length4/number_trials_length4, 2)*100) '%]']);
        
            %% dim flickering task fixation point
            
            % load data
            load([path '/sub_' subjectID '/sess_' sessionID '/block_0' num2str(block_i) '/task_trial_presentation.mat']); % sequences
            load([path '/sub_' subjectID '/sess_' sessionID '/block_0' num2str(block_i) '/task_trials_dimfixation.mat']); % dim
        
            total_hits = dim.dim_task_trial + dim.response;
            hit = length(total_hits(total_hits == 2));
            total_dim = sum(dim.dim_task_trial);
            errors = length(total_hits(total_hits == 1));
            false_negative = total_dim - hit;
            false_positive = errors - false_negative;
        
            disp(['Flickering: hits = ' num2str(hit) '/' num2str(total_dim)]);
            disp(['Flickering: false negatives = ' num2str(false_negative) '/' num2str(total_dim)]);
            disp(['Flickering: false positives = ' num2str(false_positive)]);
        
        
        end

    end

end
          
