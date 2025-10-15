%% behavioral analysis of sequential working memory task

% summary of behavioral

path_results = '/path_to_local/data/behavior';
path_preproc = '/path_to_local/results/preprocessing/M3';

%% before preprocessing

behavior_localizer = nan(length(subjects)*2, 2); % columns: (1) flickering responses, (2) color responses
behavior_task = zeros(length(subjects)*2, 5); % columns: (1) flickering responses, (2) correct responses length 3, (3) correct responses length 4
                                              %          (4) total trials
                                              %          length 3
                                              %          (5) total trials
                                              %          length 4
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

        idx = sub_i+((ses_i-1)*length(subjects));
                
        %% functional localizer - counting for the delayed responses
        
        % load data
        load([path_results '/sub_' subjectID '/sess_' sessionID '/functional_localizer.mat']); % functional_localizer
        load([path_results '/sub_' subjectID '/sess_' sessionID '/responses_localizer.mat']); % responses
        
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
                
        % color task
        total_hits = functional_localizer.data(2,:) + responses;
        hit = length(total_hits(total_hits == 2));
        total_dim = sum(functional_localizer.data(2,:));

        behavior_localizer(idx, 2) = hit;
        
        % dim task
        total_hits = functional_localizer.data(3,:) + responses;
        hit = length(total_hits(total_hits == 2));
        total_dim = sum(functional_localizer.data(2,:));

        behavior_localizer(idx, 1) = hit;
        
        %% main task
        
        blocks = 3;
                        
        for block_i = 1:blocks
        
            %% behavioral responses to main task
        
            % load data
            load([path_results '/sub_' subjectID '/sess_' sessionID '/block_0' num2str(block_i) '/task_trial_presentation.mat']); % sequences
            load([path_results '/sub_' subjectID '/sess_' sessionID '/block_0' num2str(block_i) '/responses_task.mat']); % responses
        
            % number of correct sequences, depending on sequence length
            correct_length3 = 0;
            correct_length4 = 0;
            
            number_trials_length3 = 0;
            number_trials_length4 = 0;
            
            % remove empty responses
            responses = responses(find(~isnan(responses(:,1))),:);
            sequences = sequences(find(~isnan(responses(:,1))),:);
            
            for trial_i = 1:size(sequences, 1)
            
                if sequences(trial_i, 4) ~= 0 % length 4
            
                    number_trials_length4 = number_trials_length4 + 1;
            
                    if sum(sequences(trial_i, :) == responses(trial_i, :)) == 4
            
                        correct_length4 = correct_length4 + 1;
            
                    end
            
                else % length 3
            
                    number_trials_length3 = number_trials_length3 + 1;
            
                    if sum(sequences(trial_i, 1:3) == responses(trial_i, 1:3)) == 3
            
                        correct_length3 = correct_length3 + 1;
            
                    end
            
                end
            
            end

            behavior_task(idx, 2) = behavior_task(idx, 2) + correct_length3;
            behavior_task(idx, 3) = behavior_task(idx, 3) + correct_length4;

            behavior_task(idx, 4) = behavior_task(idx, 4) + number_trials_length3;
            behavior_task(idx, 5) = behavior_task(idx, 5) + number_trials_length4;
        
            %% dim flickering task fixation point
            
            % load data
            load([path_results '/sub_' subjectID '/sess_' sessionID '/block_0' num2str(block_i) '/task_trial_presentation.mat']); % sequences
            load([path_results '/sub_' subjectID '/sess_' sessionID '/block_0' num2str(block_i) '/task_trials_dimfixation.mat']); % dim
        
            total_hits = dim.dim_task_trial + dim.response;
            hit = length(total_hits(total_hits == 2));
            total_dim = sum(dim.dim_task_trial);
            errors = length(total_hits(total_hits == 1));
            false_negative = total_dim - hit;
            false_positive = errors - false_negative;

            behavior_task(idx, 1) = behavior_task(idx, 1) + hit;
        
        end

    end

end
          
% save

summary_behavioral = [];
summary_behavioral.sample = subjects;
summary_behavioral.info = {'behavior_localizer: columns: (1) flickering responses, (2) color responses',...
                              'behavior_task: (1) flickering responses, (2) correct responses length 3, (3) correct responses length 4, (4) total trials length 3, (5) total trials length 4',...
                              'rows: session 1 of all subjects, then session 2 of all subjects'};
summary_behavioral.behavior_localizer = behavior_localizer;
summary_behavioral.behavior_task = behavior_task;

if ~isfile([path_preproc '/summary_behavioral_n17.mat'])

    save([path_preproc '/summary_behavioral_n17.mat'], 'summary_behavioral');

else

    disp('ERROR: existing file, change filename');

end



%% extract information - behavioral before preprocessing

load([path_preproc '/summary_behavioral_n17.mat']); % summary_behavioral
summary_behavioral.info'

% localizer fixation task (total 12)
tmp = summary_behavioral.behavior_localizer(:,1);
tmp(tmp == 0) = 12;
summary_behavioral.behavior_localizer(:,1) = tmp;

mean(summary_behavioral.behavior_localizer(:,1))
std(summary_behavioral.behavior_localizer(:,1))
min(summary_behavioral.behavior_localizer(:,1))
max(summary_behavioral.behavior_localizer(:,1))

% localizer color task (total 12)
tmp = summary_behavioral.behavior_localizer(:,2);
tmp(tmp == 0) = 12;
summary_behavioral.behavior_localizer(:,2) = tmp;

mean(summary_behavioral.behavior_localizer(:,2))
std(summary_behavioral.behavior_localizer(:,2))
min(summary_behavioral.behavior_localizer(:,2))
max(summary_behavioral.behavior_localizer(:,2))

% WM task fixation (total 10)
mean(summary_behavioral.behavior_task(:,1))
std(summary_behavioral.behavior_task(:,1))
min(summary_behavioral.behavior_task(:,1))
max(summary_behavioral.behavior_task(:,1))

% WM performance - sequence length 3 (session 1)
mean(summary_behavioral.behavior_task(1:length(subjects),2))
std(summary_behavioral.behavior_task(1:length(subjects),2))
min(summary_behavioral.behavior_task(1:length(subjects),2))
max(summary_behavioral.behavior_task(1:length(subjects),2))

% WM performance - sequence length 3 (session 2)
mean(summary_behavioral.behavior_task((length(subjects)+1):end,2))
std(summary_behavioral.behavior_task((length(subjects)+1):end,2))
min(summary_behavioral.behavior_task((length(subjects)+1):end,2))
max(summary_behavioral.behavior_task((length(subjects)+1):end,2))

% WM performance - sequence length 4 (session 1)
mean(summary_behavioral.behavior_task(1:length(subjects),3))
std(summary_behavioral.behavior_task(1:length(subjects),3))
min(summary_behavioral.behavior_task(1:length(subjects),3))
max(summary_behavioral.behavior_task(1:length(subjects),3))

% WM performance - sequence length 4 (session 2)
mean(summary_behavioral.behavior_task((length(subjects)+1):end,3))
std(summary_behavioral.behavior_task((length(subjects)+1):end,3))
min(summary_behavioral.behavior_task((length(subjects)+1):end,3))
max(summary_behavioral.behavior_task((length(subjects)+1):end,3))




% WM performance - sequence all lengths (session 1)
tmp = [summary_behavioral.behavior_task(1:length(subjects),2) summary_behavioral.behavior_task(1:length(subjects),3)];
tmp = tmp(:,1) + tmp(:, 2);
mean(tmp)
std(tmp)
min(tmp)
max(tmp)

% WM performance - sequence all lengths (session 2)
tmp = [summary_behavioral.behavior_task((length(subjects)+1):end,2) summary_behavioral.behavior_task((length(subjects)+1):end,3)];
tmp = tmp(:,1) + tmp(:, 2);
mean(tmp)
std(tmp)
min(tmp)
max(tmp)

% WM performance - sequence length 3 (all sessions)
tmp = [summary_behavioral.behavior_task(:,2)];
mean(tmp)
std(tmp)
min(tmp)
max(tmp)

% WM performance - sequence length 4 (all sessions)
tmp = [summary_behavioral.behavior_task(:,3)];
mean(tmp)
std(tmp)
min(tmp)
max(tmp)



%% linear regression

responses_length3_session1 = summary_behavioral.behavior_task(1:17,2)/115;
responses_length3_session2 = summary_behavioral.behavior_task(18:34,2)/110;
responses_length4_session1 = summary_behavioral.behavior_task(1:17,3)/95;
responses_length4_session2 = summary_behavioral.behavior_task(18:34,3)/100;

responses_var = [responses_length3_session1; responses_length3_session2; responses_length4_session1; responses_length4_session2];
length_var = [3*ones(1,length(subjects)*2) 4*ones(1,length(subjects)*2)]';
session_var = [ones(1,length(subjects)) 2*ones(1,length(subjects)) ones(1,length(subjects)) 2*ones(1,length(subjects))]';

% Combine into a table
dataTable = table(length_var, session_var, responses_var);

% Fit a linear model with interaction terms
lm = fitlm(dataTable, 'responses_var ~ length_var * session_var');
