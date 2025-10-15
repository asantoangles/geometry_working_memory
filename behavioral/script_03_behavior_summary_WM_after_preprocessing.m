%% behavioral analysis of sequential working memory task

path_results = '/path_to_local/data/behavior';
path_preproc = '/path_to_local/results/preprocessing/M3';

zthr = 25;



%% after preprocessing

behavior_task = zeros(length(subjects)*2, 5); % columns: (1) empty column, (2) correct responses length 3, (3) correct responses length 4, (4) total length 3, (5) total length 4

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

        disp(['sub_' subjectID ' - sess_' sessionID]);

        idx = sub_i+((ses_i-1)*length(subjects));
                        
        %% task

        load([path_preproc '/sub_' subjectID '/sess_' sessionID '/GERE_task_allblocks_wholetrial_baseline1sec_validtrialszthr10.mat']); % data_task

        for trial_i = 1:size(data_task.trial, 2)

            % if correct response
            if data_task.trialinfo(trial_i, 3) == 1

                if data_task.trialinfo(trial_i, 7) == 0                 % sequence length 3
                    behavior_task(idx, 2) = behavior_task(idx, 2) + 1;
                else                                                    % sequence length 4
                    behavior_task(idx, 3) = behavior_task(idx, 3) + 1;
                end

            end

            % count valid trials
            if data_task.trialinfo(trial_i, 7) == 0                 % sequence length 3
                behavior_task(idx, 4) = behavior_task(idx, 4) + 1;
            else                                                    % sequence length 4
                behavior_task(idx, 5) = behavior_task(idx, 5) + 1;
            end

        end

    end

end
          
% save
summary_behav_preproc = [];
summary_behav_preproc.sample = subjects;
summary_behav_preproc.info = {'behavior_task: (1) empty, (2) correct responses length 3, (3) correct responses length 4, (4) total length 3, (5) total length 4',...
                              'rows: session 1 of all subjects, then session 2 of all subjects'};
summary_behav_preproc.behavior_task = behavior_task;

if ~isfile([path_preproc '/summary_behav_preproc_n17_zthr' num2str(zthr) '.mat'])

    save([path_preproc '/summary_behav_preproc_n17_zthr' num2str(zthr) '.mat'], 'summary_behav_preproc');

else

    disp('ERROR: existing file, change filename');

end


%% extract information - behavioral after preprocessing

load([path_preproc '/summary_behav_preproc_n17_zthr' num2str(zthr) '.mat']); % summary_behav_preproc
summary_behav_preproc.info'

% WM performance - sequence length 3 (session 1)
mean(summary_behav_preproc.behavior_task(1:length(subjects),2))
std(summary_behav_preproc.behavior_task(1:length(subjects),2))
min(summary_behav_preproc.behavior_task(1:length(subjects),2))
max(summary_behav_preproc.behavior_task(1:length(subjects),2))

% WM performance - sequence length 3 (session 2)
mean(summary_behav_preproc.behavior_task((length(subjects)+1):end,2))
std(summary_behav_preproc.behavior_task((length(subjects)+1):end,2))
min(summary_behav_preproc.behavior_task((length(subjects)+1):end,2))
max(summary_behav_preproc.behavior_task((length(subjects)+1):end,2))

% WM performance - sequence length 4 (session 1)
mean(summary_behav_preproc.behavior_task(1:length(subjects),3))
std(summary_behav_preproc.behavior_task(1:length(subjects),3))
min(summary_behav_preproc.behavior_task(1:length(subjects),3))
max(summary_behav_preproc.behavior_task(1:length(subjects),3))

% WM performance - sequence length 4 (session 2)
mean(summary_behav_preproc.behavior_task((length(subjects)+1):end,3))
std(summary_behav_preproc.behavior_task((length(subjects)+1):end,3))
min(summary_behav_preproc.behavior_task((length(subjects)+1):end,3))
max(summary_behav_preproc.behavior_task((length(subjects)+1):end,3))


% WM performance - sequence all lengths (session 1)
tmp = [summary_behav_preproc.behavior_task(1:length(subjects),2) summary_behav_preproc.behavior_task(1:length(subjects),3)];
tmp = tmp(:,1) + tmp(:, 2);
mean(tmp)
std(tmp)
min(tmp)
max(tmp)

% WM performance - sequence all lengths (session 2)
tmp = [summary_behav_preproc.behavior_task((length(subjects)+1):end,2) summary_behav_preproc.behavior_task((length(subjects)+1):end,3)];
tmp = tmp(:,1) + tmp(:, 2);
mean(tmp)
std(tmp)
min(tmp)
max(tmp)


% WM performance - sequence length 3 (all sessions)
tmp = [summary_behav_preproc.behavior_task(:,2)];
mean(tmp)
std(tmp)
min(tmp)
max(tmp)

% WM performance - sequence length 4 (all sessions)
tmp = [summary_behav_preproc.behavior_task(:,3)];
mean(tmp)
std(tmp)
min(tmp)
max(tmp)

%%% pool sessions

correct_length3 = [summary_behav_preproc.behavior_task(1:17,2) summary_behav_preproc.behavior_task(18:34,2)];
correct_length3 = correct_length3(:,1) + correct_length3(:, 2);
correct_length4 = [summary_behav_preproc.behavior_task(1:17,3) summary_behav_preproc.behavior_task(18:34,3)];
correct_length4 = correct_length4(:,1) + correct_length4(:, 2);
correct_lengthall = correct_length3 + correct_length4;

total_length3 = [summary_behav_preproc.behavior_task(1:17,4) summary_behav_preproc.behavior_task(18:34,4)];
total_length3 = total_length3(:,1) + total_length3(:, 2);
total_length4 = [summary_behav_preproc.behavior_task(1:17,5) summary_behav_preproc.behavior_task(18:34,5)];
total_length4 = total_length4(:,1) + total_length4(:, 2);
total_lengthall = total_length3 + total_length4;

% correct trials - sequence all lengths (all sessions)
mean(correct_lengthall)
std(correct_lengthall)
min(correct_lengthall)
max(correct_lengthall)

% all trials - sequence all lengths (all sessions)
mean(total_lengthall)
std(total_lengthall)
min(total_lengthall)
max(total_lengthall)


%% linear regression

responses_length3_session1 = summary_behav_preproc.behavior_task(1:17,2)/115;
responses_length3_session2 = summary_behav_preproc.behavior_task(18:34,2)/110;
responses_length4_session1 = summary_behav_preproc.behavior_task(1:17,3)/95;
responses_length4_session2 = summary_behav_preproc.behavior_task(18:34,3)/100;

responses_var = [responses_length3_session1; responses_length3_session2; responses_length4_session1; responses_length4_session2];
length_var = [3*ones(1,length(subjects)*2) 4*ones(1,length(subjects)*2)]';
session_var = [ones(1,length(subjects)) 2*ones(1,length(subjects)) ones(1,length(subjects)) 2*ones(1,length(subjects))]';

% Combine into a table
dataTable = table(length_var, session_var, responses_var);

% Fit a linear model with interaction terms
lm = fitlm(dataTable, 'responses_var ~ length_var * session_var');




%% binomial test - above chance responses

number_locations = 8;
p_chance_3 = 1/(number_locations^3);
p_chance_4 = 1/(number_locations^4);

for sub_i = 1:17

    disp(' ')
    disp(['subject ' num2str(sub_i)])

    % Inputs
    N_3 = summary_behav_preproc.behavior_task(sub_i,4);  % Total number of 3-item sequence trials
    N_4 = summary_behav_preproc.behavior_task(sub_i,5);  % Total number of 4-item sequence trials
    k_3 = summary_behav_preproc.behavior_task(sub_i,2);  % Number of correct trials for 3-item sequences
    k_4 = summary_behav_preproc.behavior_task(sub_i,3);  % Number of correct trials for 4-item sequences
        
    % Binomial tests
    % For 3-item sequences
    p_value_3 = 1 - binocdf(k_3 - 1, N_3, p_chance_3);
    
    % For 4-item sequences
    p_value_4 = 1 - binocdf(k_4 - 1, N_4, p_chance_4);
    
    % Display results
    disp([num2str(p_value_3) ' ' num2str(p_value_4)]);

end

