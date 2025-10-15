%% count outputs

% path to outputs
path_results = ['/path_to_local/results/preprocessing/' folder];

%% summary valid trials

zthr = 25;

summary_table = zeros(length(subjects), 8); 
rownames = {};
colnames = {'localizer sess1' 'localizer sess2' 'task sess1' 'correct sess1' 'task sess2' 'correct sess2' 'task total' 'correct total'};

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

        if session == 1

            if isempty(rownames) 
    
                rownames{1} = subject_ID;
    
            else
    
                rownames{end+1} = subject_ID;
                
            end

        end


        %% localizer

        if isfile([path_results '/' subject_ID '/' session_ID ...
                '/GERE_localizer_validtrialszthr' num2str(zthr) '.mat'])

            load([path_results '/' subject_ID '/' session_ID ...
                '/GERE_localizer_validtrialszthr' num2str(zthr) '.mat']); % data_validtrials

            if session == 1
                summary_table(sub_i, 1) = size(data_validtrials.time, 2);
            else
                summary_table(sub_i, 2) = size(data_validtrials.time, 2);
            end


        end

        %% task

        for block_i = 1:3

            if isfile([path_results '/' subject_ID '/' session_ID ...
                        '/GERE_task_block' num2str(block_i) '_wholetrial_baseline1sec_validtrialszthr'...
                        num2str(zthr) '.mat'])

                load([path_results '/' subject_ID '/' session_ID ...
                    '/GERE_task_block' num2str(block_i) '_wholetrial_baseline1sec_validtrialszthr'...
                    num2str(zthr) '.mat']);

                if session == 1
                    summary_table(sub_i, 3) = summary_table(sub_i, 3) + size(data_validtrials.time, 2);
                    summary_table(sub_i, 4) = summary_table(sub_i, 4) + sum(data_validtrials.trialinfo(:,3));
                else
                    summary_table(sub_i, 5) = summary_table(sub_i, 5) + size(data_validtrials.time, 2);
                    summary_table(sub_i, 6) = summary_table(sub_i, 6) + sum(data_validtrials.trialinfo(:,3));
                end

            end

        end

    end

end

% valid trials in two session
summary_table(:,7) = summary_table(:, 3) + summary_table(:, 5);
summary_table(:,8) = summary_table(:, 4) + summary_table(:, 6);

% mean and std
summary_table(end+1,:) = mean(summary_table);
summary_table(end+1,:) = std(summary_table);

rownames{end+1} = 'mean';
rownames{end+1} = 'std';

% add row/column names
summary_table = array2table(summary_table, 'RowNames', rownames, 'VariableNames', colnames);

% save
save([path_results '/summary_valid_trials_zthr' num2str(zthr) '.mat'], 'summary_table');

