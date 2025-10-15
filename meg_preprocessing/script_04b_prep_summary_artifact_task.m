%% automatic artifact detection (trial-wise)
% https://www.fieldtriptoolbox.org/tutorial/automatic_artifact_rejection/

% path to outputs
path_results = ['/path_to_local/results/preprocessing/' folder];

%% inspect results - survival trials

threshold_muscle_blinks = [25];

% emtpy matrix for results
rownames        = {'wholetrial',...
                   'wholetrial_baseline1sec',...
                   'delay',...    
                   'baseline'};
colnames = {'25'};

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

        if isfile([path_results '/' subject_ID '/' session_ID '/GERE_task_artifact_thr_20to30.mat'])

            load([path_results '/' subject_ID '/' session_ID '/GERE_task_artifact_thr_20to30.mat']);
                
            survivor_trials_table = table2array(survivor_trials_table);

            summary_trials = nan(4,3);
            for col_i = 1:3
                
               summary_trials(1,col_i) = sum(survivor_trials_table(1:3,col_i), 'omitnan');
               summary_trials(2,col_i) = sum(survivor_trials_table(4:6,col_i), 'omitnan');
               summary_trials(3,col_i) = sum(survivor_trials_table(7:9,col_i), 'omitnan');
               summary_trials(4,col_i) = sum(survivor_trials_table(10:12,col_i), 'omitnan');

            end

            summary_trials = array2table(summary_trials, 'RowNames', rownames, 'VariableNames', colnames)

        end

    end

end
