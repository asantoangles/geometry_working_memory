%% automatic artifact detection (trial-wise)
% https://www.fieldtriptoolbox.org/tutorial/automatic_artifact_rejection/

% path to outputs
path_results = ['/path_to_local/results/preprocessing/' folder];

%% inspect results for all subjects

threshold_muscle_blinks = [25];

output = nan(1, length(threshold_muscle_blinks));
rownames = {};
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

        if isfile([path_results '/' subject_ID '/' session_ID '/GERE_localizer_artifact_thr_20to30.mat'])

            load([path_results '/' subject_ID '/' session_ID '/GERE_localizer_artifact_thr_20to30.mat']);
            
            if sub_i == 1 && ses_i == 1

                output(1,:) = table2array(survivor_trials_table);
                rownames{1} = ['sub' num2str(subject) ' - sess' num2str(session)];

            else

                output(end+1,:) = table2array(survivor_trials_table);
                rownames{end+1} = ['sub' num2str(subject) ' - sess' num2str(session)];
                
            end

        end

    end

end

% add row/column names
output_allsubjects = array2table(output, 'RowNames', rownames, 'VariableNames', colnames)

