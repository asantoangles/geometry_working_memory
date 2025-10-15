%% preprocessing summary

% path to data
path_results = ['/path_to_local/results/preprocessing/' folder];

zthr = 25;

rejected_trials = nan(length(subjects)*2, 2); % columns: localizer, task
rejected_components = nan(length(subjects)*2, 2); % columns: localizer, task

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

        idx = sub_i+((ses_i-1)*length(subjects));
        
        %% rejected trials 

        % localizer
        load([path_results '/' subject_ID '/' session_ID ...
            '/GERE_localizer_validtrialszthr' num2str(zthr) '.mat']); % data_validtrials

        total_trials = 560;
        rejected_trials(idx, 1) = total_trials - size(data_validtrials.trial, 2);

        % task
        load([path_results '/' subject_ID '/' session_ID ...
            '/GERE_task_allblocks_wholetrial_baseline1sec_validtrialszthr' num2str(zthr) '.mat']); % data_task

        total_trials = 210;
        rejected_trials(idx, 2) = total_trials - size(data_task.trial, 2);

        %% rejected components

        % localizer
        load([path_results '/' subject_ID '/' session_ID ...
            '/GERE_localizer_ICA.mat']); % data_comp

        rejected_components(idx, 1) = length(data_comp.rejected_components);

        % check rank of data
        load([path_results '/' subject_ID '/' session_ID ...
            '/GERE_localizer_validtrialszthr' num2str(zthr) '_postICA.mat']); % data
        data_rank = rank(data.trial{1}*data.trial{1}');

        if (207 - data_rank) ~= length(data_comp.rejected_components)
            disp('Attention here! missmatch between rejected components and data rank');
        end

        % task
        load([path_results '/' subject_ID '/' session_ID ...
            '/GERE_task_ICA_wholetrial_baseline1sec.mat']); % data_comp

        rejected_components(idx, 2) = length(data_comp.rejected_components);

        % check rank of data
        load([path_results '/' subject_ID '/' session_ID ...
            '/GERE_task_zthr' num2str(zthr) '_wholetrial_baseline1sec_postICA.mat']); % data
        data_rank = rank(data.trial{1}*data.trial{1}');

        if (207 - data_rank) ~= length(data_comp.rejected_components)
            disp('Attention here! missmatch between rejected components and data rank');
        end


    end

end

% save

summary_preprocessing = [];
summary_preprocessing.sample = subjects;
summary_preprocessing.info = {'rejected_trials: columns (1) localizer, (2) task',...
                              'rejected_components: columns (1) localizer, (2) task',...
                              'rows: session 1 of all subjects, then session 2 of all subjects'};
summary_preprocessing.rejected_trials = rejected_trials;
summary_preprocessing.rejected_components = rejected_components;

if ~isfile([path_results '/summary_preprocessing_n17_zthr' num2str(zthr) '.mat'])

    save([path_results '/summary_preprocessing_n17_zthr' num2str(zthr) '.mat'], 'summary_preprocessing');

else

    disp('ERROR: existing file, change filename');

end

%% extract information - localizer

load([path_results '/summary_preprocessing_n17_zthr' num2str(zthr) '.mat']); % summary_preprocessing
summary_preprocessing.info'

% rejected trials - localizer
mean(summary_preprocessing.rejected_trials(:,1))
std(summary_preprocessing.rejected_trials(:,1))
min(summary_preprocessing.rejected_trials(:,1))
max(summary_preprocessing.rejected_trials(:,1))

% rejected trials - localizer (session 1)
mean(summary_preprocessing.rejected_trials(1:length(subjects),1))
std(summary_preprocessing.rejected_trials(1:length(subjects),1))
min(summary_preprocessing.rejected_trials(1:length(subjects),1))
max(summary_preprocessing.rejected_trials(1:length(subjects),1))

% rejected trials - localizer (session 2)
mean(summary_preprocessing.rejected_trials((length(subjects)+1):end,1))
std(summary_preprocessing.rejected_trials((length(subjects)+1):end,1))
min(summary_preprocessing.rejected_trials((length(subjects)+1):end,1))
max(summary_preprocessing.rejected_trials((length(subjects)+1):end,1))

% rejected ICs - localizer
mean(summary_preprocessing.rejected_components(:,1))
std(summary_preprocessing.rejected_components(:,1))
min(summary_preprocessing.rejected_components(:,1))
max(summary_preprocessing.rejected_components(:,1))

% rejected ICs - localizer (session 1)
mean(summary_preprocessing.rejected_components(1:length(subjects),1))
std(summary_preprocessing.rejected_components(1:length(subjects),1))
min(summary_preprocessing.rejected_components(1:length(subjects),1))
max(summary_preprocessing.rejected_components(1:length(subjects),1))

% rejected ICs - localizer (session 2)
mean(summary_preprocessing.rejected_components((length(subjects)+1):end,1))
std(summary_preprocessing.rejected_components((length(subjects)+1):end,1))
min(summary_preprocessing.rejected_components((length(subjects)+1):end,1))
max(summary_preprocessing.rejected_components((length(subjects)+1):end,1))


% total number of trials
trials_preproc = abs(summary_preprocessing.rejected_trials(:,1) - 560);
trials_preproc = trials_preproc(1:length(subjects)) + trials_preproc((length(subjects)+1):end);
mean(trials_preproc)
std(trials_preproc)
min(trials_preproc)
max(trials_preproc)

% trials - session 1
trials_preproc = abs(summary_preprocessing.rejected_trials(:,1) - 560);
trials_preproc = trials_preproc(1:length(subjects));
mean(trials_preproc)
std(trials_preproc)
min(trials_preproc)
max(trials_preproc)

% trials - session 2
trials_preproc = abs(summary_preprocessing.rejected_trials(:,1) - 560);
trials_preproc = trials_preproc((length(subjects)+1):end);
mean(trials_preproc)
std(trials_preproc)
min(trials_preproc)
max(trials_preproc)






%%% summary across sessions

% valid trials - localizer
rejected_trials = summary_preprocessing.rejected_trials(1:17,1) + summary_preprocessing.rejected_trials(18:34,1);
valid_trials = abs(rejected_trials - 560*2);
mean(valid_trials)
std(valid_trials)
min(valid_trials)
max(valid_trials)

% valid trials - task
rejected_trials = summary_preprocessing.rejected_trials(1:17,2) + summary_preprocessing.rejected_trials(18:34,2);
valid_trials = abs(rejected_trials - 210*2);
mean(valid_trials)
std(valid_trials)
min(valid_trials)
max(valid_trials)


%% extract information - task

load([path_results '/summary_preprocessing_n17_zthr' num2str(zthr) '.mat']); % summary_preprocessing
summary_preprocessing.info'

% rejected trials - localizer
mean(summary_preprocessing.rejected_trials(:,2))
std(summary_preprocessing.rejected_trials(:,2))
min(summary_preprocessing.rejected_trials(:,2))
max(summary_preprocessing.rejected_trials(:,2))

% rejected trials - localizer (session 1)
mean(summary_preprocessing.rejected_trials(1:length(subjects),2))
std(summary_preprocessing.rejected_trials(1:length(subjects),2))
min(summary_preprocessing.rejected_trials(1:length(subjects),2))
max(summary_preprocessing.rejected_trials(1:length(subjects),2))

% rejected trials - localizer (session 2)
mean(summary_preprocessing.rejected_trials((length(subjects)+1):end,2))
std(summary_preprocessing.rejected_trials((length(subjects)+1):end,2))
min(summary_preprocessing.rejected_trials((length(subjects)+1):end,2))
max(summary_preprocessing.rejected_trials((length(subjects)+1):end,2))

% rejected ICs - localizer
mean(summary_preprocessing.rejected_components(:,2))
std(summary_preprocessing.rejected_components(:,2))
min(summary_preprocessing.rejected_components(:,2))
max(summary_preprocessing.rejected_components(:,2))

% rejected ICs - localizer (session 1)
mean(summary_preprocessing.rejected_components(1:length(subjects),2))
std(summary_preprocessing.rejected_components(1:length(subjects),2))
min(summary_preprocessing.rejected_components(1:length(subjects),2))
max(summary_preprocessing.rejected_components(1:length(subjects),2))

% rejected ICs - localizer (session 2)
mean(summary_preprocessing.rejected_components((length(subjects)+1):end,2))
std(summary_preprocessing.rejected_components((length(subjects)+1):end,2))
min(summary_preprocessing.rejected_components((length(subjects)+1):end,2))
max(summary_preprocessing.rejected_components((length(subjects)+1):end,2))


% total number of trials
trials_preproc = abs(summary_preprocessing.rejected_trials(:,2) - 210);
trials_preproc = trials_preproc(1:length(subjects)) + trials_preproc((length(subjects)+1):end);
mean(trials_preproc)
std(trials_preproc)
min(trials_preproc)
max(trials_preproc)

% trials - session 1
trials_preproc = abs(summary_preprocessing.rejected_trials(:,2) - 210);
trials_preproc = trials_preproc(1:length(subjects));
mean(trials_preproc)
std(trials_preproc)
min(trials_preproc)
max(trials_preproc)

% trials - session 2
trials_preproc = abs(summary_preprocessing.rejected_trials(:,2) - 210);
trials_preproc = trials_preproc((length(subjects)+1):end);
mean(trials_preproc)
std(trials_preproc)
min(trials_preproc)
max(trials_preproc)







