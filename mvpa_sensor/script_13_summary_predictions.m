%% GERE project

% path to data
if isfolder('/path_to_local')
    path_preprocessed_data = '/path_to_local/results/preprocessing/M3';
    path_predictions = '/path_to_local/github/results/replay/data/predictions';
    path_results = '/path_to_local/results/mvpa/P1';
end

number_stimuli = 8;
classifiers = {'l2_e-1'}; classifier = classifiers{1};
inputs = {'bins'};
time_delay = 791; % bins
subjects = 1:17;
sessions = 1:2;

folders = {'task_stim2delay', 'localizer2delay',...
           'task_stim2task_stim', 'localizer2task_stim'};

%% summary of decoded state space - prob

if 1 == 2

    for folder_i = 1:length(folders)
    
        folder_space = folders{folder_i};
    
        fprintf('\n==================\n');
        disp(['prob ' folder_space]);
        fprintf('==================\n');
    
        if ~isfile([path_predictions '/' folder_space '/prob_in_out.mat'])
    
            % preallocate cells
            prob_l3_in_all  = cell(length(subjects), length(sessions));
            prob_l3_out_all = cell(length(subjects), length(sessions));
            prob_l4_in_all  = cell(length(subjects), length(sessions));
            prob_l4_out_all = cell(length(subjects), length(sessions));
            
            for sub_i = 1:length(subjects)
            
                subject = subjects(sub_i);
            
                % loop over sessions
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
                
                    % load data
                    load([path_preprocessed_data '/' subject_ID '/' session_ID...
                            '/GERE_task_zthr25_wholetrial_baseline1sec_postICA.mat']); % data
    
    
                    % load predictions
                    load([path_predictions '/' subject_ID '/' session_ID '/' folder_space '/' inputs{1} '/' classifiers{1} '/predictions.mat']);
            
                    % dval to prob
                    predictions.prob = 1 ./ (1 + exp(-predictions.dval));
                    
                    % --- use number of trials from predictions ---
                    nTrials = size(predictions.prob,1);
                    nTime   = size(predictions.prob,2);
                    nStim   = size(predictions.prob,3);
                    
                    %% -------------------------
                    % Extract sequences (aligned to predictions)
                    % -------------------------
                    seq_all = cell(nTrials,1);
                    seq_len = zeros(nTrials,1);
                    
                    for trial_i = 1:nTrials
                        seq = data.trialinfo(trial_i,4:7);
                        seq(seq == 0) = [];
                        seq_all{trial_i} = seq;
                        seq_len(trial_i) = length(seq);
                    end
                    
                    % indices for each condition
                    idx_l3 = find(seq_len == 3);
                    idx_l4 = find(seq_len == 4);
                    
                    n_l3 = length(idx_l3);
                    n_l4 = length(idx_l4);
                    
                    %% -------------------------
                    % Preallocate tensors
                    % -------------------------
                    prob_l3_in  = zeros(n_l3, nTime, 3);
                    prob_l3_out = zeros(n_l3, nTime, nStim-3);
                    
                    prob_l4_in  = zeros(n_l4, nTime, 4);
                    prob_l4_out = zeros(n_l4, nTime, nStim-4);
                    
                    %% -------------------------
                    % Fill tensors
                    % -------------------------
                    for k = 1:n_l3
                    
                        trial_i = idx_l3(k);
                        sequence = seq_all{trial_i};
                    
                        % --- IN (keep repeats, preserve order) ---
                        stim_in = sequence;
                    
                        % --- OUT (from missing unique stimuli, trimmed) ---
                        stim_unique  = unique(sequence);
                        stim_missing = setdiff(1:nStim, stim_unique);
                    
                        n_out = nStim - length(stim_in);
                        stim_out = stim_missing(1:n_out);
                    
                        % assign
                        prob_l3_in(k,:,:)  = predictions.prob(trial_i,:,stim_in);
                        prob_l3_out(k,:,:) = predictions.prob(trial_i,:,stim_out);
                    
                    end
                    
                    for k = 1:n_l4
                    
                        trial_i = idx_l4(k);
                        sequence = seq_all{trial_i};
                    
                        stim_in = sequence;
                    
                        stim_unique  = unique(sequence);
                        stim_missing = setdiff(1:nStim, stim_unique);
                    
                        n_out = nStim - length(stim_in);
                        stim_out = stim_missing(1:n_out);
                    
                        prob_l4_in(k,:,:)  = predictions.prob(trial_i,:,stim_in);
                        prob_l4_out(k,:,:) = predictions.prob(trial_i,:,stim_out);
                    
                    end
            
                    % -------------------------
                    % STORE (subject × session)
                    % -------------------------
                    prob_l3_in_all{sub_i, ses_i}  = prob_l3_in;
                    prob_l3_out_all{sub_i, ses_i} = prob_l3_out;
            
                    prob_l4_in_all{sub_i, ses_i}  = prob_l4_in;
                    prob_l4_out_all{sub_i, ses_i} = prob_l4_out;        
                            
                end
            
            end
            
            prob_in_out = [];
            prob_in_out.prob_l3_in_all = prob_l3_in_all;
            prob_in_out.prob_l3_out_all = prob_l3_out_all;
            prob_in_out.prob_l4_in_all = prob_l4_in_all;
            prob_in_out.prob_l4_out_all = prob_l4_out_all;
            
            if ~isfolder([path_predictions '/' folder_space])
                mkdir([path_predictions '/' folder_space])
            end
            
            save([path_predictions '/' folder_space '/prob_in_out.mat'], 'prob_in_out');
        
        else
    
            load([path_predictions '/' folder_space '/prob_in_out.mat']); % prob_in_out
            prob_l3_in_all = prob_in_out.prob_l3_in_all;
            prob_l3_out_all = prob_in_out.prob_l3_out_all;
            prob_l4_in_all = prob_in_out.prob_l4_in_all;
            prob_l4_out_all = prob_in_out.prob_l4_out_all;
            
        end
    
        %% model 01
        
        % all probabilities
        
        %% =========================================================
        % HIERARCHICAL MODEL: P(in) > P(out)
        % Level 1: within-subject (pooled L3 + L4, and pooled sessions)
        % Level 2: across-subject t-test on subject t-values
        %% =========================================================
        
        fprintf('\n model 01 - all probabilities\n');
        
        nSub  = size(prob_l3_in_all,1);
        nSess = size(prob_l3_in_all,2);
        
        sub_t = nan(nSub,1);
        
        for sub_i = 1:nSub
        
            in_all  = [];
            out_all = [];
        
            for ses_i = 1:size(prob_l3_in_all,2)
        
                %% -------------------------
                % Skip empty sessions
                % -------------------------
                if isempty(prob_l3_in_all{sub_i,ses_i})
                    continue;
                end
        
                %% -------------------------
                % L3
                % -------------------------
                in_l3  = prob_l3_in_all{sub_i,ses_i};
                out_l3 = prob_l3_out_all{sub_i,ses_i};
        
                in_all  = [in_all;  in_l3(:)];
                out_all = [out_all; out_l3(:)];
        
                %% -------------------------
                % L4
                % -------------------------
                in_l4  = prob_l4_in_all{sub_i,ses_i};
                out_l4 = prob_l4_out_all{sub_i,ses_i};
        
                in_all  = [in_all;  in_l4(:)];
                out_all = [out_all; out_l4(:)];
        
            end
        
            %% -------------------------
            % LEVEL 1: within-subject test
            % -------------------------
            if numel(in_all) > 1 && numel(out_all) > 1
        
                [~,~,~,stats] = ttest2(in_all, out_all);
                sub_t(sub_i) = stats.tstat;
        
            end
        
        end
        
        %% =========================================================
        % LEVEL 2: across-subject inference
        %% =========================================================
        
        % group-level test
        [h,p,ci,stats] = ttest(sub_t);
        
        
        fprintf('     t = %.4f, p = %.6f\n', stats.tstat, p);
        
        
        
        
        %% model 02
        
        % probabilities above threshold
        
        %% =========================================================
        % HIERARCHICAL MODEL: P(in) > P(out) - WITH probability thresholding (p > 0.7)
        % Level 1: within-subject (pooled L3 + L4, and pooled sessions)
        % Level 2: across-subject t-test on subject t-values
        %% =========================================================
        
        fprintf('\n model 02 - probabilities above threshold\n');
            
        nSub  = size(prob_l3_in_all,1);
        sub_t = nan(nSub,1);
        p_thr = 0.5;
        
        for sub_i = 1:nSub
        
            in_all  = [];
            out_all = [];
        
            for ses_i = 1:size(prob_l3_in_all,2)
        
                %% -------------------------
                % L3
                % -------------------------
                in_l3  = prob_l3_in_all{sub_i,ses_i};
                out_l3 = prob_l3_out_all{sub_i,ses_i};
        
                in_l3  = in_l3(:);
                out_l3 = out_l3(:);
        
                % thresholding
                in_l3  = in_l3(in_l3  > p_thr);
                out_l3 = out_l3(out_l3 > p_thr);
        
                %% -------------------------
                % L4
                % -------------------------
                in_l4  = prob_l4_in_all{sub_i,ses_i};
                out_l4 = prob_l4_out_all{sub_i,ses_i};
        
                in_l4  = in_l4(:);
                out_l4 = out_l4(:);
        
                % thresholding
                in_l4  = in_l4(in_l4  > p_thr);
                out_l4 = out_l4(out_l4 > p_thr);
        
                %% -------------------------
                % Pool across sessions + lengths
                % -------------------------
                in_all  = [in_all;  in_l3; in_l4];
                out_all = [out_all; out_l3; out_l4];
        
            end
        
            %% -------------------------
            % LEVEL 1: within-subject test
            % -------------------------
            if numel(in_all) > 1 && numel(out_all) > 1
        
                [~,~,~,stats] = ttest2(in_all, out_all);
                sub_t(sub_i) = stats.tstat;
        
            end
        
        end
        
        %% =========================================================
        % LEVEL 2: across-subject inference
        %% =========================================================
        
        % group-level test
        [h,p,ci,stats] = ttest(sub_t);
        
        
        fprintf('     t = %.4f, p = %.6f\n', stats.tstat, p);
        
        
        
        
        %% model 03
        
        % probabilities below threshold
        
        %% =========================================================
        % HIERARCHICAL MODEL: P(in) > P(out) - WITH probability thresholding (p < 0.7)
        % Level 1: within-subject (pooled L3 + L4, and pooled sessions)
        % Level 2: across-subject t-test on subject t-values
        %% =========================================================
        
        fprintf('\n model 03 - probabilities below threshold\n');
        
        nSub  = size(prob_l3_in_all,1);
        sub_t = nan(nSub,1);
        p_thr = 0.5;
        
        for sub_i = 1:nSub
        
            in_all  = [];
            out_all = [];
        
            for ses_i = 1:size(prob_l3_in_all,2)
        
                %% -------------------------
                % L3
                % -------------------------
                in_l3  = prob_l3_in_all{sub_i,ses_i};
                out_l3 = prob_l3_out_all{sub_i,ses_i};
        
                in_l3  = in_l3(:);
                out_l3 = out_l3(:);
        
                % thresholding
                in_l3  = in_l3(in_l3  < p_thr);
                out_l3 = out_l3(out_l3 < p_thr);
        
                %% -------------------------
                % L4
                % -------------------------
                in_l4  = prob_l4_in_all{sub_i,ses_i};
                out_l4 = prob_l4_out_all{sub_i,ses_i};
        
                in_l4  = in_l4(:);
                out_l4 = out_l4(:);
        
                % thresholding
                in_l4  = in_l4(in_l4  < p_thr);
                out_l4 = out_l4(out_l4 < p_thr);
        
                %% -------------------------
                % Pool across sessions + lengths
                % -------------------------
                in_all  = [in_all;  in_l3; in_l4];
                out_all = [out_all; out_l3; out_l4];
        
            end
        
            %% -------------------------
            % LEVEL 1: within-subject test
            % -------------------------
            if numel(in_all) > 1 && numel(out_all) > 1
        
                [~,~,~,stats] = ttest2(in_all, out_all);
                sub_t(sub_i) = stats.tstat;
        
            end
        
        end
        
        %% =========================================================
        % LEVEL 2: across-subject inference
        %% =========================================================
        
        % group-level test
        [h,p,ci,stats] = ttest(sub_t);
        
        
        fprintf('     t = %.4f, p = %.6f\n', stats.tstat, p);
        
        
        
        
        
        
        
        
        %% model 04
        
        % probabilities above threshold + time segments
        
        %% =========================================================
        % HIERARCHICAL MODEL: P(in) > P(out) - WITH probability thresholding (p >
        % 0.7) & time segments of more than 10 above threshold timepoints
        % Level 1: within-subject (pooled L3 + L4, and pooled sessions)
        % Level 2: across-subject t-test on subject t-values
        %% =========================================================
        
        fprintf('\n model 04 - probabilities above threshold + time segments\n');
        
        nSub  = size(prob_l3_in_all,1);
        
        sub_t = nan(nSub,1);
        
        p_thr = 0.5;
        min_seg = 20;
        
        for sub_i = 1:nSub
        
            in_all  = [];
            out_all = [];
        
            for ses_i = 1:size(prob_l3_in_all,2)
        
                %% -------------------------
                % L3
                % -------------------------
                in_l3  = prob_l3_in_all{sub_i,ses_i};   % (trial × time × stim)
                out_l3 = prob_l3_out_all{sub_i,ses_i};
        
                in_all  = [in_all;  get_valid_segments(in_l3,  p_thr, min_seg)];
                out_all = [out_all; get_valid_segments(out_l3, p_thr, min_seg)];
        
                %% -------------------------
                % L4
                % -------------------------
                in_l4  = prob_l4_in_all{sub_i,ses_i};
                out_l4 = prob_l4_out_all{sub_i,ses_i};
        
                in_all  = [in_all;  get_valid_segments(in_l4,  p_thr, min_seg)];
                out_all = [out_all; get_valid_segments(out_l4, p_thr, min_seg)];
        
            end
        
            %% -------------------------
            % LEVEL 1: within-subject test
            %% -------------------------
            if numel(in_all) > 1 && numel(out_all) > 1
        
                [~,~,~,stats] = ttest2(in_all, out_all);
                sub_t(sub_i) = stats.tstat;
        
            end
        
        end
        
        %% =========================================================
        % LEVEL 2: across-subject inference
        %% =========================================================
        
        % group-level test
        [h,p,ci,stats] = ttest(sub_t);
        
        
        fprintf('     t = %.4f, p = %.6f\n', stats.tstat, p);
        
        
        
    end

end

%% summary of decoded state space - dval

for folder_i = 1:length(folders)

    folder_space = folders{folder_i};
    
    fprintf('\n==================\n');
    disp(['dval ' folder_space]);
    fprintf('==================\n');

    if ~isfile([path_predictions '/' folder_space '/dval_in_out.mat'])

        % preallocate cells
        dval_l3_in_all  = cell(length(subjects), length(sessions));
        dval_l3_out_all = cell(length(subjects), length(sessions));
        dval_l4_in_all  = cell(length(subjects), length(sessions));
        dval_l4_out_all = cell(length(subjects), length(sessions));
        
        for sub_i = 1:length(subjects)
        
            subject = subjects(sub_i);
        
            % loop over sessions
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
            
                % load data
                load([path_preprocessed_data '/' subject_ID '/' session_ID...
                        '/GERE_task_zthr25_wholetrial_baseline1sec_postICA.mat']); % data
        
                % load predictions
                load([path_predictions '/' subject_ID '/' session_ID '/' folder_space '/' inputs{1} '/' classifiers{1} '/predictions.mat']);
                
                % --- use number of trials from predictions ---
                nTrials = size(predictions.dval,1);
                nTime   = size(predictions.dval,2);
                nStim   = size(predictions.dval,3);
                
                %% -------------------------
                % Extract sequences (aligned to predictions)
                % -------------------------
                seq_all = cell(nTrials,1);
                seq_len = zeros(nTrials,1);
                
                for trial_i = 1:nTrials
                    seq = data.trialinfo(trial_i,4:7);
                    seq(seq == 0) = [];
                    seq_all{trial_i} = seq;
                    seq_len(trial_i) = length(seq);
                end
                
                % indices for each condition
                idx_l3 = find(seq_len == 3);
                idx_l4 = find(seq_len == 4);
                
                n_l3 = length(idx_l3);
                n_l4 = length(idx_l4);
                
                %% -------------------------
                % Preallocate tensors
                % -------------------------
                dval_l3_in  = zeros(n_l3, nTime, 3);
                dval_l3_out = zeros(n_l3, nTime, nStim-3);
                
                dval_l4_in  = zeros(n_l4, nTime, 4);
                dval_l4_out = zeros(n_l4, nTime, nStim-4);
                
                %% -------------------------
                % Fill tensors
                % -------------------------
                for k = 1:n_l3
                
                    trial_i = idx_l3(k);
                    sequence = seq_all{trial_i};
                
                    % --- IN (keep repeats, preserve order) ---
                    stim_in = sequence;
                
                    % --- OUT (from missing unique stimuli, trimmed) ---
                    stim_unique  = unique(sequence);
                    stim_missing = setdiff(1:nStim, stim_unique);
                
                    n_out = nStim - length(stim_in);
                    stim_out = stim_missing(1:n_out);
                
                    % assign
                    dval_l3_in(k,:,:)  = predictions.dval(trial_i,:,stim_in);
                    dval_l3_out(k,:,:) = predictions.dval(trial_i,:,stim_out);
                
                end
                
                for k = 1:n_l4
                
                    trial_i = idx_l4(k);
                    sequence = seq_all{trial_i};
                
                    stim_in = sequence;
                
                    stim_unique  = unique(sequence);
                    stim_missing = setdiff(1:nStim, stim_unique);
                
                    n_out = nStim - length(stim_in);
                    stim_out = stim_missing(1:n_out);
                
                    dval_l4_in(k,:,:)  = predictions.dval(trial_i,:,stim_in);
                    dval_l4_out(k,:,:) = predictions.dval(trial_i,:,stim_out);
                
                end
        
                % -------------------------
                % STORE (subject × session)
                % -------------------------
                dval_l3_in_all{sub_i, ses_i}  = dval_l3_in;
                dval_l3_out_all{sub_i, ses_i} = dval_l3_out;
        
                dval_l4_in_all{sub_i, ses_i}  = dval_l4_in;
                dval_l4_out_all{sub_i, ses_i} = dval_l4_out;        
                        
            end
        
        end
        
        dval_in_out = [];
        dval_in_out.dval_l3_in_all = dval_l3_in_all;
        dval_in_out.dval_l3_out_all = dval_l3_out_all;
        dval_in_out.dval_l4_in_all = dval_l4_in_all;
        dval_in_out.dval_l4_out_all = dval_l4_out_all;
        
        if ~isfolder([path_predictions '/' folder_space])
            mkdir([path_predictions '/' folder_space])
        end
        
        save([path_predictions '/' folder_space '/dval_in_out.mat'], 'dval_in_out');
    
    else

        load([path_predictions '/' folder_space '/dval_in_out.mat']); % prob_in_out
        dval_l3_in_all = dval_in_out.dval_l3_in_all;
        dval_l3_out_all = dval_in_out.dval_l3_out_all;
        dval_l4_in_all = dval_in_out.dval_l4_in_all;
        dval_l4_out_all = dval_in_out.dval_l4_out_all;

    end


    %% model 01
    
    % all probabilities
    
    %% =========================================================
    % HIERARCHICAL MODEL: P(in) > P(out)
    % Level 1: within-subject (pooled L3 + L4, and pooled sessions)
    % Level 2: across-subject t-test on subject t-values
    %% =========================================================
    
    fprintf('\n model 01 - all dvals\n');

    nSub  = size(dval_l3_in_all,1);
    nSess = size(dval_l3_in_all,2);
    
    sub_t = nan(nSub,1);
    
    for sub_i = 1:nSub
    
        in_all  = [];
        out_all = [];
    
        for ses_i = 1:size(dval_l3_in_all,2)
    
            %% -------------------------
            % Skip empty sessions
            % -------------------------
            if isempty(dval_l3_in_all{sub_i,ses_i})
                continue;
            end
    
            %% -------------------------
            % L3
            % -------------------------
            in_l3  = dval_l3_in_all{sub_i,ses_i};
            out_l3 = dval_l3_out_all{sub_i,ses_i};
    
            in_all  = [in_all;  in_l3(:)];
            out_all = [out_all; out_l3(:)];
    
            %% -------------------------
            % L4
            % -------------------------
            in_l4  = dval_l4_in_all{sub_i,ses_i};
            out_l4 = dval_l4_out_all{sub_i,ses_i};
    
            in_all  = [in_all;  in_l4(:)];
            out_all = [out_all; out_l4(:)];
    
        end
    
        %% -------------------------
        % LEVEL 1: within-subject test
        % -------------------------
        if numel(in_all) > 1 && numel(out_all) > 1
    
            [~,~,~,stats] = ttest2(in_all, out_all);
            sub_t(sub_i) = stats.tstat;
    
        end
    
    end
    
    %% =========================================================
    % LEVEL 2: across-subject inference
    %% =========================================================
    
    % group-level test
    [h,p,ci,stats] = ttest(sub_t);
    
    
    fprintf('     t = %.4f, p = %.6f\n', stats.tstat, p);
    
    
    
    
    %% model 02
    
    % probabilities above threshold
    
    %% =========================================================
    % HIERARCHICAL MODEL: P(in) > P(out) - WITH probability thresholding (p > 0.7)
    % Level 1: within-subject (pooled L3 + L4, and pooled sessions)
    % Level 2: across-subject t-test on subject t-values
    %% =========================================================
    
    fprintf('\n model 02 - dvals above threshold\n');
    
    nSub  = size(dval_l3_in_all,1);
    sub_t = nan(nSub,1);
    p_thr = 0;
    
    for sub_i = 1:nSub
    
        in_all  = [];
        out_all = [];
    
        for ses_i = 1:size(dval_l3_in_all,2)
    
            %% -------------------------
            % L3
            % -------------------------
            in_l3  = dval_l3_in_all{sub_i,ses_i};
            out_l3 = dval_l3_out_all{sub_i,ses_i};
    
            in_l3  = in_l3(:);
            out_l3 = out_l3(:);
    
            % thresholding
            in_l3  = in_l3(in_l3  > p_thr);
            out_l3 = out_l3(out_l3 > p_thr);
    
            %% -------------------------
            % L4
            % -------------------------
            in_l4  = dval_l4_in_all{sub_i,ses_i};
            out_l4 = dval_l4_out_all{sub_i,ses_i};
    
            in_l4  = in_l4(:);
            out_l4 = out_l4(:);
    
            % thresholding
            in_l4  = in_l4(in_l4  > p_thr);
            out_l4 = out_l4(out_l4 > p_thr);
    
            %% -------------------------
            % Pool across sessions + lengths
            % -------------------------
            in_all  = [in_all;  in_l3; in_l4];
            out_all = [out_all; out_l3; out_l4];
    
        end
    
        %% -------------------------
        % LEVEL 1: within-subject test
        % -------------------------
        if numel(in_all) > 1 && numel(out_all) > 1
    
            [~,~,~,stats] = ttest2(in_all, out_all);
            sub_t(sub_i) = stats.tstat;
    
        end
    
    end
    
    %% =========================================================
    % LEVEL 2: across-subject inference
    %% =========================================================
    
    % group-level test
    [h,p,ci,stats] = ttest(sub_t);
    
    
    fprintf('     t = %.4f, p = %.6f\n', stats.tstat, p);
    
    
    
    
    %% model 03
    
    % probabilities below threshold
    
    %% =========================================================
    % HIERARCHICAL MODEL: P(in) > P(out) - WITH probability thresholding (p < 0.7)
    % Level 1: within-subject (pooled L3 + L4, and pooled sessions)
    % Level 2: across-subject t-test on subject t-values
    %% =========================================================
    
    fprintf('\n model 03 - dvals below threshold\n');
    
    nSub  = size(dval_l3_in_all,1);
    sub_t = nan(nSub,1);
    p_thr = 0;
    
    for sub_i = 1:nSub
    
        in_all  = [];
        out_all = [];
    
        for ses_i = 1:size(dval_l3_in_all,2)
    
            %% -------------------------
            % L3
            % -------------------------
            in_l3  = dval_l3_in_all{sub_i,ses_i};
            out_l3 = dval_l3_out_all{sub_i,ses_i};
    
            in_l3  = in_l3(:);
            out_l3 = out_l3(:);
    
            % thresholding
            in_l3  = in_l3(in_l3  < p_thr);
            out_l3 = out_l3(out_l3 < p_thr);
    
            %% -------------------------
            % L4
            % -------------------------
            in_l4  = dval_l4_in_all{sub_i,ses_i};
            out_l4 = dval_l4_out_all{sub_i,ses_i};
    
            in_l4  = in_l4(:);
            out_l4 = out_l4(:);
    
            % thresholding
            in_l4  = in_l4(in_l4  < p_thr);
            out_l4 = out_l4(out_l4 < p_thr);
    
            %% -------------------------
            % Pool across sessions + lengths
            % -------------------------
            in_all  = [in_all;  in_l3; in_l4];
            out_all = [out_all; out_l3; out_l4];
    
        end
    
        %% -------------------------
        % LEVEL 1: within-subject test
        % -------------------------
        if numel(in_all) > 1 && numel(out_all) > 1
    
            [~,~,~,stats] = ttest2(in_all, out_all);
            sub_t(sub_i) = stats.tstat;
    
        end
    
    end
    
    %% =========================================================
    % LEVEL 2: across-subject inference
    %% =========================================================
    
    % group-level test
    [h,p,ci,stats] = ttest(sub_t);
    
    
    fprintf('     t = %.4f, p = %.6f\n', stats.tstat, p);
    
    
    
    
    
    
    
    
    %% model 04
    
    % probabilities above threshold + time segments
    
    %% =========================================================
    % HIERARCHICAL MODEL: P(in) > P(out) - WITH probability thresholding (p >
    % 0.7) & time segments of more than 10 above threshold timepoints
    % Level 1: within-subject (pooled L3 + L4, and pooled sessions)
    % Level 2: across-subject t-test on subject t-values
    %% =========================================================
    
    fprintf('\n model 04 - dvals above threshold + time segments\n');
    
    nSub  = size(dval_l3_in_all,1);
    
    sub_t = nan(nSub,1);
    
    p_thr = 0;
    min_seg = 20;
    
    for sub_i = 1:nSub
    
        in_all  = [];
        out_all = [];
    
        for ses_i = 1:size(dval_l3_in_all,2)
    
            %% -------------------------
            % L3
            % -------------------------
            in_l3  = dval_l3_in_all{sub_i,ses_i};   % (trial × time × stim)
            out_l3 = dval_l3_out_all{sub_i,ses_i};
    
            in_all  = [in_all;  get_valid_segments(in_l3,  p_thr, min_seg)];
            out_all = [out_all; get_valid_segments(out_l3, p_thr, min_seg)];
    
            %% -------------------------
            % L4
            % -------------------------
            in_l4  = dval_l4_in_all{sub_i,ses_i};
            out_l4 = dval_l4_out_all{sub_i,ses_i};
    
            in_all  = [in_all;  get_valid_segments(in_l4,  p_thr, min_seg)];
            out_all = [out_all; get_valid_segments(out_l4, p_thr, min_seg)];
    
        end
    
        %% -------------------------
        % LEVEL 1: within-subject test
        %% -------------------------
        if numel(in_all) > 1 && numel(out_all) > 1
    
            [~,~,~,stats] = ttest2(in_all, out_all);
            sub_t(sub_i) = stats.tstat;
    
        end
    
    end
    
    %% =========================================================
    % LEVEL 2: across-subject inference
    %% =========================================================
    
    % group-level test
    [h,p,ci,stats] = ttest(sub_t);
    
    
    fprintf('     t = %.4f, p = %.6f\n', stats.tstat, p);
    
    


end



disp('Analysis done');


%% helper functions

function vals_out = get_valid_segments(data, p_thr, min_seg)

% data: (trial × time × stim)

vals_out = [];

[nTrial, nTime, nStim] = size(data);

for tr = 1:nTrial
    for st = 1:nStim

        ts = squeeze(data(tr,:,st));  % (time)

        % threshold mask
        mask = ts > p_thr;

        % find contiguous segments
        d = diff([0 mask 0]);
        start_idx = find(d == 1);
        end_idx   = find(d == -1) - 1;

        for i = 1:length(start_idx)

            seg_len = end_idx(i) - start_idx(i) + 1;

            if seg_len >= min_seg

                seg = ts(start_idx(i):end_idx(i));
                vals_out = [vals_out; seg(:)];

            end

        end

    end
end

end