%% GERE project

% path to data
if isfolder('/path_to_local')
    path_root = '/path_to_local';
end

path_trialinfo = [path_root '/results/preprocessing/trialinfo'];
path_predictions = [path_root '/github/results/neural_sequences/data/predictions'];
path_sensitivity = [path_root '/github/results/neural_sequences/sensitivity_analysis'];
addpath([path_root '/github/scripts/neural_sequences/utilities']);

number_stimuli = 8;
classifiers = {'l2_e-1'}; classifier = classifiers{1};
time_delay = 791; % bins
subjects = 1:17;
sessions = 1:2;
number_iterations = 50;
nlags = 200;

folder_replay = 'sequence_trial_avg_delay';
folders = {'task_stim2delay', 'localizer2delay'};
amplitudes = [1.5 1.4 1.3 1.2 1.1 1 0.9 0.8 0.7 0.6 0.5];
lags = [100, 200, 300]; % milliseconds
percentiles = [85, 90, 95];

if strcmp(path_root, '/scratch/as14864/GERE')
    folders = folders(folder_ii);
    amplitudes = amplitudes(amplitude_ii);
    percentiles = percentiles(percentile_ii);
    lags = lags(lag_ii);
end

for percentile_i = 1:length(percentiles)

    percentile = percentiles(percentile_i);

    for folder_i = 1:length(folders)
    
        folder_space = folders{folder_i};
        
        for amplitude_i = 1:length(amplitudes)
    
            injected_percentage = amplitudes(amplitude_i);
            
            for lag_i = 1:length(lags)
    
                lag_ms = lags(lag_i);
    
                for iter_i = 1:number_iterations
    
                    disp(['percentile ' num2str(percentile) ' - ' folder_space ' - injected ' num2str(injected_percentage *100) ' - lag ' num2str(lag_ms) ' - iter ' num2str(iter_i)])
           
                    disp(datetime);

                    outdir = [path_sensitivity '/' folder_space '/percentile_' num2str(percentile) '/injected_' num2str(injected_percentage * 100) '/lag_' num2str(lag_ms) '/iter_' num2str(iter_i)];
                         
                    if ~isfolder(outdir)
                        mkdir(outdir);
                    end
                    
                    if ~isfile([outdir '/injected_significance_length3.mat'])
                    
                        %% --------------------
                        % find idx from the temporal lag
                        %  --------------------
        
                        disp('  find idx from the temporal lag');
        
                        disp(datetime);
                        
                        % Parameters
                        timepoints = 4000;   % delay duration in ms
                        bin_length = 50;     % window length in ms
                        bin_step = 5;        % step size in ms
                        
                        bin_number = round(timepoints / bin_step) - (bin_step * 2) + 1;
                        
                        % Create bins
                        bins = zeros(bin_number, bin_length);
                        
                        for bin_i = 1:bin_number
                            if bin_i == 1
                                bins(bin_i,:) = 1:bin_length;
                            else
                                bins(bin_i,:) = bins(bin_i-1,:) + bin_step;
                            end
                        end
                        
                        % Bin centers in milliseconds
                        time_test = mean(bins, 2);
                        
                        % Find bin corresponding to a desired time                        
                        [~, lag_idx] = min(abs(time_test - lag_ms));
                        
                        %% --------------------
                        % pool predictions
                        %  --------------------
        
                        disp('  pool predictions');
                        
                        disp(datetime);
                        
                        % preallocate cells
                        dval_l3_in  = cell(length(subjects), length(sessions));
                        dval_l4_in  = cell(length(subjects), length(sessions));
                        dval_l3_out  = cell(length(subjects), length(sessions));
                        dval_l4_out  = cell(length(subjects), length(sessions));
        
                        % loop over subjects
                        for sub_i = 1:length(subjects)
                        
                            subject = subjects(sub_i);
                        
                            for ses_i = 1:length(sessions)
                        
                                session = sessions(ses_i);
                        
                                %% load data
                        
                                % set paths
                                if subject < 10
                                    subject_ID = ['sub_0' num2str(subject)];
                                else
                                    subject_ID = ['sub_' num2str(subject)];
                                end
                                                
                                session_ID = ['sess_0' num2str(session)];
                                    
                                % load data
                                load([path_trialinfo '/' subject_ID '/' session_ID...
                                        '/data_trialinfo.mat']); % data_trialinfo
        
                                data = data_trialinfo; clear data_trialinfo
                
                                % load predictions
                                load([path_predictions '/' subject_ID '/' session_ID '/' folder_space '/' classifiers{1} '/predictions.mat']);
                                        
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
                                dval_l3_in_trials  = zeros(n_l3, nTime, 3);                
                                dval_l4_in_trials  = zeros(n_l4, nTime, 4);
                                
                                dval_l3_out_trials  = zeros(n_l3, nTime, 5);                
                                dval_l4_out_trials  = zeros(n_l4, nTime, 4);
                                
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
                                    dval_l3_in_trials(k,:,:)  = predictions.dval(trial_i,:,stim_in);
                                    dval_l3_out_trials(k,:,:)  = predictions.dval(trial_i,:,stim_out);
                                    
                                end
                                
                                for k = 1:n_l4
                                
                                    trial_i = idx_l4(k);
                                    sequence = seq_all{trial_i};
                                
                                    stim_in = sequence;
                                
                                    stim_unique  = unique(sequence);
                                    stim_missing = setdiff(1:nStim, stim_unique);
                                
                                    n_out = nStim - length(stim_in);
                                    stim_out = stim_missing(1:n_out);
                                
                                    dval_l4_in_trials(k,:,:)  = predictions.dval(trial_i,:,stim_in);
                                    dval_l4_out_trials(k,:,:)  = predictions.dval(trial_i,:,stim_out);
                                    
                                end
                        
                                % -------------------------
                                % STORE (subject × session)
                                % -------------------------
                                dval_l3_in{sub_i, ses_i}  = dval_l3_in_trials;
                                dval_l4_in{sub_i, ses_i}  = dval_l4_in_trials;
                                     
                                dval_l3_out{sub_i, ses_i}  = dval_l3_out_trials;
                                dval_l4_out{sub_i, ses_i}  = dval_l4_out_trials;
                                
                            end
                        
                        end
        
        
                        %% --------------------
                        % in vs out dval - time-resolved
                        % percentile of IN-vs-OUT differences
                        % --------------------
                        
                        disp('  in vs out dval - time-resolved');
                        disp(datetime);
                        
                        dval_in_vs_out_empirical = nan(length(subjects), length(sessions));
                        
                        % Loop over subjects
                        for sub_i = 1:length(subjects)
                        
                            % Loop over sessions
                            for ses_i = 1:length(sessions)
                        
                                %% -------------------------
                                % L3
                                % -------------------------
                        
                                in_l3  = dval_l3_in{sub_i,ses_i};
                                out_l3 = dval_l3_out{sub_i,ses_i};
                        
                                % Average OUT stimuli within each trial and time bin
                                % 40 x 791 x 5 --> 40 x 791
                                out_l3_mean = mean(out_l3, 3);
                        
                                % IN minus average OUT
                                % 40 x 791 x 3
                                dval_diff_l3 = in_l3 - out_l3_mean;
                        
                                % Pool trials, time bins and IN stimuli
                                dval_diff_l3_all = dval_diff_l3(:);
                        
                                % percentile
                                dval_effect_l3 = prctile(dval_diff_l3_all, percentile);
                        
                        
                                %% -------------------------
                                % L4
                                % -------------------------
                        
                                in_l4  = dval_l4_in{sub_i,ses_i};
                                out_l4 = dval_l4_out{sub_i,ses_i};
                        
                                % Average OUT stimuli within each trial and time bin
                                % 40 x 791 x 4 --> 40 x 791
                                out_l4_mean = mean(out_l4, 3);
                        
                                % IN minus average OUT
                                % 40 x 791 x 4
                                dval_diff_l4 = in_l4 - out_l4_mean;
                        
                                % Pool trials, time bins and IN stimuli
                                dval_diff_l4_all = dval_diff_l4(:);
                        
                                % percentile
                                dval_effect_l4 = prctile(dval_diff_l4_all, percentile);
                        
                        
                                %% -------------------------
                                % Combine L3 and L4
                                % -------------------------
                        
                                dval_in_vs_out_empirical(sub_i,ses_i) = ...
                                    mean([dval_effect_l3 dval_effect_l4]);
                        
                            end
                        end
                        
                        
        
        
                        %% --------------------
                        % inject sequences
                        %  --------------------
        
                        disp('  inject sequences');
                        
                        disp(datetime);
                        
                        % Loop over subjects
                        for sub_i = 1:length(subjects)
                        
                            subject = subjects(sub_i);
                                        
                            % Loop over sessions
                            for ses_i = 1:length(sessions)
                
                                session = sessions(ses_i);
                                        
                                % Injection amplitude
                                injected_amplitude = dval_in_vs_out_empirical(sub_i) * injected_percentage;


                                %%% length 3
                                
                                sequence_length = 3;
                                
                                % Sequence duration in bins
                                sequence_duration = lag_idx * (sequence_length - 1);
                                
                                % Valid onset range
                                first_onset = 5;
                                last_onset = 791 - sequence_duration - 3;
                                
                                % Number of sequences to inject per trial
                                n_sequences = 20;
                                
                                % Loop over trials
                                for trial_i = 1:size(dval_l3_in{sub_i,ses_i}, 1)
                                
                                    % ---------------------------------------------------------
                                    % Randomize sequence onsets independently for this trial
                                    % ---------------------------------------------------------
                                
                                    % Minimum separation between sequence onsets
                                    % Choose a value that is not a multiple of lag_idx
                                    min_onset_separation = (lag_idx * sequence_length) + randi([-10 10]);
                                    
                                    candidate_onsets = first_onset:last_onset;
                                    candidate_onsets = candidate_onsets(randperm(length(candidate_onsets)));
                                
                                    onsets = [];
                                
                                    for onset_idx = candidate_onsets
                                
                                        if isempty(onsets) || ...
                                                all(abs(onset_idx - onsets) >= min_onset_separation)
                                
                                            onsets(end+1) = onset_idx;
                                
                                        end
                                
                                        % Stop once we have the required number
                                        if length(onsets) == n_sequences
                                            break
                                        end
                                
                                    end
                                
                                    % Sort chronologically
                                    onsets = sort(onsets);
                                
                                    % ---------------------------------------------------------
                                    % Inject sequences
                                    % ---------------------------------------------------------
                                
                                    for onset_idx = onsets
                                
                                        % Inject L3 sequence: 1 -> 2 -> 3
                                        for state_i = 1:sequence_length
                                
                                            time_idx = onset_idx + lag_idx * (state_i - 1);
                                
                                            time_window = (time_idx-3):(time_idx+3);
                                
                                            dval_l3_in{sub_i,ses_i}(trial_i,time_window,state_i) = ...
                                                dval_l3_in{sub_i,ses_i}(trial_i,time_window,state_i) + ...
                                                injected_amplitude;
                                
                                        end
                                
                                    end
                                
                                end



                                %%% length 4
                                
                                sequence_length = 4;
                                
                                % Sequence duration in bins
                                sequence_duration = lag_idx * (sequence_length - 1);
                                
                                % Valid onset range
                                first_onset = 5;
                                last_onset = 791 - sequence_duration - 3;
                                
                                % Number of sequences to inject per trial
                                n_sequences = 20;
                                
                                % Loop over trials
                                for trial_i = 1:size(dval_l4_in{sub_i,ses_i}, 1)
                                
                                    % ---------------------------------------------------------
                                    % Randomize sequence onsets independently for this trial
                                    % ---------------------------------------------------------
                                
                                    % Minimum separation between sequence onsets
                                    % Choose a value that is not a multiple of lag_idx
                                    min_onset_separation = (lag_idx * sequence_length) + randi([-10 10]);
                                    
                                    candidate_onsets = first_onset:last_onset;
                                    candidate_onsets = candidate_onsets(randperm(length(candidate_onsets)));
                                
                                    onsets = [];
                                
                                    for onset_idx = candidate_onsets
                                
                                        if isempty(onsets) || ...
                                                all(abs(onset_idx - onsets) >= min_onset_separation)
                                
                                            onsets(end+1) = onset_idx;
                                
                                        end
                                
                                        % Stop once we have the required number
                                        if length(onsets) == n_sequences
                                            break
                                        end
                                
                                    end
                                
                                    % Sort chronologically
                                    onsets = sort(onsets);
                                
                                    % ---------------------------------------------------------
                                    % Inject sequences
                                    % ---------------------------------------------------------
                                
                                    for onset_idx = onsets
                                
                                        % Inject L3 sequence: 1 -> 2 -> 3 -> 4
                                        for state_i = 1:sequence_length
                                
                                            time_idx = onset_idx + lag_idx * (state_i - 1);
                                
                                            time_window = (time_idx-3):(time_idx+3);
                                
                                            dval_l4_in{sub_i,ses_i}(trial_i,time_window,state_i) = ...
                                                dval_l4_in{sub_i,ses_i}(trial_i,time_window,state_i) + ...
                                                injected_amplitude;
                                
                                        end
                                
                                    end
                                
                                end

                            end
                
                        end
                        
                        %% --------------------
                        % TDLM
                        % --------------------
        
                        disp('  TDLM');
                        
                        disp(datetime);
                        
                        % Loop over subjects
                        for sub_i = 1:length(subjects)
                        
                            subject = subjects(sub_i);
                                        
                            % Loop over sessions
                            for ses_i = 1:length(sessions)
                
                                session = sessions(ses_i);
        
                                % set paths
                                if subject < 10
                                    subject_ID = ['sub_0' num2str(subject)];
                                else
                                    subject_ID = ['sub_' num2str(subject)];
                                end
                                        
                                session_ID = ['sess_0' num2str(session)];
        
                                % loop over sequence lengths
                                sequence_lengths = [3, 4];
                                for seq_i = 1:length(sequence_lengths)
        
                                    if sequence_lengths(seq_i) == 3
                                        predictions.dval = dval_l3_in{sub_i, ses_i};
                                    else
                                        predictions.dval = dval_l4_in{sub_i, ses_i};
                                    end
                        
                                    % sigmoid function
                                    predictions.prob = (1 ./ (1 + exp(-predictions.dval)));
                                        
                                    %% sequenceness and null distributions
                                    
                                    Msf_alltrials = cell(0);
                                    Msb_alltrials = cell(0);
                                    
                                    for trial_i = 1:size(predictions.prob, 1)
                                                        
                                        %% TDLM - sequence
                        
                                        % settings TDLM
                                        nstates = sequence_lengths(seq_i);
                                        nbins=nlags+1;
                        
                                        % input - decoded state space
                                        X = squeeze(predictions.prob(trial_i,:,:)); % time by states
                        
                                        % compute shifted version of state space X(t + Dt)
                                        warning off
                        
                                        dm=[toeplitz(X(:,1),[zeros(nbins,1)])];
                                        dm=dm(:,2:end);                         
                                        
                                        for kk=2:nstates
                                            temp=toeplitz(X(:,kk),[zeros(nbins,1)]);
                                            temp=temp(:,2:end);
                                            dm=[dm temp]; 
                                        end
                                        % TOEPLITZ(C,R) is non-symmetric Toeplitz matrix having C as its first column and R as its first row
                                        % matrix dm contains the shifted versions of the original state columns of
                                        % X, arranged in a way that reflects the Toeplitz structure,
                                        % dm is X(t + Dt) with dimensions (time by nlags *
                                        % states), shifting columns down and adding zeros on
                                        % top (shift up time)
                        
                                        % dm dimensions: time by states*nlags, where
                                        % columns are organized as follow: all lags of
                                        % state 1, then all lags of state 2,...
                        
                                        warning on
                        
                                        %% TDLM
                                        
                                        % X = state space
                        
                                        % transition matrix from sequence (deleting
                                        % empty rows/columns, and ordering by
                                        % appearence in sequence, e.g. [2 5 4] is
                                        % encoded as [1 2 3], where [1==2, 2==5,
                                        % 3==4]
                                        T = zeros(nstates,nstates);
                                        for item_i = 2:nstates
                                            T(item_i-1, item_i) = 1;
                                        end
                                        
                                        % number shuffle for null distribution
                                        nshuf = 24; % maximum permutations of 4 states
                        
                                        [Msf, Msb] = TDLM_trial(X,T,nstates,nlags,nshuf);
                                
                                        Msf_alltrials{end+1} = Msf;
                                        Msb_alltrials{end+1} = Msb;
                        
                                    end
        
                        
                                    %% average across permutations for each time lag
        
                                    %%% average across trials
        
                                    % empty matrix for permutations
                                    Msf = nan(size(Msf_alltrials{1}, 2), size(Msf_alltrials{1}, 3), size(Msf_alltrials, 2)); % permutation, time lag, trial
                                    Msb = nan(size(Msb_alltrials{1}, 2), size(Msb_alltrials{1}, 3), size(Msb_alltrials, 2)); % permutation, time lag, trial
        
                                    % pool trials
                                    for trial_i = 1:size(Msf_alltrials, 2)
        
                                        Msf(:,:,trial_i) = squeeze(Msf_alltrials{trial_i});
                                        Msb(:,:,trial_i) = squeeze(Msb_alltrials{trial_i});
        
                                    end
        
                                    % average across trials
                                    Msf = median(Msf, 3);
                                    Msb = median(Msb, 3);
                        
                                    % empirical
                                    Msf_empirical = Msf(1,:);
                                    Msb_empirical = Msb(1,:);
        
                                    % null
                                    Msf_null = Msf(2:end,:);
                                    Msb_null = Msb(2:end,:);
                                    if sequence_lengths(seq_i) == 3
                                        Msf_null = Msf_null(1:5,:);
                                        Msb_null = Msb_null(1:5,:);
                                    end
        
                                    if sequence_lengths(seq_i) == 3
        
                                        Msf_empirical_length3 = Msf_empirical;
                                        Msb_empirical_length3 = Msb_empirical;
                                        Msf_null_length3 = Msf_null;
                                        Msb_null_length3 = Msb_null;
                                    else
        
                                        Msf_empirical_length4 = Msf_empirical;
                                        Msb_empirical_length4 = Msb_empirical;
                                        Msf_null_length4 = Msf_null;
                                        Msb_null_length4 = Msb_null;
        
                                    end
        
                                end
        
                                %% save
                    
                                outdir_subj = [path_sensitivity '/' folder_space '/percentile_' num2str(percentile) '/injected_' num2str(injected_percentage * 100) '/lag_' num2str(lag_ms) '/iter_' num2str(iter_i) '/tmp/' subject_ID '/' session_ID];
                                
                                if ~isfolder(outdir_subj)
                                    mkdir(outdir_subj);
                                end
                                            
                                save([outdir_subj '/Msf_empirical_length3.mat'], 'Msf_empirical_length3', '-v7.3');
                                save([outdir_subj '/Msb_empirical_length3.mat'], 'Msb_empirical_length3', '-v7.3');
                        
                                save([outdir_subj '/Msf_empirical_length4.mat'], 'Msf_empirical_length4', '-v7.3');
                                save([outdir_subj '/Msb_empirical_length4.mat'], 'Msb_empirical_length4', '-v7.3');
                                
                                save([outdir_subj '/Msf_null_length3.mat'], 'Msf_null_length3', '-v7.3');
                                save([outdir_subj '/Msb_null_length3.mat'], 'Msb_null_length3', '-v7.3');
                    
                                save([outdir_subj '/Msf_null_length4.mat'], 'Msf_null_length4', '-v7.3');
                                save([outdir_subj '/Msb_null_length4.mat'], 'Msb_null_length4', '-v7.3');
        
                            end
        
                        end
        
                        %% group statistics on sequenceness
        
                        disp('  TDLM - stats');
        
                        disp(datetime);
                        
                        % loop over sequence lengths
                        sequence_lengths = [3, 4];
                        for seq_i = 1:length(sequence_lengths)
        
                            % outputs
                            null_Msf_thr_pos_allsubjects = [];
                            null_Msf_thr_neg_allsubjects = [];
                            null_Msb_thr_pos_allsubjects = [];
                            null_Msb_thr_neg_allsubjects = [];
                            
                            % Loop over subjects
                            for sub_i = 1:length(subjects)
                            
                                subject = subjects(sub_i);
                                            
                                % Loop over sessions
                                for ses_i = 1:length(sessions)
                    
                                    session = sessions(ses_i);
        
                                    % set paths
                                    if subject < 10
                                        subject_ID = ['sub_0' num2str(subject)];
                                    else
                                        subject_ID = ['sub_' num2str(subject)];
                                    end
                                            
                                    session_ID = ['sess_0' num2str(session)];
        
                                    if sequence_lengths(seq_i) == 3
                                        
                                        load([outdir_subj '/Msf_null_length3.mat']);
                                        load([outdir_subj '/Msb_null_length3.mat']);
                                                            
                                        Msf_null = Msf_null_length3;
                                        Msb_null = Msb_null_length3;
        
                                    else
        
                                        load([outdir_subj '/Msf_null_length4.mat']);
                                        load([outdir_subj '/Msb_null_length4.mat']);
                                        
                                        Msf_null = Msf_null_length4;
                                        Msb_null = Msb_null_length4;
        
                                    end
        
                                    % maximum (first column) and minimum (second column) across time lags for each permutation
                                    null_Msf = [max(Msf_null, [],2) min(Msf_null, [],2)];
                                    null_Msb = [max(Msb_null, [],2) min(Msb_null, [],2)];
        
                                    % threshold (P_fwe = 0.05)
                                    null_Msf_thr_pos = prctile(null_Msf(:,1), 97.5);
                                    null_Msf_thr_neg = prctile(null_Msf(:,2), 2.5);
                                    null_Msb_thr_pos = prctile(null_Msb(:,1), 97.5);
                                    null_Msb_thr_neg = prctile(null_Msb(:,2), 2.5);
        
                                    null_Msf_thr_pos_allsubjects = [null_Msf_thr_pos_allsubjects null_Msf_thr_pos];
                                    null_Msf_thr_neg_allsubjects = [null_Msf_thr_neg_allsubjects null_Msf_thr_neg];
                                    null_Msb_thr_pos_allsubjects = [null_Msb_thr_pos_allsubjects null_Msb_thr_pos];
                                    null_Msb_thr_neg_allsubjects = [null_Msb_thr_neg_allsubjects null_Msb_thr_neg];
        
                                end
        
                            end
        
                            if sequence_lengths(seq_i) == 3
                                null_Msf_thr_pos_l3 = median(null_Msf_thr_pos_allsubjects);
                                null_Msf_thr_neg_l3 = median(null_Msf_thr_neg_allsubjects);
                                null_Msb_thr_pos_l3 = median(null_Msb_thr_pos_allsubjects);
                                null_Msb_thr_neg_l3 = median(null_Msb_thr_neg_allsubjects);
                            else
                                null_Msf_thr_pos_l4 = median(null_Msf_thr_pos_allsubjects);
                                null_Msf_thr_neg_l4 = median(null_Msf_thr_neg_allsubjects);
                                null_Msb_thr_pos_l4 = median(null_Msb_thr_pos_allsubjects);
                                null_Msb_thr_neg_l4 = median(null_Msb_thr_neg_allsubjects);
                            end
        
                        end
        
                        %%% significance of injected sequences
        
                        if Msf_empirical_length3(lag_idx) > null_Msf_thr_pos_l3
                            injected_significance_length3 = 1;
                        else
                            injected_significance_length3 = 0;
                        end
        
                        if Msf_empirical_length4(lag_idx) > null_Msf_thr_pos_l4
                            injected_significance_length4 = 1;
                        else
                            injected_significance_length4 = 0;
                        end
        
                        %% save                                    
                        save([outdir '/injected_significance_length3.mat'], 'injected_significance_length3', '-v7.3');
                        save([outdir '/injected_significance_length4.mat'], 'injected_significance_length4', '-v7.3');

                        save([outdir '/Msf_empirical_length3.mat'], 'Msf_empirical_length3', '-v7.3');
                        save([outdir '/Msf_empirical_length4.mat'], 'Msf_empirical_length4', '-v7.3');

                        save([outdir '/null_Msf_thr_pos_l3.mat'], 'null_Msf_thr_pos_l3', '-v7.3');
                        save([outdir '/null_Msf_thr_pos_l4.mat'], 'null_Msf_thr_pos_l4', '-v7.3');
                        
                        % delete intermediate files
                        outdir_tmp = [path_sensitivity '/' folder_space '/percentile_' num2str(percentile) '/injected_' num2str(injected_percentage * 100) '/lag_' num2str(lag_ms) '/iter_' num2str(iter_i) '/tmp'];
        
                        system(['rm -r ' outdir_tmp]);

                    end
    
                end
    
            end
    
        end
    
    end
    
    disp(datetime);

end

disp('  Analysis done');
