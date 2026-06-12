%% GERE project
% replay

folder_replay = 'sequence_trial_avg_stim';

if isfolder('/path_to_local')
    path_root = '/path_to_local';
    subjects = 1:17;
end

path_data = [path_root, '/github/results/replay/data'];
path_results = [path_root, '/github/results/replay/' folder_replay];
addpath([path_root, '/github/scripts/replay/utilities'])

classifier = 'l2_e-1';
nlags = 200;
number_stimuli = 8;
folders = {'task_stim2task_stim'...
           'localizer2task_stim'};
sessions = 1:2;

%% replay

for folder_i = 1:length(folders)

    disp(folders{folder_i});
            
    for sub_i = 1:length(subjects)
    
        subject = subjects(sub_i);
    
        % set paths
        if subject < 10
            subject_ID = ['sub_0' num2str(subject)];
            subjectID = ['sub0' num2str(subject)];
        else
            subject_ID = ['sub_' num2str(subject)];
            subjectID = ['sub' num2str(subject)];
        end

        disp(subject_ID);
                    
        for ses_i = 1:length(sessions)
    
            session = sessions(ses_i);
                                                            
            session_ID = ['sess_0' num2str(session)];
        
            outdir = [path_results '/' subject_ID '/' session_ID '/' classifier '/' folders{folder_i}];
        
            if ~isfolder(outdir)
                mkdir(outdir);
            end
            
            %% load sequences
    
            indir = [path_data '/sequences/' subject_ID '/' session_ID];
    
            load([indir '/sequences.mat']); % sequences

            %% load predictions

            indir = [path_data '/predictions/' subject_ID '/' session_ID '/' folders{folder_i} '/' classifier];
            load([indir '/predictions.mat']); % predictions

            % sigmoid function
            predictions.prob = (1 ./ (1 + exp(-predictions.dval)));
                
            %% sequenceness and null distributions
            
            Msf_alltrials = cell(0);
            Msb_alltrials = cell(0);
            
            for trial_i = 1:size(predictions.prob, 1)

                disp([subject_ID ' - ' session_ID ' - trial ' num2str(trial_i)])

                % subset predictions for trial in the sequence
                sequence = sequences.trialinfo(trial_i, 4:7);
                sequence = sequence(sequence > 0);

                %% TDLM - sequence

                % settings TDLM
                nstates = length(sequence);
                nbins=nlags+1;

                % input - decoded state space
                X = squeeze(predictions.prob(trial_i,:,sequence)); % time by states

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
                for item_i = 2:length(sequence)
                    T(item_i-1, item_i) = 1;
                end

                % number of states (for this trial)
                nstates = length(sequence);

                % number shuffle for null distribution
                nshuf = 24; % maximum permutations of 4 states

                [Msf, Msb] = TDLM_trial(X,T,nstates,nlags,nshuf);
        
                Msf_alltrials{end+1} = Msf;
                Msb_alltrials{end+1} = Msb;

            end

            %% average across permutations for each time lag

            % empirical sequenceness
            Msf_empirical = Msf_alltrials;
            Msb_empirical = Msb_alltrials;

            for trial_i = 1:size(Msf_alltrials, 2)
                Msf_empirical{trial_i} = squeeze(Msf_empirical{trial_i}(1,1,:));
                Msb_empirical{trial_i} = squeeze(Msb_empirical{trial_i}(1,1,:));
            end

            % empty matrix for permutations
            Msf = nan(size(Msf_alltrials{1}, 2), size(Msf_alltrials{1}, 3), size(Msf_alltrials, 2)); % permutation, time lag, trial
            Msb = nan(size(Msb_alltrials{1}, 2), size(Msb_alltrials{1}, 3), size(Msb_alltrials, 2)); % permutation, time lag, trial
            
            % pool trials
            for trial_i = 1:size(Msf_alltrials, 2)
                
                Msf(:,:,trial_i) = squeeze(Msf_alltrials{trial_i});
                Msb(:,:,trial_i) = squeeze(Msb_alltrials{trial_i});
            
            end

            % subset by sequence length (and remove first
            % permutation, that is, empirical)
            Msf_length3 = nan(5, size(Msf, 2), size(Msf, 3)); % permutations by lags by trials
            Msb_length3 = nan(5, size(Msf, 2), size(Msf, 3));

            Msf_length4 = nan(23, size(Msf, 2), size(Msf, 3));
            Msb_length4 = nan(23, size(Msf, 2), size(Msf, 3));

            trials_length3 = [];
            trials_length4 = [];
            for trial_i = 1:size(Msf_alltrials, 2)
                if isnan(Msf(10,1,trial_i)) % length 3
                    trials_length3 = [trials_length3 trial_i];
                    Msf_length3(:,:,trial_i) = Msf(2:6,:,trial_i);
                    Msb_length3(:,:,trial_i) = Msb(2:6,:,trial_i);
                else
                    trials_length4 = [trials_length4 trial_i];
                    Msf_length4(:,:,trial_i) = Msf(2:24,:,trial_i);
                    Msb_length4(:,:,trial_i) = Msb(2:24,:,trial_i);
                end
            end
    
            Msf_length3 = Msf_length3(:,:,trials_length3);
            Msb_length3 = Msb_length3(:,:,trials_length3);
            
            Msf_length4 = Msf_length4(:,:,trials_length4);
            Msb_length4 = Msb_length4(:,:,trials_length4);
            
            % average permutations across trials
            Msf_length3 = median(Msf_length3, 3);
            Msb_length3 = median(Msb_length3, 3);
            
            Msf_length4 = median(Msf_length4, 3);
            Msb_length4 = median(Msb_length4, 3);

            %% save outputs

            outdir = [path_results '/' subject_ID '/' session_ID '/' classifier '/' folders{folder_i}];

            save([outdir '/Msf_empirical.mat'], 'Msf_empirical', '-v7.3');
            save([outdir '/Msb_empirical.mat'], 'Msb_empirical', '-v7.3');
    
            save([outdir '/Msf_length3.mat'], 'Msf_length3', '-v7.3');
            save([outdir '/Msb_length3.mat'], 'Msb_length3', '-v7.3');

            save([outdir '/Msf_length4.mat'], 'Msf_length4', '-v7.3');
            save([outdir '/Msb_length4.mat'], 'Msb_length4', '-v7.3');
        
        end
    
    end
        
end






%% pool across subjects

for folder_i = 1:length(folders)

    disp(folders{folder_i});
                    
    % empty structures for data
    Msf_l3_allsub = nan(nlags, 0);
    Msb_l3_allsub = nan(nlags, 0);
    Msf_l4_allsub = nan(nlags, 0);
    Msb_l4_allsub = nan(nlags, 0);

    null_Msf_thr_pos_l3_allsub = [];
    null_Msf_thr_neg_l3_allsub = [];
    null_Msb_thr_pos_l3_allsub = [];
    null_Msb_thr_neg_l3_allsub = [];

    null_Msf_thr_pos_l4_allsub = [];
    null_Msf_thr_neg_l4_allsub = [];
    null_Msb_thr_pos_l4_allsub = [];
    null_Msb_thr_neg_l4_allsub = [];

    % loop over trials
    for sub_i = 1:length(subjects)
    
        subject = subjects(sub_i);
    
        for ses_i = 1:length(sessions)
    
            session = sessions(ses_i);
    
            %% load data
    
            % set paths
            if subject < 10
                subject_ID = ['sub_0' num2str(subject)];
                subjectID = ['sub0' num2str(subject)];
            else
                subject_ID = ['sub_' num2str(subject)];
                subjectID = ['sub' num2str(subject)];
            end
            
            session_ID = ['sess_0' num2str(session)];
        
            %% load sequences
    
            indir = [path_data '/sequences/' subject_ID '/' session_ID];
    
            load([indir '/sequences.mat']); % sequences

            %% null distribution computed at trial level, and averaged across trials for each time lag,
            % then, for each permutation, take the
            % minimum/maximum across time lags, to get two null
            % distributions (minimum and maximum), and
            % threshold at 95% to get the statistical bound
            % see Schwartenbeck 2023

            %% length 3

            outdir = [path_results '/' subject_ID '/' session_ID '/' classifier '/' folders{folder_i}];
            
            % load sequenceness - empirical
            load([outdir '/Msf_empirical.mat']); % Msf_empirical
            load([outdir '/Msb_empirical.mat']); % Msb_empirical

            % find idx of trials of length 3 and 4
            trials_length3 = [];
            trials_length4 = [];
            for trial_i = 1:size(sequences.trialinfo, 1)

                if sequences.trialinfo(trial_i, 7) == 0

                    trials_length3 = [trials_length3 trial_i];

                else

                    trials_length4 = [trials_length4 trial_i];

                end
            end

            % subset trials - length 3
            Msf_empirical = Msf_empirical(1, trials_length3);
            Msb_empirical = Msb_empirical(1, trials_length3);

            % average sequences across trials (for each time lag)
            Msf_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
            Msb_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));

            for trial_i = 1:size(Msf_tmp,2)
        
                Msf_tmp(:,trial_i) = Msf_empirical{trial_i};
                Msb_tmp(:,trial_i) = Msb_empirical{trial_i};

            end

            Msf = median(Msf_tmp, 2);
            Msb = median(Msb_tmp, 2);

            % pool across subjects
            Msf_l3_allsub = [Msf_l3_allsub Msf];
            Msb_l3_allsub = [Msb_l3_allsub Msb];

            %%% load null distribution 

            load([outdir '/Msf_length3.mat']); % Msf_length3
            load([outdir '/Msb_length3.mat']); % Msb_length3

            % maximum (first column) and minimum (second column) across time lags for each permutation
            null_Msf = [max(Msf_length3, [],2) min(Msf_length3, [],2)];
            null_Msb = [max(Msb_length3, [],2) min(Msb_length3, [],2)];

            % threshold (P_fwe = 0.05)
            null_Msf_thr_pos = percentile(null_Msf(:,1), 97.5);
            null_Msf_thr_neg = percentile(null_Msf(:,2), 2.5);
            null_Msb_thr_pos = percentile(null_Msb(:,1), 97.5);
            null_Msb_thr_neg = percentile(null_Msb(:,2), 2.5);

            null_Msf_thr_pos_l3_allsub = [null_Msf_thr_pos_l3_allsub null_Msf_thr_pos];
            null_Msf_thr_neg_l3_allsub = [null_Msf_thr_neg_l3_allsub null_Msf_thr_neg];
            null_Msb_thr_pos_l3_allsub = [null_Msb_thr_pos_l3_allsub null_Msb_thr_pos];
            null_Msb_thr_neg_l3_allsub = [null_Msb_thr_neg_l3_allsub null_Msb_thr_neg];

            %% length 4

            % load sequenceness - empirical
            load([outdir '/Msf_empirical.mat']); % Msf_empirical
            load([outdir '/Msb_empirical.mat']); % Msb_empirical

            % find idx of trials of length 3 and 4
            trials_length3 = [];
            trials_length4 = [];
            for trial_i = 1:size(sequences.trialinfo, 1)

                if sequences.trialinfo(trial_i, 7) == 0

                    trials_length3 = [trials_length3 trial_i];

                else

                    trials_length4 = [trials_length4 trial_i];

                end
            end

            % subset trials - length 3
            Msf_empirical = Msf_empirical(1, trials_length4);
            Msb_empirical = Msb_empirical(1, trials_length4);

            % average sequences across trials (for each time lag)
            Msf_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));
            Msb_tmp = nan(size(Msf_empirical{1}, 1), size(Msf_empirical,2));

            for trial_i = 1:size(Msf_tmp,2)
        
                Msf_tmp(:,trial_i) = Msf_empirical{trial_i};
                Msb_tmp(:,trial_i) = Msb_empirical{trial_i};

            end

            Msf = median(Msf_tmp, 2);
            Msb = median(Msb_tmp, 2);

            % pool across subjects
            Msf_l4_allsub = [Msf_l4_allsub Msf];
            Msb_l4_allsub = [Msb_l4_allsub Msb];

            %%% load null distribution 

            load([outdir '/Msf_length4.mat']); % Msf_length4
            load([outdir '/Msb_length4.mat']); % Msb_length4

            % maximum (first column) and minimum (second column) across time lags for each permutation
            null_Msf = [max(Msf_length4, [],2) min(Msf_length4, [],2)];
            null_Msb = [max(Msb_length4, [],2) min(Msb_length4, [],2)];

            % threshold (P_fwe = 0.05)
            null_Msf_thr_pos = percentile(null_Msf(:,1), 97.5);
            null_Msf_thr_neg = percentile(null_Msf(:,2), 2.5);
            null_Msb_thr_pos = percentile(null_Msb(:,1), 97.5);
            null_Msb_thr_neg = percentile(null_Msb(:,2), 2.5);

            null_Msf_thr_pos_l4_allsub = [null_Msf_thr_pos_l4_allsub null_Msf_thr_pos];
            null_Msf_thr_neg_l4_allsub = [null_Msf_thr_neg_l4_allsub null_Msf_thr_neg];
            null_Msb_thr_pos_l4_allsub = [null_Msb_thr_pos_l4_allsub null_Msb_thr_pos];
            null_Msb_thr_neg_l4_allsub = [null_Msb_thr_neg_l4_allsub null_Msb_thr_neg];



        end

    end

    outdir_fig = [path_results '/group_results/' classifier '/' folders{folder_i} '/figures'];

    if ~isfolder(outdir_fig)
        mkdir(outdir_fig);
    end

    % save
    save([outdir_fig '/data_individual_plots.mat'],...
        'Msf_l3_allsub', 'Msb_l3_allsub',...
        'null_Msf_thr_pos_l3_allsub', 'null_Msb_thr_pos_l3_allsub',...
        'null_Msf_thr_neg_l3_allsub', 'null_Msb_thr_neg_l3_allsub',...
        'Msf_l4_allsub', 'Msb_l4_allsub',...
        'null_Msf_thr_pos_l4_allsub', 'null_Msb_thr_pos_l4_allsub',...
        'null_Msf_thr_neg_l4_allsub', 'null_Msb_thr_neg_l4_allsub',...
        '-v7.3');

end
            

disp('Analysis done');
