%% GERE project
% replay

folder_replay = 'sequence_concat_delay';

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
folders = {'task_stim2delay'...
           'localizer2delay'};
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
        
        outdir = [path_results '/' subject_ID '/' classifier '/' folders{folder_i}];
    
        if ~isfolder(outdir)
            mkdir(outdir);
        end
    
        X_l3 = nan(0,3);
        X_l4 = nan(0,4);
    
        X_l3_trials = zeros(1,2); % trials per session
        X_l4_trials = zeros(1,2);
        
        for ses_i = 1:length(sessions)
    
            session = sessions(ses_i);
                                                            
            session_ID = ['sess_0' num2str(session)];
        
            %% load sequences
    
            indir = [path_data '/sequences/' subject_ID '/' session_ID];
    
            load([indir '/sequences.mat']); % sequences

            %% load predictions

            indir = [path_data '/predictions/' subject_ID '/' session_ID '/' folders{folder_i} '/' classifier];
            load([indir '/predictions.mat']); % predictions

            % sigmoid function
            predictions.prob = (1 ./ (1 + exp(-predictions.dval)));
                    
            %% concatenate trials
            
            for trial_i = 1:size(predictions.prob, 1)
    
                % subset predictions for trial in the sequence
                sequence = sequences.trialinfo(trial_i, 4:7);
                sequence = sequence(sequence > 0);

                %% TDLM - sequence

                % settings TDLM
                nstates = length(sequence);
                nbins=nlags+1;

                % input - decoded state space
                X = squeeze(predictions.prob(trial_i,:,sequence)); % time by states

                % append
                if size(X,2) == 3
                    X_l3 = [X; X_l3];
                    X_l3_trials(1,ses_i) = X_l3_trials(1, ses_i) + 1;
                else
                    X_l4 = [X; X_l4];
                    X_l4_trials(1, ses_i) = X_l4_trials(1, ses_i) + 1;
                end

            end

        end

        %% TDLM
                    
        % X = state space

        % transition matrix from sequence (deleting
        % empty rows/columns, and ordering by
        % appearence in sequence, e.g. [2 5 4] is
        % encoded as [1 2 3], where [1==2, 2==5,
        % 3==4]
        nstates = 3;
        T_l3 = zeros(nstates,nstates);
        for item_i = 2:nstates
            T_l3(item_i-1, item_i) = 1;
        end

        nstates = 4;
        T_l4 = zeros(nstates,nstates);
        for item_i = 2:nstates
            T_l4(item_i-1, item_i) = 1;
        end

        % number shuffle for null distribution
        nshuf_l3 = 6; % maximum permutations of 4 states
        nshuf_l4 = 24; % maximum permutations of 3 states

        % TDLM
        if ~isfile([outdir '/Msb_l3.mat'])

            [Msf_l3, Msb_l3] = TDLM_concatenate(X_l3,T_l3,3,nlags,nshuf_l3, X_l3_trials);
            save([outdir  '/Msf_l3.mat'], 'Msf_l3', '-v7.3');
            save([outdir  '/Msb_l3.mat'], 'Msb_l3', '-v7.3');

        end

        if ~isfile([path_results '/' subject_ID '/' classifier '/' folders{folder_i}  '/Msb_l4.mat'])

            [Msf_l4, Msb_l4] = TDLM_concatenate(X_l4,T_l4,4,nlags,nshuf_l4, X_l4_trials);
            save([outdir '/Msf_l4.mat'], 'Msf_l4', '-v7.3');
            save([outdir '/Msb_l4.mat'], 'Msb_l4', '-v7.3');

        end

    end
                        
end

disp('Analysis done');

