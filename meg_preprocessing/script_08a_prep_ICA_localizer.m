%% preprocessing MEG data for WM task [GERE project]

% path to outputs
path_results = ['/path_to_local/results/preprocessing/' folder];
path_preprocessing = '/path_to_local/results/preprocessing/triggers';

%% loop over trials

zthr = 25;

layout = [];

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
                        
        meg_data = [path_preprocessing '/' subject_ID '/sess_0' num2str(session)...
            '/GERE_' subjectID '_sess' num2str(session) '_localizer.con'];
        
        if isfile(meg_data)
        
            %% compute layout

            if isempty(layout)
            
                % prepare 2D layout
                hdr            = ft_read_header(meg_data);
                cfg            = [];
                cfg.grad       = hdr.grad;
                layout         = ft_prepare_layout(cfg);

            end
            
            %% attenuate artifacts with ICA
            % https://www.fieldtriptoolbox.org/tutorial/ica_artifact_cleaning/
            % https://www.fieldtriptoolbox.org/example/ica_ecg/
            % https://www.fieldtriptoolbox.org/example/ica_eog/
            
            % localizer
            data_locked_file = [path_results '/' subject_ID '/' session_ID ...
                '/GERE_localizer_validtrialszthr' num2str(zthr) '.mat'];
    
            if isfile(data_locked_file)
    
                if ~isfile([path_results '/' subject_ID '/sess_0' num2str(session)...
                    '/GERE_localizer_ICA_validtrialszthr'...
                            num2str(zthr) '.mat'])
    
                    load(data_locked_file); % data_validtrials
                            
                    % [flux] Then we apply a 1 - 40 Hz  bandpass filter. 
                    % The 1 Hz highpass is important for removing slow drifts 
                    % which otherwise would make the ICA decomposition less efficient. 
                    % We also apply padding around our trialdata 
                    % (3 sec on each side of each trial) to reduce filtering artifacts. 
                    
                    cfg                             = [];
                    cfg.hpfreq                      = 1;
                    cfg.padding                     = 10.5;
                    cfg.padtype                     = 'zero';
                    data_no_artifacts_downsampled_filt = ft_preprocessing(cfg, data_validtrials);
                    
                    % ICA Step 1 - ICA Decomposition
                            
                    % rank of data (first trial) - if dimensionality reduction technique has been applied
                    data_rank = rank(data_no_artifacts_downsampled_filt.trial{1}*data_no_artifacts_downsampled_filt.trial{1}');
                    
                    % run ica
                    cfg                             = [];  
                    cfg.method                      = 'fastica';  
                    cfg.numcomponent                = data_rank;
                    data_comp = ft_componentanalysis(cfg, data_no_artifacts_downsampled_filt);   
                    
                    % save
                    save([path_results '/' subject_ID '/sess_0' num2str(session)...
                    '/GERE_localizer_ICA_validtrialszthr'...
                            num2str(zthr) '.mat'], 'data_comp');
    
                end
    
            end


           

        else

            disp('No MEG data');

        end

    end

end
