%% preprocessing MEG data for WM task [GERE project]

% path to data
path_results = ['/path_to_local/results/preprocessing/' folder];
path_preprocessing = '/path_to_local/results/preprocessing/triggers';

%% loop over trials

zthr = 25;

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

        % localizer
                
        if ~isfile([path_results '/' subject_ID '/sess_0' num2str(session)...
                '/GERE_localizer_validtrialszthr' num2str(zthr) '_postICA.mat'])
                
            %% compute layout
    
            meg_data = [path_preprocessing '/' subject_ID '/' session_ID...
                '/GERE_' subjectID '_sess' num2str(session) '_localizer.con'];

            if isfile(meg_data)
            
                % prepare 2D layout
                hdr            = ft_read_header(meg_data);
                cfg            = [];
                cfg.grad       = hdr.grad;
                layout         = ft_prepare_layout(cfg);
                                
                %% ICA step 2 - Identifying components
        
                % load ICA
                load([path_results '/' subject_ID '/sess_0' num2str(session)...
                        '/GERE_localizer_ICA_validtrialszthr' num2str(zthr) '.mat']); % data_comp
    
                if ~isfield(data_comp,'rejected_components')
    
                    components2plot = [1:20; 21:40];
            
                    % loop over bins of 20 components for plotting
                    for comp_i = 1:size(components2plot, 1)
                        
                        % topographic plot
                        cfg           = [];
                        cfg.component = components2plot(comp_i, :);
                        cfg.layout    = layout;
                        cfg.comment   = 'no';
                        ft_topoplotIC(cfg, data_comp);
                        saveas(gcf,[path_results '/' subject_ID '/sess_0' num2str(session) '/GERE_localizer_ICA_' num2str(components2plot(comp_i, 1)) 'to' num2str(components2plot(comp_i, end)) 'ICs.fig']);
                        close all;
                        openfig([path_results '/' subject_ID '/sess_0' num2str(session) '/GERE_localizer_ICA_' num2str(components2plot(comp_i, 1)) 'to' num2str(components2plot(comp_i, end)) 'ICs.fig']);
            
                        % display response
                        if comp_i == 1
                            rejected_components = [1000 ];
                        end
        
                        while rejected_components(end) ~= 0
                            rejected_components = [ rejected_components input(sprintf('\n Components to reject (enter 0 when done): '))];
                        end
        
                        % reset rejected components
                        rejected_components = rejected_components(1:(end-1));
                
                        close all;
        
                    end
            
                    % reset rejected components (after loop)
                    rejected_components = rejected_components(2:end);
                    rejected_components = unique(rejected_components);
            
                    % save
                    data_comp.rejected_components = rejected_components;
                    save([path_results '/' subject_ID '/sess_0' num2str(session) '/GERE_localizer_ICA.mat'], 'data_comp'); % data_comp
                
                end
                
                %% ICA Step 3 - Rejecting components
    
                % load data
                data_locked_file = [path_results '/' subject_ID '/' session_ID ...
                    '/GERE_localizer_validtrialszthr' num2str(zthr) '.mat'];
                load(data_locked_file); % data_validtrials
                
                if ~isempty(rejected_components)
                        
                    % reject components of the original data in its sampling rate
                    cfg = [];  
                    cfg.component = data_comp.rejected_components;
                    data = ft_rejectcomponent(cfg, data_comp, data_validtrials);
            
                else
    
                    data = data_validtrials;
    
                end
    
                % save
                save([path_results '/' subject_ID '/sess_0' num2str(session)...
                    '/GERE_localizer_validtrialszthr' num2str(zthr) '_postICA.mat'], 'data');
    
            end

        end



    end

end
