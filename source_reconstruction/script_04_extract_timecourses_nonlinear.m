%% source reconstruction fieldtrip
% https://www.youtube.com/watch?v=Ez72OFjSABs&ab_channel=fieldtriptoolbox

% extract timecourse of virtual channels (sources)
% https://www.fieldtriptoolbox.org/workshop/paris2019/handson_sourceanalysis/

clear variables; close all;

path_utilities = '/path_to_local/scripts/source_reconstruction/utilities';
addpath(path_utilities);
addpath('/path_to_software/fieldtrip-lite-2020083');
path_mri = '/path_to_local/data/MRI';
path_head = '/path_to_local/data/head_shape';
path_meg = '/path_to_local/data/MEG';
path_preproc_data = '/path_to_local/results/preprocessing/M3';
path_results = '/path_to_local/results/source_reconstruction';
path_scripts = '/path_to_local/scripts/source_reconstruction';
path_atlas = '/path_to_software/atlas/Schaefer_2018/MNI';

atlas_parcels = [200];

% settings
ft_defaults;

% n = 15 - exclude 2 subjects without MRI (14 and 41)
subjects = [5 7 18 23 25 31 34 37 40 45 47 53 61 201 202];
sessions = 1:2;

%% localizer

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

        disp(['localizer ' subject_ID ' ' session_ID]);

        if ~isfile([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(end)) '_localizer_zthr25.mat'])

            % load forward and head models
            load([path_results '/' subject_ID '/' session_ID '/headmodel_localizer.mat']); % headmodel
            load([path_results '/' subject_ID '/' session_ID '/forwardmodel_localizer.mat']); % forwardmodel
    
            % load preprocessed data
            load([path_preproc_data '/' subject_ID '/' session_ID '/GERE_localizer_validtrialszthr25_postICA.mat']); % data
    
            % compute data covariance matrix
            cfg                   = [];
            cfg.covariance        = 'yes';
            cfg.covariancewindow  = 'all';
            data_tlock = ft_timelockanalysis(cfg, data);
    
            % beamforming
            cfg = [];
            cfg.method              = 'lcmv';
            cfg.sourcemodel         = forwardmodel;
            cfg.headmodel           = headmodel;
            cfg.lcmv.keepfilter     = 'yes';
            cfg.lcmv.lambda         = '5%'; % 
            cfg.lcmv.projectmom     = 'yes'; % project dipole time series in direction of maximal power (see below)
            source_data = ft_sourceanalysis(cfg, data_tlock);
    
            for number_parcels_i = 1:length(atlas_parcels)
    
                number_parcels = atlas_parcels(number_parcels_i);
    
                if ~isfile([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(number_parcels) '_localizer_zthr25.mat'])
    
                    %% find parcel for each source
            
                    % import atlas
                    atlas_individual = [path_results '/' subject_ID '/' session_ID '/atlas_nonlinear_localizer_Schaefer2018_' num2str(number_parcels) 'Parcels_7Networks_order_FSLMNI152_1mm.nii.gz'];
                    atlas = ft_read_atlas(atlas_individual);
            
                    number_sources = size(forwardmodel.pos,1);
                    parcellation = zeros(number_sources,1);
            
                    for source_i = 1:number_sources
            
                        if forwardmodel.inside(source_i) == 1
            
                            atlas_coordinates = forwardmodel.pos(source_i,:);
                            atlas_coordinates = [(atlas_coordinates -1) 1];
                            voxel_coordinates = round((atlas.transform) \ (atlas_coordinates'));
            
                            try
            
                                parcellation(source_i,1) = atlas.parcellation(voxel_coordinates(1), voxel_coordinates(2), voxel_coordinates(3));
            
                            catch
                            end
            
                        end
            
                    end
            
                    %% extract timeseries in sources
            
                    data_atlas = data;
                    data_atlas.sources_per_parcel = data_atlas.trial;
            
                    for trial_i = 1:size(data.trial,2)

                        sources_per_parcel = [];
                
                        number_timepoints = size(data.trial{trial_i},2);
                        source_timecourses = [];
                        source_timecourses.timecourses = nan(number_sources, number_timepoints);
                        source_timecourses.pos = forwardmodel.pos;
                        source_timecourses.inside = forwardmodel.inside;
                        source_timecourses.atlas = nan(number_parcels, number_timepoints);
            
                        for source_i = 1:number_sources
            
                            if forwardmodel.inside(source_i) == 1
            
                                spatial_filter = source_data.avg.filter{source_i};
                    
                                if ~isempty(spatial_filter)
                    
                                    source_timecourses.timecourses(source_i,:) = abs(spatial_filter * data.trial{trial_i});
            
                                    % take absolute value to deal with sign ambiguity
                                    % https://mailman.science.ru.nl/pipermail/fieldtrip/2011-August/017056.html
            
                                end
            
                            end
            
                        end
            
                        %% eigenvariate of sources within each parcel
                        
                        % temporal standard deviation for all sources
                        temporalSTD = max(std(source_timecourses.timecourses, [], 2), eps);                
            
                        % eigenvariate of sources within each atlas parcel
                        for parcel_i = 1:number_parcels
            
                            [idx, ~] = find(parcellation == parcel_i);
                            output = source_timecourses.timecourses(idx,:);
            
                            % remove NaN
                            idx_no_nan = ~isnan(output);
                            idx_no_nan = idx_no_nan(:,1);
                            output = output(idx_no_nan,:);
                                        
                            %%% PCA - roinets toolbox
            
                            % if parcellation is empty
                            if size(output,1)==0

                                sources_per_parcel = [sources_per_parcel 0];
            
                                % store
                                source_timecourses.atlas(parcel_i,:) = mean(source_timecourses.timecourses, 'omitnan');
            
                            % if parcellation has only one source
                            elseif size(output,1)==1

                                sources_per_parcel = [sources_per_parcel 1];
            
                                source_timecourses.atlas(parcel_i,:) = output;
            
                            else

                                sources_per_parcel = [sources_per_parcel size(output,1)];
                                
                                % SVD
                                [U, S, V]  = ROInets.fast_svds(output, 1);
                                PCAscores  = S * V';
            
                                % restore sign and scaling of parcel time-series
                                % U indicates the weight with which each voxel in the
                                % parcel contributes to the 1st PC
                                TSsign          = sign(mean(U));
                                relVoxelWeights = abs(U) ./ sum(abs(U)); % normalise the linear combination
            
                                % weight the temporal STDs from the ROI by the proportion used in 1st PC
                                tmp_STD = temporalSTD(idx);
                                tmp_STD = tmp_STD(idx_no_nan);
                                TSscale         = dot(relVoxelWeights, tmp_STD); 
                                nodeTS          = TSsign .*                               ...
                                                  (TSscale / max(std(PCAscores), eps)) .* ... 
                                                  PCAscores;
            
                                % store
                                source_timecourses.atlas(parcel_i,:) = nodeTS;
            
                            end
            
                        end
            
                        % store
                        data_atlas.trial{trial_i} = source_timecourses.atlas;
                        data_atlas.sources_per_parcel{trial_i} = sources_per_parcel;

                    end
                    
                    save([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(number_parcels) '_localizer_zthr25.mat'], 'data_atlas', '-v7.3');
    
                end
    
            end

        end

    end

end

%% task

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

        if ~isfile([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(atlas_parcels(end)) '_task_zthr25.mat'])

            %% data task
    
            % load data_task
            load([path_preproc_data '/' subject_ID '/' session_ID...
                    '/GERE_task_zthr25_wholetrial_baseline1sec_postICA.mat']);
    
            data_task = data;
    
            clear data
    
            % set runs in data_task
            data_task.trialinfo(:,end+1) = 0;
            data_task.trialinfo(1,end) = 1;
            run_count = 1;
            for trial_i = 2:size(data_task.trial,2)
                if data_task.trialinfo(trial_i,1) > data_task.trialinfo(trial_i-1,1)
                    data_task.trialinfo(trial_i,end) = run_count;
                else
                    run_count = run_count + 1;
                    data_task.trialinfo(trial_i,end) = run_count;
                end
            end

            for run_i = 1:3
    
                filename = ['task' num2str(run_i)];
    
                disp([subject_ID ' ' session_ID ' task ' num2str(run_i)]);
    
                % load forward and head models
                load([path_results '/' subject_ID '/' session_ID '/headmodel_' filename '.mat']); % headmodel
                load([path_results '/' subject_ID '/' session_ID '/forwardmodel_' filename '.mat']); % forwardmodel
    
                % load preprocessed data
                load([path_preproc_data '/' subject_ID '/' session_ID '/GERE_task_zthr25_wholetrial_baseline1sec_postICA.mat']); % data
    
                % compute data covariance matrix
                cfg                   = [];
                cfg.covariance        = 'yes';
                cfg.covariancewindow  = 'all';
                data_tlock = ft_timelockanalysis(cfg, data);
    
                % beamforming
                cfg = [];
                cfg.method              = 'lcmv';
                cfg.sourcemodel         = forwardmodel;
                cfg.headmodel           = headmodel;
                cfg.lcmv.keepfilter     = 'yes';
                cfg.lcmv.lambda         = '5%'; % 
                cfg.lcmv.projectmom     = 'yes'; % project dipole time series in direction of maximal power (see below)
                source_data = ft_sourceanalysis(cfg, data_tlock);
    
                eval(['source_data_' filename ' = source_data;'])
    
            end
    
            for number_parcels_i = 1:length(atlas_parcels)
    
                number_parcels = atlas_parcels(number_parcels_i);
                
                if ~isfile([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(number_parcels) '_task_zthr25.mat'])
    
                    %% find parcel for each source
            
                    % import atlas
                    atlas_individual = [path_results '/' subject_ID '/' session_ID '/atlas_nonlinear_localizer_Schaefer2018_' num2str(number_parcels) 'Parcels_7Networks_order_FSLMNI152_1mm.nii.gz'];
                    atlas = ft_read_atlas(atlas_individual);
    
                    for run_i = 1:3

                        filename = ['task' num2str(run_i)];
                        eval(['source_data = source_data_' filename ';'])
            
                        number_sources = size(source_data.pos,1);
                        parcellation = zeros(number_sources,1);
                
                        for source_i = 1:number_sources
                
                            if source_data.inside(source_i) == 1
                
                                atlas_coordinates = source_data.pos(source_i,:);
                                atlas_coordinates = [(atlas_coordinates -1) 1];
                                voxel_coordinates = round((atlas.transform) \ (atlas_coordinates'));
                
                                try
                
                                    parcellation(source_i,1) = atlas.parcellation(voxel_coordinates(1), voxel_coordinates(2), voxel_coordinates(3));
                
                                catch
                                end
                
                            end
                
                        end
    
                        eval(['parcellation_' filename ' = parcellation;'])
            
                    end
    
                    %% extract timeseries in sources
            
                    data_atlas = data_task;
                    data_atlas.sources_per_parcel = data_atlas.trial;
            
                    for trial_i = 1:size(data_task.trial,2)

                        sources_per_parcel = [];
                
                        eval(['source_data = source_data_task' num2str(data_task.trialinfo(trial_i,end)) ';']);
                        eval(['parcellation = parcellation_task' num2str(data_task.trialinfo(trial_i,end)) ';']);

                        number_sources = size(source_data.pos,1);
            
                        number_timepoints = size(data_task.trial{trial_i},2);
                        source_timecourses = [];
                        source_timecourses.timecourses = nan(number_sources, number_timepoints);
                        source_timecourses.pos = source_data.pos;
                        source_timecourses.inside = source_data.inside;
                        source_timecourses.atlas = nan(number_parcels, number_timepoints);
            
                        for source_i = 1:number_sources
            
                            if source_data.inside(source_i) == 1
            
                                spatial_filter = source_data.avg.filter{source_i};
        
                                if ~isempty(spatial_filter)
                    
                                    source_timecourses.timecourses(source_i,:) = abs(spatial_filter * data_task.trial{trial_i});
            
                                    % take absolute value to deal with sign ambiguity
                                    % https://mailman.science.ru.nl/pipermail/fieldtrip/2011-August/017056.html
            
                                end
        
                            end
            
                        end
            
                        %% eigenvariate of sources within each parcel
                        
                        % temporal standard deviation for all sources
                        temporalSTD = max(std(source_timecourses.timecourses, [], 2), eps);                
            
                        % eigenvariate of sources within each atlas parcel
                        for parcel_i = 1:number_parcels
            
                            [idx, ~] = find(parcellation == parcel_i);
                            output = source_timecourses.timecourses(idx,:);
            
                            % remove NaN
                            idx_no_nan = ~isnan(output);
                            idx_no_nan = idx_no_nan(:,1);
                            output = output(idx_no_nan,:);
                                        
                            %%% PCA - roinets toolbox
            
                            % if parcellation is empty
                            if size(output,1)==0

                                sources_per_parcel = [sources_per_parcel 0];
            
                                % store
                                source_timecourses.atlas(parcel_i,:) = mean(source_timecourses.timecourses, 'omitnan');
            
                            % if parcellation has only one source
                            elseif size(output,1)==1

                                sources_per_parcel = [sources_per_parcel 1];
            
                                source_timecourses.atlas(parcel_i,:) = output;
            
                            else

                                sources_per_parcel = [sources_per_parcel size(output,1)];
                                
                                % SVD
                                [U, S, V]  = ROInets.fast_svds(output, 1);
                                PCAscores  = S * V';
            
                                % restore sign and scaling of parcel time-series
                                % U indicates the weight with which each voxel in the
                                % parcel contributes to the 1st PC
                                TSsign          = sign(mean(U));
                                relVoxelWeights = abs(U) ./ sum(abs(U)); % normalise the linear combination
            
                                % weight the temporal STDs from the ROI by the proportion used in 1st PC
                                tmp_STD = temporalSTD(idx);
                                tmp_STD = tmp_STD(idx_no_nan);
                                TSscale         = dot(relVoxelWeights, tmp_STD); 
                                nodeTS          = TSsign .*                               ...
                                                  (TSscale / max(std(PCAscores), eps)) .* ... 
                                                  PCAscores;
            
                                % store
                                source_timecourses.atlas(parcel_i,:) = nodeTS;
            
                            end
            
                        end
            
                        % store
                        data_atlas.trial{trial_i} = source_timecourses.atlas;
                        data_atlas.sources_per_parcel{trial_i} = sources_per_parcel;
            
                    end
                    
                    save([path_results '/' subject_ID '/' session_ID '/data_nonlinear_atlas' num2str(number_parcels) '_task_zthr25.mat'], 'data_atlas', '-v7.3');
    
                end
    
            end

        end

    end

end


disp('Analysis done');
