%% source reconstruction fieldtrip
% https://www.youtube.com/watch?v=Ez72OFjSABs&ab_channel=fieldtriptoolbox

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

% settings
ft_defaults;

% n = 15 - exclude 2 subjects without MRI (14 and 41)
subjects = [5 7 18 23 25 31 34 37 40 45 47 53 61 201 202];
sessions = 1:2;

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

        if ~isfolder([path_results '/' subject_ID '/' session_ID])
            mkdir([path_results '/' subject_ID '/' session_ID]);
        end

        disp([subject_ID ' ' session_ID]);
            
        % load MRI data
        mri = ft_read_mri([path_mri '/' subject_ID '/t1.nii.gz']);
        mri.coordsys = 'neuromag';

        %% 1. Compute the forward model
        % X(t) = h(r) * s(r, t)
        % X(t) - matrix of timecourses per sensor
        % s(r, t) - matrix of timecourses of r sources 
        % h(r) - leadmatrix, matrix of contribution of each sensor on r sources, 
        % see 1.1 to 1.4 below
                    
        % loop over runs (run 0 is localizer, while runs 1 to 3 are the WM task)
        for run_i = 0:3

            if run_i == 0

                filename = 'localizer';

            else

                filename = ['task' num2str(run_i)];

            end

            disp([subject_ID ' ' session_ID ' run ' num2str(run_i)]);

            if ~isfile([path_results '/' subject_ID '/' session_ID '/figure_alignment_' filename '.fig'])
                                    
                %% 1.1. sensor location (where is the brain with respect to the sensors?)
                % Align sensors to headhape 
    
                %%% points/basic files
    
                cd([path_head '/' subject_ID '/' session_ID]);
                system('cp *basic.txt basic.txt');
                system('cp *points.txt points.txt');
    
                cd([path_results '/' subject_ID '/' session_ID]);
    
                point_file = [path_head '/' subject_ID '/' session_ID '/points.txt'];
                basic_file = [path_head '/' subject_ID '/' session_ID '/basic.txt'];
    
                convertLaserFiles_automatic_aniol(point_file, basic_file, 'headscan');
                close all;
    
                % remove intermediary
                system(['rm ' path_head '/' subject_ID '/' session_ID '/basic.txt']);
                system(['rm ' path_head '/' subject_ID '/' session_ID '/points.txt']);
    
                % read .hsp .elp (fiducials in headspace)- output 'shape.mat'
                parsePolhemus('headscan');
    
                % load shape.mat = headshape and fiducial in head coordinate system
                x = load('shape'); % in mm, headshape and fiducial data from the polhemus
                shape = x.shape;

                % remove intermediary file
                system('rm shape.mat');

                % apply rotation of headshape to mri (rotate 180 degrees in
                % vertical axis), that will make further registrations
                % easier
                rotation_matrix   = [-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1];
                cfg               = [];
                cfg.m             = rotation_matrix;
                shape       = ft_transform_geometry(cfg.m, shape);
    
                %%% markers
    
                cd([path_meg '/' subject_ID '/' session_ID]);

                % system('cp *-1.mrk localizer.mrk');
                try
                    system(['cp *-' num2str(run_i+1) '.mrk task' num2str(run_i) '.mrk']);
    
                    % markers in dewar coordinate system
                    marker = [path_meg '/' subject_ID '/' session_ID '/task' num2str(run_i) '.mrk'];
        
                    %%% meg
                    
                    cd([path_meg '/' subject_ID '/' session_ID]);
                    system(['cp *0' num2str(run_i+1) '_analysis_01.con meg_data.con']);
            
                    sensor                      = ft_read_sens([path_meg '/' subject_ID '/' session_ID '/meg_data.con'],...
                        'senstype', 'meg','fileformat', 'yokogawa_con'); % sensor positions in head coordinate system, same than hdr.grad
                    sensor                      = ft_convert_units(sensor,'mm');
        
                    % get marker file, read t in and convert to cm
                    mrk                         = ft_read_headshape(marker, 'format', 'yokogawa_mrk');
                    mrk                         = ft_convert_units(mrk, 'mm');

                % take the post-run markers if the pre-run markers do not
                % work properly
                catch

                    system(['cp *-' num2str(run_i+2) '.mrk task' num2str(run_i) '.mrk']);
    
                    % markers in dewar coordinate system
                    marker = [path_meg '/' subject_ID '/' session_ID '/task' num2str(run_i) '.mrk'];
        
                    %%% meg
                    
                    cd([path_meg '/' subject_ID '/' session_ID]);
                    system(['cp *0' num2str(run_i+1) '_analysis_01.con meg_data.con']);
            
                    sensor                      = ft_read_sens([path_meg '/' subject_ID '/' session_ID '/meg_data.con'],...
                        'senstype', 'meg','fileformat', 'yokogawa_con'); % sensor positions in head coordinate system, same than hdr.grad
                    sensor                      = ft_convert_units(sensor,'mm');
        
                    % get marker file, read t in and convert to cm
                    mrk                         = ft_read_headshape(marker, 'format', 'yokogawa_mrk');
                    mrk                         = ft_convert_units(mrk, 'mm');

                end
    
                % match between labels in shape.fid.label and mrk.fid.label
                shape.fid.label;
                %     {'NASION'  }
                %     {'LPA'     }
                %     {'RPA'     }
                %     {'LPAred'  }
                %     {'RPAyel'  }
                %     {'PFblue'  }
                %     {'LPFwh'   }
                %     {'RPFblack'}
                mrk.fid.label;
                %     {'nas'    } = CPF
                %     {'lpa'    } = LPA
                %     {'rpa'    } = RPA
                %     {'Marker4'} = LPF
                %     {'Marker5'} = RPF
                markers_dewarcoordsys       = mrk.fid.pos; % five fiducials
                markers_headcoordsys        = shape.fid.pnt([1 2 3 7 8],:); % select the fiducials that match THIS IS DIFFERENT FROM ORIGINAL SCRIPT)
    
                % get transformation matrix for markers (dewar coordsys) to points on head shape (head coordsys) 
                [R,T,Yf,Err]                = rot3dfit(markers_dewarcoordsys, markers_headcoordsys); % calc rotation transform
                meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix
    
                % use transformation matrix to align sensors location to headshape. 
                sensor_aln                  = ft_transform_geometry(meg2head_transm, sensor);
                        
                save([path_results '/' subject_ID '/' session_ID '/sensor_aln_' filename '.mat'], 'sensor_aln');
    
                % figure; hold on
                % ft_plot_sens(sensor_aln);
                % ft_plot_headshape(shape);

                % remove intermediary files
                system(['rm ' path_meg '/' subject_ID '/' session_ID '/task*.mrk']);
                system(['rm ' path_meg '/' subject_ID '/' session_ID '/meg_data.con']);
    
                cd(path_scripts);
                    
                %% 1.2. head model or volume conduction model (with MRI data)
    
                %%% Coregistration of anatomical MRI image with MEG sensor array
                % https://www.fieldtriptoolbox.org/workshop/paris2019/handson_anatomy/
    
                % remove other fiducials than nas, lpa and rpa from shape,
                % and rename as headshape
                headshape = shape;
                headshape.fid.pos = headshape.fid.pnt(1:3,:);
                headshape.fid.label = headshape.fid.label(1:3,:);
                headshape.pos = headshape.pnt;
                headshape.fid = rmfield(headshape.fid, 'pnt');
                headshape = rmfield(headshape, 'pnt');
                headshape = ft_convert_units(headshape,'mm');
    
                % %%%%% plot sensors and shape in head coordinate system
                % figure; hold on;
                % ft_plot_sens(sensor_aln);
                % ft_plot_headshape(headshape);

                if ~isfile([path_results '/' subject_ID '/' session_ID '/t1_aln_' filename '.nii.gz'])
                                                
                    % automated alignment of MRI to headshape
                    cfg = [];
                    cfg.method = 'headshape';
                    cfg.headshape.headshape = headshape;
                    cfg.headshape.interactive = 'no';
                    t1_aln = ft_volumerealign(cfg, mri);
                                
                    % manual alignment of MRI to headshape (when manual alignment done, press QUIT)
                    cfg = [];
                    cfg.method = 'headshape';
                    cfg.headshape.headshape = headshape;
                    cfg.headshape.icp = 'no';
                    t1_aln = ft_volumerealign(cfg, t1_aln);
        
                    % write volume
                    cfg = [];
                    cfg.parameter = 'anatomy';
                    cfg.filename = [path_results '/' subject_ID '/' session_ID '/t1_aln_' filename '.nii.gz'];
                    cfg.filetype = 'nifti';  
                    ft_volumewrite(cfg, t1_aln)
        
                    % save
                    save([path_results '/' subject_ID '/' session_ID '/t1_aln_' filename '.mat'], 't1_aln');

                else

                    load([path_results '/' subject_ID '/' session_ID '/t1_aln_' filename '.mat']);

                end
                
                % % check coregistration
                % cfg = [];
                % cfg.method = 'headshape';
                % cfg.headshape.headshape = headshape_rotated;
                % check_alignment = ft_volumerealign(cfg, t1_aln);
    
                % segment the original mri
                cfg                = [];
                cfg.output         = 'brain';
                seg_t1_aln = ft_volumesegment(cfg, t1_aln);
                
                % construct the volume conductor model (i.e. head model) for each subject
                cfg        = [];
                cfg.method = 'singleshell';
                headmodel  = ft_prepare_headmodel(cfg, seg_t1_aln);
    
                % save outputs
                save([path_results '/' subject_ID '/' session_ID '/headmodel_' filename '.mat'], 'headmodel');
    
                %% 1.3. source model
                % which locations do you want to 'scan'?
                % source model describes a set of positions (and possibly orientations) of 
                % equivalent current dipoles that are taken into consideration when doing the source reconstruction
                % https://www.fieldtriptoolbox.org/tutorial/sourcemodel/
                
                % create template grid
                cfg = [];
                cfg.method = 'basedonresolution';
                cfg.resolution  = 0.8;
                cfg.unit   = 'cm';
                cfg.inwardshift = -0.5;
                cfg.headmodel   = headmodel;
                sourcemodel   = ft_prepare_sourcemodel(cfg);
                sourcemodel  = ft_convert_units(sourcemodel,'mm');
    
                %% 1.4. leadfield matrix
                % leadfield [h(r)] coefficients that weights the timecourse of each sensor to compute the
                % timecourse of a source
                % https://www.fieldtriptoolbox.org/example/make_leadfields_using_different_headmodels/
    
                % prepare the leadfield for the single-shell model
                cfg                 = [];
                cfg.grad            = sensor_aln;
                cfg.headmodel       = headmodel;
                cfg.sourcemodel     = sourcemodel;
                cfg.unit            = 'mm';
                forwardmodel  = ft_prepare_leadfield(cfg);
    
                save([path_results '/' subject_ID '/' session_ID '/forwardmodel_' filename '.mat'], 'forwardmodel');
    
                cd(path_scripts);
    
                % make a figure of the single subject head model, source
                % model and sensors
                figure; hold on;
                ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
                ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
                ft_plot_sens(sensor_aln);
                ft_plot_headshape(headshape);

                % Save the plot as a PNG file
                saveas(gcf, [path_results '/' subject_ID '/' session_ID '/figure_alignment_' filename '.png']);
                saveas(gcf, [path_results '/' subject_ID '/' session_ID '/figure_alignment_' filename '.fig']);

                close all

            end

        end

    end
    
end

disp('Analysis done');
