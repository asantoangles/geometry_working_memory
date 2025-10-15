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
path_atlas = '/path_to_software/atlas/Schaefer_2018/MNI';

atlas_parcels = 200;

% settings
ft_defaults;

% n = 15 (zthr 10) - exclude 2 subjects without MRI (14 and 41)
subjects = [5 7 18 23 25 31 34 37 40 45 47 53 61 201 202];
sessions = 1:2;


%% linear + nonlinear registration of structural to mni

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

        for run_i = 0:3
    
            if run_i == 0
    
                filename = 'localizer';
    
            else
    
                filename = ['task' num2str(run_i)];
    
            end
        
            for number_parcels_i = 1:length(atlas_parcels)

                number_parcels = atlas_parcels(number_parcels_i);

                disp(' '); disp(' '); disp(' '); 
                disp([subject_ID ' ' session_ID]);                
                disp(datetime); disp(filename); disp(['atlas ' num2str(number_parcels)]);
        
                atlas = [path_atlas '/Schaefer2018_' num2str(number_parcels) 'Parcels_7Networks_order_FSLMNI152_1mm.nii.gz'];
                t1 = [path_results '/' subject_ID '/' session_ID '/t1_aln_' filename '.nii.gz'];
                mni = '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz';
                structural2mni = [path_results '/' subject_ID '/' session_ID '/structural2mni_' filename '.mat'];
                mni2structural = [path_results '/' subject_ID '/' session_ID '/mni2structural_' filename '.mat'];
                atlas_individual_nonlinear = [path_results '/' subject_ID '/' session_ID '/atlas_nonlinear_' filename '_Schaefer2018_' num2str(number_parcels) 'Parcels_7Networks_order_FSLMNI152_1mm'];
                warp_structural2mni = [path_results '/' subject_ID '/' session_ID '/warp_structural2mni_' filename];
                warp_mni2structural = [path_results '/' subject_ID '/' session_ID '/warp_mni2structural_' filename];

                disp('linear registration');

                if ~isfile(structural2mni)

                    if strcmp(subject_ID, 'sub_40')

                        % registration structural to mni (linear) with
                        % brain extracted image
                        t1_mni = [path_results '/' subject_ID '/' session_ID '/t1_aln_' filename '_mni.nii.gz'];
                        mni_brain = '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz';
                        t1_brain = [path_results '/' subject_ID '/' session_ID '/t1_aln_' filename '_brain.nii.gz'];
                        system(['bet ' t1 ' ' t1_brain]);
                        system(['flirt -in ' t1_brain ' -ref ' mni_brain ' -omat ' structural2mni ' -o ' t1_mni]);
                        system(['rm ' t1_brain]);
                        system(['rm ' t1_mni]);

                    else

                        % registration structural to mni (linear)
                        system(['flirt -in ' t1 ' -ref ' mni ' -omat ' structural2mni]);

                    end

                end

                disp('nonlinear registration');

                if ~isfile([warp_structural2mni '.nii.gz'])
                        
                    % nonlinear transformation mni to structural
                    system(['fnirt --in=' t1 ' --ref=' mni ' --aff=' structural2mni ' --cout=' warp_structural2mni]);
                    
                end

                disp('apply to atlas');

                if ~isfile([atlas_individual_nonlinear '.nii.gz'])
        
                    % compute inverse matrices (from mni to structural), to apply later on the atlas
                    system(['convert_xfm -omat ' structural2mni ' -inverse ' mni2structural]);
                    system(['invwarp --ref=' t1 ' --warp=' warp_structural2mni ' --out=' warp_mni2structural]);
    
                    % atlas to individual
                    system(['applywarp --ref=' t1 ' --in=' atlas ' --warp=' warp_mni2structural ' --out=' atlas_individual_nonlinear]);
    
                end

            end
    
        end

    end

end

disp('Analysis done');
