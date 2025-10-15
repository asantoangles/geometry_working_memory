function X = compute_X_lm_baseline_corrected(data_task, data_baseline, sequence_length, number_locations, shifted_locations)

% compute_X_lm_baseline_corrected computes a baseline-corrected population 
% activity matrix (X) using ridge regression across trials.
%
% This function estimates the contribution of each location and sequence
% rank to neural activity after baseline correction. It builds a design
% matrix describing the spatial and temporal structure of the task 
% (locations × ranks) and fits a linear model to predict parcel activity 
% from these design variables using ridge regression, separately for each
% neural feature (meg sensor, source parcel, neural unit).
%
% The resulting matrix X contains the regression coefficients for each 
% feature (cortical parcel in source reconstructed meg data), 
% reflecting baseline-corrected neural responses associated with 
% specific locations and sequence ranks.
%
% INPUTS:
%   data_task        - FieldTrip-like data structure containing task-period 
%                      trials (fields: trial, time, trialinfo, label)
%   data_baseline    - FieldTrip-like data structure containing baseline 
%                      trials matched to data_task
%   sequence_length  - Length of the stimulus sequence (e.g., 3, 4 or 34)
%   number_locations - Number of possible spatial locations (typically 4 or 8)
%   shifted_locations (optional) - Logical flag (0 or 1). If set to 1, 
%                      location pooling is shifted circularly across trials. 
%                      Default = 0.
%
% OUTPUT:
%   X - [number_locations * ranks] × [number of parcels] matrix of 
%       ridge regression coefficients, representing baseline-corrected
%       neural activity across conditions and ranks.
%
% DETAILS:
% - For sequence_length = 34, sequence lengths of 3 and 4 are included, and the number of ranks is set to 4.
% - Baseline correction is applied at the trial level as the difference 
%   between the mean task activity and the mean baseline activity 
%   per parcel and trial.
% - Ridge regression is performed separately for each parcel using λ = 1.
%
% NOTE:
%   The resulting X matrix is not column-demeaned.
%
% REFERENCES:
%   Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. 
%   Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. 
%   bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

if nargin < 5
    shifted_locations = 0;
end


%% subset data_task (trials depending on sequence length)
if sequence_length == 34

    ranks = 4;

else

    ranks = sequence_length;

    idx_length = [];

    % subset data_task
    for trial_i = 1:size(data_task.trial, 2)
    
        sequence = data_task.trialinfo(trial_i,4:7);
        sequence = sequence(sequence > 0);

        if length(sequence) == ranks

            idx_length = [idx_length trial_i];

        end

    end

    data_task.time = data_task.time(idx_length);
    data_task.trial = data_task.trial(idx_length);
    data_task.trialinfo = data_task.trialinfo(idx_length,:);

    data_baseline.time = data_baseline.time(idx_length);
    data_baseline.trial = data_baseline.trial(idx_length);
    data_baseline.trialinfo = data_baseline.trialinfo(idx_length,:);

end

%% compute matrix design for linear model

if number_locations == 8

    matrix_design = zeros(size(data_task.trial, 2), number_locations*ranks);

    for trial_i = 1:size(data_task.trial, 2)
        for rank_i = 1:ranks
            for loc_i = 1:number_locations
                if data_task.trialinfo(trial_i,3+rank_i) == loc_i
                    matrix_design(trial_i, loc_i+(number_locations*(rank_i-1))) = 1;
                end
            end
        end
    end


elseif number_locations == 4

    matrix_design = zeros(size(data_task.trial, 2), number_locations*ranks);

    for trial_i = 1:size(data_task.trial, 2)

        for rank_i = 1:ranks

            if shifted_locations == 0
        
                if data_task.trialinfo(trial_i,3+rank_i) == 1
                    pooled_location = 1;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 2
                    pooled_location = 1;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 3
                    pooled_location = 2;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 4
                    pooled_location = 2;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 5
                    pooled_location = 3;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 6
                    pooled_location = 3;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 7
                    pooled_location = 4;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 8
                    pooled_location = 4;
                end

            elseif shifted_locations == 1
        
                if data_task.trialinfo(trial_i,3+rank_i) == 1
                    pooled_location = 1;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 2
                    pooled_location = 2;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 3
                    pooled_location = 2;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 4
                    pooled_location = 3;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 5
                    pooled_location = 3;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 6
                    pooled_location = 4;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 7
                    pooled_location = 4;
                elseif data_task.trialinfo(trial_i,3+rank_i) == 8
                    pooled_location = 1;
                end
            end

            matrix_design(trial_i, pooled_location+(number_locations*(rank_i-1))) = 1;

        end

    end

end

%% linear model - ridge regression

X = nan(number_locations*ranks,size(data_task.label, 1));

for parcel_i = 1:size(data_task.label, 1)

    % response variable
    y = [];

    for trial_i = 1:size(data_task.trial,2)

        y_baseline_corrected = mean(data_task.trial{trial_i}(parcel_i,:)) - mean(data_baseline.trial{trial_i}(parcel_i,:));

        y = [y y_baseline_corrected];

    end

    y = y';

    % ridge regression (with centering
    lambda = 1;
    X(:, parcel_i) = ridge(y, matrix_design, lambda, 1);

end

    


end
