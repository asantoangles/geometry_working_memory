function X = compute_X_matrix_lm(data_task, data_baseline, sequence_length, number_locations, shifted_locations)

% COMPUTE_X_MATRIX_LM  Estimate population tuning vectors using ridge regression.
%
%   X = compute_X_matrix_lm(data_task, data_baseline, ...
%                           sequence_length, number_locations, shifted_locations)
%
%   This function computes the *X matrix* used in subsequent neural
%   geometry analyses. Each row of X corresponds to a (location × rank)
%   condition, and each column corresponds to a parcel/channel. The values
%   in X are ridge-regression estimates of baseline-corrected mean activity.
%
%   The function performs the following steps:
%
%       1) **Subset trials** to a specific sequence length (if applicable).
%
%       2) **Construct a design matrix** encoding which location belongs
%          to each rank. For 4-location analyses, the 8 original locations
%          are pooled according to task geometry, with optional circular
%          shifts via `shifted_locations`.
%
%       3) **Compute baseline-corrected responses**:
%
%                   y(trial) = mean(task_signal) − mean(baseline_signal)
%
%       4) **Fit a ridge regression model** for each parcel:
%
%                   y = matrix_design * β   (λ = 1)
%
%          The resulting coefficients β form the rows of X.
%
%   NOTE:
%       • The returned X matrix is **not column-demeaned**.
%
%   INPUTS:
%       data_task        FieldTrip-like struct containing task-evoked data:
%                           • trial{t}   – channels × time samples
%                           • trialinfo  – includes rank × location sequence
%
%       data_baseline    Same structure as data_task, providing baseline
%                       periods matched trial-by-trial.
%
%       sequence_length  Required sequence length. 
%
%       number_locations Number of stimulus locations for each rank (4 or 8).
%
%       shifted_locations Optional (0 or 1). Controls how the 8 original
%                         spatial locations are pooled into 4:
%                           0 – default pooling
%                           1 – circularly shifted pooling
%                         Default: 0.
%
%   OUTPUTS:
%       X                (number_locations × ranks) × n_parcels matrix of
%                        ridge regression coefficients. Each row corresponds
%                        to a (location, rank) condition; each column is a parcel.

if nargin < 5
    shifted_locations = 0;
end


%% subset data_task (trials depending on sequence length)

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
