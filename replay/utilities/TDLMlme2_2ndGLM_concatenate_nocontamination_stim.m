function [Msf, Msb] = TDLMlme2_2ndGLM_concatenate_nocontamination_stim(X,T,nstates,nlags,nshuf,X_trials)

% TDLMlme2_2ndGLM_concatenate_nocontamination_stim is a version of 
% TDLMlme2_2ndGLM_concatenate_nocontamination adapted for the stimulus period. 
% For full documentation, see TDLMlme2_2ndGLM_concatenate_nocontamination.m
%
%   Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. 
%   Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. 
%   bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

if size(X,2) == 3
    uniquePerms = [...
   1   2   3;...
   1   3   2;...
   3   1   2;...
   2   1   3;...
   3   2   1;...
   2   3   1];
elseif size(X,2) == 4
    uniquePerms = [...
   1   2   3   4;...
   3   1   4   2;...
   4   2   3   1;...
   3   1   2   4;...
   2   4   1   3;...
   1   4   2   3;...
   4   1   3   2;...
   2   3   4   1;...
   1   3   2   4;...
   1   4   3   2;...
   4   2   1   3;...
   3   4   1   2;...
   2   4   3   1;...
   1   2   4   3;...
   4   3   1   2;...
   2   1   3   4;...
   4   3   2   1;...
   4   1   2   3;...
   2   1   4   3;...
   3   2   4   1;...
   2   3   1   4;...
   1   3   4   2;...
   3   4   2   1;...
   3   2   1   4];
end

nbins=nlags+1;

if nlags == 500
    trial_tp = 1200; % timepoints
else
    trial_tp = 241;  % bins
end

warning off
dm=[toeplitz(X(:,1),[zeros(nbins,1)])];
dm=dm(:,2:end);

for kk=2:nstates
    temp=toeplitz(X(:,kk),[zeros(nbins,1)]);
    temp=temp(:,2:end);
    dm=[dm temp]; 
end

warning on
Y=X;

Msf=nan(1,nshuf,nlags);
Msb=nan(1,nshuf,nlags);

if size(uniquePerms,1)<nshuf
    nshuf=size(uniquePerms,1);
end

%%% find contaminated regions
number_trials = size(dm,1) / trial_tp;
contaminated_rows = [];
for trial_i = 1:number_trials
    contaminated_rows = [contaminated_rows (1+(trial_tp*(trial_i-1))):(trial_tp*(trial_i-1)+nlags)];
end

%%% set to zero contaminated regions in shifted state space

% split dm
if size(X,2) == 3
    dm1 = dm(:,1:nlags);
    dm2 = dm(:,(1:nlags)+nlags);
    dm3 = dm(:,(1:nlags)+nlags*2);
else
    dm1 = dm(:,1:nlags);
    dm2 = dm(:,(1:nlags)+nlags);
    dm3 = dm(:,(1:nlags)+nlags*2);
    dm4 = dm(:,(1:nlags)+nlags*3);
end

for dm_i = 1:size(X,2)

    eval(['dm = dm' num2str(dm_i) ';'])

    last_column = 0;
    
    % for each row in dm
    for row_i = 1:size(dm)
        
        % if the row is contaminated
        if ismember(row_i,contaminated_rows) == 1
    
            idx = ismember(contaminated_rows,row_i);
            row_cont_i = contaminated_rows(idx);
        
            % if consecutive contaminated rows (meaning, same trial)
            if ismember(row_cont_i-1,contaminated_rows)
    
                last_column = last_column + 1;
    
            else
    
                last_column = 1;
    
            end
    
            % set to zero
            dm(row_i,last_column:end) = 0;
        
        end
    
    end

    eval(['dm' num2str(dm_i) ' = dm;']);

end

% merge dm's
if size(X,2) == 3
    dm = [dm1 dm2 dm3];
else
    dm = [dm1 dm2 dm3 dm4];
end


%% GLM: state regression, with other lages 

betas = nan(nlags, nstates, nstates); 
vars = nan(nlags, nstates, nstates);

% columns betas: lags, Y, shifted_Y
% in transition matrix, shifted_Y are rows, Y are columns
% in Tr, the vector is concatenation of columns

bins  = nlags;

% prepare data matrix for lme
subject_variable = nan(size(Y,1),1);
session_variable = nan(size(Y,1),1);
last_idx = 0;
for sub_i = 1:size(X_trials,1)
    for ses_i = 1:size(X_trials,2)
        total_tp = trial_tp * X_trials(sub_i, ses_i);
        subject_variable((last_idx+1):(last_idx+total_tp),1) = sub_i;
        session_variable((last_idx+1):(last_idx+total_tp),1) = ses_i;
        last_idx = last_idx + total_tp;
    end
end

total_trials = reshape(X_trials.', [], 1); % vectorize trials per session
trial_variable = ones(size(subject_variable));
last_tp = 0;
for session_i = 1:length(total_trials)

    number_trials = total_trials(session_i);

    for trial_i = 1:number_trials
        trial_variable(((1:trial_tp)+(trial_tp*(trial_i-1)))+last_tp) = trial_i;
    end

    last_tp = last_tp + (trial_tp*number_trials);

end

% Convert the predictor variable to categorical
subject_variable = categorical(subject_variable);
session_variable = categorical(session_variable);
trial_variable = categorical(trial_variable);

for ilag=1:bins

    temp_zinds = (1:bins:nstates*nlags) + ilag - 1; 
    shifted_Y=dm(:,temp_zinds);   
        
    %% GLM
    for istate=1:size(Y,2)
        
        % lme
        if size(Y,2) == 3
            varnames_lme = {'shifted_X', 'X1', 'X2', 'X3', 'subject', 'session', 'trial'};
            table_lme = table(shifted_Y(:,istate), Y(:,1), Y(:,2), Y(:,3), subject_variable, session_variable, trial_variable, 'VariableNames', varnames_lme);
            lme = fitlme(table_lme, 'shifted_X  ~ X1 + X2 + X3 + (1|subject) + (1|session) + (1|trial)');
        else
            varnames_lme = {'shifted_X', 'X1', 'X2', 'X3', 'X4', 'subject', 'session', 'trial'};
            table_lme = table(shifted_Y(:,istate), Y(:,1), Y(:,2), Y(:,3), Y(:,4), subject_variable, session_variable, trial_variable, 'VariableNames', varnames_lme);
            lme = fitlme(table_lme, 'shifted_X  ~ X1 + X2 + X3 + X4 + (1|subject) + (1|session) + (1|trial)');
        end

        betas(ilag,:,istate) = lme.Coefficients.Estimate(2:end);
        vars(ilag,:,istate) = lme.Coefficients.SE(2:end);
        
    end
        
end

% reshape to lags by nstates*nstates
% columns of betaM: 1-3 are Y 1, 4-6 are Y 2
% columns betas: lags, Y, shifted_Y
betaM = nan(nlags, nstates*nstates);

if nstates == 3
    betaM(:,1:3) = [betas(:,1,1) betas(:,1,2) betas(:,1,3)];
    betaM(:,4:6) = [betas(:,2,1) betas(:,2,2) betas(:,2,3)];
    betaM(:,7:9) = [betas(:,3,1) betas(:,3,2) betas(:,3,3)];
elseif nstates == 4
    betaM(:,1:4)    = [betas(:,1,1) betas(:,1,2) betas(:,1,3) betas(:,1,4)];
    betaM(:,5:8)    = [betas(:,2,1) betas(:,2,2) betas(:,2,3) betas(:,2,4)];
    betaM(:,9:12)   = [betas(:,3,1) betas(:,3,2) betas(:,3,3) betas(:,3,4)];
    betaM(:,13:16)  = [betas(:,4,1) betas(:,4,2) betas(:,4,3) betas(:,4,4)];
end

%% GLM2

% Define template matrices (Tr)
Tauto = eye(nstates);  % Identity matrix for self-transitions
Tconst = ones(nstates);  % Constant matrix for average of all transitions
Tconf = Tconst - Tauto;

for iShuf = 1:nshuf
    rp = uniquePerms(iShuf,:);
    T1 = T(rp,rp); 
    T2 = T1'; % backwards is transpose of forwards 

    if ~issymmetric(T1)
        Tr = [Tconf(:), T1(:), T2(:)];
    else
        Tr = [Tconf(:), T1(:)];
    end
    
    for ilag=1:nlags
       
        % Perform regression
        b_lag = betaM(ilag, :);
        Z = Tr \ b_lag';
                                         
        Msf(1,iShuf,ilag) = Z(2); % forward (unique contribution)
        Msb(1,iShuf,ilag) = Z(end); % backward (unique contribution)
        
    end
    
end
