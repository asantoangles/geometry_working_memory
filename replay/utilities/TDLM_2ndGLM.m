function [Msf, Msb] = TDLM_2ndGLM(X,T,nstates,nlags,nshuf)

% TDLM_2ndGLM estimates forward and backward sequenceness in a decoded neural state space
% using a two-step GLM approach, based on Liu et al. (2019)
%
% [Msf, Msb] = TDLM_2ndGLM(X, T, nstates, nlags, nshuf)
%
% Inputs:
%   X       - (time × nstates) decoded state probabilities over time
%   T       - (nstates × nstates) hypothesized transition matrix
%   nstates - number of states
%   nlags   - number of time lags to consider
%   nshuf   - number of permutations for shuffled transition matrices
%
% Outputs:
%   Msf     - (1 × nshuf × nlags) forward sequenceness
%   Msb     - (1 × nshuf × nlags) backward sequenceness
%
% Sequenceness is computed by:
%   1) Constructing lagged predictors for each state (first GLM)
%   2) Regressing first-level betas onto template transition matrices 
%      to estimate forward and backward sequence contributions (second GLM)
%
% Notes:
%   - First permutation corresponds to the empirical transition matrix
%   - Supports 3- or 4-state sequences with automatic permutation handling
%
% REFERENCES:
%
%   Liu Y, Dolan RJ, Higgins C, Penagos H, Woolrich MW, Ólafsdóttir HF, Barry C, Kurth-Nelson Z, Behrens TE. 
%   Temporally delayed linear modelling (TDLM) measures replay in both animals and humans. Elife. 2021 Jun 7;10:e66917. 
%   doi: 10.7554/eLife.66917. PMID: 34096501; PMCID: PMC8318595.
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
    
%% GLM: state regression, with other lages 

betas = nan(nlags, nstates, nstates); 
vars = nan(nlags, nstates, nstates); 
% columns betas: lags, Y, shifted_Y
% in transition matrix, shifted_Y are rows, Y are columns
% in Tr, the vector is concatenation of columns

bins  = nlags;

for ilag=1:bins

    temp_zinds = (1:bins:nstates*nlags) + ilag - 1; 
    shifted_Y=dm(:,temp_zinds);   
        
    %% GLM
    for istate=1:size(Y,2)
        
        [cope,varcope,~]=ols_TDLM(shifted_Y(:,istate),[ones(length(Y),1),Y]);

        betas(ilag,:,istate) = cope(2:end);
        vars(ilag,:,istate) = varcope(2:end);
        
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
    betaM(:,13:16) = [betas(:,4,1) betas(:,4,2) betas(:,4,3) betas(:,4,4)];
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
