function [dFro, dProc] = subspace_metrics(X1, X2)
% SUBSPACE_METRICS Compute complementary distances between two matrices of
% PC scores (conditions by PCs)
%
% Inputs:
%   X1, X2 : conditions by PCs matrices (e.g., PC score matrices)
%
% Outputs:
%   dFro   : Frobenius norm distance (||X1 - X2||_F)
%   dProc  : Procrustes distance without scaling

%% 1) Frobenius norm distance
Xn = X1 / norm(X1, 'fro');
Yn = X2 / norm(X2, 'fro');
dFro = norm(Xn - Yn, 'fro');

%% 2) Procrustes distance 
% procrustes allows specifying scaling; here we disable it (scaling=1)
[dProc, ~, ~] = procrustes(X2, X1, 'scaling', true, 'reflection', false);

end