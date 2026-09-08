function eigvecs = compute_eigvecs(X)

%% decomposition of X with SVD (PCA)

warning off

% Calculate the covariance matrix of the data
covariance_matrix = cov(X);

% Perform Singular Value Decomposition (SVD) on the covariance matrix
[eigvecs, ~, ~] = svd(covariance_matrix);

warning on

end