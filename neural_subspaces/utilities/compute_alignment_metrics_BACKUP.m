function [angle, vaf] = compute_alignment_metrics(planes, ranks)

% COMPUTE_ALIGNMENT_METRICS  Compute subspace alignment and variance sharing across ranks.
%
%   [angle, vaf] = compute_alignment_metrics(X, planes, ranks)
%
%   This function quantifies how similar the low-dimensional neural
%   subspaces are across different item-ranks in the task. It computes:
%
%       (1) PRINCIPAL ANGLE: the geometric angle between 2-D planes estimated for 
%           each rank, using the normal vectors of those planes.
%
%       (2) VAF (Variance Accounted For): how much variance in the neural
%           population activity of one rank can be explained by the
%           principal components defining the subspace of another rank.
%
%   These metrics together characterize the degree to which rank-specific
%   neural representations lie in shared or distinct subspaces.
%
%   INPUTS:
%       planes   Struct returned by compute_subspaces, containing:
%                   • plane_r<k> – 2D basis vectors for rank k
%                   • score_r<k> – projected rank-specific scores
%                   • Y_r<k>     – demeaned population activity
%
%       ranks    Number of ranks (3 or 4).
%
%   OUTPUTS:
%       angle    ranks × ranks matrix of pairwise angles (degrees)
%                between rank-specific 2-D planes.
%                Values are bounded to [0, 90] degrees.
%
%       vaf      ranks × ranks matrix of Variance Accounted For (VAF):
%                       vaf(i,j) = variance in rank i activity
%                                  explained by rank j’s subspace.
%                Values are bounded to [0, 1].
%
%   ------------------------------------------------------------------
%   REFERENCES:
%       Xie Y, Hu P, Li J, Chen J, Song W, Wang XJ, Yang T, Dehaene S, Tang S, Min B, Wang L. 
%       Geometry of sequence working memory in macaque prefrontal cortex. 
%       Science. 2022 Feb 11;375(6581):632-639. doi: 10.1126/science.abm0204. 
%       Epub 2022 Feb 10. PMID: 35143322.
% 
%       Panichello MF, Buschman TJ. Shared mechanisms underlie the control of working memory and attention. 
%       Nature. 2021 Apr;592(7855):601-605. doi: 10.1038/s41586-021-03390-w. 
%       Epub 2021 Mar 31. PMID: 33790467; PMCID: PMC8223505.

%% principal angle between rank subspaces

angle = zeros(ranks,ranks);

for rank_i = 1:ranks

    eval(['plane1 = planes.plane_r' num2str(rank_i) ';']);

    for rank_ii = 1:ranks

        eval(['plane2 = planes.plane_r' num2str(rank_ii) ';']);

        cross_product1 = cross(plane1(:,1), plane1(:,2));
        cross_product2 = cross(plane2(:,1), plane2(:,2));

        % Calculate the dot product of the cross products
        dot_product = dot(cross_product1, cross_product2);
        
        % Calculate the magnitudes of the cross products
        magnitude1 = norm(cross_product1);
        magnitude2 = norm(cross_product2);
        
        % Calculate the cosine of the angle between the two planes
        angle(rank_i, rank_ii) = dot_product / (magnitude1 * magnitude2);

    end

end

% convert to degrees
angle = real(acosd(angle));

% reduce to range [0 90]
% if angle > 90, 
angle(angle > 90) = abs(180 - angle(angle > 90));

%% variance accounted for (VAF)

vaf = zeros(ranks,ranks);

for rank_i = 1:ranks

    eval(['plane1 = planes.plane_r' num2str(rank_i) ';']);
    eval(['score1 = planes.score_r' num2str(rank_i) ';']);
    eval(['data1 = planes.Y_r' num2str(rank_i) ';']);

    for rank_ii = 1:ranks

        eval(['plane2 = planes.plane_r' num2str(rank_ii) ';']);
        eval(['score2 = planes.score_r' num2str(rank_ii) ';']);

        vaf(rank_i, rank_ii) = (norm(plane2*plane2'*plane1*score1','fro')/norm(plane1*score1','fro')).^2;

    end

end

end

