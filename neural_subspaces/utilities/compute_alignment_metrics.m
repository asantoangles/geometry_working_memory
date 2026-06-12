function [angle, vaf, angle_min] = compute_alignment_metrics(planes, ranks)

% COMPUTE_ALIGNMENT_METRICS  
% Quantify subspace alignment between subspaces defined by PC scores.
%
%   [angle, vaf, angle_min] = compute_alignment_metrics(planes, ranks)
%
%   This function measures how similar the low-dimensional neural 
%   subspaces associated with different item ranks are. Two complementary 
%   metrics are computed:
%
%   (1) PRINCIPAL ANGLES: For 2-D subspaces embedded in R^3, the angle 
%       between planes is computed from the cross-product normal vectors. 
%       For subspaces in R^n with n > 3, where no unique normal exists, 
%       principal angles are computed via a QR–SVD approach: each subspace 
%       is orthonormalized using QR decomposition, the inner-product 
%       matrix between bases is formed, and its singular values give the 
%       cosines of the principal angles. Because each subspace is 2-D, 
%       two principal angles are returned (angle_min, angle).
%
%   (2) VARIANCE ACCOUNTED FOR (VAF): VAF(i,j) quantifies how much of the 
%       variance in rank i activity is explained by projecting it onto the 
%       subspace of rank j, expressed as a normalized Frobenius norm.
%
%   INPUTS:
%       planes   Struct from compute_subspaces containing:
%                   • plane_r<k>  – basis vectors for rank k
%                   • score_r<k>  – PC scores for rank k
%                   • Y_r<k>      – demeaned population activity
%
%       ranks    Number of item ranks (typically 3 or 4).
%
%   OUTPUTS:
%       angle       ranks × ranks matrix of principal angles (degrees).
%       angle_min   ranks × ranks matrix of the smaller principal angle 
%                   (only for n > 3 case).
%       vaf         ranks × ranks matrix of VAF values in [0, 1].
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
% 
%       Santo-Angles A, Yang J, Zhou Y, Chu W K H, Lindsay G W, Sreenivasan K K
%       (2025). Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not
%       bioRxiv 2025.09.05.674385; doi: https://doi.org/10.1101/2025.09.05.674385


%% principal angle between rank subspaces

% subspaces defined by 3 PCs, it returns one principal angle
if size(planes.plane_r1, 1) == 3

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

    angle_min = 0;

   
% subspaces defined by more than 3 PCs, it returns two principal angles
else

    angle_min = zeros(ranks,ranks);
    angle = zeros(ranks,ranks);
    
    for rank_i = 1:ranks
    
        eval(['plane1 = planes.plane_r' num2str(rank_i) ';']);
    
        for rank_ii = 1:ranks
    
            eval(['plane2 = planes.plane_r' num2str(rank_ii) ';']);
    
            [U1,~] = qr(plane1,0);
            [U2,~] = qr(plane2,0);
    
            C = U1' * U2;
            sigma = svd(C);
    
            angle_min(rank_i, rank_ii) = min(acosd(min(max(sigma, -1), 1)));
            angle(rank_i, rank_ii) = max(acosd(min(max(sigma, -1), 1)));
    
        end
    
    end

end

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

