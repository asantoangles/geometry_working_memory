function [volume, distance_by_separation] = compute_distance_accum(points)

% compute_distance_accum computes geometric measures from a set of points
% in high-dimensional space, including the convex hull volume and the mean
% pairwise distance between points as a function of their separation.
%
% SYNTAX:
%   [volume, distance_by_separation] = compute_distance_accum(points)
%
% INPUTS:
%   points : [N × D] matrix
%       Coordinates of N points in D-dimensional space (D ≥ 3). Each row
%       represents one point, and each column corresponds to a spatial
%       dimension. This function is designed to handle cases where the
%       dimensionality exceeds three or when the points are not full-rank.
%
% OUTPUTS:
%   volume : scalar
%       The volume enclosed by the convex hull of the input points, computed
%       using MATLAB’s convhulln algorithm. If the convex hull cannot be
%       computed directly (e.g., when the data are rank-deficient or the
%       number of points equals the number of dimensions), the function
%       attempts to recover a valid subset of linearly independent points
%       using QR decomposition before recomputing the hull volume.
%       If all attempts fail, volume is set to 0.
%
%   distance_by_separation : [1 × (N−1)] vector
%       The average Euclidean distance between points as a function of their
%       separation index. For each separation i (from 1 to N−1), the function
%       computes the mean distance between all pairs of points that are i
%       indices apart, wrapping around cyclically to preserve symmetry.
%
% DESCRIPTION:
%   This function generalizes compute_distance.m to support points in spaces
%   of dimensionality greater than 3. It first attempts to compute the convex
%   hull volume directly via convhulln. If this fails—typically due to
%   degeneracy or full-rank limitations—it uses QR decomposition with column
%   pivoting to identify a linearly independent subset of points, augments
%   them with the origin, and recomputes the hull volume on this reduced set.
%
%   The second part computes how Euclidean distance between points varies
%   with their separation index, providing a measure of geometric smoothness
%   or spatial continuity across ordered points in high-dimensional space.
%
%
% EXAMPLE:
%   % Generate random 6D points
%   points = rand(20,6);
%
%   % Compute convex hull volume and distance profile
%   [volume, dist_sep] = compute_distance_accum(points);
%
% REFERENCES:
%   Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. 
%   Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. 
%   bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

%% volume with convex full algorithm

try

    % Compute the convex hull of the points    
    [~, volume] = convhulln(points);

catch
    try

        %% when points are not full-rank, or the number of points is the same than the number of dimensions

        % Transpose to apply QR on rows
        [~, R, E] = qr(points', 'vector');
    
        % Estimate rank with tolerance
        tol = max(size(R)) * eps(norm(R, 'fro'));
        rank_points = sum(abs(diag(R)) > tol);
    
        % Extract the linearly independent rows
        points = points(E(1:rank_points), E(1:rank_points));

        % add zero point
        points(end+1, :) = 0;

        [~, volume] = convhulln(points);

    catch
        volume = 0;
    end
end

%% distance as a function of separation of points

% Compute distances between points for all separations
num_points = size(points, 1);
distance_by_separation = zeros(1, num_points-1);

% loop over separations between points
for i = 1:num_points-1

    % Compute distances for the given separation
    distances = zeros(0, 1);
    for j = 1:num_points

        % Calculate distance with point separated on the right
        k = j + i;
        if k > size(points,1)
            k = k - size(points,1);
        elseif k < 1
            k = k + size(points,1);
        end
        distances(end+1) = norm(points(k, :) - points(j, :));  % Euclidean distance

        % Calculate distance with point separated on the left
        k = j - i;
        if k > size(points,1)
            k = k - size(points,1);
        elseif k < 1
            k = k + size(points,1);
        end
        distances(end+1) = norm(points(k, :) - points(j, :));  % Euclidean distance

    end
    
    % Calculate the average distance for the given separation
    distance_by_separation(i) = mean(distances);

end


end