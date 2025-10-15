function [volume, distance_by_separation] = compute_distance(points)

% compute_distance computes two geometric measures from a set of points 
% defining PC scores (rows as samples, columns as PCs):
% the volume enclosed by their convex hull and the mean pairwise distance
% between points as a function of their separation.
%
% SYNTAX:
%   [volume, distance_by_separation] = compute_distance(points)
%
% INPUTS:
%   points : [N × D] matrix
%       Coordinates of N points in D-dimensional space. Each row represents
%       one point, and each column corresponds to a spatial dimension (PC)
%
% OUTPUTS:
%   volume : scalar
%       The volume enclosed by the convex hull of the input points,
%       computed using a full convex hull triangulation approach. 
%       If the convex hull cannot be computed (e.g., due to degenerate 
%       configurations), volume is returned as an empty array [].
%
%   distance_by_separation : [1 × (N−1)] vector
%       The average Euclidean distance between points as a function of their
%       separation in the ordered list of points. For each possible 
%       separation i (from 1 to N−1), the function computes the mean 
%       distance between all pairs of points that are i indices apart, 
%       wrapping around cyclically to ensure symmetry.
%
% DESCRIPTION:
%   The function first attempts to compute the convex hull of the point 
%   cloud using MATLAB’s convhulln algorithm. The convex hull defines the 
%   smallest convex polytope enclosing all points. Its volume is estimated 
%   as the sum of the volumes of tetrahedra formed between each triangular 
%   face of the hull and the origin.
%
%   In the second part, the function computes how the Euclidean distance 
%   between points depends on their index-based separation, providing a 
%   measure of spatial organization or smoothness of the point set.
%
% EXAMPLE:
%   % Generate random 3D points
%   points = rand(10,3);
%
%   % Compute convex hull volume and distance profile
%   [volume, dist_sep] = compute_distance(points);
%
% REFERENCES:
%   Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. 
%   Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. 
%   bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

%% volume with convex full algorithm

try

    % Compute the convex hull of the points
    convex_hull = convhulln(points);
    
    % Calculate the volume of the convex hull
    volume = 0;
    for i = 1:size(convex_hull, 1)
    
        % Extract vertices of a triangular face
        v1 = points(convex_hull(i, 1), :);
        v2 = points(convex_hull(i, 2), :);
        v3 = points(convex_hull(i, 3), :);
        
        % Calculate the volume of the tetrahedron formed by the face and the origin
        tetrahedron_volume = abs(dot(v1, cross(v2, v3))) / 6;
        
        % Accumulate the volume
        vsolume = volume + tetrahedron_volume;
    end

catch
    volume = [];
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