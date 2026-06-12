function [volume, distance_by_separation] = compute_separability_metrics(points)

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
        volume = volume + tetrahedron_volume;
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