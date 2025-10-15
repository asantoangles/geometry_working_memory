%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: convertLaserFiles.m
%
% Converts the output of the laser to .elp and .hsp files.
%   Prints the final headshape at the end
%
% Inputs:
%   * electrodes_file: The .txt file with the laser fiducial and sensor
%       points. 
%       * NOTE: For now, we're using the same order to get the points as
%           before.
%   * headshape_file: The .txt file with the laser headshape points.
%   * out_file_name: The name of the .elp and .hsp files. They will both
%       have this name, with different file extensions.
%   * cluster_distance: The distance between clusters for outliers in mm
%       (optional)
%   * num_clusters: Number of clusters to keep (optional)
%
% Outputs: None, but will create both .hsp and .elp files.
%
% Usage: convertLaserFiles('Laser_All_Points.txt','Laser_All_Headshape.txt','Laser_All')
%
% Author: Doug Bemis
% Date: 3/1/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convertLaserFiles(electrodes_file, headshape_file, out_file_name, cluster_distance, num_clusters)

% The distance between clusters
if ~exist('cluster_distance','var')
    cluster_distance = 10;
end

% The number of clusters to keep
if ~exist('num_clusters','var')
    num_clusters = 1;
end

% Read the electrodes in 
disp('Getting transformation...');
fid = fopen(electrodes_file);

% These are just headers
for s = 1:3
    fgetl(fid);
end

% Then, the fiducials
laser_fids = zeros(3,3);
for f = 1:3
    line = fgetl(fid);
    [laser_fids(f,1) laser_fids(f,2) laser_fids(f,3)] = strread(line,'%f%f%f'); 
end

% And the sensors
laser_sensors = zeros(5,3);
for s = 1:5
    line = fgetl(fid);
    [laser_sensors(s,1) laser_sensors(s,2) laser_sensors(s,3)] = strread(line,'%f%f%f'); 
end
fclose(fid);

% Now, get the transformation
transform = getTransformation(laser_fids);

% Convert the fiducials and sensors
new_fids = transformPoints(laser_fids, transform);
new_sensors = transformPoints(laser_sensors, transform);

% And write the elp file
disp('Writing elp file...');
writeELPFile([out_file_name '.elp'], new_fids, new_sensors);

% Read the headshape in
fid = fopen(headshape_file);

% Ignore the headers for now
for h = 1:3
    fgetl(fid);
end

% Get the points
disp('Reading headshape points...');
line = fgetl(fid);
[perc line] = strtok(line);
num_points = strtok(line);
num_points = str2double(num_points);
laser_points = fscanf(fid,'%f%f%f',[3,num_points]);

% Transform them
new_points = transformPoints(laser_points', transform);

% And write the hsp file
writeHSPFile([out_file_name '.hsp'], new_fids, new_points, cluster_distance, num_clusters);
disp('Done.');


% Helper to write out the elp file
function writeHSPFile(file_name, fids, points,cluster_distance, num_clusters)

% Good for later
num_points = size(points,1);

% Only allowed 5000, so pare it down here
if num_points > 5000
    tmp = points;
    points = zeros(0,3);
    spacing = num_points / 5000;
    next_point = spacing;
    for p = 1:num_points
        if p >= next_point
            points(end+1,:) = tmp(p,:); %#ok<AGROW>
            next_point = next_point+spacing;
        end
    end
    num_points = size(points,1);
    if num_points > 5000
        error('Still too many points. Exiting...');
    end
end

% Need to find the outliers
points = removeOutliers(points,cluster_distance,num_clusters);
num_points = length(points);

% Write the header
disp('Writing hsp file...');
fid = fopen(file_name,'w');
fprintf(fid,'3\t200\n');
fprintf(fid,'//Shape file\n');
fprintf(fid,'//Minor revision number\n');
fprintf(fid,'2\n');
fprintf(fid,'//Subject Name\n');
fprintf(fid,'%%N\tName    \n');
fprintf(fid,'////Shape code, number of digitized points\n');
fprintf(fid,['0\t' num2str(num_points) '\n']);
fprintf(fid,'//Position of fiducials X+, Y+, Y- on the subject\n');

% Now the fiducials
for f = 1:3
    fprintf(fid,'%%F\t');
    for p = 1:3
        fprintf(fid,[num2str(fids(f,p)) '\t']);
    end
    fprintf(fid,'\n');
end


% And the points
fprintf(fid,'//No of rows, no of columns; position of digitized points\n');
fprintf(fid,[num2str(num_points) '\t3\n']);
fprintf(fid,'%f %f %f\n',points');
fclose(fid);

% And show the output
close(gcf);
plot3(points(:,1),points(:,2),points(:,3),'b.');



% Helper to write out the elp file
function writeELPFile(file_name, fids, sensors)

% Write the header
fid = fopen(file_name,'w');
fprintf(fid,'3\t2\n');
fprintf(fid,'//Probe file\n');
fprintf(fid,'//Minor revision number\n');
fprintf(fid,'1\n');
fprintf(fid,'//ProbeName\n');
fprintf(fid,'%%N\tName    \n');
fprintf(fid,'//Probe type, number of sensors\n');
fprintf(fid,'0\t5\n');
fprintf(fid,'//Position of fiducials X+, Y+, Y- on the subject\n');

% Now the fiducials
for f = 1:3
    fprintf(fid,'%%F\t');
    for p = 1:3
        fprintf(fid,[num2str(fids(f,p)) '\t']);
    end
    fprintf(fid,'\n');
end

% And the sensors
sensor_names = {'RED','YELLOW','BLUE','WHITE','BLACK'};
for s = 1:5
    fprintf(fid,'//Sensor type\n');
    fprintf(fid,'%%S\t4000\n');
    fprintf(fid,['//Sensor name and data for sensor # ' num2str(s) '\n']);
    fprintf(fid,['%%N\t' num2str(s-1) '-' sensor_names{s} '   \n']);
    for p = 1:3
        fprintf(fid,[num2str(sensors(s,p)) '\t']);
    end
    fprintf(fid,'\n');
end
fclose(fid);


% Helper to transform some points
% Expect points to be rows
function points = transformPoints(points, transform)

% Scale first
points = points*transform.scale;

% Then translate
points = points - repmat(transform.origin, size(points,1),1);

% And rotate
points = points*transform.rotation;


% Get the transformation from laser coordinates to digitizer coordinates
function transform = getTransformation(fids)

% This appears to be right, but we should probably test a little more
transform.scale = 0.001;

% First scale
fids = transform.scale*fids;

% Then move the origin
% It's always midway between the left and right fiducial
transform.origin = (fids(2,:)+fids(3,:))/2;
fids = fids - repmat(transform.origin,3,1);

% Get the axes - (the l is for laser)
l_axes = zeros(3,3);

% x-axis - should be straight through the first point
l_axes(1,:) = fids(1,:); l_axes(1,:) = l_axes(1,:) / norm(l_axes(1,:));

% z-axis - perpendicular to x-axis and two other fiducials
l_axes(3,:) = cross(l_axes(1,:),fids(2,:)-fids(3,:)); l_axes(3,:) = l_axes(3,:) / norm(l_axes(3,:));

% y-axis - perpendicular to x and z axes
l_axes(2,:) = cross(l_axes(3,:), l_axes(1,:));  l_axes(2,:) = l_axes(2,:) / norm(l_axes(2,:));

% Get the rotation to upright by solving the rotation equation
transform.rotation = [1 0 0; 0 1 0; 0 0 1] / l_axes;


% Quick clustering function to remove outliers
function points = removeOutliers(points, cluster_distance, num_clusters)

% Get all the distances first
disp(['Clustering outliers [' num2str(cluster_distance) 'mm; ' num2str(num_clusters) ' cluster(s)]...']);
dist = zeros(length(points(:,1)),length(points(:,1)));
for p = 1:length(points)
    for d = 1:3
        dist(:,p) = dist(:,p) + (points(:,d) - points(p,d)) .^ 2; 
    end
end
dist = sqrt(dist);

% Store which cluster each is in
clusters = zeros(1,length(points));

% Make the clusters
for p = 1:length(points)
    
    % Check all that have been assigned
    p_clusters = clusters(dist(p,clusters>0) < cluster_distance/1000);
    
    % If none are close, make a new one
    if isempty(p_clusters)
        clusters(p) = max(clusters)+1;
    else
        
        % Set to the first one
        clusters(p) = p_clusters(1);
        
        % Set all the others as well
        for c = 2:length(p_clusters)
            clusters(clusters == p_clusters(c)) = p_clusters(1);
        end
    end    
end

% Find cluster sizes
c_ind = unique(clusters);
cluster_sizes = zeros(1,length(c_ind));
for c = 1:length(c_ind)
    cluster_sizes(c) = sum(clusters == c_ind(c));
end
[val ind] = sort(cluster_sizes,'descend');

% Put all the clusters to reject in one 
max_cluster_ind = ind(1);
max_cluster = c_ind(max_cluster_ind);
for c = 1:num_clusters
    clusters(clusters == c_ind(ind(c))) = max_cluster; 
end

% Now, plot to see and get the sizes
colors = {'g','r','m','b','c','y'};
for c = 1:length(c_ind)
    p = points(clusters == c_ind(c),:);
    if c == max_cluster_ind
        plot3(p(:,1),p(:,2),p(:,3),'k.');
    else
        hold on;
        plot3(p(:,1),p(:,2),p(:,3),[colors{mod(c,length(colors))+1} 'o']);        
    end
end

% See if we want to reject them
c = '';
while ~strcmpi(c,'y') && ~strcmpi(c,'n')
    % c = input('Reject points (y/n)? ','s');
    c = 'y'
end

% If not, then see why...
if strcmpi(c,'n')

    % See if we want to try again
    c = '';
    while ~strcmpi(c,'y') && ~strcmpi(c,'n')
        c = input('Try again with different parameters (y/n)? ','s');
    end

    % If so, get the new distance
    if strcmpi(c,'y')
        cl_d = '';
        while isnan(str2double(cl_d))
            cl_d = input('Enter new distance (mm): ','s');
        end
        cl_n = '';
        while isnan(str2double(cl_n))
            cl_n = input('Enter new number of clusters to keep: ','s');
        end
        close(gcf);
        clear dist;
        points = removeOutliers(points, str2double(cl_d), str2double(cl_n));
    end
    
% Otherwise, only keep non-outliers
else
    points = points(clusters == max_cluster,:);
end




