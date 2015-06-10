function [X, num_zones] = load_grid_dat(direct)
% LOAD_GRID_DAT load grid data an apply length scale for old jet data
%
%   [x, y] = load_grid_dat(direct, l_scale)
%#ok<*AGROW>

% convert from ft to meters
ft2meters = 0.3048;

% Currently hard coded may change
grid_file = fopen([direct filesep 'Raw Data' filesep 'grid.dat']);

% Grid files are organized as one column
grid = fscanf(grid_file, '%g');
fclose(grid_file);

% File entry is number of zones
num_zones = grid(1);    

% Get the size of each zone
zones_dims = reshape(grid(2:2*num_zones+1), 2, num_zones)';  

if isequal(zones_dims, circshift(zones_dims, 1, 1))
    dimensions = zones_dims(1,:);
else
    error('load_grid_dat assumes that zones are of equal size');
end

% Clip off grid information
grid = grid(2*num_zones+2:end);
offset = 0;

x = [];
y = [];
z = [];

for zone = 1:num_zones
    x_temp = reshape(grid((1:prod(dimensions))+offset), dimensions);
    offset = offset + prod(dimensions); % offset x
    
    y_temp = reshape(grid((1:prod(dimensions))+offset), dimensions);
    offset = offset + prod(dimensions); % offset y
    
    z_temp = reshape(grid((1:prod(dimensions))+offset), dimensions);
    offset = offset + prod(dimensions); % offset z

    %  concatinate matrices
    if zone ~= 1
        x = [fliplr(x_temp), x];
        y = [fliplr(y_temp), y];
        z = [fliplr(z_temp), z];
    else
        x = [x_temp, x];
        y = [y_temp, y];
        z = [z_temp, z];
    end
end

X.x = x*ft2meters;
X.y = y*ft2meters;
X.z = z*ft2meters;
end