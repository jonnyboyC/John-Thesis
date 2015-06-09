function [x, y] = load_grid_dat(direct, l_scale)


% Currently hard coded may change
grid_file = fopen([direct filesep 'Raw Data' filesep 'grid.dat']);

% Grid files are organized as one column
grid = fscanf(grid_file, '%g');

% File entry is number of zones
num_zones = grid(1);    

% Get the size of each zone
zones_dims = reshape(grid(2:2*num_zones+1), 2, num_zones);  

if isequal(zones_dims, circshift(zones_dims, 1, 1))
    dimensions = zones_dims(:,1);
else
    error('load_grid_dat assumes that zones are of equal size');
end

% Clip off grid information
grid = grid(2*num_zones+2:end);
offset = 0;

x = [];
y = [];

for i = 1:num_zones
    x_temp = reshape(grid([1:prod(dimensions)+offset), dimensions);
    offset = offset + prod(dimensions); % offset x
    
    y_temp = reshape(grid([1:prod(dimensions)+offset), dimensions);
    offset = offset + prod(dimensions); % offset y
    offset = offset + prod(dimensions); % offset z

    % Add together concatinate matrices
    x = [x; x_temp];
    y = [y; y_temp];
end

x = x*l_scale;
y = y*l_scale;
end