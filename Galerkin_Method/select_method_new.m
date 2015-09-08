function methods_X = select_method_new(bnd_idx, bnd_X, dimensions, xi)
% SELECT_METHOD Derivative preprocessing step deteriminng which 
% finite difference method to use at each grid point, should provide 
% large speedup if more than a 2 or 3 variables are calculated
%
%   methods_X = select_method(bnd_idx, bnd_X, dimensions, xi) , select the
%       appropriate finite difference method based on the flow boundaries
%
% Key 
% method 0: no derivative inside boundary
% method 1: 4th order central difference
% method 2: 1st order forward difference
% method 3: 2nd order forward difference
% method 4: 3rd order forward difference
% method 5: 4th order forward difference
% method 6: 3rd order stencil forward biased difference
% method 7: 4th order stencil forward biased difference
% method 8: no derivative 1 point stencil
% method 9: 2nd order central difference
% method 10: 3rd order stencil backwards biased difference
% method 11: 4th order stencil backwards biased difference
% method 12: 1st order backwards difference
% method 13: 2nd order backwards difference
% method 14: 3rd order backwards difference
% method 15: 4th order backwards difference


comps = flow_ncomps(bnd_X);
x = flow_comps(bnd_X);

% Preallocate method matrics
for i = 1:comps
    methods_X.(xi{i}) = zeros(dimensions);
end

start = 1;
index = {1, 1, 1};

if any(strcmp(xi, 'x'))
    methods_X.x = outer_loop(1, index, bnd_idx, dimensions, methods_X.x, start);
end

if any(strcmp(xi, 'xi'))
    methods_X.xi = outer_loop(1, index, bnd_idx, dimensions, methods_X.xi, start);
end

if any(strcmp(xi, 'y'))
    methods_X.y = outer_loop(2, index, bnd_idx, dimensions, methods_X.y, start);
end

if any(strcmp(xi, 'eta'))
    methods_X.eta = outer_loop(2, index, bnd_idx, dimensions, methods_X.eta, start);
end

if any(strcmp(xi, 'z'))
    methods_X.z = outer_loop(3, index, bnd_idx, dimensions, methods_X.z, start);
end

if any(strcmp(xi, 'phi'))
    methods_X.phi = outer_loop(3, index, bnd_idx, dimensions, methods_X.phi, start);
end

% Find all points on a boundary
solid_bnd = double(bnd_idx == 0);
range = [-2, -1, 1, 2];

% Remove those that are within 2 indices of the edge
for i = 1:comps
    idx_front = flow_index([1 2], i, solid_bnd);
    idx_back = flow_index([-1 0], i, solid_bnd);
    solid_bnd(idx_front{:}) = 0;
    solid_bnd(idx_back{:}) = 0;
end

% Remove all portions of that are imaginary boundaries
for i = 1:comps
    solid_bnd(bnd_X.(x{i}) ~= 0) = 0;
end

% Change all phyiscal boundaries finite different to 4th order central
for i = 1:comps
    central = zeros(dimensions);

    shift = zeros(1, comps);
    shift(i) = 1;

    for j = 1:length(range)
        central(logical(circshift(solid_bnd, range(j)*shift))) = 1;
    end

    central(bnd_idx <= 0) = 0;
    idx_front = flow_index([1 2], i, solid_bnd);
    idx_back = flow_index([-1 0], i, solid_bnd);

    central(idx_front{:}) = 0;
    central(idx_back{:}) = 0;
    methods_X.(xi{i})(logical(central)) = 9;
end
end

function methods = outer_loop(target, index, bnd_idx, dimensions, methods, start)
% OUTER_LOOP loop through the other dimensions step in for additional loops
% when necessary
%
%   methods = outer_loop(target, index, bnd_idx, dimensions, methods,
%       start)  


% Edge case when target is in the last position
if length(dimensions) == target
    stop = length(dimensions) - 1;
else
    stop = length(dimensions);
end

% Skip target dimension in outer loops
if start == target
    start = start + 1;
end

% set index
i = 1;
index{start} = i;

% Loop and step in
if start ~= stop
    while i <= dimensions(start)
        methods = outer_loop(target, index, bnd_idx, dimensions, methods, start + 1);
        i = i + 1;
        index{start} = i;
    end
else
    while i <= dimensions(start)
        methods = inner_loop(target, index, bnd_idx, dimensions, methods);
        i = i + 1;
        index{start} = i;
    end
end
end


function methods = inner_loop(target, index, bnd_idx, dimensions, methods)
% INNER_LOOP actually select methods to be used at each point along the
% mesh used by outer_loop
%
%   methods = INNER_LOOP(target, index, bnd_idx, dimensions, methods)

% Get index to create temporary vector
linear_index = index;
linear_index{target} = 1:dimensions(target);

% Create tempoary values
bnd_temp = bnd_idx(linear_index{:});
methods_temp = methods(linear_index{:});

% set index
i = 1;

% Loop through elements oooof this particular slice
while i < dimensions(target)
    
    % If inside/on boundary continue
    if bnd_temp(i) <= 0
        i = i + 1;
        continue;
    end
    
    % Get current index 
    i_temp = i;
    
    % Find the consecutive points inside a boundary
    while i <= dimensions(target) && bnd_temp(i) > 0
        i = i + 1;
    end
    
    % Determine size of gap
    gap = i - i_temp;
    
    % Fill row of methods with current method number
    switch gap
        case 1
            methods_temp(i_temp) = 8;
        case 2
            methods_temp(i_temp) = 2;
            methods_temp(i_temp + 1) = 12;
        case 3
            methods_temp(i_temp) = 3;
            methods_temp(i_temp + 1) = 9;
            methods_temp(i_temp + 2) = 13;
        case 4
            methods_temp(i_temp) = 4;
            methods_temp(i_temp + 1) = 6;
            methods_temp(i_temp + 2) = 10;
            methods_temp(i_temp + 3) = 14;
        otherwise
            methods_temp(i_temp) = 5;
            methods_temp(i_temp + 1) = 7;
            for j = 1:gap-4
                methods_temp(i_temp + 1 + j) = 1;
            end
            methods_temp(i_temp + gap - 2) = 11;
            methods_temp(i_temp + gap - 1) = 15;
    end
end

methods(linear_index{:}) = methods_temp;
end