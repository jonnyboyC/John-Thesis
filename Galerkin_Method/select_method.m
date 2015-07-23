function [methods_X] = select_method(bnd_idx, bnd_X, dimensions, xi, closed_bnd)
% SELECT_METHOD Derivative preprocessing step in deterimine which 
% finite difference method to use at each grid point, should provide 
% large speedup if more than a 2 or 3 variables are calculated
%
% [methodsX, methodsY] = select_method(bnd_idx, bnd_x, bnd_y, dimensions,
% closed_bnd), select the appropriate finite difference method based on the
% flow boundaries
%
% Key 
% method 0: no derivative inside boundary
% method 1: 1st order forward difference
% method 2: 2nd order forward difference
% method 3: 3rd order forward difference
% method 4: 4th order forward difference
% method 5: 3rd order stencil forward biased difference
% method 6: 4th order stencil forward biased difference
% method 7: no derivative 1 point stencil
% method 8: 2nd order central difference
% method 9: 4th order central difference
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

i = 1;
j = 1;
gap = 0;

% Main loop y direction derivatives
while j <= dimensions(2)
    while i <= dimensions(1) 
        
        % If in boundary leave method as 0 and continue
        if bnd_idx(i,j) <= 0
            i = i + 1;
            continue;
        end
        
        % Get current index
        i_temp = i;
        
        % Determine number of spaces until a boundary or edge of image
        while i <= dimensions(1) && bnd_idx(i,j) > 0
            gap = gap + 1;
            i = i + 1;
        end
        
        % Assign appropriate methods
        switch gap
            case 1
                methods_X.(xi{1})(i_temp, j)   = 7;
            case 2
                methods_X.(xi{1})(i_temp, j)   = 1;
                methods_X.(xi{1})(i_temp+1, j) = 12;
            case 3
                methods_X.(xi{1})(i_temp, j)   = 2;
                methods_X.(xi{1})(i_temp+1, j) = 8;
                methods_X.(xi{1})(i_temp+2, j) = 13;
            case 4
                methods_X.(xi{1})(i_temp, j)   = 3;
                methods_X.(xi{1})(i_temp+1, j) = 5;
                methods_X.(xi{1})(i_temp+2, j) = 10;
                methods_X.(xi{1})(i_temp+3, j) = 14;
            otherwise 
                methods_X.(xi{1})(i_temp, j)   = 4;
                methods_X.(xi{1})(i_temp+1, j) = 6;
                for k = 1:gap-4
                    methods_X.(xi{1})(i_temp+1+k, j) = 9;
                end
                methods_X.(xi{1})(i_temp+gap-4+2, j) = 11;
                methods_X.(xi{1})(i_temp+gap-4+3, j) = 15; 
        end
        gap = 0;
    end
    j = j + 1;
    i = 1;
end

i = 1;
j = 1;
gap = 0;

% Main loop x direction derivatives
while i <= dimensions(1) 
    while j <= dimensions(2)
        % If in boundary leave method as 0 and continue
        if bnd_idx(i,j) <= 0
            j = j + 1;
            continue;
        end
        
        % Get current index
        j_temp = j;
        
        % Determine number of spaces until a boundary or edge of image
        while j <= dimensions(2) && bnd_idx(i,j) > 0 
            gap = gap + 1;
            j = j + 1;
        end
        
        % Assign appropriate methods
        switch gap
            case 1
                methods_X.(xi{2})(i, j_temp)   = 7;
            case 2
                methods_X.(xi{2})(i, j_temp)   = 1;
                methods_X.(xi{2})(i, j_temp+1) = 12;
            case 3
                methods_X.(xi{2})(i, j_temp)   = 2;
                methods_X.(xi{2})(i, j_temp+1) = 8;
                methods_X.(xi{2})(i, j_temp+2) = 13;
            case 4
                methods_X.(xi{2})(i, j_temp)   = 3;
                methods_X.(xi{2})(i, j_temp+1) = 5;
                methods_X.(xi{2})(i, j_temp+2) = 10;
                methods_X.(xi{2})(i, j_temp+3) = 14;
            otherwise 
                methods_X.(xi{2})(i, j_temp)   = 4;
                methods_X.(xi{2})(i, j_temp+1) = 6;
                for k = 1:gap-4
                    methods_X.(xi{2})(i, j_temp+1+k) = 9;
                end
                methods_X.(xi{2})(i, j_temp+gap-4+2) = 11;
                methods_X.(xi{2})(i, j_temp+gap-4+3) = 15; 
        end
        gap = 0;
    end
    i = i + 1;
    j = 1;
end

if closed_bnd
    % Determine portion of boundary represented by solid surface
    solid_bnd = double(bnd_idx == 0);
    
    for i = 1:comps
        solid_bnd(bnd_X.(x{i}) ~= 0) = 0;
    end
    
    centralX = zeros(dimensions);
    centralY = zeros(dimensions);
    
    % Index matrix of two indexes to the left and right of solid boundary
    centralX(logical(circshift(solid_bnd, [1 0]) + circshift(solid_bnd, [2 0]) +...
        circshift(solid_bnd, [-1 0]) + circshift(solid_bnd, [-2 0]))) = 1;
    
    % Index matrix of two indexes to the above and below of solid boundary
    centralY(logical(circshift(solid_bnd, [0 1]) + circshift(solid_bnd, [0 2]) +...
        circshift(solid_bnd, [0 -1]) + circshift(solid_bnd, [0 -2]))) = 1;
    
    % remove points in the boundary and by the edge of image
    centralX(bnd_idx <=0) = 0;
    centralX(1:2, :) = 0;
    centralX(end-1:end, :) = 0;
    
    centralY(bnd_idx <=0) = 0;
    centralY(:, 1:2) = 0;
    centralY(:, end-1:end) = 0;

    % Set these methods to 4th order central
    methods_X.(xi{2})(logical(centralX)) = 9;
    methods_X.(xi{1})(logical(centralY)) = 9;
end
end