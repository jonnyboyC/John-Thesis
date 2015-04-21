function [methodsX, methodsY] = select_method(bnd_idx, dimensions)
% Derivative preprocessing step in deterimine which finite difference
% method to use at each grid point, should provide large speedup

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


% Preallocate method matrics
methodsX = zeros(dimensions);
methodsY = zeros(dimensions);

i = 1;
j = 1;
gap = 0;

% Main loop x direction derivatives
while i <= dimensions(1) 
    while j <= dimensions(2)
        % If in boundary leave method as 0 and continue
        if bnd_idx(i,j) == -1
            j = j + 1;
            continue;
        end
        
        % Get current index
        j_temp = j;
        
        % Determine number of spaces until a boundary or edge of image
        while j <= dimensions(2) && bnd_idx(i,j) ~= -1 
            gap = gap + 1;
            j = j + 1;
        end
        
        % Assign appropriate methods
        switch gap
            case 1
                methodsX(i, j_temp)   = 7;
            case 2
                methodsX(i, j_temp)   = 1;
                methodsX(i, j_temp+1) = 12;
            case 3
                methodsX(i, j_temp)   = 2;
                methodsX(i, j_temp+1) = 8;
                methodsX(i, j_temp+2) = 13;
            case 4
                methodsX(i, j_temp)   = 3;
                methodsX(i, j_temp+1) = 5;
                methodsX(i, j_temp+2) = 10;
                methodsX(i, j_temp+3) = 14;
            otherwise 
                methodsX(i, j_temp)   = 4;
                methodsX(i, j_temp+1) = 6;
                for k = 1:gap-4
                    methodsX(i, j_temp+1+k) = 9;
                end
                methodsX(i, j_temp+gap-4+2) = 11;
                methodsX(i, j_temp+gap-4+3) = 15; 
        end
        gap = 0;
    end
    i = i + 1;
    j = 1;
end

i = 1;
j = 1;
gap = 0;

% Main loop y direction derivatives
while j <= dimensions(2)
    while i <= dimensions(1) 
        
        % If in boundary leave method as 0 and continue
        if bnd_idx(i,j) == -1
            i = i + 1;
            continue;
        end
        
        % Get current index
        i_temp = i;
        
        % Determine number of spaces until a boundary or edge of image
        while i <= dimensions(1) && bnd_idx(i,j) ~= -1
            gap = gap + 1;
            i = i + 1;
        end
        
        % Assign appropriate methods
        switch gap
            case 1
                methodsY(i_temp, j)   = 7;
            case 2
                methodsY(i_temp, j)   = 1;
                methodsY(i_temp+1, j) = 12;
            case 3
                methodsY(i_temp, j)   = 2;
                methodsY(i_temp+1, j) = 8;
                methodsY(i_temp+2, j) = 13;
            case 4
                methodsY(i_temp, j)   = 3;
                methodsY(i_temp+1, j) = 5;
                methodsY(i_temp+2, j) = 10;
                methodsY(i_temp+3, j) = 14;
            otherwise 
                methodsY(i_temp, j)   = 4;
                methodsY(i_temp+1, j) = 6;
                for k = 1:gap-4
                    methodsY(i_temp+1+k, j) = 9;
                end
                methodsY(i_temp+gap-4+2, j) = 11;
                methodsY(i_temp+gap-4+3, j) = 15; 
        end
        gap = 0;
    end
    j = j + 1;
    i = 1;
end
end