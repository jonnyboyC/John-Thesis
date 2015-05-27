function theta = momentum_thickness(mean_u, low, scale, bnd_idx)
dimensions = size(bnd_idx);
theta = zeros(dimensions(1),1);

i = 1;
j = 1;
gap = 0;
high = 1-low;

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
        
        uy = mean_u(i,j_temp:j_temp+gap-1)-low;
        theta(i) = sum(uy.*(high-uy)*scale) + theta(i);
        gap = 0;
    end
    i = i + 1;
    j = 1;        
end
end