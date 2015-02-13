function [bnd_x, bnd_y] = open_boundaries(bnd_x, bnd_y, mean_u, sides, mask)
% Alter the automatic boundary detection 

% Perform base masking based on which sides are specified
if ~ismember(sides, 'left')
    bnd_x(bnd_x > 0) = 0;
end
if ~ismember(sides, 'right')
    bnd_x(bnd_x < 0) = 0;
end
if ~ismember(sides, 'top')
    bnd_y(bnd_y > 0) = 0;
end
if ~ismember(sides, 'bottom')
    bnd_y(bnd_y < 0) = 0;
end

% check if a mask is requested
if mask == true
    mask_mat = generate_mask(bnd_x, bnd_y, mean_u);
else
    mask_mat = false(size(mean_u));
end

% Mask anything that has been requested 
bnd_x(mask_mat) = 0;
bnd_y(mask_mat) = 0;
end

function mask_mat = generate_mask(bnd_x, bnd_y, mean_u)


end

