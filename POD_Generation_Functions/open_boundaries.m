function [bnd_x, bnd_y] = open_boundaries(bnd_x, bnd_y, sides)
% TODO include more stuff
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
end