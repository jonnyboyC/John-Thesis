function [xn, yn, un, vn] = image_rotation(x, y, u, v, flip_x, flip_y)
%% Rotation PIV images, if data is originally presented in a flipped format

% prefill for speed
xn = zeros(size(x));
yn = zeros(size(y));
un = zeros(size(u));
vn = zeros(size(v));

% TODO come up with more elegant solution
if all(x <= 0)
    x = x*-1;
end

% create vector of indexes to avoid loops
x_idx = 1:size(x,1);
y_idx = 1:size(y,2);

% Flip x and y if requested
if flip_x
    x = -x; 
    u = -u;
end
if flip_y
    y = -y;
    v = -y;
end

% perform flips
if x(1,1) > x(end,1) && y(1,1) > y(1,end)
    xn(x_idx, y_idx) = x(size(x,1) - x_idx + 1, size(x,2) - y_idx + 1);
    yn(x_idx, y_idx) = y(size(y,1) - x_idx + 1, size(y,2) - y_idx + 1);
    un(x_idx, y_idx) = u(size(u,1) - x_idx + 1, size(u,2) - y_idx + 1);
    vn(x_idx, y_idx) = v(size(v,1) - x_idx + 1, size(v,2) - y_idx + 1);
elseif x(1,1) > x(end,1) && y(1,1) < y(1, end)
    xn(x_idx,:) = x(size(x,1)-x_idx+1,:);
    yn(x_idx,:) = y(size(y,1)-x_idx+1,:);
    un(x_idx,:) = u(size(u,1)-x_idx+1,:);
    vn(x_idx,:) = v(size(v,1)-x_idx+1,:);
elseif x(1,1) < x(end,1) && y(1,1) > y(1,end)
    xn(:,y_idx) = x(:,size(x,2)-y_idx+1);
    yn(:,y_idx) = y(:,size(y,2)-y_idx+1);
    un(:,y_idx) = u(:,size(u,2)-y_idx+1);
    vn(:,y_idx) = v(:,size(v,2)-y_idx+1);
end
end