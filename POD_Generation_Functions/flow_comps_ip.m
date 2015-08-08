function [x_ip, u_ip] = flow_comps_ip(X, U)
% FLOW_COMPS_IP returns the velocity vectors that are in the plane of the
% non-singleton mesh dimensions assumes that x -> u, y -> v, z -> w
%
%   x = FLOW_COMPS_IP(X)
%   [x, u] = FLOW_COMPS_IP(X, U)

if nargin == 1 && nargout == 1
    x_ip = flow_comps_ns(X);
else
    x = flow_comps_ns(X);
    u = flow_comps(U);

    u_ip = cell(3,1);
    if any(strcmp('x', x)) && any(strcmp('u', u));
        u_ip{1} = 'u';
    end 
    if any(strcmp('y', x)) && any(strcmp('v', u));
        u_ip{2} = 'v';
    end
    if any(strcmp('z', x)) && any(strcmp('w', u));
        u_ip{3} = 'w';
    end

    x_ip = x;
    u_ip(cellfun('isempty', u_ip)) = [];
end
end