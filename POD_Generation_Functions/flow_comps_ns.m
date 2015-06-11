function x = flow_comps_ns(X)
% FLOW_COMPS_NS returns the components that are not singulton i.e there are
% actually used
%
%   flow_comps = FLOW_COMPS_NS(X) 

x = flow_comps(X);
for i = length(x):-1:1;
    if all(X.(x{i}) == X.(x{i})(1))
        x(i) = [];
    end
end
