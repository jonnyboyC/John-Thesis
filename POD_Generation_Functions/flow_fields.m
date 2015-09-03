function fields = flow_fields(X)
% FLOW_FIELDS return the number of total fields in a structure
%
% fields = FLOW_FIELDS(X)

fields = 0;

% Check that a structure was passed in
if isstruct(X) 
    fields = step_in(X);
end
end

function fields = step_in(X)
% STEP_IN Recursively step into a structure to count nested struture fields
%
%   fields = STEP_IN(X)

if isstruct(X)
    x = flow_comps(X);
    count = 0;
    for i = 1:length(x)
        count = count + step_in(X.(x{i}));
    end
    fields = count;
else
    fields = 1;
end

end