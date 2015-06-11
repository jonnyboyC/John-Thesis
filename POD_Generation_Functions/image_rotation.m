function [X, U] = image_rotation(X, U, flow_flip)
%% Rotation PIV images, if data is originally presented in a flipped format

% Get fields
[x, u] = flow_comps(X, U);

% Flip all variables request, assums flip if a cell of strings with strings
% such as 'x', 'u' etc
for i = 1:length(flow_flip)
    if isfield(X, flow_flip{i})
        X.(flow_flip{i}) = -X.(flow_flip{i});
    end
    if isfield(U, flow_flip{i})
        U.(flow_flip{i}) = -U.(flow_flip{i});
    end
    
    % if it is not found in either display to prompt
    if ~isfield(X, flow_flip{i}) && ~isfield(U, flow_flip{i})
        fprintf('%s was not found in data, check to make sure it way typed correctly', flow_flip{i});
    end
end

% TODO need to confirm this works

% perform flips
for i = 1:flow_ncomps(X)
    if X.(x{i})(1) > X.(x{i})(end)
        % Flip all components about inversed component
        for j = 1:flow_ncomps(X)
            X.(x{j}) = flip(X.(x{j}), i);
            U.(u{j}) = flip(U.(u{j}), i);
        end
    end
end

end