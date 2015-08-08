function path = get_mode_folder(path, num_modes, custom)
% GET_MODE_FOLDER given a path to the data directory of choise select the
% subfolder coresponding tot he correct number of modes or custom mode
% combination
%
%   path = get_mode_folder(path, num_modes, custom)


% Get files in folder
files = dir(path);

% Find matches to pattern
if custom
    match = regexp({files.name}, ['^c\w*' num2str(num_modes) ]);
else
    match = regexp({files.name}, ['^m\w*' num2str(num_modes) ]);
end

% Report if not match is found
if all(cellfun(@isempty, match))
    error('Galerkin coefficients for %d modes have not been produced', ...
        num_modes);
end

% Ensure our pattern match goes to a folder
for i = 1:length(files)
    if ~isempty(match{i}) && files(i).isdir
        break;
    end
end

% Extend and return path
path = [path filesep files(i).name];
end