function wildcard = get_wildcard(run_num, direct, data_folder, num_modes)

% Find newest file
path = [direct filesep data_folder];

% if num_modes specify move down the file tree
path = sort_modes(path, num_modes);

if ischar(run_num) && strcmp(run_num, 'first') 
        
    % short files by order
    get_sorted(path)

    matches_idx = 0;

    % Return newest file not folder
    for i = 1:length(idx)
        if ~files(i).isdir 
            if ~matches_idx(1)
                wildcard = files(i).name;
                return;
            else
            if i == matches_idx(1)
                wildcard = files(i).name;
                return;
            end
            end
        end
    end

    % if not direct match return generic wild card
    wildcard = '*.mat';
    
elseif isscalar(run_num) && run_num ~= 0
    wildcard = ['*' num2str(run_num) '*.mat'];
else
    wildcard = '*.mat';
end
end

function path = sort_modes(path, num_modes)

if num_modes ~= 0
    path = [path filesep 'modes_' num2str(num_modes)];
    if ~isdir(path)
        % exit if folder not found
       error('Galerkin coefficients for %d modes have not been produced', ...
           num_modes);
    end
end
end

function sorted_files = get_sorted(path)
% Sort files to newest
files = dir(path);
[~, idx] = sort([files.datenum], 2, 'descend');
sorted_files = files(idx);
end