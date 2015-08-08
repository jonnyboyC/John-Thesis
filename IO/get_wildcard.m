function wildcard = get_wildcard(run_num, path)
% GET_WILDCARD find specific mat file in question
%
%   wildcard = GET_WILDCARD(run_num, path) if 'first' will select newest

% If first is passed get newest
if ischar(run_num) && strcmp(run_num, 'first') 
        
    % sort files by order
    files = dir(path);
    [~, idx] = sort([files.datenum], 2, 'descend');
    files = files(idx);

    % Return newest file not folder
    for i = 1:length(files)
        if ~files(i).isdir 
            wildcard = files(i).name;
            return;
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