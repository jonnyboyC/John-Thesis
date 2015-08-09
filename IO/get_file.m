function file = get_file(run_num, path)
% GET_WILDCARD find specific mat file in question
%
%   wildcard = GET_WILDCARD(run_num, path) if 'first' will select newest

files = dir(path);

% If first is passed get newest
if ischar(run_num) && strcmp(run_num, 'first') 
        
    % sort files by order
    [~, idx] = sort([files.datenum], 2, 'descend');
    files = files(idx);

    % Return newest file not folder 
    for i = 1:length(files)
        if ~files(i).isdir 
            file = files(i).name;
            return;
        end
    end

    % if not direct match return generic wild card
    expr = '*.mat';
elseif isscalar(run_num) && run_num ~= 0
    expr = ['\w*' num2str(run_num) '\w*'];    
else
    expr = '\w*';
end

% If not match or didn't use "first" then match with generic wildcard
match = regexp({files.name}, expr);

for i = 1:length(files)
    if ~isempty(match{i}) && ~files(i).isdir
        file = files(i).name;
        return;
    end
end

file = 'not found';

end