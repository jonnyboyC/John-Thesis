function stoc_matrix = gen_stochastic_matrix(num_clusters, groups, varargin)
% Generate the approximate stochastic matrix from the group data

extras = length(varargin);

% Fill out optional parameters
switch extras
    case 1
        multiplier = varargin{1};
    case 2
        multiplier = varargin{1};
        classify = varargin{2};
    case 3
        multiplier = varargin{1};
        classify = varargin{2};
        outlier_mode = varargin{3};
end

% If we are using a outlier mode add extra row/column
if ~isempty(outlier_mode.km)
    stoc_matrix = zeros(num_clusters + 1);
    num_clusters = num_clusters + 1;
else
    stoc_matrix = zeros(num_clusters);
end

% Add values to each entry
for i = 1:multiplier
    % Break up into indivdual chains
    sub_group = groups(i:multiplier:end);
    
    for j = 1:length(sub_group) - 1
        % Add elements
        stoc_matrix(sub_group(j), sub_group(j+1)) = ...
            stoc_matrix(sub_group(j), sub_group(j+1)) + 1;
    end
end

% Generate MLE of the stochastic matrix
for i = 1:num_clusters
    row = stoc_matrix(i,:);
    row_col = sum(row);
    
    if row_col == 0;
        row_col = 1;
    end
    row = row/row_col;
    
    % Modifiy stochastic matrix so all transitions are classify
    if classify
        blanks = sum(row == 0);
        row = (1 - 0.001*blanks)*row;
        row(row == 0) = 0.001;
    end
    stoc_matrix(i,:) = row;
end
