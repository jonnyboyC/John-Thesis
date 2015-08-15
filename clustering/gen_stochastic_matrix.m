function model = gen_stochastic_matrix(model, num_clusters, varargin)
% Generate the approximate stochastic matrix from the group data

extras = length(varargin);
groups = model.groups;

% Fill out optional parameters
switch extras
    case 1
        multiplier = varargin{1};
    case 2
        multiplier = varargin{1};
        classify = varargin{2};
end

% Add an extra outlier mode
model.stoch = zeros(num_clusters + 1);
num_clusters = num_clusters + 1;

% Add values to each entry
for i = 1:multiplier
    % Break up into indivdual chains
    sub_group = groups(i:multiplier:end);
    
    for j = 1:length(sub_group) - 1
        % Add elements
        model.stoch(sub_group(j), sub_group(j+1)) = ...
            model.stoch(sub_group(j), sub_group(j+1)) + 1;
    end
end

% Generate MLE of the stochastic matrix
for i = 1:num_clusters
    row = model.stoch(i,:);
    row_col = sum(row);
    
    if row_col == 0;
        row_col = 1;
    end
    row = row/row_col;
    model.stoch(i,:) = row;
end

% Before any modifcation determine if model is transient
model = stationary(model);
if any(model.stat > .99)
    model.steady = false;
else
    model.steady = true;
end

% Modifiy stochastic matrix so all transitions are valid if we are
% classifying
if classify
    for i = 1:num_clusters
        row = model.stoch(i,:);
        blanks = sum(row == 0);
        row = (1 - 0.001*blanks)*row;
        row(row == 0) = 0.001;
        model.stoch(i,:) = row;
    end
end
