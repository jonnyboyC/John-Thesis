function stoc_matrix = gen_stochastic_matrix(num_clusters, groups, multiplier, valid)
% Generate the approximate stochastic matrix from the group data

stoc_matrix = zeros(num_clusters);

% Add values to each entry
for i = 1:multiplier
    sub_group = groups(i:multiplier:end);
    for j = 1:length(sub_group)-1
        stoc_matrix(sub_group(j), sub_group(j+1)) = ...
            stoc_matrix(sub_group(j), sub_group(j+1)) + 1;
    end
end

% determine probability as ratio
for i = 1:num_clusters
    row = stoc_matrix(i,:);
    row_col = sum(row);
    if row_col == 0;
        row_col = 1;
    end
    row = row/row_col;
    % Modifiy stochastic matrix so all transitions are valid
    if valid
        blanks = sum(row == 0);
        row = (1 - 0.001*blanks)*row;
        row(row == 0) = 0.001;
    end
    stoc_matrix(i,:) = row;
end
