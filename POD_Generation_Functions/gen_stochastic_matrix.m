function stoc_matrix = gen_stochastic_matrix(num_clusters, groups)
% Generate the approximate stochastic matrix from the group data

stoc_matrix = zeros(num_clusters);

% Add values to each entry
for i = 1:length(groups)-1;
    stoc_matrix(groups(i), groups(i+1)) = stoc_matrix(groups(i), groups(i+1)) + 1;
end

% determine probability as ratio
for i = 1:num_clusters
    row = stoc_matrix(i,:);
    row_col = sum(row);
    if row_col == 0;
        row_col = 1;
    end
    row = row/row_col;
%     blanks = sum(row == 0);
%     row = (1 - 0.001*blanks)*row;
%     row(row == 0) = 0.001;
    stoc_matrix(i,:) = row;
end
