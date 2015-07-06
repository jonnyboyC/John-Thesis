function stoc_matrix = gen_stochastic_matrix(num_clusters, groups)
% Generate the approximate stochastic matrix from the group data

stoc_matrix = zeros(num_clusters);

% Add values to each entry
for i = 1:length(groups)-1;
    stoc_matrix(groups(i), groups(i+1)) = stoc_matrix(groups(i), groups(i+1)) + 1;
end

% determine probability as ratio
for i = 1:num_clusters
    sum_col = sum(stoc_matrix(i,:));
    if sum_col == 0;
        sum_col = 1;
    end
    stoc_matrix(i,:) = stoc_matrix(i,:)/sum_col;
end