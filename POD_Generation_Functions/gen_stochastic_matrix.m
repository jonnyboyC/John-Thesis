function [stoc_matrix] = gen_stochastic_matrix(h, groups)
% Generate the approximate stochastic matrix from the group data

figure(h);

num_cluster = max(groups);
stoc_matrix = zeros(num_cluster);

% Add values to each entry
for i = 1:length(groups)-1;
    stoc_matrix(groups(i), groups(i+1)) = stoc_matrix(groups(i), groups(i+1)) + 1;
end

for i = 1:num_cluster
   stoc_matrix(i,:) = stoc_matrix(i,:)/sum(stoc_matrix(i,:));
end
pcolor(stoc_matrix);
drawnow;