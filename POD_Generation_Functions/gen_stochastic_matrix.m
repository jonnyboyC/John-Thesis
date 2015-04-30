function [stoc_matrix] = gen_stochastic_matrix(h, groups, direct, save_figures)
% Generate the approximate stochastic matrix from the group data

figure(h);

num_cluster = max(groups);
stoc_matrix = zeros(num_cluster);

% Add values to each entry
for i = 1:length(groups)-1;
    stoc_matrix(groups(i), groups(i+1)) = stoc_matrix(groups(i), groups(i+1)) + 1;
end

% determine probability as ratio
for i = 1:num_cluster
   stoc_matrix(i,:) = stoc_matrix(i,:)/sum(stoc_matrix(i,:));
end

% plot transition matrix
pcolor(stoc_matrix);
colorbar;
xlabel('Next State');
ylabel('Current State');

if ~isempty(save_figures)
    for i = 1:length(save_figures)
        saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
            'Clusters' filesep 'Stochastic_matirx_modes', num2str(size(groups,2))], save_figures{i});
    end
end

drawnow;