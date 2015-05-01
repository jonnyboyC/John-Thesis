function [stoc_matrix] = gen_stochastic_matrix(h, num_clusters, groups, direct, make_plot, save_figures)
% Generate the approximate stochastic matrix from the group data

figure(h);
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

if make_plot == true
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
end