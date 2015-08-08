function plot_finite_diff_method(X, U)
% PLOT_FINITE_DIFF_METHOD create a plot of the finite different methods
% used for some flow field
%
%   PLOT_FINITE_DIFF_METHOD(X, U) 

% If using real data get mean
ensemble_dim = flow_dims(U);
mean_U = mean_comps(U, ensemble_dim);

% get components
x = flow_comps_ip(X);

% Find boundaries and select method
[bnd_X, bnd_idx] = refine_bounds(X, U, mean_U, pwd, false, true);
methods = select_method(bnd_idx, bnd_X, dimensions, x);

% Info for plotting
m = flow_comps(methods);
comps = flow_ncomps(methods);

% Plot methods in each direction
for k = 1:length(comps)
    figure;
    ax = gca;

    % Convert symbols to string
    symbols = finite_diff_code(methods.(m{k}));

    % Setup axis
    ax.XLim = [1 size(bnd_idx,1) + 1];
    ax.XTick = 1:size(bnd_idx,1) + 1;
    ax.YLim = [1 size(bnd_idx,2) + 1];
    ax.YTick = 1:size(bnd_idx,2) + 1;
    hold(ax, 'on');
    grid(ax, 'on');
    

    for i = 1:size(bnd_idx,1)
        for j = 1:size(bnd_idx,2)
            if bnd_idx(i,j) == -1
                patch = fill([i, i, i+1, i+1], [j, j+1, j+1, j], 'k');
                patch.FaceAlpha = 0.3;
            end
            if bnd_idx(i,j) == 0 && bnd_X.x(i,j) == 0 && bnd_X.y(i,j) == 0
                patch = fill([i, i, i+1, i+1], [j, j+1, j+1, j], 'r');
                patch.FaceAlpha = 0.15;
            end
            if bnd_idx(i,j) == 0 && (~bnd_X.x(i,j) == 0 || ~bnd_X.y(i,j) == 0)
                patch = fill([i, i, i+1, i+1], [j, j+1, j+1, j], 'k');
                patch.FaceAlpha = 0.15;
            end
            text(i+0.1,j+0.5, symbols{i,j});
        end
    end
end