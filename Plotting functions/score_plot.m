function [h, R2, R, p] = score_plot(full_name, sub_name, score, var, xaxis, yaxis, ...
    top_title, direct, num_clusters, save_name, scale)

if length(full_name) <= 1
    R2 = 0;
    R = 0;
    p = 0;
    return;
end

% Set up figure
h = figure;
ax = newplot;
hold(ax, 'on');
grid(ax, 'on');
models = {'GM_', 'GM1', 'GM2', 'GM3'};
model_color = {'m', 'k', 'g', 'r'};
sub_models = {'NL_', 'til_', ''};
sub_shape = {'x', 'o', 's'};
$ 
try 
for i = 1:length(models)
    for j = 1:length(sub_models)
        match1 = strncmp(full_name, models{i}, 3);
        if j ~= 3
            match2 = strncmp(sub_name, sub_models{j}, 3);
        else
            match2 = ~(strncmp(sub_name, sub_models{1}, 3) | ...
                       strncmp(sub_name, sub_models{2}, 3));
        end
        if any(match1 & match2)
            lax = plot(ax, score(match1 & match2), (var(match1 & match2)), ...
                [model_color{i}, sub_shape{j}]);
            lax.MarkerSize = 12;
        end
    end
end
catch 
    R2 = 0;
    R = 0;
    p = 0;
    return;
end 

direct_ext = [direct filesep 'Figures' filesep 'Scores' filesep 'clusters_' ...
    num2str(num_clusters)];

if ~exist(direct_ext, 'dir') 
    mkdir(direct_ext);
end

% Linear regression fit
fitted_model = fitlm(score, var, 'linear');
fh = fitted_model.plot;
delete(fh(1));

[R, p] = corrcoef(score, var);
R2 = fitted_model.Rsquared.Ordinary;

lh = legend;
lh.String{1} = [lh.String{1} '   R2: ' num2str(R2,3)];

% Format plot
ax.YScale = scale;
ax.XLabel = xlabel(xaxis, 'Interpreter', 'tex');
ax.YLabel = ylabel(yaxis, 'Interpreter', 'tex');
ax.Title = title([top_title ', Correlation: ' num2str(R(2,1))], 'Interpreter', 'tex');
ax.Title.FontWeight = 'normal';
ax.FontSize = 12;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Title.FontSize = 18;

% Save figure
file_name = [direct_ext filesep save_name];

% Save figure in Figure\Galerkin folder
saveas(h, file_name, 'fig');
saveas(h, file_name, 'png');

end