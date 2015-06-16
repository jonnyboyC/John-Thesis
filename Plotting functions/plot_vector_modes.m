function plot_vector_modes(data, U, num_modes, dimensions, nondim, varname, lambda, direct, save_figures, streamlines)
% PLOT_VECTOR_MODES plot a set of vector fields modes
%
% PLOT_VECTOR_MODES(data, u, v, num_modes, dimensions, varname, lambda,
% direct, save_figures, streamlines)

% Determine energy content of each mode
energy = lambda(1:num_modes)./sum(lambda)*100;

% get the components that are in the active plane
dims = flow_dims(data.X);
[x, u] = flow_comps_ip(data.X, U);

% Calculate vector magnitude
magnitude = sqrt(U.(u{1}).^2 + U.(u{2}).^2);

% Determine maximum values of each mode
cmax = max(abs(magnitude));

% Scale all images to max magnitude of one
for i = 1:num_modes
    magnitude(:,i) = magnitude(:,i)/cmax(i);
    U.(u{1})(:,i) = U.(u{1})(:,i)/cmax(i);
    U.(u{2})(:,i) = U.(u{2})(:,i)/cmax(i);
end

% Get color scaling
cmax = 1;
cmin = 0;

% Preallocate figure subplot handles
h_mag_sub = gobjects(4, 1);
if streamlines
    h_dir_sub = cell(4, 1);
else
    h_dir_sub = gobjects(4,1);
end

ax_sub = gobjects(4, 1);
plot_img_num = 1;

% Generate plot handles
h = figure('Name',['  Variable: ' varname ',  (' num2str(sum(energy),4) '%)'],'color','w');
h.Position = [400, 400, 1100, 650];
movegui(h,'center')


for i = 1:num_modes
    % Update pod data
    data.pod = reshape(magnitude(:,i),dimensions);
    for j = 1:dims 
        data.U.(u{j}) = reshape(U.(u{j})(:,i), dimensions);
    end
    
    % After Originally generating 4 plots, update values
    if plot_img_num > 4
        % Save images in requested format
        if any(ismember({'fig', 'jpg', 'png'}, save_figures))
            for j = 1:size(save_figures,2)
                saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
                    'Modes' filesep varname '_modes_' num2str(i-4) '_' num2str(i-1)], save_figures{j});
            end
        end
        plot_img_num = 1;
    end
    % Make plot active to avoid writing to wrong axis
    figure(h);
    subplot(2,2,plot_img_num)
    
    % plot individual plots
    if i <= 4
        [h_mag_sub(i), h_dir_sub(i), ax_sub(i)] = plot_vector_field(data, streamlines);
        if nondim
            ax_sub(plot_img_num).XLabel = xlabel([x{1} '/L'], 'fontname','times new roman','fontsize',14);
            ax_sub(plot_img_num).YLabel = ylabel([x{2} '/L'], 'fontname','times new roman','fontsize',14);
        else
            ax_sub(plot_img_num).XLabel = xlabel(x{1}, 'fontname','times new roman','fontsize',14);
            ax_sub(plot_img_num).YLabel = ylabel(x{2}, 'fontname','times new roman','fontsize',14);
        end
        ax_sub(plot_img_num).CLim = [cmin cmax];
        cax = colorbar('peer', ax_sub(plot_img_num));
    else
        [h_mag_sub(plot_img_num), h_dir_sub(plot_img_num)] ...
            = plot_vector_field(data, streamlines, h_mag_sub(plot_img_num), h_dir_sub(plot_img_num));
    end
    
    % update plot title
    ax_sub(plot_img_num).Title = title([varname ': mode ' num2str(i) ' (' num2str((energy(i)),3) ' %)'], 'fontname','times new roman','fontsize', 16);
    
    % update plot position
    plot_img_num = plot_img_num + 1;
    
    % If last iteration save files
    if i == num_modes
        if any(ismember({'fig', 'jpg', 'png'}, save_figures))
            for j = 1:size(save_figures,2)
                saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
                    'Modes' filesep varname '_modes' num2str(j-4) '_' num2str(j-1)], save_figures{j});
            end
        end
    end
end
end