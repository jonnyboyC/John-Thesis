function plot_vector_modes(data, u, v, num_modes, dimensions, varname, lambda, direct, save_figures)
% Plot vector modes of the velocity field

% Determine energy content of each mode
energy = lambda(1:num_modes)./sum(lambda)*100;

% Calculate vector magnitude 
magnitude = sqrt(u.^2 + v.^2);

% Determine maximum values of each mode
cmax=max(abs(magnitude));

% Scale all images to max magnitude of one
for i = 1:num_modes
    magnitude(:,i) = magnitude(:,i)/cmax(i);
    u(:,i) = u(:,i)/cmax(i);
    v(:,i) = v(:,i)/cmax(i);
end

% Get color scaling
cmax = 1;
cmin = 0;

% Preallocate figure subplot handles
h_surf_sub = gobjects(4, 1);
h_quiver_sub = gobjects(4, 1);
ax_sub = gobjects(4, 1);
plot_img_num = 1;

% Generate plot handles
h = figure('Name',['  Variable: ' varname ',  (' num2str(sum(energy),4) '%)'],'color','w');
h.Position = [400, 400, 1100, 650];
movegui(h,'center')


for i = 1:num_modes
    % Update pod data
    data.pod = reshape(magnitude(:,i),dimensions);
    data.u   = reshape(u(:,i), dimensions);
    data.v   = reshape(v(:,i), dimensions);
    
    % After Originally generating 4 plots, update values
    if plot_img_num > 4 
        % Save images in requested format
        if any(ismember({'fig', 'jpg'}, save_figures))
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
        [h_surf_sub(i), h_quiver_sub(i), ax_sub(i)] = plot_vector_field(data);
        ax_sub(plot_img_num).XLabel = xlabel('x/D', 'fontname','times new roman','fontsize',12);
        ax_sub(plot_img_num).YLabel = ylabel('y/D', 'fontname','times new roman','fontsize',12);
        ax_sub(plot_img_num).CLim = [cmin cmax];
        colorbar;
    else
        [h_surf_sub(plot_img_num), h_quiver_sub(plot_img_num)] ...
            = plot_vector_field(data, h_surf_sub(plot_img_num), h_quiver_sub(plot_img_num));
    end
    
    % update plot title
    ax_sub(plot_img_num).Title = title([varname ': mode ' num2str(i) ' (' num2str((energy(i)),3) ' %)'], 'fontname','times new roman','fontsize', 14);

    % update plot position
    plot_img_num=plot_img_num+1;
    
    % If last iteration save files
    if i == num_modes   
        if any(ismember({'fig', 'jpg'}, save_figures))
            for j = 1:size(save_figures,2)
                saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
                    'Modes' filesep varname '_modes' num2str(i-4) '_' num2str(i-1)], save_figures{j});            
            end
        end
    end
    
end 