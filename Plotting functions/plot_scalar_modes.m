function plot_scalar_modes(data, U, num_modes, dimensions, non_dim, lambda, direct, save_figures)
% PLOT_SCALAR_MODES plot a set of scalar fields modes
%
% PLOT_SCALAR_MODES(data, u, v, num_modes, dimensions, varname, lambda,
% direct, save_figures, streamlines)

% Determine energy content of each mode
energy = lambda(1:num_modes)./sum(lambda)*100;

comps = flow_ncomps(U);
[x, u] = flow_comps(data.X, U);

for i = 1:comps
    % Determine maximum values of each POD mode
    cmax = max(abs(U.(u{i})));
    
    % Scale all images to max magnitude of one
    for j = 1:num_modes
        U.(u{i})(:,j) = U.(u{i})(:,j)/cmax(j);
    end
    
    % Get color scaling
    cmax = 1;
    cmin = -1;
    
    % Preallocate figure subplot handles
    h_sub = gobjects(4, 1);
    ax_sub = gobjects(4, 1);
    plot_img_num = 1;
    
    % Generate plot handles
    h = figure('Name',['  Variable: ' u{i} ',  (' num2str(sum(energy),4) '%)'],'color','w');
    h.Position = [400, 400, 1300, 600];
    movegui(h,'center')
    
    for j = 1:num_modes
        % Update pod data
        data.field = reshape(U.(u{i})(:,j),dimensions);
        
        % After Originally generating 4 plots, update values
        if plot_img_num > 4
            % Save images in requested format
            if any(ismember({'fig', 'jpg', 'png'}, save_figures))
                for k = 1:size(save_figures,2)
                    saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
                        'Modes' filesep u{i} '_modes_' num2str(j-4) '_' num2str(j-1)], save_figures{k});
                end
            end
            plot_img_num = 1;
        end
        % Make plot active to avoid writing to wrong axis
        figure(h);
        subplot(2,2,plot_img_num)
        
        % plot individual plots
        if j <= 4
            [h_sub(j), ax_sub(j)] = plot_scalar_field(data);        
            if non_dim
                ax_sub(plot_img_num).XLabel = xlabel([x{1} '/L'], 'fontname','times new roman','fontsize',14);
                ax_sub(plot_img_num).YLabel = ylabel([x{2} '/L'], 'fontname','times new roman','fontsize',14);
            else
                ax_sub(plot_img_num).XLabel = xlabel(x{1}, 'fontname','times new roman','fontsize',14);
                ax_sub(plot_img_num).YLabel = ylabel(x{2}, 'fontname','times new roman','fontsize',14);
            end
            ax_sub(plot_img_num).CLim = [cmin cmax];
            colorbar;
        else
            h_sub(plot_img_num) = plot_scalar_field(data, h_sub(plot_img_num));
        end
        
        % update plot title
        ax_sub(plot_img_num).Title = title([u{i} ': mode ' num2str(j) ' (' num2str((energy(j)),3) ' %)'], 'fontname','times new roman','fontsize', 14);
        
        % update plot position
        plot_img_num=plot_img_num+1;
        
        % If last iteration save files
        if j == num_modes
            if any(ismember({'fig', 'jpg'}, save_figures))
                for k = 1:size(save_figures,2)
                    saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
                        'Modes' filesep u{i} '_modes' num2str(j-4) '_' num2str(j-1)], save_figures{k});
                end
            end
        end
    end
end
end