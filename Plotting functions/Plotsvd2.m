function Plotsvd2(data,pod,dimensions,varname,sigma, bnd_idx, direct, save_figures)
% Plot plot pod modes for u v and vorticity

% Determine number of modes
[~,num_modes] = size(pod);

% Determine energy content of each mode
energy = sigma(1:num_modes)./sum(sigma)*100;

plot_img_num = 1;
cma=max(abs(pod));

% Scale all images to max magnitude of one
for i = 1:num_modes
    pod(:,i)=pod(:,i)/cma(i);
end

% Believe this is passed to scale colors
cma=1;
cmi=-cma;
data.cmax = cma;

% Preallocate figure plot handles
h_sub = gobjects(4, 1);
ax_sub = gobjects(4, 1);

h = figure('Name',['  Variable: ' varname ',  (' num2str(sum(energy),4) '%)'],'color','w');
h.Position = [500, 500, 800, 400];
for i = 1:num_modes
    data.pod = reshape(pod(:,i),dimensions(1),dimensions(2));
    
    % After 4 plots create a new figure
    if plot_img_num > 4 
        % Save images in requested format
        plot_img_num = 1;
        if any(ismember({'fig', 'jpg'}, save_figures))
            for j = 1:size(save_figures,2)
                saveas(h, [direct '\Figures\POD\POD_' varname '_' num2str(i-4) '_' num2str(i-1)], save_figures{j});
            end
        end
        h.Name = ['  Variable: ' varname ',  (' num2str(sum(energy),4) '%)'];
    end
    
    subplot(2,2,plot_img_num)
    % plot individual plots
    if i <= 4
        [h_sub(i), ax_sub(i)] = Plottec2(data, 0, bnd_idx);
        colorbar;
    else
        h_sub(plot_img_num) = Plottec2(data,h_sub(plot_img_num));
    end
    % update plot quantities
    ax_sub(plot_img_num).Title = title([varname ', m = ' num2str(i) ', (' num2str((energy(i)),3) ' %)'], 'fontname','times new roman','fontsize', 14);
    ax_sub(plot_img_num).XLabel = xlabel('x/D', 'fontname','times new roman','fontsize',12);
    ax_sub(plot_img_num).YLabel = ylabel('y/D', 'fontname','times new roman','fontsize',12);
    ax_sub(plot_img_num).CLim = [cmi cma];

    % update plot position
    plot_img_num=plot_img_num+1;
    
    % If last iteration save files
    if i == num_modes   
        if any(ismember({'fig', 'jpg'}, save_figures))
            for j = 1:size(save_figures,2)
                saveas(h, [direct '\Figures\POD\POD_' varname '_' num2str(i-4) '_' num2str(i-1)], save_figures{j});
            end
        end
    end
    
end 
hold off