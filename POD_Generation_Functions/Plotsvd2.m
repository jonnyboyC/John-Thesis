function Plotsvd2(data,pod,sz,varname,sigma, direct, save_figures)
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

% Preallocate figure plot handles
h_sub = gobjects(4, 1);
ax_sub = gobjects(4, 1);

h = figure('Name',['  Variable: ' varname ',  (' num2str(sum(energy),4) '%)'],'color','w');
for i = 1:num_modes
    data.pod = reshape(pod(:,i),sz(1),sz(2));
    
    % After 4 plots create a new figure
    if plot_img_num > 4 
        % Save images in requested format
        plot_img_num = 1;
        if strcmp(save_figures, 'fig') || strcmp(save_figures, 'jpg') 
           saveas(h, [direct '\Figures\POD\POD_' varname '_' num2str(i-4) '_' num2str(i-1)], save_figures);
        end
        h.Name = ['  Variable: ' varname ',  (' num2str(sum(energy),4) '%)'];
    end
    
    % plot individual plots
    subplot(2,2,plot_img_num)
    if i <= 4
        [h_sub(i), ax_sub(i)] = Plottec2(data);
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
        if strcmp(save_figures, 'fig') || strcmp(save_figures, 'jpg') 
           saveas(h, [direct '\Figures\POD\POD_' varname '_' num2str(i-4) '_' num2str(i-1)], save_figures);
        end
    end
    
end 
hold off