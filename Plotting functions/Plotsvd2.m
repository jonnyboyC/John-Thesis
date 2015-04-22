function Plotsvd2(x, y, pod, dimensions,varname,lambda, bnd_idx, direct, save_figures)
% Plot pod modes for u v and vorticity

% Determine number of modes
[~,num_modes] = size(pod);

% Determine energy content of each mode
energy = lambda(1:num_modes)./sum(lambda)*100;

% Determine maximum values of each POD mode
cmax=max(abs(pod));

% Scale all images to max magnitude of one
for i = 1:num_modes
    pod(:,i) = pod(:,i)/cmax(i);
end

% Get color scaling
cmax = 1;
cmin = -cmax;

% set plot data
data.x = x;
data.y = y;
data.bnd_idx = bnd_idx;

% Preallocate figure subplot handles
h_sub = gobjects(4, 1);
ax_sub = gobjects(4, 1);
plot_img_num = 1;

% Generate plot handles
h = figure('Name',['  Variable: ' varname ',  (' num2str(sum(energy),4) '%)'],'color','w');
h.Position = [500, 500, 800, 500];

for i = 1:num_modes
    % Update pod data
    data.pod = reshape(pod(:,i),dimensions);
    
    % After Originally generating 4 plots, update values
    if plot_img_num > 4 
        % Save images in requested format
        if any(ismember({'fig', 'jpg'}, save_figures))
            for j = 1:size(save_figures,2)
                saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
                    'Modes' filesep 'POD_' varname '_' num2str(i-4) '_' num2str(i-1)], save_figures{j});
            end
        end
        plot_img_num = 1;
    end
    % Make plot active to avoid writing to wrong axis
    figure(h);
    subplot(2,2,plot_img_num)
    
    % plot individual plots
    if i <= 4
        [h_sub(i), ax_sub(i)] = Plottec2(data);
        ax_sub(plot_img_num).XLabel = xlabel('x/D', 'fontname','times new roman','fontsize',12);
        ax_sub(plot_img_num).YLabel = ylabel('y/D', 'fontname','times new roman','fontsize',12);
        ax_sub(plot_img_num).CLim = [cmin cmax];
        colorbar;
    else
        h_sub(plot_img_num) = Plottec2(data, h_sub(plot_img_num));
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
                    'Modes' filesep 'POD_' varname '_' num2str(i-4) '_' num2str(i-1)], save_figures{j});            
            end
        end
    end
    
end 