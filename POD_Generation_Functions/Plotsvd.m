function Plotsvd(data,pod,sz,varname,sigma, direct, save_figures)


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

h = figure('Name',['  Variable: ' varname ',  (' num2str(sum(energy),4) '%)'],'color','w');
for i = 1:num_modes
    data.pod = reshape(pod(:,i),sz(1),sz(2));
    
    % After 4 plots create a new figure
    if plot_img_num > 4
        % May want to include possibility of saving images  
        plot_img_num = 1;
        if strcmp(save_figures, 'fig') || strcmp(save_figures, 'jpg') 
           saveas(h, [direct '\Figures\POD\POD_' varname '_' num2str(i-4) '_' num2str(i-1)], save_figures);
        end
        h = figure('Name',['  Variable: ' varname ',  (' num2str(sum(energy),4) '%)'],'color','w');
    end
    
    % plot indivdual plots
    subplot(2,2,plot_img_num)   
    Plottec(data);
    title([varname ', m = ' num2str(i) ', (' num2str((energy(i)),3) ' %)'],'fontname','times new roman','fontsize',14);%,round'FontWeight','demi'
    xlabel('x/D','fontname','times new roman','fontsize',12);
    ylabel('y/D','fontname','times new roman','fontsize',12);
    hold on 
    plot_img_num=plot_img_num+1;
    set(gca,'clim',[cmi cma])
    colorbar
end 
hold off