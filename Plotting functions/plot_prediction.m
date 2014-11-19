function plot_prediction(pod_u, pod_v, x, y, modal_amp, t, num_pods, dimensions, direct, MOD)

% Check to match sure requested image instead too large
if size(modal_amp,1) > 10000
    disp(['Size of prediction is too large until further work is preformed' ...
        'on the plot_prediction function']);
    return;
end

% TODO add vorticity

% Begin fill all data structures
data_u.pod = [];
data_u.xg = x;
data_u.yg = y;

data_v.pod = [];
data_v.xg = x;
data_v.yg = y;

% data_vor.pod = [];
% data_vor.xg = x;
% data_vor.yg = y;

data_m.pod = [];
data_m.xg = x;
data_m.yg = y;

% Fill all plots with blank images to set renderer to opengl
dummie = zeros(2,2);
h = figure;
subplot(3,1,1)
pcolor(dummie);
axis tight
set(gca, 'nextplot', 'replacechildren');
set(h, 'Renderer', 'opengl');

subplot(3,1,2)
pcolor(dummie);
axis tight
set(gca, 'nextplot', 'replacechildren');
set(h, 'Renderer', 'opengl');

subplot(3,1,3)
pcolor(dummie);
axis tight
set(gca, 'nextplot', 'replacechildren');
set(h, 'Renderer', 'opengl');

Hz = 1/(t(2) - t(1));

% Intialize Video creator
ext = '\Figures\Movies\POD_';
if nargin == 10
    ext = [ext MOD];
end
writer = VideoWriter([direct ext num2str(num_pods) '_' ...
   num2str(t(1)) '_' num2str(t(end)) 's_' num2str(ceil(Hz)) 'Hz_Galerkin.avi']);
writer.Quality = 100;
writer.FrameRate = 60;
open(writer);

%% TODO May need to relook at this to make it more memory efficient
data_temp.xg = x;
data_temp.yg = y;
data_temp.pod = zeros(dimensions(1), dimensions(2));

data_u = zeros(dimensions(1), dimensions(2), size(modal_amp,1));
data_v = zeros(dimensions(1), dimensions(2), size(modal_amp,1));
data_m = zeros(dimensions(1), dimensions(2), size(modal_amp,1));

% calculate predicted images
for i = 1:size(modal_amp,1)
    for j = 1:num_pods
        data_u(:,:,i) =  data_u(:,:,i)...
            + reshape(pod_u(:,j)*modal_amp(i,j),dimensions(1), dimensions(2));
        data_v(:,:,i) =  data_v(:,:,i)...
            + reshape(pod_v(:,j)*modal_amp(i,j),dimensions(1), dimensions(2));
    end
    data_m(:,:,i) = sqrt(data_u(:,:,i).^2+data_v(:,:,i).^2);
end
data = {data_u, data_v, data_m};

% Determine max values for u v and magnitude
cmax = zeros(3,1);
cmin = zeros(3,1);
for i = 1:length(data)
    cmax(i) = max(data{i}(:));
    cmin(i) = -cmax(i);
end

% Preallocated figure handles and axes handles
h_sub = gobjects(3,1);
ax_sub = gobjects(3,1);

% Plot results, print current image number, and save images to .avi video
for i = 1:length(modal_amp(:,1))  
    fprintf('image %d of %d\n', i, length(modal_amp(:,1)));
    for j = 1:length(data);
        data_temp.pod = squeeze(data{j}(:,:,i));
        data_temp.cmax = cmax(j);
        subplot(3,1,j);
        if i == 1
            [h_sub(j), ax_sub(j)] = Plottec2(data_temp);
            colorbar;
        else
            h_sub(j) = Plottec2(data_temp, h_sub(j));
        end
        ax_sub(j).ZLim = [cmin(j), cmax(j)];
        ax_sub(j).CLim = [cmin(j), cmax(j)];
    end
    frame = getframe(gcf);
    writeVideo(writer, frame)
end
end