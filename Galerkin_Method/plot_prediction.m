function plot_prediction(pod_u, pod_v, x, y, modal_amp, num_pods, dimensions, direct)

% Begin fill all data structures
data_u.pod = [];
data_u.xg = x;
data_u.yg = y;

data_v.pod = [];
data_v.xg = x;
data_v.yg = y;

data_m.pod = [];
data_m.xg = x;
data_m.yg = y;

% Fill all plots with blank images to set renderer to opengl
h = figure(1);
dummie = zeros(2,2);
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


% Preallocate size for predicated flow
data_u.pod = zeros(dimensions(1), dimensions(2), length(modal_amp(:,1)));
data_v.pod = zeros(dimensions(1), dimensions(2), length(modal_amp(:,1)));
data_m.pod = zeros(dimensions(1), dimensions(2), length(modal_amp(:,1)));

% Calculate predication by summing modes
for i = 1:length(modal_amp(:,1))
    for j = 1:num_pods
        data_u.pod(:,:,i) =  data_u.pod(:,:,i)...
            + reshape(pod_u(:,j)*modal_amp(i,j),dimensions(1), dimensions(2));
        data_v.pod(:,:,i) =  data_v.pod(:,:,i)...
            + reshape(pod_v(:,j)*modal_amp(i,j),dimensions(1), dimensions(2));
    end
    data_m.pod(:,:,i) = sqrt(data_u.pod(:,:,i).^2+data_v.pod(:,:,i).^2);
end

% Intialize Video creator
writer = VideoWriter([direct '\Figures\Movies\POD' num2str(num_pods) '_Galerkin.avi']);
open(writer);

%% TODO May need to relook at this to make it more memory efficient
data_temp.xg = x;
data_temp.yg = y;
data_temp.pod = zeros(dimensions(1), dimensions(2));
data = {data_u, data_v, data_m};

% Determine max values for u v and magnitude
cmax = zeros(3,1);
cmin = zeros(3,1);
for i = 1:length(data)
    cmax(i) = max(max(max(abs(data{i}.pod))));
    cmin(i) = -cmax(i);
end

% Preallocated figure handles and axes handles
h_sub = gobjects(3,1);
ax_sub = gobjects(3,1);

% Plot results, print current image number, and save images to .avi video
for i = 1:length(modal_amp(:,1))
    fprintf('image %d of %d\n', i, length(modal_amp(:,1)));
    for j = 1:length(data);
        data_temp.pod = squeeze(data{j}.pod(:,:,i));
        subplot(3,1,j);
        if i == 1
            [h_sub(j), ax_sub(j)] = Plottec2(data_temp);
            colorbar;
        else
            h_sub(j) = Plottec2(data_temp, h_sub(j));
        end
        ax_sub(j).ZLim = [cmin(j) cmax(j)];
        ax_sub(j).CLim = [cmin(j) cmax(j)];
    end
    frame = getframe(gcf);
    writeVideo(writer, frame)
end
end