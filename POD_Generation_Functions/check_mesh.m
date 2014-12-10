function uniform= check_mesh(x, y)
% Check to ensure that the mesh is consistently spaced
mean_step_x = mean(mean(x(1:end-1,:),2) - mean(x(2:end,:),2));
mean_step_y = mean(mean(y(:,1:end-1),1) - mean(y(:,2:end),1));

std_step_x = std(mean(x(1:end-1,:),2) - mean(x(2:end,:),2));
std_step_y = std(mean(y(:,1:end-1),1) - mean(y(:,2:end),1));

if std_step_x/mean_step_x < 1e-6 && std_step_y/mean_step_y < 1e-6
    uniform = true;
else
    uniform = false;
end
end