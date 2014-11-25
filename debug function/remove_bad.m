load('E:\thesis\PIVData\Forced\Processed Data\Processed.mat');
mean_u = mean(ui, 3);
mean_v = mean(vi, 3);

error_u = zeros(size(ui,3), 1);
error_v = zeros(size(vi,3), 1);
for i = 1:size(ui,3)
     error_u(i) = sum(sum((mean_u-ui(:,:,i)).^2));
     error_v(i) = sum(sum((mean_v-vi(:,:,i)).^2));
end

figure(10);
plot(error_u);
figure(11);
plot(error_v);