function [u, v, mean_u, mean_v] = filter_images(u, v, bnd_idx, bnd_x, bnd_y)
for i = 1:size(u,3)
    u_temp = u(:,:,i);
    v_temp = v(:,:,i);
    u_gauss = imgaussfilt(u_temp);
    v_gauss = imgaussfilt(v_temp);
    u_gauss = imgaussfilt(u_gauss);
    v_gauss = imgaussfilt(v_gauss);
    u_gauss(bnd_idx == 1) = u_temp(bnd_idx == 1);
    v_gauss(bnd_idx == 1) = v_temp(bnd_idx == 1);
    u(:,:,i) = imguidedfilter(u_gauss);
    v(:,:,i) = imguidedfilter(v_gauss);
end

mean_u = mean(u,3);
mean_v = mean(v,3);