%% This script simply changes the names and saves saves the variables in a 
% .mat file in order to be used by Gal_Proj
mean_u = muz{1};
mean_v = mvz{1};
pod_u1 = phiutc{1};
pod_v1 = phivtc{1};
vol_frac = vt{1};
bnd_idx = bi;
lambda2 = sgtc2;
eig_func = psitc;
x = x{1};
y = y{1};
dimensions = sa{1};
save('D:\shear layer\PIVData\Old Data\POD Data\POD.mat', 'mean_u', 'mean_v',...
    'pod_u1', 'pod_v1', 'vol_frac', 'bnd_idx', 'lambda2', 'x', 'y', 'dimensions', ...
    'eig_func');
