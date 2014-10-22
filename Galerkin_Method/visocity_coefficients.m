function [l_dot, l, q_2dot, q_dot, q] = visocity_coefficients(mean_u, ...
    mean_v, x, y, pod_u, pod_v, dimensions, vol_frac, bnd_idx, z)
   
%% TODO understand what these do
at = 0;
clt = 0;
cqt = 0;
cbt = 0;
cct = 0;
cdt = 0;

num_elem = numel(x);
num_modes = size(pod_u, 2);

% TODO figure out what is really being calculated here
[xxi, yxi, xet, yet, aj] = metric(x, y);

% Calculate coefficients for u's & v's derivatives
[udx, ud2x, udy, ud2y] = derivatives(mean_u, dimensions, z, xxi, yxi,...
    xet, yet, aj, bnd_idx);
[vdx, vd2x, vdy, vd2y] = derivatives(mean_v, dimensions, z, xxi, yxi,...
    xet, yet, aj, bnd_idx);

% Calculate coefficients for for pod_u's & pod_v's derivatives
[pod_udx, pod_ud2x, pod_udy, pod_ud2y] = derivatives(pod_u, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);
[pod_vdx, pod_vd2x, pod_vdy, pod_vd2y] = derivatives(pod_v, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);

% Convert all matrix quantities to vectors
udx     = reshape(udx, num_elem, 1);
udy     = reshape(udy, num_elem, 1);
ud2x    = reshape(ud2x, num_elem, 1);
ud2y    = reshape(ud2y, num_elem, 1);
vdx     = reshape(vdx, num_elem, 1);
vdy     = reshape(vdy, num_elem, 1);
vd2x    = reshape(vd2x, num_elem, 1);
vd2y    = reshape(vd2y, num_elem, 1);

% Convert all matrix quantities to vectors
pod_udx     = reshape(pod_udx, num_elem, num_modes);
pod_udy     = reshape(pod_udy, num_elem, num_modes);
pod_ud2x    = reshape(pod_ud2x, num_elem, num_modes);
pod_ud2y    = reshape(pod_ud2y, num_elem, num_modes);
pod_vdx     = reshape(pod_vdx, num_elem, num_modes);
pod_vdy     = reshape(pod_vdy, num_elem, num_modes);
pod_vd2x    = reshape(pod_vd2x, num_elem, num_modes);
pod_vd2y    = reshape(pod_vd2y, num_elem, num_modes);

% Need to check this out, should it be distance formula?
d2u = (ud2x + ud2y);
d2v = (vd2x + vd2y);
d2pod_u = (pod_ud2x + pod_ud2y);
d2pod_v = (pod_vd2x + pod_vd2y);

% clear ud2x ud2y vd2x vd2y pod_ud2x pod_ud2y pod_vd2x pod_vd2y

% Find inner produce with itself for pod modes 
% TODO find if at above is the actual inner product, currently allowing for
% at to have some sort of linear offset. Fudge factor?
au = inner_prod(pod_u, pod_u, vol_frac);
av = inner_prod(pod_v, pod_v, vol_frac);
at = at + au + av; 

% Find inner product between 2nd for mean velocity and pod modes
clu = inner_prod(d2u, pod_u, vol_frac);
clv = inner_prod(d2v, pod_v, vol_frac);
clt = clt + clu + clv;

% Find inner product between 2nd derivative of pod modes and pod modes
cbu = inner_prod(d2pod_u, pod_u, vol_frac);
cbv = inner_prod(d2pod_v, pod_v, vol_frac);
cbt = cbt + cbu + cbv;

% TODO Look at Navier Stokes Equation
qu = mean_u.*udx+mean_v.*udy;
qv = mean_u.*vdx+mean_v.*vdy;

% Inner product of pod modes with calculated quantity. Similar situation
% with cqt we are allowing for a linear offset, may be fudge factor
cqu = inner_prod(qu, pod_u, vol_frac);
cqv = inner_prod(qv, pod_v, vol_frac);
cqt = cqt + cqu + cqv;

% TODO look for better names for all these variables 
% Linear terms
% Group 1
u_pod_ux = mean_u*ones(1,num_modes).*pod_udx;
u_pod_vx = mean_u*ones(1,num_modes).*pod_vdx;

% Group 2
v_pod_uy = mean_v*ones(1,num_modes).*pod_udy;
v_pod_vy = mean_v*ones(1,num_modes).*pod_vdy;

% Group 3
ux_pod_u = pod_u.*(udx*ones(1,num_modes));
vx_pod_u = pod_u.*(vdx*ones(1,num_modes));

% Group 4
uy_pod_v = pod_v.*(udy*ones(1,num_modes));
vy_pod_v = pod_v.*(vdy*ones(1,num_modes));

% Sum of terms 
cu = u_pod_ux + v_pod_uy + ux_pod_u + uy_pod_v;
cv = u_pod_vx + v_pod_vy + vx_pod_u + vy_pod_v;

ccu = inner_prod(cu, pod_u, vol_frac);
ccv = inner_prod(cv, pod_v, vol_frac);
cct = cct + ccu + ccv;

% Quadractic Terms
% Appears to create 

pod_u_pod_u_x = zeros(size(pod_u,1), num_modes, num_modes);
pod_v_pod_u_y = zeros(size(pod_u,1), num_modes, num_modes);
pod_u_pod_v_x = zeros(size(pod_u,1), num_modes, num_modes);
pod_v_pod_v_y = zeros(size(pod_u,1), num_modes, num_modes);
cdu = zeros(num_modes, num_modes, num_modes);
cdv = zeros(num_modes, num_modes, num_modes);


for k = 1:num_modes
    pod_u_pod_u_x(:,:,k) = (pod_u(:,k)*ones(1,num_modes)).*pod_udx;
    pod_v_pod_u_y(:,:,k) = (pod_v(:,k)*ones(1,num_modes)).*pod_udy;
    pod_u_pod_v_x(:,:,k) = (pod_u(:,k)*ones(1,num_modes)).*pod_vdx;
    pod_v_pod_v_y(:,:,k) = (pod_v(:,k)*ones(1,num_modes)).*pod_vdy;
end

for k = 1:num_modes
    cdu(:,:,k) = inner_prod(pod_u_pod_u_x(:,:,k) + pod_v_pod_u_y(:,:,k), pod_u, vol_frac);
    cdv(:,:,k) = inner_prod(pod_v_pod_v_y(:,:,k) + pod_u_pod_v_x(:,:,k), pod_v, vol_frac);
end

cdu = reshape(cdu, num_modes, num_modes^2);
cdv = reshape(cdv, num_modes, num_modes^2);
cdt = cdt + cdu + cdv;

l_dot = clt;
l = cbt;
q_2dot = -cqt;
q_dot = -cct;
q = -cdt;


end