function [l_dot, l, q_2dot, q_dot, q] = visocity_coefficients_fast(mean_u, ...
    mean_v, x, y, pod_u, pod_v, dimensions, vol_frac, run_num, over_coef, direct)

% max modes in memeory at once
memory_limit = 250;

clt = 0;
cqt = 0;
cbt = 0;
cct = 0;
cdt = 0;

num_elem = numel(x);
num_modes = size(pod_u, 2);

% If we are using the same number of cutoff modes and overwrite is set to
% false look for previous data
if nargin == 11 && over_coef == false;
   saved_files = dir([direct '\Viscous Coeff\*.mat']);
   if size(saved_files,1) == 1
       data = load([direct '\Viscous Coeff\Coeff.mat'], 'cutoff', 'run_num');
       if data.cutoff == num_modes && data.run_num == run_num
           data = load([direct '\Viscous Coeff\Coeff.mat']);
           l_dot =  data.l_dot;
           l        = data.l;
           q_2dot   = data.q_2dot;
           q_dot    = data.q_dot;
           q        = data.q;
           return;
       end
   end
end

hu = mean(mean(x(1:end-1,:),2) - mean(x(2:end,:),2));
hv = mean(mean(y(:,1:end-1),1) - mean(y(:,2:end),1));

% TODO determine why these all need to be negative to match
mean_u = reshape(-mean_u, dimensions(1), dimensions(2));
mean_v = reshape(-mean_v, dimensions(1), dimensions(2));
pod_u = regroup(-pod_u, dimensions);
pod_v = regroup(-pod_v, dimensions);

[udy,udx] = gradient(mean_u, hu, hv);   % gradient of mean u flow
[vdy,vdx] = gradient(mean_v, hu, hv);   % gradient of mean v flow
d2u = del2(mean_u, hu, hv)*4;           % laplacain of mean u flow
d2v = del2(mean_v, hu, hv)*4;           % laplacian of mean v flow

pod_udx = zeros(size(pod_u));   % x derivative pod_u
pod_udy = zeros(size(pod_u));   % y derivative pod_u
d2pod_u = zeros(size(pod_u));   % laplacain pod_u

pod_vdx = zeros(size(pod_v));   % x derivative pod_v
pod_vdy = zeros(size(pod_v));   % y derivative poy_v
d2pod_v = zeros(size(pod_v));   % laplacian pod_v

for i = 1:size(pod_u,3)
   [pod_udy(:,:,i),pod_udx(:,:,i)] = gradient(pod_u(:,:,i), hu, hv);
   [pod_vdy(:,:,i),pod_vdx(:,:,i)] = gradient(pod_v(:,:,i), hu, hv);
   d2pod_u(:,:,i) = del2(pod_u(:,:,i), hu, hv)*4;
   d2pod_v(:,:,i) = del2(pod_v(:,:,i), hu, hv)*4;
end

% Convert all matrix quantities to vectors for mean derivatives and
% laplacian
udx     = reshape(udx, num_elem, 1);
udy     = reshape(udy, num_elem, 1);
vdx     = reshape(vdx, num_elem, 1);
vdy     = reshape(vdy, num_elem, 1);
d2u     = reshape(d2u, num_elem, 1);
d2v     = reshape(d2v, num_elem, 1);

% Convert all matrix quantities to vectors for pod derivatives and
% laplacian
pod_udx     = reshape(pod_udx, num_elem, num_modes);
pod_udy     = reshape(pod_udy, num_elem, num_modes);
pod_vdx     = reshape(pod_vdx, num_elem, num_modes);
pod_vdy     = reshape(pod_vdy, num_elem, num_modes);
d2pod_u     = reshape(d2pod_u, num_elem, num_modes);
d2pod_v     = reshape(d2pod_v, num_elem, num_modes);

mean_u = reshape(mean_u, num_elem, 1);
mean_v = reshape(mean_v, num_elem, 1);
pod_u = reshape(pod_u, num_elem, num_modes);
pod_v = reshape(pod_v, num_elem, num_modes);

l_dot = inner_prod(d2u, pod_u, vol_frac) + inner_prod(d2v, pod_v, vol_frac);
l = inner_prod(d2pod_u, pod_u, vol_frac) + inner_prod(d2pod_v, pod_v, vol_frac);

% Free memory
clear d2u d2v d2pod_u d2pod_v

l_dot = l_dot + clt;
l     = l + cbt;

% Free memory 
clear clt cbt

qu = mean_u.*udx + mean_v.*udy;
qv = mean_u.*vdx + mean_v.*vdy;

q_2dot = -(inner_prod(qu, pod_u, vol_frac) + inner_prod(qv, pod_v, vol_frac));
q_2dot = q_2dot + cqt;

% Free memory 
clear cqt qu qv

% TODO look for better names for all these variables 
% Linear terms
% Group 1
u_pod_ux = mean_u*ones(1,num_modes).*pod_udx;
u_pod_vx = mean_u*ones(1,num_modes).*pod_vdx;

clear mean_u

% Group 2
v_pod_uy = mean_v*ones(1,num_modes).*pod_udy;
v_pod_vy = mean_v*ones(1,num_modes).*pod_vdy;

% Free memeory
clear mean_v

% Group 3
ux_pod_u = pod_u.*(udx*ones(1,num_modes));
vx_pod_u = pod_u.*(vdx*ones(1,num_modes));

% Free memory
clear udx vdx

% Group 4
uy_pod_v = pod_v.*(udy*ones(1,num_modes));
vy_pod_v = pod_v.*(vdy*ones(1,num_modes));

% Free memory
clear udy vdy

% Sum of terms 
cu = u_pod_ux + v_pod_uy + ux_pod_u + uy_pod_v;
cv = u_pod_vx + v_pod_vy + vx_pod_u + vy_pod_v;

clear u_pod_ux v_pod_uy ux_pod_u uy_pod_v
clear u_pod_vx v_pod_vy vx_pod_u vy_pod_v

q_dot = -(inner_prod(cu, pod_u, vol_frac) + inner_prod(cv, pod_v, vol_frac));
q_dot = q_dot + cct;

% Free memory
clear cct cu cv

% If Problem has over 400 modes need to break problem into chunks
if num_modes < memory_limit;
    
    % Quadractic terms preallocation
    cdu = zeros(num_modes, num_modes, num_modes);
    cdv = zeros(num_modes, num_modes, num_modes);
    
    % Calculate terms
    for k = 1:num_modes
        pod_u_pod_u_x = (pod_u(:,k)*ones(1,num_modes)).*pod_udx;
        pod_v_pod_u_y = (pod_v(:,k)*ones(1,num_modes)).*pod_udy;
        cdu(:,:,k) = inner_prod(pod_u_pod_u_x + pod_v_pod_u_y, pod_u, vol_frac);
        
        pod_u_pod_v_x = (pod_u(:,k)*ones(1,num_modes)).*pod_vdx;
        pod_v_pod_v_y = (pod_v(:,k)*ones(1,num_modes)).*pod_vdy;
        cdv(:,:,k) = inner_prod(pod_v_pod_v_y + pod_u_pod_v_x, pod_v, vol_frac);
    end
else
    % delete parrallel pool to free up memeory
    pool = gcp();
    delete(pool);
    
    % total to be calculted
    total = num_modes;
    
    % Break problem into 400 mode chunks
    for i = 1:ceil(total/memory_limit)
        
        % Quadractic term preallocation
        if total > memory_limit
            cdu = zeros(num_modes, num_modes, memory_limit);
            cdv = zeros(num_modes, num_modes, memory_limit);
        else
            cdu = zeros(num_modes, num_modes, memory_limit);
            cdv = zeros(num_modes, num_modes, memory_limit);
        end
        
        % Remainder of problem
        total = total - memory_limit;
        
        % Calculate terms
        for k = 1:size(cdu,3)
            pod_u_pod_u_x = (pod_u(:,k)*ones(1,num_modes)).*pod_udx;
            pod_v_pod_u_y = (pod_v(:,k)*ones(1,num_modes)).*pod_udy;
            cdu(:,:,k) = inner_prod(pod_u_pod_u_x + pod_v_pod_u_y, pod_u, vol_frac);
            
            pod_u_pod_v_x = (pod_u(:,k)*ones(1,num_modes)).*pod_vdx;
            pod_v_pod_v_y = (pod_v(:,k)*ones(1,num_modes)).*pod_vdy;
            cdv(:,:,k) = inner_prod(pod_v_pod_v_y + pod_u_pod_v_x, pod_v, vol_frac);
        end
        
        % Save chunck and clear memory
        save([direct '\Viscous Coeff\Temp' num2str(i) '.mat'], 'cdu', 'cdv', '-v7.3');
        clear cdu cdv
        fprintf('%d of %d coefficients computed', num_modes-total, num_modes);
    end
    
    % Load Individual chunk and concatenate
    for i = 1:ceil(num_modes/memory_limit)
        data = load([direct '\Viscous Coeff\Temp' num2str(i) '.mat'], 'cdu', 'cdv');
        if i == 1
            cdu = data.cdu;
            cdv = data.cdv;
        else
            cdu = cat(cdu, data.cdu, 3);
            cdv = cat(cdv, data.cdv, 3);
        end
        
        % Delete tempoary chunk
        delete([direct '\Viscous Coeff\Temp' num2str(i) '.mat']);
    end
    parpool;
end

cdu = reshape(cdu, num_modes, num_modes^2);
cdv = reshape(cdv, num_modes, num_modes^2);
q = -(cdt + cdu + cdv);

% Save data
if nargin == 11
    cutoff = num_modes;
    save([direct '\Viscous Coeff\Coeff.mat'], 'l_dot', 'l', 'q_2dot', 'q_dot', 'q', 'cutoff', 'run_num', '-v7.3'); 
end
end

