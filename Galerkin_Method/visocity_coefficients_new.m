function [l_dot, l, q_2dot, q_dot, q] = visocity_coefficients_new(coef_problem)

mean_u      = coef_problem.mean_u;
mean_v      = coef_problem.mean_v;
x           = coef_problem.x;
y           = coef_problem.y;
pod_u       = coef_problem.pod_u;
pod_v       = coef_problem.pod_v;
dimensions  = coef_problem.dimensions;
vol_frac    = coef_problem.vol_frac;
bnd_idx     = coef_problem.bnd_idx;
z           = coef_problem.z;
run_num     = coef_problem.run_num;
override_coef = coef_problem.override_coef;
direct      = coef_problem.direct;
uniform     = coef_problem.uniform;

clear coef_problem

% Offset, currently not in use
clt = 0;
cqt = 0;
cbt = 0;
cct = 0;
cdt = 0;

num_elem = numel(x);
num_modes = size(pod_u, 2);

% max modes in memeory at once
bytesPerDouble = 8;
[~, system] = memory;
memory_limit = floor(system.PhysicalMemory.Available/(bytesPerDouble*num_modes*num_modes*1.85));

% If we are using the same number of cutoff modes and overwrite is set to
% false look for previous data
if override_coef == false;
   saved_files = dir([direct '\Viscous Coeff\Coeff_*']);
   if size(saved_files,1) ~= 0
       run_nums = regexp({saved_files.name}, num2str(run_num));
       if any(cell2mat(run_nums))
           data = load([direct '\Viscous Coeff\Coeff_' num2str(run_num) '.mat']);
           l_dot =  data.l_dot;
           l        = data.l;
           q_2dot   = data.q_2dot;
           q_dot    = data.q_dot;
           q        = data.q;
           return;
       end
   end
end

%if uniform == true
    % Use build in laplcian, and gradient functions
    [udx1, udy1, vdx1, vdy1, pod_udx1, pod_udy1, pod_vdx1, pod_vdy1, l_dot1, l1] = ...
        components_fast(x, y, mean_u, mean_v, pod_u, pod_v, dimensions, vol_frac, num_modes, num_elem);
%else
    % Old method allows for non_uniform mesh, SLOW
    [udx, udy, vdx, vdy, pod_udx, pod_udy, pod_vdx, pod_vdy, l_dot, l] = ...
        components(x, y, mean_u, mean_v, pod_u, pod_v, dimensions,  vol_frac, num_modes, num_elem, bnd_idx, z);
%end
clear x y dimensions bnd_idx z

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
pool = gcp();
delete(pool);

% If Problem has over 400 modes need to break problem into chunks
if num_modes < memory_limit;
    
    % Quadractic terms preallocation
    cduv = zeros(num_modes, num_modes, num_modes);
    
    % Calculate terms
    for k = 1:num_modes
        pod_u_pod_u_x = (pod_u(:,k)*ones(1,num_modes)).*pod_udx;
        pod_v_pod_u_y = (pod_v(:,k)*ones(1,num_modes)).*pod_udy;        
        pod_u_pod_v_x = (pod_u(:,k)*ones(1,num_modes)).*pod_vdx;
        pod_v_pod_v_y = (pod_v(:,k)*ones(1,num_modes)).*pod_vdy;
        cduv(:,:,k) = inner_prod(pod_u_pod_u_x + pod_v_pod_u_y, pod_u, vol_frac) + ...
                      inner_prod(pod_u_pod_v_x + pod_v_pod_v_y, pod_v, vol_frac);
        fprintf('%d of %d coefficients computed\n', k, num_modes);
        

    end
else
    % Create one worker to save files to harddrive
    pool = parpool('local', 1);
    
    exists = dir([direct '\Viscous Coeff\cduv.mat']);
    if size(exists,1) == 1
        delete([direct '\Viscous Coeff\cduv.mat']);
    end
    
    % Create harddrive copy to help alleviate memory problem
    data = matfile([direct '\Viscous Coeff\cduv.mat'], 'Writable', true);
    data.cduv(num_modes,num_modes,num_modes) = 0;
    
    for k = 1:num_modes
        pod_u_pod_u_x = (pod_u(:,k)*ones(1,num_modes)).*pod_udx;
        pod_v_pod_u_y = (pod_v(:,k)*ones(1,num_modes)).*pod_udy;
        pod_u_pod_v_x = (pod_u(:,k)*ones(1,num_modes)).*pod_vdx;
        pod_v_pod_v_y = (pod_v(:,k)*ones(1,num_modes)).*pod_vdy;
        cduv = inner_prod(pod_u_pod_u_x + pod_v_pod_u_y, pod_u, vol_frac);
                           inner_prod(pod_v_pod_v_y + pod_u_pod_v_x, pod_v, vol_frac);
        f = parfeval(pool, @save_cduv, 0, data, cduv, k);
        fprintf('%d of %d coefficients computed\n',  k, num_modes);
        if mod(k,20) == 0
            wait(f);
        end
    end

    delete(pool)
    cduv = data.cduv;
end

clear pod_u pod_v pod_udx pod_udy pod_vdx pod_vdy 
clear pod_u_pod_u_x pod_u_pod_v_x  pod_v_pod_u_y pod_v_pod_v_y

cduv = reshape(cduv, num_modes, num_modes^2);
q = -(cdt + cduv);

cutoff = num_modes;
save([direct '\Viscous Coeff\Coeff_' num2str(run_num) '.mat'], ...
    'l_dot', 'l', 'q_2dot', 'q_dot', 'q', 'cutoff', 'run_num', '-v7.3'); 

parpool('local', 3);
end

% allow writing to disk asynchronously
function save_cduv(data, cduv, k)
    data.cduv(:,:,k) = cduv;
end

