function [l_dot, l, q_2dot, q_dot, q] = visocity_coefficients_ws(coef_problem)

% TODO really comb this over 

mean_u      = coef_problem.mean_u;
mean_v      = coef_problem.mean_v;
x           = coef_problem.x;
y           = coef_problem.y;
pod_u       = coef_problem.pod_u;
pod_v       = coef_problem.pod_v;
dimensions  = coef_problem.dimensions;
vol_frac    = coef_problem.vol_frac;
run_num     = coef_problem.run_num;
override_coef = coef_problem.override_coef;
direct      = coef_problem.direct;
uniform     = coef_problem.uniform;
bnd_x       = coef_problem.bnd_x;
bnd_y       = coef_problem.bnd_y;

clear coef_problem

% Offset, currently not in use
clt = 0;
cqt = 0;
cbt = 0;
cct = 0;
cdt = 0;

num_elem = numel(x);
num_modes = size(pod_u, 2);

% If we are using the same number of cutoff modes and overwrite is set to
% false look for previous data
if override_coef == false;
   saved_files = dir([direct '\Viscous Coeff\Coeff_*']);
   if size(saved_files,1) ~= 0
       match_run = regexp({saved_files.name}, num2str(run_num));
       match_modes = regexp({saved_files.name}, ['m' num2str(num_modes) '\.']);
       weak = regexp({saved_files.name}, 'wk');
       if any(~cellfun(@isempty, match_run) & ~cellfun(@isempty, match_modes) & ~cellfun(@isempty, weak))
           data = load([direct '\Viscous Coeff\Coeff_' num2str(run_num) '_wk_m' num2str(num_modes) '.mat']);
           l_dot =  data.l_dot;
           l        = data.l;
           q_2dot   = data.q_2dot;
           q_dot    = data.q_dot;
           q        = data.q;
           return;
       end
   end
end

if uniform == true
    % Use build in laplcian, and gradient functions
    [udx, udy, vdx, vdy, pod_udx, pod_udy, pod_vdx, pod_vdy] = ...
        components_ws_fast(x, y, mean_u, mean_v, pod_u, pod_v, dimensions, num_modes, num_elem);
else
    [udx, udy, vdx, vdy, pod_udx, pod_udy, pod_vdx, pod_vdy] = ...
        components_ws(x, y, mean_u, mean_v, pod_u, pod_v, dimensions, num_modes);
end

qu = mean_u.*udx + mean_v.*udy;
qv = mean_u.*vdx + mean_v.*vdy;

cqu = inner_prod(qu, pod_u, vol_frac);
cqv = inner_prod(qv, pod_v, vol_frac);

cqux = inner_prod(udx, pod_udx, vol_frac);
cqvx = inner_prod(vdx, pod_vdx, vol_frac);

cquy = inner_prod(udy, pod_udy, vol_frac);
cqvy = inner_prod(vdy, pod_vdy, vol_frac);

if  min(min(y)) < -0.0001
    cqu1= inner_prods(udx,   pod_u,x,y,1);
    cqu2= inner_prods(udy,   pod_u,x,y,2);
    cqu3= inner_prods(udx,   pod_u,x,y,3);
    cqu4= inner_prods(mean_u,pod_u,x,y,4);
    cqv1= inner_prods(vdx,   pod_v,x,y,1);
    cqv2= inner_prods(vdy,   pod_v,x,y,2);
    cqv3= inner_prods(vdx,   pod_v,x,y,3);
    cqv4= inner_prods(mean_v,pod_v,x,y,4);
else
    cqu1=0;
    cqu2=0;
    cqu3=0;
    cqu4=0;
    cqv1=0;
    cqv2=0;
    cqv3=0;
    cqv4=0;
end

cqnt = -(cqux + cquy + cqu1 - cqu2 - cqu3 - cqu4)...
       -(cqvx + cqvy + cqv1 - cqv2 - cqv3 - cqv4);
   
l_dot = clt + cqnt;
q_2dot = cqt + cqu + cqv;

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

ccux = inner_prod(pod_udx, pod_udx, vol_frac);
ccvx = inner_prod(pod_vdx, pod_vdx, vol_frac);

ccuy = inner_prod(pod_udy, pod_udy, vol_frac);
ccvy = inner_prod(pod_vdy, pod_vdy, vol_frac);

if  true %min(min(y)) < -0.0001
    ccu1= inerpros(pod_udx, pod_u,x,y,1);
    ccu2= inerpros(pod_udy, pod_u,x,y,2);
    ccu3= inerpros(pod_udx, pod_u,x,y,3);
 %   ccu4= inerpros(pod_u,   pod_u,x,y,4);
    ccv1= inerpros(pod_vdx, pod_v,x,y,1);
    ccv2= inerpros(pod_vdy, pod_v,x,y,2);
    ccv3= inerpros(pod_vdx, pod_v,x,y,3);
 %   ccv4= inerpros(pod_v,   pod_v,x,y,4);
else
    ccu1=0;
    ccu2=0;
    ccu3=0;
    ccu4=0;
    ccv1=0;
    ccv2=0;
    ccv3=0;
    ccv4=0;
end

surf_u = surf_inner_prod(pod_udx, pod_u, vol_frac, bnd_x) + ...
         surf_inner_prod(pod_vdv, pod_u, vol_frac, bnd_y);
surf_v = surf_inner_prod(pod_udx, pod_u, vol_frac, bnd_x) + ...
         surf_inner_prod(pod_vdv, pod_u, vol_frac, bnd_y);


ccnt = (ccux + ccuy + ccvx + ccvy);
l = cbt + ccnt;
   
q_dot = cct -(inner_prod(cu, pod_u, vol_frac) + inner_prod(cv, pod_v, vol_frac));
q_dot = cct + q_dot;


% max modes in memeory at once
bytesPerDouble = 8;
[~, system] = memory;
memory_limit = floor(system.PhysicalMemory.Available/(bytesPerDouble*num_modes*num_modes*1.2));

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

cduv = reshape(cduv, num_modes, num_modes*num_modes);
q = -(cdt + cduv);

cutoff = num_modes;
save([direct '\Viscous Coeff\Coeff_' num2str(run_num) '_wk_m' num2str(num_modes) '.mat'], ...
    'l_dot', 'l', 'q_2dot', 'q_dot', 'q', 'cutoff', 'run_num', '-v7.3'); 

end

% allow writing to disk asynchronously
function save_cduv(data, cduv, k)
    data.cduv(:,:,k) = cduv;
end
