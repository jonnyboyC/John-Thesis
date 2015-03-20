function [l, q] = visocity_coefficients(coef_problem)
% Unpack Variables
x           = coef_problem.x;
y           = coef_problem.y;
pod_u       = coef_problem.pod_u;
pod_v       = coef_problem.pod_v;
dimensions  = coef_problem.dimensions;
vol_frac    = coef_problem.vol_frac;
z           = coef_problem.z;
run_num     = coef_problem.run_num;
direct      = coef_problem.direct;
override_coef = coef_problem.override_coef;
uniform     = coef_problem.uniform;
use_chunks  = coef_problem.use_chunks;
bnd_idx     = coef_problem.bnd_idx;


clear coef_problem

num_elem = numel(x);
num_modes = size(pod_u, 2);

% If we are using the same number of cutoff modes and overwrite is set to
% false look for previous data
if override_coef == false;
   saved_files = dir([direct '\Viscous Coeff\Coeff_*']);
   if size(saved_files,1) ~= 0
       match_run = regexp({saved_files.name}, num2str(run_num));
       match_modes = regexp({saved_files.name}, ['m' num2str(num_modes) '\.']);
       orig = regexp({saved_files.name}, 'og');
       if any(~cellfun(@isempty, match_run) & ~cellfun(@isempty, match_modes) & ~cellfun(@isempty, orig))
           data = load([direct '\Viscous Coeff\Coeff_' num2str(run_num) '_og_m' num2str(num_modes) '.mat']);
           l        = data.l;
           q        = data.q;
           return;
       end
   end
end


if ~uniform 
    % Use build in laplcian, and gradient functions
    [pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac, l] = ...
        components_fast(x, y, pod_u, pod_v, dimensions, vol_frac, num_modes, num_elem, bnd_idx);
else
    % Old method allows for non_uniform mesh, SLOW
    [pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac, l] = ...
        components(x, y, pod_u, pod_v, dimensions, vol_frac, num_modes, num_elem, bnd_idx, z);
end

% Free memory
clear x y dimensions bnd_idx z mean_u mean_v
h = waitbar(0, 'Calculating quadratic terms');

% If Problem has over 400 modes need to break problem into chunks
if use_chunks == false
    
    if num_modes > 700
       delete(gcp('nocreate')); 
    end
    
    % Quadractic terms preallocation
    q = zeros(num_modes, num_modes, num_modes);

    % Calculate terms
    for k = 1:num_modes
        pod_u_pod_u_x = (repmat(pod_u(:,k),1,num_modes)).*pod_udx;
        pod_v_pod_u_y = (repmat(pod_v(:,k),1,num_modes)).*pod_udy;        
        pod_u_pod_v_x = (repmat(pod_u(:,k),1,num_modes)).*pod_vdx;
        pod_v_pod_v_y = (repmat(pod_v(:,k),1,num_modes)).*pod_vdy;
        q(:,:,k) = -inner_prod(pod_u_pod_u_x + pod_v_pod_u_y, pod_u, vol_frac) ...
                   -inner_prod(pod_u_pod_v_x + pod_v_pod_v_y, pod_v, vol_frac);
               
        % Update wait bar
        waitbar(k/num_modes, h, sprintf('%d of %d coefficients computed\n', k, num_modes))
    end
else
    % Create one worker to save files to harddrive
    delete(gcp);
    pool = parpool('local', 1);
    
    exists = dir([direct '\Viscous Coeff\cduv.mat']);
    if size(exists,1) == 1
        delete([direct '\Viscous Coeff\cduv.mat']);
    end
    
    % Create harddrive copy to help alleviate memory problem
    data = matfile([direct '\Viscous Coeff\cduv.mat'], 'Writable', true);
    data.q(num_modes,num_modes,num_modes) = 0;
    
    for k = 1:num_modes
        pod_u_pod_u_x = (repmat(pod_u(:,k),1,num_modes)).*pod_udx;
        pod_v_pod_u_y = (repmat(pod_v(:,k),1,num_modes)).*pod_udy;        
        pod_u_pod_v_x = (repmat(pod_u(:,k),1,num_modes)).*pod_vdx;
        pod_v_pod_v_y = (repmat(pod_v(:,k),1,num_modes)).*pod_vdy;
        q = -inner_prod(pod_u_pod_u_x + pod_v_pod_u_y, pod_u, vol_frac) ...
            -inner_prod(pod_v_pod_v_y + pod_u_pod_v_x, pod_v, vol_frac);
        f = parfeval(pool, @save_q, 0, data, q, k);
        
        % Update wait bar
        waitbar(k/num_modes, h, sprintf('%d of %d coefficients computed\n', k, num_modes))
        
        % Allow writing to caught up
        if mod(k,20) == 0
            wait(f);
        end
    end

    delete(pool)
    q = data.q;
end

% close waitbar
close(h);

% Free memory
clear pod_u pod_v pod_udx pod_udy pod_vdx pod_vdy f
clear pod_u_pod_u_x pod_u_pod_v_x pod_v_pod_u_y pod_v_pod_v_y

q = reshape(q, [], num_modes*num_modes);

cutoff = num_modes;
save([direct '\Viscous Coeff\Coeff_' num2str(run_num) '_og_m' num2str(num_modes) '.mat'], ...
     'l', 'q', 'cutoff', 'run_num', '-v7.3'); 

if isempty(gcp('nocreate'));
    parpool('local', 4);
end
end

% allow writing to disk asynchronously
function save_q(data, q, k)
    data.q(:,:,k) = q;
end

