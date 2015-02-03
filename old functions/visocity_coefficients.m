function [l_dot, l, q_2dot, q_dot, q] = visocity_coefficients(mean_u, ...
    mean_v, x, y, pod_u, pod_v, dimensions, vol_frac, bnd_idx, z, over_coef, direct)
% Determine the Galkerin coefficients using the method proposed in
% Caraballo's dissertation
at = 0;
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
memory_limit = system.PhysicalMemory.Available/(bytesPerDouble*num_mode*num_modes*8);

% If we are using the same number of cutoff modes and overwrite is set to
% false look for previous data
if nargin == 11 && over_coef == false;
   saved_files = dir([direct '\Viscous Coeff\*.mat']);
   if size(saved_files,1) == 1
       data = load([direct '\Viscous Coeff\*.mat'], 'cutoff');
       if data.cutoff == num_modes && data.run_num == run_num
           data = load([direct 'Viscous Coeff\*.mat']);
           l_dot =  data.l_dot;
           l        = data.l;
           q_2dot   = data.q_2dot;
           q_dot    = data.q_dot;
           q        = data.q;
           return;
       end
   end
end

% TODO figure out what is really being calculated here
[xxi, yxi, xet, yet, aj] = metric2(x, y);

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

clear ud2x ud2y vd2x vd2y pod_ud2x pod_ud2y pod_vd2x pod_vd2y

% Find inner produce with itself for pod modes 
% TODO find if at above is the actual inner product, currently allowing for
% at to have some sort of linear offset. Fudge factor?
au = inner_prod(pod_u, pod_u, vol_frac);
av = inner_prod(pod_v, pod_v, vol_frac);
at = at + au + av; 

% Find inner product between 2nd for mean velocity and pod modes
clu = inner_prod(d2u, pod_u, vol_frac);
clv = inner_prod(d2v, pod_v, vol_frac);
l_dot = clt + clu + clv;

% Free memory
clear clt clu clv

% Find inner product between 2nd derivative of pod modes and pod modes
cbu = inner_prod(d2pod_u, pod_u, vol_frac);
cbv = inner_prod(d2pod_v, pod_v, vol_frac);
l = cbt + cbu + cbv;

% Free memory 
clear cbt cbu cbv

% TODO Look at Navier Stokes Equation
qu = mean_u.*udx+mean_v.*udy;
qv = mean_u.*vdx+mean_v.*vdy;

% Inner product of pod modes with calculated quantity. Similar situation
% with cqt we are allowing for a linear offset, may be fudge factor
cqu = inner_prod(qu, pod_u, vol_frac);
cqv = inner_prod(qv, pod_v, vol_frac);
q_2dot = -(cqt + cqu + cqv);

% Free memory 
clear cqt cqu cqv

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

ccu = inner_prod(cu, pod_u, vol_frac);
ccv = inner_prod(cv, pod_v, vol_frac);
q_dot = -(cct + ccu + ccv);

% Free memory
clear cct ccu ccv

% Quadractic Terms
cduv = zeros(num_modes, num_modes, num_modes);
cdv = zeros(num_modes, num_modes, num_modes);

% If Problem has over 400 modes need to break problem into chunks
if num_modes < memory_limit;
    
    % Quadractic terms preallocation
    cduv = zeros(num_modes, num_modes, num_modes);
    cdv = zeros(num_modes, num_modes, num_modes);
    
    % Calculate terms
    for k = 1:num_modes
        pod_u_pod_u_x = (pod_u(:,k)*ones(1,num_modes)).*pod_udx;
        pod_v_pod_u_y = (pod_v(:,k)*ones(1,num_modes)).*pod_udy;
        cduv(:,:,k) = inner_prod(pod_u_pod_u_x + pod_v_pod_u_y, pod_u, vol_frac);
        
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
            cduv = zeros(num_modes, num_modes, memory_limit);
        else
            cduv = zeros(num_modes, num_modes, total);
        end
        
        % Remainder of problem
        total = total - memory_limit;
        
        % Calculate terms
        for k = 1:size(cduv,3)
            pod_u_pod_u_x = (pod_u(:,k)*ones(1,num_modes)).*pod_udx;
            pod_v_pod_u_y = (pod_v(:,k)*ones(1,num_modes)).*pod_udy;
            cduv(:,:,k) = inner_prod(pod_u_pod_u_x + pod_v_pod_u_y, pod_u, vol_frac);
            
            pod_u_pod_v_x = (pod_u(:,k)*ones(1,num_modes)).*pod_vdx;
            pod_v_pod_v_y = (pod_v(:,k)*ones(1,num_modes)).*pod_vdy;
            cduv(:,:,k) = cduv(:,:,k) + inner_prod(pod_v_pod_v_y + pod_u_pod_v_x, pod_v, vol_frac);
        end
        
        % Save chunck and clear memory
        save([direct '\Viscous Coeff\Temp' num2str(i) '.mat'], 'cduv', '-v7.3');
        clear cdu
        fprintf('%d of %d coefficients computed', num_modes-total, num_modes);
    end
    
    % Load Individual chunk and concatenate
    for i = 1:ceil(num_modes/400)
        data = load([direct '\Viscous Coeff\Temp' num2str(i) '.mat'], 'cduv');
        if i == 1
            cduv = data.cdu;
        else
            cduv = cat(cduv, data.cdu, 3);
        end
        
        % Delete tempoary chunk
        delete([direct '\Viscous Coeff\Temp' num2str(i) '.mat']);
    end
    parpool;
end


cduv = reshape(cduv, num_modes, num_modes^2);
q = -(cdt + cduv);

% Save data
if nargin == 11
    cutoff = num_modes;
    save([direct '\Viscous Coeff\Coeff.mat'], 'l_dot', 'l', 'q_2dot', 'q_dot', 'q', 'cutoff', 'run_num'); 
end
end