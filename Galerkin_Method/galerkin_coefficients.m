function [l, q] = galerkin_coefficients(coef_problem)
% GALERKIN_COEFFICIENT calculate the linear terms resultant from the
% visocity term and the quadractic terms resultant from the convection term
% in the Navier Stoke equations 
%
% [l, q]= GALERKIN_COEFFICENT(coef_problem)
 
% Unpack Variables
X           = coef_problem.X;
pod_U       = coef_problem.pod_U;
dimensions  = coef_problem.dimensions;
vol_frac    = coef_problem.vol_frac;
run_num     = coef_problem.run_num;
direct      = coef_problem.direct;
override_coef = coef_problem.override_coef;
use_chunks  = coef_problem.use_chunks;
bnd_idx     = coef_problem.bnd_idx;
bnd_X       = coef_problem.bnd_X;
custom      = coef_problem.custom;
uniform     = coef_problem.uniform;


clear coef_problem

[x, u] = flow_comps_ip(X, pod_U);
dims = flow_dims(X);

num_elem = numel(X.(x{1}));
num_modes = size(pod_U.(u{1}), 2);

% If we are using the same number of cutoff modes and overwrite is set to
% false look for previous data
if override_coef == false && custom == false
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

% Calculate terms, allows for nonuniform mesh
[pod_UdX, pod_U, vol_frac, l] = ...
    components(X, pod_U, dimensions, vol_frac, num_modes, num_elem, uniform, bnd_idx, bnd_X);

% Free memory
clear X bnd_idx bnd_X
h = waitbar(0, 'Calculating quadratic terms');

% If memory is an issue pass flag use chunks to write partitial results to
% disk
if use_chunks == false
    
    % If using more than 700 modes clear matlab workers for memory
    if num_modes > 700
       delete(gcp('nocreate')); 
    end
    
    % Quadractic terms preallocation
    q = zeros(num_modes, num_modes, num_modes);

    % Calculate quadractic term ((uk * grad(uj)),ui)
    for k = 1:num_modes
        for i = 1:dims;
            
            % Calculate directional derivative component
            pod_UDX = 0;
            for j = 1:dims;
                pod_UDX = pod_UDX + ...
                    repmat(pod_U.(u{j})(:,k), 1, num_modes).*pod_UdX.(u{i}).(x{j});
            end
            
            % Add components inner product to quadractic term
            q(:,:,k) = q(:,:,k) - inner_prod(pod_UDX, pod_U.(u{i}), vol_frac);
        end
        
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
    
    % Calculate quadractic term ((uk * grad(uj)),ui)
    for k = 1:num_modes
        for i = 1:dims;
            
            % Calculate directional derivative component
            pod_UDX = 0;
            for j = 1:dims;
                pod_UDX = pod_UDX + ...
                    repmat(pod_U.(u{j})(:,k), 1, num_modes).*pod_UdX.(u{i}).(x{j});
            end
            
            % Add components inner product to quadractic term
            q = q - inner_prod(pod_UDX, pod_U.(u{i}), vol_frac);
        end
        
        % Write to disk on a separate thread
        f = parfeval(pool, @save_q, 0, data, q, k);
        
        % Update wait bar
        waitbar(k/num_modes, h, sprintf('%d of %d coefficients computed\n', k, num_modes))
        
        % Allow time for writing to caught up so less data is in memory
        if mod(k,20) == 0
            wait(f);
        end
    end

    % Delete 1 matlab worker
    delete(pool)
    
    % Pull for matrix into memory
    q = data.q;
end

% Close waitbar
close(h);

% Free memory before starting workers again
clear pod_U pod_UdX

q = reshape(q, [], num_modes*num_modes);

cutoff = num_modes; %#ok<NASGU>
if ~custom
    save([direct '\Viscous Coeff\Coeff_' num2str(run_num) '_og_m' num2str(num_modes) '.mat'], ...
         'l', 'q', 'cutoff', 'run_num', '-v7.3'); 
end

if isempty(gcp('nocreate'));
    parpool('local');
end
end

% Allow writing to disk asynchronously
function save_q(data, q, k)
    data.q(:,:,k) = q;
end

