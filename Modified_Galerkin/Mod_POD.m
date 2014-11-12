function Mod_POD(varargin)
close all;
clc;
% TODO really fill out varargin/help
switch nargin
    case 0
        RD_nm = 8;
        plot_type = 'none';
        save_mod = true;
        init = 1;
        line_range = 200;
    case 1
        RD_nm = varargin{1};
        plot_type = 'none';
        save_mod = true;
        init = 1;
        line_range = 200;
    case 2
        RD_nm = varargin{1};
        plot_type = varargin{2};
        save_mod = true;
        init = 1;
        line_range = 200;
    case 3
        RD_nm = varargin{1};
        plot_type = varargin{2};
        save_mod = varargin{3};
        init = 1;
        line_range = 200;
    case 4
        RD_nm = varargin{1};
        plot_type = varargin{2};
        save_mod = varargin{3};
        init = varargin{4};
        line_range = 200;
    case 5
        RD_nm = varargin{1};
        plot_type = varargin{2};
        save_mod = varargin{3};
        init = varargin{4};
        line_range = varargin{5};
    otherwise
        error('Too many input arguments');
end

% Start a parallel pool for use in brute force line search
pool = gcp();

[data, direct] = prompt_folder({'POD', 'Galerkin'});
update_folders(direct);
data1 = load(data{1}, 'eig_func', 'lambda2', 'pod_u1', 'pod_v1', 'dimensions', 'x', 'y');
data2 = load(data{2}, 'num_pods', 'q', 'ci', 'li', 't', 'modal_amp');

% Need to explictly declare all the loaded variables for parfor loop
eig_func    = data1.eig_func;
lambda2     = data1.lambda2;
pod_u1      = data1.pod_u1;
pod_v1      = data1.pod_v1;
dimensions  = data1.dimensions;
x           = data1.x;
y           = data1.y;

num_pods    = data2.num_pods;
q           = data2.q;
ci          = data2.ci;
li          = data2.li;
t           = data2.t;
modal_amp   = data2.modal_amp;

clear data1 data2

OG_nm = num_pods; % Original number of modes

% If less modes exist than are requested, reduced to fewer modes
if RD_nm > OG_nm
    RD_nm = OG_nm;
end

% Truncate POD
pod_ut = pod_u1(:,1:OG_nm);
pod_vt = pod_v1(:,1:OG_nm);

epsilon_0 = sum(li*lambda2(1:OG_nm)); % intial guess for epsilon

% Find initial value for line search and plot
transfer_0 = optimal_rotation(epsilon_0, ci, li, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 3000);
epsilon_i(1) = epsilon_0;
transfer_i(1) = transfer_0;


% Brute force Line Search, rough pass parrallel loop. Idea is to quickly
% sweep a large area to look for sign changes to use for a finer pass in
% fzero

search_space = 1:20:line_range;
exclusion_list = [];
exit = false;
for i = 1:length(search_space)-1
    range = search_space(i+1)-search_space(i);
    search_range = search_space(i):search_space(i+1)-1;
    epsilon_temp = zeros(range, 1);
    transfer_temp = zeros(range, 1);
    parfor j = search_range
        if mod(j, 2)
            flip = 1;
        else 
            flip = -1;
        end
        amp = floor((j)/2)/4;
        epsilon_temp(j) = epsilon_0*(1-amp*flip);
        transfer_temp(j) = optimal_rotation(epsilon_temp(j), ci, li, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 3000);
    end

    epsilon_i       = [epsilon_i; epsilon_temp];
    transfer_i      = [transfer_i; transfer_temp];

    % Get sorted list of epsilon and unresolved transfer terms to try to find
    % and interval of sign change
    [epsilon_i, idx] = sort(epsilon_i);
    transfer_i = transfer_i(idx);

    % Plot results of brute force sweep to see if range needs to be increased
    ax = newplot;
    plot(ax, epsilon_i, transfer_i);
    ax.Title.String = 'Line Search for Sign Change';
    ax.XLabel.String = 'Epsilon Value';
    ax.YLabel.String = 'Unresolved Transfer Term';
    drawnow;

    sign_t = sign(transfer_i);
    sign_flip = 0;

    % Find first instances of a sign flip
    for k = 1:length(sign_t)-1
        % See if coarse pass found a flip
        if sign_t(k)*sign_t(k+1) < 0 && ~any(exclusion_list == k)
            % If so run values using more minimization steps then check
            % again
            t1 = optimal_rotation(epsilon_i(k), ci, li, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000);
            t2 = optimal_rotation(epsilon_i(k+1), ci, li, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000);
            if t1*t2 < 0 
                sign_flip = k;
                exit = true;
                break;
            else
                exclusion_list = [exclusion_list, k];
            end
        end
    end
    if exit == true;
        break;
    end
end

% define bounded range for fzero or exit with error, 
% TODO in the future may want to either take the closest value, or simply
% extend the region that is searched 
if sign_flip ~= 0
    epsilon_range = [epsilon_i(k) epsilon_i(k+1)];
    options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', ...
    'FunValCheck', 'on');

    [epsilon_final, ~, ~, OUTPUT] = fzero(@(epsilon) optimal_rotation...
        (epsilon, ci, li, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000), epsilon_range, options);
    close all;
    disp(OUTPUT);
else
    disp('no sign flip detected');
    [~, idx] = min(abs(transfer_i));
    options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', ...
    'FunValCheck', 'on');

    [epsilon_final, ~, ~, OUTPUT] = fminbnd(@(epsilon) abs(optimal_rotation...
        (epsilon, ci, li, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000)), epsilon_i(idx-1), epsilon_i(idx+1), options);
    disp(OUTPUT);
end

% Final calculation of transformation matrix and new constant linear and
% quadratic terms
[~, X, C_til, L_til, Q_til] = ...
    optimal_rotation(epsilon_final, ci, li, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000);

% Transform pod modes and modal amplitudes 
pod_u_til = zeros(size(pod_vt,1), RD_nm);
pod_v_til = zeros(size(pod_vt,1), RD_nm);
modal_amp_til = zeros(size(modal_amp,1), RD_nm);
for i = 1:RD_nm
    for j = 1:size(pod_ut,1)
        pod_u_til(j,i) = sum(X(:,i)'.*pod_ut(j,:));
        pod_v_til(j,i) = sum(X(:,i)'.*pod_vt(j,:));
    end
    for j = 1:size(modal_amp,1)
        modal_amp_til(j,i) = sum(X(:,i)'.*modal_amp(j,:));
    end
end

% Plot either modal amplitude, time evolution movie, or nothing
if strcmp(plot_type, 'amp')
    plot_amp(modal_amp_til, t, direct, init, true);
elseif strcmp(plot_type, 'video')
    plot_prediction(pod_u_til, pod_v_til, x, y, modal_amp_til, RD_nm, dimensions, direct);
elseif strcmp(plot_type, 'both');
    plot_prediction(pod_u_til, pod_v_til, x, y, modal_amp_til, RD_nm, dimensions, direct);
    plot_amp(modal_amp, t, direct, init, true);
elseif strcmp(plot_pred, 'none')
else
    epi_error('When specifying plot type, choose either amp, video, both or none');
end

% Save important coefficients
if save_mod == true
    save([direct '\Mod Galerkin Coeff\Coeff_m' num2str(RD_nm) 'i' num2str(init) '.mat'], ...
        'X', 'C_til', 'L_til', 'Q_til', 'pod_u_til', 'pod_v_til', 'modal_amp_til');
end
end


% Use this function to minimally rotate
function [rep, X, C_til, L_til, Q_til] = ...
    optimal_rotation(epsilon, ci, li, q, OG_nm, RD_nm, lambda2, eig_func, t, init, evals)
    
tic;
X = constrained_POD(eig_func',li,OG_nm,RD_nm,epsilon,evals);
%inputs of constrained_POD are the POD temporal coefficients,eig_func',
%the linear Galerkin matrix, li, the transformation dimensions N and
%n, and the transfer term parameter epsi
Lam_til = X'*diag(lambda2(1:OG_nm))*X;

% Modified reduced order model coefficients
L_til = X'*li*X;
C_til = X'*ci;
Q_til=zeros(RD_nm,RD_nm,RD_nm);
Q = reshape(q,OG_nm,OG_nm,OG_nm);
for i = 1:RD_nm
    for j = 1:RD_nm
        for k = 1:RD_nm
            for p = 1:OG_nm
                for r = 1:OG_nm
                    for q = 1:OG_nm
                        Q_til(i,j,k) = Q_til(i,j,k) + X(p,i)*Q(p,q,r)*X(q,j)*X(r,k);
                    end
                end
            end
        end
    end
end
% solve new system of equation
Q_til=reshape(Q_til,RD_nm,RD_nm*RD_nm);
fc_til=[C_til L_til Q_til]; % constant, linear and quadratic terms coefficients

reduced_model_coeff= ode_coefficients(RD_nm,RD_nm,fc_til);

%Solution of the system of equation
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
%(Control implicit in the modes)

% Look more into tspan;
tspan = t;
[~,Y_til] = ode15s(@(t, y) system_odes(t, y, -reduced_model_coeff)...
    , tspan ,eig_func(init,1:RD_nm), options);   %(base line)(f = 500 Hz)

rep = error_til(Lam_til,Y_til);
toc;
end