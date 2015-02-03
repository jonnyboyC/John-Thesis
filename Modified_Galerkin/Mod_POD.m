function res = Mod_POD(varargin)
% TODO really fill out varargin/help

format long g
close all
clc;

%List of fields that will be checked
fields = {  'RD_nm',        'plot_type',    'save_mod', ...
            'init',         'line_range',   'direct' ,...
            'run_num',      'type'};

% Parse problem structure provided to set it up correctly
if nargin == 1
    problem = parse_inputs(fields, @setdefaults_mod, varargin{1});
else
    fprintf('Provide a single structure as input, use help Galerkin_Proj for information.\n');
    fprintf('Using Defaults\n\n');
    problem = parse_inputs(fields, @setdefaults_mod);
end

RD_nm       = problem.RD_nm;
plot_type   = problem.plot_type;
save_mod    = problem.save_mod;
init        = problem.init;
line_range  = problem.line_range;
direct      = problem.direct;
run_num     = problem.run_num;
type        = problem.type;

% Check that parallel pool is ready
if isempty(gcp)
    parpool;
end

gcp();

% Handle File IO
if strcmp(direct, '');
    [data1, direct] = prompt_folder('POD', run_num);
    [data2, direct] = prompt_folder('Galerkin', run_num, direct);
else
    [data1, direct] = prompt_folder('POD', run_num, direct);
    [data2, direct] = prompt_folder('Galerkin', run_num, direct);
end

% Make sure folders are up to date and load collected data
update_folders(direct);

% TODO load variables based on type provided
vars1 = load(data1, 'eig_func', 'lambda2', 'pod_u', 'pod_v', 'dimensions', 'x', 'y');
vars2 = load(data2, 'num_pods', 'q', 'ci_c', 'li_c', 't3', 'modal_amp_vis2'); %'c', 'l', 't',

% TODO working up to here

% Need to explictly declare all the loaded variables for parfor loop
eig_func    = vars1.eig_func;
lambda2     = vars1.lambda2;
pod_u       = vars1.pod_u;
pod_v       = vars1.pod_v;
dimensions  = vars1.dimensions;
x           = vars1.x;
y           = vars1.y;

num_pods    = vars2.num_pods;
q           = vars2.q;
c           = vars2.ci_c;
l           = vars2.li_c;
t           = vars2.t3;
modal_amp   = vars2.modal_amp_vis2;

clear var1 var2

OG_nm = num_pods; % Original number of modes

% If less modes exist than are requested, reduced to fewer modes
if RD_nm > OG_nm
    RD_nm = OG_nm;
end

% Truncate POD
pod_ut = pod_u(:,1:OG_nm);
pod_vt = pod_v(:,1:OG_nm);

% intial guess for epsilon
epsilon_0 = sum(l*lambda2(1:OG_nm)); 

% Find initial value for line search and plot
transfer_0 = optimal_rotation(epsilon_0, c, l, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 3000);
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
        amp = floor((j)/2)/8;
        epsilon_temp(j) = epsilon_0*(1-amp*flip);
        transfer_temp(j) = optimal_rotation(epsilon_temp(j), c, l, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 6000);
    end
    
    epsilon_i       = [epsilon_i; epsilon_temp];
    transfer_i      = [transfer_i; transfer_temp];

    % Get sorted list of epsilon and unresolved transfer terms to try to find
    % and interval of sign change
    [epsilon_i, idx] = sort(epsilon_i);
    transfer_i = transfer_i(idx);

    % TODO creating separate plotting function
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
            t1 = optimal_rotation(epsilon_i(k), c, l, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000);
            t2 = optimal_rotation(epsilon_i(k+1), c, l, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000);
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
        (epsilon, c, l, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000), epsilon_range, options);
    close all;
    disp(OUTPUT);
else
    disp('no sign flip detected');
    [~, idx] = min(abs(transfer_i));
    options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', ...
    'FunValCheck', 'on');

    [epsilon_final, ~, ~, OUTPUT] = fminbnd(@(epsilon) abs(optimal_rotation...
        (epsilon, c, l, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000)), epsilon_i(idx-1), epsilon_i(idx+1), options);
    disp(OUTPUT);
end

% Final calculation of transformation matrix and new constant linear and
% quadratic terms
[~, X, C_til, L_til, Q_til] = ...
    optimal_rotation(epsilon_final, c, l, q, OG_nm, RD_nm, lambda2, eig_func, t, init, 18000);

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

%% Temp
Gal_coeff_til = [C_til L_til Q_til];

reduced_model_coeff_vis = -ode_coefficients(RD_nm, RD_nm, Gal_coeff_til);
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

tic1 = tic;
[t, modal_amp_til] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis), 0:0.0001:20, ...
    eig_func(init,1:RD_nm), options);
toc(tic1);
%% Temp


% Plot either modal amplitude, time evolution movie, or nothing
if any(strcmp(plot_type, 'amp'))
    plot_amp(modal_amp_til, t, direct, init, 'Mod');
end
if any(strcmp(plot_type, 'video'))
    plot_prediction(pod_ut, pod_vt, x, y, modal_amp_til, t, num_pods, dimensions, direct)
end
if any(strcmp(plot_type, 'fft'))
    modal_fft(modal_amp_til, 1:4, size(modal_amp,1)/10,...
        4096, [0 2000], direct, 'Mod')
end
% Save important coefficients
if save_mod == true
    save([direct '\Mod Galerkin Coeff\Coeff_m' num2str(RD_nm) 'i' num2str(init) '.mat'], ...
        'X', 'C_til', 'L_til', 'Q_til', 'pod_u_til', 'pod_v_til', 'modal_amp_til' ,'epsilon_final');
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
[~,Y_til] = ode113(@(t, y) system_odes(t, y, -reduced_model_coeff)...
    , tspan ,eig_func(init,1:RD_nm), options);   %(base line)(f = 500 Hz)

rep = error_til(Lam_til,Y_til);
toc;
end