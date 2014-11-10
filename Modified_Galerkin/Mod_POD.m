function Mod_POD(varargin)
% TODO really fill out varargin/help
switch nargin
    case 0
        RD_nm = 8;
        plot_type = 'none';
        init = 1;
    case 1
        RD_nm = varargin{1};
        plot_type = 'none';
        init = 1;
    case 2
        RD_nm = varargin{1};
        plot_type = varargin{2};
        init = 1;
    case 3
        RD_nm = varargin{1};
        plot_type = varargin{2};
        init = varargin{3};
    otherwise
        error('Too many input arguments');
        
end

% See if the use of parrallel minimization actually is benificial
% group = gcp();

direct = 'D:\shear layer\dummie';
[data, ~] = prompt_folder({'POD', 'Galerkin'});
load(data{1}, 'eig_func', 'lambda2', 'pod_u1', 'pod_v1', 'dimensions', 'x', 'y');
load(data{2});

OG_nm = num_pods; % Original number of modes

% If less modes exist than are requested, reduced to fewer modes
if RD_nm > OG_nm
    RD_nm = OG_nm;
end

pod_ut = pod_u1(:,1:OG_nm);
pod_vt = pod_v1(:,1:OG_nm);

epsilon_i = sum(l*lambda2(1:OG_nm)); % intial guess for epsilon

options = optimset('PlotFcns', {@optimplotx, @optimplotfval});
[epsilon_final, ~, EXITFLAG, ~] = fzero(@(epsilon) optimal_rotation...
    (epsilon, ci, l, q, OG_nm, RD_nm, lambda2, eig_func, t), epsilon_i, options, init);
close all;

disp(EXITFLAG);

[~, X] = ...
    optimal_rotation(epsilon_final, ci, l, q, OG_nm, RD_nm, lambda2, eig_func, t, init);

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

if strcmp(plot_type, 'amp')
    plot_amp(modal_amp_til, t);
elseif strcmp(plot_type, 'video')
    plot_prediction(pod_u_til, pod_v_til, x, y, modal_amp_til, RD_nm, dimensions, direct);
elseif strcmp(plot_type, 'both');
    plot_prediction(pod_u_til, pod_v_til, x, y, modal_amp_til, RD_nm, dimensions, direct);
    plot_amp(modal_amp, t);
elseif strcmp(plot_pred, 'none')
else
    error('When specifying plot type, choose either amp, video, both or none');
end
end

function [rep, X] = ... %, C_til, L_til, Q_til] = ...
    optimal_rotation(epsilon, ci, l, q, OG_nm, RD_nm, lambda2, eig_func, t, init)
    
X = constrained_POD(eig_func',l,OG_nm,RD_nm,epsilon);
%inputs of constrained_POD are the POD temporal coefficients,eig_func',
%the linear Galerkin matrix, l, the transformation dimensions N and
%n, and the transfer term parameter epsi
Lam_til = X'*diag(lambda2(1:OG_nm))*X;

% Modified reduced order model coefficients
L_til = X'*l*X;
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
end