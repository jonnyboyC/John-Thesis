function Galerkin_Proj(varargin)
% GALERKIN_PROJ will calculate time coefficients for a provided POD basis
%
% GALERKIN_PROJ(NUM_PODS) generated time coefficients for NUM_PODS number
% of pod modes
%
% Galerkin_PROJ(NUM_PODS, PLOT_PRED) generated time coefficients for
% NUM_PODS number of pod modes. If plot_pred is true a short move is saved
% to the data folder of the predicted flow.

%%%%%%%%%%%%%%% VARIABLES USED FROM CHABOT MAIN CODE %%%%%%%%%%%%%%%%%%%%%%
%
% pod_u1: pods in the u direction that have had thier sign flipped
% pod_v1: pods in the v direction that have had their sign flipped
% mean_u: mean velocity in the u direction
% mean_v: mean velocity in the v direction
% dimensions: size of the piv images
% x: matrix of x values
% y: matrix of y values
% bnd_idx: matrix of -1's , 0's, 1's defining the detected boundaries
%
%

% Set format, clear figures, and set up correct directory
format long g
close all

%%%%%%%%%%%% ALTER THIS TO MAKE PORTABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('D:\shear layer');

% set up function 
switch nargin
    case 0 
        % Default: run simulation for 10 pod modes, don't plot, run for
        % 100s, don't save coefficients
        num_pods = 10;
        plot_pred = false;
        save_coef = false;
        tspan = [0 100];
    case 1
        num_pods = varargin{1};
    	plot_pred = false;
        save_coef = false;
        tspan = [0 100];
    case 2
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = false;        
        tspan = [0 100];
    case 3 
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = [0 100];
    case 4
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = varargin{4}      
    otherwise
        error('Too many input arguments');
end

% matlabpool local 4
[data, direct] = prompt_folder('POD');
load(data{1});

%% TODO Chunk of variables need to sort them out
% List of variables in this function i'm not sure about
% dl   tn   tm   acm   plm   pqm   fr   tr
% l   q   niu   ci   ni   li   fcuhi1

mu0=1;
sc=1;
Re0=1140;                   %Reynolds number
z=ones(size(x));            %Depth of velocity field 
M0=.25;

%Solution of the time coefficient
sf=50000;                   %Sampling frequency
Ns=2048;                    %Number samples?
dl=mu0/(sc*sf);             
tn=0:dl:dl*Ns*6*40;         
tm=zeros(10,size(tn,2));
acm=zeros(10,size(tn,2),4);
plm=zeros(10,size(tn,2));
pqm=zeros(10,size(tn,2));

fr=1000*sc/mu0;
tr=tn;

% TODO if we decide we want to calculate for range of number of pod modes,
% place for loop here

% Take a truncation of calculated pods
pod_ut = pod_u1(:,1:num_pods);
pod_vt = pod_v1(:,1:num_pods);

[l_dot, l, q_2dot, q_dot, qi] = visocity_coefficients(mean_u, mean_v, ...
    x, y, pod_ut, pod_vt, dimensions, vol_frac, bnd_idx, z);

niu = viscious_dis(eig_func_norm, num_pods, lambda2, l, q_dot, qi);
ci = l_dot/Re0.*niu+q_2dot;

% matrix of niu with extra scaling on off diagonal
ni = diag(niu) + (ones(num_pods)-eye(num_pods))/Re0;
li = ni.*l+q_dot;
fcuhi1 = [ci li qi];

% Will reduce number of coefficients if we want a smaller model
reduced_model_coeff = ode_coefficients(num_pods, num_pods, fcuhi1);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

tic;
[t, modal_amp] = ode45(@(t,y) system_odes(t,y,-reduced_model_coeff), tspan, ...
    eig_func_norm(1,1:num_pods), options);
toc;

% Provide only modal flucations ie only turblent portion of modes
% TODO check the validity of this statement
modal_amp = modal_amp - ones(size(modal_amp,1), 1)*mean(modal_amp);
% plot(t, modal_amp(:,1), 'b');

if plot_pred == true
    plot_prediction(pod_ut, pod_vt, x, y, modal_amp, num_pods, dimensions, direct)
end

if save_coef == true
    save([direct '\Galerkin Coeff\Coeff.mat'], 'ci', 'li', 'qi');
end
