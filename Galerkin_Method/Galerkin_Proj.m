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

% Set format, clear figures, and set up correct directory
format long g
close all

%%%%%%%%%%%% ALTER THIS TO MAKE PORTABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('D:\shear layer');

% set up function 
switch nargin
    case 0 
        % Default: run simulation for 10 pod modes, don't plot, run for
        % 100s, save coefficients
        num_pods = 10;
        plot_pred = 'none';
        save_coef = true;
        tspan = [0 100];
        init = 1;
        direct = ''
    case 1
        num_pods = varargin{1};
    	plot_pred = 'none';
        save_coef = true;
        tspan = [0 100];
        init = 1;
        direct = '';
    case 2
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = true;        
        tspan = [0 100];
        init = 1;
        direct = '';  
    case 3 
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = [0 100];
        init = 1;
        direct = '';
    case 4
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = varargin{4};  
        init = 1;
        direct = '';
    case 5
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = varargin{4};  
        init = varargin{5};
        direct = '';
    case 6
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = varargin{4};  
        init = varargin{5};
        direct = varargin{6};
    otherwise
        error('Too many input arguments');
end

% matlabpool local 4
if strcmp(direct, '');
    [data, direct] = prompt_folder('POD');
else
    [data, direct] = prompt_folder('POD', direct);
end
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

[l_dot, l, q_2dot, q_dot, q] = visocity_coefficients(mean_u, mean_v, ...
    x, y, pod_ut, pod_vt, dimensions, vol_frac, bnd_idx, z);

niu = viscious_dis(eig_func, num_pods, lambda2, l, q_dot, q);
ci = l_dot/Re0.*niu+q_2dot;

% matrix of niu with extra scaling on off diagonal
ni = diag(niu) + (ones(num_pods)-eye(num_pods))/Re0;
li = ni.*l+q_dot;
fcuhi1 = [ci li q];

% Will reduce number of coefficients if we want a smaller model
reduced_model_coeff = -ode_coefficients(num_pods, num_pods, fcuhi1);
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

% Was previously using ode113 now using 15s, preformance up 4x, may need to
% change in the future.
tic;
[t, modal_amp] = ode15s(@(t,y) system_odes(t,y,reduced_model_coeff), tspan, ...
    eig_func(init,1:num_pods), options);
toc;

% Was previously using eig_func_norm

% Provide only modal flucations ie only turblent portion of modes
% TODO check the validity of this statement
% modal_amp = modal_amp - ones(size(modal_amp,1), 1)*mean(modal_amp);


if strcmp(plot_pred, 'amp')
    plot_amp(modal_amp(:, 1:num_pods), t, direct, init);
elseif strcmp(plot_pred, 'video')
    plot_prediction(pod_ut, pod_vt, x, y, modal_amp, t, num_pods, dimensions, direct)
elseif strcmp(plot_pred, 'both');
    plot_amp(modal_amp(:, 1:num_pods), t, direct, init);
    plot_prediction(pod_ut, pod_vt, x, y, modal_amp, t, num_pods, dimensions, direct)
elseif strcmp(plot_pred, 'none')
else
    error('When specifying plot type, choose either amp, video, both or none');
end

if save_coef == true
    save([direct '\Galerkin Coeff\Coeff_m' num2str(num_pods) 'i' num2str(init) '.mat'],...
        'ci', 'li', 'q', 'num_pods', 'modal_amp', 't', 'l_dot', 'l', 'q_2dot', 'q_dot', 'q');
end
