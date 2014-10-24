function Mod_POD(varargin)
% TODO really fill out varargin/help
switch nargin
    case 0
        RD_nm = 8;
    case 1
        RD_nm = varargin{1};
    otherwise
        error('Too many input arguments');
end
[data, direct] = prompt_folder({'POD', 'Galerkin'});
load(data{1}, 'eig_func_norm', 'lambda2', 'num_pods');
load(data{2});

OG_nm = num_pods; % Original number of modes

% If less modes exist than are requested, reduced to fewer modes
if RD_nm > OG_nm
    RD_nm = OG_nm;
end

rep = 1;                            % residual from epsilon
epsi = sum(li*lambda2(1:OG_nm));    % intial guess for epsilon
nv = 1;                             % Leaving for now will probably delete
iter = 0;                           % optimization iteratin

esc=0;                              % escape criteria
epst = epsi*1.01;                   % Don't want to rotate worse than initial guess

% TODO have figure update x and y data ins
figure
hold on
while abs(epst) > abs(epsi) %500 %(rep ~= 0 ||)ct < 1 % 
    [X, L, epsilon, lambda] = constrained_POD(eig_func_norm',li,OG_nm,RD_nm,epsi);
    %inputs of constrained_POD are the POD temporal coefficients,eig_func_norm',
    %the linear Galerkin matrix, li, the transformation dimensions N and
    %n, and the transfer term parameter epsi
    Lam_til = X'*diag(lambda2(1:OG_nm))*X;
    
    % Modified reduced order model coefficients
    L_til = X'*li*X;
    C_til = X'*ci;
    Q_til=zeros(RD_nm,RD_nm,RD_nm);
    Q = reshape(qi,OG_nm,OG_nm,OG_nm);
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
    
    in=1; % initial point (to)
    
    reduced_model_coeff= ode_coefficients(RD_nm,RD_nm,fc_til);
    
    %Solution of the system of equation
    options = odeset('RelTol',1e-6,'AbsTol',1e-8);
    %(Control implicit in the modes)
    
    % Look more into tspan;
    tspan = [0 100]; % tn(1):dl:tn(2000);
    [~,Y_til] = ode113(@(t, y) system_odes(t, y, -reduced_model_coeff)...
        , tspan ,eig_func_norm(in,1:RD_nm), options);   %(base line)(f = 500 Hz)
    
    rep = error_til(Lam_til,Y_til);
    plot(epsi,rep,'*')
    drawnow;
       % saveas(gcf, ['M:\MASTERS\Spring 2014 Thesis figures\alpha',num2str(ia),'_Mod_POD_bl_',num2str(N),' to ',num2str(n),' a POSITIVE EPS CORRECTION'],'fig');
    
    if nv == 1
        rep1 = rep;
        eps = epsi;
        epsi = epsi-epsi/100;
        nv = nv+1;
        repm = rep;
        epsm = epsi;
    else
       if abs(rep) < abs(repm)
            repm = rep;
            epsm=epsi;
       end
       
       if sign(rep1) == sign(rep)
           if esc > 0 && abs(rep1) < abs(rep)
               break
           end
            eps=epsi;
            epsi = epsi-epst/100;
       else
            teps = (eps+epsi)/2;
            eps = epsi;
            epsi = teps;
            esc=1;
       end
       
       rep1 = rep;
        
    end
    
    iter = iter+1;
    if iter> 200
        break
    end
end

%calculation with minimum error
epsi = epsm;
[X, L, epsilon, lambda] = constrained_POD(eig_func_norm',li,OG_nm,RD_nm,epsi);
Lam_til = X'*diag(lambda2(1:OG_nm))*X;

% Modified redued order model coefficients
L_til = X'*li*X;
C_til = X'*ci;
Q_til=zeros(RD_nm,RD_nm,RD_nm);
Q = reshape(qi,OG_nm,OG_nm,OG_nm);
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

in=1; % initial point (to)

reduced_model_coeff= ode_coefficients(RD_nm,RD_nm,fc_til);

%Solution of the system of equation
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
%(Control implicit in the modes)

% Again look into tspan
tspan = [0 100];
[T1,Y1] = ode113(@(t, y) system_odes(t, y, -reduced_model_coeff)...
    , tspan ,eig_func_norm(in,1:RD_nm), options);   %(base line)(f = 500 Hz)

rep = error_til(Lam_til,Y1);
  
figure
  
for i = 1:2
        subplot(2,2,i)%(211)%
        if RD_nm>4
            ep=4;
        else
            ep =RD_nm;
        end
        hold on
%         if vd > Nsamples
            plot((0:size(eig_func_norm(in:end,1))-1)*sc/mu0,eig_func_norm(in:end,i),'. k')
            %                 plot((0:size(psitc1(in:end,1))-1)*sc/mu0,psitc1(in:end,1:ep),'.g')
%         else
%             plot((0:vd-1)*sc/mu0,eig_func_norm(in:vd,i),'. k')
%             %                 plot((0:size(psitc1,1)-1)*sc/mu0,psitc1(:,1:ep),'.g')
%         end
        
        plot(T1*sc/mu0,Y1(:,i),'m','MarkerSize',3);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)
        
        ylim([-1 1]);
        
        title(['a^',num2str(i), '(t); ' num2str(RD_nm) 'POD modes'],'fontname','times new roman','fontsize',11)%,'FontWeight','demi'
        xlabel('Time (s)','fontname','times new roman','fontsize',11)%,'units', 'norm','position',[.5,-.105])
        ylabel('Amplitude','fontname','times new roman','fontsize',11)
        
        if tn(vd)*sc/mu0 > 0.035
            xlim([0.03 0.035]);
        else
            xlim([0 tn(vd)*sc/mu0]);
        end

              subplot(2,2,2+i)%(212)%(223)%
              hold on
        for mn = 1
            run  edgarfouriercoeff;%_bl frequency
        end
        grid on
end
     %   saveas(gcf, ['M:\MASTERS\Spring 2014 Thesis figures\alpha',num2str(ia),'_Mod_POD_',num2str(N),' to POSITIVE EPS CORRECTION',num2str(n)],'fig');
end