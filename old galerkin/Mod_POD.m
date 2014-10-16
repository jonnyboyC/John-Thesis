% script to canculate new POD modes
% global L epsilon lambda
clear global L  lambda epsilon
tic
%L;L = li_i;

N = 16; %counter; %from plmtc_3_comp_piv_bl
n = 8;  %counter2;

rep = 1;
epsi = sum(li_i*sgtc2(1:N));%-.000001; %.005;%
nv = 1;
ct = 0;

esc=0;
epst = epsi*1.01;
figure
hold on
while abs(epst) > abs(epsi) %500 %(rep ~= 0 ||)ct < 1 % 
    [X] = constrained_POD(psitc',li_i,N,n,epsi);
    %inputs of constrained_POD are the POD temporal coefficients,psitc',
    %the linear Galerkin matrix, li_i, the transformation dimensions N and
    %n, and the transfer term parameter epsi
    Lam_til = X'*diag(sgtc2(1:N))*X;
    
    % Modified reduced order model coefficients
    L_til = X'*li_i*X;
    C_til = X'*ci_i;
    Q_til=zeros(n,n,n);
    Q = reshape(q_i,N,N,N);
    for i = 1:n
        for j = 1:n
            for k = 1:n
                for p = 1:N
                    for r = 1:N
                        for q = 1:N
                            Q_til(i,j,k) = Q_til(i,j,k) + X(p,i)*Q(p,q,r)*X(q,j)*X(r,k);
                        end
                    end
                end
            end
        end
    end
    % solve new system of equation
    Q_til=reshape(Q_til,n,n*n);
    fc_til=[C_til L_til Q_til]; % constant, linear and quadratic terms coefficients
    
    in=1; % initial point (to)
    
    rmct=-odecoe(n,n,fc_til);
    
    %Solution of the system of equation
    options = odeset('RelTol',1e-6,'AbsTol',1e-8);
    %(Control implicit in the modes)
    
    [T_til,Y_til] = ode113(@eqdifvecm,[tn(1):dl:tn(2000)],psitc(in,1:n),options,-rmct);   %(base line)(f = 500 Hz)
    
    [rep] = error_til(Lam_til,Y_til);
    plot(epsi,rep,'*')
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
    
    ct = ct+1;
    if ct> 200
        break
    end
end

%calculation with minimum error
epsi = epsm;
[X] = constrained_POD(psitc',li_i,N,n,epsi);
Lam_til = X'*diag(sgtc2(1:N))*X;

% Modified redued order model coefficients
L_til = X'*li_i*X;
C_til = X'*ci_i;
Q_til=zeros(n,n,n);
Q = reshape(q_i,N,N,N);
for i = 1:n
    for j = 1:n
        for k = 1:n
            for p = 1:N
                for r = 1:N
                    for q = 1:N
                        Q_til(i,j,k) = Q_til(i,j,k) + X(p,i)*Q(p,q,r)*X(q,j)*X(r,k);
                    end
                end
            end
        end
    end
end

% solve new system of equation
Q_til=reshape(Q_til,n,n*n);
fc_til=[C_til L_til Q_til]; % constant, linear and quadratic terms coefficients

in=1; % initial point (to)

rmct=-odecoe(n,n,fc_til);

%Solution of the system of equation
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
%(Control implicit in the modes)

[T1,Y1] = ode113(@eqdifvecm,[tn(1):dl:tn(4*2048)],psitc(in,1:n),options,-rmct);   %(base line)(f = 500 Hz)

[rep] = error_til(Lam_til,Y1);
  
figure
  
for i = 1:2
        subplot(2,2,i)%(211)%
        if n>4
            ep=4;
        else
            ep =n;
        end
        hold on
%         if vd > Nsamples
            plot((0:size(psitc(in:end,1))-1)*sc/mu0,psitc(in:end,i),'. k')
            %                 plot((0:size(psitc1(in:end,1))-1)*sc/mu0,psitc1(in:end,1:ep),'.g')
%         else
%             plot((0:vd-1)*sc/mu0,psitc(in:vd,i),'. k')
%             %                 plot((0:size(psitc1,1)-1)*sc/mu0,psitc1(:,1:ep),'.g')
%         end
        
        plot(T1*sc/mu0,Y1(:,i),'m','MarkerSize',3);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)
        
        ylim([-1 1]);
        
        title(['a^',num2str(i), '(t); ' num2str(n) 'POD modes'],'fontname','times new roman','fontsize',11)%,'FontWeight','demi'
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
toc