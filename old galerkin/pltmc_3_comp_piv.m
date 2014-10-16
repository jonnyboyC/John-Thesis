% Routine to solve the system of equation forplotsvd1(data,phiutc(:,1:4),s,'u',sgtc,div,sec) the low dimensional model
% clear
% pack
%define normalization values, need to adjust
mu0=1;
sc=1;
re0=1140;
z{1}=ones(size(x{1}));
M0=.25;

tic
warning off MATLAB:ode45:IntegrationTolNotMet

%Solution of the time coefficient
sf=50000;
Ns=2048;
dl=mu0/(sc*sf);
tn=0:dl:dl*Ns*6*40;%size(psitc(:,1),1);%t0(in:end)-t0(in);
tm=zeros(10,size(tn,2));
acm=zeros(10,size(tn,2),4);
plm=zeros(10,size(tn,2));
pqm=zeros(10,size(tn,2));

fr=1000*sc/mu0;
tr=tn;

%  V=1.8*(sin(2*pi*fr/2*tr-pi/6)+0.5*sin(2*pi*fr*tr))*0.01239920704897;%+dt
%plot(tn(1),V(1),'+')
%  [cl,cq]=lqse_tc_sp_piv_mt(V,Aa1_pc11,Ba1_pc11,4);
% plot(tn(1:end-3),cq)
%  options = odeset('OutputFcn','odephas3');
% ci_i=lo_i.*niu_i'+qoo_i;%/re0
% ni_i=diag(niu_i)*ones(nm,nm);%;
% li_i=ni_i.*l_i+qo_i;
% fcuhi1=[ci_i li_i q_i];
%
% ci_i1=lo_i1.*niu_i1'+qoo_i1;%/re0
% ni_i1=diag(niu_i1)*ones(nm,nm);%;
% li_i1=ni_i1.*l_i1+qo_i1;
% fcuhi11=[ci_i1 li_i1 q_i1];
%
% ci_i2=lo_i2.*niu_i2'+qoo_i2;%/re0
% ni_i2=diag(niu_i2)*ones(nm,nm);%;
% li_i2=ni_i2.*l_i2+qo_i2;
% fcuhi12=[ci_i2 li_i2 q_i2];

for n =6 % number of modes to be used
    n
%     cm=4
    %     [Aa1_v1,Ba1_v1,a1l1v,a1q1v]=lqse_sp_piv_mt(vcn1,psc11(:,1:cm),Nsamples,1,dp);
    %    [Aa2_v1,Ba2_v1,a2l1v,a2q1v]=lqse_sp_piv_mt(vcn2,psc12(:,1:cm),Nsamples,1,dp);
    %     [Aam_v1,Bam_v1,aml1v,amq1v]=lqse_sp_piv_mt(vcn,psc1m(:,1:cm),1.5*Nsamples,1,dp);
    %      run cont_mod
    phiutca{1}=phiutc{1}(:,1:n);
    phivtca{1}=phivtc{1}(:,1:n);
%     phum{1}=phucm{1}(:,1:cm);
%     phvm{1}=phvcm{1}(:,1:cm);
    [lo_i,l_i,qoo_i,qo_i,q_i]=coef_visc(muz,mvz,x,y,phiutca,phivtca,t0,sa,re0,vt,z{1},M0, bi);
    niu_i=visc_dis_n(psitc,n,sgtc2,l_i,qo_i,q_i);
    ci_i=lo_i/re0.*niu_i'+qoo_i;%
    ni_i=diag(niu_i)+(ones(n,n)-eye(n))/re0;%;%diag(niu_i)*ones(n,n);
    li_i=ni_i.*l_i+qo_i;
    fcuhi1=[ci_i li_i q_i];

%     %Rempfer's approach
%     [lo_r,l_r,qoo_r,qo_r,q_r,lc_r,qom_r,qc_r,qm_r]=coef_visc_remp1(muz,mvz,data.x,data.y,phiutca,phivtca,phum,phvm,t0,sa,re0,vt,data.zg{1},M0);
%     [niu_r,niu_rc]=visc_dis_n_cont(psitc,psitc3,psc1m,n,sgtc2,sgtc32,l_r,qo_r,q_r,lc_r,qom_r,qm_r,qc_r);
%     ci_r=lo_r.*niu_r'+qoo_r;%/re0
%     ni_r=diag(niu_r);%+(ones(n,n)-eye(n))/re0;%diag(niu_r)*ones(n,n);
%     li_r=ni_r.*l_r+qo_r;
%     ni_rc=diag(niu_rc);%+(ones(n,n)-eye(n))/re0;%diag(niu_rc)*ones(n,cm)
%     li_rc=lc_r/re0+qom_r;%ni_rc.*
%     fcuhr1=[ci_r li_r q_r];
%     fcuhrc=[li_rc qm_r qc_r];
%     c1=mean(psc1m(:,1));
% 
%     F0=ci_r;
%     G0=li_r;
%     H0=q_r;
% 
%     P0=li_rc;
%     Q0=qm_r;
%     R0=qc_r;
    %rempfer's approach adding the mean flow difference
    %     [lo_rm,l_rm,qoo_rm,qo_rm,q_rm,lc_rm,qom_rm,qc_rm,qm_rm,lco_rm,qoc_rm,qcom_rm,qco_rm]=coef_visc_remp_m(muz,mvz,muz1,mvz1,data.x,data.y,phiutca,phivtca,phuc1,phvc1,t0,sa,re0,vt,data.zg{1},M0);
    %     [niu_rm,niu_rcm]=visc_dis_n_cont(psitc,psitc1,psc11,n,sgtc2,sgtc12,l_rm,qo_rm,q_rm,lc_rm,qom_rm,qm_rm,qc_rm);
    %     ci_rm=lo_rm.*niu_rm'+qoo_rm;%/re0
    %     ni_rm=diag(niu_rm)*ones(n,n);%;
    %     li_rm=ni_rm.*l_rm+qo_rm;
    %     ni_rcm=diag(niu_rcm)*ones(n,n);
    %     li_rcm=ni_rcm.*lc_rm+qom_rm;
    %     lio_rcm=niu_rcm'.*lco_rm+qoc_rm;
    %     fcuhr1m=[ci_rm li_rm q_rm];
    %     fcuhrcm=[li_rcm qm_rm,qc_rm];
    %     fduc=[lio_rcm qcom_rm qco_rm];
    %     c1=mean(psc11(:,1));

%     % Onder's
%     [lo_e,l_e,qoo_e,qo_e,q_e,P0_e,Q0_e]=coef_visc_incomp_ond(muz,mvz,data.x,data.y,phiutca,phivtca,t0,sa,re0,vt,data.zg{1},M0);
%     niu_e=visc_dis_n(psitc,n,sgtc2,l_e,qo_e,q_e);
%     F0_e=lo_e.*niu_e'+qoo_e;%/re0
%     ni_e=diag(niu_e);%+(ones(n,n)-eye(n))/re0;%diag(niu_e)*ones(n,n)
%     G0_e=ni_e.*l_e+qo_e;
%     H0_e=q_e;

    % Weak
    %     [lo_w,l_w,qoo_w,qo_w,q_w,P0_w]=coef_visc_incomp_ws(muz,mvz,data.xg,data.yg,phiutca,phivtca,t0,sa,re0,vt,data.zg{1},M0);
    %     niu_w=visc_dis_n(psitc,n,sgtc2,l_w,qo_w,q_w);
    %     F0_w=lo_w.*niu_w'+qoo_w;%/re0
    %     ni_w=diag(niu_w)*ones(n,n);
    %     G0_w=ni_w.*l_w+qo_w;
    %     H0_w=q_w;
    %     P0_w=P0_w.*(diag(niu_w)*ones(n,2));

    %     phiutcb{1}=phiutc4{1}(:,1:n);
    %     phivtcb{1}=phivtc4{1}(:,1:n);
    %             [lo_i1,l_i1,qoo_i1,qo_i1,q_i1]=coef_visc(muz4,mvz4,data.x,data.y,phiutcb,phivtcb,t0,sa,re0,vt,data.zg{1},M0);
    %             niu_i1=visc_dis_n(psitc4,n,sgtc42,l_i1,qo_i1,q_i1);
    %             ci_i1=lo_i1/re0+qoo_i1;%.*niu_i1'
    %             ni_i1=diag(niu_i1)+(ones(n,n)-eye(n))/re0;%diag(niu_i1)*ones(n,n);%;
    %             li_i1=ni_i1.*l_i1+qo_i1;
    %             fcuhi11=[ci_i1 li_i1 q_i1];

    % Onder's
    %     [lo_e1,l_e1,qoo_e1,qo_e1,q_e1,P0_e1,Q0_e1]=coef_visc_incomp_ond(muz4,mvz4,data.x,data.y,phiutcb,phivtcb,t0,sa,re0,vt,data.zg{1},M0);
    %     niu_e1=visc_dis_n(psitc4,n,sgtc42,l_e1,qo_e1,q_e1);
    %     F0_e1=lo_e1/re0+qoo_e1;%.*niu_e1'
    %     ni_e1=diag(niu_e1)+(ones(n,n)-eye(n))/re0;%diag(niu_e1)*ones(n,n);%
    %     G0_e1=ni_e1.*l_e1+qo_e1;
    %     H0_e1=q_e1;
    % c2=mean(psitc2(:,1));
    % Weak
    %     [lo_w1,l_w1,qoo_w1,qo_w1,q_w1,P0_w1]=coef_visc_incomp_ws(muz2,mvz2,data.xg,data.yg,phiutcb,phivtcb,t0,sa,re0,vt,data.zg{1},M0);
    %     niu_w1=visc_dis_n(psitc2,n,sgtc22,l_w1,qo_w1,q_w1);
    %     F0_w1=lo_w1.*niu_w1'+qoo_w1;%/re0
    %     ni_w1=diag(niu_w1)*ones(n,n);
    %     G0_w1=ni_w1.*l_w1+qo_w1;
    %     H0_w1=q_w1;
    %     P0_w1=P0_w1.*(diag(niu_w1)*ones(n,2));

    %     phiutcc{1}=[phiutc5{1}(:,1:n)];% phiutc2{1}(:,1)];+1
    %     phivtcc{1}=[phivtc5{1}(:,1:n)];% phivtc2{1}(:,1)];+1
    %     psitc2a=[psitc5(:,1:n)];% psitc2(:,1)];+1
    %     sgtc22a=[sgtc25];%(2:n+1); sgtc22(1); sgtc22(n+2:end)];
    %         [lo_i2,l_i2,qoo_i2,qo_i2,q_i2]=coef_visc(muz3,mvz3,data.x,data.y,phiutcc,phivtcc,t0,sa,re0,vt,data.zg{1},M0);
    %         niu_i2=visc_dis_n(psitc2a,n,sgtc22a,l_i2,qo_i2,q_i2);
    %         ci_i2=lo_i2.*niu_i2'+qoo_i2;%/re0
    %         ni_i2=diag(niu_i2)*ones(n,n);%;
    %         li_i2=ni_i2.*l_i2+qo_i2;
    %         fcuhi12=[ci_i2 li_i2 q_i2];
    %
    %     % Onder
    %     [lo_e2,l_e2,qoo_e2,qo_e2,q_e2,P0_e2,Q0_e2]=coef_visc_incomp_ond(muz3,mvz3,data.x,data.y,phiutcc,phivtcc,t0,sa,re0,vt,data.zg{1},M0);
    %     niu_e2=visc_dis_n(psitc2a,n,sgtc22a,l_e2,qo_e2,q_e2);
    %     F0_e2=lo_e2.*niu_e2'+qoo_e2;%+1/re0
    %     ni_e2=diag(niu_e2)*ones(n,n);%+1+1
    %     G0_e2=ni_e2.*l_e2+qo_e2;
    %     H0_e2=q_e2;
    %
    %     % Weak
    %     [lo_w2,l_w2,qoo_w2,qo_w2,q_w2,P0_w2]=coef_visc_incomp_ws(muz3,mvz3,data.xg,data.yg,phiutcc,phivtcc,t0,sa,re0,vt,data.zg{1},M0);
    %     niu_w2=visc_dis_n(psitc2a,n,sgtc22a,l_w2,qo_w2,q_w2);
    %     F0_w2=lo_w2.*niu_w2'+qoo_w2;%/re0
    %     ni_w2=diag(niu_w2)*ones(n,n);
    %     G0_w2=ni_w2.*l_w2+qo_w2;
    %     H0_w2=q_w2;
    %     P0_w2=P0_w2.*(diag(niu_w2)*ones(n,2));

    %Selecting only the required coefficients according to the number of
    %modes ( control implicit in the POD modes)
%     rmc=-odecoe(n,n,fcuhr1);
    %     rmc1=-odecoe(n,n,fcuhi11);
    %         rmc2=-odecoe(n,n,fcuhi12);

    %     fco0=-[F0_e G0_e H0_e];
    %     fco1=-[F0_e1 G0_e1 H0_e1];
    %     fco2=-[F0_e2 G0_e2 H0_e2];

    %     fcw0=-[F0_w G0_w H0_w];
    %     fcw1=-[F0_w1 G0_w1 H0_w1];
    %     fcw2=-[F0_w2 G0_w2 H0_w2];

        rmc0=odecoe(n,n,fcuhi1);
    %     rmc1=odecoe(n,n,fcuhi11);
    %     rmc2=odecoe(n,n,fcuhi12);

    %      rmc0(:,2:end)= rmc0(:,2:end).*(ones(n,size(rmc0,2)-1)-eye(n,size(rmc0,2)-1)*.5);
    %      rmc0(:,n+1:end)= rmc0(:,n+1:end)*6;
    %        rmc0(:,1)= rmc0(:,1)*-.4;
    %      rmc(:,1)= rmc(:,1)*-1;+1+1
    %      rmco0=-odecoe(n,n,fco0);
    %      rmco1=-odecoe(n,n,fco1);
    %     rmco2=-odecoe(n,n,fco2);

    %     rmcw0=-odecoe(n,n,fcw0);
    %     rmcw1=-odecoe(n,n,fcw1);
    %     rmcw2=-odecoe(n,n,fcw2);

    in=1; % initial point (to)
    %Selecting only the required coefficients according to the number of
    %modes ( explicit control scheme), coefficients of the control independet terms
    %         frmc0=odecoe(n,n,fco0);
    %         frmc1=odecoe(n,n,fco1);
    %         frmc2=odecoe(n,n,fco2);

    %     tn1=t1(1:end)-t1(1);
    %     tn2=t2(1:end)-t2(1);
    %     dt=.19921;
    %Solution of the system of equation
    options = odeset('RelTol',1e-6,'AbsTol',1e-8);
    %(Control implicit in the modes)
    %       [T0,Y0] = ode45(@eqdifvecm,tn,psitc(in,1:n),[],rmc);   %(base line)(in):dt:tn(end)
    [T1,Y1] = ode113(@eqdifvecm,[tn(1):dl:tn(end)],psitc(in,1:n),options,-rmc0);   %(base line)(f = 500 Hz)
    %     [T2,Y2] = ode113(@eqdifvecm,[tn(1):dl:tn(end)],psitc(in,1:n),options,-rmc1);   %(base line)(f = 900 Hz)
    %         [T3,Y3] = ode45(@eqdifvecm,tn,psitc3(in,1:n),options,-rmc2);   %(base line)(in):dt:tn(end)options
    %     [T1,Y1] = ode45(@eqdifvecm,tn,psitc(in,1:n),[],rmcw0);   %(base line)(f = 500 Hz)
    %     [T5,Y5] = ode45(@eqdifvecm,tn,psitc(in,1:n),[],rmcw2);   %(base line)(f = 900 Hz)

    %     sum([diag(Y3'*Y3)/size(Y3,1) sgtc(1:n)])
    %     sum([diag(Y1'*Y1)/size(Y1,1) sgtc1(1:n)])
    %     sum([diag(Y2'*Y2)/size(Y2,1) sgtc21(1:n)])
    %(explicit Control without the control term)
    % Onder's
    %     [T4,Y4] = ode45(@eqdifvecm,tn,psitc0(in,1:n),[],-frmc0);  %(base line)
    %     [T5,Y5] = ode45(@eqdifvecm,tn,psitc1(in,1:n),[],-frmc1);  %(f = 500 Hz)
    %     [T5,Y5] = ode45(@eqdifvecm,[tn2(in):dt:tn2(end)],psitc2(in,1:n),[],-frmc2);  %(f = 900 Hz)

    %Weak solution
    %     [T6,Y6] = ode45(@eqdifvecm,[tn(1):dt:tn(end)],psitc0(in,1:n),[],rmcw0);   %(base line)(in):dt:tn(end)
    %     [T7,Y7] = ode45(@eqdifvecm,[tn(1):dt:tn(end)],psitc1(in,1:n),[],rmcw1);   %(f = 500 Hz)
    %     [T8,Y8] = ode45(@eqdifvecm,[tn(1):dt:tn(end)],psitc2(in,1:n),[],rmcw2);   %(f = 900 Hz)
    %
    %(explicit Control with the control term included)
    % Onder's
    ec=0; %0 quadratic estimation, 1 linear estimation
    %       [T2,Y2] = ode45(@eqdifveci,tn,psitc(in,1:n),[],rmco0,P0_e(1:n,:),Q0_e(1:n,:,1:n),0,0*sc/mu0,0,1,mu0);  %(Base line)
    %     [T2,Y2] = ode45(@eqdifveci,tn,psitc(in,1:n),[],rmco1,P0_e1(1:n,:),Q0_e1(1:n,:,1:n),3.7/mu0,0*sc/mu0,0,1,mu0);  %(Base line)
    %     [T3,Y3] = ode45(@eqdifveci,tn,psitc1(in,1:n),[],rmco2,P0_e2(1:n,:),Q0_e2(1:n,:,1:n),0,0,0,1,mu0);  %(Base line)
%     [T2,Y2] = ode113(@eqdifveci_cont_r,tn,psitc(in,1:n),options,-rmc0,fcuhrc,0,0*sc/mu0,Aam_v1(:,1:cm),Bam_v1(:,1:cm),dl,ec);  %(Base line)
%     [T3,Y3] = ode113(@eqdifveci_cont_r,tn,psitc(in,1:n),options,-rmc0,fcuhrc,3.7,1610*sc/mu0,Aam_v1(:,1:cm),Bam_v1(:,1:cm),dl,ec);  %(Base line)
%     %     [T4,Y4] = ode45(@eqdifveci,tn,psitc1(in,1:n),[],rmco2,P0_e2(1:n,:),Q0_e2(1:n,:,1:n),500/mu0,3920*sc/mu0,0,1,mu0);  %(Base line)
%     [T4,Y4] = ode113(@eqdifveci_cont_r,tn,psitc(in,1:n),options,-rmc0,fcuhrc,3.7,3920*sc/mu0,Aam_v1(:,1:cm),Bam_v1(:,1:cm),dl,ec);  %(Base line)
    % T3 = T3a;
    % Y3 = Y3a;
    % T4 = T4a;
    % Y4 = Y4a;

    %     [T11,Y11] = ode45(@eqdifvecc,[tn2(in):dt:tn2(end)],psitc2(in,1:n),[],-frmc2,-P2(1:n,:),-Q2(1:n,:,1:n),40/mu2,900*sc/mu2,0,ci2,mu2); %(f = 900 Hz)

    %Weak solution
    %     [T1,Y1] = ode45(@eqdifveccs,tn,psitc(in,1:n),[],rmcw0,P0_w(1:n,:),0,0*sc/mu0,0,0,mu0);  %(Base line)(in):dt:tn(end)
    %     [T2,Y2] = ode45(@eqdifveccs,tn,psitc1(in,1:n),[],rmcw1,P0_w1(1:n,:),5/mu0,3920*sc/mu0,0,0,mu0);  %(Base line)(in):dt:tn(end)
    %     [T3,Y3] = ode45(@eqdifveccs,tn,psitc3(in,1:n),[],rmcw2,P0_w2(1:n,:),0,0*sc/mu0,0,0,mu0); %(f = 500 Hz)
    %     [T4,Y4] = ode45(@eqdifveccs,tn,psitc3(in,1:n),[],rmcw2,P0_w2(1:n,:),5/mu0,3920*sc/mu0,0,0,mu0); %(f = 900 Hz)
    %     % %
    Y1=Y1-ones(size(Y1,1),1)*mean(Y1);
%     Y2=Y2-ones(size(Y2,1),1)*mean(Y2);
%     Y3=Y3-ones(size(Y3,1),1)*mean(Y3);
%     Y4=Y4-ones(size(Y4,1),1)*mean(Y4);
    
    %Plotting the results of the first coefficient using the selected
    %numbers of modes ("n")
    vd=max([size(Y1,1)]);% size(Y2,1) size(Y3,1)]);
    %     if vd > 2047 %size(tn,2)/4
    figure('Name',[num2str(n)  ' modes'],'color','w')

    %(base line)
    if n< 4
        rm=n;
    else
        rm=4;
    end
    for i = 1:2
        subplot(2,2,i)%(211)%
        if n>4
            ep=4;
        else
            ep =n;
        end
        hold on
        if vd > Nsamples
            plot((0:size(psitc(in:end,1))-1)*sc/mu0,psitc(in:end,i),'. k')
            %                 plot((0:size(psitc1(in:end,1))-1)*sc/mu0,psitc1(in:end,1:ep),'.g')
        else
            plot((0:vd-1)*sc/mu0,psitc(in:vd,i),'. k')
            %                 plot((0:size(psitc1,1)-1)*sc/mu0,psitc1(:,1:ep),'.g')
        end
        %          plot(T0*sc/mu0,Y0(:,1:ep),'b','MarkerSize',3);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)

        plot(T1*sc/mu0,Y1(:,i),'m','MarkerSize',3);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)
%         plot(T2*sc/mu0,Y2(:,i),'b','MarkerSize',3);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)
%         plot(T3*sc/mu0,Y3(:,i),'r --','MarkerSize',3);
%         plot(T4*sc/mu0,Y4(:,i),'g -.','MarkerSize',3)
        %             plot(T6*sc/mu0,Y6(:,1:ep),'b-','MarkerSize',3)
        %         plot(T12,Y12(:,1),'m','MarkerSize',3)
        ylim([-1 1]);
        %
        %         legend('Org.','ODE','ODE v_d_i_s''Ond w/o cont.','Ond w. cont.')%,'ODE','weak w. cont.','qse','weak w/o cont.'
        title(['a^',num2str(i), '(t); ' num2str(n) 'POD modes'],'fontname','times new roman','fontsize',11)%,'FontWeight','demi'
        xlabel('Time (s)','fontname','times new roman','fontsize',11)%,'units', 'norm','position',[.5,-.105])
        ylabel('Amplitude','fontname','times new roman','fontsize',11)
        if tn(vd)*sc/mu0 > 0.035
            xlim([0.03 0.035]);
        else
            xlim([0 tn(vd)*sc/mu0]);
        end

        %             clear ue1 ve1
        %     ue1=phu(:,1:n)*Y3';
        %     ve1=phv(:,1:n)*Y3';
        % Estimated pressure
        %         if n>6
        %             [Epl,pc_e]=pest_piv(Y3(:,1:6),Ap0_a,Bp0_a,size(Y3(:,1:6),2),mp0,6);
        %             [Epl,pc_e1]=pest_piv(Y2(:,1:6),Ap1_a,Bp1_a,size(Y2(:,1:6),1),mp1,6);
        %             [Epl,pc_e2]=pest_piv(Y3(:,1:6),Ap0_a,Bp0_a,size(Y3(:,1:6),1),mp2,6);
        %             [Epl,pc_e3]=pest_piv(Y4(:,1:6),Ap1_a,Bp1_a,size(Y4(:,1:6),1),mp2,6);
        %             %         pc_e=Ap0'*Y3(:,1:6)';
        %         else
        %             [Ap0_a,Bp0_a,Epl0,Epq0,mp0]=lqsep_piv_al(pcn,psitc(:,1:n),Nsamples,t0,n);
        %             [Ap1_a,Bp1_a,Epl1,Epq1,mp1]=lqsep_piv_al(pcn1,psitc1(:,1:n),Nsamples,t0,n);
        %             [Ap2_a,Bp2_a,Epl2,Epq2,mp2]=lqsep_piv_al(pcn2 ,psitc2(:,1:n),Nsamples,t0,n);
        %             [Epl,pc_e]=pest_piv(Y1,Ap0_a,Bp0_a,size(Y1,1),mp0,n);
        %             [Epl,pc_e1]=pest_piv(Y2,Ap1_a,Bp1_a,size(Y2,1),mp1,n);
        %             [Epl,pc_e2]=pest_piv(Y3,Ap0_a,Bp0_a,size(Y3,1),mp2,n);
        %             [Epl,pc_e3]=pest_piv(Y4,Ap1_a,Bp1_a,size(Y4,1),mp2,n);
        %             %         pc_e=Ap0(1:n,:)'*Y3(:,1:n)';
        %         end
        % estimated time coeffcient
        %             [Eal0,Eaq]=lqse_tc_sp_piv(pc_e,Aa_p1,Ba_p1);
        %             plot(T3*sc/mu0,Eaq(1:ep,:)','b','MarkerSize',3);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)

        subplot(2,2,2+i)%(212)%(223)%
        %             plot(T3*sc/mu0,pc_e)
        %             hold on
        %             %             plot(T3*sc/mu0,Epl,'r-')
        %             plot((0:1:255)*1/50000,pi(5,1:256)','k-')
        %             xlim([0 .01])
        %             %             ylim([-.05 .05])
        %             subplot(224)
        for mn = 1:4
            run  edgarfouriercoeff;%_bl frequency
        end
        grid on
    end
    %     figure
    %     run fft_press
    % %          hold on
    % %          plot(tn(1:end),psitc1(in:end,i),'k')
    % plot3(Y3(:,1),Y3(:,2),Y3(:,3));
    %           plot(Y3(:,1),Y3(:,2),'b','MarkerSize',2);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)
    % %          plot(T4,Y4(:,i),'r','MarkerSize',3);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)
    % %         %         plot(T5,Y5(:,i),'ro-','MarkerSize',3)
    % %         plot(T7,Y7(:,i),'cd-','MarkerSize',3)
    % %         %         plot(T13,Y13(:,1),'m','MarkerSize',3)
    %          ylim([-2.5 2.5]);
    % %
    % %         %         legend('Org.','ODE','Ond w/o cont.','Ond w. cont.')%'ODE v_d_i_s','''weak w. cont.',qse','weak w/o cont.'
    %          title('500 Hz ','fontname','times new roman','fontsize',12)%,'FontWeight','demi'
    % %         %     xlabel('Time','fontname','times new roman','fontsize',11)%,'units', 'norm','position',[.5,-.105])
    %          ylabel('a^1(t)','fontname','times new roman','fontsize',11)
    % %         % end
    %          subplot(313)%(2,2,2*i)
    %          hold on
    %          plot(tn(1:end),psitc2(in:end,i),'k')
    %          plot(T2,Y2(:,i),'b','MarkerSize',2);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)
    %          plot(T5,Y5(:,i),'r','MarkerSize',3);%'color',[ 0.000 0.502 0.000 ],'Marker','o','MarkerSize',4)
    % %         %         plot(T5,Y5(:,i),'ro-','MarkerSize',3)
    % %         plot(T8,Y8(:,i),'cd-','MarkerSize',3)
    % %         %         plot(T13,Y13(:,1),'m','MarkerSize',3)
    %          ylim([-2.5 2.5]);
    % %
    % %         %                 legend('Org.','ODE','Ond w/o cont.','Ond w. cont.')%'ODE v_d_i_s','''weak w. cont.',qse','weak w/o cont.'
    %          title('900 Hz ','fontname','times new roman','fontsize',12)%,'FontWeight','demi'
    % %         %     xlabel('Time','fontname','times new roman','fontsize',11)%,'units', 'norm','position',[.5,-.105])
    %          ylabel('a^1(t)','fontname','times new roman','fontsize',11)
    %     end
    %         clear ue1 ve1
    %         ue1=phu(:,1:n)*Y1';
    %         ve1=phv(:,1:n)*Y1';
    % tm(n-2,1:size(T3,1))=T3;
    % acm(n-2,1:size(Y3,1),1:size(Y3,2))=Y3;
    % plm(n-2,1:size(Epl,2))=Epl(5,:);
    % pqm(n-2,1:size(pc_e,2))=pc_e(5,:);
    clear T0 T1 T2 T3 T4 T5 T6 T7 T8
    %     clear Y0 Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8
    %     clear rmc0 rmc1 rmc0_1 rmc1_1% rmc2
    % clear frmc0 frmc1 frmc2
end
toc