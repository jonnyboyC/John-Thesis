function [lo,l,qoo,qo,q]=coef_visc(uu,vv,xg,yg,phuz,phvz,t,sa,re,vtz,z,M,bi)
nzones=size(vtz,2);
at=0;
clt=0;
cqt=0;
cbt=0;
cct=0;
cdt=0;

for nz=1:nzones
    s=sa{nz};
    x=xg{nz};
    y=yg{nz};
    vt=vtz{nz};
    nm=size(phuz{nz},2);
    nx=size(x,1);
    ny=size(y,2);
    
%     [xxi,yxi,xet,yet,aj]=metric1(x,y);
     [xxi,yxi,xet,yet,aj]=metric(x,y);
    %mean velocities
    mu=uu{nz};
    mv=vv{nz};
    
    phu=phuz{nz};
    phv=phvz{nz};
        
    % clear uu vv aa
    % pack

    [udx,ud2x,udy,ud2y] = derdxmt(mu,x,y,s,ones(nx,ny),xxi,yxi,xet,yet,aj,bi);
    [vdx,vd2x,vdy,vd2y] = derdxmt(mv,x,y,s,ones(nx,ny),xxi,yxi,xet,yet,aj,bi);
    
    [phudx,phud2x,phudy,phud2y] = derdxmt(phu,x,y,s,ones(nx,ny),xxi,yxi,xet,yet,aj,bi);
    [phvdx,phvd2x,phvdy,phvd2y] = derdxmt(phv,x,y,s,ones(nx,ny),xxi,yxi,xet,yet,aj,bi);
   
%     [udx,ud2x,udy,ud2y] = derdxmt1(mu,x,y,s,ones(nx,ny),xxi,yxi,xet,yet,aj);%derdx1(mu,x,y,s);
%     [vdx,vd2x,vdy,vd2y] = derdxmt1(mv,x,y,s,ones(nx,ny),xxi,yxi,xet,yet,aj);%derdx1(mv,x,y,s);
%     
%     [phudx,phud2x,phudy,phud2y] = derdxmt1(phu,x,y,s,ones(nx,ny),xxi,yxi,xet,yet,aj);%derdx1(phu,x,y,s);
%     [phvdx,phvd2x,phvdy,phvd2y] = derdxmt1(phv,x,y,s,ones(nx,ny),xxi,yxi,xet,yet,aj);%derdx1(phv,x,y,s);  

    udx=reshape(udx,nx*ny,1);
    udy=reshape(udy,nx*ny,1);
    ud2x=reshape(ud2x,nx*ny,1);
    ud2y=reshape(ud2y,nx*ny,1);
    vdx=reshape(vdx,nx*ny,1);
    vdy=reshape(vdy,nx*ny,1);
    vd2x=reshape(vd2x,nx*ny,1);
    vd2y=reshape(vd2y,nx*ny,1);
    
    phudx=reshape(phudx,nx*ny,nm);
    phudy=reshape(phudy,nx*ny,nm);
    phud2x=reshape(phud2x,nx*ny,nm);
    phud2y=reshape(phud2y,nx*ny,nm);
    phvdx=reshape(phvdx,nx*ny,nm);
    phvdy=reshape(phvdy,nx*ny,nm);
    phvd2x=reshape(phvd2x,nx*ny,nm);
    phvd2y=reshape(phvd2y,nx*ny,nm);
    
    d2u=(ud2x+ud2y);
    d2v=(vd2x+vd2y);
    d2phu=(phud2x+phud2y);
    d2phv=(phvd2x+phvd2y);
    clear ud2x ud2y vd2x vd2y phud2x phud2y phvd2x phvd2y
   
    
    %coefficient for the system of ODE.
    %left side coefficient
    au=inerpro(phu,phu,t,xg,yg,[1 1:nx-1],vt);
    av=inerpro(phv,phv,t,xg,yg,[1 1:nx-1],vt);
    at=at+au+av;
    
    % Constant Terms
    lu= d2u;
    lv=d2v;
    
    bu=d2phu;
    bv=d2phv;
    
    clu = inerpro(lu,phu,t,xg,yg,[1 1:nx-1],vt);
    clv = inerpro(lv,phv,t,xg,yg,[1 1:nx-1],vt);
    clt=clt+clu+clv;
    
    cbu = inerpro(bu,phu,t,xg,yg,[1 1:nx-1],vt);
    cbv = inerpro(bv,phv,t,xg,yg,[1 1:nx-1],vt);
    cbt=cbt+cbu+cbv;
    
    qu=mu.*udx+mv.*udy;
    qv=mu.*vdx+mv.*vdy;
    
    cqu = inerpro(qu,phu,t,xg,yg,[1 1:nx-1],vt);
    cqv = inerpro(qv,phv,t,xg,yg,[1 1:nx-1],vt);
    cqt=cqt+cqu+cqv;
    
    % Linear Terms
    uphux=mu*ones(1,nm).*phudx;
    uphvx=mu*ones(1,nm).*phvdx;
    uxphu=phu.*(udx*ones(1,nm));
    vxphu=phu.*(vdx*ones(1,nm));
    vphuy=mv*ones(1,nm).*phudy;
    vphvy=mv*ones(1,nm).*phvdy;
    uyphv=phv.*(udy*ones(1,nm));
    vyphv=phv.*(vdy*ones(1,nm));
    
    cu=uphux+uxphu+vphuy+uyphv;
    cv=uphvx+vxphu+vphvy+vyphv;
    
    ccu = inerpro(cu,phu,t,xg,yg,[1 1:nx-1],vt);
    ccv = inerpro(cv,phv,t,xg,yg,[1 1:nx-1],vt);
    cct=cct+ccu+ccv;
    
    
    % Quadratic Terms
    for k=1:nm
        phuphux(:,:,k)=(phu(:,k)*ones(1,nm)).*phudx;
        phvphuy(:,:,k)=(phv(:,k)*ones(1,nm)).*phudy;
        phuphvx(:,:,k)=(phu(:,k)*ones(1,nm)).*phvdx;
        phvphvy(:,:,k)=(phv(:,k)*ones(1,nm)).*phvdy;
    end
    
    for j=1:nm
        cdu(:,:,j) = inerpro(phuphux(:,:,j)+phvphuy(:,:,j),phu,t,xg,yg,[1 1:nx-1],vt);
        cdv(:,:,j) = inerpro(phuphvx(:,:,j)+phvphvy(:,:,j),phv,t,xg,yg,[1 1:nx-1],vt);
    end
    
    clear phuphux phvphuy phuphvx phvphvy
    
    cdu=reshape(cdu,nm,nm*nm);
    cdv=reshape(cdv,nm,nm*nm);
    cdt=cdt+cdu+cdv;
    
    clear cdu cdv
end
% need to verify the inner products and the summation
% Final Coefficient for the system of ODE
lo = clt;
l = cbt;
qoo = -cqt;
qo = -cct;
q = -cdt;

%at=[diag(au);diag(av); diag(ah)]*ones(1,size(fc,2));
%fc=fc./at;



