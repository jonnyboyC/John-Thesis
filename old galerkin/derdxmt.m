function [dx,d2x,dy,d2y] = derdxmt(p,x,y,s,z,xxi,yxi,xet,yet,aj,bi)
% [xix,xiy,etx,ety]=metric(x,y);
nm=size(p,2);
nx=size(x,1);
ny=size(y,2);
% z1=z;
% z1(97:117,73:87)=0;
% z1(:,1:3)=0;
% z1(13:22,73:87)=0;
for j=1:nm
    pp{1} =  reshape(p(:,j),nx,ny);%vect2blocks(p(:,j),s,1);
    [dxic,detc]=visder(pp{1},x,y,z,0,bi);

    [dx(:,:,j),dy(:,:,j)]=dervt(dxic,detc,xxi,yxi,xet,yet,aj,x,y);
    [dxic,detc]=visder(dx(:,:,j),x,y,z,1,bi);
    [d2x(:,:,j),dyb]=dervt(dxic,detc,xxi,yxi,xet,yet,aj,x,y);
    [dxic,detc]=visder(dy(:,:,j),x,y,z,1,bi);
    [dxb,d2y(:,:,j)]=dervt(dxic,detc,xxi,yxi,xet,yet,aj,x,y);
% dx(:,:,j)=dx(:,:,j).*z1;
% d2x(:,:,j)=d2x(:,:,j).*z1;
% dy(:,:,j)=dy(:,:,j).*z1;
% d2y(:,:,j)=d2y(:,:,j).*z1;
end
