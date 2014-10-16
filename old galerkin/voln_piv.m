function vt=voln_piv(x,y,nz,mu,bi)
% finding boundaries of image
s=[size(x,1); size(x,2)];
cnt=1;
mu=reshape(mu,s(1),s(2));
vt=ones(size(x));

% [ X, Y, xb, yb, ib, jb, bi ] = boundary_check(x, y, mu);


% % finding boundary points
l=size(x);
% b1=1;
% b2=l(1);

% return
% if min(min(y)) < 0
%     b3=l(2);
% else
%     b3=1;
% end
%
% if x(end,1)>0
%     for i=1:l(1)-1
%         if (x(i,1) <0) & (x(i+1,1) > 0)
%             b1=i;                   % Front wall of the cavity
%         elseif (x(i,1) < 4) & (x(1+i,1) > 4)
%             b2=i;               % Rear wall of the cavity
%         end
%     end
% else
%     b1=1;                   % Front wall of the cavity
%     b2=l(1);
% end
%
% for i=1:l(2)
%     if y(1,i)<0 & y(1,i+1)>0
%         b3=i;                   % Floor of the entrance section
%     end
% end

% for j=1:nz
n2=size(x,2);
np=size(x,1);
hl=zeros(np,n2);
hr=zeros(np,n2);
ht=zeros(np,n2);
hb=zeros(np,n2);

% area elements
hl(2:np,:)=(x(2:end,:)-x(1:end-1,:));
hr(1:np-1,:)=(x(2:end,:)-x(1:end-1,:));
hb(:,2:end)=(y(:,2:end)-y(:,1:end-1));
ht(:,1:end-1)=(y(:,2:end)-y(:,1:end-1));

for r=1:np-1
    for q = 1:n2-1
        xm(r,q)=mean([x(r,q) x(r,q+1) x(r+1,q) x(r+1,q+1)]);
        ym(r,q)=mean([y(r,q) y(r,q+1) y(r+1,q) y(r+1,q+1)]);
    end
end
cnt = 1;
% ld=max(max(ib));

for i = 1:size(x,1)
    for j=1:size(x,2)
        if  bi(i,j)== -1
            vt(i,j)=0;
                 
        elseif (i == 1 && (j == size(x,2)|| bi(i,j) == 0 && bi(i+1,j) == 0 && bi(i+1,j-1) > 0))  % Top Left Corner
            wd=(xm(i,j-1)-x(i,j));
            lh=y(i,j)-ym(i,j-1);
            vt(i,j)=lh*wd;
        
        elseif (i == size(x,1) && (j == size(x,2)|| bi(i,j) == 0 && bi(i-1,j) == 0 && bi(i-1,j-1) > 0)) % Top Right Corner
            wd=(x(i,j)-xm(i-1,j-1));
            lh=y(i,j)-ym(i-1,j-1);
            vt(i,j)=lh*wd;
        
        elseif ( i == 1 && (j == 1|| bi(i,j) == 0 && bi(i+1,j) == 0 && bi(i+1,j+1) > 0)) % Bottom Left Corner
            wd=(xm(i,j)-x(i,j));
            lh=ym(i,j)-y(i,j);
            vt(i,j)=lh*wd;
        
        elseif (i == size(x,1) && (j == 1 || bi(i,j) == 0 && bi(i-1,j) == 0 && bi(i-1,j+1) > 0)) % Bottom Right Corner
            wd=(x(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            vt(i,j)=lh*wd;
        
        elseif (i < size(x,1) && j > 1 && bi(i,j) == 0 &&  bi(i+1,j) == 0 && bi(i,j-1) == 0 && bi(i+1,j-1) > 0) % Upper Left Corner
            wd=(xm(i,j-1)-x(i,j));
            lh=y(i,j)-ym(i,j-1);
            vt(i,j)=lh*wd;
        
        elseif (i > 1 && j > 1 && bi(i,j) == 0 &&  bi(i-1,j) == 0 && bi(i,j-1) == 0 && bi(i-1,j-1) > 0) % Upper Right Corner
            wd=(x(i,j)-xm(i-1,j-1));
            lh=y(i,j)-ym(i-1,j-1);
            vt(i,j)=lh*wd;
        
        elseif (i < size(x,1) && j < size(x,2) && bi(i,j) == 0 &&  bi(i+1,j) == 0 && bi(i,j+1) == 0 && bi(i+1,j+1) > 0) % Lower Left Corner
            wd=(xm(i,j)-x(i,j));
            lh=ym(i,j)-y(i,j);
            vt(i,j)=lh*wd;
        
        elseif (i > 1 && j< size(x,2)  && bi(i,j)== 0 &&  bi(i-1,j)== 0 && bi(i,j+1)== 0 && bi(i-1,j+1)> 0) % lower Right Corner
            wd=(x(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            vt(i,j)=lh*wd;
            
        elseif (j == 1 || i > 1 && i < size(x,1) && j > 1 && bi(i,j) == 0 &&  bi(i-1,j) == 0 && bi(i+1,j) == 0 &&  bi(i,j+1) > 0) % Bottom Boundary
            wd=(xm(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            rh=ym(i,j)-y(i,j);
            vt(i,j)=(rh+lh)/2*wd;
        
        elseif (i == 1 || i > 1 && j > 1 && j < size(x,2)  &&  bi(i,j) == 0 &&  bi(i,j-1) == 0 && bi(i,j+1) == 0 && bi(i+1,j) > 0) % Left Boundary
            wd=(xm(i,j)-x(i,j));
            lh=ym(i,j)-ym(i,j-1);
            vt(i,j)=lh*wd;
        
        elseif (i == size(x,1) || i < size(x,1) && j > 1 && j < size(x,2) && bi(i,j) == 0 &&  bi(i,j+1) == 0 && bi(i,j-1) == 0 && bi(i-1,j)> 0) % Right Boundary
            wd=(x(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-ym(i-1,j-1);
            vt(i,j)=lh*wd;
        
        elseif (j==size(x,2) || i > 1 && i < size(x,1) && j < size(x,2) && bi(i,j) == 0 &&  bi(i-1,j) == 0 && bi(i+1,j) == 0 && bi(i,j-1)> 0) % Top Boundary
            wd=(xm(i,j-1)-xm(i-1,j-1));
            lh=y(i,j)-ym(i-1,j-1);
            rh=y(i,j)-ym(i,j-1);
            vt(i,j)=(rh+lh)/2*wd;
            
        elseif (i< size(x,1) && j > 1 && bi(i,j) == 0 &&  bi(i+1,j) == 0 && bi(i,j-1) == 0 && bi(i-1,j+1) > 0) % Upper Left Edge
            wd=(xm(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            rh=ym(i,j)-y(i,j);
            vt(i,j)=(rh+lh)/2*wd;
            wd=(x(i,j)-xm(i-1,j-1));
            lh=(y(i,j)-ym(i-1,j-1));
            vt(i,j)=vt(i,j)+lh*wd;
            
        elseif (i > 1 && j > 1 && bi(i,j) == 0 &&  bi(i-1,j) == 0 && bi(i,j-1) == 0 && bi(i+1,j+1) > 0) % Upper Right Edge
            wd=(xm(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            rh=ym(i,j)-y(i,j);
            vt(i,j)=(rh+lh)/2*wd;
            wd=(xm(i,j-1)-x(i,j));
            rh=(y(i,j)-ym(i,j-1));
            vt(i,j)=vt(i,j)+rh*wd;
            
        elseif (i< size(x,1) && j < size(x,2) && bi(i,j) == 0 &&  bi(i+1,j) == 0 && bi(i,j+1) == 0 && bi(i-1,j-1) > 0) % Lower Left Edge
            wd=(xm(i,j-1)-xm(i-1,j-1));
            lh=y(i,j)-ym(i-1,j-1);
            rh=y(i,j)-ym(i,j-1);
            vt(i,j)=(rh+lh)/2*wd;
            wd=(x(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            vt(i,j)=vt(i,j)+lh*wd;
            
        elseif (i > 1 && j< size(x,2)  && bi(i,j) == 0 &&  bi(i-1,j) == 0 && bi(i,j+1) == 0 && bi(i+1,j-1) > 0) % lower Right Edge
            wd=(xm(i,j-1)-xm(i-1,j-1));
            lh=y(i,j)-ym(i-1,j-1);
            rh=y(i,j)-ym(i,j-1);
            vt(i,j)=(rh+lh)/2*wd;
            wd=(xm(i,j)-x(i,j));
            rh=(ym(i,j)-y(i,j));
            vt(i,j)=vt(i,j)+rh*wd;
        
        else
            wd=((xm(i,j)-xm(i-1,j))+(xm(i,j-1)-xm(i-1,j-1)))/2;
            lh=ym(i-1,j)-ym(i-1,j-1);
            rh=ym(i,j)-ym(i,j-1);
            vt(i,j)=(rh+lh)/2*wd;
        end
    end
%     if i < ld
%         cnt=cnt+1;
%     end
end
% surf(x,y,vt)
% view(0,90)