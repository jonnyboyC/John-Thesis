function vt=voln_piv2(x,y,bnd_idx)
% find vol elements for portions of the velocity image that aren't
% boundaries 
sz=size(x);
vt=ones(size(sz));

% range to fill xm,ym
x_range = 1:sz(1)-1;
y_range = 1:sz(2)-1;

% Determine the interstitial averaged values between the centers of each
% point
xm(x_range, y_range)= (x(x_range,y_range)+x(x_range,y_range+1)+x(x_range+1,y_range)+x(x_range+1,y_range+1))/4;
ym(x_range, y_range)= (y(x_range,y_range)+y(x_range,y_range+1)+y(x_range+1,y_range)+y(x_range+1,y_range+1))/4;

% Determine vol elements for the whole grid
% TODO may want to look at this later, for additional speed up
for i = 1:size(x,1)
    for j=1:size(x,2)
        if  bnd_idx(i,j)== -1
            vt(i,j)=0;
                 
        elseif (i == 1 && (j == size(x,2)|| bnd_idx(i,j) == 0 && bnd_idx(i+1,j) == 0 && bnd_idx(i+1,j-1) > 0))  % Top Left Corner
            wd=(xm(i,j-1)-x(i,j));
            lh=y(i,j)-ym(i,j-1);
            vt(i,j)=lh*wd;
        
        elseif (i == size(x,1) && (j == size(x,2)|| bnd_idx(i,j) == 0 && bnd_idx(i-1,j) == 0 && bnd_idx(i-1,j-1) > 0)) % Top Right Corner
            wd=(x(i,j)-xm(i-1,j-1));
            lh=y(i,j)-ym(i-1,j-1);
            vt(i,j)=lh*wd;
        
        elseif ( i == 1 && (j == 1|| bnd_idx(i,j) == 0 && bnd_idx(i+1,j) == 0 && bnd_idx(i+1,j+1) > 0)) % Bottom Left Corner
            wd=(xm(i,j)-x(i,j));
            lh=ym(i,j)-y(i,j);
            vt(i,j)=lh*wd;
        
        elseif (i == size(x,1) && (j == 1 || bnd_idx(i,j) == 0 && bnd_idx(i-1,j) == 0 && bnd_idx(i-1,j+1) > 0)) % Bottom Right Corner
            wd=(x(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            vt(i,j)=lh*wd;
        
        elseif (i < size(x,1) && j > 1 && bnd_idx(i,j) == 0 &&  bnd_idx(i+1,j) == 0 && bnd_idx(i,j-1) == 0 && bnd_idx(i+1,j-1) > 0) % Upper Left Corner
            wd=(xm(i,j-1)-x(i,j));
            lh=y(i,j)-ym(i,j-1);
            vt(i,j)=lh*wd;
        
        elseif (i > 1 && j > 1 && bnd_idx(i,j) == 0 &&  bnd_idx(i-1,j) == 0 && bnd_idx(i,j-1) == 0 && bnd_idx(i-1,j-1) > 0) % Upper Right Corner
            wd=(x(i,j)-xm(i-1,j-1));
            lh=y(i,j)-ym(i-1,j-1);
            vt(i,j)=lh*wd;
        
        elseif (i < size(x,1) && j < size(x,2) && bnd_idx(i,j) == 0 &&  bnd_idx(i+1,j) == 0 && bnd_idx(i,j+1) == 0 && bnd_idx(i+1,j+1) > 0) % Lower Left Corner
            wd=(xm(i,j)-x(i,j));
            lh=ym(i,j)-y(i,j);
            vt(i,j)=lh*wd;
        
        elseif (i > 1 && j< size(x,2)  && bnd_idx(i,j)== 0 &&  bnd_idx(i-1,j)== 0 && bnd_idx(i,j+1)== 0 && bnd_idx(i-1,j+1)> 0) % lower Right Corner
            wd=(x(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            vt(i,j)=lh*wd;
            
        elseif (j == 1 || i > 1 && i < size(x,1) && j > 1 && bnd_idx(i,j) == 0 &&  bnd_idx(i-1,j) == 0 && bnd_idx(i+1,j) == 0 &&  bnd_idx(i,j+1) > 0) % Bottom Boundary
            wd=(xm(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            rh=ym(i,j)-y(i,j);
            vt(i,j)=(rh+lh)/2*wd;
        
        elseif (i == 1 || i > 1 && j > 1 && j < size(x,2)  &&  bnd_idx(i,j) == 0 &&  bnd_idx(i,j-1) == 0 && bnd_idx(i,j+1) == 0 && bnd_idx(i+1,j) > 0) % Left Boundary
            wd=(xm(i,j)-x(i,j));
            lh=ym(i,j)-ym(i,j-1);
            vt(i,j)=lh*wd;
        
        elseif (i == size(x,1) || i < size(x,1) && j > 1 && j < size(x,2) && bnd_idx(i,j) == 0 &&  bnd_idx(i,j+1) == 0 && bnd_idx(i,j-1) == 0 && bnd_idx(i-1,j)> 0) % Right Boundary
            wd=(x(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-ym(i-1,j-1);
            vt(i,j)=lh*wd;
        
        elseif (j==size(x,2) || i > 1 && i < size(x,1) && j < size(x,2) && bnd_idx(i,j) == 0 &&  bnd_idx(i-1,j) == 0 && bnd_idx(i+1,j) == 0 && bnd_idx(i,j-1)> 0) % Top Boundary
            wd=(xm(i,j-1)-xm(i-1,j-1));
            lh=y(i,j)-ym(i-1,j-1);
            rh=y(i,j)-ym(i,j-1);
            vt(i,j)=(rh+lh)/2*wd;
            
        elseif (i< size(x,1) && j > 1 && bnd_idx(i,j) == 0 &&  bnd_idx(i+1,j) == 0 && bnd_idx(i,j-1) == 0 && bnd_idx(i-1,j+1) > 0) % Upper Left Edge
            wd=(xm(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            rh=ym(i,j)-y(i,j);
            vt(i,j)=(rh+lh)/2*wd;
            wd=(x(i,j)-xm(i-1,j-1));
            lh=(y(i,j)-ym(i-1,j-1));
            vt(i,j)=vt(i,j)+lh*wd;
            
        elseif (i > 1 && j > 1 && bnd_idx(i,j) == 0 &&  bnd_idx(i-1,j) == 0 && bnd_idx(i,j-1) == 0 && bnd_idx(i+1,j+1) > 0) % Upper Right Edge
            wd=(xm(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            rh=ym(i,j)-y(i,j);
            vt(i,j)=(rh+lh)/2*wd;
            wd=(xm(i,j-1)-x(i,j));
            rh=(y(i,j)-ym(i,j-1));
            vt(i,j)=vt(i,j)+rh*wd;
            
        elseif (i< size(x,1) && j < size(x,2) && bnd_idx(i,j) == 0 &&  bnd_idx(i+1,j) == 0 && bnd_idx(i,j+1) == 0 && bnd_idx(i-1,j-1) > 0) % Lower Left Edge
            wd=(xm(i,j-1)-xm(i-1,j-1));
            lh=y(i,j)-ym(i-1,j-1);
            rh=y(i,j)-ym(i,j-1);
            vt(i,j)=(rh+lh)/2*wd;
            wd=(x(i,j)-xm(i-1,j));
            lh=ym(i-1,j)-y(i,j);
            vt(i,j)=vt(i,j)+lh*wd;
            
        elseif (i > 1 && j< size(x,2)  && bnd_idx(i,j) == 0 &&  bnd_idx(i-1,j) == 0 && bnd_idx(i,j+1) == 0 && bnd_idx(i+1,j-1) > 0) % lower Right Edge
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
end
