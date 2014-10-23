function [dprdxic,dprdetc]=visder(pr,dimensions,z,bi)
dxi=1;%./real(l(1)-1);
det=1;%./real(l(2)-1);
rdxi=1./dxi;
rdet=1./det;
rdxi12=rdxi/12;
rdet12=rdet/12;
nx=dimensions(1);
ny=dimensions(2);
dprdxic=zeros(dimensions);
dprdetc=zeros(dimensions);

for i = 1:nx
    for j = 1:ny
        
        if  (bi(i,j) < 1 || i > 1 && i < nx && bi(i,j) == 0 && bi(i-1,j) < 1 && bi(i+1,j) < 1) %== -1
            dprdxic(i,j) = 0;
            
        elseif (i == 1 || ( i < nx && bi(i-1,j) == 0 && bi(i,j) > 0)) % left boundary 
            if (bi(i,j) < 1 && bi(i+1,j) <1)
                dprdxic(i,j) = 0;
            elseif (i == nx-1 || bi(i+2,j) == 0 )
                dprdxic(i,j) = (pr(i+1,j)-pr(i,j))*rdxi12*12;
            elseif (i == nx-2 || bi(i+3,j) == 0 )
                dprdxic(i,j) = (-pr(i+2,j)+4*pr(i+1,j)-3*pr(i,j))*rdxi12*6;
            elseif (i == nx-3 || bi(i+4,j) == 0 )
                dprdxic(i,j) = (2*pr(i+3,j)-9*pr(i+2,j)+18*pr(i+1,j)-11*pr(i,j))*rdxi12*2;
            else
                dprdxic(i,j) = (-3*pr(i+4,j)+16*pr(i+3,j)-36*pr(i+2,j)+48*pr(i+1,j)-25*pr(i,j))*rdxi12;
            end
            
        elseif (i == 2 || i < nx && (bi(i-2,j) == 0 && bi(i-1,j) > 0)) % left boundary plus one
            if (bi(i,j) == 0 && bi(i+1,j) <1)
                dprdxic(i,j)=(pr(i,j)-pr(i-1,j))*rdxi12*12;
            elseif (i == nx-1 || bi(i+2,j) == 0 )
                dprdxic(i,j)=(pr(i+1,j)-pr(i-1,j))*rdxi12*6;
            elseif (i == nx-2 || bi(i+3,j) == 0 )
                dprdxic(i,j)=(-pr(i+2,j)+6*pr(i+1,j)-3*pr(i,j)-2*pr(i-1,j))*rdxi12*2;
            else
                dprdxic(i,j)=(pr(i+3,j)-6*pr(i+2,j)+18*pr(i+1,j)-10*pr(i,j)-3*pr(i-1,j))*rdxi12;
            end
            
        elseif (i == nx-1 || (i < nx && bi(i+2,j)== 0 && bi(i+1,j)> 0)) % Right boundary minus one
            if (bi(i,j) == 0 && bi(i-1,j) <1)
                dprdxic(i,j)=(pr(i+1,j)-pr(i,j))*rdxi12*12;
             elseif (i == 2 || bi(i-2,j) == 0)
                dprdxic(i,j)=(pr(i+1,j)-pr(i-1,j))*rdxi12*6;
            elseif (i == 3 || bi(i-3,j) == 0)
                dprdxic(i,j)=(pr(i-2,j)-6*pr(i-1,j)+3*pr(i,j)+2*pr(i+1,j))*rdxi12*2;
            else
                dprdxic(i,j)=(-pr(i-3,j)+6*pr(i-2,j)-18*pr(i-1,j)+10*pr(i,j)+3*pr(i+1,j))*rdxi12;
            end
            
        elseif (i == nx || (bi(i+1,j) == 0 && bi(i,j)> 0)) % Right boundary
            if (bi(i,j) < 1 && bi(i-1,j) <1)
                dprdxic(i,j)=0;
            elseif (i == 2  ||bi(i-2,j) == 0  )
                dprdxic(i,j)=(pr(i,j)-pr(i-1,j))*rdxi12*12;
            elseif ( i == 3 || bi(i-3,j) == 0)
                dprdxic(i,j)=(pr(i-2,j)-4*pr(i-1,j)+3*pr(i,j))*rdxi12*6;
            elseif (i == 4 || bi(i-4,j) == 0)
                dprdxic(i,j)=(-2*pr(i-3,j)+9*pr(i-2,j)-18*pr(i-1,j)+11*pr(i,j))*rdxi12*2;
            else
                dprdxic(i,j)=(3*pr(i-4,j)-16*pr(i-3,j)+36*pr(i-2,j)-48*pr(i-1,j)+25*pr(i,j))*rdxi12;
            end
            
        else
            dprdxic(i,j)=(-pr(i+2,j)+8*pr(i+1,j)-8*pr(i-1,j)+pr(i-2,j))*rdxi12;
            
        end
        
        if  (bi(i,j) < 1 || j > 1 && j < ny && bi(i,j) == 0 && bi(i,j-1) < 1 && bi(i,j+1) < 1)%== -1 
            dprdetc(i,j) = 0;
            
        elseif (j == 1 || ( j < ny && bi(i,j-1) == 0 && bi(i,j)> 0)) % Bottom boundary 
            if (bi(i,j) < 1 && bi(i,j+1) <1)
                dprdetc(i,j)=0;
            elseif (j == ny-1 || bi(i,j+2) == 0)
                dprdetc(i,j)=(pr(i,j+1)-pr(i,j))*rdet12*12;
            elseif (j == ny-2 || bi(i,j+3) == 0)
                dprdetc(i,j)=(-pr(i,j+2)+4*pr(i,j+1)-3*pr(i,j))*rdet12*6;
            elseif (j == ny-3 || bi(i,j+4) == 0)
                dprdetc(i,j)=(2*pr(i,j+3)-9*pr(i,j+2)+18*pr(i,j+1)-11*pr(i,j))*rdet12*2;
            else
                dprdetc(i,j)=(-3*pr(i,j+4)+16*pr(i,j+3)-36*pr(i,j+2)+48*pr(i,j+1)-25*pr(i,j))*rdet12;
            end
            
        elseif (j == 2 || (bi(i,j-2) == 0 && bi(i,j-1)> 0)) % Bottom boundary plus one
            if (bi(i,j) == 0 && bi(i,j+1) <1)
                dprdetc(i,j)=(pr(i,j)-pr(i,j-1))*rdet12*12;
            elseif (j== ny-1 || bi(i,j+2) == 0)
                dprdetc(i,j)=(pr(i,j+1)-pr(i,j-1))*rdet12*6;
            elseif (j == ny-2 || bi(i,j+3) == 0)
                dprdetc(i,j)=(-pr(i,j+2)+6*pr(i,j+1)-3*pr(i,j)-2*pr(i,j-1))*rdet12*2;
            else
                dprdetc(i,j)=(pr(i,j+3)-6*pr(i,j+2)+18*pr(i,j+1)-10*pr(i,j)-3*pr(i,j-1))*rdet12;
            end
            
        elseif (j == ny-1 || (j < ny && bi(i,j+2) == 0 && bi(i,j+1)> 0)) % Top boundary minus one
            if (bi(i,j) == 0 && bi(i,j-1) <1)
                dprdetc(i,j)=(pr(i,j+1)-pr(i,j))*rdet12*12;
             elseif (j==2 || bi(i,j-2) == 0)
                dprdetc(i,j)=(pr(i,j+1)-pr(i,j-1))*rdet12*6;
            elseif (j == 3 || bi(i,j-3) == 0)
                dprdetc(i,j)=(pr(i,j-2)-6*pr(i,j-1)+3*pr(i,j)+2*pr(i,j+1))*rdet12*2;
            else
                dprdetc(i,j)=(-pr(i,j-3)+6*pr(i,j-2)-18*pr(i,j-1)+10*pr(i,j)+3*pr(i,j+1))*rdet12;
            end
            
        elseif (j == ny || (bi(i,j+1) == 0 && bi(i,j)> 0)) % Top Boundary
            if (bi(i,j) < 1 && bi(i,j-1) <1)
                dprdetc(i,j)=0;
            elseif (j == 2  || bi(i,j-2) == 0)
                dprdetc(i,j)=(pr(i,j)-pr(i,j-1))*rdet12*12;
            elseif ( j == 3 || bi(i,j-3) == 0)
                dprdetc(i,j)=(pr(i,j-2)-4*pr(i,j-1)+3*pr(i,j))*rdet12*6;
            elseif (j == 4 || bi(i,j-4) == 0)
                dprdetc(i,j)=(-2*pr(i,j-3)+9*pr(i,j-2)-18*pr(i,j-1)+11*pr(i,j))*rdet12*2;
            else
                dprdetc(i,j)=(3*pr(i,j-4)-16*pr(i,j-3)+36*pr(i,j-2)-48*pr(i,j-1)+25*pr(i,j))*rdet12;
            end
            
        else
            dprdetc(i,j)=(-pr(i,j+2)+8*pr(i,j+1)-8*pr(i,j-1)+pr(i,j-2))*rdet12;
            
        end
        
    end
end

dprdxic=dprdxic.*z;
dprdetc=dprdetc.*z;
