function [xxi,yxi,xet,yet,aj] = metric2(x,y)
% This function Appears to simply take an advanced approach to determining
% the spacing between nodes of the pvi image, currently should come out
% exactly the same as taking the difference between neighbors.

dxi = 1;            % May be scaling constants
det = 1;            % "     "       "
dxi12 = 12.*dxi;    
det12 = 12.*det;

xxi = zeros(size(x));
yxi = zeros(size(x));
xet = zeros(size(x));
yet = zeros(size(x));
aj = zeros(size(x));


% Loop through each x and y and generated xxi yxi xet yet for each x/y,
% currently calculates for values even inside the boundary. 
for j = 1:size(x,2)
    for i=1:size(x,1)
        if (i==1)           
            xxi(i,j)=(-3*x(i+4,j)+16*x(i+3,j)-36*x(i+2,j)+48*x(i+1,j)-25*x(i,j))/dxi12;
            yxi(i,j)=(-3*y(i+4,j)+16*y(i+3,j)-36*y(i+2,j)+48*y(i+1,j)-25*y(i,j))/dxi12;
        elseif (i==2)       
            xxi(i,j)=(x(i+3,j)-6*x(i+2,j)+18*x(i+1,j)-10*x(i,j)-3*x(i-1,j))/dxi12;
            yxi(i,j)=(y(i+3,j)-6*y(i+2,j)+18*y(i+1,j)-10*y(i,j)-3*y(i-1,j))/dxi12; 
        elseif (i == size(x,1)-1)
            xxi(i,j)=(-x(i-3,j)+6*x(i-2,j)-18*x(i-1,j)+10*x(i,j)+3*x(i+1,j))/dxi12;
            yxi(i,j)=(-y(i-3,j)+6*y(i-2,j)-18*y(i-1,j)+10*y(i,j)+3*y(i+1,j))/dxi12;
        elseif (i == size(x,1))
            xxi(i,j)=(3*x(i-4,j)-16*x(i-3,j)+36*x(i-2,j)-48*x(i-1,j)+25*x(i,j))/dxi12;
            yxi(i,j)=(3*y(i-4,j)-16*y(i-3,j)+36*y(i-2,j)-48*y(i-1,j)+25*y(i,j))/dxi12;
        else
            xxi(i,j)=(-x(i+2,j)+8*x(i+1,j)-8*x(i-1,j)+x(i-2,j))/dxi12;
            yxi(i,j)=(-y(i+2,j)+8*y(i+1,j)-8*y(i-1,j)+y(i-2,j))/dxi12;
        end
        
        if (j==1)
            xet(i,j)=(-3*x(i,j+4)+16*x(i,j+3)-36*x(i,j+2)+48*x(i,j+1)-25*x(i,j))/det12;
            yet(i,j)=(-3*y(i,j+4)+16*y(i,j+3)-36*y(i,j+2)+48*y(i,j+1)-25*y(i,j))/det12;
        elseif (j==2)
            xet(i,j)=(x(i,j+3)-6*x(i,j+2)+18*x(i,j+1)-10*x(i,j)-3*x(i,j-1))/det12;
            yet(i,j)=(y(i,j+3)-6*y(i,j+2)+18*y(i,j+1)-10*y(i,j)-3*y(i,j-1))/det12; 
        elseif j == size(x,2)-1
            xet(i,j)=(-x(i,j-3)+6*x(i,j-2)-18*x(i,j-1)+10*x(i,j)+3*x(i,j+1))/det12;
            yet(i,j)=(-y(i,j-3)+6*y(i,j-2)-18*y(i,j-1)+10*y(i,j)+3*y(i,j+1))/det12;
        elseif j == size(x,2)
            xet(i,j)=(3*x(i,j-4)-16*x(i,j-3)+36*x(i,j-2)-48*x(i,j-1)+25*x(i,j))/det12;
            yet(i,j)=(3*y(i,j-4)-16*y(i,j-3)+36*y(i,j-2)-48*y(i,j-1)+25*y(i,j))/det12;
        else
            xet(i,j)=(-x(i,j+2)+8*x(i,j+1)-8*x(i,j-1)+x(i,j-2))/det12;
            yet(i,j)=(-y(i,j+2)+8*y(i,j+1)-8*y(i,j-1)+y(i,j-2))/det12;
        end    
        % Would think this would be a very large number since xxi*yet and
        % yxi*xet should be very close to each other
        aj(i,j)=1./(xxi(i,j)*yet(i,j)-yxi(i,j)*xet(i,j));
    end
end
