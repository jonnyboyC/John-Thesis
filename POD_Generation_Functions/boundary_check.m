function [ X, Y, xb, yb, ib, jb, bi ] = boundary_check(x, y, muz)
%Defining the boundary and checking
%processing velocity data
%muz{1}=mean(u{1},3);
%mvz{1}=mean(v{1},3);
s=[size(x,1); size(x,2)];
s1=s;
bi=ones(size(x));
xh=[];
xv=[];
% finding boundaries of image
cnt=1;
for ii = 1:s(1)
    for jj = 1:s(2)
        if muz(ii,jj) == 0
            if (ii > 2 && ii < s(1)-1 && jj > 2 && jj < s(2)-1)
                if  (muz(ii-2,jj) == 0 && muz(ii-1,jj) == 0 && muz(ii+1,jj) == 0 && muz(ii+2,jj) == 0 && ...
                        muz(ii,jj-2) == 0 && muz(ii,jj-1) == 0 && muz(ii,jj+1) == 0 && muz(ii,jj+2) == 0)
                    if (abs(muz(ii-1,jj-1)) > 0 || abs(muz(ii-1,jj+1)) > 0 ||abs(muz(ii+1,jj-1)) > 0 ||abs(muz(ii+1,jj+1)) > 0)
                        xh(cnt)=x(ii,jj);
                        yh(cnt)=y(ii,jj);
                        indhx(cnt)=ii;
                        indhy(cnt)=jj;
                        bi(ii, jj)=0;
                        cnt=cnt+1;
                    else
                        bi(ii,jj)=-1;
                    end
                    
                elseif (muz(ii-2,jj) == 0 && muz(ii-1,jj) == 0 && abs(muz(ii+1,jj))> 0 ...
                        || muz(ii+2,jj) == 0 && muz(ii+1,jj) == 0 && abs(muz(ii-1,jj))> 0 ...
                        || muz(ii,jj-2) == 0 && muz(ii,jj-1) == 0 && abs(muz(ii,jj+1))> 0 ...
                        || muz(ii,jj+2) == 0 && muz(ii,jj+1) == 0 && abs(muz(ii,jj-1))> 0)
                    if (muz(ii-1,jj-1) == 0 || muz(ii-1,jj+1) == 0 || muz(ii+1,jj-1) == 0 || muz(ii+1,jj+1) == 0)
                        xh(cnt)=x(ii,jj);
                        yh(cnt)=y(ii,jj);
                        indhx(cnt)=ii;
                        indhy(cnt)=jj;
                        bi(ii,jj)=0;
                        cnt=cnt+1;
                    else
                        bi(ii,jj)=-1;
                    end
                    
                elseif  (muz(ii-1,jj) == 0 && muz(ii+1,jj) == 0 && muz(ii,jj-1) == 0 && muz(ii,jj+1) == 0)
                    if (abs(muz(ii-1,jj-1)) > 0 || abs(muz(ii-1,jj+1)) > 0 ||abs(muz(ii+1,jj-1)) > 0 ||abs(muz(ii+1,jj+1)) > 0)
                        xh(cnt)=x(ii,jj);
                        yh(cnt)=y(ii,jj);
                        indhx(cnt)=ii;
                        indhy(cnt)=jj;
                        bi(ii, jj)=0;
                        cnt=cnt+1;
                    else
                        bi(ii,jj)=-1;
                    end
                end
                
            elseif (ii > 1 && ii < s(1) && jj > 1 && jj < s(2))
                
                if  (muz(ii-1,jj) == 0 && muz(ii+1,jj) == 0 && muz(ii,jj-1) == 0 && muz(ii,jj+1) == 0)
                    if (abs(muz(ii-1,jj-1)) > 0 || abs(muz(ii-1,jj+1)) > 0 ||abs(muz(ii+1,jj-1)) > 0 ||abs(muz(ii+1,jj+1)) > 0)
                        xh(cnt)=x(ii,jj);
                        yh(cnt)=y(ii,jj);
                        indhx(cnt)=ii;
                        indhy(cnt)=jj;
                        bi(ii, jj)=0;
                        cnt=cnt+1;
                    else
                        bi(ii,jj)=-1;
                    end
                    
                elseif (muz(ii-1,jj) == 0 && abs(muz(ii+1,jj))> 0 ...
                        ||  muz(ii+1,jj) == 0 && abs(muz(ii-1,jj))> 0 ...
                        ||  muz(ii,jj-1) == 0 && abs(muz(ii,jj+1))> 0 ...
                        ||  muz(ii,jj+1) == 0 && abs(muz(ii,jj-1))> 0)
                    
                    if (muz(ii-1,jj-1) == 0 || muz(ii-1,jj+1) == 0 || muz(ii+1,jj-1) == 0 || muz(ii+1,jj+1) == 0)
                        xh(cnt)=x(ii,jj);
                        yh(cnt)=y(ii,jj);
                        indhx(cnt)=ii;
                        indhy(cnt)=jj;
                        bi(ii,jj)=0;
                        cnt=cnt+1;
                    else
                        bi(ii,jj)=-1;
                    end
                end
            else
                if ii == 1
                    if jj ==1
                        if (muz(ii+1,jj) ==0 && abs(muz(ii,jj+1)) > 0 || muz(ii,jj+1) == 0 && abs(muz(ii+1,jj)) > 0 || muz(ii,jj+1) == 0 && muz(ii+1,jj) == 0 && abs(muz(ii+1,jj+1))> 0)
                            xh(cnt)=x(ii,jj);
                            yh(cnt)=y(ii,jj);
                            indhx(cnt)=ii;
                            indhy(cnt)=jj;
                            bi(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (muz(ii,jj+1) == 0 && muz(ii+1,jj) == 0 && muz(ii+1,jj+1) == 0)
                            bi(ii,jj)=-1;
                        end
                    elseif jj == s(2)
                        if (muz(ii+1,jj) ==0 && abs(muz(ii,jj-1)) > 0 || muz(ii,jj-1) == 0 && abs(muz(ii+1,jj)) > 0 || muz(ii,jj-1) == 0 && muz(ii+1,jj) == 0 && abs(muz(ii+1,jj-1))> 0)
                            xh(cnt)=x(ii,jj);
                            yh(cnt)=y(ii,jj);
                            indhx(cnt)=ii;
                            indhy(cnt)=jj;
                            bi(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (muz(ii,jj-1) == 0 && muz(ii+1,jj) == 0 && muz(ii+1,jj-1) == 0)
                            bi(ii,jj)=-1;
                        end
                    else
                        if ((muz(ii,jj+1) == 0 && muz(ii+1,jj) == 0 && abs(muz(ii,jj-1)) > 0) ...
                                ||(muz(ii,jj-1) == 0 && muz(ii+1,jj) == 0 && abs(muz(ii,jj+1)) > 0))
                            xh(cnt)=x(ii,jj);
                            yh(cnt)=y(ii,jj);
                            indhx(cnt)=ii;
                            indhy(cnt)=jj;
                            bi(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (muz(ii,jj-1) == 0 && muz(ii,jj+1) == 0 && muz(ii+1,jj) == 0)
                            bi(ii,jj)=-1;
                        end
                    end
                    
                elseif ii == s(1)
                    if jj == 1
                        if (muz(ii-1,jj) ==0 && abs(muz(ii,jj+1)) > 0 || muz(ii,jj+1) == 0 && abs(muz(ii-1,jj)) > 0 || muz(ii,jj+1) == 0 && muz(ii-1,jj) == 0 && abs(muz(ii-1,jj+1))> 0)
                            xh(cnt)=x(ii,jj);
                            yh(cnt)=y(ii,jj);
                            indhx(cnt)=ii;
                            indhy(cnt)=jj;
                            bi(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (muz(ii,jj+1) == 0 && muz(ii-1,jj) == 0 && muz(ii-1,jj+1) == 0)
                            bi(ii,jj)=-1;
                        end
                        
                    elseif jj == s(2)
                        if (muz(ii-1,jj) ==0 && abs(muz(ii,jj-1)) > 0 || muz(ii,jj-1) == 0 && abs(muz(ii-1,jj)) > 0 || muz(ii,jj-1) == 0 && muz(ii-1,jj) == 0 && abs(muz(ii-1,jj-1))> 0)
                            xh(cnt)=x(ii,jj);
                            yh(cnt)=y(ii,jj);
                            indhx(cnt)=ii;
                            indhy(cnt)=jj;
                            bi(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (muz(ii,jj-1) == 0 && muz(ii-1,jj) == 0 && muz(ii-1,jj-1) == 0)
                            bi(ii,jj)=-1;
                        end
                    else
                        if ((muz(ii,jj+1) == 0 && muz(ii-1,jj) == 0 && abs(muz(ii,jj-1)) > 0) ...
                                || (muz(ii,jj-1) == 0 && muz(ii-1,jj) == 0 && abs(muz(ii,jj+1)) > 0))
                            xh(cnt)=x(ii,jj);
                            yh(cnt)=y(ii,jj);
                            indhx(cnt)=ii;
                            indhy(cnt)=jj;
                            bi(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (muz(ii,jj-1) == 0 && muz(ii,jj+1) == 0 && muz(ii-1,jj) == 0)
                            bi(ii,jj)=-1;
                        end
                    end
                else
                    if jj == 1
                        if ((muz(ii+1,jj) == 0 || muz(ii-1,jj) == 0) && abs(muz(ii,jj+1)) > 0 ...
                            || muz(ii-1,jj) == 0 && muz(ii+1,jj) == 0 && muz(ii,jj+1) == 0 ...
                            && (abs(muz(ii-1,jj+1)) > 0 || abs(muz(ii+1,jj+1)) > 0))
                            xh(cnt)=x(ii,jj);
                            yh(cnt)=y(ii,jj);
                            indhx(cnt)=ii;
                            indhy(cnt)=jj;
                            bi(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (muz(ii-1,jj) == 0 && muz(ii+1,jj) == 0 && muz(ii,jj+1) == 0)
                            bi(ii,jj)=-1;
                        end
                    elseif jj == s(2)
                        if ((muz(ii+1,jj) == 0 || muz(ii-1,jj) == 0) && abs(muz(ii,jj-1)) > 0 ...
                           || muz(ii-1,jj) == 0 && muz(ii+1,jj) == 0 && muz(ii,jj-1) == 0 ...
                            && (abs(muz(ii-1,jj-1)) > 0 || abs(muz(ii+1,jj-1) > 0)))
                            xh(cnt)=x(ii,jj);
                            yh(cnt)=y(ii,jj);
                            indhx(cnt)=ii;
                            indhy(cnt)=jj;
                            bi(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (muz(ii-1,jj) == 0 && muz(ii+1,jj) == 0 && muz(ii,jj-1) == 0)
                            bi(ii,jj)=-1;
                        end
                    end
                end
            end
        end
    end
end

if size(xh,2)>1
    [Xh,IX] = sort(xh);
    Yh=yh(IX);
    ih=indhx(IX);
    jh=indhy(IX);
end

if size(xh,2)>1 
    Xn = Xh;
    Yn = Yh;
    in = ih;
    jn = jh;
    xn = xh;
    yn = yh;
else
    X = [];
    Y = [];
    ib=[];
    jb=[];
    xb=[];
    yb=[];
end
Rp=[];
ct = 1;
if size(xh,2)>1 
    for i = 1 : size(Xn,2)-1
        for j = i+1:size(Xn,2)
            if (Xn(1,i) == Xn(1,j) && Yn(1,i) == Yn(1,j))
                Rp(ct)= j;
                ct=ct+1;
            end
        end
    end
    ct=1;
    tv=1;
    for i = 1 : size(Xn,2)
        if (ct <= size(Rp,2) && Rp(ct) == i)
            ct=ct+1;
        else
            X(tv)=Xn(1,i);
            Y(tv)=Yn(1,i);
            ib(tv)=in(1,i);
            jb(tv)=jn(1,i);
            xb(tv)=xn(1,i);
            yb(tv)=yn(1,i);
            tv=tv+1;
        end
    end
end

% plot(X,Y)

