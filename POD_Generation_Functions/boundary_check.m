function [ bnd_x, bnd_y, bnd_idx ] = boundary_check(x, mean_u)
%Defining the boundary and checking
%processing velocity data
%muz{1}=mean(u{1},3);
%mvz{1}=mean(v{1},3);
s=[size(x,1); size(x,2)];
bnd_idx=ones(size(x));
% finding boundaries of image
cnt=1;
for ii = 1:s(1)
    for jj = 1:s(2)
        if mean_u(ii,jj) == 0
            if (ii > 2 && ii < s(1)-1 && jj > 2 && jj < s(2)-1)
                if  (mean_u(ii-2,jj) == 0 && mean_u(ii-1,jj) == 0 && mean_u(ii+1,jj) == 0 && mean_u(ii+2,jj) == 0 && ...
                        mean_u(ii,jj-2) == 0 && mean_u(ii,jj-1) == 0 && mean_u(ii,jj+1) == 0 && mean_u(ii,jj+2) == 0)
                    if (abs(mean_u(ii-1,jj-1)) > 0 || abs(mean_u(ii-1,jj+1)) > 0 ||abs(mean_u(ii+1,jj-1)) > 0 ||abs(mean_u(ii+1,jj+1)) > 0)
                        bnd_idx(ii, jj)=0;
                        cnt=cnt+1;
                    else
                        bnd_idx(ii,jj)=-1;
                    end
                elseif (mean_u(ii-2,jj) == 0 && mean_u(ii-1,jj) == 0 && abs(mean_u(ii+1,jj))> 0 ...
                        || mean_u(ii+2,jj) == 0 && mean_u(ii+1,jj) == 0 && abs(mean_u(ii-1,jj))> 0 ...
                        || mean_u(ii,jj-2) == 0 && mean_u(ii,jj-1) == 0 && abs(mean_u(ii,jj+1))> 0 ...
                        || mean_u(ii,jj+2) == 0 && mean_u(ii,jj+1) == 0 && abs(mean_u(ii,jj-1))> 0)
                    if (mean_u(ii-1,jj-1) == 0 || mean_u(ii-1,jj+1) == 0 || mean_u(ii+1,jj-1) == 0 || mean_u(ii+1,jj+1) == 0)
                        bnd_idx(ii,jj)=0;
                        cnt=cnt+1;
                    else
                        bnd_idx(ii,jj)=-1;
                    end
                elseif  (mean_u(ii-1,jj) == 0 && mean_u(ii+1,jj) == 0 && mean_u(ii,jj-1) == 0 && mean_u(ii,jj+1) == 0)
                    if (abs(mean_u(ii-1,jj-1)) > 0 || abs(mean_u(ii-1,jj+1)) > 0 ||abs(mean_u(ii+1,jj-1)) > 0 ||abs(mean_u(ii+1,jj+1)) > 0)
                        bnd_idx(ii, jj)=0;
                        cnt=cnt+1;
                    else
                        bnd_idx(ii,jj)=-1;
                    end
                end
            elseif (ii > 1 && ii < s(1) && jj > 1 && jj < s(2))
                
                if  (mean_u(ii-1,jj) == 0 && mean_u(ii+1,jj) == 0 && mean_u(ii,jj-1) == 0 && mean_u(ii,jj+1) == 0)
                    if (abs(mean_u(ii-1,jj-1)) > 0 || abs(mean_u(ii-1,jj+1)) > 0 ||abs(mean_u(ii+1,jj-1)) > 0 ||abs(mean_u(ii+1,jj+1)) > 0)
                        bnd_idx(ii, jj)=0;
                        cnt=cnt+1;
                    else
                        bnd_idx(ii,jj)=-1;
                    end
                    
                elseif (mean_u(ii-1,jj) == 0 && abs(mean_u(ii+1,jj))> 0 ...
                        ||  mean_u(ii+1,jj) == 0 && abs(mean_u(ii-1,jj))> 0 ...
                        ||  mean_u(ii,jj-1) == 0 && abs(mean_u(ii,jj+1))> 0 ...
                        ||  mean_u(ii,jj+1) == 0 && abs(mean_u(ii,jj-1))> 0)
                    
                    if (mean_u(ii-1,jj-1) == 0 || mean_u(ii-1,jj+1) == 0 || mean_u(ii+1,jj-1) == 0 || mean_u(ii+1,jj+1) == 0)
                        bnd_idx(ii,jj)=0;
                        cnt=cnt+1;
                    else
                        bnd_idx(ii,jj)=-1;
                    end
                end
            else
                if ii == 1
                    if jj ==1
                        if (mean_u(ii+1,jj) ==0 && abs(mean_u(ii,jj+1)) > 0 || mean_u(ii,jj+1) == 0 && abs(mean_u(ii+1,jj)) > 0 || mean_u(ii,jj+1) == 0 && mean_u(ii+1,jj) == 0 && abs(mean_u(ii+1,jj+1))> 0)
                            bnd_idx(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (mean_u(ii,jj+1) == 0 && mean_u(ii+1,jj) == 0 && mean_u(ii+1,jj+1) == 0)
                            bnd_idx(ii,jj)=-1;
                        end
                    elseif jj == s(2)
                        if (mean_u(ii+1,jj) ==0 && abs(mean_u(ii,jj-1)) > 0 || mean_u(ii,jj-1) == 0 && abs(mean_u(ii+1,jj)) > 0 || mean_u(ii,jj-1) == 0 && mean_u(ii+1,jj) == 0 && abs(mean_u(ii+1,jj-1))> 0)
                            bnd_idx(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (mean_u(ii,jj-1) == 0 && mean_u(ii+1,jj) == 0 && mean_u(ii+1,jj-1) == 0)
                            bnd_idx(ii,jj)=-1;
                        end
                    else
                        if ((mean_u(ii,jj+1) == 0 && mean_u(ii+1,jj) == 0 && abs(mean_u(ii,jj-1)) > 0) ...
                                ||(mean_u(ii,jj-1) == 0 && mean_u(ii+1,jj) == 0 && abs(mean_u(ii,jj+1)) > 0))
                            bnd_idx(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (mean_u(ii,jj-1) == 0 && mean_u(ii,jj+1) == 0 && mean_u(ii+1,jj) == 0)
                            bnd_idx(ii,jj)=-1;
                        end
                    end
                elseif ii == s(1)
                    if jj == 1
                        if (mean_u(ii-1,jj) ==0 && abs(mean_u(ii,jj+1)) > 0 || mean_u(ii,jj+1) == 0 && abs(mean_u(ii-1,jj)) > 0 || mean_u(ii,jj+1) == 0 && mean_u(ii-1,jj) == 0 && abs(mean_u(ii-1,jj+1))> 0)
                            bnd_idx(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (mean_u(ii,jj+1) == 0 && mean_u(ii-1,jj) == 0 && mean_u(ii-1,jj+1) == 0)
                            bnd_idx(ii,jj)=-1;
                        end
                        
                    elseif jj == s(2)
                        if (mean_u(ii-1,jj) ==0 && abs(mean_u(ii,jj-1)) > 0 || mean_u(ii,jj-1) == 0 && abs(mean_u(ii-1,jj)) > 0 || mean_u(ii,jj-1) == 0 && mean_u(ii-1,jj) == 0 && abs(mean_u(ii-1,jj-1))> 0)
                            bnd_idx(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (mean_u(ii,jj-1) == 0 && mean_u(ii-1,jj) == 0 && mean_u(ii-1,jj-1) == 0)
                            bnd_idx(ii,jj)=-1;
                        end
                    else
                        if ((mean_u(ii,jj+1) == 0 && mean_u(ii-1,jj) == 0 && abs(mean_u(ii,jj-1)) > 0) ...
                                || (mean_u(ii,jj-1) == 0 && mean_u(ii-1,jj) == 0 && abs(mean_u(ii,jj+1)) > 0))
                            bnd_idx(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (mean_u(ii,jj-1) == 0 && mean_u(ii,jj+1) == 0 && mean_u(ii-1,jj) == 0)
                            bnd_idx(ii,jj)=-1;
                        end
                    end
                else
                    if jj == 1
                        if ((mean_u(ii+1,jj) == 0 || mean_u(ii-1,jj) == 0) && abs(mean_u(ii,jj+1)) > 0 ...
                            || mean_u(ii-1,jj) == 0 && mean_u(ii+1,jj) == 0 && mean_u(ii,jj+1) == 0 ...
                            && (abs(mean_u(ii-1,jj+1)) > 0 || abs(mean_u(ii+1,jj+1)) > 0))
                            bnd_idx(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (mean_u(ii-1,jj) == 0 && mean_u(ii+1,jj) == 0 && mean_u(ii,jj+1) == 0)
                            bnd_idx(ii,jj)=-1;
                        end
                    elseif jj == s(2)
                        if ((mean_u(ii+1,jj) == 0 || mean_u(ii-1,jj) == 0) && abs(mean_u(ii,jj-1)) > 0 ...
                           || mean_u(ii-1,jj) == 0 && mean_u(ii+1,jj) == 0 && mean_u(ii,jj-1) == 0 ...
                            && (abs(mean_u(ii-1,jj-1)) > 0 || abs(mean_u(ii+1,jj-1) > 0)))
                            bnd_idx(ii,jj)=0;
                            cnt=cnt+1;
                        elseif   (mean_u(ii-1,jj) == 0 && mean_u(ii+1,jj) == 0 && mean_u(ii,jj-1) == 0)
                            bnd_idx(ii,jj)=-1;
                        end
                    end
                end
            end
        end
    end
end

[bnd_x, bnd_y] = graident(bnd_idx);
