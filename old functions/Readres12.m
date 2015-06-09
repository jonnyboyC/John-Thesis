function [r,u,v,w,t,ar] = readres(fname,nzone,s,np)
%echo off
%fname = 'css2.dat';
%nzone = 1;
%s=[65 37];
cn=2;
cni=cn;
rb=0;
while cn > 0
    format long
    fid = fopen(fname);
    %fid1 = fopen('swc.dat','a');
    i = 'integer*4';
    f = 'single';
    
    R = cell(nzone,1);
    U = cell(nzone,1);
    V = cell(nzone,1);
    W = cell(nzone,1);
    t = [];
    ar=[];
    kh=0;
    k1 = 0;
    k2 = 0;
    k3=0;
    k4=12;
    k5=24;
    k6=36;
    k7=49;
    k8=62;
    k9=74;
    k10=87;
    %stop = 0;
    %while ~stop
    for q= 1:1000
        time=fread(fid,1,f);
        %if abs(round((q-1)/8)-((q-1)/8))> 0
        fn=rand*10;
        if fn > 1.32
            %if q == 0
            k1 = k1+1;
            for j = 1:nzone
                kmax=fread(fid,1,i);
                jmax=fread(fid,1,i);
                offset = 3;
                len = kmax*jmax;
                r1 = reshape(fread(fid,len,f),kmax,jmax);
                u1 = reshape(fread(fid,len,f),kmax,jmax);
                v1 = reshape(fread(fid,len,f),kmax,jmax);
                w1 = reshape(fread(fid,len,f),kmax,jmax);
                clear r1 u1 v1 w1
            end
        else
            %k2 = k2+1;
            er=(rand*10);
            
            if er< 1.25
                if k3 > 11
                    if k4 > 23
                        if k5 > 35
                            if k6 > 48
                                if k7 > 61
                                    if k8 > 73
                                        if k9 > 86
                                            if k10 < 100
                                                k10=k10+1;
                                                k2=k10;
                                            end
                                        else
                                            k9=k9+1;
                                            k2=k9;
                                        end
                                    else
                                        k8=k8+1;
                                        k2=k8;
                                    end
                                else
                                    k7=k7+1;
                                    k2=k7;
                                end
                            else
                                k6=k6+1;
                                k2=k6;
                            end
                        else
                            k5=k5+1;
                            k2=k5;
                        end
                    else
                        k4=k4+1;
                        k2=k4;
                    end
                else
                    k3=k3+1;
                    k2=k3;
                end
                
            elseif er < 2.5
                if k4 > 23
                    if k5 > 35
                        if k6 > 48
                            if k7 > 61
                                if k8 > 73
                                    if k9 > 86
                                        if k10 >99
                                            if k3<12
                                                k3=k3+1;
                                                k2=k3;
                                            end
                                        else
                                            k10=k10+1;
                                            k2=k10;
                                        end
                                    else
                                        k9=k9+1;
                                        k2=k9;
                                    end
                                else
                                    k8=k8+1;
                                    k2=k8;
                                end
                            else
                                k7=k7+1;
                                k2=k7;
                            end
                        else
                            k6=k6+1;
                            k2=k6;
                        end
                    else
                        k5=k5+1;
                        k2=k5;
                    end
                else
                    k4=k4+1;
                    k2=k4;
                end
                
            elseif er < 3.75
                if k5 > 35
                    if k6 > 48
                        if k7 > 61
                            if k8 > 73
                                if k9 > 86
                                    if k10 >99
                                        if k3 > 11
                                            if k4 < 24
                                                k4=k4+1;
                                                k2=k4;
                                            end
                                        else
                                            k3=k3+1;
                                            k2=k3;
                                        end
                                    else
                                        k10=k10+1;
                                        k2=k10;
                                    end
                                else
                                    k9=k9+1;
                                    k2=k9;
                                end
                            else
                                k8=k8+1;
                                k2=k8;
                            end
                        else
                            k7=k7+1;
                            k2=k7;
                        end
                    else
                        k6=k6+1;
                        k2=k6;
                    end
                else
                    k5=k5+1;
                    k2=k5;
                end
                
            elseif er < 5
                if k6 > 48
                    if k7 > 61
                        if k8 > 73
                            if k9 > 86
                                if k10 > 99
                                    if k3 > 11
                                        if k4 > 23
                                            if k5<36
                                                k5=k5+1;
                                                k2=k5;
                                            end
                                        else
                                            k4=k4+1;
                                            k2=k4;
                                        end
                                    else
                                        k3=k3+1;
                                        k2=k3;
                                    end
                                else
                                    k10=k10+1;
                                    k2=k10;
                                end
                            else
                                k9=k9+1;
                                k2=k9;
                            end
                        else
                            k8=k8+1;
                            k2=k8;
                        end
                    else
                        k7=k7+1;
                        k2=k7;
                    end
                else
                    k6=k6+1;
                    k2=k6;
                end 
                
            elseif er < 6.25
                if k7 > 61
                    if k8 > 73
                        if k9 > 86
                            if k10 > 99
                                if k3 > 11
                                    if k4 >23
                                        if k5 > 35
                                            if k6<49
                                                k6=k6+1;
                                                k2=k6;
                                            end
                                        else
                                            k5=k5+1;
                                            k2=k5;
                                        end
                                    else
                                        k4=k4+1;
                                        k2=k4;
                                    end
                                else
                                    k3=k3+1;
                                    k2=k3;
                                end
                            else
                                k10=k10+1;
                                k2=k10;
                            end
                        else
                            k9=k9+1;
                            k2=k9;
                        end
                    else
                        k8=k8+1;
                        k2=k8;
                    end
                else
                    k7=k7+1;
                    k2=k7;
                end
                
            elseif er < 7.5
                if k8 > 73
                    if k9 > 86
                        if k10 > 99
                            if k3 > 11
                                if k4 > 23
                                    if k5 >35
                                        if k6 > 48
                                            if k7<62
                                                k7=k7+1;
                                                k2=k7;
                                            end
                                        else
                                            k6=k6+1;
                                            k2=k6;
                                        end
                                    else
                                        k5=k5+1;
                                        k2=k5;
                                    end
                                else
                                    k4=k4+1;
                                    k2=k4;
                                end
                            else
                                k3=k3+1;
                                k2=k3;
                            end
                        else
                            k10=k10+1;
                            k2=k10;
                        end
                    else
                        k9=k9+1;
                        k2=k9;
                    end
                else
                    k8=k8+1;
                    k2=k8;
                end
                
            elseif er < 8.75
                if k9 > 86
                    if k10 > 99
                        if k3 > 11
                            if k4 > 23
                                if k5 > 35
                                    if k6 >48
                                        if k7 > 61
                                            if k8<74
                                                k8=k8+1;
                                                k2=k8;
                                            end
                                        else
                                            k7=k7+1;
                                            k2=k7;
                                        end
                                    else
                                        k6=k6+1;
                                        k2=k6;
                                    end
                                else
                                    k5=k5+1;
                                    k2=k5;
                                end
                            else
                                k4=k4+1;
                                k2=k4;
                            end
                        else
                            k3=k3+1;
                            k2=k3;
                        end
                    else
                        k10=k10+1;
                        k2=k10;
                    end
                else
                    k9=k9+1;
                    k2=k9;
                end
            else
                if k10 > 99
                    if k3 > 11
                        if k4 > 23
                            if k5 > 35
                                if k6 > 48
                                    if k7 > 61
                                        if k8 > 73
                                            if k9<87
                                                k9=k9+1;
                                                k2=k9;
                                            end
                                        else
                                            k8=k8+1;
                                            k2=k8;
                                        end
                                    else
                                        k7=k7+1;
                                        k2=k7;
                                    end
                                else
                                    k6=k6+1;
                                    k2=k6;
                                end
                            else
                                k5=k5+1;
                                k2=k5;
                            end
                        else
                            k4=k4+1;
                            k2=k4;
                        end
                    else
                        k3=k3+1;
                        k2=k3;
                    end
                else
                    k10=k10+1;
                    k2=k10;
                end
            end
            %kr = [kr; k2];
            kh=kh+1;
            for j = 1:nzone
                kmax=fread(fid,1,i);
                jmax=fread(fid,1,i);
                offset = 3;
                len = kmax*jmax;
                r2 = reshape(fread(fid,len,f),kmax,jmax);
                u2 = reshape(fread(fid,len,f),kmax,jmax)/3.280839;
                v2 = reshape(fread(fid,len,f),kmax,jmax)/3.280839;
                w2 = reshape(fread(fid,len,f),kmax,jmax)/3.280839;
                if kh > 100
                    k2 = abs(round(rand*100));
                    if k2 < 1 
                        k2=1
                    end
                end
                if cn<cni
                    kk=k2+rb*100;
                else
                    kk=k2;
                end
                mk(kk)=kk;
                r{j}(:,:,kk) = r2(np,:);
                u{j}(:,:,kk) = u2(np,:);
                v{j}(:,:,kk) = v2(np,:);
                w{j}(:,:,kk) = w2(np,:);
                clear r2 v2 w2 u2
                
            end
            top=[q kk];
            ar(k2,1)=q;
            %ar(kh,2)=k2;
            %t = [t time];
            t1(kk) = time;
        end
    end
    if size(ar,1) == 100
        if cn >1
           cn=cn-1;
           rb=rb+1;
           status=fclose(fid);
        else
            cn=0;
        end
    else
        status=fclose(fid);
    end
end
t=t1;
%nm=length(kr)