function [dprdxic,dprdetc]=visder_par(pr,nx,ny,z,bi)
% Preliminary step to calculating finite difference

% Check that parallel pool is ready
if isempty(gcp)
    parpool;
end

pool = gcp;
workers = pool.NumWorkers;

% Find ranges
range = cell(workers,1);
start = zeros(workers,1);
finish = zeros(workers,1);

% Split image into equal parts for each worker
for i = 1:workers;
    start(i) = (i-1)*floor((nx/workers))+1;
    if i ~= workers
        finish(i) = i*floor((nx/workers));
    else
        finish(i) = nx;
    end
    range{i} = start(i):finish(i);
end

% preallocate results
dprdxic=zeros(nx, ny);
dprdetc=zeros(nx, ny);

% Distrubute jobs to workers
for i = 1:workers
    par_handles(i) = parfeval(pool, @par_worker, 2, pr, range{i}, nx, ny, bi);
end

% Place results in global matrix when they finish
for i = 1:workers
    [idx, ic_temp, tc_temp] = fetchNext(par_handles);
    dprdxic(start(idx):finish(idx),:) = ic_temp;
    dprdetc(start(idx):finish(idx),:) = tc_temp;
end

% Multiply by depth field
dprdxic=dprdxic.*z;
dprdetc=dprdetc.*z;
end

function [dprdxic, dprdetc] = par_worker(pr,range,nx,ny,bi)

dxi=1;
det=1;
rdxi=1./dxi;
rdet=1./det;
rdxi12=rdxi/12;
rdet12=rdet/12;
dprdxic=zeros(length(range), ny);
dprdetc=zeros(length(range), ny);

for i = range
    for j = 1:ny
        k = i-range(1)+1;
        % If within a boundary or on it's edge
        if  (bi(i,j) < 1 || i > 1 && i < nx && bi(i,j) == 0 && bi(i-1,j) < 1 && bi(i+1,j) < 1)
            dprdxic(k,j) = 0;
            
        % Either left side of image or just right of a boundary    
        elseif (i == 1 || ( i < nx && bi(i-1,j) == 0 && bi(i,j) > 0)) 
            if (i == nx-1 || bi(i+2,j) == 0 )
                dprdxic(k,j) = (pr(i+1,j)-pr(i,j))*rdxi12*12;
            elseif (i == nx-2 || bi(i+3,j) == 0 )
                dprdxic(k,j) = (-pr(i+2,j)+4*pr(i+1,j)-3*pr(i,j))*rdxi12*6;
            elseif (i == nx-3 || bi(i+4,j) == 0 )
                dprdxic(k,j) = (2*pr(i+3,j)-9*pr(i+2,j)+18*pr(i+1,j)-11*pr(i,j))*rdxi12*2;
            else
                dprdxic(k,j) = (-3*pr(i+4,j)+16*pr(i+3,j)-36*pr(i+2,j)+48*pr(i+1,j)-25*pr(i,j))*rdxi12;
            end
            
            
        % Either one right of left side of image, or just one right of a
        % boundary
        elseif (i == 2 || i < nx && (bi(i-2,j) == 0 && bi(i-1,j) > 0)) 
            if (i == nx-1 || bi(i+2,j) == 0 )
                dprdxic(k,j)=(pr(i+1,j)-pr(i-1,j))*rdxi12*6;
            elseif (i == nx-2 || bi(i+3,j) == 0 )
                dprdxic(k,j)=(-pr(i+2,j)+6*pr(i+1,j)-3*pr(i,j)-2*pr(i-1,j))*rdxi12*2;
            else
                dprdxic(k,j)=(pr(i+3,j)-6*pr(i+2,j)+18*pr(i+1,j)-10*pr(i,j)-3*pr(i-1,j))*rdxi12;
            end
            
        % Either one left of right side of image, or just one left of a
        % boundary
        elseif (i == nx-1 || (i > 1 && i < nx && bi(i+2,j)== 0 && bi(i+1,j)> 0)) 
            if (i == 2 || bi(i-2,j) == 0)
                dprdxic(k,j)=(pr(i+1,j)-pr(i-1,j))*rdxi12*6;
            elseif (i == 3 || bi(i-3,j) == 0)
                dprdxic(k,j)=(pr(i-2,j)-6*pr(i-1,j)+3*pr(i,j)+2*pr(i+1,j))*rdxi12*2;
            else
                dprdxic(k,j)=(-pr(i-3,j)+6*pr(i-2,j)-18*pr(i-1,j)+10*pr(i,j)+3*pr(i+1,j))*rdxi12;
            end
            
        % Either right side of image, or just left of a boundary
        elseif (i == nx || (i > 1 && bi(i+1,j) == 0 && bi(i,j)> 0)) 
            if (i == 2  ||bi(i-2,j) == 0  )
                dprdxic(k,j)=(pr(i,j)-pr(i-1,j))*rdxi12*12;
            elseif ( i == 3 || bi(i-3,j) == 0)
                dprdxic(k,j)=(pr(i-2,j)-4*pr(i-1,j)+3*pr(i,j))*rdxi12*6;
            elseif (i == 4 || bi(i-4,j) == 0)
                dprdxic(k,j)=(-2*pr(i-3,j)+9*pr(i-2,j)-18*pr(i-1,j)+11*pr(i,j))*rdxi12*2;
            else
                dprdxic(k,j)=(3*pr(i-4,j)-16*pr(i-3,j)+36*pr(i-2,j)-48*pr(i-1,j)+25*pr(i,j))*rdxi12;
            end
            
        else
            
            % If not near a boundary or the edge of the image
            dprdxic(k,j)=(-pr(i+2,j)+8*pr(i+1,j)-8*pr(i-1,j)+pr(i-2,j))*rdxi12; 
        end
        
        % If within a boundary or on it's edge
        if  (bi(i,j) < 1 || (j > 1 && j < ny && bi(i,j) == 0 && bi(i,j-1) < 1 && bi(i,j+1) < 1))
            dprdetc(k,j) = 0;
            
        % Either bottom of image, or just above a boundary
        elseif (j == 1 || ( j < ny && bi(i,j-1) == 0 && bi(i,j)> 0))  
            if (j == ny-1 || bi(i,j+2) == 0)
                dprdetc(k,j)=(pr(i,j+1)-pr(i,j))*rdet12*12;
            elseif (j == ny-2 || bi(i,j+3) == 0)
                dprdetc(k,j)=(-pr(i,j+2)+4*pr(i,j+1)-3*pr(i,j))*rdet12*6;
            elseif (j == ny-3 || bi(i,j+4) == 0)
                dprdetc(k,j)=(2*pr(i,j+3)-9*pr(i,j+2)+18*pr(i,j+1)-11*pr(i,j))*rdet12*2;
            else
                dprdetc(k,j)=(-3*pr(i,j+4)+16*pr(i,j+3)-36*pr(i,j+2)+48*pr(i,j+1)-25*pr(i,j))*rdet12;
            end
            
        % Either one above the bottom of image, or just one above a
        % boundary
        elseif (j == 2 || (j < ny && bi(i,j-2) == 0 && bi(i,j-1)> 0)) 
            if (bi(i,j) == 0 && bi(i,j+1) <1)
                dprdetc(k,j)=(pr(i,j)-pr(i,j-1))*rdet12*12;
            elseif (j == ny-1 || bi(i,j+2) == 0)
                dprdetc(k,j)=(pr(i,j+1)-pr(i,j-1))*rdet12*6;
            elseif (j == ny-2 || bi(i,j+3) == 0)
                dprdetc(k,j)=(-pr(i,j+2)+6*pr(i,j+1)-3*pr(i,j)-2*pr(i,j-1))*rdet12*2;
            else
                dprdetc(k,j)=(pr(i,j+3)-6*pr(i,j+2)+18*pr(i,j+1)-10*pr(i,j)-3*pr(i,j-1))*rdet12;
            end
            
            
        % Either one below the top of image, or just one below a boundary    
        elseif (j == ny-1 || (j > 1 && j < ny && bi(i,j+2) == 0 && bi(i,j+1)> 0)) 
            if (j == 2 || bi(i,j-2) == 0)
                dprdetc(k,j)=(pr(i,j+1)-pr(i,j-1))*rdet12*6;
            elseif (j == 3 || bi(i,j-3) == 0)
                dprdetc(k,j)=(pr(i,j-2)-6*pr(i,j-1)+3*pr(i,j)+2*pr(i,j+1))*rdet12*2;
            else
                dprdetc(k,j)=(-pr(i,j-3)+6*pr(i,j-2)-18*pr(i,j-1)+10*pr(i,j)+3*pr(i,j+1))*rdet12;
            end
            
        % Either top of image, or just below a boundary    
        elseif (j == ny || (j > 1 && bi(i,j+1) == 0 && bi(i,j)> 0)) 
            if (bi(i,j) < 1 && bi(i,j-1) <1)
                dprdetc(k,j)=0;
            elseif (j == 2  || bi(i,j-2) == 0)
                dprdetc(k,j)=(pr(i,j)-pr(i,j-1))*rdet12*12;
            elseif ( j == 3 || bi(i,j-3) == 0)
                dprdetc(k,j)=(pr(i,j-2)-4*pr(i,j-1)+3*pr(i,j))*rdet12*6;
            elseif (j == 4 || bi(i,j-4) == 0)
                dprdetc(k,j)=(-2*pr(i,j-3)+9*pr(i,j-2)-18*pr(i,j-1)+11*pr(i,j))*rdet12*2;
            else
                dprdetc(k,j)=(3*pr(i,j-4)-16*pr(i,j-3)+36*pr(i,j-2)-48*pr(i,j-1)+25*pr(i,j))*rdet12;
            end
            
        % If not near a boundary or the edge of the image
        else
            dprdetc(k,j)=(-pr(i,j+2)+8*pr(i,j+1)-8*pr(i,j-1)+pr(i,j-2))*rdet12; 
        end
    end
end
end
