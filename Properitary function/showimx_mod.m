function [lhs1,lhs2,lhs3,lhs4,lhs5,lhs6,lhs7]=showimx(A, frame)
% CALL:      [lhs1,lhs2,lhs3,lhs4,lhs5,lhs6,lhs7]=showimx(A, frame);
%
% FUNCTION:  Displaying data of LaVision's IMX structure
%            (one vector field, all image frames or only single image frame)
%
% ARGUMENTS: A  = IMX-structure created by READIMX/READIM7 function
%
% RETURN:    in case of images (image type=0):
%               lhs1 = scaled x-coordinates
%               lhs2 = scaled y-coordinates
%               lhs3 = scaled image intensities
%            in case of 2D vector fields (A.IType = 1,2 or 3):
%               lhs1 = scaled x-coordinates
%               lhs2 = scaled y-coordinates
%               lhs3 = scaled vx-components of vectors
%               lhs4 = scaled vy-components of vectors
%               lhs5 = vector choice field 
%            in case of 3D vector fields (A.IType = 4 or 5):
%               lhs1 = scaled x-coordinates
%               lhs2 = scaled y-coordinates
%               lhs3 = scaled z-coordinates
%               lhs4 = scaled vx-components of vectors
%               lhs5 = scaled vy-components of vectors
%               lhs6 = scaled vz-components of vectors
%               lhs7 = vector choice field
if nargin==0,
	help showimx, return
elseif nargin<2
    frame = -1;
elseif ~isscalar(frame)
    frame = -1;
else
    frame = floor(frame);
end
if ~isfield(A,'DaVis'),
	help showimx, return
end
%Check image type and data format
nx = size(A.Data,1);
nz = A.Nz;
ny = A.Ny;
%set data range
baseRangeX = 1:nx; 
baseRangeY = 1:ny;
baseRangeZ = 1:nz;
%initialize left handside values
lhs1 = double( (baseRangeX-0.5)*A.Grid*A.ScaleX(1)+A.ScaleX(2) ); % x-range
lhs2 = double( (baseRangeY-0.5)*A.Grid*A.ScaleY(1)+A.ScaleY(2) ); % y-range
lhs3 = double(0);
lhs4 = double(0);
lhs5 = double(0);
lhs6 = double(0);
lhs7 = double(0);
if A.IType<=0, % grayvalue image format
    % Calculate frame range
    if frame>0 && frame<A.Nf,
        baseRangeY = baseRangeY + frame*ny;
    end
    lhs3 = double(A.Data(:,baseRangeY)');
    % Display image
	imagesc(lhs1,lhs2,lhs3);
elseif A.IType==2, % simple 2D vector format: (vx,vy)
	% Calculate vector position and components
	[lhs1,lhs2] = ndgrid(lhs1,lhs2);
	lhs3 = double(A.Data(:,baseRangeY   ))*A.ScaleI(1)+A.ScaleI(2);
	lhs4 = double(A.Data(:,baseRangeY+ny))*A.ScaleI(1)+A.ScaleI(2);
    if A.ScaleY(1)<0.0,
        lhs4 = -lhs4;
    end
	% quiver(lhs1,lhs2,lhs3,lhs4);
elseif (A.IType==3 || A.IType==1) , % normal 2D vector format + peak: sel+4*(vx,vy) (+peak)
	% Calculate vector position and components
    [lhs1,lhs2] = ndgrid(lhs1,lhs2);
    lhs3 = lhs1*0;
    lhs4 = lhs2*0;
    % Get choice
    maskData = double(A.Data(:,baseRangeY));
	% Build best vectors from choice field
	for i = 0:5,
		mask = ( maskData ==(i+1) );
		if (i<4) % get best vectors
			dat = double(A.Data(:,baseRangeY+(2*i+1)*ny));
			lhs3(mask) = dat(mask);
			dat = double(A.Data(:,baseRangeY+(2*i+2)*ny));
			lhs4(mask) = dat(mask);
		else    % get interpolated vectors
			dat = double(A.Data(:,baseRangeY+7*ny));
			lhs3(mask) = dat(mask);
			dat = double(A.Data(:,baseRangeY+8*ny));
			lhs4(mask) = dat(mask);
		end
    end
	lhs3 = lhs3*A.ScaleI(1)+A.ScaleI(2);
	lhs4 = lhs4*A.ScaleI(1)+A.ScaleI(2);
    %Display vector field
    if A.ScaleY(1)<0.0,
        lhs4 = -lhs4;
    end
	% quiver(lhs1,lhs2,lhs3,lhs4);
elseif A.IType==4,
	% Calculate vector position and components
    lhs3 = double((baseRangeZ-0.5)*A.Grid*A.ScaleZ(1)+A.ScaleZ(2));
    lhs4 = double(zeros([nx ny nz]));
    lhs5 = lhs4;
    lhs6 = lhs4;
    for iz=1:nz,
        lhs4(:,:,iz)=double( A.Data(:,baseRangeY+(iz-1)*ny) );
        lhs5(:,:,iz)=double( A.Data(:,baseRangeY+(iz+0)*ny) );
        lhs6(:,:,iz)=double( A.Data(:,baseRangeY+(iz+1)*ny) );
    end
	[lhs1,lhs2,lhs3] = ndgrid(lhs1,lhs2,lhs3);
	lhs4 = lhs4*A.ScaleI(1)+A.ScaleI(2);
	lhs5 = lhs5*A.ScaleI(1)+A.ScaleI(2);
	lhs6 = lhs6*A.ScaleI(1)+A.ScaleI(2);
elseif (A.IType==5 || A.IType==6),
	% Calculate vector position and components
    lhs3 = double((baseRangeZ-0.5)*A.Grid*A.ScaleZ(1)+A.ScaleZ(2));
    lhs4 = double(zeros([nx ny nz]));
    lhs5 = lhs4;
    lhs6 = lhs4;
    lhs7 = lhs4;
    blockSize = size(A.Data,2)/A.Ny;
    for iz=1:nz,
        px = double(zeros([nx ny]));
        py = double(zeros([nx ny]));
        pz = double(zeros([nx ny]));
        prange = baseRangeY + ((iz-1)*blockSize*ny);
        % Build best vectors from best choice field
        maskData = double( A.Data(:,prange));
        maskData = bitand( maskData, ones(size(maskData))*255);
        for i = 0:5,
            mask = (maskData==(i+1));
            if (i<4) % get best vectors
                dat = double(A.Data(:,prange+(3*i+1)*ny));
                px(mask) = dat(mask);
                dat = double(A.Data(:,prange+(3*i+2)*ny));
                py(mask) = dat(mask);
                dat = double(A.Data(:,prange+(3*i+3)*ny));
                pz(mask)=dat(mask);
            else    % get interpolated vectors
                dat = double(A.Data(:,prange+10*ny));
                px(mask) = dat(mask);
                dat = double(A.Data(:,prange+11*ny));
                py(mask) = dat(mask);
                dat = double(A.Data(:,prange+12*ny));
                pz(mask) = dat(mask);
            end
        end
        lhs4(:,:,iz) = px;
        lhs5(:,:,iz) = py;
        lhs6(:,:,iz) = pz;
        lhs7(:,:,iz) = maskData;
    end
	[lhs1,lhs2,lhs3] = ndgrid(lhs1,lhs2,lhs3);
	lhs4 = lhs4*A.ScaleI(1)+A.ScaleI(2);
	lhs5 = lhs5*A.ScaleI(1)+A.ScaleI(2);
	lhs6 = lhs6*A.ScaleI(1)+A.ScaleI(2);
end
end

