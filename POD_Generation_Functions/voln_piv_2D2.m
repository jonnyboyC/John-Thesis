function vol_frac = voln_piv_2D(x,y,bnd_idx,bnd_x, bnd_y)
% VOLN_PIV_2D determine the volumne contained in each flow pixel
%
% vol_frac = VOLN_PIV_2D(x,y,bnd_idx, bnd_x, bnd_y) determine the volumn contained in
% each flow pixel taking into account the boundary


vol_frac = zeros(size(x));

xlim = size(x,1);
ylim = size(x,2);

x_range = 1:size(x,1)-1;
y_range = 1:size(x,2)-1;

x_median = (x(x_range,y_range)+x(x_range,y_range+1)+x(x_range+1,y_range)+x(x_range+1,y_range+1))/4;
y_median = (y(x_range,y_range)+y(x_range,y_range+1)+y(x_range+1,y_range)+y(x_range+1,y_range+1))/4;

for i = 1:size(x,1)
    for j = 1:size(x,2)
        % In flow
        if bnd_idx(i,j) == 1 && bnd_x(i,j) == 0 && bnd_y(i,j) == 0
            width = ((x_median(i,j) - x_median(i-1,j)) + (x_median(i,j-1) - x_median(i-1,j-1)))/2;
            left_hand = y_median(i-1,j) - y_median(i-1,j-1);
            right_hand = y_median(i,j) - y_median(i,j-1);
            vol_frac(i,j) = width*(right_hand + left_hand)/2;
            
        % out of flow
        elseif bnd_idx(i,j) == -1
            vol_frac(i,j) = 0;
            
        % Top boundary
        elseif bnd_y(i,j) == -1 && j > 1 && i > 1 && i < xlim && ...
                bnd_y(i-1,j) == -1 && bnd_y(i+1,j) == -1 && bnd_idx(i,j-1) == 1
            width = x_median(i,j-1)-x_median(i-1,j-1);
            left_hand = y(i,j)-y_median(i-1,j-1);
            right_hand = y(i,j)-y_median(i,j-1);
            vol_frac(i,j) = width*(right_hand + left_hand)/2;
            
        % Bottom boundary
        elseif bnd_y(i,j) == 1 && j < ylim && i > 1 && i < xlim && ...
                bnd_y(i-1,j) == 1 && bnd_y(i+1,j) == 1 && bnd_idx(i,j+1) == 1
            width = x_median(i,j)-x_median(i-1,j);
            left_hand = y_median(i-1,j)-y(i,j);
            right_hand = y_median(i,j)-y(i,j);
            vol_frac(i,j) = (right_hand+left_hand)/2*width;
            
        % Left boundary
        elseif bnd_x(i,j) == 1 && i < xlim && j > 1 && j < ylim && ...
                bnd_x(i,j-1) == 1 && bnd_x(i,j+1) == 1 && bnd_idx(i+1,j) == 1
            width = x_median(i,j) - x(i,j);
            left_hand = y_median(i,j) - y_median(i,j-1);
            vol_frac(i,j) = left_hand*width;
            
        % Right boundary
        elseif bnd_x(i,j) == -1 && i > 1 && j > 1 && j < ylim && ...
                bnd_x(i,j-1) == -1 && bnd_x(i,j+1) == -1 && bnd_idx(i-1,j) == 1
            width = x(i,j) - x_median(i-1,j);
            left_hand = y_median(i-1,j) - y_median(i-1,j-1);
            vol_frac(i,j) = left_hand*width;
            
        % Upper Left Inside Corner
        elseif bnd_x(i,j) == 1 && bnd_y(i,j) == -1 && i < xlim && j > 1 && ...
                bnd_x(i,j-1) == 1 && bnd_y(i+1,j) == -1 && bnd_idx(i+1,j-1) == 1
            width = x_median(i,j-1) - x(i,j);
            left_hand = y(i,j) - y_median(i,j-1);
            vol_frac(i,j) = left_hand*width;
            
        % Upper Right Inside Corner
        elseif bnd_x(i,j) == -1 && bnd_y(i,j) == -1 && i > 1 && j > 1 && ...
                bnd_x(i,j-1) == -1 && bnd_y(i-1,j) == -1 && bnd_idx(i-1,j-1) == 1
            width = x(i,j) - x_median(i-1,j-1);
            left_hand = y(i,j) - y_median(i-1,j-1);
            vol_frac(i,j) = left_hand*width;
            
        % Lower Left Inside Corner
        elseif bnd_x(i,j) == 1 && bnd_y(i,j) == 1 && i < xlim && j < ylim && ...
                bnd_x(i,j+1) == 1 && bnd_y(i+1,j) == 1  && bnd_idx(i+1,j+1) == 1
            width = x_median(i,j) - x(i,j);
            left_hand = y_median(i,j) - y(i,j);
            vol_frac(i,j) = left_hand*width;
        
        % Lower Right Inside Corner
        elseif bnd_x(i,j) == -1 && bnd_y(i,j) == 1 && i > xlim && j > 1 && ...
                bnd_x(i,j+1) == -1 && bnd_y(i-1,j)== -1 &&  bnd_idx(i-1,j+1) == 1
            width = x(i,j) - x_median(i-1,j);
            left_hand = y_median(i-1,j) - y(i,j);
            vol_frac(i,j) = left_hand*width;
            
        % Upper Left Outside Corner
        elseif bnd_x(i,j) == 1 && bnd_y(i,j) == -1 && i < xlim && j > 1 && i > 1 && j < ylim &&...
                bnd_x(i,j+1) == 1 && bnd_y(i-1,j) == -1 && bnd_idx(i+1,j-1) == 1
            width = x_median(i,j) - x_median(i-1,j);
            left_hand = y_median(i-1,j) - y(i,j);
            right_hand = y_median(i,j) - y(i,j);
            vol_frac(i,j) = width*(right_hand + left_hand)/2;
            width = x(i,j) - x_median(i-1,j-1);
            left_hand = (y(i,j) - y_median(i-1,j-1));
            vol_frac(i,j) = vol_frac(i,j) + left_hand*width;
            
        % Upper Right Outside Corner
        elseif bnd_x(i,j) == -1 && bnd_y(i,j) == -1 && i > 1 && j > 1 && i < xlim && j < ylim &&...
                bnd_x(i,j+1) == -1 && bnd_y(i+1,j) == -1 && bnd_idx(i-1,j-1) == 1
            width = x_median(i,j) - x_median(i-1,j);
            left_hand = y_median(i-1,j) - y(i,j);
            right_hand = y_median(i,j) - y(i,j);
            vol_frac(i,j) = width*(right_hand + left_hand)/2;
            width = x_median(i,j-1) - x(i,j);
            right_hand = y(i,j) - y_median(i,j-1);
            vol_frac(i,j) = vol_frac(i,j) + right_hand*width;
            
        % Lower Left Outside Corner
        elseif bnd_x(i,j) == 1 && bnd_y(i,j) == 1 && i < xlim && j < ylim && j > 1 && i > 1 &&...
                bnd_x(i,j-1) == 1 && bnd_y(i-1,j) == 1  && bnd_idx(i+1,j+1) == 1
            width = x_median(i,j-1) - x_median(i-1,j-1);
            left_hand = y(i,j) - y_median(i-1,j-1);
            right_hand = y(i,j) - y_median(i,j-1);
            vol_frac(i,j) = width*(right_hand + left_hand)/2;
            width = (x(i,j) - x_median(i-1,j));
            left_hand = y_median(i-1,j) - y(i,j);
            vol_frac(i,j) = vol_frac(i,j) + left_hand*width;
            
        % Lower Right Outside Corner
        elseif bnd_x(i,j) == -1 && bnd_y(i,j) == 1 && i > xlim && j > 1 && i > 1 && j < ylim && ...
                bnd_x(i,j-1)== -1 && bnd_y(i+1,j)== 1 &&  bnd_idx(i-1,j+1) == 1
            width = x_median(i,j-1) - x_median(i-1,j-1);
            left_hand = y(i,j) - y_median(i-1,j-1);
            right_hand = y(i,j) - y_median(i,j-1);
            vol_frac(i,j) = width*(right_hand + left_hand)/2;
            width = x_median(i,j)-x(i,j);
            right_hand = y_median(i,j)-y(i,j);
            vol_frac(i,j) = vol_frac(i,j)+right_hand*width;
        else
            vol_frac(i,j) = 0;
        end
    end
end

end