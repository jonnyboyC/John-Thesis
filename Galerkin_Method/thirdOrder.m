function [s1, s2, s3, s4, s5, s6, s7, s8, s9] = thirdOrder(x1, x2, x3, x4, position)
% THIRDORDER third order finite difference methods for GENERATE_STENCIL
if position == 1
    x  = x1;
    s1 = 0;
    s2 =  (x2*x3 - 2*x*x3 - 2*x*x4 - 2*x*x2 + x2*x4 + x3*x4 + 3*x^2)...
        /((x1 - x2)*(x1 - x3)*(x1 - x4));
    s3 = -(x1*x3 - 2*x*x3 - 2*x*x4 - 2*x*x1 + x1*x4 + x3*x4 + 3*x^2)...
        /((x1 - x2)*(x2 - x3)*(x2 - x4));
    s4 =  (x1*x2 - 2*x*x2 - 2*x*x1 - 2*x*x4 + x1*x4 + x2*x4 + 3*x^2)...
        /((x1 - x3)*(x2 - x3)*(x3 - x4));
    s5 = -(x1*x2 - 2*x*x2 - 2*x*x3 - 2*x*x1 + x1*x3 + x2*x3 + 3*x^2)...
        /((x1 - x4)*(x2 - x4)*(x3 - x4));
    s6 = 0;
    s7 = 0;
    s8 = 0;
    s9 = 0;
elseif position == 2
    x  = x2;
    s1 = 0;
    s2 = 0;
    s3 =  (x2*x3 - 2*x*x3 - 2*x*x4 - 2*x*x2 + x2*x4 + x3*x4 + 3*x^2)...
        /((x1 - x2)*(x1 - x3)*(x1 - x4));
    s4 = -(x1*x3 - 2*x*x3 - 2*x*x4 - 2*x*x1 + x1*x4 + x3*x4 + 3*x^2)...
        /((x1 - x2)*(x2 - x3)*(x2 - x4));
    s5 =  (x1*x2 - 2*x*x2 - 2*x*x1 - 2*x*x4 + x1*x4 + x2*x4 + 3*x^2)...
        /((x1 - x3)*(x2 - x3)*(x3 - x4));
    s6 = -(x1*x2 - 2*x*x2 - 2*x*x3 - 2*x*x1 + x1*x3 + x2*x3 + 3*x^2)...
        /((x1 - x4)*(x2 - x4)*(x3 - x4));
    s7 = 0;
    s8 = 0;
    s9 = 0;
elseif position == 3
    x  = x3;
    s1 = 0;
    s2 = 0;
    s3 = 0;
    s4 =  (x2*x3 - 2*x*x3 - 2*x*x4 - 2*x*x2 + x2*x4 + x3*x4 + 3*x^2)...
        /((x1 - x2)*(x1 - x3)*(x1 - x4));
    s5 = -(x1*x3 - 2*x*x3 - 2*x*x4 - 2*x*x1 + x1*x4 + x3*x4 + 3*x^2)...
        /((x1 - x2)*(x2 - x3)*(x2 - x4));
    s6 =  (x1*x2 - 2*x*x2 - 2*x*x1 - 2*x*x4 + x1*x4 + x2*x4 + 3*x^2)...
        /((x1 - x3)*(x2 - x3)*(x3 - x4));
    s7 = -(x1*x2 - 2*x*x2 - 2*x*x3 - 2*x*x1 + x1*x3 + x2*x3 + 3*x^2)...
        /((x1 - x4)*(x2 - x4)*(x3 - x4));
    s8 = 0;
    s9 = 0;
else
    x  = x4;
    s1 = 0;
    s2 = 0;
    s3 = 0;
    s4 = 0;
    s5 =  (x2*x3 - 2*x*x3 - 2*x*x4 - 2*x*x2 + x2*x4 + x3*x4 + 3*x^2)...
        /((x1 - x2)*(x1 - x3)*(x1 - x4));
    s6 = -(x1*x3 - 2*x*x3 - 2*x*x4 - 2*x*x1 + x1*x4 + x3*x4 + 3*x^2)...
        /((x1 - x2)*(x2 - x3)*(x2 - x4));
    s7 =  (x1*x2 - 2*x*x2 - 2*x*x1 - 2*x*x4 + x1*x4 + x2*x4 + 3*x^2)...
        /((x1 - x3)*(x2 - x3)*(x3 - x4));
    s8 = -(x1*x2 - 2*x*x2 - 2*x*x3 - 2*x*x1 + x1*x3 + x2*x3 + 3*x^2)...
        /((x1 - x4)*(x2 - x4)*(x3 - x4));
    s9 = 0;
end
