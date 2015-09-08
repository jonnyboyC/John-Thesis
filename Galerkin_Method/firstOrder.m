function [s1, s2, s3, s4, s5, s6, s7, s8, s9] = firstOrder(x1, x2, position)
% FIRSTORDER used in generate_stencile to calculate coefficients for first
% order finite difference methods, using analytical solution to Lagrange's
% equation of unequal spacing

if position == 1
    s1 = 0;
    s2 = 0;
    s3 = 0;
    s4 = 1/(x1-x2);
    s5 = -1/(x1-x2);
    s6 = 0;
    s7 = 0;
    s8 = 0;
    s9 = 0;
else
    s1 = 0;
    s2 = 0;
    s3 = 0;
    s4 = 0;
    s5 = 1/(x1-x2);
    s6 = -1/(x1-x2);
    s7 = 0;
    s8 = 0;
    s9 = 0;
end