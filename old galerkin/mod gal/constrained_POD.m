function [X] = constrained_POD(a,li_i,N,n,epsi)
 global lambda L epsilon

 L=li_i;
 epsilon = epsi;
for i=1:N;
 lambda(i) = mean(a(i,:).*a(i,:));
 end

 x0 = eye(n,n);
 if N > n
 x0(N,:) = 0;
 end
options = optimset('Algorithm','interior-point');
%  problem = createOptimProblem('fmincon', ...
%  'objective', @objective, ...
%  'nonlcon', @constraint, ...
%  'x0',x0);
%  [x,fval,EXITFLAG,OUTPUT,LAMBDA] = fmincon(problem)
[x,fval,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@objective,x0,[], [], [],[],[],[],@constraint, options);% 
OUTPUT.message

 X = x*(x'*x)^(-1/2);
 end