function rmc=odecoe(nm,dm,fcuh)
% Routine to arrenge the coefficients of the system of ODEs
%dm = input('desired number of modes = ');
mc=modcoef(nm,fcuh);
rmc=mc(1:dm,1:dm+1);
np=dm+1;
for i =1:dm
    for j = i:dm
        np=np+1;
        rmc(:,np)=mc(1:dm,i/2*(2*nm-i+3)+1+(j-i));
    end
end
% [T3,Y3] = ode45(@eqdifvecm,t1(1:end),psitc(1,1:dm),[],-rmc);
% figure
% plot(t1,psitc(:,1))
% hold on
% plot(T3,Y3(:,1),'ro-')