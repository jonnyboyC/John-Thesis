function niu=visc_dis_n_new(psi,nm,landa,D0,Dk,C0,Ck,Ckk,re)
nt = size(psi,1);
nl=size(C0,1);
d=zeros(nm,nt);
P=zeros(nm,nt);

% for i = 1: nl
%     psi(:,i)=psi(:,i)/landa(i);
% end

    for l = 1:nt
    for i = 1:nm
        d(i,l)=D0(i);
        for j = 1:nm
            d(i,l)=d(i,l)+Dk(i,j)*psi(l,j);
        end
    end
end



for l = 1:nt
    for i = 1:nm
        ckk=reshape(Ckk(i,:),nl,nl);
        for j = nm+1:nl
            P(i,l)=P(i,l)+Ck(i,j)*psi(l,j)+Dk(i,j)*psi(l,j)/re;
            for j1 = nm+1:nl
                P(i,l)=P(i,l)+ckk(j,j1)*psi(l,j)*psi(l,j1);
            end
        end
    end
end

vs1=mean(P.*d,2)./mean(d.^2,2);
niu=mean(psi(:,1:nm)'.*P,2)./mean(psi(:,1:nm)'.*d,2);
vs2=mean((psi(:,1:nm)').^2.*P.*d,2)./mean((psi(:,1:nm)'.*d).^2,2);
niu=niu';
end
 
