function niu=viscious_dis_old(psi,nm,landa,b,d,e)
% for i = 1:nm
%     for j = 1:nm
%         ac(i,j)=mean(psi(:,i).*psi(:,j));
%     end
% end

for k = 1:nm
    nv=1;   
    for i = 1:nm
        for j = 1:nm
            nc(k,i,j)=e(k,nv)*mean(psi(:,i).*psi(:,j).*psi(:,k));
            nv=nv+1;
        end
    end
    mf_e(k)=d(k,k)*landa(k);
    
    bc(k)=b(k,k)*landa(k);
    qc(k)=sum(sum(nc(k,:,:)));
    uc(k)=qc(k)+mf_e(k);
    niu(k) = -uc(k)/bc(k);
end