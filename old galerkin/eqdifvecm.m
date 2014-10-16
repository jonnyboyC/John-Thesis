function dy = eqdifvecm(t,y,fc)
nm=size(fc,1);
nt=size(fc,2);
dy = zeros(nm,1);    % a column vector
for i = 1:nm
    l=1;
    dy(i) = fc(i,1);
    for j= 1:nm
        l=l+1;
        dy(i) = dy(i)+fc(i,l)*y(j);
    end
    for j= 1:nm
        for k=j:nm
            l=l+1;
            dy(i) = dy(i)+fc(i,l)*y(j)*y(k);
        end
    end
end