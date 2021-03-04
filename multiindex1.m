function ind=multiindex1(n,d) % max n in each dimension

ind=(1:n)';
if d~=1
for k=1:d-1
    tmp3=[];
    for i=1:n
        tmp=i*ones(n^k,1);
        tmp2=[tmp ind];
        tmp3=[tmp3;tmp2];
    end
    ind=tmp3;
end
end

end