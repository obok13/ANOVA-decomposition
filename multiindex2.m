function ind=multiindex2(n,d) % max of sum is n
if d==1
    ind=(1:n)';
else
    ind=ones(1,d);
    for i=0:n-d
        dividers = nchoosek(1:(i+d-1), d-1);
        ndividers = size(dividers, 1);
        b = cat(2, zeros(ndividers, 1), dividers, (i+d)*ones(ndividers, 1));
        ind = [ind;diff(b, 1, 2)];
    end
    ind(1,:)=[];
end
