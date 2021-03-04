clear;

d = 4; % dimension of input
od = 4; % order of gpc
q = 2; % cutoff dimension
thres = 0.001; % threshold in adaptive criterion
N_test = 1000; % # of test samples

[x,w] = knots_Legendre(od+1,0,1); % domain: [0,1]^d, measure: lesbesgue
qx = repmat(x',1,d);
qw = repmat(w',1,d);

% basiscoeff of gpc
syms x;
basiscoeff = cell(d,od);
for idxd = 1:d
for idxod = 1:od
    basiscoeff{idxd,idxod} = sym2poly(legendreP(idxod,2*x-1));
    alp=sum(polyval(basiscoeff{idxd,idxod},qx(:,idxd)).^2.*qw(:,idxd));
    basiscoeff{idxd,idxod} = basiscoeff{idxd,idxod}/sqrt(alp); % normalization
end
end
clear x;

for idxd=1:q
    Set{idxd}=nchoosek(1:d,idxd);
end
anc = 0.5*ones(1,d); % anchor pt
fu.set{1} = []; % u
fu.qx{1} = anc; % quad pts of [0,1]^u
fu.qw{1} = 1; % quad weights of [0,1]^u
fu.val{1} = func(fu.qx{1},d); % f_u value on fu.qx
fu.mean{1} = fu.val{1}; % Ef_u
fu.var{1}= 0; % Varf_u
sizeSet(1)=1;
sumvar=0;
fval = 1; % # of function call

for sizeu=1:q
    
    idxfu = sizeSet(sizeu)+1;
    nodidx_pt=multiindex1(od+1,sizeu); % list of quadrature points
    for idxu=1:size(Set{sizeu},1)
        u=Set{sizeu}(idxu,:);
        fu.set{idxfu,1}=u;
        fu.qx{idxfu,1}=ones((od+1)^sizeu,1)*anc;
        clear tmpw;
        for idxsizeu=1:sizeu
             fu.qx{idxfu}(:,u(idxsizeu))=qx(nodidx_pt(:,idxsizeu),u(idxsizeu));
             tmpw(:,idxsizeu)=qw(nodidx_pt(:,idxsizeu),u(idxsizeu));
        end
        fu.qw{idxfu,1} = prod(tmpw,2);
        fu.val{idxfu,1} = func(fu.qx{idxfu},d);
        fval = fval + size(fu.qx{idxfu,1},1);
        fu.val{idxfu,1}=fu.val{idxfu} - fu.val{1};
        if sizeu > 1
            for sizev = 1:sizeu-1
                for idxfu2 = sizeSet(sizev)+1:sizeSet(sizev+1)
                    if all(ismember(fu.set{idxfu2},u)) % if fu.set{idxfu2} \subset u
                        tmpidx = nodidx_pt(:,ismember(u,fu.set{idxfu2}))-1;
                        tmpidx = sum(tmpidx.*((od+1).^flip(0:sizev-1)),2)+1;
                        fu.val{idxfu,1}=fu.val{idxfu,1} - fu.val{idxfu2,1}(tmpidx,:);
                    end
                end
            end
        end
        fu.mean{idxfu,1} = fu.qw{idxfu}'*fu.val{idxfu};
        fu.var{idxfu,1} = fu.qw{idxfu}'*(fu.val{idxfu}.^2);
        fu.var{idxfu,1} = fu.var{idxfu} - fu.mean{idxfu}.^2;
        fu.var{idxfu,1} = norm(fu.var{idxfu,1});
        sumvar = sumvar + fu.var{idxfu};
        idxfu = idxfu + 1;
    end
    sizeSet(sizeu+1) = idxfu - 1;
    
    if sizeu<q && size(Set{sizeu},1)~=0
        for idxfu=sizeSet(sizeu)+1:sizeSet(sizeu+1)
            if fu.var{idxfu}/sumvar < thres
                for idxd=sizeu+1:q
                    idxset=1;
                    while idxset <= size(Set{idxd},1)
                        if all(ismember(fu.set{idxfu},Set{idxd}(idxset,:))) % if fu.set{idxfu} \subset Set{idxd}(idxset,:)
                            Set{idxd}(idxset,:)=[];
                        else
                            idxset=idxset+1;
                        end
                    end
                end
            end
        end
    end
    
end

fprintf('fval: %d\n',fval);

% gpc coefficients
idxcoeff=1;
gpccoeff{idxcoeff,1} = 0;
for idxfu = 1:size(fu.set,1)
    gpccoeff{idxcoeff} = gpccoeff{idxcoeff} + fu.mean{idxfu,1};
end
idxcoeff = 2;
for sizeu = 1:q
    nodidx_basis = multiindex2(od,sizeu);
    for idxu = 1:size(Set{sizeu},1)
        u = Set{sizeu}(idxu,:);
        for idxnod = 1:size(nodidx_basis,1)
            gpccoeff{idxcoeff,1} = 0;
            for idxfu = sizeSet(sizeu)+1:size(fu.set,1)
                u2 = fu.set{idxfu};
                sizeu2 = size(u2,2);
                nodidx_pt = multiindex1(od+1,sizeu2); % list of quadrature pts
                if all(ismember(u,u2)) % if u \subset u2
                    clear tmpw;
                    idxsizeu = 1;
                    for idxsizeu2 = 1:sizeu2
                        if ~ismember(u2(idxsizeu2),u) % if u2(idxsizeu2) \not\in u
                            tmpw(:,idxsizeu2)=qw(nodidx_pt(:,idxsizeu2),u2(idxsizeu2));
                        else
                            tmpw(:,idxsizeu2)=polyval(basiscoeff{u2(idxsizeu2),nodidx_basis(idxnod,idxsizeu)},fu.qx{idxfu}(:,u2(idxsizeu2))).*qw(nodidx_pt(:,idxsizeu2),u2(idxsizeu2));
                            idxsizeu = idxsizeu + 1;
                        end
                    end
                    polyw = prod(tmpw,2);
                    gpccoeff{idxcoeff} = gpccoeff{idxcoeff} + polyw'*fu.val{idxfu};
                end
            end
            idxcoeff = idxcoeff + 1;
        end
    end
end

% campare results
X_test = rand(N_test,d);
approx = ones(size(X_test,1),1)*gpccoeff{1};
approxmean = gpccoeff{1};
approxvar = 0;
idxcoeff = 2;
for sizeu = 1:q
    nodidx_basis = multiindex2(od,sizeu);
    for idxu = 1:size(Set{sizeu},1)
        u = Set{sizeu}(idxu,:);
        for idxnod = 1:size(nodidx_basis,1)
            Phi = ones(size(X_test,1),1);
            for idxsizeu = 1:sizeu
                Phi = Phi.*polyval(basiscoeff{u(idxsizeu),nodidx_basis(idxnod,idxsizeu)},X_test(:,u(idxsizeu)));
            end
            approx = approx + Phi*gpccoeff{idxcoeff};
            approxvar = approxvar + gpccoeff{idxcoeff}.^2;
            idxcoeff = idxcoeff + 1;
        end
    end
end
exact = func(X_test,d);
exactmean = mean(exact,1);
exactvar = mean(exact.^2,1);
exactvar = exactvar - exactmean.^2;

errormean = norm(approxmean-exactmean)/norm(exactmean);
errorvar = norm(approxvar-exactvar)/norm(exactvar);
error = norm(approx-exact)/norm(exact);

fprintf('error: %.3e\n',error);
fprintf('error_mean: %.3e\n',errormean);
fprintf('error_var: %.3e\n',errorvar);

function output = func(x,d)

output = (1+1/d)^d*prod(x.^(1/d),2);

end