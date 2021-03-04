clear

% parameters
d = 3;
cutdim=2;
p = 3; % polynomial order

% GMM model
nmix = 2;
w = ones(nmix,1)/nmix;
mu(1,:) = ones(1,d); mu(2,:) = -ones(1,d);
Sigma = zeros(d,d,nmix);
for idxmix = 1:nmix 
    Sigma(:,:,idxmix) = eye(d)/10 + ones(d)/100;
end
GMobj.mu = mu;
GMobj.sigma = Sigma;
GMobj.w = w;
GMobj.nmix = size(mu,1);
GMobj.dim = size(mu,2);
GMobj.type = 'Gaussian Mixture';

% train and test data
rng(0);
nsamples = 300;
X_train = random(gmdistribution(GMobj.mu,GMobj.sigma,GMobj.w), nsamples);
X_test = random(gmdistribution(GMobj.mu,GMobj.sigma,GMobj.w), nsamples);
realf = func(X_test,d);
realf_train = func(X_train,d);

% construct basis function
Muidx = (0:p)+(1:p)';
for idxd=1:d
    mu = mean(X_train(:,idxd).^(0:2*p-1),1);
    Mu = mu(Muidx);
    for idxod=1:p
        Mur = zeros(idxod+1);
        Mur(1:idxod,1:idxod+1) = Mu(1:idxod,1:idxod+1);
        Mur(idxod+1,idxod+1) = 1;
        e1 = zeros(idxod+1,1);
        e1(idxod+1) = 1;
        polycoeff = flipud(Mur\e1);
        phi_train{idxd}(:,idxod) = polyval(polycoeff,X_train(:,idxd));
        phi_test{idxd}(:,idxod) = polyval(polycoeff,X_test(:,idxd));
        phi_train{idxd}(:,idxod) = phi_train{idxd}(:,idxod)/sqrt(mean(phi_train{idxd}(:,idxod).^2));
        phi_test{idxd}(:,idxod) = phi_test{idxd}(:,idxod)/sqrt(mean(phi_test{idxd}(:,idxod).^2));
    end
end

for idxd=1:cutdim
    S{idxd}=nchoosek((1:d),idxd);
end

DMORPHf=mean(realf);
resf=realf_train-mean(realf_train);

count_basis=0; % # of basis
count_cp=1;

for idxd=1:cutdim

for idxset=1:size(S{idxd},1)
    clear ho Sr
    u=S{idxd}(idxset,:); % 크기가 idxd인 subset 하나 선택
    count_ho=1;
    ho(:,count_ho) = ones(nsamples,1);
    
    % size of SS is 1
    for idxelt = 1:idxd
        for idxod = 1:p
            count_basis=count_basis+1;
            Phi_train(:,count_basis)=phi_train{u(idxelt)}(:,idxod);
            Phi_test(:,count_basis)=phi_test{u(idxelt)}(:,idxod);
            if idxd > 1
                count_ho=count_ho+1;
                ho(:,count_ho) = Phi_train(:,count_basis);
            end
        end
    end

    % size of SS is 2
    for idxd2=2:idxd
        SS=nchoosek(u,idxd2);
        for idxset2=1:size(SS,1)
            uu=SS(idxset2,:);
%             mulind=multiindex1(p,idxd2);
            mulind=multiindex2(p,idxd2);
            for idxmul = 1:size(mulind,1)
                count_basis=count_basis+1;
                Phi_train(:,count_basis)=ones(nsamples,1);
                Phi_test(:,count_basis)=ones(nsamples,1);
                for idxelt = 1:idxd2
                    Phi_train(:,count_basis)=Phi_train(:,count_basis).*phi_train{uu(idxelt)}(:,mulind(idxmul,idxelt));
                    Phi_test(:,count_basis)=Phi_test(:,count_basis).*phi_test{uu(idxelt)}(:,mulind(idxmul,idxelt));
                end
                Phi_train(:,count_basis)=Phi_train(:,count_basis)/norm(Phi_train(:,count_basis));
                Phi_test(:,count_basis)=Phi_test(:,count_basis)/norm(Phi_test(:,count_basis));
                if idxd2 < idxd
                    count_ho=count_ho+1;
                    ho(:,count_ho) = Phi_train(:,count_basis);
                end
            end
        end
    end
    count_cp=count_cp+1;
    cp(count_cp) = count_basis;
    
    % construct B
    for idxho=1:count_ho
        Sr(:,idxho)=mean(Phi_train(:,cp(count_cp-1)+1:cp(count_cp)).*ho(:,idxho),1)';
    end
    if count_cp == 2
        B=Sr*Sr';
    else
        B=blkdiag(B,Sr*Sr');
    end
end
end

% construct P
A=Phi_train'*Phi_train;
t=size(A,1);
c=pinv(Phi_train)*resf;
P=eye(t)-pinv(A)*A;

% construct c_inf
[U,D,V]=svd(P*B);
r=rank(D,1e-8);
if r<t
    c_inf=V(:,r+1:t)*inv(U(:,r+1:t)'*V(:,r+1:t))*U(:,r+1:t)'*c;
else
    c_inf=c;
end

DMORPHf = DMORPHf + Phi_test*c_inf;

error = norm(realf - DMORPHf)/norm(realf)

figure;
scatter(realf,DMORPHf,'.'); hold on;
MAX=max(max([realf DMORPHf]));
MIN=min(min([realf DMORPHf]));
plot([MIN MAX],[MIN MAX])
xlabel('true f')
ylabel('DMORPH f')

function output = func(x,d)

mu1=0.5; mu2=0.5; mu3=0.5;
a0=1; a1=2; b0=2; b1=3; c0=3; c1=1; c2=2;
d0=1; d1=2; d2=2; d3=3;

output = (a1*(x(:,1)-mu1)+a0).*(b1*(x(:,2)-mu2)+b0)...
    +c2*(x(:,2)-mu2).^2+c1*(x(:,2)-mu2)+c0...
    +d3*(x(:,3)-mu3).^3+d2*(x(:,3)-mu3).^2+d1*(x(:,3)-mu3)+d0;

end