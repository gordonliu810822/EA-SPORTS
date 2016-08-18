function obj = sparseLDApath(X,Z,Xte,Zte,clus,lamseq, opts)
% Pko: that is orthogonal to (Y'Y)^{-1/2}Y'Xbi, i<k
if nargin < 7
    opts = [];
    opts.verbose = 1;    
    opts.maxIters = 30; 
    opts.tol = 1e-10;
    %opts.gamma = 2.7;
    opts.standardize = 0; % standardized: 1 or not: 0.
end

[n,K] = size(Z); %K: # of clusters.
p = size(X,2);   %p: # of feature.
nlam = length(lamseq);

Zclass = zeros(n,1);
for k = 1:K
    Zclass(Z(:,k)==1) =k;
end

mu = zeros(K,p);
sumZvec = sum(Z,1);

%% precomputing for standardization of X, Xte
if opts.standardize ~= 1
    tmp = zeros(p,p);
    for k = 1:K
        for j = 1:p
            mu(k,j) = sum(Z(:,k).*X(:,j))/sumZvec(k);
        end
        tmp = tmp + ((X - repmat(mu(k,:),n,1)).*repmat(Z(:,k),1,p).^0.5)'*((X - repmat(mu(k,:),n,1)).*repmat(Z(:,k),1,p).^0.5);
    end
    SigmaW = tmp/n;
    
    sigma2j = diag(SigmaW);
    if min(sigma2j) < 1e-10
       error('Some features have 0 within-class standard deviation.');
    end
    Xm = mean(X,1);
    X1 = X - repmat(Xm,n,1);
    X = X1./repmat(sigma2j.^0.5',n,1);
    if isempty(Xte) ~= 1
        ntest = size(Xte,1);
        tmp = Xte - repmat(Xm,ntest,1);
        Xte = tmp./repmat(sigma2j.^0.5',ntest,1);
    end
    if isempty(Zte) ~= 1
        ntest = size(Zte,1);
        SigmaBt = Xte'*Zte*((Zte'*Zte)\Zte'*Xte)/ntest;
        tmp = zeros(p,p);
        sumZtevec = sum(Zte,1);
        mute = zeros(K,p);
        for k = 1:K
            for j = 1:p
                mute(k,j) = sum(Zte(:,k).*Xte(:,j))/sumZtevec(k);
            end
            tmp = tmp + ((Xte - repmat(mute(k,:),ntest,1)).*repmat(Zte(:,k),1,p).^0.5)'*((Xte - repmat(mute(k,:),ntest,1)).*repmat(Zte(:,k),1,p).^0.5);
        end
        SigmaWt = tmp/ntest;
        % SigmaWt = Xte'*Xte/ntest;
    end
end

A = X'*Z*((Z'*Z)^0.5\eye(K))/sqrt(n);
 
  
%% set initial cell
funcVal = cell(nlam,1);
discr  = cell(nlam,1);
Xproj  = cell(nlam,1);
Xteproj= cell(nlam,1);
if isempty(Xte) ~= 1
    Zpred  = cell(nlam,1);
end

crit = zeros(nlam,1);
%% update path
for i = 1:nlam
    funcVal{i} = zeros(opts.maxIters,clus);
    funcVal{i}(1,:) = -inf;
    discr{i} = zeros(p,clus);
    lam = lamseq(i);
    
    for k = 1:clus
        
        if k == 1
            Pko = eye(K);
        elseif k > 1
            [U,S,~] = svd(A'*discr{i});
            u_vec = U(:,diag(S)>1e-10);
            Pko = eye(K) - u_vec*u_vec';
        end
        [U,S,~] = svd(A*Pko);
        d = S(1,1)^2;
        discr{i}(:,k) = U(:,1);
        funcVal{i}(1,k) = discr{i}(:,k)'*A*Pko*A'*discr{i}(:,k) - lam*d*sum(abs(discr{i}(:,k)));
        
        
        for iter = 2:opts.maxIters
            tmp = A*Pko*(A'*discr{i}(:,k));
            dif = abs(tmp) - lam*d/2;
            dif(dif < 0) = 0;
            S = sign(tmp).*dif;
            %S = S/(1-1/opts.gamma);
            %S(abs(tmp) > opts.gamma*lam*d/2) = tmp(abs(tmp) > opts.gamma*lam*d/2);
            discr{i}(:,k) = S;
            delta = sqrt(discr{i}(:,k)'*discr{i}(:,k));
            if delta < 1e-5
                discr{i}(:,k) = zeros(p,1);
            else
                discr{i}(:,k) = discr{i}(:,k)/delta;
            end
            funcVal{i}(iter,k) = discr{i}(:,k)'*A*Pko*A'*discr{i}(:,k) - lam*d*sum(abs(discr{i}(:,k)));
            
            if opts.verbose
                fprintf('Iter %d:   funcVal = %f \n', iter, funcVal{i}(iter,k));
            end
            
            if iter > 1
                if abs(funcVal{i}(iter,k)-funcVal{i}(iter-1,k)) <= opts.tol
                    break
                end
            end
        end 
    end
    Xproj{i} = X*discr{i}; % dim: n-by-clus 
    if isempty(Xte) ~= 1
        Xteproj{i} = Xte*discr{i};% dim: ntest-by-clus 
        Zpred{i} = zeros(ntest,clus);
        for k = 1:clus
            Zpred{i}(:,k) = classify(Xproj{i}(:,1:k),Xteproj{i}(:,1:k),Zclass);
        end
        
        crit(i)=trace(discr{i}'*SigmaBt*discr{i});
    end

end

obj.discr = discr;
obj.funcVal = funcVal;
obj.crit = crit; 
obj.Xproj = Xproj;
if isempty(Xte) ~= 1
    obj.Zpred = Zpred;
end
end

function pred = classify(xtr,xte,ytr,equalpriors)
if nargin < 4
    equalpriors = false;
end
p = size(xtr,2);
ntest = size(xte,1);
clus = length(unique(ytr));
prior = repmat(1/clus,clus,1);
if ~equalpriors
    for k = 1:clus
        prior(k) = mean(ytr==k);
    end
end
%% classify obs to nearest centroid in training
mu = zeros(p,clus);
for k = 1:clus
    mu(:,k) = mean(xtr(ytr==k,:),1);
end
dd = xte*mu;
dd0= ones(1,p)*mu.^2/2 - log(prior');
negdists = dd - repmat(dd0,ntest,1);
[~,pred] = max(negdists,[],2);
end


