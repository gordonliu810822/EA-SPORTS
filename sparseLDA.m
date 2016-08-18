function obj = sparseLDA(X,Z,clus,lam, opts)
% Pko: that is orthogonal to (Y'Y)^{-1/2}Y'Xbi, i<k
if nargin < 5
    opts = [];
    opts.verbose = 1;    
    opts.maxIters = 100; 
    opts.tol = 1e-10;
    opts.gamma = 10000;
    opts.standardize = 0; % standardized: 1 or not: 0.
end

[n,K] = size(Z); %K: # of clusters.
p = size(X,2);   %p: # of feature.

mu = zeros(K,p);
sumZvec = sum(Z,1);

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
    X1 = X - repmat(mean(X,1),n,1);
    X = X1./repmat(sigma2j.^0.5',n,1);
end

A = X'*Z*((Z'*Z)^0.5\eye(K))/sqrt(n);
 
discr = zeros(p,clus);
  
%% update vector
funcVal = zeros(opts.maxIters,clus);
funcVal(1,:) = -inf;

for k = 1:clus

    if k == 1
        Pko = eye(K);
    elseif k > 1
        [U,S,~] = svd(A'*discr);
        u_vec = U(:,diag(S)>1e-10);
        Pko = eye(K) - u_vec*u_vec';
    end
    [U,S,~] = svd(A*Pko);
    d = S(1,1)^2;
    discr(:,k) = U(:,1);
    funcVal(1,k) = discr(:,k)'*A*Pko*A'*discr(:,k) - lam*d*sum(abs(discr(:,k)));
        
    
    for iter = 2:opts.maxIters
        tmp = A*Pko*(A'*discr(:,k));
        dif = abs(tmp) - lam*d/2;
        dif(dif < 0) = 0;
        discr(:,k) = sign(tmp).*dif;
        delta = sqrt(discr(:,k)'*discr(:,k));
        if delta < 1e-5
            discr(:,k) = zeros(p,1);
        else 
            discr(:,k) = discr(:,k)/delta;
        end
        funcVal(iter,k) = discr(:,k)'*A*Pko*A'*discr(:,k) - lam*d*sum(abs(discr(:,k)));
        
        if opts.verbose
            fprintf('Iter %d:   funcVal = %f \n', iter, funcVal(iter,k));
        end
        
        if iter > 1
            if abs(funcVal(iter,k)-funcVal(iter-1,k)) <= opts.tol
                break
            end
        end
    end
    
end

obj.discr = discr;
obj.funcVal = funcVal;


    