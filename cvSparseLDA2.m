function out = cvSparseLDA2(X,Z, opts)
% Pko: that is orthogonal to (Y'Y)^{-1/2}Y'Xbi, i<k
if nargin < 3
    opts = [];
    opts.maxIters = 100; 
    opts.tol = 1e-10;
    opts.epsilon = 0.001;
    opts.nlam = 100;
    opts.maxLam = 2;
    opts.nfold = 5;
end
clus = size(Z,2)-1;

lamseq = logspace(log10(opts.epsilon*opts.maxLam),log10(opts.maxLam), opts.nlam);

% mySeed = 10;
% rng(mySeed);  
indices = crossvalind('Kfold', size(Z,1), opts.nfold);
error = zeros(opts.nfold,opts.nlam,clus);
% save('indices.dat', 'indices','-ascii','-tabs');
for fold = 1: opts.nfold
    Xtrain = X(indices ~= fold,:);
    Xtest = X(indices == fold,:);
    Ztrain = Z(indices ~= fold,:);
    Ztest = Z(indices == fold,:);
    ZclassTe = zeros(size(Ztest,1),1);
    for k = 1:size(Z,2)
        ZclassTe(Ztest(:,k)==1) =k;
    end
    
    obj = sparseLDApath(Xtrain,Ztrain,Xtest,Ztest,clus,lamseq);
    for i = 1:opts.nlam
        for k = 1:clus            
            % error(fold,i,k)=sum(obj.Zpred{i}(:,k) ~= ZclassTe);
             error(fold,i,k) = -obj.crit(i);
        end
    end
    
end
out.err = error;
out.errmean = reshape(mean(error,1),opts.nlam,clus);
out.lamseq = lamseq;
out.onese = reshape(var(error).^0.5/sqrt(opts.nfold),opts.nlam,clus);
[C,~]=min(out.errmean);
[~,out.bestK]=min(C);
out.bestLam1se = max(lamseq(out.errmean(:,out.bestK) <= min(out.errmean(:,out.bestK)+out.onese(:,out.bestK))));
end


    