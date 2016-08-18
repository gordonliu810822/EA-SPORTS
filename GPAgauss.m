function obj = GPAgauss(gwasPval, annMat, opts)

%% 
% GPA -- Gwas prioritization using Pleiotropy information and tissue-specific
% expressions.
% input
% pvalue: a ngene-by-nGWAS pvalue matrix
% ann: a ngene-by-nAnn matrix, whose elements are normalized in tissues.

% opts: inclues initialization and some other parameters.
% output
% obj
% obj.loglik: the loglike during the iteractions
% obj.pi_vec: final Pi vector, e.g., pi_00, pi_10, pi_01, pi_11.
% obj.piMat: all the Pi vector during the iteractions
% obj.betaAlpha: final alpha parameter for the beta distr
% obj.betaAlphaMat: all alpha parameter during iteractions
% obj.Z: The posterior of the the hidden variables, e.g., Z_00, Z_10, Z_01,
% Z_11.
% obj.maxDiff: max diff of parameters change.
% obj.mu: mu an nComb-by-nAnn matrix
% obj.muMat: all mu during iteractions.
%
%
%%
if nargin < 3
    opts = [];
    opts.verbose = 1;
    
    opts.maxIters = 2000; 
    opts.epsStopLogLik = 1e-10;
    opts.epsStopPE = 1e-10;
    opts.initBetaMean = 0.1;
    opts.initPi = 0.1;
    opts.initQ1 = 0.75;
    opts.lbPi1 = 0.001;
    opts.lbBetaAlpha = 0.001;
    opts.lbQ1 = 0.001;
    opts.hom =0; %indicator whehter the homogeneous variance is used.
    opts.constraintPi = 0; %indicator whehter the constrain model on pi is used.
    opts.constraintMu = 0; %indicator whehter the constrain model on mu is used.
    opts.empiricalNull = 0;%if 1, will use beta(alpha0,1) for null distribution.
end


%%
%check pvalue matrix 
if ( (sum(sum(gwasPval<=0))>0) || (sum(sum(gwasPval>1))>0))
    error('InValidate pvalue matrix: There exists some p-values: p<=0 or p>1 !');
end


%% Get nsnp, nGWAS, nAnn and check compatibility
[nsnp1, nGWAS] = size(gwasPval);



%%
% check annMat
if (~isempty(annMat))
   [nsnp2, nAnn] = size(annMat);
   if nsnp1 ~=nsnp2
       error('the length of pvalue does not equal to the length of annotation');
   else
       nsnp = nsnp1;
       fprintf('Info: Annotation data is provided.\n');
       fprintf('Info: Number of GWAS: %d; Number of Annotation data: %d. \n', nGWAS, nAnn);
   end
    %index0 = (annMat==0);
    %index1 = (annMat==1);
    %if (sum(index0(:)) + sum(index1(:)) ~= nsnp2*nAnn)
    %    error('Invalidate annMat: there are some entries in annMat do not belong to {0,1}.');
    %end
    %clear index0 index1 nsnp1 nsnp2;
    clear nsnp1 nsnp2;
else
    nsnp = nsnp1;
    clear nsnp1;
    
    fprintf('Info: no Annotation data is provided.\n');
    fprintf('Info: Number of GWAS: %d. \n', nGWAS);
end

%%define a binary matrix for multiple GWAS data sets
%
binaryMat = comb_state(nGWAS);

nComb = 2^nGWAS;
%% Pre-computation

logPval = log(gwasPval);

% initialization of parameters

betaAlpha = opts.initBetaMean/(1-opts.initBetaMean)*ones(1,nGWAS);

pi_vec = opts.initPi.^sum(binaryMat,2);
pi_vec(pi_vec<opts.lbPi1) = opts.lbPi1;
pi_vec(1) = 1-sum(pi_vec(2:end));

if (~isempty(annMat) )
    if opts.constraintMu == 1
        mu = opts.initQ1*ones(1,nAnn);
        %mu: nComb-by-nAnn
        tmp = (annMat - repmat(mu,nsnp,1))'*(annMat - repmat(mu,nsnp,1));
        if opts.hom == 1 
            %sigCov = zeros(nAnn,nAnn);
            sigCov = tmp/nsnp;
        elseif opts.hom ==0
            sigCov = zeros(nAnn,nAnn,nComb);
            for g=1:nComb
                sigCov(:,:,g) = tmp/nsnp;
            end
        end
        clear tmp;
    elseif opts.constraintMu == 0
        mu = opts.initQ1*ones(nComb,nAnn);
        tmp = 0;
        for g=1:nComb          %mu: nComb-by-nAnn
            tmp = tmp + (annMat - repmat(mu(g,:),nsnp,1))'*(annMat - repmat(mu(g,:),nsnp,1));
        end
        if opts.hom == 1 
            %sigCov = zeros(nAnn,nAnn);
            sigCov = tmp/nsnp;
        elseif opts.hom ==0
            sigCov = zeros(nAnn,nAnn,nComb);
            for g=1:nComb
                sigCov(:,:,g) = tmp/nsnp;
            end
        end
        clear tmp;
    end
end

fprintf('EM starts...\n');
%% EM step

piMat = zeros(opts.maxIters, nComb); % prepare a matrix to track the pi vector;
piMat(1,:) = pi_vec;
betaAlphaMat = zeros(opts.maxIters,nGWAS); %prepare a matrix to track the alpha value of beta distribution;
betaAlphaMat(1,:) = betaAlpha;

if (opts.empiricalNull == 1)
    betaAlpha0Mat = zeros(opts.maxIters,nGWAS);
    betaAlpha0 = ones(1,nGWAS);
    betaAlpha0Mat(1,:) = betaAlpha0;
end

loglik = zeros(opts.maxIters,1);
loglik(1) = -inf;

maxDiff_vec = zeros(opts.maxIters,1);% recode changes in parameter estimate.
maxDiff_vec(1) = inf;

% if opts.verbose>=1
%     pi
% end

if (~isempty(annMat))
    if opts.constraintMu == 1
        muMat = zeros(1, nAnn, opts.maxIters) ;
        muMat(:,:,1) = mu;
        if  opts.hom == 1 
            sigMat = zeros(nAnn,nAnn, opts.maxIters);
            sigMat(:,:,1) = sigCov;
        elseif opts.hom ==0
            sigMat = zeros(nAnn,nAnn, nComb,opts.maxIters);
            sigMat(:,:,:,1) = sigCov;
        end
    elseif opts.constraintMu == 0
        muMat = zeros(nComb, nAnn, opts.maxIters) ;
        muMat(:,:,1) = mu;
        if opts.hom == 1 
            sigMat = zeros(nAnn,nAnn, opts.maxIters);
            sigMat(:,:,1) = sigCov;
        elseif opts.hom == 0
            sigMat = zeros(nAnn,nAnn, nComb,opts.maxIters);
            sigMat(:,:,:,1) = sigCov;
        end
    end
end

%pre-allocated matrice for effecient computation
llmat = zeros(nsnp,nGWAS); %for track loglik
Z = zeros(nsnp,nComb);  %
betaDistr = zeros(nsnp,nGWAS);
beta0Distr = zeros(nsnp,nGWAS);
gaussianDistr = zeros(nsnp,nComb);

Zmarg = cell(nGWAS,1);
sumZmargList = cell(nGWAS,1);
for k = 1:nGWAS
    Zmarg{k} = zeros(nsnp,2);
end

if (opts.empiricalNull == 1)
    Zmarg0 = cell(nGWAS,1);
    sumZmarg0List = cell(nGWAS,1);
    for k = 1:nGWAS
        Zmarg0{k} = zeros(nsnp,2);
    end
end
%% main EM  iteration

for iter = 2: opts.maxIters
    %pre-calculate Beta dist and gaussian dist
    for k = 1:nGWAS
        betaDistr(:,k) = betaAlpha(k).*gwasPval(:,k).^( betaAlpha(k)-1); %betapdf(gwasPval(:,k), betaAlpha(k), 1);% 
        if (opts.empiricalNull == 1)
            beta0Distr(:,k)= betaAlpha0(k).*gwasPval(:,k).^( betaAlpha0(k)-1);
        end
    end
    
    if (~isempty(annMat))
        if opts.hom == 1 
            for g = 1:nComb
                if opts.constraintMu == 1
                    gaussianDistr(:,g) = mvnpdf(annMat,mu,sigCov);
                elseif opts.constraintMu == 0 
                    gaussianDistr(:,g) = mvnpdf(annMat,mu(g,:),sigCov);
                end
            end
        else
            for g = 1:nComb
                if opts.constraintMu == 1
                    gaussianDistr(:,g) = mvnpdf(annMat,mu,sigCov(:,:,g));
                elseif opts.constraintMu == 0
                    gaussianDistr(:,g) = mvnpdf(annMat,mu(g,:),sigCov(:,:,g));
                end
            end            
        end
    end

    % E step
    for g = 1:nComb
        
        %mixing proportion
        Z(:,g) = pi_vec(g);
        for k=1:nGWAS
            if binaryMat(g,k)==1
                Z(:,g)=Z(:,g).*betaDistr(:,k);
            end
            if binaryMat(g,k)==0
                if (opts.empiricalNull == 1)
                    Z(:,g)=Z(:,g).*beta0Distr(:,k);
                end
            end
        end
        
        %annotation
        if (~isempty(annMat))
           Z(:,g) = Z(:,g).* gaussianDistr(:,g);
        end
    end
    
    %normalize Z
    denom = sum(Z,2);
    for g=1:nComb
        Z(:,g) = Z(:,g)./denom;
    end
    
    %marginal Z
    for k=1:nGWAS
        Zmarg{k}(:,2) = sum(Z(:,binaryMat(:,k)==1),2);
        Zmarg{k}(:,1) = 1-Zmarg{k}(:,2);
    end
    
    if (opts.empiricalNull == 1)
        for k=1:nGWAS
            Zmarg0{k}(:,2) = sum(Z(:,binaryMat(:,k)==0),2);
            Zmarg0{k}(:,1) = 1-Zmarg0{k}(:,2);
        end
    end
    %M step:
    
    %pre-computation
    sumZvec = sum(Z,1);
    for k = 1:nGWAS
        sumZmargList{k} = sum(Zmarg{k});
    end

    if (opts.empiricalNull == 1)
        for k = 1:nGWAS
            sumZmarg0List{k} = sum(Zmarg0{k});
        end
    end
    %M-step: update Beta parameters
    for k = 1:nGWAS
        betaAlpha(k) = sumZmargList{k}(2)/sum(Zmarg{k}(:,2).*(-logPval(:,k)));
        
        if (betaAlpha(k)<=opts.lbBetaAlpha)
            betaAlpha(k)=opts.lbBetaAlpha;
        end
        
        
        if (betaAlpha(k)>= 1-opts.lbBetaAlpha)
            betaAlpha(k)=1-opts.lbBetaAlpha;
        end
        
        if (opts.empiricalNull == 1)
            betaAlpha0(k) = sumZmarg0List{k}(2)/sum(Zmarg0{k}(:,2).*(-logPval(:,k)));
        end
    end
    
    %M-Step: update pi_vec
    pi_vec = sumZvec/nsnp;
    if opts.constraintPi == 1
        tmp_pi = zeros(nGWAS,1);
        for k = 1:nGWAS
            tmp_pi(k) = sum(pi_vec(binaryMat(:,k) == 1));
        end
        for g = 1:nComb
            marg_pi = zeros(nGWAS,1);
            for k = 1:nGWAS
                if ( binaryMat(g,k) ==1 )
                    marg_pi(k) = tmp_pi(k);
                else
                    marg_pi(k) = 1 - tmp_pi(k);
                end
            end
            pi_vec(g) = prod(marg_pi);
        end
        clear marg_pi tmp_pi;
    end
    
    
    
    %M-Step: update gaussian parameters
    if (~isempty(annMat))
        if opts.constraintMu == 1
            for d = 1:nAnn
                mu(1,d) = sum(annMat(:,d))./nsnp;
            end
            if opts.hom == 1
                tmp = 0;
                for g=1:nComb          %mu: nComb-by-nAnn
                    tmp = tmp + ((annMat - repmat(mu,nsnp,1)).*repmat(Z(:,g),1,nAnn).^0.5)'*((annMat - repmat(mu,nsnp,1)).*repmat(Z(:,g),1,nAnn).^0.5);
                end
                sigCov = tmp/nsnp;
                clear tmp;
            elseif opts.hom == 0
                for g=1:nComb          %mu: nComb-by-nAnn
                    sigCov(:,:,g) = ((annMat - repmat(mu,nsnp,1)).*repmat(Z(:,g),1,nAnn).^0.5)'*((annMat - repmat(mu,nsnp,1)).*repmat(Z(:,g),1,nAnn).^0.5)/sumZvec(g);
                end
            end

        elseif opts.constraintMu == 0
            if opts.hom == 1
                tmp = 0;
                for g=1:nComb          %mu: nComb-by-nAnn
                    for d = 1:nAnn
                        mu(g,d) = sum(Z(:,g).*annMat(:,d))/sumZvec(g);
                    end
                    tmp = tmp + ((annMat - repmat(mu(g,:),nsnp,1)).*repmat(Z(:,g),1,nAnn).^0.5)'*((annMat - repmat(mu(g,:),nsnp,1)).*repmat(Z(:,g),1,nAnn).^0.5);
                end
                sigCov = tmp/nsnp;
                clear tmp;
            elseif opts.hom == 0
                for g=1:nComb          %mu: nComb-by-nAnn
                    for d = 1:nAnn
                        mu(g,d) = sum(Z(:,g).*annMat(:,d))./sumZvec(g);
                    end
                    sigCov(:,:,g) = ((annMat - repmat(mu(g,:),nsnp,1)).*repmat(Z(:,g),1,nAnn).^0.5)'*((annMat - repmat(mu(g,:),nsnp,1)).*repmat(Z(:,g),1,nAnn).^0.5)/sumZvec(g);
                end
            end
        end
    end
    
    %check the lower bound: pi
    if any(pi_vec<opts.lbPi1)
        pi_vec(pi_vec<opts.lbPi1) = opts.lbPi1;
        pi_vec = pi_vec/sum(pi_vec);
    end
    
    %check the lower bound for q1
    %if (~isempty(annMat))
    %    q1(q1<opts.lbQ1) = opts.lbQ1;
    %    q1(q1>1-opts.lbQ1) = 1-opts.lbQ1;
    %end
    
    % recode estimate in each iter
    piMat(iter,:) = pi_vec;
    betaAlphaMat(iter,:) = betaAlpha;
    
    if (opts.empiricalNull == 1)
        betaAlpha0Mat(iter,:) = betaAlpha0;
    end
    
    if (~isempty(annMat) )
        muMat(:,:,iter) = mu;
        if (opts.hom==1)
            sigMat(:,:,iter) = sigCov;
        else 
            sigMat(:,:,:,iter) = sigCov;
        end
    end
    
    %track loglik and estimates
    %pre-calculate Beta distr
    for k = 1:nGWAS
        betaDistr(:,k) = betaAlpha(k).*gwasPval(:,k).^( betaAlpha(k)-1); %betapdf(gwasPval(:,k), betaAlpha(k), 1);% 
        if (opts.empiricalNull == 1)
            beta0Distr(:,k)= betaAlpha0(k).*gwasPval(:,k).^( betaAlpha0(k)-1);
        end
    end

    if (~isempty(annMat))
        if opts.constraintMu == 1
            if  opts.hom == 1
                for g = 1:nComb
                    gaussianDistr(:,g) = mvnpdf(annMat,mu,sigCov);
                end
            elseif opts.hom == 0
                for g = 1:nComb
                    gaussianDistr(:,g) = mvnpdf(annMat,mu,sigCov(:,:,g));
                end
            end
        elseif opts.constraintMu == 0
            if opts.hom == 1
                for g = 1:nComb
                    gaussianDistr(:,g) = mvnpdf(annMat,mu(g,:),sigCov);
                end
            elseif opts.hom == 0
                for g = 1:nComb
                    gaussianDistr(:,g) = mvnpdf(annMat,mu(g,:),sigCov(:,:,g));
                end
            end
        end
    end

    % track complete loglik
    for g = 1:nComb
        llmat(:,g) = pi_vec(g);
        
        % emission for GWAS
        for k = 1:nGWAS
            if (binaryMat(g,k) == 1)
                llmat(:,g) = llmat(:,g).*betaDistr(:,k); 

            end
            if (opts.empiricalNull == 1 && binaryMat(g,k) == 0)
                    llmat(:,g) = llmat(:,g).*beta0Distr(:,k); 
            end 
        end
        
        % if annotation exist
        
        if (~isempty(annMat))
            llmat(:,g) = llmat(:,g).*gaussianDistr(:,g);
            %for d=1:nAnn
            %    tmp = [1-q1(d,g); q1(d,g)];
            %    llmat(:,g) = llmat(:,g).*tmp(annMat(:,d)+1);
            %end
        end
    end
    
    loglik(iter) = sum(log(sum(llmat,2)));
    
    if opts.verbose >=1
        if mod(iter,10)==0
            fprintf('%d-th iter: loglik %f\n', iter, loglik(iter));
        end
    end
    
    % track parameter estimates
    
    piDiff = abs(piMat(iter,:) - piMat(iter-1,:))./piMat(iter,:);
    betaAlphaDiff = abs( betaAlphaMat(iter,:) - betaAlphaMat(iter-1,:))./ betaAlphaMat(iter,:);
    
    maxDiff = max([piDiff(:); betaAlphaDiff(:)]);
    
    if (~isempty(annMat))
        muDiff = (muMat(:,:,iter) - muMat(:,:,iter-1))./muMat(:,:, iter);
        maxDiff = max([maxDiff; muDiff(:)]);
    end
    
    maxDiff_vec(iter) = maxDiff;
    if (opts.verbose>=1)
        if mod(iter,10)==0
            fprintf('max changes in parameters: %f\n', maxDiff);
        end
    end
    
    % track loglik and parameter estimates
    
    if (loglik(iter) < loglik(iter-1))
%         fprintf('Info: LogLik deseases! Estimates are jittered.\n');
%         
%         pi_vec = jitter(piMat(iter-1,:));
%         pi_vec(pi_vec<=0) = 0.01;
%         pi_vec(pi_vec>=1) = 0.99;
%         pi_vec = pi_vec/sum(pi_vec);
%         piMat(iter, :) = pi_vec;
%         
%         betaAlpha = jitter(betaAlphaMat(iter-1,:));
%         betaAlpha(betaAlpha<=0) = 0.01;
%         betaAlpha(betaAlpha>=1) = 0.99;
%         betaAlphaMat(iter,:) = betaAlpha;
%         
%         if (~isempty(annMat))
%             mu = jitter(muMat(:,:,iter-1));
%             %q1 = jitter(q1Mat(:,:,iter-1));
%             %q1(q1<=0) = 0.01;
%             %q1(q1>=1) = 0.99;
%             %q1Mat(:,:,iter) = q1;
%         end
        
    else
         % stop criterion
         if (loglik(iter) - loglik(iter-1) < opts.epsStopLogLik)
                piMat((iter+1):end,:) = [];
                betaAlphaMat((iter+1):end,:) = [];                
                loglik(iter+1:end) = [];
                
                if(~isempty(annMat))
                    muMat(:,:,(iter+1):end) = [];
                end
                fprintf('Info: Algorithm converges in %d iters because there is no improvements in Log likelihood.\n', iter);
                break;
             
         else
             
         
             if (maxDiff < opts.epsStopPE)
                fprintf('Info: Algorithm converges in %d iters because there is no improvements in parameter estimates.\n',iter);
                piMat((iter+1):end,:) = [];
                betaAlphaMat((iter+1):end,:) = [];
                loglik(iter+1:end) = [];
                if(~isempty(annMat))
                    muMat(:,:,(iter+1):end) = [];
                end
                break;
             end
         end
    end
    
     
    
end

if iter == opts.maxIters
    fprintf('Info: Algorithm stops because it reaches maxIters = %d.\n', opts.maxIters);
    fprintf('Max changes in parameters: %f\n', maxDiff);
    fprintf('Changes in loglik: %f\n', loglik(iter) - loglik(iter-1));
end
    obj.loglik = loglik;

    obj.pi_vec = pi_vec;
    obj.piMat = piMat;
    
    obj.betaAlpha = betaAlpha;
    obj.betaAlphaMat = betaAlphaMat;
    
    if (opts.empiricalNull == 1)
        obj.betaAlpha0 = betaAlpha0;
        obj.betaAlpha0Mat = betaAlpha0Mat;
    end
    
    obj.Z = Z;
    obj.Zmarg = Zmarg;
    
    obj.maxDiff = maxDiff_vec;
    
    if (~isempty(annMat))
        obj.mu = mu;
        obj.muMat = muMat;
        obj.sigCov = sigCov;
        obj.sigMat = sigMat;
    end
   


