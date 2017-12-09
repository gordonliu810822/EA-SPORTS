clear;
% ngene : number of genes
% ntiss : number of tissues
% nnon  : number of nonzero coefficients in tissue
% snr   : signal to noise ratio
% alpha : simulated alpha value for non-null distribution
% pi0   : proportion of pvalues from null distriubtion
% rho   : covariance among two studies.
ngene = 20000; ntiss = 100; nnon = 5; snr = 2; alpha = 0.4; pi0 = 0.9; rho = 0;


[pvalue,z,Anno] = generativeModel2(ngene,ntiss,nnon, snr,rho, pi0, alpha);

z1 = z(:,1);
z2 = z(:,2);
%indicator for true status (null or non-null) for each pvalue in each study

opts1.hom = 0;
options1 = GPAgaussSet(opts1);
fm11_nSum = GPAgauss(pvalue,[],options1); %joint analysis using no tissues
fm10_nSum = GPAgauss(pvalue(:,1),[],options1);%separate analysis using no tissues
fm01_nSum = GPAgauss(pvalue(:,2),[],options1);%separate analysis using no tissues

% using posterior prob as response to select tissues
opts.nfold = 5;
opts.maxLam = 2;
opts.nlam = 100;
options = cvSparseLDASet(opts);
out = cvSparseLDA2(Anno,fm11_nSum.Z,options);
lamopt = out.bestLam1se;
obj = sparseLDA(Anno,fm11_nSum.Z, out.bestK, lamopt);

% joint analysis using tissues
opts2.hom = 0;
options2 = GPAgaussSet(opts2);
summ = Anno*obj.discr;
if ( sum(sum(summ) ~= 0) )
    summ = summ(:,sum(summ) ~= 0);
    
    fm11_wSum = GPAgauss(pvalue,summ,options2);
else
    fm11_wSum = fm11_nSum;
end

% seperate analysis using tissues
out10 = cvSparseLDA2(Anno,fm10_nSum.Z,options);
lamopt = out10.bestLam1se;
obj10 = sparseLDA(Anno,fm10_nSum.Z, 1, lamopt);

summ10 = Anno*obj10.discr;
if ( sum(summ10)~=0 )
    fm10_wSum = GPAgauss(pvalue(:,1),summ10,options2);
else
    fm10_wSum = fm10_nSum;
end

% seperate analysis using tissues
out01 = cvSparseLDA2(Anno,fm01_nSum.Z,options);
lamopt = out01.bestLam1se;
obj01 = sparseLDA(Anno,fm01_nSum.Z, 1, lamopt);

summ01 = Anno*obj01.discr;
if ( sum(summ01)~=0 )
    fm01_wSum = GPAgauss(pvalue(:,2),summ01,options2);
else
    fm01_wSum = fm01_nSum;
end

