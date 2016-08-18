function [pvalue,z,Anno] = generativeModel2(ngene,ntiss,nnon, snr, rho, pi0,alpha)

Anno = randn(ngene, ntiss);
Anno = Anno - ones(ngene,1)*mean(Anno);
Anno = Anno./(ones(ngene,1)*std(Anno))/sqrt(ntiss);

SIGMA = [1,rho; rho,1];
betas = mvnrnd([0,0],SIGMA,ntiss);
b00 = 0.5; e00 = b00/snr;

%study 1
nonzero_ind = 1:nnon;
beta1 = zeros(ntiss,1);
beta1(nonzero_ind,:) = b00*sqrt(ntiss/nnon)*betas(nonzero_ind,1);
res1 = e00*randn(ngene,1);
signal1 = Anno*beta1;
linpred1 = signal1+ res1; % snr = sqrt((var(signal))/var(res));

qt08 = quantile(linpred1,pi0);

z1 = zeros(ngene, 1);
z1(linpred1>= qt08) = 1;

%study 2

z2 = zeros(ngene,1);
indx = 1:ngene;
pi1 = 1-pi0;

matched = int64(ngene*(pi1)*(pi1+ rho*pi0));
nmatched= sum(z1) - matched;
z2(randsample(indx(z1==1),matched)) = 1;
z2(randsample(indx(z1==0),nmatched)) = 1;

pvalue = zeros(ngene,2);
pvalue(z1 ~= 1,1) = unifrnd(0,1,length(z1(z1~=1)),1);
pvalue(z1 == 1,1) = betarnd(alpha,1,length(z1(z1==1)),1);

pvalue(z2 ~= 1,2) = unifrnd(0,1,length(z2(z2~=1)),1);
pvalue(z2 == 1,2) = betarnd(alpha,1,length(z2(z2==1)),1);

z = [z1,z2];