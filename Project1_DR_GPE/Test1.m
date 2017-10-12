% TEST 1
%% Test on GPML toolbox On 1-in-1out data 

clear all

%% Generate toy data
n=10;
limt_up=50;
limt_down=-50;
% x=linspace(limt_down, limt_up, n)';
x = rand(n,1)*(limt_up-limt_down)+limt_down;

z=linspace(limt_down, limt_up, 100)';
%  plot(x,y,'+');

% y=exp(x)+sin(x);
y=2.*x+x.^2+sin(x)*50;
yz=2.*z+z.^2+sin(z)*50;

%% GP
% %Structure 1
% meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]); % could be working EXTREMELY well
% likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

%Structure 2
covfunc = @covSEiso; 
hyp.cov = [0; 0];

likfunc = @likGauss; 
sn = 0.1;
hyp.lik = log(sn);

meanfunc=[];
hyp.mean=[];




% TRAIN & TEST EXACT
% hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, y);
% exp(hyp.lik)
% nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y)
% [m s] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, z);



% MCMC
cov = {@covSum,{@covSEiso,@covNoise}}; hyp.cov = [0; 2; -Inf];       % ell,sf,sn
lik =  {'likLogistic'};  hyp.lik  = [];                  % likLogistic or likErf
mn  = {'meanZero'};      hyp.mean = [];
inf = 'infEP';

% 3) predict using approximated inference
% [nlZ,dnlZ,post] = gp(hyp, inf, mn, cov, lik, xtr, ytr);
% [ymu,ys2,fmu,fs2,junk,post] = gp(hyp,inf,mn,cov,lik,xtr,post,xte);

% 4) prediction using MCMC sampling
% set MCMC parameters, see some more details in inf/infMCMC.m
% We have two samplers implemented, namely
%  hmc - Hybrid Monte Carlo, and
%  ess - Elliptical Slice Sampling.
% par.sampler = 'hmc'; par.Nsample = 20;
par.sampler = 'ess'; % par.Nais = 5; par.Nsample = 200; par.Nskip = 100;
tic
[posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,x,y,par);
toc
[m,s,fmus,fs2s,junk,posts] = gp(hyp,@infMCMC,mn,cov,lik,x,posts,z);

% TRAIN & TEST MCMC
% hyp = minimize(hyp, @gp, -100, @infMCMC, meanfunc, covfunc, likfunc, x, y);
% exp(hyp.lik)
% nlml2 = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, x, y)
% [m s] = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, x, y, z);





figure(1)
f = [m+2*sqrt(s); flipdim(m-2*sqrt(s),1)];
fill([z; flipdim(z,1)], f, [7 7 7]/8)

hold on; plot(z, m);  plot(x, y, '+');hold off;

figure(2)
plot(x,y,'+');
hold on;
plot(z, m,'--');
plot(z, yz,'-');
hold off