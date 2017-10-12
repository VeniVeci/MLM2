% TEST 2
%% Test on GPML toolbox On 2-in-1-out data 

clear all

%% Generate toy data
n=30;
limt_up=5;
limt_down=-5;
step=0.5;
% x=linspace(limt_down, limt_up, n)';
x = rand(n,2)*(limt_up-limt_down)+limt_down;

[GridX,GridY] = meshgrid(limt_down:step:limt_up,limt_down:step:limt_up);
[lx,ly]=size(GridX);

z(:,1) = reshape(GridX,[],1)';
z(:,2) = reshape(GridY,[],1)';

% z(:,1)=linspace(limt_down, limt_up, 20)';
% z(:,2)=linspace(limt_down, limt_up, 20)';

%  plot(x,y,'+');

% y=exp(x)+sin(x);
% y=sin(x(:,1))*10+sin(x(:,2))*10+x(:,1).*x(:,2);
% yz=sin(z(:,1))*10+sin(z(:,2))*10+z(:,1).*z(:,2);
% GridZ=sin(GridX)*10+sin(GridY)*10+GridX.*GridY;

%option 1
y = x(:,1) .* exp(-x(:,1).^2 - x(:,2).^2);
GridZ = GridX .* exp(-GridX.^2 - GridY.^2);

%option 1
y=sin(x(:,1))*10+sin(x(:,2))*10+x(:,1).*x(:,2);
GridZ=sin(GridX)*10+sin(GridY)*10+GridX.*GridY;


% PLOT
% surf(GridX,GridY,GridZ);
% hold on
% scatter3(x(:,1),x(:,2),y)
% hold off



% surf(z(:,1),z(:,1),yz);

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

% TRAIN & TEST MCMC
% hyp = minimize(hyp, @gp, -100, @infMCMC, meanfunc, covfunc, likfunc, x, y);
% exp(hyp.lik)
% nlml2 = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, x, y)
% [m s] = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, x, y, z);


%% TRY MCMC
% 2) set GP parameters
cov = {@covSum,{@covSEiso,@covNoise}}; hyp.cov = [0; 2; -Inf];       % ell,sf,sn
% cov = @covSEiso; hyp.cov = [0; 0; ];       % ell,sf,sn

% lik =  {'likLogistic'};  hyp.lik  = [];                  % likLogistic or likErf
lik =  {'likGauss'};                 % likLogistic or likErf
sn = 0.1;
hyp.lik = log(sn);
mn  = {'meanZero'};      hyp.mean = [];
% inf = 'infEP';

% 3) predict using approximated inference
% [nlZ,dnlZ,post] = gp(hyp, inf, mn, cov, lik, xtr, ytr);
% [ymu,ys2,fmu,fs2,junk,post] = gp(hyp,inf,mn,cov,lik,xtr,post,xte);

% 4) prediction using MCMC sampling
% set MCMC parameters, see some more details in inf/infMCMC.m
% We have two samplers implemented, namely
%  hmc - Hybrid Monte Carlo, and
%  ess - Elliptical Slice Sampling.
par.sampler = 'hmc'; 
% par.sampler = 'ess'; 
%  par.Nsample = 20;
% par.sampler = 'ess'; % par.Nais = 5; par.Nsample = 200; par.Nskip = 100;
tic
[posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,x,y,par);
% [posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,x,y);
toc
[m,s,fmus,fs2s,junk,posts] = gp(hyp,@infMCMC,mn,cov,lik,x,posts,z);



% MCMC
% par.sampler = 'hmc'; par.Nsample = 20;
% mn  = {'meanZero'};      hyp.mean = [];
% tic
% [posts,nlZs,dnlZs] = infMCMC(hyp,mn,covfunc,likfunc,x,y,par);
% toc
% [m,s,fmus,fs2s,junk,posts] = gp(hyp,@infMCMC,mn,covfunc,likfunc,x,posts,z);



% TRAIN & TEST EXACT
% hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, y);
% exp(hyp.lik)
% nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y)
% [m s] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, z);

mup=m+2*sqrt(s);
mdown=m-2*sqrt(s);
figure(1)
plot(m,'-');
hold on
plot(mup,'--');
plot(mdown,'--');
hold off

figure(2)
Gridm = reshape(m,[lx,ly]);
surf(GridX,GridY,Gridm);
hold on
scatter3(x(:,1),x(:,2),y,'or')
hold off

figure(3)
surf(GridX,GridY,GridZ);
hold on
scatter3(x(:,1),x(:,2),y,'or')
hold off


% figure(1)
% f = [m+2*sqrt(s); flipdim(m-2*sqrt(s),1)];
% fill([z; flipdim(z,1)], f, [7 7 7]/8)
% 
% hold on; plot(z, m);  plot(x, y, '+');hold off;
% 
% figure(2)
% plot(x,y,'+');
% hold on;
% plot(z, m,'--');
% plot(z, yz,'-');
% hold off