%% DEIM_vs_GP_Demo
clear
% Solve y(x,mu)=(1-x)cos(3*pi*mu*(x+1))*e^(-(1+x)*mu);

X=linspace(-1,1,100)';
Mu=linspace(1,3*pi,80)';

for i=1:length(Mu)
    Y_train(:,i)=(1-X).*cos(3*pi*Mu(i)*(X+1)).*exp(-(1+X)*Mu(i));   
end

%% DEIM 
[U,S,V]=svd(Y_train);

[phi,Uu,P] = DEIM(U);

Mu_test=linspace(1,3*pi,1000)';
for i=1:length(Mu_test)
    Y(:,i)=(1-X).*cos(3*pi*Mu_test(i)*(X+1)).*exp(-(1+X)*Mu_test(i));   
end

Dim_deim=20:20:100;

% for j =1:length(Dim_deim)
%     Pr=P(:,1:Dim_deim(j));
%     Ur=U(:,1:Dim_deim(j));
%     Dr=Ur*inv(Pr'*Ur);
%     Y_approx= Dr*(Pr'*Y);
%     SE=(Y-Y_approx).^2;
%     MSSE_DEIM(j)=mean(SE(:));   
% end

i=1;
for j =Dim_deim
    Pr=P(:,1:j);
    Ur=U(:,1:j);
    Dr=Ur*inv(Pr'*Ur);
    Y_approx= Dr*(Pr'*Y);
    SE=(Y-Y_approx).^2;
    MSSE_DEIM(i)=mean(SE(:));   
    i=i+1;
end



%% GP Direct prediction 
% Z=S*V';
Z=U'*Y_train;
Z=Z';

% GP Structure
Dim_X=1; % only a temporary value as reminder. Would be detect later by the script.
covfunc = {@covSum,{@covSEard,@covNoise}}; 
hyp.cov = [zeros(Dim_X+1,1);0];

likfunc = @likGauss; 
sn = 0.1;
hyp.lik = log(sn);

meanfunc=[];
hyp.mean=[];


InfMethod=@infExact;


[num_Z,dim_Z]=size(Z);     
for i=1:dim_Z 
    %Clean memory of hyperparameter
%     [num_X,dim_X]=size(X);
    hyp.cov = [zeros(Dim_X+1,1);0];
    sn = 0.1;
    hyp.lik = log(sn);
    hyp.mean=[];

    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, Mu, Z(:,i));
    exp(hyp.lik);
    nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, Mu, Z(:,i));        
    [m(:,i) s(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, Mu, Z(:,i), Mu_test);
end

for i=1:dim_Z 
    Y_approx= U(:,1:i)*m(:,1:i)';
    SE=(Y-Y_approx).^2;
    MSSE_GP(i)=mean(SE(:));   
    
end


figure
plot(Dim_deim,log(MSSE_DEIM),'*-')
hold on 
plot([1:100],log(MSSE_GP),'.-')
hold off
title('MSSE VS Dimension of DEIM')

%% Indirect GP




