% Test_GPR_KL_v4
%
% This script is written to track down the fail prombel when PCA_GPR is
% applied to the data set: "5050PoleFlow"
% 20-Nov-2013 

clear

%Assign parameter

new_dim=3;
num_test=300;

num_train_star=50;
num_train_step=50;
num_train_end=200;

%Assign memory for recoder
SSE=[];

%% PCA-GPR test
%----------------------
i=1;

for num_train = num_train_star:num_train_step:num_train_end
    
    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    %------------------------------------------------------------------------------------------
    % Testing------------------------------
    
    [Y_star_svd,Yvar_star_svd,t_svd]=GPR_SVD(X,Y,X_star,new_dim);
    options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
    [Y_star_kpca,Yvar_star_kpca,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
     

    SquErr=(Y_starorig-Y_star_kpca).^2;
    means=mean(Y_starorig,2);
    SSErr=sum(SquErr,2);
    RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
    RateErr=abs(Y_starorig-Y_star_kpca)./Y_starorig;
    RateErr=mean(RateErr,2);

    RecSSErr_kpca(:,i)   =SSErr;
    RecRateSsErr_kpca(:,i)=RatSsErr;
    RecRateErr_kpca(:,i)=RateErr;
    RecTime_kpca(:,i)=t_kpca;

    %SVD part
    SquErr=(Y_starorig-Y_star_svd).^2;
    means=mean(Y_starorig,2);
    SSErr=sum(SquErr,2);
    RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
    RateErr=abs(Y_starorig-Y_star_svd)./Y_starorig;
    RateErr=mean(RateErr,2);

    RecSSErr_svd(:,i)   =SSErr;
    RecRateSsErr_svd(:,i)=RatSsErr;
    RecRateErr_svd(:,i)=RateErr;
    RecTime_svd(:,i)=t_svd;

    i=i+1;
    
end

figure(1)
boxplot(RecRateSsErr_kpca(:,1:1:i-1),{1:1:i-1});
title(sprintf('Square Sum Error Rate of Each pixel of KPCA-GPR'));
figure(2)
boxplot(RecRateSsErr_svd(:,1:1:i-1),{1:1:i-1});
title(sprintf('Square Sum Error Rate of Each pixel of LPCA-GPR'));

RecRateSsErr_comb=[];
for j=1:i-1
    RecRateSsErr_comb=[RecRateSsErr_comb,RecRateSsErr_kpca(:,j),RecRateSsErr_svd(:,j)]; 
end
figure(3)
boxplot(RecRateSsErr_comb);
title(sprintf('Square Sum Error Rate of Each pixel of KLPCA-GPR'));




figure(4)
boxplot(RecRateErr_kpca(:,1:1:i-1),{1:1:i-1});
title(sprintf('Square Sum Error Rate of Each pixel of KPCA-GPR'));
figure(5)
boxplot(RecRateErr_svd(:,1:1:i-1),{1:1:i-1});
title(sprintf('Square Sum Error Rate of Each pixel of LPCA-GPR'));

RecRateErr_comb=[];
for j=1:i-1
    RecRateErr_comb=[RecRateErr_comb,RecRateErr_kpca(:,j),RecRateErr_svd(:,j)]; 
end
figure(6)
boxplot(RecRateErr_comb);
title(sprintf('Square Sum Error Rate of Each pixel of KLPCA-GPR'));


