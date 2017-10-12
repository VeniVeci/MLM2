%Test_GPR_KIkI

clear

%Assign parameter
new_dim_start=1;    
new_dim_step =1;
new_dim_end  =3;

num_train=150;                  %Defining the input parameters
num_test=300;

%------------------------------------------------------------------------------------
%Defining the input

[X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);

% Y=Y*10000;
% Y_starorig=Y_starorig*10000;

%------------------------------------------------------------------------------------------
% Testing------------------------------

    for new_dim=new_dim_start:new_dim_step:new_dim_end
        
%       options = struct('ker','poly','arg',[2,0],'new_dim',new_dim); 
        options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
%       options = struct('ker','rbf','arg',50,'new_dim',new_dim); 
%       options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);
        
        [Y_star_kpca,Yvar_star_kpca,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
%       [Y_star_svd,Yvar_star_svd,t_svd]=GPR_SVD(X,Y,X_star,new_dim);
        
        % Gaussian process with isomap and kernel isomap
        options = struct('dim_new',new_dim,'neighbor',50,'d2p_method','Dw', 'd2p_Dwpara',3,'d2p_points',50);      
        [Y_star_isomap,Yvar_star_isomap,t_isomap]=GPR_Isomap(X,Y,X_star,options);
        [Y_star_kisomap,Yvar_star_kisomap,t_kisomap]=GPR_kIsomap(X,Y,X_star,options);

        %------------------------------------------------------------ 
        %Record
             
%       Y_star_kpca=real(Y_star_kpca);
        
        %KPCA part
        SquErr=(Y_starorig-Y_star_kpca).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
        RateErr=abs(Y_starorig-Y_star_kpca)./Y_starorig;
        RateErr=mean(RateErr,2);
        
        RecSSErr_kpca(:,new_dim)   =SSErr;
        RecRateSsErr_kpca(:,new_dim)=RatSsErr;
        RecRateErr_kpca(:,new_dim)=RateErr;
        RecTime_kpca(:,new_dim)=t_kpca;
  
        %Isomap part
        SquErr=(Y_starorig-Y_star_isomap).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
        RateErr=abs(Y_starorig-Y_star_isomap)./Y_starorig;
        RateErr=mean(RateErr,2);
        
        RecSSErr_isomap(:,new_dim)   =SSErr;
        RecRateSsErr_isomap(:,new_dim)=RatSsErr;
        RecRateErr_isomap(:,new_dim)=RateErr;
        RecTime_isomap(:,new_dim)=t_isomap;
        
        %kIsomap part
        SquErr=(Y_starorig-Y_star_kisomap).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
        RateErr=abs(Y_starorig-Y_star_kisomap)./Y_starorig;
        RateErr=mean(RateErr,2);
        
        RecSSErr_kisomap(:,new_dim)   =SSErr;
        RecRateSsErr_kisomap(:,new_dim)=RatSsErr;
        RecRateErr_kisomap(:,new_dim)=RateErr;
        RecTime_kisomap(:,new_dim)=t_kisomap;
                     
    end

figure(1)
boxplot(RecRateSsErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of KPCA-GPR'));

figure(2)
boxplot(RecRateSsErr_isomap(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of Isomap-GPR'));
  
figure(3)
boxplot(RecRateSsErr_kisomap(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of Kernel Isomap-GPR'));
    
temp=[RecRateSsErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),RecRateSsErr_isomap(:,new_dim_start:new_dim_step:new_dim_end)];
figure(4)
boxplot(temp);
