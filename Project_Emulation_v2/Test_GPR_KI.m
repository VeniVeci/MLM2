%Test_GPR_KI

clear

%Assign parameter
new_dim_start=1;    
new_dim_step =1;
new_dim_end  =3;

num_train=60;                  %Defining the input parameters
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
        
        % Gaussian process with isomap
        options = struct('dim_new',new_dim,'neighbor',10,'d2p_method','Dw', 'd2p_Dwpara',10,'d2p_points',10);      
        [Y_star_isomap,Yvar_star_isomap,t_isomap]=GPR_Isomap(X,Y,X_star,options);

        %------------------------------------------------------------ 
        %Record
             
%         Y_star_kpca=real(Y_star_kpca);
        
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
                     
    end

figure(1)
boxplot(RecRateSsErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of KPCA-GPR'));

figure(2)
boxplot(RecRateSsErr_isomap(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of Isomap-GPR'));
    
    




% 
% 
% 
% clear
% 
% %Assign parameter
% new_dim_start=1;    
% new_dim_end  =5;
% 
% num_train=100;                  %Defining the input         %80
% num_test=1;
% 
% iteration=10;                   %Defining test times
% 
% %-----------------------------------------------------------------------
% %temp variable
% Rec=zeros(1,6);
% 
% %-----------------------------------------------------------------------
% 
% for i = 1:iteration  
%     
%     %------------------------------------------------------------------------------------
%     %Defining the input
% 
%     [X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);
% 
%         %X=load('Inputs2.txt');
%         %Y=load('Outputs2.txt');
%         %X_star=load('InputTest14');
%         %Y_starorig=load('OutputTest14');
% 
%     [np_train,Dim_X]=size(X);
%     [np_test, Dim_Y]=size(Y_starorig);
%     %------------------------------------------------------------------------------------------
%     % Testing------------------------------
% 
%     for new_dim=new_dim_start:new_dim_end
%         
%         %-------------------------------------
%         % Gaussian process with kernel PCA
%         
%         %options = struct('ker','poly','arg',[2,0],'new_dim',new_dim); 
%         options = struct('ker','rbf','arg',50,'new_dim',new_dim); 
%         %options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);
%         
%         [Y_star_kpca,Yvar_star_kpca,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
%         
%         
%         %-------------------------
%         % Gaussian process with isomap
%         options = struct('dim_new',new_dim,'neighbor',25,'d2p_method','Dw', 'd2p_Dwpara',10,'d2p_points',10);
%         
%         [Y_star_isomap,Yvar_star_isomap,t_isomap]=GPR_Isomap(X,Y,X_star,options);
% 
%         %------------------------------------------------------------
%         % Record
%         % KPCA
%         SquErr_kpca=(Y_starorig-Y_star_kpca).^2;
%         RateErr_kpca=(abs(Y_starorig-Y_star_kpca))./Y_starorig;
%         RecSsErr(new_dim,1)=sum(SquErr_kpca(:));
%         RecARsErr(new_dim,1)=sum(RateErr_kpca(:))/(np_test*Dim_Y); %Relative sum err
%         RecTime(new_dim,1)=t_kpca;
%         
%         % Isomap
%         SquErr_isomap=(Y_star_isomap-Y_starorig).^2;
%         RateErr_isomap=(abs(Y_star_isomap-Y_starorig))./Y_starorig;
%         RecSsErr(new_dim,2)=sum(SquErr_isomap(:));
%         RecARsErr(new_dim,2)=sum(RateErr_isomap(:))/(np_test*Dim_Y); %Relative sum err
%         RecTime(new_dim,2)=t_isomap;
% 
%         %Plot
%         figure(new_dim)
%         title(sprintf('NonRescaled New dim= %g', new_dim))
%         plot(Y_star_kpca,'-ob')
%         hold on
%         plot(Y_star_isomap,'-+r')
%         plot(Y_starorig,'-dk')
%         hold off
%         legend('Y^*--KPCA.','Y^*--isomap','Y^*--real.','Location','northeast')
%         
%     end
%     
%     Record=[RecSsErr,RecARsErr,RecTime];
%     
%     Rec=[Rec;Record];
%     
% end
%     
% Rec=removerows(Rec,'ind',1);    
% 
% %-----------------------------------------------
% % Deleting blank line. It help showing a brief result
%  temp=(sum(Rec'))';  
%  temp=temp';
%  index_void=find(temp==0);
% 
%  Rec2=removerows(Rec,'ind',index_void);    
