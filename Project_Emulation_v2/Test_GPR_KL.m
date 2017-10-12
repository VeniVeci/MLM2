%Test_GPR_KL

clear

%Assign parameter
new_dim_start=3;    
new_dim_step =1;
new_dim_end  =3;

num_train=50;                  %Defining the input parameters
num_test=100;

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
        
        %options = struct('ker','poly','arg',[2,0],'new_dim',new_dim); 
        options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
%         %options = struct('ker','rbf','arg',50,'new_dim',new_dim); 
%         options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);
        
        [Y_star_kpca,Yvar_star_kpca,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
        [Y_star_svd,Yvar_star_svd,t_svd]=GPR_SVD(X,Y,X_star,new_dim);

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
  
        %SVD part
        SquErr=(Y_starorig-Y_star_svd).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
        RateErr=abs(Y_starorig-Y_star_svd)./Y_starorig;
        RateErr=mean(RateErr,2);
        
        RecSSErr_svd(:,new_dim)   =SSErr;
        RecRateSsErr_svd(:,new_dim)=RatSsErr;
        RecRateErr_svd(:,new_dim)=RateErr;
        RecTime_svd(:,new_dim)=t_svd;
                     
    end
    
%     save('RecSSErr_kpca_nptrain100.txt','RecSSErr_kpca','-ascii')
%     save('RecSSErr_svd_nptrain100.txt','RecSSErr_svd','-ascii')
    
figure(1)
boxplot(RecRateSsErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of KPCA-GPR'));

figure(2)
boxplot(RecRateSsErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of LPCA-GPR'));

figure(3)
boxplot(RecRateErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Average Error Rate of Each pixel of KPCA-GPR'));

figure(4)
boxplot(RecRateErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Average Error Rate of Each pixel of LPCA-GPR'));



% n=99;
% 
% Mass=reshape(Y_starorig(n,:),50,50);
% figure(3)
% contourf(Mass)
% title(sprintf('Real'));
% 
% Mass=reshape(Y_star_kpca(n,:),50,50);
% figure(4)
% contourf(Mass)
% title(sprintf('KPCA-GPR'));
% 
% Mass=reshape(Y_star_svd(n,:),50,50);
% figure(5)
% contourf(Mass)
% title(sprintf('LPCA-GPR'));

% figure(3)
% boxplot(RecRateErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
% title(sprintf('Average Error Rate of Each pixel of KPCA-GPR'));
% 
% figure(4)
% boxplot(RecRateErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
% title(sprintf('Average Error Rate of Each pixel of LPCA-GPR'));
    
    % A boxplot example
%     data = rand(20,24)
%     month = repmat({'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec'},1,2);
%     simobs = [repmat({'sim'},1,12),repmat({'obs'},1,12)];
%     boxplot(data,{month,simobs},'colors',repmat('rb',1,12),'factorgap',[5 2],'labelverbosity','minor');
    
%         %Plot
%         figure(new_dim)
%         title(sprintf('NonRescaled New dim= %g', new_dim))
%         plot(Y_star_kpca,'-ob')
%         hold on
%         plot(Y_star_svd,'-+r')
%         plot(Y_starorig,'-dk')
%         hold off
%         legend('Y^*-KPCA.','Y^*--LPCA','Y^*--real.','Location','northeast')


%         Mass=reshape(Y_star_svd,50,50);
%         Mass=Mass';
%         figure21=figure('InvertHardcopy','off','Color',[1 1 1]);
%         axes1 = axes('Parent',figure21,'FontSize',28,'FontName','Times');        
%         box(axes1,'on');
%         hold(axes1,'all');
%         contourf(x2,x1,Mass,'Parent',axes1)
%         % Create xlabel
%         xlabel('x / m','FontSize',28,'FontName','Times');
%         % Create ylabel
%         ylabel('y / m','FontSize',28,'FontName','Times');
%         title(sprintf('(b)'),'FontSize',28,'FontName','Times')
%         colorbar('peer',axes1,'FontSize',28,'FontName','Times')